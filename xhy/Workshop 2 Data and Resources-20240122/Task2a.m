
%a
[x_est,P_matrix] = Initialise_GNSS_KF;
% 计算状态转移矩阵
tau_s = 1; % 秒
I3 = eye(3); % 3x3 单位矩阵
O3_1 = zeros(3, 1); % 3x1 零矩阵
O3 = zeros(3, 3); % 3x3 零矩阵

Phi_k_minus_1 = [I3, tau_s * I3, O3_1, O3_1;
                 O3, I3, O3_1, O3_1;
                 O3_1', O3_1', 1, tau_s;
                 O3_1', O3_1', 0, 1];

% 计算系统噪声协方差矩阵
S_a = 5; % 加速度功率谱密度 (PSD)
S_ca = 0.01; % 钟差加速度 PSD
S_caf = 0.04; % 钟差频率 PSD

Q_k_minus_1 = [(1/3) * tau_s^3 * S_a * I3, (1/2) * tau_s^2 * S_a * I3, O3_1, O3_1;
               (1/2) * tau_s^2 * S_a * I3, tau_s * S_a * I3, O3_1, O3_1;
               O3_1', O3_1', S_ca * tau_s + (1/3) * S_caf * tau_s^3 , (1/2) * S_caf * tau_s^2;
               O3_1', O3_1', (1/2) * S_caf * tau_s^2, S_caf * tau_s];
x_k = Phi_k_minus_1 * x_est;
disp(x_k)

P_k_minus = Phi_k_minus_1 * P_matrix * Phi_k_minus_1' + Q_k_minus_1;



% 终于对咯
% 第几个时间点
timeIdx = 1;

% 读取伪距数据文件
data = csvread('Workshop2_Pseudo_ranges.csv'); % 从第二行第一列开始读取数据
data_rates = csvread('Workshop2_Pseudo_range_rates.csv'); % 从第二行第一列开始读取数据

% 获取卫星编号
satellite_numbers = data(1, 2:end); % 假设第一行是卫星编号
num_satellites = length(satellite_numbers); % 卫星数量

% 初始化矩阵
satellite_positions = zeros(3, num_satellites); % 卫星位置
satellite_rates = zeros(3, num_satellites); % 卫星位置
satellite_ranges = zeros(1, num_satellites); % 卫星距离
line_of_sight_vectors = zeros(3, num_satellites); % 视线单位向量
H_G = zeros(num_satellites, 4); % 测量矩阵

% 初始位置假设在地球中心
estimated_position = x_k(1:3);
pseudo_range_rates = data_rates(timeIdx + 1, 2:end); % 假设第二行是伪距率 
predicted_receiver_velocity = x_k(4:6); % 从状态向量中获取用户预测速度
satellite_range_rates = zeros(num_satellites, 1);

% 在指定时间计算每颗卫星的ECEF位置
for idx = 1:num_satellites
    j = satellite_numbers(idx);
    [sat_r_es_e, sat_v_es_e] = Satellite_position_and_velocity(60 * (timeIdx - 1), j);
    satellite_positions(:, idx) = sat_r_es_e;
    satellite_rates(:, idx) = sat_v_es_e;
    % 不需要打印每颗卫星的ECEF位置，除非需要进行调试
end

% 迭代条件
convergence_threshold = 0.1; % 10cm
max_iterations = 1000; % 防止无限循环
iteration = 0; % 迭代计数器
solution_difference = inf; % 初始化差异为无限大，用于迭代
% 计算每颗卫星的距离
for idx = 1:num_satellites
    j = satellite_numbers(idx);
    % 计算卫星到用户的相对位置向量
    r_ea_e = satellite_positions(:, idx) - estimated_position;
    % 初始距离估计
    range_old = 0;
    range_new = sqrt(r_ea_e' * r_ea_e);
    % 递归计算距离，直到新旧距离之差的绝对值小于某个阈值（例如1e-6米）
    while abs(range_new - range_old) > 1e-6
        range_old = range_new;
        % 更新Sagnac效应补偿矩阵C_e_prime
        C_e_prime = [1, omega_ie * range_old / c, 0;
                     -omega_ie * range_old / c, 1, 0;
                      0, 0, 1];
        r_ea_e_corrected = C_e_prime * satellite_positions(:, idx)- estimated_position;
        range_new = sqrt(r_ea_e_corrected' * r_ea_e_corrected);
    end
    satellite_ranges(idx) = range_new; % 存储修正后的距离
    % 打印每个卫星的预测距离
    %fprintf('卫星%d的预测距离：%f m\n', j, satellite_ranges(idx));
    line_of_sight_vectors(:, idx) = r_ea_e_corrected / range_new; % 计算视线单位向量

    r_a = satellite_positions(:, idx) - estimated_position; % 卫星到接收机的相对位置向量
    v_a = satellite_rates(:, idx) - predicted_receiver_velocity; % 卫星到接收机的相对速度向量
    first = C_e_prime * (satellite_rates(:, idx) + Omega_ie * satellite_positions(:, idx));
    second = predicted_receiver_velocity + Omega_ie * x_k(1 : 3);
    satellite_range_rates(idx) = line_of_sight_vectors(:, idx)' * (first - second); % 公式(22)
end

% i
% 计算测量矩阵H_k
H_k = zeros(num_satellites * 2, 8); % 8列矩阵，因为状态向量x有8个分量
for idx = 1:num_satellites
    u_as_e = line_of_sight_vectors(:, idx);
    H_k(idx, :) = [-u_as_e', 0, 0, 0, 1, 0]; % 按照公式(23)构造H_k的每一行
    u_as_e = line_of_sight_vectors(:, idx);
    H_k(idx + num_satellites, :) = [0, 0, 0, -u_as_e', 0, 1]; % 按照公式(23)构造H_k的每一行
end

% j
% 计算测量噪声协方差矩阵R_k
sigma_p = 10; % 伪距测量误差的标准差
sigma_r = 0.05; % 伪距率测量误差的标准差
R_k = diag([sigma_p^2 * ones(1, num_satellites), sigma_r^2 * ones(1, num_satellites)]); % 公式(24)
% k
% 计算卡尔曼增益矩阵K_k
K_k = P_k_minus * H_k' / (H_k * P_k_minus * H_k' + R_k); % 公式(25


% 预测状态向量和测量创新向量初始化
measurement_innovation = zeros(num_satellites * 2, 1);
% 假设初始预测接收机时钟偏移为0
predicted_receiver_clock_offset = 0;
% 计算预测状态向量 x̂ 和测量创新向量 δz̃
idx = 1;
for j = satellite_numbers
    measurement_innovation(idx) = data(timeIdx + 1, idx + 1)- satellite_ranges(idx) - x_k(7);
    measurement_innovation(idx + num_satellites) = data_rates(timeIdx + 1, idx + 1)- satellite_range_rates(idx) - x_k(8);
    idx = idx + 1;
end

% woc尼玛呀终于作对了，就是不能乘上那个括号里的那个r_aj那个(怪题
%disp(measurement_innovation)

x_k = x_k + K_k * measurement_innovation;
P_k = (eye(8) - K_k * H_k) * P_k_minus;


