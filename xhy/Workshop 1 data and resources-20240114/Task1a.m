% Task1a_a
timeIdx = 1;

% 定义经纬度和高度
latitude_deg = -33.821075;   % 纬度，单位为度
longitude_deg = 151.188496; % 经度，单位为度
height = 120;               % 高度，单位为米

% 使用 deg_to_rad 变量将经纬度从度转换为弧度
latitude_rad = latitude_deg * deg_to_rad;
longitude_rad = longitude_deg * deg_to_rad;

% 设定速度为零
velocity_ned = [0; 0; 0];

% 调用 pv_NED_to_ECEF 函数
[r_eb_e, v_eb_e] = pv_NED_to_ECEF(latitude_rad, longitude_rad, height, velocity_ned);

% 显示结果
disp('ECEF 坐标：');
disp(r_eb_e);

% Task1a_b
% 读取伪距数据文件
data = csvread('Workshop2_Pseudo_ranges.csv'); % 忽略第一行和第一列的标题

% 获取时间列和卫星编号行
satellite_numbers = data(1, 2:end);% 2 17 18 22...
% 获取卫星数量
num_satellites = size(satellite_numbers, 2);% 8
satellite_positions = zeros(3, num_satellites); % 8
% 在时间0时刻计算每颗卫星的ECEF位置
idx = 1;
for j = satellite_numbers
    [sat_r_es_e, ~] = Satellite_position_and_velocity(60 * (timeIdx - 1), j);
    satellite_positions(:, idx) = sat_r_es_e;
    fprintf('卫星%d的ECEF位置：X=%f m, Y=%f m, Z=%f m\n', j, sat_r_es_e(1), sat_r_es_e(2), sat_r_es_e(3));
    idx = idx + 1;
end
% Task1a_c
% 初始化距离数组
satellite_ranges = zeros(1, num_satellites);

% 计算每颗卫星的距离
idx = 1;
for j = satellite_numbers
    % 计算卫星到用户的相对位置向量
    % 初始化Sagnac效应补偿矩阵为单位矩阵
    C_e_prime = eye(3);
    % r_ea_e = C_e_prime * satellite_positions(:, idx) - r_eb_e;
    r_ea_e = satellite_positions(:, idx) - r_eb_e;
  
    % 初始距离估计
    range_old = 0;
    range_new = sqrt(r_ea_e' * r_ea_e);
    
    % 递归计算，直到新旧距离之差的绝对值小于某个阈值（例如1e-6米）
    while abs(range_new - range_old) > 1e-6
        % 更新距离估计
        range_old = range_new;
    
        % 计算Sagnac效应补偿矩阵C'_e，根据方程（2）
        C_e_prime = [1, omega_ie * range_old / c, 0;
        -omega_ie * range_old / c, 1, 0;
        0, 0, 1];
            % 应用Sagnac效应补偿矩阵来修正相对位置向量
        % r_ea_e_corrected = C_e_prime * range_old * satellite_positions(:, idx)- r_eb_e;
        r_ea_e_corrected = C_e_prime * satellite_positions(:, idx)- r_eb_e;
    
        % 根据修正后的向量重新计算距离
        range_new = sqrt(r_ea_e_corrected' * r_ea_e_corrected);
    end

    % 存储收敛后的距离
    satellite_ranges(idx) = range_new;  
    % 打印每个卫星的预测距离
    fprintf('卫星%d的预测距离：%f m\n', j, satellite_ranges(idx));
    idx = idx + 1;
end

% Task1a_d
% 计算每颗卫星的视线单位向量
line_of_sight_vectors = zeros(3, num_satellites);
idx = 1;
for j = satellite_numbers
    r_ea_e_corrected = C_e_prime * satellite_positions(:, idx) - r_eb_e;
    line_of_sight_vectors(:, idx) = r_ea_e_corrected / satellite_ranges(idx);
    fprintf('卫星%d的视线单位向量：[%f, %f, %f]\n', j, line_of_sight_vectors(1, idx), line_of_sight_vectors(2, idx), line_of_sight_vectors(3, idx));
    idx = idx + 1;
end


% Task1a_e
% 预测状态向量和测量创新向量初始化
predicted_state_vector = zeros(4, 1);
measurement_innovation = zeros(num_satellites, 1);
H_G = zeros(num_satellites, 4);

% 假设初始预测接收机时钟偏移为0
predicted_receiver_clock_offset = 0;

% 计算预测状态向量 x̂ 和测量创新向量 δz̃
idx = 1;
for j = satellite_numbers
    predicted_state_vector(1:3) = r_eb_e;
    predicted_state_vector(4) = predicted_receiver_clock_offset;

    % 预测的伪距
    % rho_j_hat = norm(satellite_positions(:, idx) - r_eb_e);  
    % 因为题目写的是measured所以其使用的是下面的
    rho_j_hat = data(timeIdx + 1, idx + 1);
    disp(rho_j_hat);
    % 测量创新
    measurement_innovation(idx) = rho_j_hat - satellite_ranges(idx) - predicted_receiver_clock_offset;

    % 测量矩阵 H^e_G 的构造
    H_G(idx, 1:3) = -line_of_sight_vectors(:, idx)';
    H_G(idx, 4) = 1;
    idx = idx + 1;
end


% Task1a_f
% 利用无加权最小二乘计算位置和接收机时钟偏移
% 最小二乘解

delta_x = inv(H_G' * H_G) * H_G' * measurement_innovation;
estimated_position_and_clock_offset = predicted_state_vector + delta_x;

% 打印计算得到的位置和时钟偏移
disp('估计的位置和接收机时钟偏移：');
disp(estimated_position_and_clock_offset);


% Task1a_f
[L_b, lambda_b, h_b, v_eb_n] = pv_ECEF_to_NED(estimated_position_and_clock_offset(1:3), estimated_position_and_clock_offset(4));


% 显示结果
fprintf('g题结果：[%f, %f, %f]\n',  L_b / deg_to_rad, lambda_b / deg_to_rad, h_b);