% 读取CSV文件
data = csvread('Workshop3_Speed_Heading.csv'); % 跳过标题行
time = data(:, 1); % 时间
speed = data(:, 2); % 速度
heading = data(:, 3); % 航向

% 初始化位置
L = zeros(length(time), 1);
lambda = zeros(length(time), 1);
L(1) = deg2rad(50.4249580); % 初始纬度，转换为弧度
lambda(1) = deg2rad(-3.5957974); % 初始经度，转换为弧度

% 初始化平均速度
v_N_ave = zeros(length(time), 1);
v_E_ave = zeros(length(time), 1);

% 循环计算每个时间点的位置
for k = 2:length(time)
    % 航向角度转弧度
    psi_k = heading(k) * deg_to_rad;
    psi_k_minus_1 = heading(k-1) * deg_to_rad;
    
    % 计算平均速度
    v_N_ave(k) = (cos(psi_k) + cos(psi_k_minus_1)) / 2 * speed(k);
    v_E_ave(k) = (sin(psi_k) + sin(psi_k_minus_1)) / 2 * speed(k);
    
    % 计算位置
    delta_t = time(k) - time(k-1);
    [R_N, R_E] = Radii_of_curvature(L(k-1)); % 子午圈曲率半径
    h = 37.4; % 地理高度
    
    L(k) = L(k-1) + v_N_ave(k) * delta_t / (R_N + h);
    lambda(k) = lambda(k-1) + v_E_ave(k) * delta_t / ((R_E + h) * cos(L(k)));
end

% 初始化平均速度
v_N = zeros(length(time), 1);
v_E = zeros(length(time), 1);
v_N(1) = speed(1) * cos(heading(1) * deg_to_rad);
v_E(1) = speed(1) * sin(heading(1) * deg_to_rad);

% 计算阻尼瞬时DR速度
for k = 2:length(time)
    v_N(k) = 1.7 * v_N_ave(k) - 0.7 * v_N(k-1);
    v_E(k) = 1.7 * v_E_ave(k) - 0.7 * v_E(k-1);
end

% 转换回度以便于观察
L = L / deg_to_rad;
lambda = lambda  / deg_to_rad;

disp([L, lambda, v_N, v_E])



% Task2
% 我真的是日了狗了，哪里写错了啊西巴
% 终于nmd做对了，我不是天才是什么
% 加载GNSS位置和速度数据
GNSS_data = csvread('Workshop3_GNSS_Pos_Vel_NED.csv');
time_GNSS = GNSS_data(:, 1);
latitude_GNSS = GNSS_data(:, 2);
longitude_GNSS = GNSS_data(:, 3);
height_GNSS = GNSS_data(:, 4);
velocity_N_GNSS = GNSS_data(:, 5);
velocity_E_GNSS = GNSS_data(:, 6);

% DR数据来自任务1的结果
% 这里假设DR_data是一个包含以下内容的矩阵：
% time_DR, latitude_DR, longitude_DR, height_DR, velocity_N_DR, velocity_E_DR
DR_data = [L, lambda, v_N, v_E]; % 需要用任务1的结果填充此处

% 初始化状态向量和协方差矩阵
x_hat = zeros(4, 1);
[R_N, R_E] = Radii_of_curvature(latitude_GNSS(1));
P_hat = diag([(0.1)^2, (0.1)^2, (10)^2 / (R_N + height_GNSS(1))^2, (10)^2 / ((R_E + height_GNSS(1)) * cos(latitude_GNSS(1)))^2]);


% 系统噪声PSD
S_DR = 0.2;

% 测量噪声
sigma_r = 5; % GNSS位置测量误差标准差，单位：米
sigma_v = 0.02; % GNSS速度测量误差标准差，单位：米/秒

% 卡尔曼滤波器的时间间隔
tau_s = 0.5;

% 5. 计算测量矩阵
H_k = [0 0 -1 0;
       0 0 0 -1;
      -1 0 0 0;
       0 -1 0 0];
% k = 1
fprintf('最终结果：time：%f°\n纬度 = %f°, 经度 = %f°, 速度 = %f米%f米\n', time_GNSS(1), latitude_GNSS(1), longitude_GNSS(1), velocity_N_GNSS(1), velocity_E_GNSS(1));

% 开始卡尔曼滤波循环
for k = 2:length(time_GNSS) 
    % 1. 计算状态转移矩阵
    [R_N, R_E] = Radii_of_curvature(latitude_GNSS(k-1));
        
    Phi_k = [1 0 0 0;
             0 1 0 0;
             tau_s/(R_N + height_GNSS(k-1)) 0 1 0;
             0 tau_s/((R_E + height_GNSS(k-1))*cos(latitude_GNSS(k-1))) 0 1];
    
    % 2. 计算系统噪声协方差矩阵
    Q_k = [S_DR*tau_s, 0, (S_DR*tau_s^2)/(2*(R_N + height_GNSS(k-1))), 0;
       0, S_DR*tau_s, 0, (S_DR*tau_s^2)/(2*(R_E + height_GNSS(k-1)) * cos(latitude_GNSS(k-1)));
       (S_DR*tau_s^2)/(2*(R_N + height_GNSS(k-1))), 0, (S_DR*tau_s^3)/(3*(R_N + height_GNSS(k-1))^2), 0;
       0, (S_DR*tau_s^2)/(2*((R_E + height_GNSS(k-1)) * cos(latitude_GNSS(k-1)))), 0, (S_DR*tau_s^3)/(3*((R_E + height_GNSS(k-1))^2 * cos(latitude_GNSS(k-1))^2))];

    
    % 3. 传播状态估计
    x_hat_minus = Phi_k * x_hat;
    
    % 4. 传播误差协方差矩阵
    P_hat_minus = Phi_k * P_hat * Phi_k' + Q_k;
    
    % 6. 计算测量噪声协方差矩阵
    R_k = [sigma_r^2/(R_N + height_GNSS(k))^2 0 0 0;
           0 sigma_r^2/((R_E + height_GNSS(k))*cos(latitude_GNSS(k)))^2 0 0;
           0 0 sigma_v^2 0;
           0 0 0 sigma_v^2];
    
    % 7. 计算卡尔曼增益矩阵
    K_k = P_hat_minus * H_k' / (H_k * P_hat_minus * H_k' + R_k);
    
    % 8. 构造测量创新向量
    delta_z_k = [latitude_GNSS(k) - DR_data(k,1);
             longitude_GNSS(k) - DR_data(k,2);
             velocity_N_GNSS(k) - DR_data(k,3);
             velocity_E_GNSS(k) - DR_data(k,4)];
    %disp(delta_z_k)
    delta_z_k = delta_z_k - H_k * x_hat_minus;
    %disp(delta_z_k)

    x_hat = x_hat_minus + K_k * delta_z_k;
    % disp(K_k * delta_z_k)
    P_hat = (eye(4) - K_k*H_k)*P_hat_minus;

    %woc我瞎改的后来发现就应该这么改，x作为变换，前两维是速度后两维才是经纬度...
    L_b = DR_data(k,1) - x_hat(3);
    lambda_b = DR_data(k,2) - x_hat(4);
    v_n = DR_data(k,3) - x_hat(1);
    v_e = DR_data(k,4) - x_hat(2);
    fprintf('最终结果：time：%f°\n纬度 = %f°, 经度 = %f°, 速度 = %f米%f米\n', time_GNSS(k), L_b, lambda_b, v_n, v_e);

    
end
