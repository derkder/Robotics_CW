clear;

%% dr 数据读取预处理与初始化
%Dead Reckoning
%Load Data
run('Define_Constants.m');

DR_data_path = "Dead_reckoning.csv";
data = table2array(readtable(DR_data_path));
time = data(:, 1);
num_epoch = length(time);
delta_t = 0.5; %time interval
v_ave_wheel = data(:, 2:5); 
gyro = data(:, 6); %rad/s
heading_mag = data(:, 7); %degree
heading_mag = deg2rad(heading_mag);%rad

% 加载GNSS位置和速度数据
GNSS_data = csvread('GNSS_Pos_Vel_NED.csv');
time_GNSS = GNSS_data(:, 1);
latitude_GNSS = GNSS_data(:, 2) * deg_to_rad;
longitude_GNSS = GNSS_data(:, 3) * deg_to_rad;
height_GNSS = GNSS_data(:, 4);
velocity_N_GNSS = GNSS_data(:, 5);
velocity_E_GNSS = GNSS_data(:, 6);

%Estimate the outputs of sensors
%Wheel odmetry
%Wheel odometry carlibration
scale_ODO = 0.03; %scale factor
v_ave_wheel = v_ave_wheel * (1 - scale_ODO);

%Notice that when the lawnmower movement direction remains almost constant, 
% the speed changes relatively smoothly. 
% There is a more noticeable change in speed when turning.

%So when the direction change slightly
%We filter out data from the wheel 
%that has a significantly different speed compared to the others.
%By observing the images, we manually selected 1.12 and 1.32 as the threshold values.
%%And when turning,we believe in the speed of the outer wheel.
v_ave_wheel(v_ave_wheel > 1.32) = 1.32;
v_ave_wheel(v_ave_wheel > 1 & v_ave_wheel < 1.12) = 1.12;

%Note that the back tires are driving tires
%And the given GNSS antenna position is not the center of the wheel
%positions. We need to compute the velosity of the antenna
a = 0.3;
b = 0.2;
prop1 = 0.2 * 0.2 / (0.5 * 0.4);
prop2 = 0.3 * 0.2 / (0.5 * 0.4);

v_ave_antenna = v_ave_wheel(:, 1) * prop1 + v_ave_wheel(:, 2) * prop1 + ...
    v_ave_wheel(:, 3) * prop2 + v_ave_wheel(:, 4) * prop2;

%% heading_gyro计算
%Use Gyro-Magnetometer Integration to estimate
%heading.
%Preprocess gyro data (Carlibration)
bia_gyro = 1; %deg/s
bia_gyro = deg2rad(bia_gyro);
scale_gyro = 0.01;
gyro = gyro * (1 - scale_gyro);

%Compute heading using gyro data 
heading_gyro = zeros(num_epoch, 1);
heading_gyro(1) = heading_mag(1);
for k = 2 : num_epoch
    heading_gyro(k) = heading_gyro(k-1) + delta_t * gyro(k);
end

%Kalman filtering based integration
%State contains gyro-derived heading error and gyro bias
x_initial = [0; 0];
sigma_gyro = 1e-4;
P_initial = diag([sigma_gyro^2 * delta_t^2, bia_gyro^2]);

%Initialisation
xk_est_prior = x_initial;
Pk_prior = P_initial;
gyro_PSD = 3 * 1e-6; %gyrp random noise PSD
%Measurement matrix
H_k = [-1, 0];
%Measurement noise matrix
sigma_mag = deg2rad(4);
R_k = sigma_mag^2;

heading_gyro_correct = zeros(num_epoch, 1);
heading_gyro_correct(1) = heading_gyro(1);

for k = 2 : num_epoch
    %Transition matrix
    Phi_k_minus = [1, delta_t;
                   0, 1];
    %System noise covariance matrix
    Q_k_minus = [gyro_PSD * delta_t, 0;
                 0, 0];
    %Prediction stage
    %Propagate the state estimates
    xk_est_prior = Phi_k_minus * xk_est_prior;
    %Propagate the error covariance matrix
    Pk_prior = Phi_k_minus * Pk_prior * Phi_k_minus' + Q_k_minus;
    %Compute Kalman gain
    K_k = Pk_prior * H_k' / (H_k * Pk_prior * H_k' + R_k);
    %Formulate measurement innovation vector
    delta_z = heading_mag(k) - heading_gyro(k) + xk_est_prior(1);

    %Correction stage
    %Update state estimate
    xk_est_posterior = xk_est_prior + K_k * delta_z;
    %Update error covariance matrix
    Pk_posterior = (eye(2) - K_k * H_k) * Pk_prior;

    heading_gyro_correct(k) = heading_gyro(k) - xk_est_posterior(1);
end

% 初始化位置
L = zeros(length(time), 1);
lambda = zeros(length(time), 1);
% 改成从GNSS读
L(1) = latitude_GNSS(1); % 初始纬度，转换为弧度
lambda(1) = longitude_GNSS(1); % 初始经度，转换为弧度

% 初始化平均速度
v_N_ave = zeros(length(time), 1);
v_E_ave = zeros(length(time), 1);


% 初始化平均速度
v_N = zeros(length(time), 1);
v_E = zeros(length(time), 1);
v_N(1) = v_ave_antenna(1) * cos(heading_gyro_correct(1));
v_E(1) = v_ave_antenna(1) * sin(heading_gyro_correct(1));

% Task2
% DR数据来自任务1的结果
% 这里假设DR_data是一个包含以下内容的矩阵：
% time_DR, latitude_DR, longitude_DR, height_DR, velocity_N_DR, velocity_E_DR
% DR_data = [L, lambda, v_N, v_E]; % 需要用任务1的结果填充此处

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
fprintf('最终结果：time：%f°\n纬度 = %f°, 经度 = %f°, 速度 = %f米%f米\n', time_GNSS(1), latitude_GNSS(1) / deg_to_rad, longitude_GNSS(1) / deg_to_rad, velocity_N_GNSS(1), velocity_E_GNSS(1));

% 开始卡尔曼滤波循环
for k = 2:length(time_GNSS) 
    % 航向角度转弧度
    psi_k = heading_gyro_correct(k);
    psi_k_minus_1 = heading_gyro_correct(k-1);
    % 计算平均速度
    v_N_ave(k) = (cos(psi_k) + cos(psi_k_minus_1)) / 2 * v_ave_antenna(k);
    v_E_ave(k) = (sin(psi_k) + sin(psi_k_minus_1)) / 2 * v_ave_antenna(k);   
    % 计算位置
    delta_t = time(k) - time(k-1);
    [R_N, R_E] = Radii_of_curvature(L(k - 1)); % 子午圈曲率半径
    h = height_GNSS(k); % 地理高度   
    L(k) = L(k-1) + v_N_ave(k) * delta_t / (R_N + h);
    lambda(k) = lambda(k-1) + v_E_ave(k) * delta_t / ((R_E + h) * cos(L(k)));
    v_N(k) = 1.7 * v_N_ave(k) - 0.7 * v_N(k-1);
    v_E(k) = 1.7 * v_E_ave(k) - 0.7 * v_E(k-1);
    DR_data(k, :) =  [L(k) * rad_to_deg, lambda(k) * rad_to_deg, v_N(k), v_E(k)];
    x_hat = zeros(4, 1);


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
    delta_z_k = [latitude_GNSS(k) - DR_data(k,1) * deg_to_rad;
             longitude_GNSS(k) - DR_data(k,2) * deg_to_rad;
             velocity_N_GNSS(k) - DR_data(k,3);
             velocity_E_GNSS(k) - DR_data(k,4)];
    %disp(delta_z_k)
    delta_z_k = delta_z_k - H_k * x_hat_minus;
    %disp(delta_z_k)

    x_hat = x_hat_minus + K_k * delta_z_k;
    % disp(K_k * delta_z_k)
    P_hat = (eye(4) - K_k*H_k)*P_hat_minus;

    %woc我瞎改的后来发现就应该这么改，x作为变换，前两维是速度后两维才是经纬度...
    L_b = DR_data(k,1) - x_hat(3) * rad_to_deg;
    lambda_b = DR_data(k,2) - x_hat(4) * rad_to_deg;
    v_n = DR_data(k,3) - x_hat(1);
    v_e = DR_data(k,4) - x_hat(2);
    fprintf('最终结果：time：%f°\n纬度 = %f°, 经度 = %f°, 速度 = %f米%f米\n', time_GNSS(k), L_b, lambda_b, v_n, v_e);

    L(k) = L_b * deg_to_rad;
    lambda(k) = lambda_b * deg_to_rad;
    v_N(k) = v_n;
    v_E(k) = v_e;

    heading_gyro_correct(k) = abs(atan2(v_E(k), v_N(k)));
end

figure;
% 绘制轨迹
plot(lambda * rad_to_deg, L * rad_to_deg, '-o'); % 经度为x轴，纬度为y轴，'-'表示线形，'o'表示数据点
% 添加标签和标题
xlabel('Longtitude (degrees)');
ylabel('Latitude (degrees)');
title('Integerate');
% 添加网格线以便更好地查看数据点
grid on;
% 可以选择设置坐标轴的限制范围，使得轨迹更清晰
xlim([min(lambda * rad_to_deg) max(lambda * rad_to_deg)]);
ylim([min(L * rad_to_deg) max(L * rad_to_deg)]);
