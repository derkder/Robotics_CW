% 主函数
% Task0
ts = 0.5; % 传播间隔为0.5秒
acceleration_north = 4; % 北向加速度为4 m/s²

% 调用函数并打印结果
pos_error = position_error(ts, acceleration_north);
fprintf('位置误差为: %.2f 米\n', pos_error);

H = measurement_matrix();
disp('测量矩阵 H:');
disp(H);



% Task1
%a
% 估计的初始位置和速度
estimated_position = [2447019; -5884199; -284783]; % (m)
estimated_velocity = [184; 77; 0]; % (m/s)

% 初始不确定性值
position_uncertainty = 10; % (m)
velocity_uncertainty = 5; % (m/s)

% 初始状态向量 X0
X0 = [estimated_position; estimated_velocity];

% 初始误差协方差矩阵 P0
P0 = diag([position_uncertainty^2; position_uncertainty^2; position_uncertainty^2;
            velocity_uncertainty^2; velocity_uncertainty^2; velocity_uncertainty^2]);

% 输出初始状态向量和初始误差协方差矩阵
disp('初始状态向量 X0:');
disp(X0);
disp('初始误差协方差矩阵 P0:');
disp(P0);

%b
ts = 1;
Phi = [eye(3), ts * eye(3);
        zeros(3,3), eye(3)];

% c
s_ae = 5;
Q =  [(1 / 3) * s_ae * (ts ^ 3) * eye(3), (1 / 2) * s_ae * (ts ^ 2) * eye(3);
        (1 / 2) * s_ae * (ts ^ 2) * eye(3), s_ae * ts * eye(3)];

%d
x_k = Phi * X0;
disp(x_k);
%e
p_k = Phi * P0 * Phi' + Q;
disp(p_k);

%f
H = [1, 0, 0, 0, 0, 0;
        0, 1, 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0];
%g
phi_x = 2.5;
R_k = [phi_x ^ 2, 0, 0;
        0, phi_x ^ 2, 0;
        0, 0, phi_x ^ 2];
%h
K_k= p_k * H' * inv(H * p_k * H' + R_k);
disp(K_k);
%i
data = csvread('Workshop2_GNSS_Pos_ECEF.csv');
r_ea_k = data(1, 2 : end);
z_k = r_ea_k' - x_k(1 : 3);
disp(z_k);
%j
x_k = x_k + K_k * z_k;
disp(x_k);
%k
freedom = 6;
p_k = (eye(freedom) - K_k * H) * p_k;
disp(p_k);
%i
[L_b, lambda_b, h_b, v_eb_n] = pv_ECEF_to_NED(x_k(1 : 3), x_k(4 : 6));
fprintf('最终结果：纬度 = %f°, 经度 = %f°, 高度 = %f米\n, 速度 = %f米%f米%f米\n', L_b / deg_to_rad, lambda_b / deg_to_rad, h_b, v_eb_n(1), v_eb_n(2), v_eb_n(3));

% 状态转移矩阵函数
function Phi = state_transition_matrix(ts)
    Phi = [1, 0, ts, 0, 0.5*ts^2, 0;
           0, 1, 0, ts, 0, 0.5*ts^2;
           0, 0, 1, 0, ts, 0;
           0, 0, 0, 1, 0, ts;
           0, 0, 0, 0, 1, 0;
           0, 0, 0, 0, 0, 1];
end

% 位置误差函数
function pos_error = position_error(ts, acceleration_north)
    % Define the state transition matrix without acceleration
    Phi_no_accel = state_transition_matrix(ts);
    % disp(Phi_no_accel)
    % Initial state vector (assuming initial velocities and accelerations are zero)
    X0 = [0; 0; 0; 0; acceleration_north; 0];
    
    % Compute the final state by multiplying the initial state with the transition matrix
    X_final = Phi_no_accel * X0;
    
    % The position error is the difference in the north position
    pos_error = X_final(1);
end

% 测量矩阵函数
function H = measurement_matrix()
    H = [1, 0, 0, 0, 0, 0; % North position
         0, 1, 0, 0, 0, 0; % East position
         0, 0, 1, 0, 0, 0; % North velocity
         0, 0, 0, 1, 0, 0]; % East velocity
end
