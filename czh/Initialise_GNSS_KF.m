function [x_est,err_cov_matrix] = Initialise_GNSS_KF
%Initialise GNSS EKF state estimates and error covariance matrix 
%
% Outputs:
%   x_est                 Kalman filter estimates:
%     Rows 1-3            estimated ECEF user position (m)
%     Rows 4-6            estimated ECEF user velocity (m/s)
%     Row 7               estimated receiver clock offset (m) 
%     Row 8               estimated receiver clock drift (m/s)
%   P_matrix              state estimation error covariance matrix


% Initialize state estimates  of the machine at time 0 using the least square method
pseudo_range = csvread("Pseudo_ranges.csv");
range_rate = csvread("Pseudo_range_rates.csv");
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
Omega_ie = Skew_symmetric([0,0,omega_ie]);
c = 299792458; % Speed of light in m/s

[rows , columns] = size(pseudo_range);
Sat_num_array = pseudo_range(1,2 : columns);
Sat_number = length(Sat_num_array);
user_position = [0;0;0];
user_velocity = [0;0;0];
clock_offset = 0;
clock_drift = 0;
disp("Cartesian ECEF position of the satellites at time ");
positionMatrix = zeros(3, Sat_number);
velocityMatrix = zeros(3, Sat_number);
% Cartesian ECEF satellite position 
for idx = 1 : Sat_number
    [sat_position, sat_velocity] = Satellite_position_and_velocity(0, Sat_num_array(idx)); % ECEF position
    positionMatrix(:, idx) = sat_position;
    velocityMatrix(:, idx) = sat_velocity;
end  
% Position determination
positionDiff = 100;   
predicted_Pos_state_vector = zeros(4, 1); 
while norm(positionDiff) > 0.1
    % Predict the ranges from the approximate user position to each satellite 
    identity_matrix = eye(3);
    Range_matrix = zeros(Sat_number, 1);
    Sagnac_effect = identity_matrix * positionMatrix - user_position; % 3 ✖ 8 
    initial_range = diag(sqrt(Sagnac_effect' * Sagnac_effect)); 
    Sagnac_matrix = zeros(Sat_number, 3, 3);
    for idx = 1 : Sat_number
        Sagnac_matrix(idx, : , :) = eye(3);
        Sagnac_matrix(idx, 1, 2) = omega_ie * initial_range(idx)/ c;
        Sagnac_matrix(idx, 2, 1) = -omega_ie * initial_range(idx)/ c;
        Sagnac_effect = reshape(Sagnac_matrix(idx, : , :), 3, 3) * positionMatrix(:, idx) - user_position;
        Range_matrix(idx,:) = sqrt(Sagnac_effect' * Sagnac_effect); 
    end 
    
    % Compute the line-of-sight unit vector from the approximate user position to each
    lineOfSight_vecto_matrix = zeros(Sat_number, 3);
    for idx = 1 : Sat_number
        lineOfSight_vecto_matrix(idx,:) = (reshape(Sagnac_matrix(idx, : , :), 3, 3) * positionMatrix(:, idx) - user_position) / initial_range(idx); 
    end 
   
    % Formulate the predicted state vector, measurement innovation vector, measurement matrix
    
    predicted_Pos_state_vector(1 : 3, 1) = user_position;
    predicted_Pos_state_vector(4, 1) = clock_offset; 


    range_innovation_vector = pseudo_range(2 , 2 : columns)' - Range_matrix - clock_offset;
   
    measuremen_matrix = ones(Sat_number,4);
    measuremen_matrix(:, 1 : 3) = -lineOfSight_vecto_matrix;
 
    % Compute the position and receiver clock offset
    Updated_Pos_state_vector = predicted_Pos_state_vector + (measuremen_matrix' * measuremen_matrix) \ measuremen_matrix' * range_innovation_vector; % 4 * 1
    
    % Convert this Cartesian ECEF user position solution to latitude, longitude and height
    Updated_user_position = Updated_Pos_state_vector(1 : 3, :);
    clock_offset = Updated_Pos_state_vector(4, :);

    positionDiff = Updated_user_position - user_position;
    user_position = Updated_user_position;
end    

% Velocity determination
predicted_Vel_state_vector = zeros(4, 1); 
predicted_Vel_state_vector(1 : 3, 1) = user_velocity;
predictd_range_rates = zeros(Sat_number, 1);
for idx = 1 : Sat_number
    predictd_range_rates(idx,:) = lineOfSight_vecto_matrix(idx,:) * (reshape(Sagnac_matrix(idx, : , :), 3, 3) * (velocityMatrix(:, idx) + Omega_ie * positionMatrix(:, idx)) - (user_velocity + Omega_ie * user_position));
end 
rate_innovation_vector = range_rate(2 , 2 : columns)' - predictd_range_rates - clock_drift;
Updated_velocity_vector = predicted_Vel_state_vector + (measuremen_matrix' * measuremen_matrix) \ measuremen_matrix' * rate_innovation_vector; % 4 * 1
Updated_user_velocity = Updated_velocity_vector(1 : 3, :);
clock_drift = Updated_velocity_vector(4, :);
[latitude, longtitude, height, velocity] = pv_ECEF_to_NED(Updated_user_position,Updated_user_velocity);
disp([rad2deg(latitude), rad2deg(longtitude), height]);

x_est = [  Updated_user_position(1); Updated_user_position(2); Updated_user_position(3);
    velocity(1); velocity(2);  velocity(3);
    clock_offset; clock_drift]; 

%Initialise error covariance matrix
clock_offset_std_dev = 100000; % Initial receiver clock offset standard deviation
clock_drift_std_dev = 200; % Receiver clock drift standard deviation
err_cov_matrix =  zeros(8);
err_cov_matrix(1,1) = 100;
err_cov_matrix(2,2) = 100;
err_cov_matrix(3,3) = 100;
err_cov_matrix(4,4) = 0.01;
err_cov_matrix(5,5) = 0.01;
err_cov_matrix(6,6) = 0.01;
err_cov_matrix(7,7) = clock_offset_std_dev^2; 
err_cov_matrix(8,8) = clock_drift_std_dev^2;

% Ends