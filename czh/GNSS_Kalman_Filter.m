%  GNSS Kalman Filter Multiple Epochs
format bank
% a) Initialise the Kalman filter state vector estimate and error covariance matrix 
[estimated_state, error_covariance] = Initialise_GNSS_KF;
% b) Compute the transition matrix
interval = 0.5;
I_3 = eye(3); 
Zeros_3 = zeros(3); 
Zeors_3_1 = zeros(3, 1);  
transition_matrix = [I_3, interval * I_3, Zeors_3_1, Zeors_3_1;
                 Zeros_3, I_3, Zeors_3_1, Zeors_3_1;
                 Zeors_3_1', Zeors_3_1', 1, interval;
                 Zeors_3_1', Zeors_3_1', 0, 1]; 

% c) Compute the system noise covariance matrix  
PSD = 0.01;
Clock_phase= 0.01;
Clock_frequency = 0.04;
noise_cov_matrix = [1/3 * PSD * interval^3 * I_3, 1/2 * PSD * interval^2 * I_3, Zeors_3_1, Zeors_3_1;
               1/2 * PSD * interval^2 * I_3, PSD * interval * I_3, Zeors_3_1, Zeors_3_1;
               Zeors_3_1', Zeors_3_1', Clock_phase * interval + 1/3 * Clock_frequency * interval^3, 1/2 * Clock_frequency * interval^2;
               Zeors_3_1', Zeors_3_1', 1/2 * Clock_frequency * interval^2, Clock_frequency * interval];

pseudo_range = csvread("Pseudo_ranges.csv");
pseudo_range_rate = csvread("Pseudo_range_rates.csv");
[rows , columns] = size(pseudo_range);
Sat_num_array = pseudo_range(1,2 : columns);
time_array = pseudo_range(2 : rows, 1);
time_length = length(time_array);
latitude_list = zeros(1,time_length); 
longtitude_list = zeros(1,time_length); 

filename = 'GNSS_Pos_Vel1.csv';
csv = fopen(filename, 'w');

for time_idx = 1 :  time_length
    positionMatrix = zeros(3, Sat_number);
    velocityMatrix = zeros(3, Sat_number);
    Sat_number = length(Sat_num_array);
    disp("Time " + time_array(time_idx));
    % d)  Use the transition matrix to propagate the state estimates
    propagate_state = transition_matrix * estimated_state;

    % e) propagate the error covariance matrix
    propagate_cov = transition_matrix * error_covariance * transition_matrix' + noise_cov_matrix;

    % f) Predict the ranges from the approximate user position to each satellite 

    for idx = 1 : Sat_number
        [sat_position, sat_velocity] = Satellite_position_and_velocity(time_array(time_idx), Sat_num_array(idx));
        positionMatrix(:, idx) = sat_position;
        velocityMatrix(:, idx) = sat_velocity;
    end 
    
    Range_matrix = zeros(Sat_number, 1);
    Sagnac_effect = I_3 * positionMatrix - propagate_state(1:3); % 3 ✖ 8 
    initial_range = diag(sqrt(Sagnac_effect' * Sagnac_effect)); 
    Sagnac_matrix = zeros(Sat_number, 3, 3);
    for idx = 1 : Sat_number
        Sagnac_matrix(idx, : , :) = eye(3);
        Sagnac_matrix(idx, 1, 2) = omega_ie * initial_range(idx)/ c;
        Sagnac_matrix(idx, 2, 1) = -omega_ie * initial_range(idx)/ c;
        Sagnac_effect = reshape(Sagnac_matrix(idx, : , :), 3, 3) * positionMatrix(:, idx) - propagate_state(1:3);
        Range_matrix(idx,:) = sqrt(Sagnac_effect' * Sagnac_effect); 
    end 

    % g) Compute the line-of-sight unit vector from the approximate user position to each
    lineOfSight_vecto_matrix = zeros(Sat_number, 3);
    for idx = 1 : Sat_number
        lineOfSight_vecto_matrix(idx,:) = (reshape(Sagnac_matrix(idx, : , :), 3, 3) * positionMatrix(:, idx) - propagate_state(1:3)) / initial_range(idx); 
    end 
    
    % h) Predict the range rates from the approximate user position to each satellite
    range_rates = zeros(Sat_number, 1);
    for idx = 1 : Sat_number
        range_rates(idx,:) = lineOfSight_vecto_matrix(idx,:) * (reshape(Sagnac_matrix(idx, : , :), 3, 3) * (velocityMatrix(:, idx) + Omega_ie * positionMatrix(:, idx)) - (propagate_state(4:6,1) + Omega_ie * propagate_state(1:3,1)));
    end 

    % i) Compute the measurement matrix
    measurement_matrix = zeros(Sat_number * 2, 8);
    measurement_matrix(1:Sat_number, 1:3) = -lineOfSight_vecto_matrix; 
    measurement_matrix(Sat_number+1:Sat_number*2, 4:6) = -lineOfSight_vecto_matrix; 
    measurement_matrix(1:Sat_number, 7) = 1;
    measurement_matrix(Sat_number+1:Sat_number*2, 8) = 1;
    
    % j) Compute the measurement noise covariance matrix
    % GNSS error deviation 
    signal_error = 1; % Signal in space error standard deviation
    residual_iono = 2; % Residual ionosphere error standard deviation at zenith
    residual_tropo = 0.2; % Residual troposphere error standard deviation at zenith
    codeTrack_error = 2; % Code tracking and multipath error standard deviation
    range_rate_error = 0.02; % Range rate tracking and multipath error standard deviation
    pseudo_range_err = signal_error^2 + residual_iono^2 + residual_tropo^2 + codeTrack_error^2;
    measurement_noise_matrix = zeros(Sat_number * 2, Sat_number * 2);
    measurement_noise_matrix(1:Sat_number, 1:Sat_number) = diag(power(10,2)*ones(1,Sat_number));
    measurement_noise_matrix(Sat_number + 1 : Sat_number * 2, Sat_number + 1 : Sat_number * 2) = diag(power(0.05,2)*ones(1,Sat_number));
    
    % k) Compute the Kalman gain matrix
    Kalman_gain_matrix = propagate_cov * measurement_matrix' / (measurement_matrix * propagate_cov * measurement_matrix' + measurement_noise_matrix);
    
    % l) Formulate the measurement innovation vector
    innovation_vector = zeros(Sat_number*2 , 1);
    innovation_vector(1 : Sat_number, 1) = pseudo_range(time_idx + 1, 2:columns)' - Range_matrix(:,1) - propagate_state(7);
    innovation_vector(Sat_number + 1 : Sat_number * 2, 1) = pseudo_range_rate(time_idx + 1, 2:columns)' - range_rates(:,1) - propagate_state(8);
    
    % m) Update the state estimates
    Updated_state = propagate_state + Kalman_gain_matrix * innovation_vector;
    
    % n) Update the error covariance matrix
    Updated_error_cov = (eye(8) - Kalman_gain_matrix * measurement_matrix) * propagate_cov;

    
    % o) Convert cartesian ECEF position to latitude, longitude and height
    Updated_user_position = Updated_state(1 : 3, :);
    Updated_velocity = Updated_state(4 : 6, :);
    [new_latitude, new_longitude, new_height, new_velovcity] = pv_ECEF_to_NED(Updated_user_position,Updated_velocity);
    fprintf('   %.6f    %.6f    %.6f    %.6f    %.6f    %.6f\n', rad2deg(new_latitude), rad2deg(new_longitude), new_height, new_velovcity(1), new_velovcity(2), new_velovcity(3));

   
    % Calculate the residual covariance
    innovation_covariance = measurement_matrix * error_covariance * measurement_matrix' + measurement_noise_matrix;  

    % Outlier detection threshold
    outlier_threshold = 6; % based on chi-square distribution
    longtitude_list(time_idx) = rad2deg(new_longitude);
    latitude_list(time_idx) = rad2deg(new_latitude);
    
    reduce_innovation_vector = innovation_vector;
    while(1)
        %find the maximum residual vector 
        [maxValue, outliner_idx]  = max(abs(reduce_innovation_vector(1 : Sat_number, 1))); 
        
        sta_idx = find(abs(innovation_vector) == maxValue, 1);
        
        residual_covariance = innovation_covariance(sta_idx, sta_idx);
        standardized_residual = maxValue^2 / residual_covariance;
    
        %if the maximum residual vector is outlier  
        if standardized_residual > outlier_threshold
            % Then remove the information that provided by the 
            disp("The outlier appears on satellite " + Sat_num_array(sta_idx) + " at " + time_array(time_idx));
            positionMatrix(:, outliner_idx) = []; % remove the position of the outlier satellite
            velocityMatrix(:, outliner_idx) = []; % remove the velocity of the outlier satellite
            Sat_number = Sat_number - 1;
            reduced_columns = columns - 1;
    
            lineOfSight_vecto_matrix(outliner_idx,:) = [];
            range_rates = zeros(Sat_number, 1);
            for idx = 1 : Sat_number
                range_rates(idx,:) = lineOfSight_vecto_matrix(idx,:) * (reshape(Sagnac_matrix(idx, : , :), 3, 3) * (velocityMatrix(:, idx) + Omega_ie * positionMatrix(:, idx)) - (propagate_state(4:6,1) + Omega_ie * propagate_state(1:3,1)));
            end 
            measurement_matrix(outliner_idx, :) = [];
            measurement_matrix(Sat_number + outliner_idx, :) = [];
    
            measurement_noise_matrix = zeros(Sat_number * 2, Sat_number * 2);
            measurement_noise_matrix(1:Sat_number, 1:Sat_number) = diag(power(10,2)*ones(1,Sat_number));
            measurement_noise_matrix(Sat_number + 1 : Sat_number * 2, Sat_number + 1 : Sat_number * 2) = diag(power(0.05,2)*ones(1,Sat_number));
            Kalman_gain_matrix = propagate_cov * measurement_matrix' / (measurement_matrix * propagate_cov * measurement_matrix' + measurement_noise_matrix);
    
            reduce_innovation_vector(outliner_idx, :) = [];
            reduce_innovation_vector(Sat_number + outliner_idx, :) = [];
    
            Updated_state = propagate_state + Kalman_gain_matrix * reduce_innovation_vector;
            Updated_error_cov = (eye(8) - Kalman_gain_matrix * measurement_matrix) * propagate_cov;
    
            Updated_user_position = Updated_state(1 : 3, :);
            Updated_velocity = Updated_state(4 : 6, :);
            [new_latitude, new_longitude, new_height, new_velovcity] = pv_ECEF_to_NED(Updated_user_position,Updated_velocity);
            %fprintf(csv, '%f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', time_array(time_idx),rad2deg(new_latitude), rad2deg(new_longitude), new_height, new_velovcity(1), new_velovcity(2), new_velovcity(3));
         
            longtitude_list(time_idx) = rad2deg(new_longitude);
            latitude_list(time_idx) = rad2deg(new_latitude);
        else
            disp('No outliers detected');
            longtitude_list(time_idx) = rad2deg(new_longitude);
            latitude_list(time_idx) = rad2deg(new_latitude);
            fprintf(csv, '%f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', time_array(time_idx), rad2deg(new_latitude), rad2deg(new_longitude), new_height, new_velovcity(1), new_velovcity(2), new_velovcity(3));
            break;
        end
    end 
    estimated_state = Updated_state; 
    error_covariance = Updated_error_cov; 
end
fclose(csv);
scatter(longtitude_list, latitude_list, 5, 'MarkerFaceColor', 'b');   %  draw scatter diagram 
title('Lawmower position (PSD = 0.01)'); 
xlabel('latitude');       
ylabel('longitude');   

