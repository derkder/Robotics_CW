clc;
clear;
%Workshop 3.
%Task1 Car Dead Reckoning

%Load Data
%Car odometry and Magnetic compass data
speed_heading_path = "Workshop3_Speed_Heading.csv";
speed_heading_data = table2array(readtable(speed_heading_path));
time = speed_heading_data(:, 1);
v_ave = speed_heading_data(:, 2);
heading_inst = speed_heading_data(:, 3);

num_epochs = length(time);

%Define Constants
h = 37.4; %Geodetic height, meter

%Convert headings to radians value
heading_inst = deg2rad(heading_inst);

%Start position latitude and longitude
lat_0 = 50.4249580; %degree
lon_0 = -3.5957974;

%Compute the average northen and easten velosity between epochs k and k-1
%Initialisation
v_ave_NE = zeros(2, num_epochs);
%Epoch 0
v_ave_NE(:, 1) = [v_ave(1) * cos(heading_inst(1)); v_ave(1) * sin(heading_inst(1))];
for k = 2 : num_epochs
    psi_k = heading_inst(k);
    psi_k_minus = heading_inst(k-1);

    temp = [cos(psi_k) + cos(psi_k_minus); sin(psi_k) + sin(psi_k_minus)];
    v_ave_NE(:, k) = 0.5 * temp * v_ave(k);
end

%Update latitude and Longitude at each epoch
%Inisialisation
lat_DR = zeros(num_epochs, 1);
lon_DR = zeros(num_epochs, 1);
lat_DR(1) = deg2rad(lat_0);
lon_DR(1) = deg2rad(lon_0);

for k = 2 : num_epochs
    %Compute radii of curvature
    [R_Nk, R_Ek] = Radii_of_curvature(lat_DR(k-1));
    %Time interval
    delta_t = time(k) - time(k-1);
    %Update latitude
        lat_DR(k) = lat_DR(k-1) + v_ave_NE(1, k) * delta_t / (R_Nk + h);
    %Update longitude
    lon_DR(k) = lon_DR(k-1) + v_ave_NE(2, k) * delta_t / ((R_Ek + h) * cos(lat_DR(k)));
end

%Convert results back to degree
% Lat_DR = rad2deg(Lat_DR);
% Lon_DR = rad2deg(Lon_DR);

%Compute the damped instantaneous DR velosity at each epoch
%Inisialisation
v_inst_N = zeros(num_epochs, 1);
v_inst_E = zeros(num_epochs, 1);
v_inst_N(1) = v_ave_NE(1, 1);
v_inst_E(1) = v_ave_NE(2, 1);

for k = 2 : num_epochs
    v_inst_N(k) = 1.7 * v_ave_NE(1, k) - 0.7 * v_inst_N(k-1);
    v_inst_E(k) = 1.7 * v_ave_NE(2, k) - 0.7 * v_inst_E(k-1);
end

%Check
DR_results = [time, lat_DR, lon_DR, v_inst_N, v_inst_E];

%Task 2 Car DR/GNSS Integration

%Load Data
gnss_pos_vel_NED_path = "Workshop3_GNSS_Pos_Vel_NED.csv";
gnss_pos_vel_NED_data = table2array(readtable(gnss_pos_vel_NED_path));
lat_gnss = gnss_pos_vel_NED_data(:, 2); %degree
lon_gnss = gnss_pos_vel_NED_data(:, 3); %degree
h_gnss = gnss_pos_vel_NED_data(:, 4);
v_gnss_N = gnss_pos_vel_NED_data(:, 5);
v_gnss_E = gnss_pos_vel_NED_data(:, 6);
v_gnss_d = gnss_pos_vel_NED_data(:, 7);

%Convert GNSS latitude and longitude to radiant
lat_gnss = deg2rad(lat_gnss);
lon_gnss = deg2rad(lon_gnss);

%Time interval
tau_s = 0.5;

%Initialise states
x_initial = zeros(4, 1);

sigma_v = 0.1;
sigma_pos = 10;
[R_N0, R_E0] = Radii_of_curvature(lat_DR(1));
%Initialise state estimation error covariance matrix
P_initial = zeros(4);
P_initial(1, 1) = sigma_v^2;
P_initial(2, 2) = sigma_v^2;
P_initial(3, 3) = sigma_pos^2 / (R_N0 + h_gnss(1))^2;
P_initial(4, 4) = sigma_pos^2 / (R_E0 + h_gnss(1))^2 / (cos(lat_DR(1))^2);

%Compute measurement matrix
H_k = [0 0 -1 0;
       0 0 0 -1;
      -1 0 0 0;
       0 -1 0 0];

%DR velosity PSD
S_DR = 0.2;
%GNSS position measurement error
sigma_pos_gnss = 5;
%GNSS velosity measurement error
sigma_v_gnss = 0.02;

%Initialise save matrix
lat_correct = zeros(num_epochs, 1);
lon_correct = zeros(num_epochs, 1);
v_N_correct = zeros(num_epochs, 1);
v_E_correct = zeros(num_epochs, 1);

%The first epoch
xk_est_prior = x_initial;
Pk_prior = P_initial;
%Kalman Filtering Update
for k = 2 : num_epochs
    %Compute radii of curvature
    [R_N, R_E] = Radii_of_curvature(lat_DR(k-1));
    %Transition matrix
     Phi_k_minus = [1 0 0 0;
             0 1 0 0;
             tau_s / (R_N + h_gnss(k-1)) 0 1 0;
             0 tau_s/((R_E + h_gnss(k-1)) * cos(lat_DR(k-1))) 0 1];
     %System noise covariance matrix
     %velosity error
     var_v = S_DR * tau_s;
     %position error
     var_pos1 = S_DR * tau_s^3 / (3 * (R_N + h_gnss(k-1))^2); %Latitude
     var_pos2 = S_DR * tau_s^3 / (3 * (R_E + h_gnss(k-1))^2 * cos(lat_DR(k-1))^2); %Longitude
     %correlation error of velosity and position
     covar1 = S_DR * tau_s^2 / (2 * (R_N + h_gnss(k-1)));
     covar2 = S_DR * tau_s^2 / (2 * (R_E + h_gnss(k-1)) * cos(lat_DR(k-1)));

     Q_k_minus = [var_v 0 covar1 0;
                  0 var_v 0 covar2;
                  covar1 0 var_pos1 0;
                  0 covar2 0 var_pos2];
     %Predict stage
     %Propagate the state estimates:
     xk_est_prior = Phi_k_minus * xk_est_prior;
     %Propagate the error covariance matrix
     Pk_prior = Phi_k_minus * Pk_prior * Phi_k_minus' + Q_k_minus;

     %Compute the measurement noise covariance matrix
     var_pos_gnss1 = sigma_pos_gnss^2 / (R_N + h_gnss(k))^2; %latitude
     var_pos_gnss2 = sigma_pos_gnss^2 / (R_E + h_gnss(k))^2 / cos(lat_gnss(k)); %longitude
     var_v_gnss = sigma_v_gnss^2;
     R_k = diag([var_pos_gnss1, var_pos_gnss2, var_v_gnss, var_v_gnss]);

     %Compute the Kalman gain
     K_k = Pk_prior * H_k' / (H_k * Pk_prior * H_k' + R_k);

     %Formulate the measurement innovation vector
     delta_z = [lat_gnss(k) - lat_DR(k) + xk_est_prior(3);
                lon_gnss(k) - lon_DR(k) + xk_est_prior(4);
                v_gnss_N(k) - v_inst_N(k) + xk_est_prior(1);
                v_gnss_E(k) - v_inst_E(k) + xk_est_prior(2)];

     %Correction stage
     %Update state estimate
     xk_est_posterior = xk_est_prior + K_k * delta_z;
     %Update error covariance matrix
     Pk_posterior = (eye(4) - K_k * H_k) * Pk_prior;

     %Use the KF estimates to correct the DR solution
     latk_correct = lat_DR(k) - xk_est_posterior(3);
     lonk_correct = lon_DR(k) - xk_est_posterior(4);
     vk_N_correct = v_inst_N(k) - xk_est_posterior(1);
     vk_E_correct = v_inst_E(k) - xk_est_posterior(2);
     
     %Save data
     lat_correct(k) = latk_correct;
     lon_correct(k) = lonk_correct;
     v_N_correct(k) = vk_N_correct;
     v_E_correct(k) = vk_E_correct;
     %Update for next epoch
     xk_est_prior = xk_est_posterior;
     Pk_prior = Pk_posterior;     
end

%Convert corrected results to degree
lat_correct = rad2deg(lat_correct);
lon_correct = rad2deg(lon_correct);

%Check the results
corrected_results = [time, lat_correct, lon_correct, v_N_correct, v_E_correct];