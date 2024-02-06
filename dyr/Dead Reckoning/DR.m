clc;
clear;
%Dead Reckoning

%Load Data
DR_data_path = "Dead_reckoning.csv";
DR_data = table2array(readtable(DR_data_path));

time = DR_data(:, 1);
num_epoch = length(time);
delta_t = 0.5; %time interval
v_ave_wheel = DR_data(:, 2:5); 

gyro = DR_data(:, 6); %rad/s
heading_mag = DR_data(:, 7); %degree
%Convert magnetic heading to radiant
heading_mag = deg2rad(heading_mag);

%Estimate the outputs of sensors
%Wheel odmetry
%Wheel odometry carlibration
scale_ODO = 0.03; %scale factor
v_ave_wheel = v_ave_wheel * (1 - scale_ODO);

%The wheels may experience slipping and turning issues.
%We plot the velosities to see if there is any outliar.
figure;
plot(time, v_ave_wheel(:, 1),'r');%Front left wheel
hold on;
plot(time, v_ave_wheel(:, 2), 'g');%Front right wheel
hold on;
plot(time, v_ave_wheel(:, 3), 'b');%Back left
hold on;
plot(time, v_ave_wheel(:, 4), 'y');%Back right
hold on;

%Gyro
%Plot the values of the heading angle to observe changes in direction.
plot(time, gyro);
hold off;

%Notice that when the lawnmower movement direction remains almost constant, 
% the speed changes relatively smoothly. 
% There is a more noticeable change in speed when turning.

%So when the direction change slightly
%We filter out data from the wheel 
%that has a significantly different speed compared to the others.
%By observing the images, we manually selected 1.12 and 1.32 as the threshold values.
%%And when turning,we believe in the speed of the outer wheel.
v_ave_wheel(v_ave_wheel > 1 & v_ave_wheel > 1.32) = 1.32;
v_ave_wheel(v_ave_wheel > 1 & v_ave_wheel < 1.12) = 1.12;

figure;
plot(time, v_ave_wheel(:, 1),'r');%Front left wheel
hold on;
plot(time, v_ave_wheel(:, 2), 'g');%Front right wheel
hold on;
plot(time, v_ave_wheel(:, 3), 'b');%Back left
hold on;
plot(time, v_ave_wheel(:, 4), 'y');%Back right
hold off;

%Note that the back tires are driving tires
%And the given GNSS antenna position is not the center of the wheel
%positions. We need to compute the velosity of the antenna
a = 0.3;
b = 0.2;
prop1 = 0.2 * 0.2 / (0.5 * 0.4);
prop2 = 0.3 * 0.2 / (0.5 * 0.4);

v_ave_antenna = v_ave_wheel(:, 1) * prop1 + v_ave_wheel(:, 2) * prop1 + ...
    v_ave_wheel(:, 3) * prop2 + v_ave_wheel(:, 4) * prop2;
figure;
plot(time, v_ave_antenna);
hold off;

%Use Gyro-Magnetometer Integration to estimate
%heading.
%Preprocess gyro data (Carlibration)
bia_gyro = 1; %deg/s
bia_gyro = deg2rad(bia_gyro);
scale_gyro = 0.01;
gyro = (gyro + bia_gyro) * (1 - scale_gyro);

%Compute heading using gyro data 
heading_gyro = zeros(num_epoch, 1);
heading_gyro(1) = heading_mag(1);
for k = 2 : num_epoch
    heading_gyro(k) = heading_gyro(k-1) + delta_t * gyro(k);
end
figure;
plot(time, heading_mag, "DisplayName","Magnetic heading");
hold on;
plot(time, heading_gyro, "DisplayName","Gyro heading");
hold off;
legend("show")
xlim([0 450])
ylim([-0.5 6])

%Kalman filtering based integration
%State contains gyro-derived heading error and gyro bias
x_initial = [0; 0];
sigma_gyro = 1e-4;
P_initial = diag([sigma_gyro^2 * delta_t^2, sigma_gyro^2]);

%Initialisation
xk_est_prior = x_initial;
Pk_prior = P_initial;

gyro_PSD = 3 * 1e-6; %gyrp random noise PSD

%Measurement matrix
H_k = [-1, 0];
%Measurement noise matrix
sigma_mag = deg2rad(4);
R_k = sigma_mag^2;

%
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

figure;
plot(time,heading_mag, "DisplayName","Magnetic heading");
hold on;
plot(time,heading_gyro, "DisplayName","Gyro heading");
hold on;
plot(time,heading_gyro_correct, "DisplayName","Mag-Gyro heading");
hold off;
legend("show")
xlim([0 450])
ylim([-0.5 6])