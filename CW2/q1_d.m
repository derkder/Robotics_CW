% This script runs Q1(d)

% Create the configuration object.
configuration = drivebot.SimulatorConfiguration();

% Only enable the GPS
configuration.enableGPS = true;

% If you set this parameter to false, the simulator generates measurements
% with no noise in them. You might find this useful for debugging.
% However, unless specified otherwise, any submitted results for any questions or subparts of questions
% must have this value set to true.
configuration.perturbWithNoise = true;

% Set up the simulator
simulator = drivebot.DriveBotSimulator(configuration, 'q1_d');

% Create the localization system
drivebotSLAMSystem = drivebot.DriveBotSLAMSystem(configuration);

% Force the optimizer to run with this frequency. This lets you see what's
% happening in greater detail, but slows everything down.

%drivebotSLAMSystem.setRecommendOptimizationPeriod(20);
% This tells the SLAM system to do a very detailed check that the input
% appears to be correct but can make the code run slowly. Once you are
% confident your code is working safely, you can set this to false.
drivebotSLAMSystem.setValidateGraph(false);

% Run the main loop and correct results
results = minislam.mainLoop(simulator, drivebotSLAMSystem);


% Minimal output plots. For your answers, please provide titles and label
% the axes.

% Plot optimisation times
minislam.graphics.FigureManager.getFigure('Optimization times');
clf
plot(results{1}.vehicleStateTime, results{1}.optimizationTimes, '*')
hold on

% Plot the error curves
minislam.graphics.FigureManager.getFigure('Errors');
clf
plot(results{1}.vehicleStateTime, results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory')
xlabel('Time');
ylabel('State error');
legend('X','Y','Angles');

% Plot covariance
minislam.graphics.FigureManager.getFigure('Vehicle Covariances');
clf
plot(results{1}.vehicleStateTime, results{1}.vehicleCovarianceHistory');
xlabel('Time');
ylabel('Vehicle Covariance');
legend('X covariances','Y covariances','Angles covariances');
hold on

% Plot chi2 values
minislam.graphics.FigureManager.getFigure('chi2 values');
clf
plot(results{1}.chi2Time, results{1}.chi2History)
xlabel('chi2Time');
ylabel('chi2');
legend('chi2');
hold on

