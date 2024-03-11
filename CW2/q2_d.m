% This script runs Q2(d)

% Q2d
% Set this value to truncate the run at a specified timestep rather than
% run through the whole simulation to its end.
maxStepNumbers = [1205, 1210];
% Initialize arrays to store determinants for each iteration
determinantsAll = zeros(length(maxStepNumbers), 7);



for stepIndex = 1:length(maxStepNumbers)
    % Create the configuration object.
    configuration = drivebot.SimulatorConfiguration();
    % Enable the laser to support pure SLAM
    configuration.enableLaser = true;    
    % For this part of the coursework, this should be set to false.
    configuration.perturbWithNoise = false;
    configuration.maximumStepNumber = maxStepNumbers(stepIndex);
    % Set up the simulator
    simulator = drivebot.DriveBotSimulator(configuration, 'q2_d');

    % Create the localization system
    drivebotSLAMSystem = drivebot.DriveBotSLAMSystem(configuration);
    drivebotSLAMSystem.setRecommendOptimizationPeriod(inf);
    
    % This tells the SLAM system to do a very detailed check that the input
    % appears to be correct but can make the code run slowly. Once you are
    % confident your code is working safely, you can set this to false.
    drivebotSLAMSystem.setValidateGraph(false);
    
    % Run the main loop and correct results
    results = minislam.mainLoop(simulator, drivebotSLAMSystem);
    
    % Minimal output plots. For your answers, please provide titles and label
    % the axes.
    
    % Plot optimisation times.
    minislam.graphics.FigureManager.getFigure('Optimization times');
    clf
    plot(results{1}.optimizationTimes, '*')
    hold on
    
    % Plot the error curves.
    minislam.graphics.FigureManager.getFigure('Errors');
    clf
    plot(results{1}.vehicleStateHistory'-results{1}.vehicleStateHistory')
    
    % Plot covariance.
    minislam.graphics.FigureManager.getFigure('Vehicle Covariances');
    clf
    plot(results{1}.vehicleCovarianceHistory')
    hold on
    
    % Plot errors.
    minislam.graphics.FigureManager.getFigure('Errors');
    clf
    plot(results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory')
    hold on
    
    % Plot chi2 values.
    minislam.graphics.FigureManager.getFigure('chi2 values');
    clf
    plot(results{1}.chi2Time, results{1}.chi2History)
    hold on
    
    % Get landmark estimates after the simulation
    [x, P, landmarkIds] = drivebotSLAMSystem.landmarkEstimates();

    % Calculate the determinants of the covariance matrices for each landmark
    for i = 1:length(landmarkIds)
        determinantsAll(stepIndex, i) = det(P(:,:,i));
    end

    % Print the differences in determinants for each landmark
    for i = 1:length(landmarkIds)
        fprintf('Landmark %d: Determinant = %f\n', landmarkIds(i), determinantsAll(stepIndex, i));
    end
    
end

% Calculate the differences in determinants between the two iterations for each landmark
determinantsDifferences = diff(determinantsAll, 1);

% Print the differences in determinants for each landmark
for i = 1:length(landmarkIds)
    fprintf('Landmark %d: Determinant Difference = %f\n', landmarkIds(i), determinantsDifferences(i));
end


