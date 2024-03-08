% This script runs Q2(c)

% Create the configuration object.
configuration = drivebot.SimulatorConfiguration();

% Enable the laser to support pure SLAM
configuration.enableLaser = true;

% If you set this parameter to false, the simulator generates measurements
% with no noise in them. You might find this useful for debugging.
% However, unless specified otherwise, any submitted results must have this
% value set to true.
configuration.perturbWithNoise = true;

% Set up the simulator
simulator = drivebot.DriveBotSimulator(configuration, 'q2_c');

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

% Plot optimisation times
minislam.graphics.FigureManager.getFigure('Optimization times');
clf
plot(results{1}.optimizationTimes, '*')
hold on

% Plot the error curves
minislam.graphics.FigureManager.getFigure('Errors');
clf
plot(results{1}.vehicleStateHistory'-results{1}.vehicleStateHistory')

% Plot covariance
minislam.graphics.FigureManager.getFigure('Vehicle Covariances');
clf
plot(results{1}.vehicleCovarianceHistory')
hold on

% Plot errors
minislam.graphics.FigureManager.getFigure('Errors');
clf
plot(results{1}.vehicleStateHistory'-results{1}.vehicleTrueStateHistory')
hold on

% Plot chi2 values
minislam.graphics.FigureManager.getFigure('chi2 values');
clf
plot(results{1}.chi2Time, results{1}.chi2History)
hold on


% This is how to extract the graph from the optimizer
graph = drivebotSLAMSystem.optimizer();

% This is how to extract cell arrays of the vertices and edges from the
% graph
allVertices = graph.vertices();
%disp(length(allVertices));
allEdges = graph.edges();

% Work out the number of vehicle poses and landmarks. 
numVehicleVertices = 0;
numLandmarks = 0;

landmarkObservationsPerVehicleVertex = 0;
observationsPerLandmarkVertex = 0;

% Q2c:
% Finish implementing the code to capture information about the graph
% structure.
% 初始化统计变量
totalLandmarkObservations = 0;
landmarkObservationCounts = containers.Map('KeyType', 'int64', 'ValueType', 'int32');

% 遍历所有顶点来分类统计
for vertex = allVertices
    if isa(vertex{1}, 'drivebot.graph.VehicleStateVertex')  % 假设你有一个表示车辆姿态的类
        numVehicleVertices = numVehicleVertices + 1;
        %totalLandmarkObservations = totalLandmarkObservations + vertex{1}.numberOfEdges();
    elseif isa(vertex{1}, 'drivebot.graph.LandmarkStateVertex')  % 同样，假设你有一个表示地标的类
        numLandmarks = numLandmarks + 1;
        totalLandmarkObservations = totalLandmarkObservations + vertex{1}.numberOfEdges();
        landmarkObservationCounts(vertex{1}.id()) = vertex{1}.numberOfEdges();
    end
end

% 计算平均观测数
% 可能因为一个timestamp对应一个process vertice吧
averageLandmarkObservationsPerVehicle = totalLandmarkObservations / numVehicleVertices;
averageObservationsPerLandmark = mean(cell2mat(values(landmarkObservationCounts)));

% 输出结果
fprintf('Number of vehicle poses: %d\n', numVehicleVertices);
fprintf('Number of landmarks initialized: %d\n', numLandmarks);
fprintf('Average number of observations made by a robot at each timestep: %f\n', averageLandmarkObservationsPerVehicle);
fprintf('Average number of observations received by each landmark: %f\n', averageObservationsPerLandmark);
