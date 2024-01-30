% 第几个时间点
timeIdx = 7;

% 异常值检测阈值
T = 6;

% 测量误差标准差
sigma_rho = 5; % 5米

% 读取伪距数据文件
data = csvread('Workshop1_Pseudo_ranges.csv'); % 从第二行第一列开始读取数据

% 获取卫星编号
satellite_numbers = data(1, 2:end); % 假设第一行是卫星编号
pseudo_ranges = data(timeIdx + 1, 2:end);
num_satellites = length(satellite_numbers); % 卫星数量

% 初始化矩阵
satellite_positions = zeros(3, num_satellites); % 卫星位置
satellite_ranges = zeros(1, num_satellites); % 卫星距离
line_of_sight_vectors = zeros(3, num_satellites); % 视线单位向量
H_G = zeros(num_satellites, 4); % 测量矩阵
measurement_innovation = zeros(num_satellites, 1); % 测量创新

% 初始位置假设在地球中心
estimated_position = [0; 0; 0];
predicted_receiver_clock_offset = 0; % 初始预测的接收机时钟偏移

% 在指定时间计算每颗卫星的ECEF位置
for idx = 1:num_satellites
    j = satellite_numbers(idx);
    [sat_r_es_e, ~] = Satellite_position_and_velocity(60 * (timeIdx - 1), j);
    satellite_positions(:, idx) = sat_r_es_e;
    % 不需要打印每颗卫星的ECEF位置，除非需要进行调试
end

% 迭代条件
convergence_threshold = 0.1; % 10cm
max_iterations = 1000; % 防止无限循环
iteration = 0; % 迭代计数器
solution_difference = inf; % 初始化差异为无限大，用于迭代
hasDele = 0;

% 迭代开始
while solution_difference > convergence_threshold && iteration < max_iterations
    iteration = iteration + 1;

    % 计算每颗卫星的距离
    for idx = 1:num_satellites
        j = satellite_numbers(idx);
        % 计算卫星到用户的相对位置向量
        r_ea_e = satellite_positions(:, idx) - estimated_position;
        % 初始距离估计
        range_old = 0;
        range_new = sqrt(r_ea_e' * r_ea_e);
        % 递归计算距离，直到新旧距离之差的绝对值小于某个阈值（例如1e-6米）
        while abs(range_new - range_old) > 1e-6
            range_old = range_new;
            % 更新Sagnac效应补偿矩阵C_e_prime
            C_e_prime = [1, omega_ie * range_old / c, 0;
                         -omega_ie * range_old / c, 1, 0;
                          0, 0, 1];
            r_ea_e_corrected = C_e_prime * satellite_positions(:, idx) - estimated_position;
            range_new = sqrt(r_ea_e_corrected' * r_ea_e_corrected);
        end
        satellite_ranges(idx) = range_new; % 存储修正后的距离
        % 打印每个卫星的预测距离
        fprintf('卫星%d的预测距离：%f m\n', j, satellite_ranges(idx));
        line_of_sight_vectors(:, idx) = r_ea_e_corrected / range_new; % 计算视线单位向量
    end

    % 构建测量创新向量和测量矩阵H_G
    for idx = 1:num_satellites
        j = satellite_numbers(idx);
        % 从数据中获取伪距
        % 好像是这里错了
        rho_j_hat = pseudo_ranges(idx);
        % 测量创新
        measurement_innovation(idx) = rho_j_hat - satellite_ranges(idx) - predicted_receiver_clock_offset;
        % 构建测量矩阵H_G
        H_G(idx, :) = [-line_of_sight_vectors(:, idx)', 1];
    end

    % 最小二乘解
    delta_x = (H_G' * H_G) \ (H_G' * measurement_innovation);
    
    % 更新位置和时钟偏移
    estimated_position = estimated_position + delta_x(1:3);
    predicted_receiver_clock_offset = predicted_receiver_clock_offset + delta_x(4);

    % 计算残差向量v
    I_m = eye(num_satellites); % 创建m x m的单位矩阵
    H_G_transpose = H_G'; % H_G的转置
    %P = I_m - H_G * inv(H_G_transpose * H_G) * H_G_transpose; % 投影矩阵
    P = H_G * inv(H_G_transpose * H_G) * H_G_transpose - I_m; % 投影矩阵
    v = P * measurement_innovation; % 残差向量
    
    % 计算残差协方差矩阵C_v
    C_v = (-P) * (sigma_rho^2);
    
    % 检测异常值
    outliers = abs(v) > sqrt(diag(C_v)) * T;

    % 如果有异常值
    if any(outliers) && iteration > 2
        T = 5;
        % 找到残差最大的卫星
        [~, outlier_idx] = max(abs(v));
        % 从H矩阵和测量创新中去除这个卫星的数据

        % 打印异常值信息
        fprintf('异常值检测：卫星 %d 在迭代 %d 被移除\n', satellite_numbers(outlier_idx), iteration);

        if rank(H_G) < 4
            fprintf('警告: 矩阵H_G在工作精度内为奇异的，无法计算最小二乘解。\n');
            break; % 跳出循环
        else
            H_G(outlier_idx, :) = [];
            measurement_innovation(outlier_idx) = [];
            num_satellites = num_satellites - 1;
            satellite_numbers(outlier_idx) = [];
            satellite_positions(:, outlier_idx) = [];
            pseudo_ranges(outlier_idx) = [];
            satellite_ranges(outlier_idx) = [];
            line_of_sight_vectors(:, outlier_idx) = [];
            %H_G = zeros(num_satellites, 4); % 测量矩阵
            %measurement_innovation = zeros(num_satellites, 1); % 测量创新
            %satellite_ranges(outlier_idx) = [];
            %line_of_sight_vectors(:, outlier_idx) = [];
            % 重新计算最小二乘解
            %delta_x = (H_G' * H_G) \ (H_G' * measurement_innovation);
            %estimated_position = estimated_position + delta_x(1:3);
            %predicted_receiver_clock_offset = predicted_receiver_clock_offset + delta_x(4);
        end
    end

    
    % 打印当前迭代的位置和时钟偏移
    fprintf('迭代 %d - 估计的位置和接收机时钟偏移：X=%f m, Y=%f m, Z=%f m, Clock Offset=%f m\n', iteration, estimated_position(1), estimated_position(2), estimated_position(3), predicted_receiver_clock_offset);

    % 计算与上一次迭代解的差异
    if iteration > 1 % 从第二次迭代开始计算差异
        solution_difference = norm(delta_x(1:3));
        fprintf('迭代 %d - 位置差异：%f m\n', iteration, solution_difference);
    end

end

% 转换为NED坐标系下的经纬度和高度
% 注意：需要有pv_ECEF_to_NED函数的实现
[L_b, lambda_b, h_b, v_eb_n] = pv_ECEF_to_NED(estimated_position, predicted_receiver_clock_offset);

% 显示最终结果
fprintf('最终结果：纬度 = %f°, 经度 = %f°, 高度 = %f米\n', L_b / deg_to_rad, lambda_b / deg_to_rad, h_b);