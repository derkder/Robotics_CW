% 定义 Sagnac 效应补偿矩阵的函数
function C_e = sagnac_matrix(r_ea_e, omega_ie, c)
    % 根据方程（2）计算 Sagnac 效应补偿矩阵
    C_e = [1, omega_ie * r_ea_e(2) / c, -omega_ie * r_ea_e(1) / c;
          -omega_ie * r_ea_e(2) / c, 1, 0;
           omega_ie * r_ea_e(1) / c, 0, 1];
end

% 计算近似范围的函数
function ranges = predict_ranges(r_eb_e, satellite_positions, omega_ie, c)
    num_satellites = size(satellite_positions, 2);
    ranges = zeros(1, num_satellites);
    for j = 1:num_satellites
        % 计算用户位置到卫星的矢量
        r_ea_e = satellite_positions(:, j) - r_eb_e;
        % 首先使用单位矩阵计算初始距离估计
        range = sqrt(r_ea_e' * r_ea_e);
        % 计算 Sagnac 效应补偿矩阵
        C_e = sagnac_matrix(r_ea_e, omega_ie, c);
        % 修正矢量
        r_ea_e_corrected = C_e * r_ea_e;
        % 重新计算修正后的距离
        range_corrected = sqrt(r_ea_e_corrected' * r_ea_e_corrected);
        ranges(j) = range_corrected;
    end
end