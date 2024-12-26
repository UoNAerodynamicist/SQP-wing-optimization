function gamma = findgamma(b_span, c_root, c_tip, alpha0, S)
    % findgamma: 使用二分法计算满足升力条件的安装角 gamma
    % 常数定义
    W = 100 * 9.8; % 重力 (N)
    rho = 1.225;   % 空气密度 (kg/m^3)
    V = 20;        % 巡航速度 (m/s)

    % 二分法搜索范围
    gamma_min = -20; % 最小角度
    gamma_max = 20;  % 最大角度
    tol = 1e-7;      % 收敛精度
    max_iter = 1000; % 最大迭代次数
    
    % 计算初始升力误差
    error_min = calculateLiftError(b_span, c_root, c_tip, alpha0, gamma_min, S, W, rho, V);
    error_max = calculateLiftError(b_span, c_root, c_tip, alpha0, gamma_max, S, W, rho, V);

    % 验证初始区间是否包含解
    if error_min * error_max > 0
        warning('Unable to find a valid gamma: No root in the interval [%f, %f]', gamma_min, gamma_max);
        gamma = NaN;
        return;
    end

    % 开始二分法迭代
    iter = 0;
    while (gamma_max - gamma_min) > tol && iter < max_iter
        % 中间点
        gamma_mid = (gamma_min + gamma_max) / 2;

        % 计算中间点误差
        error_mid = calculateLiftError(b_span, c_root, c_tip, alpha0, gamma_mid, S, W, rho, V);

        % 更新区间
        if error_min * error_mid < 0
            gamma_max = gamma_mid;
            error_max = error_mid;
        else
            gamma_min = gamma_mid;
            error_min = error_mid;
        end

        iter = iter + 1;
    end

    % 计算最终解
    gamma = (gamma_min + gamma_max) / 2;

    % 检查是否达到最大迭代次数
    if iter == max_iter
        warning('Max iterations reached in findgamma.');
    end
end

function error = calculateLiftError(b_span, c_root, c_tip, alpha0, gamma, S, W, rho, V)
    % 计算升力与重力的差值

    % 调用升力线理论计算升力系数
    [C_l, ~] = sample_analysis(b_span, c_root, c_tip, alpha0, gamma);

    % 计算升力
    L = 0.5 * rho * C_l * S * V^2;

    % 计算误差
    error = L - W;
end