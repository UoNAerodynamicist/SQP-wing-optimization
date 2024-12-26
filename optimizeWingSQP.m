function optimizeWingSQP()
    % 设置变量的上下界
    lb = [2, 0.1, 0.1, -8];  % [b_min, c_root_min, c_tip_min, alpha0_min]
    ub = [4, 1, 1, 0];       % [b_max, c_root_max, c_tip_max, alpha0_max]

    % 初始解
    x0 = [3.5, 1, 0.8, -7];  % 初始猜测 [b, c_root, c_tip, alpha0]
    
    % 定义全局变量，用于记录迭代数据
    global iterData;
    iterData = struct('iterations', [], 'objective', [], 'gradNorm', [], ...
                     'constraintViolation', [], 'variableChange', [], ...
                     'gamma', []);

    % SQP优化选项
    options = optimoptions('fmincon', ...
                          'Algorithm', 'sqp', ...
                          'Display', 'iter', ...
                          'MaxIterations', 2000, ...
                          'OptimalityTolerance', 1e-6, ...
                          'StepTolerance', 1e-8, ...
                          'ConstraintTolerance', 1e-6, ...
                          'OutputFcn', @outfun);
    
    % 调用SQP优化器
    [x_opt, fval, exitFlag, output] = fmincon(@objective, x0, [], [], [], [], lb, ub, @constraints, options);
    
    % 计算最优解对应的gamma值
    S_opt = x_opt(1)  0.5  (x_opt(2) + x_opt(3));
    gamma_opt = findGamma(x_opt(1), x_opt(2), x_opt(3), x_opt(4), S_opt);
    
    % 输出最终结果
    fprintf('-------------------n');
    fprintf('Optimization Completedn');
    fprintf('-------------------n');
    fprintf('Optimal designn');
    fprintf('  b (Span) %.2fn', x_opt(1));
    fprintf('  c_root (Root chord) %.2fn', x_opt(2));
    fprintf('  c_tip (Tip chord) %.2fn', x_opt(3));
    fprintf('  alpha0 (Zero lift AoA) %.2fn', x_opt(4));
    fprintf('  Gamma %.2f degreesn', gamma_opt);
    [~, L2D_ratio] = sample_analysis(x_opt(1), x_opt(2), x_opt(3), x_opt(4), gamma_opt);
    fprintf('Optimal Lift-to-Drag Ratio %.2fn', L2D_ratio);
    fprintf('Total iterations %dn', output.iterations);

    % 根据 exitFlag 输出优化状态
    fprintf('-------------------n');
    fprintf('Exit Flag %dn', exitFlag);
    switch exitFlag
        case 1
            fprintf('Optimization terminated successfully.n');
        case 0
            fprintf('Maximum number of iterations reached.n');
        case -1
            fprintf('Optimization was terminated by the user.n');
        case -2
            fprintf('No feasible solution found.n');
        otherwise
            fprintf('Unknown termination condition.n');
    end

    % 绘制收敛图表
    plotConvergence();
end

% 目标函数
function f = objective(x)
    b_span = x(1);
    c_root = x(2);
    c_tip = x(3);
    alpha0 = x(4);
    
    % 计算机翼面积
    S = b_span  0.5  (c_root + c_tip);
    
    % 计算安装角 gamma
    gamma = findGamma(b_span, c_root, c_tip, alpha0, S);
    if isnan(gamma)
        f = 1e6; % 如果 gamma 不存在或无效，给出惩罚
        return;
    end

    % 调用 LLTanalysis 计算升阻比
    [~, L2D_ratio] = sample_analysis(b_span, c_root, c_tip, alpha0, gamma);
    
    % 权重系数
    w1 = 0.8;  % 升阻比的权重
    w2 = 0.2;  % gamma的权重
    
    % 归一化处理
    L2D_norm = L2D_ratio  100;  % 假设最大升阻比约为100
    gamma_norm = abs(gamma)  10; % 假设gamma的范围是±10度
    
    % 多目标函数：最大化LD，最小化gamma
    f = -(w1  L2D_norm - w2  gamma_norm);
end

% 非线性约束
function [c, ceq] = constraints(x)
    b_span = x(1);
    c_root = x(2);
    c_tip = x(3);
    alpha0 = x(4);
    
    % 计算机翼面积
    S = b_span  0.5  (c_root + c_tip);
    
    % 计算gamma
    gamma = findGamma(b_span, c_root, c_tip, alpha0, S);

    % 非线性不等式约束
    c = [];
    c(1) = 0.1 - S;            % S = 0.1
    c(2) = S - 4;              % S = 4
    c(3) = 2 - (b_span^2  S); % b^2  S  2
    c(4) = c_tip - c_root;     % c_tip = c_root
    c(5) = abs(gamma) - 10;    % gamma = 10度

    % 等式约束 (升力等于重力)
    if isnan(gamma)  gamma  -10  gamma  10
        ceq = 1e6; % 无效时强制违反约束
    else
        W = 100  9.8; % 重力 (N)
        rho = 1.225;   % 空气密度 (kgm^3)
        V = 20;        % 巡航速度 (ms)

        % 计算升力系数 C_l
        [C_l, ~] = sample_analysis(b_span, c_root, c_tip, alpha0, gamma);

        % 计算升力
        L = 0.5  rho  C_l  S  V^2;

        % 等式约束 升力等于重力
        ceq = L - W;
    end
end

% 输出函数用于记录每次迭代的数据
function stop = outfun(x, optimValues, state)
    global iterData;

    switch state
        case 'iter'
            % 提取设计变量
            b_span = x(1);
            c_root = x(2);
            c_tip = x(3);
            alpha0 = x(4);
            S = b_span  0.5  (c_root + c_tip);
            
            % 计算当前迭代的 gamma
            gamma = findGamma(b_span, c_root, c_tip, alpha0, S);

            % 记录迭代数据
            iterData.iterations = [iterData.iterations; optimValues.iteration];
            iterData.objective = [iterData.objective; optimValues.fval];
            iterData.gradNorm = [iterData.gradNorm; norm(optimValues.gradient)];
            iterData.variableChange = [iterData.variableChange; optimValues.stepsize];
            iterData.constraintViolation = [iterData.constraintViolation; optimValues.constrviolation];
            iterData.gamma = [iterData.gamma; gamma];
    end
    stop = false;
end

% 绘制收敛图表
function plotConvergence()
    global iterData;

    % 创建2x3网格布局的图表
    figure('Position', [100, 100, 1200, 800]);
    
    % 图1：目标函数值 vs 迭代次数
    subplot(2, 3, 1);
    plot(iterData.iterations, iterData.objective, '-o', 'LineWidth', 1.5);
    xlabel('Iteration');
    ylabel('Objective Function Value');
    title('Convergence of Objective Function');
    grid on;

    % 图2：梯度范数 vs 迭代次数
    subplot(2, 3, 2);
    plot(iterData.iterations, iterData.gradNorm, '-o', 'LineWidth', 1.5);
    xlabel('Iteration');
    ylabel('Gradient Norm');
    title('Convergence of Gradient Norm');
    grid on;

    % 图3：约束违背量 vs 迭代次数
    subplot(2, 3, 3);
    plot(iterData.iterations, iterData.constraintViolation, '-o', 'LineWidth', 1.5);
    xlabel('Iteration');
    ylabel('Constraint Violation');
    title('Convergence of Constraint Violation');
    grid on;

    % 图4：变量变化量 vs 迭代次数
    subplot(2, 3, 4);
    plot(iterData.iterations, iterData.variableChange, '-o', 'LineWidth', 1.5);
    xlabel('Iteration');
    ylabel('Change in Variables');
    title('Convergence of Variables');
    grid on;

    % 图5：gamma 值 vs 迭代次数
    subplot(2, 3, 5);
    plot(iterData.iterations, iterData.gamma, '-o', 'LineWidth', 1.5);
    xlabel('Iteration');
    ylabel('Gamma Value (degrees)');
    title('Convergence of Gamma');
    grid on;

    % 调整子图间距
    set(gcf, 'DefaultAxesFontSize', 10);
    sgtitle('Optimization Convergence Analysis', 'FontSize', 14);
end