function [s, info] = pgd_1f2(g, Phi, lambda, max_iter)
    [~, N] = size(Phi);
    s = zeros(N, 1);
    g = g(:);

    g_norm = norm(g);
    L = compute_lipschitz_constant(Phi);
    t = 1 / L;  % 计算步长

    % 停止参数：
    % 1) residual_tol: 达到目标残差比例直接停止
    % 2) plateau_tol + plateau_patience: 连续多次改进极小，判定进入地板
    residual_tol = 1e-2;
    plateau_tol = 1e-6;
    plateau_patience = 20;
    plateau_count = 0;

    % 迭代更新
    residual = Phi * s - g;  % 初始残差
    residual_norm = norm(residual);
    stop_reason = "max_iter";
    iter_done = 0;

    for iter = 1:max_iter
        grad = Phi' * residual;    % 计算梯度
        z = s - t * grad;
        s = proximal_operator(z, lambda * t); % 应用$L_\frac12$范数的近端算子

        prev_residual_norm = residual_norm;
        residual = Phi * s - g;  % 更新残差
        residual_norm = norm(residual);
        iter_done = iter;

        if residual_norm < residual_tol * g_norm
            stop_reason = "residual_tol";
            break;
        end

        improvement = (prev_residual_norm - residual_norm) / max(prev_residual_norm, eps);
        if improvement < plateau_tol
            plateau_count = plateau_count + 1;
        else
            plateau_count = 0;
        end

        if plateau_count >= plateau_patience
            stop_reason = "plateau";
            break;
        end
    end

    if nargout > 1
        info = struct();
        info.iter = iter_done;
        info.residual_norm = residual_norm;
        info.relative_residual = residual_norm / max(g_norm, eps);
        info.stop_reason = stop_reason;
        info.step_size = t;
        info.L = L;
        info.plateau_count = plateau_count;
    end
end

% $L_\frac12$ 范数的近端算子
function r = proximal_operator(x, lambda)
    r = zeros(size(x));
    for i = 1:length(x)
        if abs(x(i)) > (3/2) * lambda^(2/3)
            r(i) = x(i) - (lambda * sign(x(i))) / (abs(x(i))^(1/3));
        else
            r(i) = 0;
        end
    end
end

function L = compute_lipschitz_constant(Phi)
    ATA = Phi' * Phi;
    L = max(eig(ATA));
end
