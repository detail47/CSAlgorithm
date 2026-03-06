function [s_est, info] = iht(g, Phi, max_iter, K)
    [~, N] = size(Phi);
    s_est = zeros(N, 1);
    g = g(:);

    g_norm = norm(g);
    L = compute_lipschitz_constant(Phi);
    t = 1 / L;  % 计算步长
    
    residual_tol = 1e-2;
    plateau_tol = 1e-6;
    plateau_patience = 20;
    plateau_count = 0;

    residual = Phi * s_est - g;  % 初始残差
    residual_norm = norm(residual);
    stop_reason = "max_iter";
    iter_done = 0;

    for iter = 1:max_iter
        z = s_est - Phi' * residual * t;    % 计算梯度步长更新
        s_est = hard_thresholding(z, K); % 应用硬阈值处理

        prev_residual_norm = residual_norm;
        residual = Phi * s_est - g;  % 更新残差
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

function r = hard_thresholding(x, K)
    % 对向量x进行硬阈值处理，保留最大的K个元素
    r = zeros(size(x));
    [~, idx] = sort(abs(x), 'descend');
    r(idx(1:K)) = x(idx(1:K));
end

function L = compute_lipschitz_constant(Phi)
    ATA = Phi' * Phi;
    L = max(eig(ATA));
end
