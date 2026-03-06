function [s, info] = omp(g, Phi, K)
    [M, N] = size(Phi);
    s = zeros(N, 1);
    g = g(:);

    % 初始化
    g_norm = norm(g);
    residual = g;
    residual_norm = norm(residual);
    coef = zeros(K, 1);
    support = zeros(K, 1);
    selected_atoms = zeros(M, K);
    residual_tol = 1e-2;
    stop_reason = "max_iter";
    iter_done = 0;

    for iter = 1:K
        corr = abs(Phi' * residual);     % 计算相关性
        corr(support(1:iter-1)) = 0;    % 排除已选

        [~, idx] = max(corr);           % 选择最大相关性
        support(iter) = idx;
        selected_atoms(:, iter) = Phi(:, idx); % 选取对应列

        % coef = selected_atoms(:, 1:iter) \ g; % 最小二乘求解系数
        coef = pinv(selected_atoms(:, 1:iter)) * g; % 最小二乘求解系数，使用广义逆避免奇异
        residual = g - selected_atoms(:, 1:iter) * coef; % 更新
        residual_norm = norm(residual);
        iter_done = iter;
        if residual_norm < residual_tol * g_norm
            stop_reason = "residual_tol";
            break;
        end
    end
    s(support(1:iter)) = coef;

    if nargout > 1
        info = struct();
        info.iter = iter_done;
        info.residual_norm = residual_norm;
        info.relative_residual = residual_norm / max(g_norm, eps);
        info.stop_reason = stop_reason;
        info.support = support(1:iter);
    end
end
