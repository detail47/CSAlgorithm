function s = omp(g, Phi, K)
    [M, N] = size(Phi);
    s = zeros(N, 1);
    g = g(:);

    % 初始化
    residual = g;
    coef = zeros(K, 1);
    support = zeros(K, 1);
    selected_atoms = zeros(M, K);

    for iter = 1:K
        corr = abs(Phi' * residual);     % 计算相关性
        corr(support(1:iter-1)) = 0;    % 排除已选

        [~, idx] = max(corr);           % 选择最大相关性
        support(iter) = idx;
        selected_atoms(:, iter) = Phi(:, idx); % 选取对应列

        % coef = selected_atoms(:, 1:iter) \ g; % 最小二乘求解系数
        coef = pinv(selected_atoms(:, 1:iter)) * g; % 最小二乘求解系数，使用广义逆避免奇异
        residual = g - selected_atoms(:, 1:iter) * coef; % 更新
        if norm(residual) < 1e-10
            break;
        end
    end
    s(support(1:iter)) = coef;
end
