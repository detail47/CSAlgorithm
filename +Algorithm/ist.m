function s = ist(g, Phi, lambda, max_iter)
    [~, N] = size(Phi);
    s = zeros(N, 1);
    g = g(:);

    g_norm = norm(g);
    L = compute_lipschitz_constant(Phi);
    t = 1 / L;  % 计算步长

    % 迭代更新
    residual = Phi * s - g;  % 初始残差
    for iter = 1:max_iter
        grad = Phi' * residual;    % 计算梯度
        z = s - t * grad;
        s = l1_proximal_operator(z, lambda * t); % 应用L1范数的近端算子
        residual = Phi * s - g;  % 更新残差
        if norm(residual) < 1e-2 * g_norm
            break;
        end
    end
end

function r = l1_proximal_operator(x, param)
    r = max(0, x - param) - max(0, -x - param);
end

function L = compute_lipschitz_constant(Phi)
    ATA = Phi' * Phi;
    L = max(eig(ATA));
end
