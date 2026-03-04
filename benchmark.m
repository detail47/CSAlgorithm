% 生成稀疏信号
K = 3;
N = 80;
M = 16;

EPOCH = 1e5;
residual_normal = zeros(EPOCH, 2);
for i = 1:EPOCH
    % 生成稀疏信号
    s_true = zeros(N, 1);
    support = randperm(N, K);
    s_true(support) = randn(K, 1);
    Phi = randn(M, N);

    % 生成测量值
    g = Phi * s_true;
    g = g + 0.1 * randn(M, 1); % 添加少量噪声

    % 使用IST算法恢复信号
    lambda = 1; % 正则化参数
    max_iter = 1e4; % 最大迭代次数
    s_est = Algorithm.ist(g, Phi, lambda, max_iter);
    % 计算IST恢复误差
    error = norm(s_est - s_true);
    sum = norm(s_true);
    error_rate = error / sum;
    residual_normal(i, 1) = error_rate;

    % 使用OMP算法恢复信号
    s_est_omp = Algorithm.omp(g, Phi, K);
    % 计算OMP恢复误差
    error_omp = norm(s_est_omp - s_true);
    error_rate_omp = error_omp / sum;
    residual_normal(i, 2) = error_rate_omp;
    
    % 显示当前迭代的误差率
    if mod(i, 100) ~= 0
        continue;
    end
    fprintf('Epoch %d: IST Residual = %.4f, OMP Residual = %.4f\n', i, error_rate, error_rate_omp);
end

% 绘制误差率分布图
figure;
hold on;
edges = linspace(min(residual_normal(:)), max(residual_normal(:)), 501);
histogram(residual_normal(:, 1), edges, 'DisplayName', 'IST', 'FaceAlpha', 0.6);
histogram(residual_normal(:, 2), edges, 'DisplayName', 'OMP', 'FaceAlpha', 0.6);
title('Residual Distribution');
xlabel('Residual');
ylabel('Frequency');
legend;
hold off;

% 计算统计指标
mean_ist = mean(residual_normal(:, 1));
mean_omp = mean(residual_normal(:, 2));
std_ist = std(residual_normal(:, 1));
std_omp = std(residual_normal(:, 2));
fprintf('IST Mean Residual: %.4f, Std Dev: %.4f\n', mean_ist, std_ist);
fprintf('OMP Mean Residual: %.4f, Std Dev: %.4f\n', mean_omp, std_omp);

