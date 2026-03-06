ALGORITHMS = {'IST', 'OMP', 'IHT', 'PGD L_{1/2}'};

% 生成稀疏信号
K = 2;
N = 200;
M = 16;
max_iter = 1e4;
lambda = 1;
SNR = 3;

s_true = complex(zeros(N, 1)); % 初始化为复数零向量
support = randperm(N, K);
s_true(support) = 1 * exp(1j * 2 * pi * rand(K, 1)) + 0.1 * (randn(K, 1) + 1j * randn(K, 1)); % 稀疏信号，包含随机相位和小幅噪声

% 生成测量矩阵
Phi = (randn(M, N) + 1j * randn(M, N)) / sqrt(2); % 复数高斯随机矩阵，归一化

% 生成测量值
g_raw = Phi * s_true;
g = awgn(g_raw, SNR, 'measured'); % 添加高斯白噪声

% 评估算法性能
for i = 1:length(ALGORITHMS)
    alg_name = ALGORITHMS{i};
    params = struct();
    params.Phi = Phi;
    params.g = g;
    params.K = K;
    params.max_iter = max_iter;
    params.lambda = lambda;
    [s_est, info] = Algorithm.run_algorithm(alg_name, params);

    % 计算恢复误差
    SNR_res = Shared.compute_SNR(s_true, s_est);
    fprintf('%s: SNR = %.2f dB\n', alg_name, SNR_res);

    % 绘制恢复结果
    figure;
    stem(1:N, abs(s_true), 'r', 'filled', 'DisplayName', 'True Signal');
    hold on;
    stem(1:N, abs(s_est), 'b', 'filled', 'DisplayName', sprintf('%s Estimate', alg_name));
    title(sprintf('%s Recovery', alg_name));
    xlabel('Index');
    ylabel('Magnitude');
    legend;
    grid on;
end
