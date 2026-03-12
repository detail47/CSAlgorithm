% 生成稀疏信号
K_MAX = 3;
N = 200;
M = 15;
SNR = 10;
EPOCH = 1e4;
VIEWS = 10;
MAX_ITER = 1e4;
LAMBDA = 1;
LAMBDA_L12_IST = 5;
LAMBDA_L12_PGD = 0.5;

ALGORITHMS = {'IST', 'PGD L_{1/2}', 'PGD L_{1/2} after IST'};
SNR_Matrix = zeros(EPOCH, length(ALGORITHMS));
TOTAL = EPOCH * length(ALGORITHMS);
elapsed_time = zeros(length(ALGORITHMS), 1);
disp('Starting benchmark...');
for alg_idx = 1:length(ALGORITHMS)
    alg_name = ALGORITHMS{alg_idx};
    fprintf('Running %s...\n', alg_name);
    tich = tic;
    for epoch = 1:EPOCH
        s_true = complex(zeros(N, 1));
        K = randi([1, K_MAX]);
        non_zero_indices = randperm(N, K);
        s_true(non_zero_indices) = 1 * exp(1j * 2 * pi * rand(K, 1)) + 0.1 * (randn(K, 1) + 1j * randn(K, 1));
        Phi = (randn(M, N) + 1j * randn(M, N)) / sqrt(2);
        g_raw = Phi * s_true;
        g = awgn(g_raw, SNR, 'measured');
        params = struct();
        params.K = K;
        params.lambda = LAMBDA;
        params.max_iter = MAX_ITER;
        params.Phi = Phi;
        params.g = g;
        params.lambda_ist = LAMBDA_L12_IST;
        params.lambda_pgd = LAMBDA_L12_PGD;
        [s_est, ~] = Algorithm.run_algorithm(alg_name, params);
        SNR_Matrix(epoch, alg_idx) = Shared.compute_SNR(s_true, s_est);
    end
    elapsed_time(alg_idx) = toc(tich);
end
fid = fopen('Log/Benchmark_Results.log', 'w');
for alg_idx = 1:length(ALGORITHMS)
    alg_name = ALGORITHMS{alg_idx};
    single_time = elapsed_time(alg_idx) / EPOCH;
    str = sprintf('%s: Total Time = %.2f s, Average Time per Run = %.4f ms\n', ...
        alg_name, elapsed_time(alg_idx), single_time * 1e3);
    fprintf(fid, '%s', str);
    fprintf('%s', str);
end
fclose(fid);

% 计算统计指标
mean_SNR = zeros(length(ALGORITHMS), 1);
median_SNR = zeros(length(ALGORITHMS), 1);
std_SNR = zeros(length(ALGORITHMS), 1);
for alg_idx = 1:length(ALGORITHMS)
    alg_name = ALGORITHMS{alg_idx};
    SNR_values = SNR_Matrix(:, alg_idx);
    mean_SNR(alg_idx) = mean(SNR_values);
    median_SNR(alg_idx) = median(SNR_values);
    std_SNR(alg_idx) = std(SNR_values);
end
fid = fopen('Log/Benchmark_Statistics.log', 'w');
for alg_idx = 1:length(ALGORITHMS)
    alg_name = ALGORITHMS{alg_idx};
    str = sprintf( '%s: Mean SNR = %.2f dB, Median SNR = %.2f dB, Std SNR = %.2f dB\n', ...
        alg_name, mean_SNR(alg_idx), median_SNR(alg_idx), std_SNR(alg_idx));
    fprintf(fid, '%s', str);
    fprintf('%s', str);
end
fclose(fid);

% 绘制恢复样例
disp('Generating recovery examples...');
for v = 1:VIEWS
    s_true = complex(zeros(N, 1));
    K = randi([1, K_MAX]);
    non_zero_indices = randperm(N, K);
    s_true(non_zero_indices) = 1 * exp(1j * 2 * pi * rand(K, 1)) + 0.1 * (randn(K, 1) + 1j * randn(K, 1));
    Phi = (randn(M, N) + 1j * randn(M, N)) / sqrt(2);
    g_raw = Phi * s_true;
    g = awgn(g_raw, SNR, 'measured');
    s_est = zeros(N, length(ALGORITHMS));
    info = cell(1, length(ALGORITHMS));
    for alg_idx = 1:length(ALGORITHMS)
        alg_name = ALGORITHMS{alg_idx};
        params = struct();
        params.K = K;
        params.lambda = LAMBDA;
        params.max_iter = MAX_ITER;
        params.Phi = Phi;
        params.g = g;
        params.lambda_ist = LAMBDA_L12_IST;
        params.lambda_pgd = LAMBDA_L12_PGD;
        [s_est(:, alg_idx), info{alg_idx}] = Algorithm.run_algorithm(alg_name, params);
    end
    % 绘制恢复结果
    figure("Visible", "off");
    rows = ceil((length(ALGORITHMS) + 1) / 2 + 1);
    subplot(rows, 2, 1);
    stem(abs(g_raw), 'filled');
    title('原始测量值幅度');
    subplot(rows, 2, 2);
    stem(abs(g), 'filled');
    title('带噪声的测量值幅度');
    subplot(rows, 2, 3);
    stem(abs(s_true), 'filled');
    title('真实信号幅度');
    for alg_idx = 1:length(ALGORITHMS)
        subplot(rows, 2, alg_idx + 3);
        stem(abs(s_est(:, alg_idx)), 'filled');
        title(sprintf('%s恢复的信号幅度', ALGORITHMS{alg_idx}));
    end
    saveas(gcf, sprintf('Image/Recovery_Example_%d.png', v));
    % 保存迭代信息
    info_text = cell(1, length(ALGORITHMS));
    for alg_idx = 1:length(ALGORITHMS)
        alg_name = ALGORITHMS{alg_idx};
        info_struct = info{alg_idx};
        str_struct = jsonencode(info_struct, 'PrettyPrint', true);
        info_text{alg_idx} = sprintf('Algorithm: %s\n%s', alg_name, str_struct);
    end
    fid = fopen(sprintf('Log/Recovery_Example_%d.log', v), 'w');
    fprintf(fid, '%s\n\n', info_text{:});
    fclose(fid);
end

% 绘制误差率分布图
NUM_BAR = 50;
figure;
hold on;
edges = linspace(min(SNR_Matrix(:)), max(SNR_Matrix(:)), NUM_BAR+1);
disp('Plotting SNR distribution...');
for alg_idx = 1:length(ALGORITHMS)
    [counts, binEdges] = histcounts(SNR_Matrix(:, alg_idx), edges, 'Normalization', 'pdf');
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
    plot(binCenters, counts, 'LineWidth', 1.8);
end
title('SNR分布比较');
xlabel('SNR (dB)');
ylabel('概率密度');
legend(ALGORITHMS);
hold off;
saveas(gcf, 'Image/SNR_Distribution.png');
disp('Benchmark completed. Results saved to Image/Recovery_Example_*.png and Image/SNR_Distribution.png');
