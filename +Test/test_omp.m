% 生成稀疏信号
K = 5;
N = 100;
s_true = zeros(N, 1);
support = randperm(N, K);
s_true(support) = randn(K, 1);

% 生成测量矩阵
M = 50;
Phi = randn(M, N);

% 生成测量值
g = Phi * s_true;

% 使用OMP算法恢复信号
s_est = Algorithm.omp(g, Phi, K);

% 计算恢复误差
error = norm(s_est - s_true);
sum = norm(s_true);
error_rate = error / sum;
disp(['OMP恢复误差: ' num2str(error)]);
disp(['OMP恢复误差率: ' num2str(error_rate)]);

% 显示结果
figure;
subplot(2, 1, 1);
stem(s_true, 'filled');
title('真实信号');
subplot(2, 1, 2);
stem(s_est, 'filled');
title('OMP恢复的信号');

