% 生成稀疏信号
K = 2;
N = 100;
s_true = zeros(N, 1);
support = randperm(N, K);
s_true(support) = randn(K, 1);

% 生成测量矩阵
M = 16;
Phi = randn(M, N);

% 生成测量值
g = Phi * s_true;
g = g + 0.1 * randn(M, 1); % 添加少量噪声

% 使用IST算法恢复信号
lambda = 1; % 正则化参数
max_iter = 1e4; % 最大迭代次数
s_est = Algorithm.ist(g, Phi, lambda, max_iter);

% 计算恢复误差
error = norm(s_est - s_true);
sum = norm(s_true);
error_rate = error / sum;
disp(['IST恢复误差率: ' num2str(error_rate)]);

% 显示结果
figure;
subplot(2, 1, 1);
stem(s_true, 'filled');
title('真实信号');
subplot(2, 1, 2);
stem(s_est, 'filled');
title('IST恢复的信号');

