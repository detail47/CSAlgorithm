function [s, info] = pgd_L1_2_after_L1(g, Phi, lambda_ist, lambda_pgd, max_iter)
    % 先使用L1范数进行稀疏化，然后再使用L1/2范数进行进一步优化
    % Step 1: 使用L1范数进行稀疏化
    [s0, info_ist] = Algorithm.ist(g, Phi, lambda_ist, max_iter);
    % Step 2: 使用L1/2范数进行进一步优化
    [s, info_pgdL1_2] = Algorithm.pgd_L1_2(g, Phi, lambda_pgd, max_iter, s0);
    % 合并信息
    info.ist = info_ist;
    info.pgd = info_pgdL1_2;
end
