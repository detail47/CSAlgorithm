% 所有算法的统一接口
% IST: 迭代软阈值算法
% OMP: 正交匹配追踪算法
% IHT: 迭代硬阈值算法
% PGD L_{1/2}: 基于L_{1/2}范数的近端梯度下降算法
% PGD L_{1/2} after IST: 先使用IST算法进行预处理，然后再使用PGD L_{1/2}算法进行精细优化

function [s_est, info] = run_algorithm(alg_name, param)
    g = param.g;
    Phi = param.Phi;
    switch alg_name
        case 'IST'
            lambda = param.lambda;
            max_iter = param.max_iter;
            [s_est, info] = Algorithm.ist(g, Phi, lambda, max_iter);
        case 'OMP'
            K = param.K;
            [s_est, info] = Algorithm.omp(g, Phi, K);
        case 'IHT'
            K = param.K;
            max_iter = param.max_iter;
            [s_est, info] = Algorithm.iht(g, Phi, max_iter, K);
        case 'PGD L_{1/2}'
            lambda = param.lambda;
            max_iter = param.max_iter;
            [s_est, info] = Algorithm.pgd_L1_2(g, Phi, lambda, max_iter);
        case 'PGD L_{1/2} after IST'
            if ~isfield(param, 'lambda_ist')
                lambda_ist = param.lambda;
            else
                lambda_ist = param.lambda_ist;
            end
            if ~isfield(param, 'lambda_pgd')
                lambda_pgd = param.lambda;
            else
                lambda_pgd = param.lambda_pgd;
            end
            max_iter = param.max_iter;
            [s_est, info] = Algorithm.pgd_L1_2_after_L1(g, Phi, lambda_ist, lambda_pgd, max_iter);
        otherwise
            error('未知算法: %s', alg_name);
    end
end