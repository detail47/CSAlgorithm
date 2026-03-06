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
            [s_est, info] = Algorithm.pgd_1f2(g, Phi, lambda, max_iter);
        otherwise
            error('未知算法: %s', alg_name);
    end
end