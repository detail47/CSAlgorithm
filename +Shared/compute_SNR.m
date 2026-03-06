function SNR = compute_SNR(s_true, s_est)
    error = norm(s_est - s_true);
    sum_true = norm(s_true);
    SNR = 20 * log10(sum_true / max(error, eps));
end