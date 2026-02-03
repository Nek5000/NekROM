% Benchmark: GPODE Original vs. Incremental
rng(42); % For reproducibility
N = 5000; % Number of potential locations (rows)
K = 50;   % Number of modes (columns)
M = 200;  % Total points to select

U = randn(N, K);

fprintf('Starting Benchmark...\n');

% Time Original Version
tic;
p_orig = gpode_original(U, M);
t_orig = toc;
fprintf('Original Version:    %.4f seconds\n', t_orig);

% Time Adaptive Version
tic;
p_incr = qdeim_adaptive(U, M);
t_incr = toc;
fprintf('Adaptive Version: %.4f seconds\n', t_incr);


% Time Improved Version
tic;
p_imp = gpode_improved(U, M);
t_imp = toc;
fprintf('Improved Version:    %.4f seconds\n', t_imp);


fprintf('Adaptive Speedup: %.2fx\n', t_orig / t_incr);
fprintf('Improved Speedup: %.2fx\n', t_orig / t_imp);

% Verification: Check if the selected indices match
if isequal(p_orig, p_incr)
    fprintf('Results: Adaptive Identical ✅\n');
else
    fprintf('Results: Different (Numerical noise or logic divergence) ⚠️\n');
end
if isequal(p_orig, p_imp)
    fprintf('Results: Improved Identical ✅\n');
else
    fprintf('Results: Different (Numerical noise or logic divergence) ⚠️\n');
end
