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

% Time Incremental Version
tic;
p_incr = gpode_incremental(U, M);
t_incr = toc;
fprintf('Incremental Version: %.4f seconds\n', t_incr);

% Time Improved Version
tic;
p_imp = gpode_improved(U, M);
t_imp = toc;
fprintf('Original Version:    %.4f seconds\n', t_imp);



fprintf('Incremental Speedup: %.2fx\n', t_orig / t_incr);
fprintf('Improved Speedup: %.2fx\n', t_imp / t_incr);

% Verification: Check if the selected indices match
if isequal(p_orig, p_incr)
    fprintf('Results: Identical ✅\n');
elseif isequal(p_orig, p_imp)
    fprintf('Results: Identical ✅\n');
else
    fprintf('Results: Different (Numerical noise or logic divergence) ⚠️\n');
end
