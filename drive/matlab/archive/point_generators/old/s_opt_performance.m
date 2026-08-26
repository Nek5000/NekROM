%% S-Opt Performance Benchmark
clear; clc;

% 1. Setup Parameters
N_Bm = 2000;  % Number of candidate points
N    = 50;    % Number of basis functions
M    = 40;    % Number of points to select
Vo   = randn(N_Bm, N); 

fprintf('Benchmarking Matrix: %d x %d\n', N_Bm, N);

% 2. Profile Original (ensure original code is saved as s_opt_original.m)
if exist('s_opt_original', 'file')
    tic;
    idx_orig = s_opt_original(Vo, M, [], '');
    t_orig = toc;
    fprintf('Original:  %.4f seconds\n', t_orig);
else
    t_orig = NaN;
    fprintf('s_opt_original.m not found. Skipping original profile.\n');
end

% 3. Profile Optimized
tic;
idx_opt = s_opt(Vo, M, [], '');
t_opt = toc;
fprintf('Optimized: %.4f seconds\n', t_opt);

% 4. Compare
if ~isnan(t_orig)
    fprintf('Speedup:   %.2fx\n', t_orig / t_opt);
    if isequal(idx_orig, idx_opt)
        fprintf('Result:    SUCCESS (Indices Match)\n');
    else
        fprintf('Result:    FAILURE (Indices Differ)\n');
    end
end
