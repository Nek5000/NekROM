%% Gappy POD Three-Way Performance & Validation Benchmark
clear; clc;

% --- Parameters ---
N_rows = 5000;      % Increased for better timing resolution
p_basis = 60;       % Number of POD modes
n_samples = 180;    % Oversampling (n > p)

% Generate a basis with decaying singular values (more realistic)
rng(42);
[Q, ~] = qr(randn(N_rows, p_basis), 0);
S = diag(logspace(0, -6, p_basis)); % Add spectral decay
U_full = Q * S;

%% --- Run Algorithms ---

% 1. Original Algorithm
fprintf('Running Original...\n');
tic;
idx_orig = gappy_pod_original(U_full, n_samples);
t_orig = toc;

% 2. QR-Optimized Algorithm
fprintf('Running QR-Optimized...\n');
tic;
idx_opt = gappy_pod_qr(U_full, n_samples);
t_opt = toc;

% 3. New Algorithm (The "New" candidate)
fprintf('Running Gappy POD New...\n');
tic;
idx_new = gappy_pod_new(U_full, n_samples);
t_new = toc;

%% --- Validation Logic ---
% We use the last mode (p) as a test case for reconstruction using (1:p-1)
p = p_basis;
U_train = U_full(:, 1:p-1);
u_test  = U_full(:, p);

% Helper to compute metrics
compute_metrics = @(indices) deal(...
    norm(u_test - U_train * (U_train(indices, :) \ u_test(indices))), ...
    cond(U_train(indices, :)));

[res_orig, cond_orig] = compute_metrics(idx_orig);
[res_opt,  cond_opt]  = compute_metrics(idx_opt);
[res_new,  cond_new]  = compute_metrics(idx_new);

%% --- Report Results ---
fprintf('\n--- Performance Benchmark ---\n');
data = { 'Original', t_orig, res_orig, cond_orig; ...
         'QR-Opt',   t_opt,  res_opt,  cond_opt; ...
         'New',      t_new,  res_new,  cond_new };

fprintf('%-12s | %-10s | %-12s | %-10s\n', 'Method', 'Time (s)', 'Residual', 'Condition #');
fprintf('-------------------------------------------------------------\n');
for i = 1:size(data, 1)
    fprintf('%-12s | %-10.4f | %-12.2e | %-10.2e\n', data{i, :});
end

fprintf('\n--- Comparison Analysis ---\n');
fprintf('Speedup (Opt vs Orig): %.2fx\n', t_orig / t_opt);
fprintf('Speedup (New vs Orig): %.2fx\n', t_orig / t_new);

% Logic Check
[~, best_idx] = min([res_orig, res_opt, res_new]);
methods = {'Original', 'QR-Opt', 'New'};
fprintf('Best Reconstruction: %s\n', methods{best_idx});
