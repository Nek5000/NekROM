%% Three-Way Performance Comparison: S-Optimal Evolution
clear; clc;

% --- Simulation Parameters ---
% Increase N or M to see the 'final' version truly pull away
N_Bm = 2500;  % Total candidates
N    = 300;   % Total features
M    = 250;   % Total points to select (high M rewards qrupdate)

Vo = randn(N_Bm, N);
fprintf('Benchmarking: [%d x %d] Matrix, selecting %d points.\n\n', N_Bm, N, M);

%% 1. Original Approach
fprintf('Testing: s_opt_original... ');
tic;
idx_orig = s_opt_original(Vo, M, [], []);
t_orig = toc;
fprintf('Done (%.3fs)\n', t_orig);

%% 2. Improved (Vectorized) Approach
fprintf('Testing: s_opt_improved... ');
tic;
idx_impr = s_opt_improved(Vo, M, [], []);
t_impr = toc;
fprintf('Done (%.3fs)\n', t_impr);

%% 3. Final (QR Update) Approach
fprintf('Testing: s_opt_qr_fast... ');
tic;
idx_final = s_opt_qr_fast(Vo, M, [], []);
t_final = toc;
fprintf('Done (%.3fs)\n', t_final);

%% --- Result Visualization ---
labels = {'Original', 'Improved', 'Final (QR)'};
times = [t_orig, t_impr, t_final];

figure('Color', 'w');
b = bar(times, 'FaceColor', 'flat');
b.CData(3,:) = [0 0.5 0.8]; % Highlight the winner
set(gca, 'XTickLabel', labels);
ylabel('Execution Time (seconds)');
title('S-Optimal Algorithm Performance Comparison');
grid on;

fprintf('\nSummary of Results:\n');
fprintf('---------------------------------\n');
fprintf('Improved vs Original: %.2fx speedup\n', t_orig / t_impr);
fprintf('Final vs Original:    %.2fx speedup\n', t_orig / t_final);
fprintf('Final vs Improved:    %.2fx speedup\n', t_impr / t_final);
