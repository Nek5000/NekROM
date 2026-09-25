%% GNAT Performance Comparison Benchmark
clear; clc;

% --- Parameters ---
n = 5000;          % Number of degrees of freedom (spatial points)
m_total = 100;     % Total available basis vectors
m_used = 50;       % Basis vectors to use for sampling
nsr = 150;         % Number of samples required (usually > m_used)

% Generate a synthetic orthonormal basis Q
[Q, ~] = qr(rand(n, m_total), 0);

fprintf('Benchmarking GNAT implementations with n=%d, m=%d, ns=%d...\n\n', n, m_used, nsr);

% --- Test Original Version ---
tic;
phi_orig = gnat_original(Q, m_used, nsr);
t_orig = toc;
fprintf('Original Version Time:  %.4f seconds\n', t_orig);

% --- Test Optimized Version ---
tic;
phi_opt = gnat_optimized(Q, m_used, nsr);
t_opt = toc;
fprintf('Optimized Version Time: %.4f seconds\n', t_opt);

% --- Results Analysis ---
speedup = t_orig / t_opt;
fprintf('Speedup Factor:         %.2fx\n', speedup);

% Check if selected indices are identical
% Note: sort them because the original code returned a sort order 
% while optimized returns the actual indices.
if isequal(sort(phi_orig), sort(phi_opt))
    fprintf('Validation:             SUCCESS (Both versions selected the same indices)\n');
else
    fprintf('Validation:             MISMATCH (Check logic differences)\n');
end

%% --- Original Implementation (Included for the test) ---
function [phi] = gnat_original(Q, m_used, nsr)
    phi = []; used = []; U = Q(:,1); Q_sampled = [];
    m = min(m_used, size(Q,2));
    ns = iif(nsr > 0, nsr, m);
    n = size(Q,1);
    ns_mod_nr = mod(ns,m);
    
    % First iteration
    nsi = iif(mod(ns,m) > 0, floor(ns/m)+1, floor(ns/m));
    P = [];
    for i = 1:nsi
        s_row = -1; s_row_val = -Inf;
        for j = 1:n
            if ~any(used == j)
                if s_row == -1 || s_row_val < abs(Q(j,1))
                    s_row = j; s_row_val = abs(Q(j,1));
                end
            end
        end
        used(end+1) = s_row; phi(end+1) = s_row;
        newPcol = zeros(n,1); newPcol(s_row) = 1; P = [P newPcol];
    end
    
    % Subsequent iterations
    for l = 2:m
        M = P'*U;
        inv_M = pinv(M);
        c = inv_M * (P'*Q(:,l));
        r = Q(:,l) - U*c;
        U = [U Q(:,l)];
        nsi = iif(l-1 < ns_mod_nr, floor(ns/m)+1, floor(ns/m));
        for i = 1:nsi
            s_row = -1; s_row_val = -Inf;
            for j = 1:n
                if ~any(used == j)
                    if s_row == -1 || s_row_val < abs(r(j))
                        s_row = j; s_row_val = abs(r(j));
                    end
                end
            end
            used(end+1) = s_row; phi(end+1) = s_row;
            newPcol = zeros(n,1); newPcol(s_row) = 1; P = [P newPcol];
        end
    end
end

%% --- Optimized Implementation ---
function [phi] = gnat_optimized(Q, m_used, nsr)
    [n, total_m] = size(Q);
    m = min(m_used, total_m);
    ns = iif(nsr > 0, nsr, m);
    
    phi = zeros(1, ns);
    is_used = false(n, 1);
    base_nsi = floor(ns / m);
    ns_mod_nr = mod(ns, m);
    U = []; 
    curr_phi_idx = 0;

    for l = 1:m
        if l == 1
            r = Q(:, 1);
        else
            % Vectorized sampling indexing instead of P matrix multiplication
            sampled_indices = phi(1:curr_phi_idx);
            c = U(sampled_indices, :) \ Q(sampled_indices, l);
            r = Q(:, l) - U * c;
        end
        U = [U, Q(:, l)];
        nsi = base_nsi + (l <= ns_mod_nr);
        for i = 1:nsi
            temp_r = abs(r);
            temp_r(is_used) = -1; 
            [~, s_row] = max(temp_r);
            curr_phi_idx = curr_phi_idx + 1;
            phi(curr_phi_idx) = s_row;
            is_used(s_row) = true;
        end
    end
end

function val = iif(cond, t, f), if cond, val = t; else, val = f; end, end
