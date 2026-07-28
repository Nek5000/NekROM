% The original gpode code is in `extra` folder.
% This is the result after querying Gemini for
% improvements.
% Original paper:
% Stabilizing discrete empirical interpolation via randomized and deterministic oversampling
% by Peherstorfer et. al
% (this is the GPOD+E algorithm, essentially qdeim with oversampling.)
function p = qdeim_adaptive(U, m, tol)
    [n_samples, n_vars] = size(U);
    p = zeros(m, 1);
    mask = true(n_samples, 1); 

    % 1. Initial Selection
    [~, ~, p_init] = qr(U', 'vector');
    p(1:n_vars) = p_init(1:n_vars);
    mask(p(1:n_vars)) = false;

    % 2. Initial SVD
    Ut = U'; 
    [~, S, W] = svd(Ut(:, p(1:n_vars))', 0);
    Ub = W' * Ut; 

    % --- AUTOMATIC TOLERANCE LOGIC ---
    % If no tol is provided, use machine epsilon scaled by the 
    % largest singular value squared (the energy of the first mode).
    if nargin < 3 || isempty(tol)
        tol = eps(S(1,1)^2) * n_samples; 
    end
    % ---------------------------------

    actual_m = m; 

    for i = (n_vars + 1):m
        s_prev = S(end-1, end-1)^2;
        s_curr = S(end, end)^2;
        g = s_prev - s_curr;

        Ub_sq = sum(Ub.^2, 1);
        term1 = g + Ub_sq;
        
        % Objective function
        r = term1 - sqrt(max(0, term1.^2 - 4 * g * Ub(end, :).^2));
        
        % Find best candidate among unselected indices
        r_candidates = r;
        r_candidates(~mask) = -inf;
        [max_val, new_p] = max(r_candidates);

        % 3. Adaptive Early Stopping
        % If the improvement is below our threshold, stop oversampling.
        if max_val < tol
            actual_m = i - 1;
            p = p(1:actual_m);
            break;
        end

        p(i) = new_p;
        mask(new_p) = false;

        % 4. Incremental Update
        u_new_proj = Ub(:, new_p); 
        [~, S_new, W_rot] = svd([S; u_new_proj'], 0);
        S = S_new(1:n_vars, :); 
        Ub = W_rot' * Ub;
    end
end
