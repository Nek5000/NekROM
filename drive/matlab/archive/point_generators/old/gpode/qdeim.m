function p = qdeim_final(U, m, tol)
    if nargin < 3, tol = 1e-10; end % Default tolerance
    
    [n_samples, n_vars] = size(U);
    p = zeros(m, 1);
    mask = true(n_samples, 1); 

    % 1. Initial Selection (QR with pivoting)
    [~, ~, p_init] = qr(U', 'vector');
    p(1:n_vars) = p_init(1:n_vars);
    mask(p(1:n_vars)) = false;

    % 2. Initial SVD and Projection
    Ut = U'; 
    [~, S, W] = svd(Ut(:, p(1:n_vars))', 0);
    Ub = W' * Ut; 

    actual_m = m; % Track how many points we actually pick

    for i = (n_vars + 1):m
        % Calculate 'g' and objective 'r'
        s_prev = S(end-1, end-1)^2;
        s_curr = S(end, end)^2;
        g = s_prev - s_curr;

        Ub_sq = sum(Ub.^2, 1);
        term1 = g + Ub_sq;
        r = term1 - sqrt(max(0, term1.^2 - 4 * g * Ub(end, :).^2));
        
        % Early Stopping Check: If the max gain is negligible
        [max_val, new_p] = max(r .* mask');
        if max_val < tol
            actual_m = i - 1;
            p = p(1:actual_m);
            fprintf('Convergence reached at %d points.\n', actual_m);
            break;
        end

        p(i) = new_p;
        mask(new_p) = false;

        % 3. Incremental SVD Update
        u_new_proj = Ub(:, new_p); 
        
        % We update the SVD of the chosen rows: [S; u_new_proj']
        % W_rot is (n_vars x n_vars) because it re-aligns the n_vars columns
        [~, S_new, W_rot] = svd([S; u_new_proj'], 0);
        
        % S_new is ((n_vars+1) x n_vars). Truncate to keep n_vars.
        S = S_new(1:n_vars, :); 
        
        % Rotate the projections to match the new basis alignment
        Ub = W_rot' * Ub;
    end
%{
function p = qdeim(U, m)
    [n_samples, n_vars] = size(U);
    p = zeros(m, 1);
    
    % Initial selection via QR pivoting
    [~, ~, p_init] = qr(U', 'vector');
    p(1:n_vars) = p_init(1:n_vars);
    
    % Pre-transpose for column-major efficiency
    Ut = U'; 
    
    % Initial decomposition of the first n_vars rows
    % We maintain the SVD of the currently selected subset
    [~, S, W] = svd(U(p(1:n_vars), :), 0);
    
    for i = (n_vars + 1):m
        % 1. Compute the critical values from the current SVD
        s_prev = S(end-1, end-1)^2;
        s_curr = S(end, end)^2;
        g = s_prev - s_curr;
        
        % 2. Project all candidate rows onto the current basis W
        Ub = W' * Ut; % Efficient projection
        Ub_sq = sum(Ub.^2, 1);
        
        % 3. Calculate the objective function 'r'
        term1 = g + Ub_sq;
        r = term1 - sqrt(max(0, term1.^2 - 4 * g * Ub(end, :).^2));
        
        % 4. Find the best new index (avoiding duplicates)
        [~, I] = sort(r, 'descend');
        idx = 1;
        while any(p(1:i-1) == I(idx))
            idx = idx + 1;
        end
        new_p = I(idx);
        p(i) = new_p;
        
        % 5. INCREMENTAL UPDATE: Update S and W for the next iteration
        % We append the new row u_new = U(new_p, :) to our set.
        % This is a rank-1 update to the Gramian or a row-append to the SVD.
        u_new = U(new_p, :);
        
        % We use a standard SVD "row-addition" update logic:
        % Update W and S by performing SVD on the small [S; u_new * W] matrix
        [~, S, W_new] = svd([S; u_new * W], 0);
        W = W * W_new; % Rotate the basis
    end
%}
end
