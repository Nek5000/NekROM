% The original gpode code is in `extra` folder.
% This is the result after querying Gemini for
% improvements.
function p = gpode(U, m)
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
end
