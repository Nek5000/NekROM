function p = gpode(U, m)
    % Get initial points using QR decomposition with pivoting
    [~, ~, p_init] = qr(U', 'vector');
    
    n_vars = size(U, 2);
    p = zeros(m, 1); % Pre-allocate for speed
    p(1:n_vars) = p_init(1:n_vars);

    % Pre-transpose U for faster column access
    Ut = U'; 

    for i = (n_vars + 1):m
        % Economy SVD of the currently selected rows
        [~, S, W] = svd(U(p(1:i-1), :), 0);
        
        % Square the singular values once
        s_prev = S(end-1, end-1)^2;
        s_curr = S(end, end)^2;
        g = s_prev - s_curr;
        
        % Projected coordinates
        Ub = W' * Ut;
        Ub_sq = sum(Ub.^2, 1);
        
        % Compute r efficiently
        % Logic: r = (g + Ub_sq) - sqrt((g + Ub_sq)^2 - 4*g*Ub_last^2)
        term1 = g + Ub_sq;
        r = term1 - sqrt(term1.^2 - 4 * g * Ub(end, :).^2);
        
        % Find the best candidate not already in p
        [~, I] = sort(r, 'descend');
        
        % Find the first index in I that is not already in p(1:i-1)
        % setdiff or logical indexing is cleaner than a while loop
        available_indices = I(~ismember(I, p(1:i-1)));
        p(i) = available_indices(1);
    end
end
