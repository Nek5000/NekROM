function [Me] = get_Me(x, y)
    persistent cached_x cached_y cached_Me
    
    % Check if we've already computed this for the current x and y
    if isequal(x, cached_x) && isequal(y, cached_y)
        Me = cached_Me;
        return;
    end

    % --- Original Logic ---
    nx1 = size(x, 1);
    [zi, w] = zwgll(nx1-1);
    d = deriv_mat(zi);
    
    % Optimization: Only extract 'jac' if the other outputs aren't used
    [~, ~, ~, ~, ~, ~, ~, ~, jac, ~, ~] = deriv_geo(x, y, d);
    
    % Vectorized weight multiplication
    % Using (w * w') is an outer product; ensure w is a column vector
    Me = reshape(jac .* (w * w'), [], 1); 

    % Update Cache
    cached_x = x;
    cached_y = y;
    cached_Me = Me;
end
