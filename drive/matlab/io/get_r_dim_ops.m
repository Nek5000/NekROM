function [a, a0, b, c, c0, c1, c2, c3, u0, uk, ukmin, ukmax] = get_r_dim_ops(au_full, bu_full, cu_full, u0_full, uk_full, nb)
    % Linear Dynamics extraction
    a  = au_full(2:nb+1, 2:nb+1);
    a0 = au_full(2:nb+1, 1);
    b  = bu_full(2:nb+1, 2:nb+1);

    % Bilinear/Quadratic Terms (C-tensors)
    % Extracting slices directly to save memory
    c0 = reshape(cu_full(1:nb, 1:nb+1, 1:nb+1), nb*(nb+1), nb+1);
    c1 = cu_full(1:nb, 1, 1);
    c2 = reshape(cu_full(1:nb, 1, 1:nb+1), nb, nb+1);
    c3 = reshape(cu_full(1:nb, 1:nb+1, 1), nb, nb+1);

    % Core C-tensor reshaping
    c  = reshape(cu_full(1:nb, 2:nb+1, 2:nb+1), nb*nb, nb);

    % State/Control extraction
    u0 = u0_full(1:nb+1);
    uk = uk_full(1:nb+1, :);
    
    % Boundary calculations
    ukmin = min(uk, [], 2);
    ukmax = max(uk, [], 2);
end
