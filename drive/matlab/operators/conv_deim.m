function out_coef = conv_deim(ucoef, rom_data, method)
    % APPLY_ROM_CONVECTION Unified runtime for DEIM, CLSDEIM, and MCLSDEIM
    % method: 'deim', 'cls', or 'mcls'
    
    if nargin < 3, method = 'deim'; end

    u_p = rom_data.u_p;
    v_p = rom_data.v_p;
    ux_p = rom_data.ux_p;
    uy_p = rom_data.uy_p;
    eval_weights = [];
    if isfield(rom_data, 'eval_u_p')
        u_p = rom_data.eval_u_p;
        v_p = rom_data.eval_v_p;
        ux_p = rom_data.eval_ux_p;
        uy_p = rom_data.eval_uy_p;
    end
    if isfield(rom_data, 'eval_weights')
        eval_weights = rom_data.eval_weights;
    end

    % --- 1. Evaluate Nonlinearity at DEIM points only ---
    % Convection: (u * du/dx) + (v * du/dy)
    f_raw = (u_p(:, 2:end) * ucoef(2:end)) .* (ux_p(:, 2:end) * ucoef(2:end)) + ...
            (v_p(:, 2:end) * ucoef(2:end)) .* (uy_p(:, 2:end) * ucoef(2:end));
    if isempty(eval_weights)
        f_p = f_raw;
    else
        f_p = eval_weights .* f_raw;
    end

    % --- 2. Solve for Coefficients ---
    switch lower(method)
        case 'deim'
            % Standard DEIM
            c_hat = rom_data.interp_mat * f_p;
            
        case 'clsdeim'
            % Constrained Least Squares
            c_hat = rom_data.interp_mat * f_p;
            b = rom_data.proj_mat' * ucoef(2:end);
            lambda = (b' * c_hat) / (b' * rom_data.Ainv * b);
            c_hat = c_hat - lambda * (rom_data.Ainv * b);

        case 'mclsdeim'
            % Modified Constrained LS (Regularized with Snapshot Stats)
            rhs = (rom_data.nl_bas_p_eval' * f_p) + (rom_data.alpha * rom_data.tau * rom_data.mu);
            c_hat = rom_data.A_tau_inv * rhs;
            
            % Enforce Linear Constraint
            b = rom_data.proj_mat' * ucoef(2:end);
            lambda = (b' * c_hat) / (b' * rom_data.A_tau_inv * b);
            c_hat = c_hat - lambda * (rom_data.A_tau_inv * b);
            
        otherwise
            error('Invalid method. Use "deim", "clsdeim", or "mclsdeim".');
    end

    % --- 3. Final Mapping ---
    % Add the precomputed interaction of the zeroth mode
    out_coef = rom_data.proj_mat * c_hat + (rom_data.zmc * ucoef);
end
