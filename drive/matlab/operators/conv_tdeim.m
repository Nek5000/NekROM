function out_coef = conv_tdeim(ucoef, tcoef, rom_data, method)
    % CONV_TDEIM DEIM-family runtime for the temperature convection term u·grad(T).
    %
    % Inputs:
    % - ucoef: (nb+1)x1 velocity coefficients (mode0 + nb modes)
    % - tcoef: (nb+1)x1 temperature coefficients (mode0 + nb modes)
    % - rom_data: struct from load_tdeim_artifacts
    % - method: 'deim', 'clsdeim', or 'mclsdeim'
    %
    % Output:
    % - out_coef: nbx1 coefficient vector for Pi{u·grad(T)} (excluding mode0)

    if nargin < 4 || isempty(method)
        method = 'deim';
    end

    u_p = rom_data.u_p;
    v_p = rom_data.v_p;
    w_p = [];
    if isfield(rom_data, 'w_p')
        w_p = rom_data.w_p;
    end
    tx_p = rom_data.tx_p;
    ty_p = rom_data.ty_p;
    tz_p = [];
    if isfield(rom_data, 'tz_p')
        tz_p = rom_data.tz_p;
    end

    eval_weights = [];
    if isfield(rom_data, 'eval_weights')
        eval_weights = rom_data.eval_weights;
    end

    % --- 1. Evaluate nonlinearity at TDEIM points only ---
    up = u_p(:, 2:end) * ucoef(2:end);
    vp = v_p(:, 2:end) * ucoef(2:end);
    wp = 0.0;
    if ~isempty(w_p)
        wp = w_p(:, 2:end) * ucoef(2:end);
    end

    tx = tx_p(:, 2:end) * tcoef(2:end);
    ty = ty_p(:, 2:end) * tcoef(2:end);
    tz = 0.0;
    if ~isempty(tz_p)
        tz = tz_p(:, 2:end) * tcoef(2:end);
    end

    f_raw = (up .* tx) + (vp .* ty) + (wp .* tz);
    if isempty(eval_weights)
        f_p = f_raw;
    else
        f_p = eval_weights .* f_raw;
    end

    % --- 2. Solve for coefficients ---
    switch lower(method)
        case 'deim'
            c_hat = rom_data.interp_mat * f_p;

        case 'clsdeim'
            c_hat = rom_data.interp_mat * f_p;
            b = rom_data.proj_mat' * tcoef(2:end);
            denom = (b' * (rom_data.Ainv * b));
            if abs(denom) > 0
                lambda = (b' * c_hat) / denom;
                c_hat = c_hat - lambda * (rom_data.Ainv * b);
            end

        case 'mclsdeim'
            rhs = (rom_data.nl_bas_p_eval' * f_p) + (rom_data.alpha * rom_data.tau * rom_data.mu);
            c_hat = rom_data.A_tau_inv * rhs;

            b = rom_data.proj_mat' * tcoef(2:end);
            denom = (b' * (rom_data.A_tau_inv * b));
            if abs(denom) > 0
                lambda = (b' * c_hat) / denom;
                c_hat = c_hat - lambda * (rom_data.A_tau_inv * b);
            end

        otherwise
            error('Invalid method. Use "deim", "clsdeim", or "mclsdeim".');
    end

    % --- 3. Final mapping (plus exact mode-0 coupling) ---
    out_coef = rom_data.proj_mat * c_hat + (rom_data.zmc_u * ucoef) + (rom_data.zmc_t * tcoef);
end
