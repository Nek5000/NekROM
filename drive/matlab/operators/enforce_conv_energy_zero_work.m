function [c_fixed, info] = enforce_conv_energy_zero_work(c, u_state, mass_mat, opts)
    % ENFORCE_CONV_ENERGY_ZERO_WORK Enforce u' * B * c(u) = 0 by projection.
    %
    % Given a reduced convection vector `c` (size nb) and state coefficients
    % `u_state` (size nb), project `c` so it does no work in the B-inner product:
    %   u_state' * mass_mat * c_fixed = 0.
    %
    % This is a cheap, online-only correction that improves robustness when
    % hyper-reduction breaks the skew-adjointness/energy property of convection.

    if nargin < 4
        opts = struct();
    end
    if ~isstruct(opts)
        error('enforce_conv_energy_zero_work:InvalidOpts', 'opts must be a struct.');
    end
    if ~isfield(opts, 'mass_cross')
        opts.mass_cross = [];
    end
    if ~isfield(opts, 'mean_coef')
        opts.mean_coef = 0;
    end

    u_state = u_state(:);
    c = c(:);
    nb = numel(c);
    assert(numel(u_state) == nb, 'enforce_conv_energy_zero_work:DimMismatch', ...
        'u_state must have length %d.', nb);
    assert(all(size(mass_mat) == [nb, nb]), 'enforce_conv_energy_zero_work:MassMismatch', ...
        'mass_mat must be %dx%d.', nb, nb);

    mass_cross = opts.mass_cross(:);
    if ~isempty(mass_cross)
        assert(numel(mass_cross) == nb, 'enforce_conv_energy_zero_work:CrossMismatch', ...
            'opts.mass_cross must have length %d.', nb);
    end

    Bu = mass_mat * u_state;
    denom = u_state' * Bu;
    if ~isempty(mass_cross)
        denom = denom + (opts.mean_coef * (mass_cross' * u_state));
    end

    info = struct();
    info.applied = false;
    info.alpha = 0;
    info.work_before = u_state' * (mass_mat * c);
    if ~isempty(mass_cross)
        info.work_before = info.work_before + (opts.mean_coef * (mass_cross' * c));
    end
    info.work_after = info.work_before;
    info.denom = denom;

    if ~isfinite(denom) || denom <= eps
        c_fixed = c;
        return;
    end

    alpha = info.work_before / denom;
    c_fixed = c - alpha * u_state;

    info.applied = true;
    info.alpha = alpha;
    info.work_after = u_state' * (mass_mat * c_fixed);
    if ~isempty(mass_cross)
        info.work_after = info.work_after + (opts.mean_coef * (mass_cross' * c_fixed));
    end
end
