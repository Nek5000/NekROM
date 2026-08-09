function rom_data = load_deim_artifacts(ops_dir, nb, nbnl)
% LOAD_DEIM_ARTIFACTS Load Fortran-compatible DEIM-family artifacts from ops/.
%
% This loader understands the plain-text format written by NekROM Fortran:
% - Scalars and vectors are written one value per line
% - Matrices are written in column-major order (Fortran/MATLAB linearization)
%
% The returned rom_data struct matches the fields expected by conv_deim.m.

    if nargin < 1 || isempty(ops_dir)
        error('load_deim_artifacts:MissingOpsDir', 'ops_dir is required.');
    end
    if nargin < 2 || isempty(nb)
        error('load_deim_artifacts:MissingNb', 'nb is required.');
    end
    if nargin < 3 || isempty(nbnl)
        error('load_deim_artifacts:MissingNbnl', 'nbnl is required.');
    end

    ndeim_pts = read_int_scalar(fullfile(ops_dir, 'deim_npts'));
    ndeim_pts_os = read_int_scalar(fullfile(ops_dir, 'deim_npts_os'));
    ndeim_pts_eval = read_int_scalar(fullfile(ops_dir, 'deim_npts_eval'));

    % Fortran may clamp the requested nbnl to an internal maximum and will
    % then write ops/deim_* with the *effective* nbnl. The effective size
    % is consistent with ndeim_pts.
    nbnl_eff = min(nbnl, ndeim_pts);
    if nbnl_eff ~= nbnl
        warning('load_deim_artifacts:ClampedNbnl', ...
            'Requested nbnl=%d but ops/deim_npts=%d; using nbnl_eff=%d for ops/deim_* loading.', ...
            nbnl, ndeim_pts, nbnl_eff);
    end

    rom_data = struct();
    rom_data.inds = read_int_vector(fullfile(ops_dir, 'deim_inds'), ndeim_pts);
    rom_data.inds_os = read_int_vector_optional(fullfile(ops_dir, 'deim_inds_os'), ndeim_pts_os);
    rom_data.eval_inds = read_int_vector(fullfile(ops_dir, 'deim_eval_inds'), ndeim_pts_eval);
    rom_data.eval_weights = read_real_vector(fullfile(ops_dir, 'deim_eval_weights'), ndeim_pts_eval);

    % Pointwise POD basis and gradients at evaluation points.
    rom_data.eval_u_p = read_real_matrix(fullfile(ops_dir, 'deim_u_p'), ndeim_pts_eval, nb + 1);
    rom_data.eval_v_p = read_real_matrix(fullfile(ops_dir, 'deim_v_p'), ndeim_pts_eval, nb + 1);
    rom_data.eval_ux_p = read_real_matrix(fullfile(ops_dir, 'deim_ux_p'), ndeim_pts_eval, nb + 1);
    rom_data.eval_uy_p = read_real_matrix(fullfile(ops_dir, 'deim_uy_p'), ndeim_pts_eval, nb + 1);

    % Provide non-eval aliases so conv_deim works regardless of which fields are present.
    rom_data.u_p = rom_data.eval_u_p;
    rom_data.v_p = rom_data.eval_v_p;
    rom_data.ux_p = rom_data.eval_ux_p;
    rom_data.uy_p = rom_data.eval_uy_p;

    % DEIM mapping operators.
    rom_data.nl_bas_p_eval = read_real_matrix(fullfile(ops_dir, 'deim_nl_bas_p_eval'), ndeim_pts_eval, nbnl_eff);
    rom_data.proj_mat = read_real_matrix(fullfile(ops_dir, 'deim_proj_mat'), nb, nbnl_eff);
    rom_data.zmc = read_real_matrix(fullfile(ops_dir, 'deim_zmc'), nb, nb + 1);
    rom_data.Ainv = read_real_matrix(fullfile(ops_dir, 'deim_Ainv'), nbnl_eff, nbnl_eff);
    rom_data.interp_mat = read_real_matrix(fullfile(ops_dir, 'deim_interp_mat'), nbnl_eff, ndeim_pts_eval);

    % Optional MCLSDEIM artifacts.
    mu_path = fullfile(ops_dir, 'deim_mu');
    tau_path = fullfile(ops_dir, 'deim_tau');
    a_tau_inv_path = fullfile(ops_dir, 'deim_A_tau_inv');
    alpha_path = fullfile(ops_dir, 'deim_alpha');

    if exist(mu_path, 'file')
        rom_data.mu = read_real_vector(mu_path, nbnl_eff);
    end
    if exist(tau_path, 'file')
        rom_data.tau = read_real_matrix(tau_path, nbnl_eff, nbnl_eff);
    end
    if exist(a_tau_inv_path, 'file')
        rom_data.A_tau_inv = read_real_matrix(a_tau_inv_path, nbnl_eff, nbnl_eff);
    end
    if exist(alpha_path, 'file')
        rom_data.alpha = read_real_scalar(alpha_path);
    end
end

function value = read_real_scalar(path)
    data = dlmread(path);
    if isempty(data)
        error('load_deim_artifacts:EmptyFile', 'Empty scalar file: %s', path);
    end
    value = data(1);
end

function value = read_int_scalar(path)
    data = dlmread(path);
    if isempty(data)
        error('load_deim_artifacts:EmptyFile', 'Empty int scalar file: %s', path);
    end
    value = round(data(1));
end

function vec = read_real_vector(path, n)
    data = dlmread(path);
    if numel(data) < n
        error('load_deim_artifacts:ShortVector', ...
            'Expected %d entries in %s, got %d.', n, path, numel(data));
    end
    vec = reshape(data(1:n), [n, 1]);
end

function vec = read_int_vector(path, n)
    data = dlmread(path);
    if numel(data) < n
        error('load_deim_artifacts:ShortVector', ...
            'Expected %d entries in %s, got %d.', n, path, numel(data));
    end
    vec = reshape(round(data(1:n)), [n, 1]);
end

function vec = read_int_vector_optional(path, n)
    if n <= 0 || ~exist(path, 'file')
        vec = zeros(0, 1);
        return;
    end
    vec = read_int_vector(path, n);
end

function mat = read_real_matrix(path, m1, m2)
    data = dlmread(path);
    expected = m1 * m2;
    if numel(data) < expected
        error('load_deim_artifacts:ShortMatrix', ...
            'Expected %d entries in %s (%dx%d), got %d.', expected, path, m1, m2, numel(data));
    end
    mat = reshape(data(1:expected), [m1, m2]);
end
