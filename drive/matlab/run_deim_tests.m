function run_deim_tests(snaps_path, casename, reorder, ndeim_pts, n_os_points, ps_alg)
    % RUN_DEIM_TESTS Batch sanity checks for coarse and fine-grid DEIM setup.
    %
    % This is safe to invoke with `matlab -batch "run_deim_tests"` from the
    % drive/matlab directory. When called without arguments, it loads config.m.

    if nargin ~= 0 && nargin ~= 6
        error('run_deim_tests expects either zero or six arguments.');
    end

    if nargin == 0
        clear classes;
    end
    close all;
    clc;

    original_dir = pwd;
    cwd_cleanup = onCleanup(@() cd(original_dir));

    driver_dir = fileparts(mfilename('fullpath'));
    cd(driver_dir);

    toolkit_dir = fullfile(driver_dir, '..', '..', '..', 'NekToolKit', 'matlab');
    if exist(toolkit_dir, 'dir')
        addpath(toolkit_dir, '-begin');
    end
    addpath(fullfile(driver_dir, 'io'));
    addpath(fullfile(driver_dir, 'operators'));
    addpath(fullfile(driver_dir, 'point_generators'));

    deim_dealias = false;
    deim_dealias_cquad = false;
    deim_dealias_quad = false;
    if nargin == 0
        config;
        reorder = 1;
    end
    if ~exist('deim_alpha', 'var') || isempty(deim_alpha)
        deim_alpha = 1e-12;
    end

    bas_obj = NekSnaps(fullfile(snaps_path, strcat('bas', casename)));
    [pod_u, pod_v] = get_snaps(bas_obj, reorder);
    [x_base, y_base] = get_grid(bas_obj, reorder);
    rom_nb = size(pod_u, 2) - 1;

    nl_bas_obj = NekSnaps(fullfile(snaps_path, strcat('cba', casename)));
    [nl_bas_u, nl_bas_v] = get_snaps(nl_bas_obj, reorder);
    [x_nl, y_nl] = get_grid(nl_bas_obj, reorder);

    nl_snaps_obj = NekSnaps(fullfile(snaps_path, strcat('csn', casename)));
    [nl_snaps_u, nl_snaps_v] = get_snaps(nl_snaps_obj, reorder);

    methods = {'deim', 'clsdeim', 'mclsdeim'};
    grid_modes = [false, true];

    for deim_finegrid = grid_modes
        if deim_finegrid
            [pod_u_grid, pod_v_grid] = interp_basis_to_grid(pod_u, pod_v, x_base, y_base, x_nl, y_nl);
            x_grid = x_nl;
            y_grid = y_nl;
            nl_bas_grid = [nl_bas_u; nl_bas_v];
            nl_snaps_u_grid = nl_snaps_u;
            nl_snaps_v_grid = nl_snaps_v;
        else
            pod_u_grid = pod_u;
            pod_v_grid = pod_v;
            x_grid = x_base;
            y_grid = y_base;
            [nl_bas_u_grid, nl_bas_v_grid] = interp_basis_to_grid(nl_bas_u, nl_bas_v, x_nl, y_nl, x_grid, y_grid);
            [nl_snaps_u_grid, nl_snaps_v_grid] = interp_basis_to_grid(nl_snaps_u, nl_snaps_v, x_nl, y_nl, x_grid, y_grid);
            nl_bas_grid = [nl_bas_u_grid; nl_bas_v_grid];
        end

        if deim_finegrid
            grid_name = 'fine';
        else
            grid_name = 'coarse';
        end
        fprintf('DEIM sanity checks on %s grid\n', grid_name);

        for im = 1:numel(methods)
            method = methods{im};
            if strcmp(method, 'mclsdeim')
                snaps_u = nl_snaps_u_grid;
                snaps_v = nl_snaps_v_grid;
            else
                snaps_u = [];
                snaps_v = [];
            end

            rom_data = setup_conv_deim( ...
                pod_u_grid, pod_v_grid, nl_bas_grid, snaps_u, snaps_v, ...
                x_grid, y_grid, ndeim_pts, n_os_points, ps_alg, deim_dealias, deim_dealias_quad, deim_alpha, deim_dealias_cquad);

            assert(size(rom_data.inds, 1) == ndeim_pts, 'Incorrect DEIM point count.');
            assert(numel(unique(rom_data.inds)) == ndeim_pts, 'DEIM points must be unique.');
            assert(max(rom_data.inds) <= size(nl_bas_grid, 1), 'DEIM points out of bounds.');
            assert(size(rom_data.u_p, 1) == ndeim_pts, 'Incorrect u_p row count.');
            assert(size(rom_data.u_p, 2) == rom_nb + 1, 'Incorrect u_p column count.');
            assert(size(rom_data.v_p, 1) == ndeim_pts, 'Incorrect v_p row count.');
            assert(size(rom_data.ux_p, 1) == ndeim_pts, 'Incorrect ux_p row count.');
            assert(size(rom_data.uy_p, 1) == ndeim_pts, 'Incorrect uy_p row count.');
            assert(size(rom_data.proj_mat, 1) == rom_nb, 'Incorrect projection matrix row count.');
            assert(size(rom_data.Ainv, 1) == size(rom_data.Ainv, 2), 'Ainv must be square.');

            if deim_dealias_quad
                assert(rom_data.use_full_quadrature, 'Strict quadrature path was not enabled.');
                assert(size(rom_data.eval_u_p, 1) == size(rom_data.eval_weights, 1), 'Quadrature evaluation rows must match weights.');
                assert(size(rom_data.eval_u_p, 1) == size(rom_data.nl_bas_p_eval, 1), 'Quadrature basis rows must match evaluation rows.');
                assert(all(rom_data.eval_weights > 0), 'Quadrature weights must be positive.');
            elseif deim_dealias_cquad
                assert(isfield(rom_data, 'use_compressed_quadrature') && rom_data.use_compressed_quadrature, ...
                    'Compressed quadrature path was not enabled.');
                assert(size(rom_data.eval_u_p, 1) == size(rom_data.eval_weights, 1), 'CQuad evaluation rows must match weights.');
                assert(size(rom_data.eval_u_p, 1) == size(rom_data.nl_bas_p_eval, 1), 'CQuad basis rows must match evaluation rows.');
                assert(all(rom_data.eval_weights > 0), 'CQuad weights must be positive.');
            elseif deim_dealias
                expected_pts = ndeim_pts + n_os_points;
                assert(rom_data.use_oversampled_points, 'Oversampled DEIM path was not enabled.');
                assert(size(rom_data.inds_os, 1) == expected_pts, 'Incorrect oversampled DEIM point count.');
                assert(size(rom_data.eval_u_p, 1) == expected_pts, 'Incorrect oversampled evaluation row count.');
                assert(size(rom_data.eval_ux_p, 1) == expected_pts, 'Incorrect oversampled gradient row count.');
                assert(size(rom_data.nl_bas_p_eval, 1) == expected_pts, 'Incorrect oversampled basis row count.');
            end

            if strcmp(method, 'mclsdeim')
                assert(isfield(rom_data, 'mu') && isfield(rom_data, 'tau') && isfield(rom_data, 'A_tau_inv'), ...
                    'MCLS-DEIM statistics are missing.');
                assert(size(rom_data.A_tau_inv, 1) == size(rom_data.A_tau_inv, 2), ...
                    'A_tau_inv must be square.');
            end

            if size(rom_data.proj_mat, 2) >= 1 && norm(rom_data.proj_mat(:, 1)) > 0
                test_ucoef = [1; rom_data.proj_mat(:, 1)];
            else
                test_ucoef = [1; ones(rom_nb, 1)];
            end
            out = conv_deim(test_ucoef, rom_data, method);
            assert(numel(out) == rom_nb, 'Unexpected DEIM output size.');
            assert(all(isfinite(out)), 'DEIM output contains non-finite values.');

            fprintf('PASS: %s on %s grid\n', method, grid_name);
        end
    end

    fprintf('DEIM sanity checks complete.\n');
end
