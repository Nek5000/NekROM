% driver.m - MATLAB/Octave Driver for Galerkin-based reduced order model 
function results = driver()
    close all; clc;
    wall_tic = tic;
    original_dir = pwd;
    cwd_cleanup = onCleanup(@() cd(original_dir));

    % Anchor all relative paths to this driver, not to the launch directory.
    driver_dir = fileparts(mfilename('fullpath'));
    cd(driver_dir);

    % --- OCTAVE COMPATIBILITY LAYER ---
    if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
        try
            pkg load optim;
        catch
            warning('Octave ''optim'' package not found. Optimization features (ifcopt) will fail.');
        end
    end
    % ----------------------------------

    % Setup paths relative to the driver location.
    addpath(fullfile(driver_dir, 'point_generators'));
    addpath(fullfile(driver_dir, 'io'));
    addpath(fullfile(driver_dir, 'operators'));
    clear functions;  % reset persistent helper state between runs

    % Add NekToolKit to path (required dependency)
    toolkit_dir = fullfile(driver_dir, '..', '..', '..', 'NekToolKit', 'matlab');
    if exist(toolkit_dir, 'dir')
        addpath(toolkit_dir, '-begin');
    end

    % Check for critical NekToolKit functions (all required)
    required_nektoolkit = {'NekSnaps', 'zwgll', 'deriv_mat', 'deriv_geo', 'grad', 'interp_mat'};
    missing_funcs = {};
    for i = 1:length(required_nektoolkit)
        if ~exist(required_nektoolkit{i}, 'file')
            missing_funcs{end+1} = required_nektoolkit{i};
        end
    end

    if ~isempty(missing_funcs)
        error(['NekToolKit dependency not satisfied.\n' ...
               'Missing functions: %s\n\n' ...
               'NekToolKit is required for NekROM MATLAB driver.\n' ...
               'Install from: https://github.com/kent0/NekToolKit\n\n' ...
               'Quick install:\n' ...
               '  cd %s\n' ...
               '  git clone https://github.com/kent0/NekToolKit.git\n\n' ...
               'For detailed diagnostics, run: check_dependencies()'], ...
               strjoin(missing_funcs, ', '), ...
               fullfile(driver_dir, '..', '..', '..'));
    end

    % Load simulation parameters
    config;

    % Validate case directory and required files exist
    if ~exist(case_path, 'dir')
        error(['Case directory not found: ', case_path, ...
               '\nCheck that thiscase=''%s'' is correct in config.m'], casename);
    end

    if ~exist(snaps_path, 'dir')
        error(['Snapshots directory not found: ', snaps_path, ...
               '\nRun offline phase first: cd ', case_path, ' && makerom ', casename]);
    end

    % Check for POD basis snapshots
    bas_file = fullfile(snaps_path, strcat('bas', casename, '0.f00001'));
    if ~exist(bas_file, 'file')
        error(['POD basis file not found: ', bas_file, ...
               '\nRun offline phase first: cd ', case_path, ' && makerom ', casename]);
    end

    % Check for offline operators before any expensive DEIM setup.
    ops_dir = fullfile(case_path, 'ops');
    if ~exist(ops_dir, 'dir')
        error(['Operators directory not found: ', ops_dir, ...
               '\nRun offline phase first: cd ', case_path, ' && makerom ', casename]);
    end

    required_ops = {'nb', 'au', 'bu', 'cu', 'u0', 'ns', 'uk'};
    missing_ops = {};
    for i = 1:length(required_ops)
        op_file = fullfile(ops_dir, required_ops{i});
        if ~exist(op_file, 'file')
            missing_ops{end+1} = op_file;
        end
    end
    if ~isempty(missing_ops)
        error(['Required operator file(s) not found:\n  %s\n' ...
               'Run offline phase first: cd ', case_path, ' && makerom ', casename], ...
               strjoin(missing_ops, '\n  '));
    end

    %% ROM Setup & Basis Generation
    reorder = 1; 
    cname = fullfile(snaps_path, strcat('bas', casename));
    bas_snaps = NekSnaps(cname); 
    [pod_u, pod_v] = get_snaps(bas_snaps, reorder);
    [x_fom, y_fom] = get_grid(bas_snaps, reorder);
    inde = bas_snaps.flds{1}.inde;
    pod_u_deim = pod_u;
    pod_v_deim = pod_v;
    x_deim = x_fom;
    y_deim = y_fom;
    nl_snaps_u_deim = [];
    nl_snaps_v_deim = [];

    % Define output directory
    casedir = sprintf('output/%s_%s/', casename, datestr(now, 'yyyy-mm-dd-HH-MM-SS'));
    if ~exist(casedir, 'dir'), mkdir(casedir); end
    basepath = fullfile(casedir, 'fields', casename);
    if ~exist(fullfile(casedir, 'fields'), 'dir'), mkdir(fullfile(casedir, 'fields')); end

    %% DEIM & Snapshot Initialization
    if ismember(conv_approach, {'deim', 'clsdeim', 'mclsdeim'})
        use_fortran_deim_ops = read_env_bool('NEKROM_DEIM_FROM_OPS', false);
        if use_fortran_deim_ops
            nbnl = case_meta.deim_nbnl;
            if isempty(nbnl) || ~isfinite(nbnl) || nbnl <= 0
                error(['NEKROM_DEIM_FROM_OPS=1 requires a valid [DEIM] nbnl entry in ' ...
                       case_meta.case_file, '.']);
            end
            deim_data = load_deim_artifacts(ops_dir, nb, nbnl);
        else
            % Use the generated cba* basis directly, or interpolate the ROM basis
            % onto that grid when the fine-grid DEIM path is enabled.
            matlab_pod_basis = 0;
            nl_bas_obj_nr = NekSnaps(fullfile(snaps_path, strcat('cba', casename)));
            [nl_bas_u_nr, nl_bas_v_nr] = get_snaps(nl_bas_obj_nr, reorder);
            [x_nl, y_nl] = get_grid(nl_bas_obj_nr, reorder);

            if deim_finegrid
                x_deim = x_nl;
                y_deim = y_nl;
                [pod_u_deim, pod_v_deim] = interp_basis_to_grid(pod_u, pod_v, x_fom, y_fom, x_deim, y_deim);
            else
                [nl_bas_u_nr, nl_bas_v_nr] = interp_basis_to_grid(nl_bas_u_nr, nl_bas_v_nr, x_nl, y_nl, x_deim, y_deim);
            end

            if matlab_pod_basis || strcmp(conv_approach, 'mclsdeim')
                nl_cname = fullfile(snaps_path, strcat('csn', casename));
                nl_snaps_obj = NekSnaps(nl_cname);
                [nl_snaps_u, nl_snaps_v] = get_snaps(nl_snaps_obj, reorder);
                if deim_finegrid
                    [nl_snaps_u_deim, nl_snaps_v_deim] = deal(nl_snaps_u, nl_snaps_v);
                else
                    [nl_snaps_u_deim, nl_snaps_v_deim] = interp_basis_to_grid(nl_snaps_u, nl_snaps_v, x_nl, y_nl, x_deim, y_deim);
                end
            else
                nl_snaps_u = [];
                nl_snaps_v = [];
            end

            if matlab_pod_basis
                [nl_bas, ~, ~] = get_pod_basis_from_arrays(nl_snaps_u_deim, nl_snaps_v_deim, x_deim, y_deim, ndeim_pts, 0, 0);
            else
                nl_bas = [nl_bas_u_nr; nl_bas_v_nr];
            end

            deim_data = setup_conv_deim( ...
                pod_u_deim, pod_v_deim, nl_bas, nl_snaps_u_deim, nl_snaps_v_deim, ...
                x_deim, y_deim, ndeim_pts, n_os_points, ps_alg, deim_dealias, deim_dealias_quad, deim_alpha);

            if deim_dealias_quad
                warning('NekROM:DEIMQuadNoPersist', ...
                    ['NEKROM_DEIM_DEALIAS_QUAD=1 uses an overintegrated MATLAB-only DEIM path. ' ...
                     'Skipping ops/ persistence because the current Fortran runtime cannot load it.']);
            else
                save_deim_artifacts(ops_dir, deim_data);
            end
        end
    end

    % Load FOM Operators
    [au_full, bu_full, cu_full, u0_full, uk_full, mb, ns] = load_full_ops(ops_dir);
    Me = get_Me(x_fom, y_fom);

    if nb ~= mb
        error('Configured nb (%d) does not match ops/nb (%d).', nb, mb);
    end
    if ~isempty(case_meta.ns) && ns ~= case_meta.ns
        error('Configured ns (%d) does not match ops/ns (%d).', case_meta.ns, ns);
    end

    % Call validation tests
    if if_run_tests
        ops_ips = case_meta.ips;
        if strcmpi(ops_ips, 'HLM')
            run_tests(snaps_path, casename, nb, reorder, pod_u, pod_v, au_full, bu_full, ops_ips, subtract_mean, 1 / nu, dt);
        else
            run_tests(snaps_path, casename, nb, reorder, pod_u, pod_v, au_full, bu_full, ops_ips, subtract_mean);
        end
    end

    % Get reduced dimensional operators
    [au, a0, bu, cu, c0, c1, c2, c3, u0, uk, ukmin, ukmax] = get_r_dim_ops(au_full, bu_full, cu_full, u0_full, uk_full, nb);

    %% Initialization
    time   = 0;
    rhs    = zeros(nb, 1);
    ext    = zeros(nb, 3);
    hufac  = [];
    nan_detected = false;
    aborted_step = NaN;

    % Preallocate outputs
    num_outputs = floor(nsteps / iostep);
    ucoef = zeros(num_outputs, nb+1);
    kes = zeros(num_outputs, 1);
    momentums = zeros(num_outputs, 2);
    io_idx = 1;

    if ifleray || ifefr || iftr
       dfHfac = [];
       dfHfac = set_df(au, bu, radius, 1, dfHfac);
    end

    % Set initial condition
    u = zeros(nb+1, 3); 
    u(:,1) = u0;
    [alphas, betas] = setcoef();

    u_proj = pod_u(:, 1:nb+1) * u(:, 1);
    v_proj = pod_v(:, 1:nb+1) * u(:, 1);

    field_data = struct('u', u_proj, 'v', v_proj, 'x', x_fom, 'y', y_fom, 'inde', inde, 'size', size(x_fom), 'time', 0.0, 'iostep', 0);
    output_fields(basepath, field_data, ifvort, ifwrite, ifvis); 

    %% Integrate ROM with BDFk/EXTk
    fprintf('Starting time integration: %d steps\n', nsteps);
    progress_interval = max(1, floor(nsteps / 20));  % Report every 5%

    for istep = 1:nsteps
        time = double(istep) * dt;
        ito = min(istep, 3);

        % Progress indicator (every 5%)
        if mod(istep, progress_interval) == 0
            pct = 100 * istep / nsteps;
            fprintf('Progress: %3.0f%% (step %d/%d, time=%.4f)\n', ...
                    pct, istep, nsteps, time);
        end
        
        if istep <= 3
            hufac = [];
        end

        ext(:,3) = ext(:,2);
        ext(:,2) = ext(:,1);

        if ifleray
            utmp = [1; (dfHfac \ (dfHfac' \ u(2:end, 1)))];
        else 
            utmp = u;
        end

        % Convection Logic
        switch conv_approach 
            case 'fom'
                c_coef = conv_fom(u(:,1), pod_u, pod_v, x_fom, y_fom, true);
            case 'ftensor'
                c_coef = (reshape(c0*utmp(:,1), nb, nb+1) * u(:,1));
            case 'rtensor'
                c_coef = conv_tensor_dense(u(:,1), pod_u, pod_v, x_fom, y_fom, tensor_size);
            case {'deim', 'clsdeim', 'mclsdeim'}
                c_coef = conv_deim(u(:,1), deim_data, conv_approach);
            otherwise
                error(['Unrecognized conv_approach: ', conv_approach]);
        end

        ext(:,1) = -c_coef - nu * a0;

        if iftr
            utmp_tr = [1; (dfHfac \ (dfHfac' \ u(2:end, 1)))];
            ext(:,1) = ext(:,1) - relax * (u(2:end, 1) - utmp_tr(2:end));
        end

        rhs = (ext * alphas(:, ito)) - bu * (u(2:end, :) * betas(2:end, ito)) / dt;

        % Solve Step
        if ifcopt
            [x, ~] = fmincon(@(x)rom_residual(x, au, bu, nu, betas, dt, ito, rhs), ...
                u(2:end, 1), [], [], [], [], ukmin(2:end), ukmax(2:end));
            u_new = [1; x];
        else
            if isempty(hufac)
                h = bu * betas(1, ito) / dt + au * nu;
                hfac = chol(h);
            end
            u_new = [1; (hfac \ (hfac' \ rhs))];
        end

        if ifefr
            utmp_efr = [1; (dfHfac \ (dfHfac' \ u_new(2:end)))];
            u_new = (1 - relax) * u_new + relax * utmp_efr;
        end
            
        u = shift(u, u_new, 3);

        if any(isnan(u(:,1)))
            fprintf('NaN detected at step %d. Aborting.\n', istep);
            nan_detected = true;
            aborted_step = istep;
            break;
        end

        % IO Routine
        if mod(istep, iostep) == 0
            % Removed verbose IOSTEP print (progress indicator handles this)
            ucoef(io_idx, :) = u(:, 1)';

            u_proj = pod_u(:, 1:nb+1) * u(:, 1);
            v_proj = pod_v(:, 1:nb+1) * u(:, 1);

            kes(io_idx) = 0.5 * (u_proj' * (Me .* u_proj) + v_proj' * (Me .* v_proj));
            momentums(io_idx, :) = [sum(Me .* u_proj), sum(Me .* v_proj)];

            field_data.u = u_proj;
            field_data.v = v_proj;
            field_data.time = time;
            field_data.iostep = floor(istep / iostep);

            output_fields(basepath, field_data, ifvort, ifwrite, ifvis);
            io_idx = io_idx + 1;
        end
    end

    %% Outputs & Visualization
    ucoef = ucoef(1:io_idx-1, :);
    kes = kes(1:io_idx-1);
    momentums = momentums(1:io_idx-1, :);

    fileID = fopen(fullfile(casedir, 'ucoef'), 'w');
    fprintf(fileID, '%24.15e\n', ucoef'); 
    fclose(fileID);

    if isempty(kes)
        ke_initial = NaN;
        ke_final = NaN;
        ke_rel_drift = NaN;
        ke_peak_ratio = NaN;
    else
        ke_initial = kes(1);
        ke_final = kes(end);
        ke_rel_drift = (ke_final - ke_initial) / max(abs(ke_initial), eps);
        ke_peak_ratio = max(kes) / max(abs(ke_initial), eps);
    end

    if isempty(momentums)
        max_momentum_norm = NaN;
        final_momentum_norm = NaN;
    else
        momentum_norms = sqrt(sum(momentums.^2, 2));
        max_momentum_norm = max(momentum_norms);
        final_momentum_norm = momentum_norms(end);
    end

    if isempty(ucoef)
        max_ucoef_norm = NaN;
    else
        max_ucoef_norm = max(sqrt(sum(ucoef(:, 2:end).^2, 2)));
    end

    completed = ~nan_detected && ((io_idx - 1) == num_outputs);

    results = struct();
    results.casename = casename;
    results.conv_approach = conv_approach;
    results.deim_finegrid = deim_finegrid;
    results.deim_dealias = deim_dealias;
    results.deim_dealias_quad = deim_dealias_quad;
    results.deim_dealias_mode = deim_dealias_mode;
    results.completed = completed;
    results.nan_detected = nan_detected;
    results.aborted_step = aborted_step;
    results.nsteps = nsteps;
    results.iostep = iostep;
    results.num_outputs = num_outputs;
    results.last_io_step = io_idx - 1;
    results.final_time = time;
    results.ke_initial = ke_initial;
    results.ke_final = ke_final;
    results.ke_rel_drift = ke_rel_drift;
    results.ke_peak_ratio = ke_peak_ratio;
    results.max_momentum_norm = max_momentum_norm;
    results.final_momentum_norm = final_momentum_norm;
    results.max_ucoef_norm = max_ucoef_norm;
    results.ucoef = ucoef;
    results.kes = kes;
    results.momentums = momentums;
    results.casedir = casedir;
    results.wall_time_sec = toc(wall_tic);

    save(fullfile(casedir, 'stability.mat'), 'results');

    if ifplot
        % Energy plot
        figure(2);
        plot(kes, 'LineWidth', 1.5); xlabel('Time'); ylabel('Kinetic energy');
        title('Energy Conservation');
        if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
            print(fullfile(casedir, 'ke.pdf'), '-dpdf');
        else
            exportgraphics(gca, fullfile(casedir, 'ke.pdf'), 'ContentType', 'vector');
        end

        % Momentum plot
        figure(3);
        plot(momentums(:,1), 'LineWidth', 1.5); hold on;
        plot(momentums(:,2), 'LineWidth', 1.5);
        xlabel('Time'); ylabel('Momentum components');
        title('Momentum Conservation');
        legend('U Momentum', 'V Momentum');
        if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
            print(fullfile(casedir, 'momentum.pdf'), '-dpdf');
        else
            exportgraphics(gca, fullfile(casedir, 'momentum.pdf'), 'ContentType', 'vector');
        end
    end

    disp('Simulation complete.');
end

%#####################################
%# Local Auxiliary functions
%#####################################

function hfac = set_df(a, b, dfRadius, dfOrder, hfac)
    if isempty(hfac)
        bfac = chol(b);
        h = (dfRadius^2) * (bfac \ (bfac' \ a));
        for i = 2:dfOrder
            h = h * ((dfRadius^2) * (bfac \ (bfac' \ a)));
        end
        h = h + eye(size(a));
        hfac = chol(h);
    end
end

function [alphas, betas] = setcoef()
    alphas = zeros(3, 3); betas  = zeros(4, 3);
    alphas(1,1) =  1.0;
    alphas(1,2) =  2.0; alphas(2,2) = -1.0;
    alphas(1,3) =  3.0; alphas(2,3) = -3.0; alphas(3,3) =  1.0;

    betas(1,1) =  1.0; betas(2,1) = -1.0;
    betas(1,2) =  1.5; betas(2,2) = -2.0; betas(3,2) =  0.5;
    betas(1,3) =  11.0/6; betas(2,3) = -3.0; betas(3,3) =  1.5; betas(4,3) = -1.0/3;
end

function a = shift(a, b, n)
    for i = n:-1:2 
        a(:,i) = a(:,i-1); 
    end
    a(:,1) = b;
end

function F = rom_residual(x, a, b, diff, betas, dt, ito, rhs)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    h = b * betas(1, ito) / dt + a * diff;
    F1 = h * x - rhs;
    F = norm(F1);
end
