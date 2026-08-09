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
    enforce_skew_adjoint = read_env_bool('NEKROM_CONV_ENFORCE_SKEW_ADJOINT', false);
    skew_inner = getenv('NEKROM_CONV_SKEW_INNER');
    if isempty(skew_inner)
        skew_inner = 'l2';
    else
        skew_inner = lower(strtrim(skew_inner));
    end
    skew_apply = getenv('NEKROM_CONV_SKEW_APPLY');
    if isempty(skew_apply)
        skew_apply = 'deim';
    else
        skew_apply = lower(strtrim(skew_apply));
    end

    % Validate case directory and required files exist
    if ~exist(case_path, 'dir')
        error(['Case directory not found: ', case_path, ...
               '\nCheck that thiscase=''%s'' is correct in config.m'], casename);
    end

    if ~exist(snaps_path, 'dir')
        % Some NekROM examples keep snapshots/bases in the case root or under
        % snaps_rom/ instead of snaps/. We discover snapshots via helper
        % functions later; don't hard-require snaps/ here.
        warning('NekROM:MissingSnapsDir', ...
            ['Snapshots directory not found: %s\n' ...
             'Proceeding with snapshot discovery in case root/snaps_rom/.'], snaps_path);
    end

    % Check for POD basis snapshots (support both snaps/ and case root layouts).
    bas_prefix = resolve_snapshot_prefix(case_path, casename, 'bas');
    if isempty(bas_prefix)
        error(['POD basis files not found under ', case_path, '.\n' ...
               'Expected either snaps/bas*0.f* or bas*0.f*.\n' ...
               'Run offline phase first: cd ', case_path, ' && makerom ', casename]);
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

    % Detect thermo-fluid operator bundles (temperature equation + optional TDEIM).
    thermal_ops = {'at', 'bt', 'ct', 't0', 'tk'};
    has_thermal_ops = true;
    for i = 1:length(thermal_ops)
        if ~exist(fullfile(ops_dir, thermal_ops{i}), 'file')
            has_thermal_ops = false;
            break;
        end
    end

	    %% ROM Setup & Basis Generation
	    reorder = true;
	    bas_snaps = NekSnaps(bas_prefix);
	    has_geom = isfield(bas_snaps.flds{1}, 'x') && isfield(bas_snaps.flds{1}, 'y') && ...
	               ~isempty(bas_snaps.flds{1}.x) && ~isempty(bas_snaps.flds{1}.y);

	    if ~has_geom
	        warning('NekROM:MissingCoordinates', ...
	            ['Basis snapshots under %s do not include X/Y coordinates. ' ...
	             'Disabling element reordering and loading the grid from the main session snapshots instead.'], bas_prefix);
	        reorder = false;
	    end

	    [pod_u, pod_v, ~, pod_t] = get_snaps(bas_snaps, reorder);
	    if has_geom
	        [x_fom, y_fom] = get_grid(bas_snaps, reorder);
	        inde = bas_snaps.flds{1}.inde;
	    else
	        session_prefix = resolve_case_session_prefix(case_path, casename);
	        if isempty(session_prefix)
	            error('Unable to locate %s0.f* snapshots to load the mesh coordinates.', casename);
	        end
	        session_snaps = NekSnaps(session_prefix);
	        [x_fom, y_fom] = get_grid(session_snaps, false);
	        inde = session_snaps.flds{1}.inde;
	    end
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

	    % Load offline ROM operators.
	    ops = load_full_ops_struct(ops_dir);
	    au_full = ops.au;
	    bu_full = ops.bu;
	    cu_full = ops.cu;
	    u0_full = ops.u0;
		    uk_full = ops.uk;
		    mb = ops.nb;
		    ns = ops.ns;

	    if nb ~= mb
	        error('Configured nb (%d) does not match ops/nb (%d).', nb, mb);
	    end
	    if ~isempty(case_meta.ns) && ns ~= case_meta.ns
	        error('Configured ns (%d) does not match ops/ns (%d).', case_meta.ns, ns);
	    end

		    required_cols = nb + 1;
		    if size(pod_u, 2) < required_cols || size(pod_v, 2) < required_cols
		        error(['POD basis does not contain enough modes for nb=%d.\n' ...
		               'Expected at least %d columns (for 2:nb+1 indexing), but got pod_u=%d, pod_v=%d.\n' ...
		               'Regenerate the offline basis or reduce nb in config.m.'], ...
		               nb, required_cols, size(pod_u, 2), size(pod_v, 2));
		    end
		    if has_thermal_ops
		        if isempty(pod_t) || size(pod_t, 2) < required_cols
		            got_cols = 0;
		            if ~isempty(pod_t)
		                got_cols = size(pod_t, 2);
		            end
		            error(['This case appears to include a temperature equation (ops/{at,bt,ct,t0,tk}).\n' ...
		                   'The loaded POD basis does not include enough temperature modes for nb=%d.\n' ...
		                   'Expected at least %d columns, but got pod_t=%d.\n' ...
		                   'Regenerate the offline basis or reduce nb.'], ...
		                   nb, required_cols, got_cols);
		        end
		    end

	    Me = get_Me(x_fom, y_fom);
	    % Reduced L2 inner product matrix for energy/skew enforcement.
	    % This matches the kinetic energy calculation based on Me on the ROM grid.
	    Me_stack = [Me; Me];
	    phi_energy = [pod_u(:, 2:nb+1); pod_v(:, 2:nb+1)];
	    b_l2 = phi_energy' * bsxfun(@times, Me_stack, phi_energy);

	    %% DEIM Initialization
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
            cba_prefix = resolve_snapshot_prefix(case_path, casename, 'cba');
            if isempty(cba_prefix)
                if exist(fullfile(ops_dir, 'deim_npts'), 'file')
                    error(['Nonlinear DEIM basis snapshots (cba*) were not found under ', case_path, '.\n' ...
                           'This case does have ops/deim_* artifacts; rerun with NEKROM_DEIM_FROM_OPS=1 to use them.\n' ...
                           'Otherwise regenerate snapshots with deim:dumpnls=yes and rerun makerom.']);
                end
                error(['Nonlinear DEIM basis snapshots (cba*) were not found under ', case_path, '.\n' ...
                       'Regenerate snapshots with deim:dumpnls=yes and rerun makerom.']);
            end
            nl_bas_obj_nr = NekSnaps(cba_prefix);
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
                csn_prefix = resolve_snapshot_prefix(case_path, casename, 'csn');
                if isempty(csn_prefix)
                    error(['Nonlinear DEIM snapshot data (csn*) was not found under ', case_path, '.\n' ...
                           'Regenerate snapshots with deim:dumpnls=yes and rerun makerom.']);
                end
                nl_snaps_obj = NekSnaps(csn_prefix);
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
                x_deim, y_deim, ndeim_pts, n_os_points, ps_alg, deim_dealias, deim_dealias_quad, deim_alpha, deim_dealias_cquad);

            if deim_dealias_quad
                warning('NekROM:DEIMQuadNoPersist', ...
                    ['NEKROM_DEIM_DEALIAS_QUAD=1 uses an overintegrated MATLAB-only DEIM path. ' ...
                     'Skipping ops/ persistence because the current Fortran runtime cannot load it.']);
            elseif deim_dealias_cquad
                allow_persist = read_env_bool('NEKROM_DEIM_DEALIAS_CQUAD_PERSIST', false);
                if allow_persist
                    warning('NekROM:DEIMCQuadPersist', ...
                        ['NEKROM_DEIM_DEALIAS_CQUAD=1 is experimental and may destabilize some cases. ' ...
                         'Persisting DEIM artifacts to ops/ because NEKROM_DEIM_DEALIAS_CQUAD_PERSIST=1 was set.']);
                    save_deim_artifacts(ops_dir, deim_data);
                else
                    warning('NekROM:DEIMCQuadNoPersist', ...
                        ['NEKROM_DEIM_DEALIAS_CQUAD=1 is experimental and may destabilize some cases. ' ...
                         'Skipping ops/ persistence. Set NEKROM_DEIM_DEALIAS_CQUAD_PERSIST=1 to override.']);
                end
            else
                save_deim_artifacts(ops_dir, deim_data);
            end
        end
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

	    % Get reduced dimensional operators (velocity).
	    [au, a0, bu, cu, c0, c1, c2, c3, u0, uk, ukmin, ukmax] = get_r_dim_ops(au_full, bu_full, cu_full, u0_full, uk_full, nb);

	    % Optional reduced operators (temperature + buoyancy).
	    has_temp = has_thermal_ops && isfield(ops, 'has_thermal') && ops.has_thermal;
	    use_tdeim = false;
	    tdeim_data = struct();
	    buoy_enabled = false;
	    if has_temp
	        at_full = ops.at;
	        bt_full = ops.bt;
	        ct_full = ops.ct;
	        t0_full = ops.t0;
	        tk_full = ops.tk;

	        at = at_full(2:nb+1, 2:nb+1);
	        a0t = at_full(2:nb+1, 1);
	        bt = bt_full(2:nb+1, 2:nb+1);
	        ct0 = reshape(ct_full(1:nb, 1:nb+1, 1:nb+1), nb*(nb+1), nb+1);

	        t0 = t0_full(1:nb+1);
	        tk = tk_full(1:nb+1, :);
	        tkmin = min(tk, [], 2);
	        tkmax = max(tk, [], 2);

	        disable_tdeim = read_env_bool('MOR_DISABLE_TDEIM', false);
	        use_tdeim_from_ops = read_env_bool('NEKROM_TDEIM_FROM_OPS', true);
	        have_tdeim_ops = exist(fullfile(ops_dir, 'tdeim_npts'), 'file') ~= 0;

	        if have_tdeim_ops && use_tdeim_from_ops && ~disable_tdeim && ismember(conv_approach, {'deim', 'clsdeim', 'mclsdeim'})
	            t_nbnl = case_meta.deim_nbnl;
	            if isempty(t_nbnl) || ~isfinite(t_nbnl) || t_nbnl <= 0
	                error(['TDEIM requires a valid [DEIM] nbnl entry in ' case_meta.case_file '.']);
	            end
	            tdeim_data = load_tdeim_artifacts(ops_dir, nb, t_nbnl);
	            use_tdeim = true;
	        end

	        gx = read_env_scalar('NEKROM_GX', case_meta.buoyancy.gx);
	        gy = read_env_scalar('NEKROM_GY', case_meta.buoyancy.gy);
	        gz = read_env_scalar('NEKROM_GZ', case_meta.buoyancy.gz);

	        if isfield(ops, 'has_buoyancy') && ops.has_buoyancy
	            buoy_enabled = norm([gx, gy, gz], 2) > 0;
	            if ~buoy_enabled
	                warning('NekROM:BuoyancyDisabled', ...
	                    ['Buoyancy operators are present under ops/, but the driver gravity vector is zero.\n' ...
	                     'Set NEKROM_GX/NEKROM_GY/NEKROM_GZ to enable buoyancy forcing in the MATLAB driver.']);
	            end
	        end
	    end

    %% Initialization
    time   = 0;
	    rhs    = zeros(nb, 1);
	    ext    = zeros(nb, 3);
	    hufac  = [];
	    rhs_t  = zeros(nb, 1);
	    ext_t  = zeros(nb, 3);
	    ht_fac = [];
	    nan_detected = false;
	    aborted_step = NaN;

    % Preallocate outputs
    num_outputs = floor(nsteps / iostep);
	    ucoef = zeros(num_outputs, nb+1);
	    tcoef = zeros(num_outputs, nb+1);
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
	    t = zeros(nb+1, 3);
	    if has_temp
	        t(:,1) = t0;
	    end
	    [alphas, betas] = setcoef();

    u_proj = pod_u(:, 1:nb+1) * u(:, 1);
    v_proj = pod_v(:, 1:nb+1) * u(:, 1);

	    field_data = struct('u', u_proj, 'v', v_proj, 'x', x_fom, 'y', y_fom, 'inde', inde, 'size', size(x_fom), 'time', 0.0, 'iostep', 0);
	    if has_temp
	        t_proj = pod_t(:, 1:nb+1) * t(:, 1);
	        field_data.t = t_proj;
	    end
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
	            ht_fac = [];
	        end

	        ext(:,3) = ext(:,2);
	        ext(:,2) = ext(:,1);
	        if has_temp
	            ext_t(:,3) = ext_t(:,2);
	            ext_t(:,2) = ext_t(:,1);
	        end

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

        if enforce_skew_adjoint
            % Minimal skew-adjoint correction: enforce <u, c(u)>_B = 0 by projection.
            % Default uses the reduced L2 inner product (matches kinetic energy);
            % set NEKROM_CONV_SKEW_INNER=ips to use the POD inner product (bu).
            apply = false;
            if strcmp(skew_apply, 'all')
                apply = true;
            elseif strcmp(skew_apply, 'deim')
                apply = strcmp(conv_approach, 'deim');
            elseif strcmp(skew_apply, 'deim_family')
                apply = ismember(conv_approach, {'deim', 'clsdeim', 'mclsdeim'});
            else
                error('Invalid NEKROM_CONV_SKEW_APPLY: "%s" (use "deim", "deim_family", or "all").', skew_apply);
            end

            if apply
                if strcmp(skew_inner, 'ips')
                    b_skew = bu;
                elseif strcmp(skew_inner, 'l2')
                    b_skew = b_l2;
                else
                    error('Invalid NEKROM_CONV_SKEW_INNER: "%s" (use "l2" or "ips").', skew_inner);
                end
                [c_coef, ~] = enforce_conv_energy_zero_work(c_coef, u(2:end, 1), b_skew);
            end
        end

	        ext(:,1) = -c_coef - nu * a0;
	        if has_temp && buoy_enabled
	            buoy_term = zeros(nb, 1);
	            tvec = t(:, 1);
	            if isfield(ops, 'buxt')
	                buxt = ops.buxt;
	                if size(buxt, 1) == nb + 1
	                    buxt = buxt(2:end, :);
	                end
	                buoy_term = buoy_term - gx * (buxt * tvec);
	            end
	            if isfield(ops, 'buyt')
	                buyt = ops.buyt;
	                if size(buyt, 1) == nb + 1
	                    buyt = buyt(2:end, :);
	                end
	                buoy_term = buoy_term - gy * (buyt * tvec);
	            end
	            if isfield(ops, 'buzt')
	                buzt = ops.buzt;
	                if size(buzt, 1) == nb + 1
	                    buzt = buzt(2:end, :);
	                end
	                buoy_term = buoy_term - gz * (buzt * tvec);
	            end
	            ext(:,1) = ext(:,1) + buoy_term;
	        end

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
	                hufac = chol(h);
	            end
	            u_new = [1; (hufac \ (hufac' \ rhs))];
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

	        % Temperature step (optional)
	        if has_temp
	            if use_tdeim
	                ct_coef = conv_tdeim(utmp(:, 1), t(:, 1), tdeim_data, conv_approach);
	            else
	                ct_coef = (reshape(ct0 * utmp(:, 1), nb, nb + 1) * t(:, 1));
	            end

	            ext_t(:, 1) = -ct_coef - kappa * a0t;
	            rhs_t = (ext_t * alphas(:, ito)) - bt * (t(2:end, :) * betas(2:end, ito)) / dt;

	            if isempty(ht_fac)
	                ht = bt * betas(1, ito) / dt + at * kappa;
	                ht_fac = chol(ht);
	            end
	            t_new = [1; (ht_fac \ (ht_fac' \ rhs_t))];
	            t = shift(t, t_new, 3);

	            if any(isnan(t(:, 1)))
	                fprintf('NaN detected in temperature at step %d. Aborting.\n', istep);
	                nan_detected = true;
	                aborted_step = istep;
	                break;
	            end
	        end

        % IO Routine
	        if mod(istep, iostep) == 0
	            % Removed verbose IOSTEP print (progress indicator handles this)
	            ucoef(io_idx, :) = u(:, 1)';
	            if has_temp
	                tcoef(io_idx, :) = t(:, 1)';
	            end

	            u_proj = pod_u(:, 1:nb+1) * u(:, 1);
	            v_proj = pod_v(:, 1:nb+1) * u(:, 1);

            kes(io_idx) = 0.5 * (u_proj' * (Me .* u_proj) + v_proj' * (Me .* v_proj));
            momentums(io_idx, :) = [sum(Me .* u_proj), sum(Me .* v_proj)];

	            field_data.u = u_proj;
	            field_data.v = v_proj;
	            if has_temp
	                t_proj = pod_t(:, 1:nb+1) * t(:, 1);
	                field_data.t = t_proj;
	            end
	            field_data.time = time;
	            field_data.iostep = floor(istep / iostep);

	            output_fields(basepath, field_data, ifvort, ifwrite, ifvis);
            io_idx = io_idx + 1;
        end
    end

	    %% Outputs & Visualization
	    ucoef = ucoef(1:io_idx-1, :);
	    if has_temp
	        tcoef = tcoef(1:io_idx-1, :);
	    else
	        tcoef = [];
	    end
	    kes = kes(1:io_idx-1);
	    momentums = momentums(1:io_idx-1, :);

	    fileID = fopen(fullfile(casedir, 'ucoef'), 'w');
	    fprintf(fileID, '%24.15e\n', ucoef'); 
	    fclose(fileID);
	    if has_temp
	        fileID = fopen(fullfile(casedir, 'tcoef'), 'w');
	        fprintf(fileID, '%24.15e\n', tcoef');
	        fclose(fileID);
	    end

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

	    if isempty(tcoef)
	        max_tcoef_norm = NaN;
	    else
	        max_tcoef_norm = max(sqrt(sum(tcoef(:, 2:end).^2, 2)));
	    end

    completed = ~nan_detected && ((io_idx - 1) == num_outputs);

    results = struct();
    results.casename = casename;
    results.conv_approach = conv_approach;
    results.deim_finegrid = deim_finegrid;
    results.deim_dealias = deim_dealias;
    results.deim_dealias_cquad = deim_dealias_cquad;
    results.deim_dealias_quad = deim_dealias_quad;
    results.deim_dealias_mode = deim_dealias_mode;
    results.enforce_skew_adjoint = enforce_skew_adjoint;
    results.skew_inner = skew_inner;
    results.skew_apply = skew_apply;
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
	    results.max_tcoef_norm = max_tcoef_norm;
	    results.ucoef = ucoef;
	    results.tcoef = tcoef;
	    results.has_temp = has_temp;
	    if has_temp
	        results.kappa = kappa;
	        results.use_tdeim = use_tdeim;
	        results.buoy_enabled = buoy_enabled;
	        results.gx = gx;
	        results.gy = gy;
	        results.gz = gz;
	    else
	        results.kappa = [];
	        results.use_tdeim = false;
	        results.buoy_enabled = false;
	        results.gx = 0.0;
	        results.gy = 0.0;
	        results.gz = 0.0;
	    end
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
