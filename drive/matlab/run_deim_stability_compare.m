function comparison = run_deim_stability_compare(varargin)
    % RUN_DEIM_STABILITY_COMPARE Compare coarse-grid and fine-grid DEIM stability.
    %
    % The comparison uses the currently configured case and DEIM-family
    % method, then runs the ROM twice:
    %   1. deim_finegrid = false
    %   2. deim_finegrid = true
    %
    % The driver returns a results struct for each run, and this wrapper
    % condenses the output into a short stability summary.

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

    nsteps_override = [];
    iostep_override = [];
    if nargin >= 1
        nsteps_override = varargin{1};
    end
    if nargin >= 2
        iostep_override = varargin{2};
    end

    config;

    common_env = struct( ...
        'NEKROM_CASE', casename, ...
        'NEKROM_CONV_APPROACH', conv_approach, ...
        'NEKROM_IF_RUN_TESTS', '0', ...
        'NEKROM_IFWRITE', '0', ...
        'NEKROM_IFVIS', '0', ...
        'NEKROM_IFVORT', '0', ...
        'NEKROM_IFPLOT', '0');
    if deim_dealias_quad
        common_env.NEKROM_DEIM_DEALIAS_QUAD = '1';
    else
        common_env.NEKROM_DEIM_DEALIAS = '1';
    end

    if ~isempty(nsteps_override)
        common_env.NEKROM_NSTEPS = num2str(nsteps_override);
    end
    if ~isempty(iostep_override)
        common_env.NEKROM_IOSTEP = num2str(iostep_override);
    end

    coarse = run_driver_with_env(common_env, 'NEKROM_DEIM_FINEGRID', '0');
    fine   = run_driver_with_env(common_env, 'NEKROM_DEIM_FINEGRID', '1');

    comparison = struct();
    comparison.case = casename;
    comparison.conv_approach = conv_approach;
    comparison.deim_dealias = deim_dealias;
    comparison.deim_dealias_quad = deim_dealias_quad;
    comparison.deim_dealias_mode = deim_dealias_mode;
    comparison.coarse = summarize_run(coarse);
    comparison.fine = summarize_run(fine);

    fprintf('\nDEIM stability comparison for %s (%s, mode=%s)\n', ...
        casename, conv_approach, comparison.deim_dealias_mode);
    print_summary('coarse', comparison.coarse);
    print_summary('fine', comparison.fine);
    fprintf('Comparison complete.\n');
end

function results = run_driver_with_env(env_map, fine_key, fine_value)
    keys = fieldnames(env_map);
    old_values = struct();

    for i = 1:numel(keys)
        key = keys{i};
        old_values.(key) = getenv(key);
        setenv(key, env_map.(key));
    end

    old_fine = getenv(fine_key);
    setenv(fine_key, fine_value);

    cleanup = onCleanup(@() restore_env(old_values, fine_key, old_fine));
    results = driver();
    clear cleanup;
end

function restore_env(old_values, fine_key, old_fine)
    keys = fieldnames(old_values);
    for i = 1:numel(keys)
        key = keys{i};
        if isempty(old_values.(key))
            setenv(key, '');
        else
            setenv(key, old_values.(key));
        end
    end

    if isempty(old_fine)
        setenv(fine_key, '');
    else
        setenv(fine_key, old_fine);
    end
end

function summary = summarize_run(results)
    summary.completed = results.completed;
    summary.nan_detected = results.nan_detected;
    summary.aborted_step = results.aborted_step;
    summary.wall_time_sec = results.wall_time_sec;
    summary.ke_initial = results.ke_initial;
    summary.ke_final = results.ke_final;
    summary.ke_rel_drift = results.ke_rel_drift;
    summary.ke_peak_ratio = results.ke_peak_ratio;
    summary.max_momentum_norm = results.max_momentum_norm;
    summary.final_momentum_norm = results.final_momentum_norm;
    summary.max_ucoef_norm = results.max_ucoef_norm;
end

function print_summary(label, summary)
    fprintf('%s: completed=%d, nan=%d, aborted_step=%g, wall=%.2fs, ke_drift=% .3e, ke_peak=% .3e, max_ucoef=% .3e\n', ...
        label, summary.completed, summary.nan_detected, summary.aborted_step, ...
        summary.wall_time_sec, summary.ke_rel_drift, summary.ke_peak_ratio, summary.max_ucoef_norm);
end
