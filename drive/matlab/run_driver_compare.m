function comparison = run_driver_compare(varargin)
% RUN_DRIVER_COMPARE Compare MATLAB and Fortran driver outputs case-by-case.
%
% Usage:
%   comparison = run_driver_compare()
%   comparison = run_driver_compare(cases)
%   comparison = run_driver_compare(cases, matlab_output_root, fortran_root)
%
% The harness compares the latest MATLAB output directory for each case
% against the matching Fortran driver output directory and reports snapshot,
% kinetic-energy, and momentum agreement.

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

    cases = {'ldc', 'cyl', 'shear', 't2d'};
    matlab_output_root = fullfile(driver_dir, 'output');
    fortran_root = fullfile(driver_dir, '..', '..', 'examples');

    if nargin >= 1 && ~isempty(varargin{1})
        cases = normalize_cases(varargin{1});
    end
    if nargin >= 2 && ~isempty(varargin{2})
        matlab_output_root = varargin{2};
    end
    if nargin >= 3 && ~isempty(varargin{3})
        fortran_root = varargin{3};
    end

    comparison = struct();
    comparison.generated_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    comparison.matlab_output_root = matlab_output_root;
    comparison.fortran_root = fortran_root;
    comparison.cases = cell(numel(cases), 1);

    fprintf('\nMATLAB vs Fortran driver comparison\n');
    fprintf('MATLAB output root: %s\n', matlab_output_root);
    fprintf('Fortran root:       %s\n\n', fortran_root);

    for i = 1:numel(cases)
        casename = cases{i};
        matlab_case_dir = find_latest_matlab_output(matlab_output_root, casename);
        fortran_case_dir = find_fortran_case_dir(fortran_root, casename);

        summary = compare_case_outputs(casename, matlab_case_dir, fortran_case_dir);
        comparison.cases{i} = summary;
        print_case_summary(summary);
    end

    fprintf('\nComparison complete.\n');
end

function cases = normalize_cases(input_cases)
    if ischar(input_cases)
        cases = {strtrim(input_cases)};
        return;
    end

    if iscell(input_cases)
        cases = input_cases(:);
        for i = 1:numel(cases)
            cases{i} = strtrim(cases{i});
        end
        return;
    end

    error('run_driver_compare:InvalidCases', ...
        'cases must be a char vector or a cell array of case names.');
end

function case_dir = find_latest_matlab_output(output_root, casename)
    case_dir = '';

    if ~exist(output_root, 'dir')
        return;
    end

    matches = dir(fullfile(output_root, [casename '_*']));
    matches = matches([matches.isdir]);
    if isempty(matches)
        return;
    end

    [~, idx] = max([matches.datenum]);
    case_dir = fullfile(output_root, matches(idx).name);
end

function case_dir = find_fortran_case_dir(fortran_root, casename)
    snaps_dir = fullfile(fortran_root, casename, 'snaps');
    if exist(snaps_dir, 'dir')
        case_dir = snaps_dir;
        return;
    end

    direct_case_dir = fullfile(fortran_root, casename);
    if exist(direct_case_dir, 'dir')
        case_dir = direct_case_dir;
        return;
    end

    case_dir = '';
end

function summary = compare_case_outputs(casename, matlab_case_dir, fortran_case_dir)
    summary = struct();
    summary.case = casename;
    summary.matlab_case_dir = matlab_case_dir;
    summary.fortran_case_dir = fortran_case_dir;
    summary.status = 'skip';
    summary.notes = '';
    summary.matlab_snapshot_count = 0;
    summary.fortran_snapshot_count = 0;
    summary.common_snapshot_count = 0;
    summary.grid_rel_err = NaN;
    summary.field_rel_err_fro = NaN;
    summary.field_rel_err_final = NaN;
    summary.ke_rel_err_max = NaN;
    summary.ke_rel_err_final = NaN;
    summary.momentum_rel_err_max = NaN;
    summary.momentum_rel_err_final = NaN;
    summary.matlab_completed = NaN;
    summary.matlab_nan_detected = NaN;
    summary.matlab_aborted_step = NaN;

    if isempty(matlab_case_dir) || ~exist(matlab_case_dir, 'dir')
        summary.notes = 'Missing MATLAB output directory.';
        return;
    end

    if ~exist(fortran_case_dir, 'dir')
        summary.notes = 'Missing Fortran snaps directory.';
        return;
    end

    matlab_fields_dir = fullfile(matlab_case_dir, 'fields');
    if ~exist(matlab_fields_dir, 'dir')
        summary.notes = 'Missing MATLAB fields directory.';
        return;
    end

    if ~has_primary_snapshot_stream(matlab_fields_dir, casename)
        summary.notes = 'Missing MATLAB primary case0.f* snapshots.';
        return;
    end

    if ~has_primary_snapshot_stream(fortran_case_dir, casename)
        summary.notes = 'Missing Fortran primary case0.f* snapshots.';
        return;
    end

    matlab_snaps = load_snaps_from_dir(matlab_fields_dir, casename);
    fortran_snaps = load_snaps_from_dir(fortran_case_dir, casename);

    [u_mat, v_mat] = get_snaps(matlab_snaps, false);
    [u_frt, v_frt] = get_snaps(fortran_snaps, false);
    [x_mat, y_mat] = get_grid(matlab_snaps, false);
    [x_frt, y_frt] = get_grid(fortran_snaps, false);

    summary.matlab_snapshot_count = size(u_mat, 2);
    summary.fortran_snapshot_count = size(u_frt, 2);
    summary.common_snapshot_count = min(summary.matlab_snapshot_count, summary.fortran_snapshot_count);

    if summary.common_snapshot_count == 0
        summary.notes = 'No overlapping snapshots found.';
        return;
    end

    u_mat = u_mat(:, 1:summary.common_snapshot_count);
    v_mat = v_mat(:, 1:summary.common_snapshot_count);
    u_frt = u_frt(:, 1:summary.common_snapshot_count);
    v_frt = v_frt(:, 1:summary.common_snapshot_count);

    if size(u_mat, 1) ~= size(u_frt, 1) || size(v_mat, 1) ~= size(v_frt, 1)
        summary.status = 'warn';
        summary.notes = sprintf( ...
            'Snapshot size mismatch (%dx%d vs %dx%d).', ...
            size(u_mat, 1), size(u_mat, 2), size(u_frt, 1), size(u_frt, 2));
        return;
    end

    if ~isequal(size(x_mat), size(x_frt)) || ~isequal(size(y_mat), size(y_frt))
        summary.status = 'warn';
        summary.notes = sprintf( ...
            'Grid shape mismatch (%s vs %s).', ...
            mat2str(size(x_mat)), mat2str(size(x_frt)));
        return;
    end

    summary.grid_rel_err = max( ...
        relative_error(x_frt, x_mat), ...
        relative_error(y_frt, y_mat));

    uv_mat = [u_mat; v_mat];
    uv_frt = [u_frt; v_frt];

    summary.field_rel_err_fro = relative_error(uv_frt, uv_mat);
    summary.field_rel_err_final = relative_error(uv_frt(:, end), uv_mat(:, end));

    Me_mat = get_Me(x_mat, y_mat);
    Me_frt = get_Me(x_frt, y_frt);

    [ke_mat, mom_mat] = derive_metrics(u_mat, v_mat, Me_mat);
    [ke_frt, mom_frt] = derive_metrics(u_frt, v_frt, Me_frt);

    summary.ke_rel_err_max = relative_error(ke_frt, ke_mat);
    summary.ke_rel_err_final = relative_error(ke_frt(end), ke_mat(end));
    summary.momentum_rel_err_max = relative_error(mom_frt, mom_mat);
    summary.momentum_rel_err_final = relative_error(mom_frt(end, :), mom_mat(end, :));

    stability_file = fullfile(matlab_case_dir, 'stability.mat');
    if exist(stability_file, 'file')
        data = load(stability_file, 'results');
        if isfield(data, 'results')
            summary.matlab_completed = logical(data.results.completed);
            summary.matlab_nan_detected = logical(data.results.nan_detected);
            summary.matlab_aborted_step = data.results.aborted_step;
        end
    end

    if summary.matlab_snapshot_count == summary.fortran_snapshot_count && ...
            summary.grid_rel_err < 1e-12 && ...
            summary.field_rel_err_final < 1e-8 && ...
            summary.ke_rel_err_final < 1e-8 && ...
            summary.momentum_rel_err_final < 1e-8
        summary.status = 'pass';
    else
        summary.status = 'warn';
    end

    if summary.matlab_snapshot_count ~= summary.fortran_snapshot_count
        summary.notes = sprintf('Snapshot count mismatch (%d vs %d).', ...
            summary.matlab_snapshot_count, summary.fortran_snapshot_count);
    elseif summary.grid_rel_err >= 1e-12
        summary.notes = 'Grid mismatch detected.';
    end
end

function [kes, momentums] = derive_metrics(u, v, Me)
    kes = 0.5 * (sum(u .* (Me .* u), 1) + sum(v .* (Me .* v), 1));
    momentums = [sum(Me .* u, 1).', sum(Me .* v, 1).'];
end

function snaps = load_snaps_from_dir(snapshot_dir, casename)
    if ~exist(snapshot_dir, 'dir')
        error('run_driver_compare:MissingSnapshotDir', ...
            'Snapshot directory not found: %s', snapshot_dir);
    end

    original_dir = pwd;
    cleanup = onCleanup(@() cd(original_dir));
    cd(snapshot_dir);
    snaps = NekSnaps(casename);
    clear cleanup;
end

function err = relative_error(a, b)
    denom = max(norm(b(:)), eps);
    err = norm(a(:) - b(:)) / denom;
end

function present = has_primary_snapshot_stream(snapshot_dir, casename)
    present = ~isempty(dir(fullfile(snapshot_dir, [casename '0.f*'])));
end

function print_case_summary(summary)
    fprintf('%-6s  %-4s  snaps %4d/%-4d  grid % .3e  field % .3e  KE % .3e  mom % .3e', ...
        summary.case, upper(summary.status), ...
        summary.matlab_snapshot_count, summary.fortran_snapshot_count, ...
        summary.grid_rel_err, summary.field_rel_err_final, ...
        summary.ke_rel_err_final, summary.momentum_rel_err_final);

    if ~isempty(summary.notes)
        fprintf('  %s', summary.notes);
    end

    fprintf('\n');
end
