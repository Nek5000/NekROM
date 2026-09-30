function comparison = run_nus_compare(varargin)
% RUN_NUS_COMPARE Compare MATLAB and Fortran Nusselt histories for RB-style cases.
%
% Usage:
%   comparison = run_nus_compare()
%   comparison = run_nus_compare(matlab_output_root, fortran_case_root)
%   comparison = run_nus_compare(case_prefix, matlab_output_root, fortran_case_root)
%
% The helper looks for a qoi/nus weight file in the Fortran case tree, then
% projects the MATLAB temperature coefficients onto that QoI and compares the
% resulting scalar history against every nus.dat file under the Fortran tree.
% This is meant for sweep-style examples such as rb_axi, where field snapshots
% are not the right unit of comparison.

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

    case_prefix = 'rb';
    matlab_output_root = fullfile(driver_dir, 'output');
    fortran_case_root = fullfile(driver_dir, '..', '..', 'examples', 'rb_axi');

    if nargin >= 1 && ~isempty(varargin{1})
        case_prefix = strtrim(varargin{1});
    end
    if nargin >= 2 && ~isempty(varargin{2})
        matlab_output_root = varargin{2};
    end
    if nargin >= 3 && ~isempty(varargin{3})
        fortran_case_root = varargin{3};
    end

    matlab_case_dir = find_latest_matlab_output(matlab_output_root, case_prefix);
    qoi_file = find_qoi_file(fortran_case_root);
    qoi_weights = load_qoi_weights(qoi_file);
    if isempty(qoi_weights)
        error('run_nus_compare:MissingQoI', ...
            'Could not load qoi/nus weights from %s.', qoi_file);
    end

    mat_nus = load_matlab_nus(matlab_case_dir, qoi_weights);
    if isempty(mat_nus)
        error('run_nus_compare:MissingMatlabHistory', ...
            'Could not load a MATLAB tcoef history from %s.', matlab_case_dir);
    end

    candidates = find_fortran_nus_candidates(fortran_case_root);
    if isempty(candidates)
        error('run_nus_compare:MissingFortranCandidates', ...
            'No nus.dat files found under %s.', fortran_case_root);
    end

    comparison = struct();
    comparison.case_prefix = case_prefix;
    comparison.matlab_output_root = matlab_output_root;
    comparison.fortran_case_root = fortran_case_root;
    comparison.matlab_case_dir = matlab_case_dir;
    comparison.qoi_file = qoi_file;
    comparison.candidates = cell(numel(candidates), 1);
    comparison.best = struct();

    fprintf('\nMATLAB vs Fortran Nusselt comparison\n');
    fprintf('MATLAB output root: %s\n', matlab_output_root);
    fprintf('Fortran root:       %s\n', fortran_case_root);
    fprintf('QoI weights:        %s\n\n', qoi_file);

    best_err = inf;
    best_idx = 0;
    for i = 1:numel(candidates)
        summary = compare_candidate(mat_nus, candidates(i));
        comparison.candidates{i} = summary;
        print_candidate_summary(summary);

        if summary.rel_err < best_err
            best_err = summary.rel_err;
            best_idx = i;
        end
    end

    if best_idx > 0
        comparison.best = comparison.candidates{best_idx};
        fprintf('\nBest match: %s (rel_err=% .3e)\n', ...
            comparison.best.name, comparison.best.rel_err);
    end
end

function case_dir = find_latest_matlab_output(output_root, case_prefix)
    case_dir = '';
    if ~exist(output_root, 'dir')
        return;
    end

    matches = dir(fullfile(output_root, [case_prefix '_*']));
    matches = matches([matches.isdir]);
    if isempty(matches)
        return;
    end

    [~, idx] = max([matches.datenum]);
    case_dir = fullfile(output_root, matches(idx).name);
end

function qoi_file = find_qoi_file(fortran_case_root)
    candidates = { ...
        fullfile(fortran_case_root, 'qoi', 'nus'), ...
        fullfile(fortran_case_root, 'qoi', 'nus.dat'), ...
        fullfile(fortran_case_root, 'nus.dat') ...
    };

    qoi_file = '';
    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file')
            qoi_file = candidates{i};
            return;
        end
    end
end

function qoi_weights = load_qoi_weights(qoi_file)
    qoi_weights = [];
    if isempty(qoi_file) || ~exist(qoi_file, 'file')
        return;
    end

    data = dlmread(qoi_file);
    qoi_weights = data(:);
end

function mat_nus = load_matlab_nus(matlab_case_dir, qoi_weights)
    mat_nus = [];

    if isempty(matlab_case_dir) || ~exist(matlab_case_dir, 'dir')
        return;
    end

    stability_file = fullfile(matlab_case_dir, 'stability.mat');
    if exist(stability_file, 'file')
        data = load(stability_file, 'results');
        if isfield(data, 'results') && isfield(data.results, 'tcoef') && ~isempty(data.results.tcoef)
            tcoef = data.results.tcoef;
            if size(tcoef, 2) ~= numel(qoi_weights)
                return;
            end
            mat_nus = tcoef * qoi_weights(:);
            return;
        end
    end

    tcoef_file = fullfile(matlab_case_dir, 'tcoef');
    if ~exist(tcoef_file, 'file')
        return;
    end

    raw = dlmread(tcoef_file);
    if isempty(raw) || mod(numel(raw), numel(qoi_weights)) ~= 0
        return;
    end

    mat_nus = reshape(raw, numel(qoi_weights), []).';
    mat_nus = mat_nus * qoi_weights(:);
end

function candidates = find_fortran_nus_candidates(fortran_case_root)
    candidates = struct('name', {}, 'path', {}, 'values', {}, 'step', {}, 'time', {});

    if ~exist(fortran_case_root, 'dir')
        return;
    end

    direct_nus = fullfile(fortran_case_root, 'nus.dat');
    if exist(direct_nus, 'file')
        candidates(end+1) = load_nus_candidate(fortran_case_root, direct_nus); %#ok<AGROW>
    end

    subdirs = dir(fullfile(fortran_case_root, '*'));
    rom_dirs = pick_candidate_dirs(subdirs, 'rom');
    if isempty(rom_dirs)
        rom_dirs = pick_candidate_dirs(subdirs, '');
    end

    for i = 1:numel(rom_dirs)
        nus_file = fullfile(fortran_case_root, rom_dirs{i}, 'nus.dat');
        if exist(nus_file, 'file')
            candidates(end+1) = load_nus_candidate(fullfile(fortran_case_root, rom_dirs{i}), nus_file); %#ok<AGROW>
        end
    end
end

function names = pick_candidate_dirs(subdirs, prefix)
    names = {};
    for i = 1:numel(subdirs)
        if ~subdirs(i).isdir
            continue;
        end
        name = subdirs(i).name;
        if strcmp(name, '.') || strcmp(name, '..')
            continue;
        end
        if isempty(prefix) || strncmp(name, prefix, numel(prefix))
            names{end+1} = name; %#ok<AGROW>
        end
    end
end

function candidate = load_nus_candidate(dir_path, nus_file)
    fid = fopen(nus_file, 'r');
    if fid < 0
        error('run_nus_compare:UnreadableNus', 'Unable to open %s.', nus_file);
    end
    cleanup = onCleanup(@() fclose(fid));

    parsed = textscan(fid, '%f %f %f %s', 'CollectOutput', true, 'MultipleDelimsAsOne', true);
    data = parsed{1};
    clear cleanup;

    candidate = struct();
    candidate.name = dir_path;
    candidate.path = nus_file;
    candidate.values = data;
    candidate.step = [];
    candidate.time = [];
    if size(data, 2) >= 3
        candidate.step = data(:, 1);
        candidate.time = data(:, 2);
    end
end

function summary = compare_candidate(mat_nus, candidate)
    summary = struct();
    summary.name = candidate.name;
    summary.path = candidate.path;
    summary.step = candidate.step;
    summary.time = candidate.time;

    values = candidate.values;
    if size(values, 2) >= 3
        values = values(:, 3);
        if numel(values) > 1
            values = values(2:end);
        end
    else
        values = values(:);
    end

    n = min(numel(mat_nus), numel(values));
    summary.common_count = n;
    if n == 0
        summary.rel_err = Inf;
        summary.final_err = Inf;
        summary.first_value = NaN;
        summary.first_matlab = NaN;
        return;
    end

    summary.rel_err = relative_error(mat_nus(1:n), values(1:n));
    summary.final_err = relative_error(mat_nus(n), values(n));
    summary.first_value = values(1);
    summary.first_matlab = mat_nus(1);
end

function print_candidate_summary(summary)
    fprintf('%-28s  count %4d  rel % .3e  final % .3e  first % .6f vs % .6f\n', ...
        summary.name, summary.common_count, summary.rel_err, summary.final_err, ...
        summary.first_value, summary.first_matlab);
end

function err = relative_error(a, b)
    denom = max(norm(b(:)), eps);
    err = norm(a(:) - b(:)) / denom;
end
