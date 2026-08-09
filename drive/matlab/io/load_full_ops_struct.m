function ops = load_full_ops_struct(ops_dir, varargin)
% LOAD_FULL_OPS_STRUCT Load ROM operators from ops/ into a single struct.
%
% This is a more flexible alternative to load_full_ops.m that supports both
% velocity-only and thermo-fluid operator bundles.
%
% Usage:
%   ops = load_full_ops_struct(ops_dir)
%
% Optional name/value:
%   'LoadThermal' (logical, default: auto) - force/disable thermal loads

    if nargin < 1 || isempty(ops_dir)
        error('load_full_ops_struct:MissingOpsDir', 'ops_dir is required.');
    end

    load_thermal = [];
    if ~isempty(varargin)
        for i = 1:2:numel(varargin)
            key = lower(strtrim(varargin{i}));
            value = varargin{i+1};
            switch key
                case 'loadthermal'
                    load_thermal = logical(value);
                otherwise
                    error('load_full_ops_struct:BadArg', 'Unknown option: %s', key);
            end
        end
    end

    nb = read_real_scalar(fullfile(ops_dir, 'nb'));
    nb = round(nb);
    ns = read_real_scalar(fullfile(ops_dir, 'ns'));
    ns = round(ns);

    ops = struct();
    ops.ops_dir = ops_dir;
    ops.nb = nb;
    ops.ns = ns;

    ops.au = read_real_matrix(fullfile(ops_dir, 'au'), nb + 1, nb + 1);
    ops.bu = read_real_matrix(fullfile(ops_dir, 'bu'), nb + 1, nb + 1);
    ops.cu = read_real_tensor(fullfile(ops_dir, 'cu'), nb, nb + 1, nb + 1);
    ops.u0 = read_real_vector(fullfile(ops_dir, 'u0'), nb + 1);
    ops.uk = read_real_matrix(fullfile(ops_dir, 'uk'), nb + 1, ns);

    thermal_present = exist(fullfile(ops_dir, 'at'), 'file') || ...
                     exist(fullfile(ops_dir, 'bt'), 'file') || ...
                     exist(fullfile(ops_dir, 'ct'), 'file');
    if isempty(load_thermal)
        load_thermal = thermal_present;
    end

    ops.has_thermal = false;
    if load_thermal
        if ~exist(fullfile(ops_dir, 'at'), 'file') || ...
           ~exist(fullfile(ops_dir, 'bt'), 'file') || ...
           ~exist(fullfile(ops_dir, 'ct'), 'file') || ...
           ~exist(fullfile(ops_dir, 't0'), 'file') || ...
           ~exist(fullfile(ops_dir, 'tk'), 'file')
            error('load_full_ops_struct:MissingThermalOps', ...
                'Thermal operators requested but one or more of ops/{at,bt,ct,t0,tk} is missing in %s.', ops_dir);
        end

        ops.at = read_real_matrix(fullfile(ops_dir, 'at'), nb + 1, nb + 1);
        ops.bt = read_real_matrix(fullfile(ops_dir, 'bt'), nb + 1, nb + 1);
        ops.ct = read_real_tensor(fullfile(ops_dir, 'ct'), nb, nb + 1, nb + 1);
        ops.t0 = read_real_vector(fullfile(ops_dir, 't0'), nb + 1);
        ops.tk = read_real_matrix(fullfile(ops_dir, 'tk'), nb + 1, ns);
        ops.has_thermal = true;

        % Optional buoyancy coupling operators.
        ops.has_buoyancy = false;
        if exist(fullfile(ops_dir, 'buxt'), 'file')
            ops.buxt = read_buoyancy_matrix(fullfile(ops_dir, 'buxt'), nb);
            ops.has_buoyancy = true;
        end
        if exist(fullfile(ops_dir, 'buyt'), 'file')
            ops.buyt = read_buoyancy_matrix(fullfile(ops_dir, 'buyt'), nb);
            ops.has_buoyancy = true;
        end
        if exist(fullfile(ops_dir, 'buzt'), 'file')
            ops.buzt = read_buoyancy_matrix(fullfile(ops_dir, 'buzt'), nb);
            ops.has_buoyancy = true;
        end
        if exist(fullfile(ops_dir, 'but'), 'file')
            ops.but = read_buoyancy_matrix(fullfile(ops_dir, 'but'), nb);
            ops.has_buoyancy = true;
        end
    else
        ops.has_buoyancy = false;
    end
end

function value = read_real_scalar(path)
    data = dlmread(path);
    if isempty(data)
        error('load_full_ops_struct:EmptyFile', 'Empty scalar file: %s', path);
    end
    value = data(1);
end

function vec = read_real_vector(path, n)
    data = dlmread(path);
    if numel(data) < n
        error('load_full_ops_struct:ShortVector', ...
            'Expected %d entries in %s, got %d.', n, path, numel(data));
    end
    vec = reshape(data(1:n), [n, 1]);
end

function mat = read_real_matrix(path, m1, m2)
    data = dlmread(path);
    expected = m1 * m2;
    if numel(data) < expected
        error('load_full_ops_struct:ShortMatrix', ...
            'Expected %d entries in %s (%dx%d), got %d.', expected, path, m1, m2, numel(data));
    end
    mat = reshape(data(1:expected), [m1, m2]);
end

function ten = read_real_tensor(path, n1, n2, n3)
    data = dlmread(path);
    expected = n1 * n2 * n3;
    if numel(data) < expected
        error('load_full_ops_struct:ShortTensor', ...
            'Expected %d entries in %s (%dx%dx%d), got %d.', expected, path, n1, n2, n3, numel(data));
    end
    ten = reshape(data(1:expected), [n1, n2, n3]);
end

function mat = read_buoyancy_matrix(path, nb)
    % Buoyancy coupling operators may be stored as:
    % - (nb+1)x(nb+1) (full-space)
    % - nb x (nb+1) (dynamic velocity modes only, consistent with cu's first dimension)
    data = dlmread(path);
    count = numel(data);
    full_count = (nb + 1) * (nb + 1);
    dyn_count = nb * (nb + 1);

    if count < min(full_count, dyn_count)
        error('load_full_ops_struct:ShortBuoyancyMatrix', ...
            'Expected %d or %d entries in %s, got %d.', full_count, dyn_count, path, count);
    end

    if count >= full_count
        mat = reshape(data(1:full_count), [nb + 1, nb + 1]);
        return;
    end

    if count >= dyn_count
        mat = reshape(data(1:dyn_count), [nb, nb + 1]);
        return;
    end

    error('load_full_ops_struct:BadBuoyancyMatrix', ...
        'Unable to parse buoyancy operator %s with nb=%d (got %d entries).', path, nb, count);
end
