function case_meta = load_case_metadata(case_path, casename)
% LOAD_CASE_METADATA Read stable ROM metadata from a case directory.
%
% The loader prefers authoritative offline metadata under ops/ when it is
% available, then falls back to the case .mor file, and finally to the
% built-in defaults used by the MATLAB driver. If the expected
% <casename>.mor file is absent, it will use a single .mor file in the case
% directory when that is unambiguous.

    if nargin < 1 || isempty(case_path)
        error('load_case_metadata:MissingCasePath', 'case_path is required.');
    end

    if nargin < 2 || isempty(casename)
        [~, casename] = fileparts(regexprep(case_path, '[\\/]+$', ''));
        if isempty(casename)
            casename = 'case';
        end
    end

    case_meta = struct();
    case_meta.case_path = case_path;
    case_meta.case_file = resolve_case_file(case_path, casename);
    case_meta.ops_dir = fullfile(case_path, 'ops');

    case_meta.nb = [];
    case_meta.ns = [];
    case_meta.ips = 'H10';
    case_meta.pod_mode0 = 'avg';
    case_meta.subtract_mean = true;
    case_meta.conv_approach = '';
    case_meta.deim_nbnl = [];
    case_meta.deim_alpha = [];
    case_meta.deim_dumpnls = [];
    case_meta.deim_dumpfine = [];
    case_meta.field = '';
    case_meta.buoyancy = struct('gx', [], 'gy', [], 'gz', [], 'magnitude', [], 'angle', []);

    mor_meta = parse_mor_file(case_meta.case_file);

    if isfield(mor_meta, 'nb') && ~isempty(mor_meta.nb)
        case_meta.nb = mor_meta.nb;
    end
    if isfield(mor_meta, 'ns') && ~isempty(mor_meta.ns)
        case_meta.ns = mor_meta.ns;
    end
    if isfield(mor_meta, 'ips') && ~isempty(mor_meta.ips)
        case_meta.ips = mor_meta.ips;
    end
    if isfield(mor_meta, 'pod_mode0') && ~isempty(mor_meta.pod_mode0)
        case_meta.pod_mode0 = mor_meta.pod_mode0;
    end
    if isfield(mor_meta, 'conv_approach') && ~isempty(mor_meta.conv_approach)
        case_meta.conv_approach = mor_meta.conv_approach;
    end
    if isfield(mor_meta, 'deim_nbnl') && ~isempty(mor_meta.deim_nbnl)
        case_meta.deim_nbnl = mor_meta.deim_nbnl;
    end
    if isfield(mor_meta, 'deim_alpha') && ~isempty(mor_meta.deim_alpha)
        case_meta.deim_alpha = mor_meta.deim_alpha;
    end
    if isfield(mor_meta, 'deim_dumpnls') && ~isempty(mor_meta.deim_dumpnls)
        case_meta.deim_dumpnls = mor_meta.deim_dumpnls;
    end
    if isfield(mor_meta, 'deim_dumpfine') && ~isempty(mor_meta.deim_dumpfine)
        case_meta.deim_dumpfine = mor_meta.deim_dumpfine;
    end
    if isfield(mor_meta, 'field') && ~isempty(mor_meta.field)
        case_meta.field = mor_meta.field;
    end
    if isfield(mor_meta, 'buoyancy') && ~isempty(mor_meta.buoyancy)
        case_meta.buoyancy = merge_buoyancy(case_meta.buoyancy, mor_meta.buoyancy);
    end

    ops_nb_file = fullfile(case_meta.ops_dir, 'nb');
    if exist(ops_nb_file, 'file')
        ops_nb = read_optional_scalar_file(ops_nb_file);
        if ~isempty(ops_nb) && isfinite(ops_nb)
            case_meta.nb = ops_nb;
        end
    end

    ops_ns_file = fullfile(case_meta.ops_dir, 'ns');
    if exist(ops_ns_file, 'file')
        ops_ns = read_optional_scalar_file(ops_ns_file);
        if ~isempty(ops_ns) && isfinite(ops_ns)
            case_meta.ns = ops_ns;
        end
    end

    ops_ips_file = fullfile(case_meta.ops_dir, 'ips');
    if exist(ops_ips_file, 'file')
        case_meta.ips = read_ops_ips(case_meta.ops_dir);
    end

    case_meta.pod_mode0 = lower(strtrim(case_meta.pod_mode0));
    if isempty(case_meta.pod_mode0)
        case_meta.pod_mode0 = 'avg';
    end
    case_meta.subtract_mean = strcmp(case_meta.pod_mode0, 'avg');

    if isempty(case_meta.field)
        case_meta.field = '';
    else
        case_meta.field = lower(strtrim(case_meta.field));
    end

    % Finalize buoyancy vector if magnitude/angle were provided.
    case_meta.buoyancy = finalize_buoyancy(case_meta.buoyancy);

    if isempty(case_meta.conv_approach)
        case_meta.conv_approach = '';
    else
        case_meta.conv_approach = lower(strtrim(case_meta.conv_approach));
    end
end

function mor_file = resolve_case_file(case_path, casename)
    mor_file = fullfile(case_path, [casename '.mor']);
    if exist(mor_file, 'file')
        return;
    end

    mor_files = dir(fullfile(case_path, '*.mor'));
    if isempty(mor_files)
        return;
    end

    case_dir = regexprep(case_path, '[\\/]+$', '');
    [~, case_dir_name] = fileparts(case_dir);
    preferred_names = {[casename '.mor']};
    if ~isempty(case_dir_name)
        preferred_names{end+1} = [case_dir_name '.mor'];
    end

    mor_names = {mor_files.name};
    for i = 1:numel(preferred_names)
        idx = find(strcmpi(mor_names, preferred_names{i}), 1);
        if ~isempty(idx)
            mor_file = fullfile(case_path, mor_files(idx).name);
            return;
        end
    end

    if numel(mor_files) == 1
        mor_file = fullfile(case_path, mor_files(1).name);
        return;
    end

    warning('load_case_metadata:AmbiguousMorFile', ...
        'Using %s because %s was not found in %s.', ...
        mor_files(1).name, [casename '.mor'], case_path);
    mor_file = fullfile(case_path, mor_files(1).name);
end

function meta = parse_mor_file(mor_file)
    meta = struct();

    if ~exist(mor_file, 'file')
        return;
    end

    fid = fopen(mor_file, 'r');
    if fid < 0
        return;
    end

    cleanup = onCleanup(@() fclose(fid));
    current_section = '';

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end

        line = strip_inline_comment(strtrim(line));
        if isempty(line)
            continue;
        end

        if line(1) == '[' && line(end) == ']'
            current_section = lower(strtrim(line(2:end-1)));
            continue;
        end

        tokens = regexp(line, '^([^=]+?)\s*=\s*(.+)$', 'tokens', 'once');
        if isempty(tokens) || isempty(current_section)
            continue;
        end

        key = lower(strtrim(tokens{1}));
        value = strtrim(tokens{2});
        full_key = [current_section ':' key];

        switch full_key
            case 'general:field'
                meta.field = lower(strtrim(value));
            case 'general:nb'
                parsed_value = parse_numeric_value(value);
                if ~isempty(parsed_value)
                    meta.nb = parsed_value;
                end
            case 'general:ns'
                parsed_value = parse_numeric_value(value);
                if ~isempty(parsed_value)
                    meta.ns = parsed_value;
                end
            case 'pod:type'
                parsed_value = normalize_inner_product(value);
                if ~isempty(parsed_value)
                    meta.ips = parsed_value;
                end
            case 'pod:mode0'
                parsed_value = normalize_mode0(value);
                if ~isempty(parsed_value)
                    meta.pod_mode0 = parsed_value;
                end
            case 'deim:mode'
                parsed_value = normalize_deim_mode(value);
                if ~isempty(parsed_value)
                    meta.conv_approach = parsed_value;
                end
            case 'deim:nbnl'
                parsed_value = parse_numeric_value(value);
                if ~isempty(parsed_value)
                    meta.deim_nbnl = parsed_value;
                end
            case 'deim:alpha'
                parsed_value = parse_numeric_value(value);
                if ~isempty(parsed_value)
                    meta.deim_alpha = parsed_value;
                end
            case 'deim:dumpnls'
                parsed_value = parse_boolean_value(value);
                if ~isempty(parsed_value)
                    meta.deim_dumpnls = parsed_value;
                end
            case 'deim:dumpfine'
                parsed_value = parse_boolean_value(value);
                if ~isempty(parsed_value)
                    meta.deim_dumpfine = parsed_value;
                end
            case 'buoyancy:magnitude'
                if ~isfield(meta, 'buoyancy') || isempty(meta.buoyancy)
                    meta.buoyancy = struct('gx', [], 'gy', [], 'gz', [], 'magnitude', [], 'angle', []);
                end
                parsed_value = parse_numeric_value(value);
                if ~isempty(parsed_value)
                    meta.buoyancy.magnitude = parsed_value;
                end
            case 'buoyancy:angle'
                if ~isfield(meta, 'buoyancy') || isempty(meta.buoyancy)
                    meta.buoyancy = struct('gx', [], 'gy', [], 'gz', [], 'magnitude', [], 'angle', []);
                end
                parsed_value = parse_numeric_value(value);
                if ~isempty(parsed_value)
                    meta.buoyancy.angle = parsed_value;
                end
            case 'buoyancy:gx'
                if ~isfield(meta, 'buoyancy') || isempty(meta.buoyancy)
                    meta.buoyancy = struct('gx', [], 'gy', [], 'gz', [], 'magnitude', [], 'angle', []);
                end
                parsed_value = parse_numeric_value(value);
                if ~isempty(parsed_value)
                    meta.buoyancy.gx = parsed_value;
                end
            case 'buoyancy:gy'
                if ~isfield(meta, 'buoyancy') || isempty(meta.buoyancy)
                    meta.buoyancy = struct('gx', [], 'gy', [], 'gz', [], 'magnitude', [], 'angle', []);
                end
                parsed_value = parse_numeric_value(value);
                if ~isempty(parsed_value)
                    meta.buoyancy.gy = parsed_value;
                end
            case 'buoyancy:gz'
                if ~isfield(meta, 'buoyancy') || isempty(meta.buoyancy)
                    meta.buoyancy = struct('gx', [], 'gy', [], 'gz', [], 'magnitude', [], 'angle', []);
                end
                parsed_value = parse_numeric_value(value);
                if ~isempty(parsed_value)
                    meta.buoyancy.gz = parsed_value;
                end
        end
    end
end

function out = merge_buoyancy(a, b)
    out = a;
    fields = {'gx', 'gy', 'gz', 'magnitude', 'angle'};
    for i = 1:numel(fields)
        f = fields{i};
        if isfield(b, f) && ~isempty(b.(f))
            out.(f) = b.(f);
        end
    end
end

function buoy = finalize_buoyancy(buoy)
    if ~isfield(buoy, 'gx'), buoy.gx = []; end
    if ~isfield(buoy, 'gy'), buoy.gy = []; end
    if ~isfield(buoy, 'gz'), buoy.gz = []; end
    if ~isfield(buoy, 'magnitude'), buoy.magnitude = []; end
    if ~isfield(buoy, 'angle'), buoy.angle = []; end

    if isempty(buoy.gx) && isempty(buoy.gy) && ~isempty(buoy.magnitude) && ~isempty(buoy.angle)
        % 2D convenience: angle is degrees displacement from x-axis.
        buoy.gx = buoy.magnitude * cosd(buoy.angle);
        buoy.gy = buoy.magnitude * sind(buoy.angle);
        if isempty(buoy.gz)
            buoy.gz = 0.0;
        end
    end

    if isempty(buoy.gx), buoy.gx = 0.0; end
    if isempty(buoy.gy), buoy.gy = 0.0; end
    if isempty(buoy.gz), buoy.gz = 0.0; end
end

function line = strip_inline_comment(line)
    comment_chars = ['#', ';', '%', '!'];
    comment_idx = [];
    for i = 1:numel(comment_chars)
        idx = find(line == comment_chars(i), 1);
        if ~isempty(idx) && (isempty(comment_idx) || idx < comment_idx)
            comment_idx = idx;
        end
    end

    if ~isempty(comment_idx)
        line = strtrim(line(1:comment_idx-1));
    end
end

function value = parse_numeric_value(raw)
    value = str2double(strtrim(raw));
    if isnan(value) || ~isfinite(value)
        warning('load_case_metadata:InvalidNumericValue', ...
            'Skipping invalid numeric value "%s" in case metadata.', raw);
        value = [];
    end
end

function value = parse_boolean_value(raw)
    switch lower(strtrim(raw))
        case {'1', 'true', 'yes', 'on'}
            value = true;
        case {'0', 'false', 'no', 'off'}
            value = false;
        otherwise
            warning('load_case_metadata:InvalidBooleanValue', ...
                'Skipping invalid boolean value "%s" in case metadata.', raw);
            value = [];
    end
end

function value = normalize_inner_product(raw)
    value = upper(strtrim(raw));
    if isempty(value)
        return;
    end

    if strncmp(value, 'L2', 2)
        value = 'L2';
    elseif strncmp(value, 'H10', 3)
        value = 'H10';
    elseif strncmp(value, 'HLM', 3)
        value = 'HLM';
    elseif strncmp(value, 'OFF', 3) || strncmp(value, 'NONE', 4)
        value = 'L2';
    else
        warning('load_case_metadata:UnsupportedInnerProduct', ...
            'Ignoring unsupported POD inner-product "%s".', raw);
        value = '';
    end
end

function value = normalize_mode0(raw)
    value = lower(strtrim(raw));
    if isempty(value)
        return;
    end

    if ~(strcmp(value, 'avg') || strcmp(value, 'state'))
        warning('load_case_metadata:UnsupportedMode0', ...
            'Ignoring unsupported pod:mode0 value "%s".', raw);
        value = '';
    end
end

function value = normalize_deim_mode(raw)
    value = lower(strtrim(raw));
    if isempty(value)
        return;
    end

    if strcmp(value, 'off') || strcmp(value, 'none')
        value = 'fom';
    end

    valid_modes = {'fom', 'ftensor', 'rtensor', 'deim', 'clsdeim', 'mclsdeim'};
    if ~ismember(value, valid_modes)
        warning('load_case_metadata:UnsupportedDeimMode', ...
            'Ignoring unsupported deim:mode value "%s".', raw);
        value = '';
    end
end

function value = read_optional_scalar_file(path)
    value = [];
    fid = fopen(path, 'r');
    if fid < 0
        return;
    end
    cleanup = onCleanup(@() fclose(fid));
    raw = fscanf(fid, '%f', 1);
    if isempty(raw) || ~isfinite(raw)
        return;
    end
    value = raw;
end
