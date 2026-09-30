function snaps_obj = apply_snapshot_file_list_subset(snaps_obj)
% APPLY_SNAPSHOT_FILE_LIST_SUBSET Restrict a NekSnaps object to file.list entries.
%
% NekROM offline generation reads snapshots in file.list order and may only
% consume the first ops/ns entries. NekToolKit's NekSnaps loader scans the
% directory directly, so this helper trims/reorders the in-memory snapshot
% list to match the Fortran offline path.

    if nargin < 1 || isempty(snaps_obj) || ~isprop(snaps_obj, 'cname')
        return;
    end

    case_snaps_prefix = char(snaps_obj.cname);
    case_snaps_dir = fileparts(case_snaps_prefix);
    case_dir = fileparts(case_snaps_dir);
    if isempty(case_dir)
        return;
    end

    file_list_path = fullfile(case_dir, 'file.list');
    if ~exist(file_list_path, 'file')
        return;
    end

    ops_ns = [];
    ops_ns_path = fullfile(case_dir, 'ops', 'ns');
    if exist(ops_ns_path, 'file')
        ops_ns = read_optional_scalar_file(ops_ns_path);
        if ~isempty(ops_ns)
            ops_ns = round(ops_ns);
        end
    end

    listed_steps = read_file_list_steps(file_list_path);
    if isempty(listed_steps)
        return;
    end

    if ~isempty(ops_ns) && ops_ns > 0
        listed_steps = listed_steps(1:min(ops_ns, numel(listed_steps)));
    end

    if isempty(snaps_obj.isnaps) || isempty(snaps_obj.flds)
        return;
    end

    grid_template = [];
    for i = 1:numel(snaps_obj.flds)
        if isfield(snaps_obj.flds{i}, 'x') && isfield(snaps_obj.flds{i}, 'y')
            grid_template = snaps_obj.flds{i};
            break;
        end
    end

    [is_listed, ordered_idx] = ismember(listed_steps, snaps_obj.isnaps);
    ordered_idx = ordered_idx(is_listed);
    if isempty(ordered_idx)
        return;
    end

    snaps_obj.isnaps = snaps_obj.isnaps(ordered_idx);
    snaps_obj.flds = snaps_obj.flds(ordered_idx);

    if ~isempty(grid_template)
        has_grid = isfield(snaps_obj.flds{1}, 'x') && isfield(snaps_obj.flds{1}, 'y');
        if ~has_grid
            snaps_obj.flds{1}.x = grid_template.x;
            snaps_obj.flds{1}.y = grid_template.y;
            if isfield(grid_template, 'z')
                snaps_obj.flds{1}.z = grid_template.z;
            end
        end
    end
end

function steps = read_file_list_steps(file_list_path)
    steps = [];

    fid = fopen(file_list_path, 'r');
    if fid < 0
        return;
    end
    cleanup = onCleanup(@() fclose(fid));

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end

        line = strtrim(line);
        if isempty(line) || line(1) == '#'
            continue;
        end

        tok = regexp(line, '\.f(\d+)$', 'tokens', 'once');
        if isempty(tok)
            continue;
        end

        steps(end+1, 1) = str2double(tok{1}); %#ok<AGROW>
    end
end

function value = read_optional_scalar_file(path)
    value = [];
    fid = fopen(path, 'r');
    if fid < 0
        return;
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    raw = fscanf(fid, '%f', 1);
    if isempty(raw) || ~isfinite(raw)
        return;
    end
    value = raw;
end
