function prefix = resolve_snapshot_prefix(case_path, casename, tag)
% RESOLVE_SNAPSHOT_PREFIX Find a Nek5000 snapshot prefix for NekSnaps.
%
% NekROM examples are not fully consistent about whether snapshot/basis files
% live under case_path/ or case_path/snaps/, nor whether they include the
% casename in the prefix (e.g. bas0.f00001 vs basshear0.f00001).
%
% This helper searches common locations and naming patterns and returns the
% first matching prefix (without the trailing "0.fXXXXX").
%
% Example:
%   prefix = resolve_snapshot_prefix('../../examples/ann', 'ann', 'bas');
%   snaps = NekSnaps(prefix);

    if nargin < 1 || isempty(case_path)
        error('resolve_snapshot_prefix:MissingCasePath', 'case_path is required.');
    end
    if nargin < 2 || isempty(casename)
        casename = '';
    end
    if nargin < 3 || isempty(tag)
        error('resolve_snapshot_prefix:MissingTag', 'tag is required.');
    end

    case_path = regexprep(case_path, '[\\/]+$', '');
    snaps_path = fullfile(case_path, 'snaps');
    snaps_rom_path = fullfile(case_path, 'snaps_rom');

    candidates = {};
    if ~isempty(casename)
        candidates{end+1} = fullfile(snaps_path, [tag casename]);
        candidates{end+1} = fullfile(snaps_rom_path, [tag casename]);
        candidates{end+1} = fullfile(case_path, [tag casename]);
    end
    candidates{end+1} = fullfile(snaps_path, tag);
    candidates{end+1} = fullfile(snaps_rom_path, tag);
    candidates{end+1} = fullfile(case_path, tag);

    prefix = '';
    for i = 1:numel(candidates)
        cand = candidates{i};
        if has_field_files(cand)
            prefix = cand;
            return;
        end
    end
end

function ok = has_field_files(prefix)
    ok = false;
    files = dir([prefix '0.f*']);
    if ~isempty(files)
        ok = true;
        return;
    end
end
