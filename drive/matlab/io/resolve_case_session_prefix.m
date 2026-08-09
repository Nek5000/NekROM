function prefix = resolve_case_session_prefix(case_path, casename)
% RESOLVE_CASE_SESSION_PREFIX Find the main Nek5000 session prefix for a case.
%
% Example: for casename='ann', returns something that matches ann0.f* files
% under either case_path/snaps/ or case_path/.

    if nargin < 1 || isempty(case_path)
        error('resolve_case_session_prefix:MissingCasePath', 'case_path is required.');
    end
    if nargin < 2 || isempty(casename)
        error('resolve_case_session_prefix:MissingCasename', 'casename is required.');
    end

    case_path = regexprep(case_path, '[\\/]+$', '');
    snaps_path = fullfile(case_path, 'snaps');

    candidates = {fullfile(snaps_path, casename), fullfile(case_path, casename)};
    prefix = '';
    for i = 1:numel(candidates)
        cand = candidates{i};
        files = dir([cand '0.f*']);
        if ~isempty(files)
            prefix = cand;
            return;
        end
    end
end

