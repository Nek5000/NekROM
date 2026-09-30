function mode0_ref = load_pod_mode0_reference(snaps_obj, reorder)
% LOAD_POD_MODE0_REFERENCE Load the zeroth-mode reference for state-mode POD.
%
% For cases that store pod:mode0 = state, the Fortran offline path uses the
% case's reference stream rather than the raw snapshot mean. In practice the
% closest MATLAB match comes from the avg-stream snapshot when it exists, with
% uic as a fallback and the first raw snapshot as a last resort.

    mode0_ref = [];

    if nargin < 1 || isempty(snaps_obj) || ~isprop(snaps_obj, 'cname')
        return;
    end
    if nargin < 2 || isempty(reorder)
        reorder = false;
    end

    case_prefix = char(snaps_obj.cname);
    [case_snaps_dir, casename] = fileparts(case_prefix);
    case_dir = fileparts(case_snaps_dir);
    if isempty(case_dir) || isempty(casename)
        return;
    end

    candidate_prefixes = {['avg' casename], ['uic' casename]};
    for i = 1:numel(candidate_prefixes)
        candidate = NekSnaps(fullfile(case_dir, 'snaps', candidate_prefixes{i}));
        if isempty(candidate) || isempty(candidate.isnaps)
            continue;
        end

        [u_ref, v_ref] = get_snaps(candidate, reorder);
        if ~isempty(u_ref) && ~isempty(v_ref)
            mode0_ref = [u_ref(:, 1); v_ref(:, 1)];
            return;
        end
    end

    if isempty(snaps_obj.isnaps) || isempty(snaps_obj.flds)
        return;
    end

    [u_ref, v_ref] = get_snaps(snaps_obj, reorder);
    if ~isempty(u_ref) && ~isempty(v_ref)
        mode0_ref = [u_ref(:, 1); v_ref(:, 1)];
    end
end
