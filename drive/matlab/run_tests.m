function run_tests(snaps_path, casename, nb, reorder, pod_u, pod_v, au_full, bu_full, inner_product, varargin)
% RUN_TESTS Code to generate the basis functions in Matlab to validate Fortran output
%
% This ensures MATLAB and NekROM outputs align properly.

    fprintf('--- Running Validation Tests ---\n');

    subtract_mean = true;
    conserve_momentum = 0;
    method = "snapshots";
    mode_energy_tol = 1e-3;
    pod_tol = 1e-2;

    if nargin < 9 || isempty(inner_product)
        inner_product = 'H10';
    end
    inner_product = upper(strtrim(inner_product));

    if numel(varargin) >= 1 && ~isempty(varargin{1})
        subtract_mean = logical(varargin{1});
    end

    % Load Snapshots
    snaps_obj = NekSnaps(fullfile(snaps_path, casename));
    snaps_obj = apply_snapshot_file_list_subset(snaps_obj);
    if strcmp(inner_product, 'HLM')
        if numel(varargin) == 2
            hlm_re = varargin{1};
            hlm_dt = varargin{2};
            hlm_beta1 = 11/6;
        elseif numel(varargin) >= 3
            subtract_mean = logical(varargin{1});
            hlm_re = varargin{2};
            hlm_dt = varargin{3};
            hlm_beta1 = 11/6;
            if numel(varargin) >= 4 && ~isempty(varargin{4})
                hlm_beta1 = varargin{4};
            end
        else
            error('run_tests:HLMRequiresParams', ...
                'HLM validation requires subtract_mean, Reynolds number, and dt arguments.');
        end
        [pod_ml, ~, uk_ml] = get_pod_basis(snaps_obj, nb, reorder, subtract_mean, conserve_momentum, ...
            method, inner_product, hlm_re, hlm_dt, hlm_beta1);
    else
        if numel(varargin) >= 1 && ~isempty(varargin{1})
            subtract_mean = logical(varargin{1});
        end
        [pod_ml, ~, uk_ml] = get_pod_basis(snaps_obj, nb, reorder, subtract_mean, conserve_momentum, method, inner_product);
    end

    if isempty(uk_ml) || size(uk_ml, 1) < 2
        sig_mode_count = nb;
    else
        mode_energy = sum(abs(uk_ml(2:end, :)).^2, 2);
        if isempty(mode_energy) || all(mode_energy <= 0)
            sig_mode_count = nb;
        else
            sig_mode_count = find(mode_energy >= max(mode_energy) * mode_energy_tol, 1, 'last');
            if isempty(sig_mode_count)
                sig_mode_count = min(1, nb);
            end
        end
    end

    compare_cols = [1, 2:(sig_mode_count + 1)];
    pod_u_ml = pod_ml(1:size(pod_ml,1)/2, compare_cols);
    pod_v_ml = pod_ml(size(pod_ml,1)/2 + 1:end, compare_cols);

    [x_fom_ml, y_fom_ml] = get_grid(snaps_obj, reorder);
    [au_full_ml, bu_full_ml] = gen_Au(pod_u_ml, pod_v_ml, x_fom_ml, y_fom_ml);

    %% Test 1: Basis vectors are the same (modulo sign differences)
    pod = [pod_u; pod_v];
    pod_diff = norm(abs(pod(:, compare_cols)) - abs(pod_ml(:, compare_cols))) / norm(abs(pod(:, compare_cols)));

    if pod_diff < pod_tol
        fprintf('PASS: POD Basis vectors match for %d energetic modes (diff: %e)\n', sig_mode_count, pod_diff);
    else
        warning('FAIL: POD Basis vectors mismatch! (diff: %e)', pod_diff);
    end

    if sig_mode_count < nb
        fprintf('NOTE: Ignoring %d numerically null POD mode(s) beyond the energetic subspace.\n', nb - sig_mode_count);
    end

    %% Test 2: MATLAB and Fortran Au and Bu operators match
    bu_diff = norm(abs(bu_full(compare_cols, compare_cols)) - abs(bu_full_ml(compare_cols, compare_cols))) / norm(abs(bu_full(compare_cols, compare_cols)));
    au_diff = norm(abs(au_full(compare_cols, compare_cols)) - abs(au_full_ml(compare_cols, compare_cols))) / norm(abs(au_full(compare_cols, compare_cols)));

    if bu_diff < 1e-5
        fprintf('PASS: Bu operators match for %d energetic modes (diff: %e)\n', sig_mode_count, bu_diff);
    else
        warning('FAIL: Bu operators mismatch! (diff: %e)', bu_diff);
    end

    if au_diff < 1e-5
        fprintf('PASS: Au operators match for %d energetic modes (diff: %e)\n', sig_mode_count, au_diff);
    else
        warning('FAIL: Au operators mismatch! (diff: %e)', au_diff);
    end

    fprintf('--- Validation Tests Complete ---\n\n');
end
