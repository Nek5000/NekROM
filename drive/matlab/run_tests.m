function run_tests(snaps_path, casename, nb, reorder, pod_u, pod_v, au_full, bu_full, inner_product, varargin)
% RUN_TESTS Code to generate the basis functions in Matlab to validate Fortran output
%
% This ensures MATLAB and NekROM outputs align properly.

    fprintf('--- Running Validation Tests ---\n');

    subtract_mean = true;
    conserve_momentum = 0;
    method = "snapshots";

    if nargin < 9 || isempty(inner_product)
        inner_product = 'H10';
    end
    inner_product = upper(strtrim(inner_product));

    % Load Snapshots
    snaps_obj = NekSnaps(fullfile(snaps_path, casename)); 
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
        [pod_ml, ~, ~] = get_pod_basis(snaps_obj, nb, reorder, subtract_mean, conserve_momentum, ...
            method, inner_product, hlm_re, hlm_dt, hlm_beta1);
    else
        if numel(varargin) >= 1 && ~isempty(varargin{1})
            subtract_mean = logical(varargin{1});
        end
        [pod_ml, ~, ~] = get_pod_basis(snaps_obj, nb, reorder, subtract_mean, conserve_momentum, method, inner_product);
    end

    pod_u_ml = pod_ml(1:size(pod_ml,1)/2, 1:nb+1);
    pod_v_ml = pod_ml(size(pod_ml,1)/2 + 1:end, 1:nb+1);

    [x_fom_ml, y_fom_ml] = get_grid(snaps_obj, reorder);
    [au_full_ml, bu_full_ml] = gen_Au(pod_u_ml, pod_v_ml, x_fom_ml, y_fom_ml);

    %% Test 1: Basis vectors are the same (modulo sign differences)
    pod = [pod_u; pod_v];
    pod_diff = norm(abs(pod) - abs(pod_ml)) / norm(abs(pod));
    
    if pod_diff < 1e-5
        fprintf('PASS: POD Basis vectors match (diff: %e)\n', pod_diff);
    else
        warning('FAIL: POD Basis vectors mismatch! (diff: %e)', pod_diff);
    end

    %% Test 2: MATLAB and Fortran Au and Bu operators match
    bu_diff = norm(abs(bu_full) - abs(bu_full_ml)) / norm(abs(bu_full));
    au_diff = norm(abs(au_full) - abs(au_full_ml)) / norm(abs(au_full));
    
    if bu_diff < 1e-5
        fprintf('PASS: Bu operators match (diff: %e)\n', bu_diff);
    else
        warning('FAIL: Bu operators mismatch! (diff: %e)', bu_diff);
    end
    
    if au_diff < 1e-5
        fprintf('PASS: Au operators match (diff: %e)\n', au_diff);
    else
        warning('FAIL: Au operators mismatch! (diff: %e)', au_diff);
    end

    fprintf('--- Validation Tests Complete ---\n\n');
end
