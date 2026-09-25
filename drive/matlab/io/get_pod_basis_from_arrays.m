function [bas, u0, uk] = get_pod_basis_from_arrays(u_snaps, v_snaps, x, y, nb, ...
                                                    subtract_mean, conserve_momentum, method, inner_product, varargin)
    % GET_POD_BASIS_FROM_ARRAYS Performs POD with L2, H10, or HLM weighting.
    %
    % When subtract_mean is false, the first snapshot is used as the zeroth
    % mode reference so the returned basis still matches the Fortran layout.
    %
    % Optional HLM parameters (used only when inner_product = 'HLM'):
    %   varargin{1} - Reynolds number
    %   varargin{2} - Time step size
    %   varargin{3} - BDF beta(1,3) coefficient, default 11/6
    
    if nargin < 6; subtract_mean = true; end
    if nargin < 7; conserve_momentum = false; end
    if nargin < 8; method = 'snapshots'; end
    if nargin < 9; inner_product = 'L2'; end

    inner_product = upper(strtrim(inner_product));
    valid_inner_products = {'L2', 'H10', 'HLM'};
    if ~ismember(inner_product, valid_inner_products)
        error('get_pod_basis_from_arrays:UnsupportedInnerProduct', ...
            'Unsupported inner product "%s". Valid options: %s', ...
            inner_product, strjoin(valid_inner_products, ', '));
    end

    hlm_re = [];
    hlm_dt = [];
    hlm_beta1 = 11/6;
    if strcmp(inner_product, 'HLM')
        if numel(varargin) < 2
            error('get_pod_basis_from_arrays:HLMRequiresParams', ...
                'HLM inner product requires Reynolds number and dt arguments.');
        end

        hlm_re = varargin{1};
        hlm_dt = varargin{2};
        if numel(varargin) >= 3 && ~isempty(varargin{3})
            hlm_beta1 = varargin{3};
        end

        if ~isscalar(hlm_re) || ~isscalar(hlm_dt) || hlm_re <= 0 || hlm_dt <= 0
            error('get_pod_basis_from_arrays:InvalidHLMParams', ...
                'HLM parameters must be positive scalars: reynolds=%g dt=%g.', ...
                hlm_re, hlm_dt);
        end
    end

    % --- Mesh and Quadrature Setup ---
    nx1 = size(x,1);
    [~, w] = zwgll(nx1-1);
    d_mat = deriv_mat(zwgll(nx1-1)); 
    [~,~,~,~,~,~,~,~,jac,~,~] = deriv_geo(x, y, d_mat);
    
    nL = numel(x);
    W_diag = [reshape(jac.*(w*w'), nL, 1); reshape(jac.*(w*w'), nL, 1)];
    
    % --- Snapshot Prep ---
    snaps = [u_snaps; v_snaps];
    if subtract_mean
        mode0_ref = mean(snaps, 2);
    else
        mode0_ref = snaps(:, 1);
    end
    pod_snaps = snaps - mode0_ref;

    % --- Momentum Conservation ---
    if conserve_momentum
        eu = [ones(nL, 1); zeros(nL, 1)];
        ev = [zeros(nL, 1); ones(nL, 1)];
        E = [eu, ev];
        E_norm = sqrt(sum(E .* (W_diag .* E), 1));
        E = E ./ E_norm;
        pod_snaps = pod_snaps - E * (E' * (W_diag .* pod_snaps));
    end

    % --- Core POD Calculation ---
    if (strcmp(inner_product, 'H10') || strcmp(inner_product, 'HLM')) && strcmpi(method, 'svd')
        warning('Direct SVD is computationally prohibitive for H10/HLM. Falling back to ''snapshots''.');
        method = 'snapshots';
    end

    if strcmpi(method, 'svd')
        % Direct SVD approach
        W_sqrt = sqrt(W_diag);
        % Octave/MATLAB compatible sigma: 'lm' for Largest Magnitude
        [U, ~, ~] = svds(W_sqrt .* pod_snaps, nb, 'lm');
        bas = U ./ W_sqrt; 
        
    else
        % Method of Snapshots (Gramian)
        Gramian = compute_pod_gramian(pod_snaps, x, y, W_diag, inner_product, hlm_re, hlm_dt, hlm_beta1);
        
        Gramian = 0.5 * (Gramian + Gramian'); % Ensure symmetry
        
        % FIX: Use 'lm' (Largest Magnitude) instead of 'largestabs' for Octave
        [eigvecs, eigvals] = eigs(Gramian, nb, 'lm');
        
        % Projected modes: bas = X * eigvecs
        bas = pod_snaps * eigvecs;
        
        % Normalize basis
        if strcmp(inner_product, 'H10') || strcmp(inner_product, 'HLM')
            bas_metric = compute_pod_gramian(bas, x, y, W_diag, inner_product, hlm_re, hlm_dt, hlm_beta1);
            bas_norm = sqrt(max(real(diag(bas_metric)), 0))';
            bas_norm(bas_norm == 0) = eps;
        else
            bas_norm = sqrt(sum(bas .* (W_diag .* bas), 1));
        end
        bas = bas ./ bas_norm;
    end

    % --- Post-processing ---
    if conserve_momentum
        bas = [E, bas(:, 1:end-2)];
    end

    % --- Calculate Coefficients ---
    if strcmp(inner_product, 'H10') || strcmp(inner_product, 'HLM')
        % Handle potential matrix/vector eigvals differences
        if isvector(eigvals)
            S2 = diag(eigvals);
        else
            S2 = eigvals;
        end
        uk_fluc = diag(1 ./ bas_norm) * S2 * eigvecs';
    else
        % Standard L2 projection
        uk_fluc = bas' * (W_diag .* pod_snaps);
    end

    bas = [mode0_ref, bas];
    uk = [ones(1, size(snaps, 2)); uk_fluc];

    u0 = uk(:, 1);
end

function gramian = compute_pod_gramian(stacked_modes, x, y, W_diag, inner_product, hlm_re, hlm_dt, hlm_beta1)
    nrows = size(stacked_modes, 1);
    if mod(nrows, 2) ~= 0
        error('compute_pod_gramian:InvalidStackedModes', ...
            'Stacked velocity basis must have an even number of rows.');
    end

    nL = nrows / 2;
    u_modes = stacked_modes(1:nL, :);
    v_modes = stacked_modes(nL+1:end, :);

    switch inner_product
        case 'L2'
            gramian = stacked_modes' * (W_diag .* stacked_modes);

        case 'H10'
            [gramian, ~] = gen_Au(u_modes, v_modes, x, y);

        case 'HLM'
            [au_metric, bu_metric] = gen_Au(u_modes, v_modes, x, y);
            gramian = (1 / hlm_re) * au_metric + (hlm_beta1 / hlm_dt) * bu_metric;

        otherwise
            error('compute_pod_gramian:UnsupportedInnerProduct', ...
                'Unsupported inner product "%s".', inner_product);
    end
end
