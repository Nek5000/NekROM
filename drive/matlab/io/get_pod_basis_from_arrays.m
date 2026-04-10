function [bas, u0, uk] = get_pod_basis_from_arrays(u_snaps, v_snaps, x, y, nb, ...
                                                    subtract_mean, conserve_momentum, method, inner_product)
    % GET_POD_BASIS_FROM_ARRAYS Performs POD with mass-matrix or diffusion weighting
    % Updated for Octave/MATLAB compatibility (2026)
    
    if nargin < 6; subtract_mean = true; end
    if nargin < 7; conserve_momentum = false; end
    if nargin < 8; method = 'snapshots'; end
    if nargin < 9; inner_product = 'L2'; end

    % --- Mesh and Quadrature Setup ---
    nx1 = size(x,1);
    [~, w] = zwgll(nx1-1);
    d_mat = deriv_mat(zwgll(nx1-1)); 
    [~,~,~,~,~,~,~,~,jac,~,~] = deriv_geo(x, y, d_mat);
    
    nL = numel(x);
    W_diag = [reshape(jac.*(w*w'), nL, 1); reshape(jac.*(w*w'), nL, 1)];
    
    % --- Snapshot Prep ---
    snaps = [u_snaps; v_snaps];
    avg_snaps = mean(snaps, 2);
    pod_snaps = snaps - (subtract_mean * avg_snaps);

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
    if strcmpi(inner_product, 'H10') && strcmpi(method, 'svd')
        warning('Direct SVD is computationally prohibitive for H10. Falling back to ''snapshots''.');
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
        if strcmpi(inner_product, 'H10')
            Gramian = gen_Au(pod_snaps(1:nL, :), pod_snaps(nL+1:end, :), x, y);
        else
            Gramian = pod_snaps' * (W_diag .* pod_snaps);
        end
        
        Gramian = 0.5 * (Gramian + Gramian'); % Ensure symmetry
        
        % FIX: Use 'lm' (Largest Magnitude) instead of 'largestabs' for Octave
        [eigvecs, eigvals] = eigs(Gramian, nb, 'lm');
        
        % Projected modes: bas = X * eigvecs
        bas = pod_snaps * eigvecs;
        
        % Normalize basis
        if strcmpi(inner_product, 'H10')
            bas_A = gen_Au(bas(1:nL, :), bas(nL+1:end, :), x, y);
            bas_norm = sqrt(diag(bas_A))'; 
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
    if strcmpi(inner_product, 'H10')
        % Handle potential matrix/vector eigvals differences
        if isvector(eigvals)
            S2 = diag(eigvals);
        else
            S2 = eigvals;
        end
        uk_fluc = diag(1 ./ bas_norm) * S2 * eigvecs';
    else
        % Standard L2 projection
        if subtract_mean
            uk_fluc = bas' * (W_diag .* (snaps - avg_snaps));
        else
            uk_fluc = bas' * (W_diag .* snaps);
        end
    end

    if subtract_mean
        bas = [avg_snaps, bas];
        uk = [ones(1, size(snaps, 2)); uk_fluc];
    else
        uk = uk_fluc;
    end

    u0 = uk(:, 1);
end
