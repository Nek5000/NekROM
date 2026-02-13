function [bas, u0, uk] = get_pod_basis_from_arrays(u_snaps, v_snaps, x, y, nb, ...
                                        subtract_mean, conserve_momentum, method)
    % GET_POD_BASIS_FROM_ARRAYS Performs POD with mass-matrix weighting
    % method: 'snapshots' (default) or 'svd'
    
    if nargin < 6; subtract_mean = false; end
    if nargin < 7; conserve_momentum = false; end
    if nargin < 8; method = 'snapshots'; end

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
        % Project out momentum: P = (I - E*E'*M)
        pod_snaps = pod_snaps - E * (E' * (W_diag .* pod_snaps));
    end

    % --- Core POD Calculation ---
    if strcmpi(method, 'svd')
        % Direct SVD approach: SVD(W^0.5 * X)
        W_sqrt = sqrt(W_diag);
        [U, ~, ~] = svds(W_sqrt .* pod_snaps, nb, 'largest');
        bas = U ./ W_sqrt; % Transform back to physical space
        
    else
        % Method of Snapshots (Gramian): K = X' * M * X
        % K is [n_snaps x n_snaps], much smaller than the spatial grid
        Gramian = pod_snaps' * (W_diag .* pod_snaps);
        Gramian = 0.5 * (Gramian + Gramian'); % Ensure symmetry
        
        [eigvecs, eigvals] = eigs(Gramian, nb, 'largestabs');
        
        % Projected modes: bas = X * eigvecs * inv(sqrt(eigvals))
        % This ensures the basis is orthonormal in the M-inner product
        bas = pod_snaps * eigvecs;
        bas_norm = sqrt(sum(bas .* (W_diag .* bas), 1));
        bas = bas ./ bas_norm;
    end

    % --- Post-processing ---
    if conserve_momentum
        bas = [E, bas(:, 1:end-2)];
    end

    % Calculate coefficients
    % If mean was subtracted, uk(1,:) is the mean weight (1.0)
    if subtract_mean
        uk_fluc = bas' * (W_diag .* (snaps - avg_snaps));
        bas = [avg_snaps, bas];
        uk = [ones(1, size(snaps, 2)); uk_fluc];
    else
        uk = bas' * (W_diag .* snaps);
    end

    u0 = uk(:, 1);
end
