function rom_data = setup_conv_deim(pod_u, pod_v, nl_bas, nl_snaps_u, nl_snaps_v, x, y, ndeim_pts, n_os_points, ps_alg)
    % SETUP_ROM_CONVECTION Precomputes matrices for DEIM, CLS-DEIM, and MCLS-DEIM
    
    [nL, nb] = size(pod_u);
    nx1 = size(x, 1);
    
    % --- 1. Spectral Gradient Precomputation ---
    [zi, w] = zwgll(nx1-1);
    d_mat = deriv_mat(zi);
    [~,~,~,~,rx,ry,sx,sy,~,jaci,d] = deriv_geo(x, y, d_mat);
    my_lgrad = @(u) grad(u, rx, ry, sx, sy, jaci, d, 0);
    
    % Pre-allocate and compute gradients for all POD modes
    [ux_pods, uy_pods, vx_pods, vy_pods] = deal(zeros(nL, nb));
    for i = 1:nb
        [ux, uy] = my_lgrad(reshape(pod_u(:,i), size(x)));
        [vx, vy] = my_lgrad(reshape(pod_v(:,i), size(x)));
        ux_pods(:,i) = ux(:); uy_pods(:,i) = uy(:);
        vx_pods(:,i) = vx(:); vy_pods(:,i) = vy(:);
    end

    % --- 2. Point Selection (DEIM Indices) ---
    total_pts = ndeim_pts + n_os_points;
    switch lower(ps_alg)
        case 'sopt'
            inds = s_opt(nl_bas, total_pts, [])';
        case {'gpode','qdeim'}
            inds = gpode(nl_bas, total_pts)';
        case {'gappy_pod','deim'}
            inds = gappy_pod(nl_bas, total_pts);
        case 'gnat'
            inds = gnat(nl_bas, ndeim_pts, total_pts);
        otherwise
            error('Unknown point selection algorithm: %s', ps_alg);
    end
    inds = inds(1:ndeim_pts);
    rom_data.inds = inds;
    rom_data.nl_bas_p = nl_bas(inds, :);

    % --- 3. Projection & Mean-Flow Interaction ---
    % Integration weights (Mass Matrix)
    Me = reshape(jaci.^-1 .* (w * w'), nL, 1);
    Me_pod = [Me .* pod_u(:, 2:end); Me .* pod_v(:, 2:end)];
    rom_data.proj_mat = Me_pod' * nl_bas;
    
    % Precompute Zeroth Mode (Mean Flow) interaction terms
    c2 = Me_pod' * ([pod_u; pod_u] .* [ux_pods(:,1); vx_pods(:,1)] + ...
                    [pod_v; pod_v] .* [uy_pods(:,1); vy_pods(:,1)]);
    c3 = Me_pod' * ([pod_u(:,1); pod_u(:,1)] .* [ux_pods; vx_pods] + ...
                    [pod_v(:,1); pod_v(:,1)] .* [uy_pods; vy_pods]);
    rom_data.zmc = c2 + c3;
    rom_data.zmc(:,1) = rom_data.zmc(:,1) / 2;

    % --- 4. Sparse Stacks (Only storage for DEIM points) ---
    rom_data.u_p  = [pod_u; pod_u];     rom_data.u_p  = rom_data.u_p(inds, :);
    rom_data.v_p  = [pod_v; pod_v];     rom_data.v_p  = rom_data.v_p(inds, :);
    rom_data.ux_p = [ux_pods; vx_pods]; rom_data.ux_p = rom_data.ux_p(inds, :);
    rom_data.uy_p = [uy_pods; vy_pods]; rom_data.uy_p = rom_data.uy_p(inds, :);

    % --- 5. CLS and MCLS Matrix Inversions ---
    % Standard CLS matrix
    % Using \ instead of inv() for numerical stability
    rom_data.Ainv = (rom_data.nl_bas_p' * rom_data.nl_bas_p) \ eye(size(nl_bas, 2));
    rom_data.interp_mat = rom_data.Ainv * rom_data.nl_bas_p';

    % MCLS Statistics from training data
    if numel(nl_snaps_u) > 0 && numel(nl_snaps_v) > 0;
        nl_snapshot_proj = nl_bas' * [nl_snaps_u; nl_snaps_v];
        rom_data.mu  = mean(nl_snapshot_proj, 2);
        rom_data.tau = (cov(nl_snapshot_proj') + 1e-15*eye(size(nl_bas,2))) \ eye(size(nl_bas, 2)); 
        rom_data.alpha = 1e-12; % Regularization strength
    
        % Regularized matrix for MCLS
        rom_data.A_tau_inv = (rom_data.nl_bas_p' * rom_data.nl_bas_p + rom_data.alpha * rom_data.tau) \ eye(size(nl_bas, 2));
    end;
end
