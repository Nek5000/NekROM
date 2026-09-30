function rom_data = setup_conv_deim(pod_u, pod_v, nl_bas, nl_snaps_u, nl_snaps_v, x, y, ndeim_pts, n_os_points, ps_alg, use_oversampled_points, use_full_quadrature, deim_alpha, use_compressed_quadrature)
    % SETUP_ROM_CONVECTION Precomputes matrices for DEIM, CLS-DEIM, and MCLS-DEIM
    %
    % Requires NekToolKit functions: zwgll, deriv_mat, deriv_geo, grad, interp_mat
    % Note: Dependencies are checked at driver startup

    if nargin < 11 || isempty(use_oversampled_points)
        use_oversampled_points = false;
    end
    if nargin < 12 || isempty(use_full_quadrature)
        use_full_quadrature = false;
    end
    if nargin < 13 || isempty(deim_alpha)
        deim_alpha = 1e-12;
    end
    if nargin < 14 || isempty(use_compressed_quadrature)
        use_compressed_quadrature = false;
    end

    if use_full_quadrature && use_compressed_quadrature
        warning('setup_conv_deim:QuadratureModeConflict', ...
            'Both full and compressed quadrature requested; using full quadrature.');
        use_compressed_quadrature = false;
    end
    
    nL = numel(x);
    nb = size(pod_u, 2);
    nx1 = size(x, 1);
    assert(size(pod_u, 1) == nL && size(pod_v, 1) == nL, ...
        'setup_conv_deim:podGridMismatch', 'POD basis must match the supplied grid.');
    assert(size(nl_bas, 1) == 2 * nL, ...
        'setup_conv_deim:nlGridMismatch', 'Nonlinear basis must be stacked on the supplied grid.');
    if ~isempty(nl_snaps_u) && ~isempty(nl_snaps_v)
        assert(size(nl_snaps_u, 1) == nL && size(nl_snaps_v, 1) == nL, ...
            'setup_conv_deim:snapsGridMismatch', 'Training snapshots must match the supplied grid.');
    end
    
    % --- 1. Spectral Gradient Precomputation ---
    [zi, w] = zwgll(nx1-1);
    d_mat = deriv_mat(zi);
    [~,~,~,~,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x, y, d_mat);
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
            inds = qdeim(nl_bas, total_pts)';
        case {'gappy_pod','deim'}
            inds = gappy_pod(nl_bas, total_pts);
        case 'gnat'
            inds = gnat(nl_bas, ndeim_pts, total_pts);
        otherwise
            error('Unknown point selection algorithm: %s', ps_alg);
    end
    inds = inds(:);
    assert(numel(inds) >= ndeim_pts, ...
        'setup_conv_deim:InsufficientPoints', ...
        'Point selection returned fewer than %d indices.', ndeim_pts);

    inds_primary = inds(1:ndeim_pts);
    if use_oversampled_points
        inds_eval = inds;
    else
        inds_eval = inds_primary;
    end

    rom_data.use_oversampled_points = logical(use_oversampled_points);
    rom_data.use_full_quadrature = logical(use_full_quadrature);
    rom_data.use_compressed_quadrature = logical(use_compressed_quadrature);
    rom_data.inds = inds_primary;
    rom_data.inds_os = inds;
    rom_data.eval_inds = inds_eval;
    rom_data.sample_count = numel(inds_eval);
    rom_data.eval_weights = ones(rom_data.sample_count, 1);
    rom_data.nl_bas_p = nl_bas(inds_primary, :);
    rom_data.nl_bas_p_os = nl_bas(inds, :);
    rom_data.nl_bas_p_eval = nl_bas(inds_eval, :);

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
    rom_data.u_p  = [pod_u; pod_u];     rom_data.u_p  = rom_data.u_p(inds_primary, :);
    rom_data.v_p  = [pod_v; pod_v];     rom_data.v_p  = rom_data.v_p(inds_primary, :);
    rom_data.ux_p = [ux_pods; vx_pods]; rom_data.ux_p = rom_data.ux_p(inds_primary, :);
    rom_data.uy_p = [uy_pods; vy_pods]; rom_data.uy_p = rom_data.uy_p(inds_primary, :);

    rom_data.eval_u_p  = [pod_u; pod_u];     rom_data.eval_u_p  = rom_data.eval_u_p(inds_eval, :);
    rom_data.eval_v_p  = [pod_v; pod_v];     rom_data.eval_v_p  = rom_data.eval_v_p(inds_eval, :);
    rom_data.eval_ux_p = [ux_pods; vx_pods]; rom_data.eval_ux_p = rom_data.eval_ux_p(inds_eval, :);
    rom_data.eval_uy_p = [uy_pods; vy_pods]; rom_data.eval_uy_p = rom_data.eval_uy_p(inds_eval, :);

    rom_data.u_p_os  = [pod_u; pod_u];     rom_data.u_p_os  = rom_data.u_p_os(inds, :);
    rom_data.v_p_os  = [pod_v; pod_v];     rom_data.v_p_os  = rom_data.v_p_os(inds, :);
    rom_data.ux_p_os = [ux_pods; vx_pods]; rom_data.ux_p_os = rom_data.ux_p_os(inds, :);
    rom_data.uy_p_os = [uy_pods; vy_pods]; rom_data.uy_p_os = rom_data.uy_p_os(inds, :);

    % --- 5. CLS and MCLS Matrix Inversions ---
    % Standard CLS matrix
    % Using \ instead of inv() for numerical stability
    gram = rom_data.nl_bas_p_eval' * rom_data.nl_bas_p_eval;
    rom_data.Ainv = gram \ eye(size(gram));
    rom_data.interp_mat = rom_data.Ainv * rom_data.nl_bas_p_eval';

    % MCLS Statistics from training data
    if ~isempty(nl_snaps_u) && ~isempty(nl_snaps_v)
        nl_snapshot_proj = nl_bas' * [nl_snaps_u; nl_snaps_v];
        rom_data.mu  = mean(nl_snapshot_proj, 2);
        rom_data.tau = (cov(nl_snapshot_proj') + 1e-15*eye(size(nl_bas,2))) \ eye(size(nl_bas, 2)); 
        rom_data.alpha = deim_alpha; % Regularization strength
    
        % Regularized matrix for MCLS
        gram_tau = rom_data.nl_bas_p_eval' * rom_data.nl_bas_p_eval + rom_data.alpha * rom_data.tau;
        rom_data.A_tau_inv = gram_tau \ eye(size(gram_tau));
    end

    if use_full_quadrature || use_compressed_quadrature
        nxq = ceil(1.5 * nx1);
        [zi_q, w_q] = zwgll(nxq-1);
        interp_q = interp_mat(zi_q, zi);

        x_q = zeros(nxq, nxq, size(x, 3));
        y_q = zeros(nxq, nxq, size(y, 3));
        jac_q = zeros(nxq, nxq, size(jac, 3));
        for ie = 1:size(x, 3)
            x_q(:,:,ie) = interp_q * x(:,:,ie) * interp_q';
            y_q(:,:,ie) = interp_q * y(:,:,ie) * interp_q';
            jac_q(:,:,ie) = interp_q * jac(:,:,ie) * interp_q';
        end

        [pod_u_q, pod_v_q] = interp_basis_to_grid(pod_u, pod_v, x, y, x_q, y_q);
        [ux_q, vx_q] = interp_basis_to_grid(ux_pods, vx_pods, x, y, x_q, y_q);
        [uy_q, vy_q] = interp_basis_to_grid(uy_pods, vy_pods, x, y, x_q, y_q);

        nL_q = numel(x_q);
        Me_q = reshape(jac_q .* reshape(w_q * w_q', nxq, nxq, 1), nL_q, 1);
        sqrt_Me_q = sqrt(Me_q);
        sqrt_Me_stack = [sqrt_Me_q; sqrt_Me_q];
        Me_stack = sqrt_Me_stack .^ 2;
        Me_pod_q = [Me_q .* pod_u_q(:, 2:end); Me_q .* pod_v_q(:, 2:end)];

        [nl_bas_u, nl_bas_v] = split_stacked_basis(nl_bas, nL);
        [nl_bas_u_q, nl_bas_v_q] = interp_basis_to_grid(nl_bas_u, nl_bas_v, x, y, x_q, y_q);
        nl_bas_q = [nl_bas_u_q; nl_bas_v_q];

        % Dealiased projection operators (offline cost only).
        rom_data.proj_mat = Me_pod_q' * nl_bas_q;

        c2 = Me_pod_q' * ([pod_u_q; pod_u_q] .* [ux_q(:,1); vx_q(:,1)] + ...
                          [pod_v_q; pod_v_q] .* [uy_q(:,1); vy_q(:,1)]);
        c3 = Me_pod_q' * ([pod_u_q(:,1); pod_u_q(:,1)] .* [ux_q; vx_q] + ...
                          [pod_v_q(:,1); pod_v_q(:,1)] .* [uy_q; vy_q]);
        rom_data.zmc = c2 + c3;
        rom_data.zmc(:,1) = rom_data.zmc(:,1) / 2;

        if use_full_quadrature
            nl_bas_q_w = bsxfun(@times, sqrt_Me_stack, nl_bas_q);

            rom_data.use_full_quadrature = true;
            rom_data.use_compressed_quadrature = false;
            rom_data.eval_inds = (1:size(nl_bas_q_w, 1))';
            rom_data.sample_count = size(nl_bas_q_w, 1);
            rom_data.eval_weights = sqrt_Me_stack;
            rom_data.eval_u_p = [pod_u_q; pod_u_q];
            rom_data.eval_v_p = [pod_v_q; pod_v_q];
            rom_data.eval_ux_p = [ux_q; vx_q];
            rom_data.eval_uy_p = [uy_q; vy_q];
            rom_data.nl_bas_p_eval = nl_bas_q_w;

            gram_q = nl_bas_q_w' * nl_bas_q_w;
            rom_data.Ainv = gram_q \ eye(size(gram_q));
            rom_data.interp_mat = rom_data.Ainv * nl_bas_q_w';
        else
            % Compressed quadrature on the overintegrated grid:
            % choose a small weighted point set so the online path stays close
            % to sampled DEIM while approximating the dealiased mass inner product.
            rom_data.use_full_quadrature = false;
            rom_data.use_compressed_quadrature = true;

            % Default: keep the point budget within the Fortran runtime limit.
            % ndeim_max = 3*lbnl, and nbnl is the runtime active count.
            nbnl = size(nl_bas_q, 2);
            cquad_mult = 3;
            env_mult = getenv('NEKROM_DEIM_CQUAD_MULT');
            if ~isempty(env_mult)
                tmp = str2double(env_mult);
                if isfinite(tmp) && tmp > 0
                    cquad_mult = tmp;
                end
            end
            cquad_npts = ceil(cquad_mult * nbnl);
            env_npts = getenv('NEKROM_DEIM_CQUAD_NPTS');
            if ~isempty(env_npts)
                tmp = str2double(env_npts);
                if isfinite(tmp) && tmp > 0
                    cquad_npts = round(tmp);
                end
            end
            cquad_npts = max(cquad_npts, nbnl);
            cquad_npts = min(cquad_npts, 3 * nbnl);
            cquad_npts = min(cquad_npts, size(nl_bas_q, 1));

            [cquad_inds, cquad_wts, cquad_info] = compressed_quadrature(nl_bas_q, Me_stack, cquad_npts, []);
            if numel(cquad_inds) < nbnl
                error('setup_conv_deim:CQuadTooFewPoints', ...
                    'Compressed quadrature produced only %d points; need at least %d.', ...
                    numel(cquad_inds), nbnl);
            end

            rom_data.cquad_info = cquad_info;
            rom_data.eval_inds = cquad_inds(:);
            rom_data.sample_count = numel(rom_data.eval_inds);
            rom_data.eval_weights = sqrt(cquad_wts(:));

            % Store only sampled rows (avoid materializing 2*nL_q stacks).
            rom_data.eval_u_p = stack_select(pod_u_q, pod_u_q, rom_data.eval_inds, nL_q);
            rom_data.eval_v_p = stack_select(pod_v_q, pod_v_q, rom_data.eval_inds, nL_q);
            rom_data.eval_ux_p = stack_select(ux_q, vx_q, rom_data.eval_inds, nL_q);
            rom_data.eval_uy_p = stack_select(uy_q, vy_q, rom_data.eval_inds, nL_q);
            rom_data.nl_bas_p_eval = bsxfun(@times, rom_data.eval_weights, nl_bas_q(rom_data.eval_inds, :));

            gram_cq = rom_data.nl_bas_p_eval' * rom_data.nl_bas_p_eval;
            rom_data.Ainv = gram_cq \ eye(size(gram_cq));
            rom_data.interp_mat = rom_data.Ainv * rom_data.nl_bas_p_eval';
        end

        if ~isempty(nl_snaps_u) && ~isempty(nl_snaps_v)
            % MCLS statistics are built against the dealiased (full) quadrature norm.
            [nl_snaps_u_q, nl_snaps_v_q] = interp_basis_to_grid(nl_snaps_u, nl_snaps_v, x, y, x_q, y_q);
            snap_q = [nl_snaps_u_q; nl_snaps_v_q];
            snap_q_w = bsxfun(@times, sqrt_Me_stack, snap_q);
            nl_bas_q_w_full = bsxfun(@times, sqrt_Me_stack, nl_bas_q);
            nl_snapshot_proj = nl_bas_q_w_full' * snap_q_w;
            rom_data.mu = mean(nl_snapshot_proj, 2);
            rom_data.tau = (cov(nl_snapshot_proj') + 1e-15 * eye(size(nl_bas_q, 2))) \ eye(size(nl_bas_q, 2));
            rom_data.alpha = deim_alpha;

            if rom_data.use_full_quadrature
                gram_base = nl_bas_q_w_full' * nl_bas_q_w_full;
            else
                gram_base = rom_data.nl_bas_p_eval' * rom_data.nl_bas_p_eval;
            end
            gram_tau = gram_base + rom_data.alpha * rom_data.tau;
            rom_data.A_tau_inv = gram_tau \ eye(size(gram_tau));
        end
    end
end

function [u_half, v_half] = split_stacked_basis(stacked_basis, nL)
    u_half = stacked_basis(1:nL, :);
    v_half = stacked_basis(nL+1:end, :);
end

function out = stack_select(first_half, second_half, inds, nL)
    inds = inds(:);
    out = zeros(numel(inds), size(first_half, 2));
    mask_first = inds <= nL;
    if any(mask_first)
        out(mask_first, :) = first_half(inds(mask_first), :);
    end
    mask_second = ~mask_first;
    if any(mask_second)
        out(mask_second, :) = second_half(inds(mask_second) - nL, :);
    end
end
