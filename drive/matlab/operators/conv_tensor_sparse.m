function [out_coef] = conv_tensor_sparse(ucoef, pod_u, pod_v, x, y, validator)
    % De-aliased Sparse Convection Tensor with Automated Symmetric Validation
    
    persistent tensor nb_i nb_j nb_k
    force_skew = true; 

    if isempty(tensor)
        % --- 1. Setup & De-aliasing (3/2 Rule) ---
        [nx1, ny1, n_elem] = size(x);
        nb = size(pod_u, 2);
        nb_i = nb; nb_j = nb; nb_k = nb - 1;
        
        nx_f = ceil(1.5 * nx1); 
        nL_f_total = nx_f * nx_f * n_elem;
        [zi, ~] = zwgll(nx1-1);
        [zi_f, w_f] = zwgll(nx_f-1);
        d_f = deriv_mat(zi_f);
        Interp = interp_mat(zi_f, zi);

        % --- 2. Interpolate Weights & Geometry ---
        W2D = w_f * w_f';
        Me_f = zeros(nL_f_total, 1);
        RX = zeros(nx_f, nx_f, n_elem); RY = RX; SX = RX; SY = RX; JACI = RX;
        
        for ie = 1:n_elem
            xf_e = Interp * x(:,:,ie) * Interp';
            yf_e = Interp * y(:,:,ie) * Interp';
            [~,~,~,~,rx,ry,sx,sy,~,jaci,~] = deriv_geo(xf_e, yf_e, d_f);
            
            idx_f = (ie-1)*nx_f^2+1 : ie*nx_f^2;
            Me_f(idx_f) = (1./jaci(:)) .* W2D(:); 
            RX(:,:,ie) = rx; RY(:,:,ie) = ry; SX(:,:,ie) = sx; SY(:,:,ie) = sy; JACI(:,:,ie) = jaci;
        end

        % --- 3. Interpolate POD Modes ---
        pod_u_f = zeros(nL_f_total, nb); pod_v_f = zeros(nL_f_total, nb);
        for m = 1:nb
            for ie = 1:n_elem
                idx_c = (ie-1)*nx1^2+1 : ie*nx1^2;
                idx_f = (ie-1)*nx_f^2+1 : ie*nx_f^2;
                pod_u_f(idx_f, m) = reshape(Interp * reshape(pod_u(idx_c, m), nx1, nx1) * Interp', [], 1);
                pod_v_f(idx_f, m) = reshape(Interp * reshape(pod_v(idx_c, m), nx1, nx1) * Interp', [], 1);
            end
        end

        % --- 4. Assembly with Symmetric Validation ---
        % Test functions: modes 2 to nb (mapped to k=1:nb-1)
        pod_test_f = [bsxfun(@times, Me_f, pod_u_f(:, 2:nb)); ...
                      bsxfun(@times, Me_f, pod_v_f(:, 2:nb))];

        % Temporary dense slice storage for skew-enforcement
        tensor = sptensor([nb, nb, nb-1]);

        for i = 1:nb
            % Batch gradient for trial modes
            [ux_f, uy_f] = grad(reshape(pod_u_f, nx_f, nx_f, n_elem * nb), ...
                                repmat(RX, [1,1,nb]), repmat(RY, [1,1,nb]), ...
                                repmat(SX, [1,1,nb]), repmat(SY, [1,1,nb]), ...
                                repmat(JACI, [1,1,nb]), d_f, 0);
            [vx_f, vy_f] = grad(reshape(pod_v_f, nx_f, nx_f, n_elem * nb), ...
                                repmat(RX, [1,1,nb]), repmat(RY, [1,1,nb]), ...
                                repmat(SX, [1,1,nb]), repmat(SY, [1,1,nb]), ...
                                repmat(JACI, [1,1,nb]), d_f, 0);

            % Dense slice for current i to handle symmetry logic
            C_i = zeros(nb, nb-1);

            for j = 1:nb
                conv_ij = [pod_u_f(:,i).*reshape(ux_f(:,:,:,j),[],1) + pod_v_f(:,i).*reshape(uy_f(:,:,:,j),[],1); ...
                           pod_u_f(:,i).*reshape(vx_f(:,:,:,j),[],1) + pod_v_f(:,i).*reshape(vy_f(:,:,:,j),[],1)];
                
                % All k-projections for pair (i,j)
                C_i(j, :) = conv_ij' * pod_test_f;
            end

            if force_skew
                % Fluctuating indices: j in [2:nb], k in [1:nb-1]
                % These map to the same physical modes
                fluct_block = C_i(2:end, :); 
                C_i(2:end, :) = 0.5 * (fluct_block - fluct_block');
            end

            % Final sparsification: Apply validator and its symmetric counterpart
            for j = 1:nb
                for k = 1:nb-1
                    mode_j = j;
                    mode_k = k + 1; 
                    
                    % Check if either the interaction or its partner is valid
                    if validator(i, mode_j, mode_k) || (mode_j > 1 && validator(i, mode_k, mode_j))
                        if abs(C_i(j,k)) > 1e-15 % Avoid storing pure zeros
                            tensor(i, j, k) = C_i(j, k);
                        end
                    end
                end
            end
        end
    end

    % --- 5. Online Stage ---
    out_coef = double(ttv(tensor, {ucoef, ucoef}, [1, 2]));
end
