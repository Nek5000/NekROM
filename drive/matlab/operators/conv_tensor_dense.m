function [out_coef] = conv_tensor_dense(ucoef, pod_u, pod_v, x, y, tensor_size)
    % C_ijk = < phi_k , (phi_i . grad) phi_j >
    
    persistent tensor nb_i nb_j nb_k
    force_skew = true; 

    if isempty(tensor)
        % --- 1. Grid & Size Setup ---
        [nx1, ny1, n_elem] = size(x);
        nb = size(pod_u, 2);
        
        if nargin < 6
            nb_i = nb; nb_j = nb; nb_k = nb - 1;
        else
            nb_i = tensor_size(1); nb_j = tensor_size(2); nb_k = tensor_size(3);
        end

        % --- 2. De-aliasing Operators (3/2 Rule) ---
        nx_f = ceil(1.5 * nx1);
        nL_f_total = nx_f * nx_f * n_elem;
        [zi, ~] = zwgll(nx1-1);
        [zi_f, w_f] = zwgll(nx_f-1);
        d_f = deriv_mat(zi_f);
        Interp = interp_mat(zi_f, zi);

        % --- 3. Interpolate Geometry & Setup Weights ---
        W2D = w_f * w_f';
        Me_f = zeros(nx_f, nx_f, n_elem);
        for ie = 1:n_elem
            xf_e = Interp * x(:,:,ie) * Interp';
            yf_e = Interp * y(:,:,ie) * Interp';
            [~,~,~,~,rx,ry,sx,sy,jac,jaci,~] = deriv_geo(xf_e, yf_e, d_f);
            Me_f((ie-1)*nx_f^2+1 : ie*nx_f^2) = jac(:) .* W2D(:);
            % Store geometric factors for grad (simplified for this block)
            RX(:,:,ie) = rx; RY(:,:,ie) = ry; SX(:,:,ie) = sx; SY(:,:,ie) = sy; JACI(:,:,ie) = jaci;
        end
        Me_f = Me_f(:);

        % --- 4. Interpolate POD Modes ---
        max_mode = max([nb_i, nb_j, nb_k + 1]);
        pod_u_f = zeros(nL_f_total, max_mode);
        pod_v_f = zeros(nL_f_total, max_mode);
        
        for m = 1:max_mode
            for ie = 1:n_elem
                idx_c = (ie-1)*nx1^2+1 : ie*nx1^2;
                idx_f = (ie-1)*nx_f^2+1 : ie*nx_f^2;
                pod_u_f(idx_f, m) = reshape(Interp * reshape(pod_u(idx_c, m), nx1, nx1) * Interp', [], 1);
                pod_v_f(idx_f, m) = reshape(Interp * reshape(pod_v(idx_c, m), nx1, nx1) * Interp', [], 1);
            end
        end

        % --- 5. Gradient Computation & Assembly ---
        tensor = zeros(nb_i, nb_j, nb_k);
        % Test functions are modes 2 to nb_k+1
        pod_test_f = [bsxfun(@times, Me_f, pod_u_f(:, 2:nb_k+1)); ...
                      bsxfun(@times, Me_f, pod_v_f(:, 2:nb_k+1))];

        for i = 1:nb_i
            % Gradient of all j-modes
            [ux_f, uy_f] = grad(reshape(pod_u_f(:,1:nb_j), nx_f, nx_f, n_elem * nb_j), ...
                                repmat(RX, [1,1,nb_j]), repmat(RY, [1,1,nb_j]), ...
                                repmat(SX, [1,1,nb_j]), repmat(SY, [1,1,nb_j]), ...
                                repmat(JACI, [1,1,nb_j]), d_f, 0);
            [vx_f, vy_f] = grad(reshape(pod_v_f(:,1:nb_j), nx_f, nx_f, n_elem * nb_j), ...
                                repmat(RX, [1,1,nb_j]), repmat(RY, [1,1,nb_j]), ...
                                repmat(SX, [1,1,nb_j]), repmat(SY, [1,1,nb_j]), ...
                                repmat(JACI, [1,1,nb_j]), d_f, 0);

            conv_j = [bsxfun(@times, pod_u_f(:,i), reshape(ux_f, nL_f_total, nb_j)) + ...
                      bsxfun(@times, pod_v_f(:,i), reshape(uy_f, nL_f_total, nb_j)); ...
                      bsxfun(@times, pod_u_f(:,i), reshape(vx_f, nL_f_total, nb_j)) + ...
                      bsxfun(@times, pod_v_f(:,i), reshape(vy_f, nL_f_total, nb_j))];

            % C_slice size: [nb_j, nb_k]
            C_slice = conv_j' * pod_test_f;

            if force_skew
                % Skew-symmetry is only valid for the overlapping fluctuating modes.
                % Trial modes j (starting from 2) must match Test modes k.
                % Since k=1 is mode 2, k=2 is mode 3... we look at C_slice(2:end, :)
                n_overlap = min(nb_j - 1, nb_k);
                sub_slice = C_slice(2:n_overlap+1, 1:n_overlap);
                C_slice(2:n_overlap+1, 1:n_overlap) = 0.5 * (sub_slice - sub_slice');
            end
            
            tensor(i, :, :) = C_slice;
        end
    end

    % --- 6. Online Stage (Optimized MatVec) ---
    out_coef = zeros([size(pod_u,2)-1, 1]);
    T_mat_i = reshape(tensor, nb_i, nb_j * nb_k);
    temp_jk = reshape(ucoef(1:nb_i)' * T_mat_i, nb_j, nb_k);
    out_coef(1:nb_k) = temp_jk' * ucoef(1:nb_j);
end
