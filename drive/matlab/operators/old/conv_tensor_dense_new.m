function [out_coef] = conv_tensor_dense(ucoef, pod_u, pod_v, x, y, tensor_size)
    % PERSISTENT variables to cache the heavy lifting
    persistent tensor nb_dims
    
    % Initialize if first run
    if isempty(tensor)
        [nL, nb] = size(pod_u);
        nx1 = size(x,1);
        
        % --- Pre-processing Geometry ---
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [~,~,~,~,rx,ry,sx,sy,~,jaci,d] = deriv_geo(x,y,d);
        Me = reshape(jaci.^-1 .* (w*w'), nL, 1); % Integration weights
        
        % Handle tensor sizes
        if nargin < 6
            nb_dims = [nb, nb, nb-1];
        else
            nb_dims = tensor_size;
        end
        
        % --- Gradient Computation (Vectorized) ---
        % Reshape pod_u to act on the whole basis at once if lgrad permits,
        % otherwise, keep the loop but pre-allocate.
        ux_pods = zeros(nL, nb_dims(2)); 
        uy_pods = zeros(nL, nb_dims(2));
        vx_pods = zeros(nL, nb_dims(2)); 
        vy_pods = zeros(nL, nb_dims(2));

        for j = 1:nb_dims(2)
            [ux, uy] = grad(reshape(pod_u(:,j), size(x)), rx, ry, sx, sy, jaci, d, 0);
            [vx, vy] = grad(reshape(pod_v(:,j), size(x)), rx, ry, sx, sy, jaci, d, 0);
            ux_pods(:,j) = ux(:); uy_pods(:,j) = uy(:);
            vx_pods(:,j) = vx(:); vy_pods(:,j) = vy(:);
        end

        % Pre-weight the test functions (Weak form)
        % Using modes 2:end (assuming mode 1 is mean or discarded)
        W = [Me .* pod_u(:, 2:nb_dims(3)+1); Me .* pod_v(:, 2:nb_dims(3)+1)];
        
        % --- Tensor Assembly ---
        % Logic: T_ijk = Integral( Phi_k * (Phi_i * grad(Phi_j)) )
        tensor = zeros(nb_dims(1), nb_dims(2), nb_dims(3));
        for i = 1:nb_dims(1)
            % Convection term for this i-mode across all j-modes
            conv_x = pod_u(:,i) .* ux_pods + pod_v(:,i) .* uy_pods;
            conv_y = pod_u(:,i) .* vx_pods + pod_v(:,i) .* vy_pods;
            
            % Project onto the k-basis: [2*nL x nb_j]' * [2*nL x nb_k]
            % This results in a [nb_j x nb_k] slice
            tensor(i,:,:) = [conv_x; conv_y]' * W;
        end
    end

    % --- Online Phase (Fast) ---
    % Contract: out = u' * Tensor * u
    % Flattening for speed: [nb_i] * [nb_i x (nb_j*nb_k)] -> [1 x nb_j*nb_k]
    u_i = ucoef(1:nb_dims(1));
    u_j = ucoef(1:nb_dims(2));
    
    T_flattened = reshape(tensor, nb_dims(1), []);
    reduced_mat = reshape(u_i' * T_flattened, nb_dims(2), nb_dims(3));
    
    out_coef = (u_j' * reduced_mat)';
end
