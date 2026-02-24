function [out_coef, T_mat] = conv_tensor(ucoef, pod_u, pod_v, x, y)
    % Optimized Convection Tensor for NekROM
    % C(i,j,k) = < phi_i . grad(phi_j) , phi_k >
    
    persistent tensor nb n_test
    
    % --- Offline Stage (Pre-calculation) ---
    if isempty(tensor)
        [nx1, ny1] = size(x);
        nL = nx1 * ny1;
        nb = size(pod_u, 2);
        n_test = nb - 1; % Projection onto the POD modes 2:nb
        
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        
        % Get geometric factors. jaci = 1/J, jac = J.
        [~,~,~,~,rx,ry,sx,sy,jac,jaci,~] = deriv_geo(x, y, d);
        
        % Corrected weight: Me = J * w_i * w_j
        Me = reshape(jac .* (w * w'), nL, 1);
        
        % Pre-allocate gradient storage
        ux_p = zeros(nL, nb); uy_p = zeros(nL, nb);
        vx_p = zeros(nL, nb); vy_p = zeros(nL, nb);
        
        % Calculate all gradients once
        for i = 1:nb
            ui = reshape(pod_u(:, i), size(x));
            vi = reshape(pod_v(:, i), size(x));
            
            [ux, uy] = grad(ui, rx, ry, sx, sy, jaci, d, 0);
            [vx, vy] = grad(vi, rx, ry, sx, sy, jaci, d, 0);
            
            ux_p(:, i) = ux(:); uy_p(:, i) = uy(:);
            vx_p(:, i) = vx(:); vy_p(:, i) = vy(:);
        end

        if false;        
            % Performant version
            % Pre-weight the test functions (modes 2 to nb)
            % This combines the velocity components and mass matrix weights
            phi_k_weighted = [Me .* pod_u(:, 2:nb); Me .* pod_v(:, 2:nb)];
            
            % Assemble 3D Tensor C(i,j,k)
            % We iterate over i to keep memory overhead manageable
            tensor = zeros(nb, nb, n_test);
            for i = 1:nb
                % Convection term for a fixed i across all j: (phi_i . grad(phi_j))
                % Matrix size: (2*nL) x nb
                conv_ij = [bsxfun(@times, pod_u(:,i), ux_p) + bsxfun(@times, pod_v(:,i), uy_p);
                           bsxfun(@times, pod_u(:,i), vx_p) + bsxfun(@times, pod_v(:,i), vy_p)];
                
                % Project onto test space k: result is nb x n_test
                % This uses matrix-matrix multiply (BLAS3)
                tensor(i, :, :) = (conv_ij' * phi_k_weighted);
            end
            
            % Permute for optimal online contraction: [k, i, j]
            %tensor = permute(tensor, [3, 1, 2]);
            tensor = permute(tensor, [3,2,1]);

        else
            % k: Test function index (The mode we project onto)
            for k = 1:n_test
                % Readable version
                % Get the weighted test function (mode k+1)
                % Reshaped for a dot product
                wk_u = Me .* pod_u(:, k+1);
                wk_v = Me .* pod_v(:, k+1);
                
                % i: Trial function index (The 'convecting' velocity)
                for i = 1:nb
                    ui = pod_u(:, i);
                    vi = pod_v(:, i);
                    
                    % j: Basis function index (The 'convected' field)
                    for j = 1:nb
                        % Inner product: (phi_k) dot (phi_i * grad(phi_j))
                        % C_kij = sum( wk * (ui*ux_j + vi*uy_j) )
                        conv_x = ui .* ux_p(:, j) + vi .* uy_p(:, j);
                        conv_y = ui .* vx_p(:, j) + vi .* vy_p(:, j);

                        %tensor(k, i, j) = wk_u' * conv_x + wk_v' * conv_y;
                        tensor(k, j, i) = wk_u' * conv_x + wk_v' * conv_y;
                    end
                end
            end

        end
    end
    
    % --- Online Stage (Performance Critical) ---
    % Goal: compute sum_{i,j} C(k,i,j) * u(i) * u(j)
    % We use a reshaped matrix-vector product for O(N^2) speed.
    
    % Contract over j: temp(k, i) = C(k, i, :) * ucoef
    T_mat = reshape(tensor, n_test * nb, nb);
    temp = reshape(T_mat * ucoef, n_test, nb);
   
    %out_coef = temp*ucoef;
    
    % Contract over i: out(k) = temp(k, i) * ucoef
    out_coef = sum(temp .* ucoef', 2);

    %out_coef = (reshape(tensor*ucoef(:,1),nb,nb+1)*ucoef(:,1));
end
