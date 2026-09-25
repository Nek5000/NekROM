function [out_coef, T_mat] = conv_tensor_dealiased(ucoef, pod_u, pod_v, x, y)
    % Optimized De-aliased Convection Tensor for NekROM [N, N, E]
    
    persistent tensor nb n_test
    
    if isempty(tensor)
        % --- Configuration ---
        use_vectorized = true; % Set to true for performance, false for readability
        
        [nx1, ny1, n_elem] = size(x); 
        nb = size(pod_u, 2);
        n_test = nb - 1; 
        
        % 3/2 Rule for de-aliasing (e.g., 9 -> 14)
        nx_f = ceil(1.5 * nx1); 
        nL_f_total = nx_f * nx_f * n_elem;
        
        % --- Operators ---
        [zi, ~] = zwgll(nx1-1);           
        [zi_f, w_f] = zwgll(nx_f-1);      
        d_f = deriv_mat(zi_f);            
        Interp = interp_mat(zi_f, zi);    
        
        % --- Geometry & Interpolation ---
        x_f = zeros(nx_f, nx_f, n_elem);
        y_f = zeros(nx_f, nx_f, n_elem);
        for ie = 1:n_elem
            x_f(:,:,ie) = Interp * x(:,:,ie) * Interp';
            y_f(:,:,ie) = Interp * y(:,:,ie) * Interp';
        end
        
        [~,~,~,~,rx_f,ry_f,sx_f,sy_f,jac_f,jaci_f,~] = deriv_geo(x_f, y_f, d_f);
        
        % Fine weights: [nx_f*nx_f*n_elem, 1]
        W2D = w_f * w_f';
        Me_f = zeros(nx_f, nx_f, n_elem);
        for ie = 1:n_elem
            Me_f(:,:,ie) = jac_f(:,:,ie) .* W2D;
        end
        Me_f = Me_f(:); 
        
        % Interpolate POD to Fine Grid
        pod_u_f = zeros(nL_f_total, nb);
        pod_v_f = zeros(nL_f_total, nb);
        for i = 1:nb
            ui_c = reshape(pod_u(:, i), nx1, ny1, n_elem);
            vi_c = reshape(pod_v(:, i), nx1, ny1, n_elem);
            ui_f = zeros(nx_f, nx_f, n_elem);
            vi_f = zeros(nx_f, nx_f, n_elem);
            for ie = 1:n_elem
                ui_f(:,:,ie) = Interp * ui_c(:,:,ie) * Interp';
                vi_f(:,:,ie) = Interp * vi_c(:,:,ie) * Interp';
            end
            pod_u_f(:, i) = ui_f(:);
            pod_v_f(:, i) = vi_f(:);
        end
        
        % Batch Gradients
        [ux_f_mat, uy_f_mat] = grad(reshape(pod_u_f, nx_f, nx_f, n_elem * nb), ...
                                    repmat(rx_f, [1,1,nb]), repmat(ry_f, [1,1,nb]), ...
                                    repmat(sx_f, [1,1,nb]), repmat(sy_f, [1,1,nb]), ...
                                    repmat(jaci_f, [1,1,nb]), d_f, 0);
        [vx_f_mat, vy_f_mat] = grad(reshape(pod_v_f, nx_f, nx_f, n_elem * nb), ...
                                    repmat(rx_f, [1,1,nb]), repmat(ry_f, [1,1,nb]), ...
                                    repmat(sx_f, [1,1,nb]), repmat(sy_f, [1,1,nb]), ...
                                    repmat(jaci_f, [1,1,nb]), d_f, 0);
        
        ux_f = reshape(ux_f_mat, nL_f_total, nb); 
        uy_f = reshape(uy_f_mat, nL_f_total, nb);
        vx_f = reshape(vx_f_mat, nL_f_total, nb); 
        vy_f = reshape(vy_f_mat, nL_f_total, nb);

        % --- Tensor Assembly ---
        tensor = zeros(n_test, nb, nb);
        
        if use_vectorized
            % VECTORIZED VERSION: O(n_test * nb) iterations
            % We pre-weight the test functions once
            Phi_K_weighted_u = bsxfun(@times, Me_f, pod_u_f(:, 2:nb));
            Phi_K_weighted_v = bsxfun(@times, Me_f, pod_v_f(:, 2:nb));
            
            for i = 1:nb
                % Convection operator for fixed 'i' acting on all 'j'
                % (u_i * d/dx + v_i * d/dy) applied to all modes j
                conv_j_x = bsxfun(@times, pod_u_f(:,i), ux_f) + ...
                           bsxfun(@times, pod_v_f(:,i), uy_f);
                conv_j_y = bsxfun(@times, pod_u_f(:,i), vx_f) + ...
                           bsxfun(@times, pod_v_f(:,i), vy_f);
                
                % Contract with all k test functions simultaneously
                % Result: [n_test, nb]
                res = (Phi_K_weighted_u' * conv_j_x) + (Phi_K_weighted_v' * conv_j_y);
                tensor(:, :, i) = res;
            end
        else
            % READABLE VERSION: O(n_test * nb^2) iterations
            for k = 1:n_test
                wk_u = Me_f .* pod_u_f(:, k+1);
                wk_v = Me_f .* pod_v_f(:, k+1);
                for i = 1:nb
                    ui = pod_u_f(:, i); vi = pod_v_f(:, i);
                    for j = 1:nb
                        conv_x = ui .* ux_f(:, j) + vi .* uy_f(:, j);
                        conv_y = ui .* vx_f(:, j) + vi .* vy_f(:, j);
                        tensor(k, j, i) = wk_u' * conv_x + wk_v' * conv_y;
                    end
                end
            end
        end
    end
    
    % --- Online Stage ---
    T_mat = reshape(tensor, n_test * nb, nb);
    temp = reshape(T_mat * ucoef, n_test, nb);
    out_coef = sum(temp .* ucoef', 2);
end
