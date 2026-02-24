function [out_coef] = conv_fom_2(ucoef, pod_u, pod_v, x, y, varargin)
    % Optimized Hybrid Pseudo-ROM / FOM Convection
    
    do_dealias = true;
    if nargin > 5, do_dealias = varargin{1}; end
    reset_flag = false;
    if nargin > 6, reset_flag = varargin{2}; end

    persistent pods_grad pod_test_weighted Jh_1d Jh_1dT nx1 nx_fine M_fine_weights nModes_old use_pagemtimes

    curr_nx1 = size(x, 1);
    curr_nModes = size(pod_u, 2);
    
    % --- OFFLINE PHASE ---
    if isempty(pod_test_weighted) || curr_nx1 ~= nx1 || curr_nModes ~= nModes_old || reset_flag
        nx1 = curr_nx1;
        nModes_old = curr_nModes;
        use_pagemtimes = (exist('pagemtimes', 'builtin') == 5); % 5 = builtin

        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [~,~,~,~,rx,ry,sx,sy,jac,jaci,~] = deriv_geo(x,y,d);
        
        % 1. POD Derivatives
        pods_grad.ux = zeros(numel(x), curr_nModes); pods_grad.uy = pods_grad.ux;
        pods_grad.vx = pods_grad.ux; pods_grad.vy = pods_grad.ux;
        
        for i = 1:curr_nModes
            ui_3d = reshape(pod_u(:,i), nx1, nx1, []);
            vi_3d = reshape(pod_v(:,i), nx1, nx1, []);
            [ux, uy] = grad(ui_3d, rx, ry, sx, sy, jaci, d, 0);
            [vx, vy] = grad(vi_3d, rx, ry, sx, sy, jaci, d, 0);
            pods_grad.ux(:,i) = ux(:); pods_grad.uy(:,i) = uy(:);
            pods_grad.vx(:,i) = vx(:); pods_grad.vy(:,i) = vy(:);
        end

        % 2. De-aliasing Setup (3/2 rule)
        nx_fine = ceil(nx1 * 1.5);
        [zi_fine, w_fine] = zwgll(nx_fine-1);
        Jh_1d = interp_mat(zi_fine, zi); Jh_1dT = Jh_1d';
        
        % Build weighted fine-grid Jacobian
        jac_3d = reshape(jac, nx1, nx1, []);
        if use_pagemtimes
            jac_fine = pagemtimes(pagemtimes(Jh_1d, jac_3d), Jh_1dT);
        else
            ne = size(jac_3d, 3);
            jac_fine = zeros(nx_fine, nx_fine, ne);
            for ie = 1:ne
                jac_fine(:,:,ie) = Jh_1d * jac_3d(:,:,ie) * Jh_1dT;
            end
        end
        Wf_2d = w_fine * w_fine';
        M_fine_weights = jac_fine .* reshape(Wf_2d, nx_fine, nx_fine, 1);
        
        % 3. Pre-weight POD Test Functions (Modes 2:N) for Online Projection
        % We include the coarse mass weights Me here.
        Me_vec = reshape(jac_3d .* reshape(w * w', nx1, nx1), [], 1);
        pod_test_weighted = [bsxfun(@times, pod_u(:, 2:end), Me_vec); ...
                             bsxfun(@times, pod_v(:, 2:end), Me_vec)]';
                             
        if reset_flag, return; end
    end

    % --- ONLINE PHASE ---
    % Reconstruct FOM fields (Matrix-Vector multiply is very fast)
    u_f  = reshape(pod_u * ucoef, nx1, nx1, []);
    v_f  = reshape(pod_v * ucoef, nx1, nx1, []);
    ux_f = reshape(pods_grad.ux * ucoef, nx1, nx1, []);
    uy_f = reshape(pods_grad.uy * ucoef, nx1, nx1, []);
    vx_f = reshape(pods_grad.vx * ucoef, nx1, nx1, []);
    vy_f = reshape(pods_grad.vy * ucoef, nx1, nx1, []);

    if do_dealias
        if use_pagemtimes
            u_h = pagemtimes(pagemtimes(Jh_1d, u_f), Jh_1dT);
            v_h = pagemtimes(pagemtimes(Jh_1d, v_f), Jh_1dT);
            ux_h = pagemtimes(pagemtimes(Jh_1d, ux_f), Jh_1dT);
            uy_h = pagemtimes(pagemtimes(Jh_1d, uy_f), Jh_1dT);
            vx_h = pagemtimes(pagemtimes(Jh_1d, vx_f), Jh_1dT);
            vy_h = pagemtimes(pagemtimes(Jh_1d, vy_f), Jh_1dT);
            
            % Quadratic interaction + Fine Weights (Weak Form)
            cu_h = (u_h .* ux_h + v_h .* uy_h) .* M_fine_weights;
            cv_h = (u_h .* vx_h + v_h .* vy_h) .* M_fine_weights;
            
            % Project back to coarse grid (Integration)
            c_u = pagemtimes(pagemtimes(Jh_1dT, cu_h), Jh_1d);
            c_v = pagemtimes(pagemtimes(Jh_1dT, cv_h), Jh_1d);
        else
            ne = size(u_f, 3);
            c_u = zeros(nx1, nx1, ne); c_v = c_u;
            for ie = 1:ne
                % Process and project in one loop to save memory
                uh_ie = Jh_1d * u_f(:,:,ie) * Jh_1dT;
                vh_ie = Jh_1d * v_f(:,:,ie) * Jh_1dT;
                
                cu_h = (uh_ie .* (Jh_1d * ux_f(:,:,ie) * Jh_1dT) + ...
                        vh_ie .* (Jh_1d * uy_f(:,:,ie) * Jh_1dT)) .* M_fine_weights(:,:,ie);
                cv_h = (uh_ie .* (Jh_1d * vx_f(:,:,ie) * Jh_1dT) + ...
                        vh_ie .* (Jh_1d * vy_f(:,:,ie) * Jh_1dT)) .* M_fine_weights(:,:,ie);
                
                c_u(:,:,ie) = Jh_1dT * cu_h * Jh_1d;
                c_v(:,:,ie) = Jh_1dT * cv_h * Jh_1d;
            end
        end
        conv_vec = [c_u(:); c_v(:)];
    else
        % No de-aliasing: just multiply by coarse weights Me
        % Note: pod_test_weighted already has weights, so we don't multiply here
        conv_vec = [(u_f(:) .* ux_f(:) + v_f(:) .* uy_f(:)); ...
                    (u_f(:) .* vx_f(:) + v_f(:) .* vy_f(:))];
        % Apply weights if not using the de-alias path
        Me_vec = pod_test_weighted(1,1) / pod_u(1,2); % (Pseudo-code logic: extract Me)
        % Simplified: logic below is cleaner
    end

    % Final Projection onto POD basis
    out_coef = pod_test_weighted * conv_vec;
end
