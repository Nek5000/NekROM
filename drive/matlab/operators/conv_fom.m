function [out_coef] = conv_fom(ucoef, pod_u, pod_v, x, y, varargin)
    % Hybrid Multi-element Pseudo-ROM (MATLAB/Octave)
    
    do_dealias = true;
    if nargin > 5, do_dealias = varargin{1}; end
    reset_flag = false;
    if nargin > 6, reset_flag = varargin{2}; end

    persistent pods_grad pod_uv_weighted Jh_1d Jh_1dT nx1 nx_fine M_fine_weights invMc nModes_old use_pagemtimes

    curr_nx1 = size(x, 1);
    curr_nModes = size(pod_u, 2);
    curr_nL = numel(x);
    
    if isempty(pod_uv_weighted) || curr_nx1 ~= nx1 || curr_nModes ~= nModes_old || reset_flag
        nx1 = curr_nx1;
        nModes_old = curr_nModes;
        use_pagemtimes = exist('pagemtimes', 'builtin'); % Check environment once

        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [~,~,~,~,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        
        % 1. Pre-calculate POD derivatives
        pods_grad.ux = zeros(curr_nL, curr_nModes); pods_grad.uy = zeros(curr_nL, curr_nModes);
        pods_grad.vx = zeros(curr_nL, curr_nModes); pods_grad.vy = zeros(curr_nL, curr_nModes);
        
        for i = 1:curr_nModes
            u_3d = reshape(pod_u(:,i), nx1, nx1, []);
            v_3d = reshape(pod_v(:,i), nx1, nx1, []);
            [ux, uy] = grad(u_3d, rx, ry, sx, sy, jaci, d, 0);
            [vx, vy] = grad(v_3d, rx, ry, sx, sy, jaci, d, 0);
            pods_grad.ux(:,i) = ux(:); pods_grad.uy(:,i) = uy(:);
            pods_grad.vx(:,i) = vx(:); pods_grad.vy(:,i) = vy(:);
        end

        % 2. Setup De-aliasing
        nx_fine = ceil(nx1 * 1.5);
        [zi_fine, w_fine] = zwgll(nx_fine-1);
        Jh_1d = interp_mat(zi_fine, zi); 
        Jh_1dT = Jh_1d';
        
        % 3. Build Fine Mass Matrix
        jac_3d = reshape(jac, nx1, nx1, []);
        if use_pagemtimes
            jac_fine = pagemtimes(pagemtimes(Jh_1d, jac_3d), Jh_1dT);
        else
            jac_fine = zeros(nx_fine, nx_fine, size(jac_3d, 3));
            for ie = 1:size(jac_3d, 3)
                jac_fine(:,:,ie) = Jh_1d * jac_3d(:,:,ie) * Jh_1dT;
            end
        end
        
        Wf_2d = w_fine * w_fine';
        M_fine_weights = jac_fine .* reshape(Wf_2d, nx_fine, nx_fine, 1);
        
        % 4. Coarse Weights
        Me_3d = jac_3d .* reshape(w * w', nx1, nx1, 1);
        invMc = 1 ./ Me_3d;
        pod_uv_weighted = [bsxfun(@times, pod_u(:, 2:end), Me_3d(:)); ...
                           bsxfun(@times, pod_v(:, 2:end), Me_3d(:))]';
        if reset_flag, return; end
    end

    % --- ONLINE PHASE ---
    
    % Reconstruct FOM fields
    u_f = reshape(pod_u * ucoef, nx1, nx1, []);
    v_f = reshape(pod_v * ucoef, nx1, nx1, []);
    ux_f = reshape(pods_grad.ux * ucoef, nx1, nx1, []);
    uy_f = reshape(pods_grad.uy * ucoef, nx1, nx1, []);
    vx_f = reshape(pods_grad.vx * ucoef, nx1, nx1, []);
    vy_f = reshape(pods_grad.vy * ucoef, nx1, nx1, []);

    if do_dealias
        if use_pagemtimes
            % Modern Path
            u_h = pagemtimes(pagemtimes(Jh_1d, u_f), Jh_1dT);
            v_h = pagemtimes(pagemtimes(Jh_1d, v_f), Jh_1dT);
            ux_h = pagemtimes(pagemtimes(Jh_1d, ux_f), Jh_1dT);
            uy_h = pagemtimes(pagemtimes(Jh_1d, uy_f), Jh_1dT);
            vx_h = pagemtimes(pagemtimes(Jh_1d, vx_f), Jh_1dT);
            vy_h = pagemtimes(pagemtimes(Jh_1d, vy_f), Jh_1dT);
        else
            % Fallback Path
            ne = size(u_f, 3);
            u_h = zeros(nx_fine, nx_fine, ne); v_h = u_h; ux_h = u_h; uy_h = u_h; vx_h = u_h; vy_h = u_h;
            for ie = 1:ne
                u_h(:,:,ie) = Jh_1d * u_f(:,:,ie) * Jh_1dT;
                v_h(:,:,ie) = Jh_1d * v_f(:,:,ie) * Jh_1dT;
                ux_h(:,:,ie) = Jh_1d * ux_f(:,:,ie) * Jh_1dT;
                uy_h(:,:,ie) = Jh_1d * uy_f(:,:,ie) * Jh_1dT;
                vx_h(:,:,ie) = Jh_1d * vx_f(:,:,ie) * Jh_1dT;
                vy_h(:,:,ie) = Jh_1d * vy_f(:,:,ie) * Jh_1dT;
            end
        end
        
        % Quadratic interactions
        cu_h = u_h .* ux_h + v_h .* uy_h;
        cv_h = u_h .* vx_h + v_h .* vy_h;
        
        if use_pagemtimes
            % Modern Path Projection
            c_u = invMc .* pagemtimes(pagemtimes(Jh_1dT, M_fine_weights .* cu_h), Jh_1d);
            c_v = invMc .* pagemtimes(pagemtimes(Jh_1dT, M_fine_weights .* cv_h), Jh_1d);
        else
            % Fallback Path Projection
            c_u = zeros(nx1, nx1, ne); c_v = c_u;
            for ie = 1:ne
                c_u(:,:,ie) = invMc(:,:,ie) .* (Jh_1dT * (M_fine_weights(:,:,ie) .* cu_h(:,:,ie)) * Jh_1d);
                c_v(:,:,ie) = invMc(:,:,ie) .* (Jh_1dT * (M_fine_weights(:,:,ie) .* cv_h(:,:,ie)) * Jh_1d);
            end
        end
        conv_fom_vec = [c_u(:); c_v(:)];
    else
        conv_fom_vec = [(u_f(:) .* ux_f(:) + v_f(:) .* uy_f(:)); ...
                        (u_f(:) .* vx_f(:) + v_f(:) .* vy_f(:))];
    end

    out_coef = pod_uv_weighted * conv_fom_vec;
end
