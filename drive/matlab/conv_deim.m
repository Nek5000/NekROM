function [out_coef] = conv_deim(ucoef, pod_u, pod_v, nl_snaps_obj, ndeim_pts,istep,clsdeim,n_os_points)
    persistent proj_mat Ainv inv_p_nl u_deimu v_deimu u_deimv v_deimv ux_deimu uy_deimu vx_deimv vy_deimv;
    persistent u_deim_stack v_deim_stack ux_deim_stack uy_deim_stack tau mu A_tau_inv alpha nl_bas_inds% nl_max_coef nl_min_coef;
    persistent inds;
  
%   if isempty(proj_mat)
    if istep == 1
        % Stuff to be precomputed

        % Read in the snapshots.     
        [nl_u_snaps, nl_v_snaps] = get_grid_and_pod(nl_snaps_obj);
        nl_snaps = [nl_u_snaps; nl_v_snaps];

        x=nl_snaps_obj.flds{1}.x;
        y=nl_snaps_obj.flds{1}.y;
        nx1 = size(x,1);
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
        nL = prod(size(x));
        Me = reshape(jac.*(w*w'),nL,1);

        %% Calculate the POD of the NL snapshots
        gramian = nl_snaps'*diag([Me;Me])*nl_snaps;
        gramian = 0.5*(gramian + gramian');
        [nl_eigvecs, nl_eigvals] = eig(gramian);
        [nl_eigvals_sorted, sort_inds] = sort(diag(abs(nl_eigvals)),'descend');
        fileID = fopen('nl_eigvals_sorted.txt', 'w');
        fprintf(fileID, '%24.15e\n', nl_eigvals_sorted);
        fclose(fileID);
        nl_eigvecs = nl_eigvecs(:,sort_inds);
        nl_bas = nl_snaps*nl_eigvecs; % Using all of the modes 
        nl_bas = nl_bas(:,1:ndeim_pts);
        % Maybe 500 snapshots is not enough?
        nl_bas = orth(nl_bas,1e-16); % Use all of the snapshots for now.

        % For use with Constrained DEIM
        %nl_snapshot_proj = nl_bas'*nl_snaps;
        %nl_max_coef = max(nl_snapshot_proj,[],2);
        %nl_min_coef = min(nl_snapshot_proj,[],2); 


        nl_bas_u = nl_bas(1:nL,:);
        nl_bas_v = nl_bas(nL+1:end,:);
        nl_bas_u = orth(nl_bas_u,1e-16); % Use all of the snapshots for now.
        nl_bas_v = orth(nl_bas_v,1e-16); % Use all of the snapshots for now.
        %nl_bas = orth(nl_bas(:,1:500),1e-16);
        %for i = 1:size(nl_bas,2)
        %    nl_bas(:,i) = nl_bas(:,i)/norm(nl_bas(:,i));
        %end

        separate=false

        if separate % Calculate the QDEIM points separately for u and v
            % This probably won't work.
%           [P_u, nl_u_inds] = calc_qdeim_proj_mat(nl_u_snaps);
%           [P_v, nl_v_inds] = calc_qdeim_proj_mat(nl_v_snaps);
            [P_u, nl_u_inds] = calc_qdeim_proj_mat(nl_bas_u);
            [P_v, nl_v_inds] = calc_qdeim_proj_mat(nl_bas_v); % not sure why the separate approach is not working
            % maybe because this will not satisfiyes divergence-free constraints...
            n_deim_points_u = ceil(size(nl_bas,2)/2);
            n_deim_points_v = floor(size(nl_bas,2)/2);
            nl_u_inds = nl_u_inds(1,1:n_deim_points_u); 
            nl_v_inds = nl_v_inds(1,1:n_deim_points_v); 

            inds = [nl_u_inds, nl_v_inds + size(nl_u_snaps,1)];
        else % Use the same QDEIM points for u and v.
            % Maybe this isn't right for vector quantities.
%           [P, inds] = calc_qdeim_proj_mat(nl_snaps);
            if false;
               [P, inds] = calc_qdeim_proj_mat(nl_bas); % Should select points based on the basis
            else;
            % Can oversample if desired.
               inds = s_opt_generator(nl_bas, ndeim_pts + n_os_points, [])
               inds = inds';
            end;
            % For testing, use all rows
            %inds = [1:size(nl_snaps,1)];
            %divider = size(nl_u_snaps,1);
            %nl_u_inds = inds(inds <= divider);
            %nl_v_inds = inds(inds > divider);
            %nl_u_inds
            %exit
            %inds = [nl_u_inds, nl_v_inds];
        end;
        inds = inds(:,1:ndeim_pts);

       
        % Kento's ROM approach. Calculate the gradients of the POD modes
        ux_pods = [];
        uy_pods = [];
        vx_pods = [];
        vy_pods = [];
        for i = 1:size(pod_u,2);
            [ux_pod, uy_pod] = lgrad(reshape(pod_u(:,i),size(x)),0);
            [vx_pod, vy_pod] = lgrad(reshape(pod_v(:,i),size(x)),0);
            ux_pods = [ux_pods, reshape(ux_pod, nL,1)];
            uy_pods = [uy_pods, reshape(uy_pod, nL,1)];
            vx_pods = [vx_pods, reshape(vx_pod, nL,1)];
            vy_pods = [vy_pods, reshape(vy_pod, nL,1)];
        end;

      
        % Is this problematic? What does it mean to integrate on the DEIM points?
        u_deim = Me.*pod_u;
        v_deim = Me.*pod_v;
%       u_deim = pod_u;
%       v_deim = pod_v;
    
        if separate
            u_deimu = u_deim(nl_u_inds,:);
            v_deimu = v_deim(nl_v_inds,:);
            u_deimv = u_deim(nl_u_inds,:);
            v_deimv = v_deim(nl_v_inds,:);

            ux_deimu = ux_pods(nl_u_inds,:);
            uy_deimu = uy_pods(nl_u_inds,:);
            vx_deimv = vx_pods(nl_v_inds,:);
            vy_deimv = vy_pods(nl_v_inds,:);
%           proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*[nl_bas_u; nl_bas_v]*inv([nl_bas_u(nl_u_inds,:); nl_bas_v(nl_v_inds,:)]); 
        else
            u_deim_stack = [u_deim; u_deim];
            u_deim_stack = u_deim_stack(inds,:);
            v_deim_stack = [v_deim; v_deim];
            v_deim_stack = v_deim_stack(inds,:);


            ux_deim_stack = [ux_pods; vx_pods];
            ux_deim_stack = ux_deim_stack(inds,:);
            uy_deim_stack = [uy_pods; vy_pods];
            uy_deim_stack = uy_deim_stack(inds,:);
%           proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*nl_bas*inv(nl_bas(inds,:)); 
        end;

        proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*nl_bas;

        nl_bas_inds = nl_bas(inds,:);
        if size(nl_bas,2) == size(inds,1);   
            inv_p_nl = inv(nl_bas(inds,:));
        else;
            inv_p_nl = pinv(nl_bas(inds,:));
        end;
         
%       proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*nl_bas*inv(nl_bas);
%       proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*([Me;Me].*nl_bas)*inv(nl_bas(inds,:)); 
        % For testing
        %proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]';
  
        % Matrices for CLSDEIM
        % Doesn't seem like the 2 should be necessary
        %Ainv = inv(2*nl_bas(inds,:)'*nl_bas(inds,:));
        Ainv = inv(nl_bas(inds,:)'*nl_bas(inds,:));

        % Matrices for MCLSDEIM
        nl_snapshot_proj = nl_bas'*nl_snaps;
        tau = inv(cov(nl_snapshot_proj'));
        size(nl_snapshot_proj)
        size(tau)
        size(nl_bas_inds)
        mu = mean(nl_snapshot_proj,2);
        alpha = 1e-14;
        A_tau_inv = inv(nl_bas(inds,:)'*nl_bas(inds,:) + alpha*tau); 
    end;

    separate=false;
    mclsdeim = false;

    if separate
        conv_u_deim = ((u_deimu(:,2:end)*ucoef(2:end)).*(ux_deimu(:,2:end)*ucoef(2:end)) + (v_deimu(:,2:end)*ucoef(2:end)).*(uy_deimu(:,2:end)*ucoef(2:end)));
        conv_v_deim = ((u_deimv(:,2:end)*ucoef(2:end)).*(vx_deimv(:,2:end)*ucoef(2:end)) + (v_deimv(:,2:end)*ucoef(2:end)).*(vy_deimv(:,2:end)*ucoef(2:end))); 
        out_coef = proj_mat*inv_p_nl*[conv_u_deim; conv_v_deim];
    else
        conv_deim = ((u_deim_stack(:,2:end)*ucoef(2:end)).*(ux_deim_stack(:,2:end)*ucoef(2:end)) + (v_deim_stack(:,2:end)*ucoef(2:end)).*(uy_deim_stack(:,2:end)*ucoef(2:end)));
        out_coef = inv_p_nl*conv_deim;
        if clsdeim;
            b = proj_mat'*ucoef(2:end);
            out_coef = out_coef - ((b'*out_coef)/(b'*Ainv*b))*(Ainv*b);
        elseif mclsdeim;
            %disp('here');
            b = proj_mat'*ucoef(2:end);
            out_coef = A_tau_inv*(nl_bas_inds'*conv_deim + alpha*tau*mu); 
            out_coef = out_coef - ((b'*out_coef)/(b'*A_tau_inv*b))*(A_tau_inv*b);

        %    This does not work. It keeps the problem from blowing up for longer, but it still eventually blows up.
        %    disp('HERE');
        %    options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
        %    cmin_func = @(x) inv_p_nl@x;
        %    [out_coef,fval,exitflag,output] = fmincon(cmin_func,conv_deim,[],[],[],[],nl_max_coef,nl_min_coef,[],options);
        end;
        out_coef = proj_mat*out_coef;
    end;
end

