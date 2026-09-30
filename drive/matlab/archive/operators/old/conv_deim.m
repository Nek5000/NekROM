function [out_coef] = conv_deim(ucoef, pod_u, pod_v, nl_bas, nl_snaps_u, nl_snaps_v, x, y, ndeim_pts,istep,clsdeim,n_os_points,ps_alg)

    % Convection operator that uses DEIM points
    % TODO: Separate DEIM and clsdeim stuff. Just have the CLS DEIM function call the DEIM function

    persistent proj_mat Ainv interp_mat u_deimu v_deimu u_deimv v_deimv ux_deimu uy_deimu vx_deimv vy_deimv nb;
    persistent u_deim_stack v_deim_stack ux_deim_stack uy_deim_stack tau mu A_tau_inv alpha nl_bas_inds% nl_max_coef nl_min_coef;
    persistent inds proj_and_interp_mat zeroth_mode_contribution validator
  
%   if isempty(proj_mat)
    if istep == 1
        disp("STARTING NL POD");
        % Stuff to be precomputed

        % Read in the snapshots.     
        %[nl_u_snaps, nl_v_snaps] = get_snaps(nl_snaps_obj);
        %nl_snaps = [nl_u_snaps; nl_v_snaps];

        % Calculate mass matrix with Jacobian
        %x=nl_snaps_obj.flds{1}.x;
        %y=nl_snaps_obj.flds{1}.y;

        % Should probably make all of these functions persistent
        nx1 = size(x,1);
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        my_lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
        nL = prod(size(x));
        nb = size(pod_u,2);

        validator = @(i,j,k) i == 1 || j == 1;

        Me = reshape(jac.*(w*w'),nL,1);
        %if 0;
            %[nl_bas, ~, ~] = get_pod_basis_from_arrays(nl_snaps_u, nl_snaps_v, x, y, ndeim_pts, 0, 0);
            %size(nl_bas)
            %size(nl_bas_nr)
            %norm(nl_bas_nr(:,1:31) - nl_bas)/norm(nl_bas)
            %exit;
        %else;
            %nl_bas = nl_bas_nr;
            %Me = Me_in;
        %end;
        %nl_bas = nl_bas_nr
        %{
        size(nl_bas)
        size(nl_bas_nr)
        nl_bas'*([Me;Me].*nl_bas)
        nl_bas_nr'*([Me;Me].*nl_bas_nr)
        
        mean(abs(nl_bas),1)
        mean(abs(nl_bas_nr),1)
        %exit;
        %}
        % For use with Constrained DEIM
        %nl_snapshot_proj = nl_bas'*nl_snaps;
        %nl_max_coef = max(nl_snapshot_proj,[],2);
        %nl_min_coef = min(nl_snapshot_proj,[],2); 
    

        %nl_bas_u = nl_bas(1:nL,:);
        %nl_bas_v = nl_bas(nL+1:end,:);
        %nl_bas_u = orth(nl_bas_u,1e-16); % Use all of the snapshots for now.
        %nl_bas_v = orth(nl_bas_v,1e-16); % Use all of the snapshots for now.
        %nl_bas = orth(nl_bas(:,1:500),1e-16);
        %for i = 1:size(nl_bas,2)
        %    nl_bas(:,i) = nl_bas(:,i)/norm(nl_bas(:,i));
        %end

        separate=false

        if separate % Calculate the QDEIM points separately for u and v
        % Delete this
        exit;
        else % Use the same QDEIM points for u and v.
            % Maybe this isn't right for vector quantities.
            % Will this satisfy the divergence-free constraints?
            % Yes, because it is projected onto the divergence-free subspace
            % of the snapshots
%           [P, inds] = calc_qdeim_proj_mat(nl_snaps);
            if strcmp(ps_alg, 'sopt')
                % Can oversample if desired.
                inds = s_opt(nl_bas, ndeim_pts + n_os_points, []);
                inds = inds';
            elseif strcmp(ps_alg, 'gpode') || strcmp(ps_alg, 'qdeim')
                inds = gpode(nl_bas, ndeim_pts + n_os_points);
                inds = inds';
            elseif strcmp(ps_alg, 'gappy_pod') || strcmp(ps_alg, 'deim');
                inds = gappy_pod(nl_bas, ndeim_pts + n_os_points);
            elseif strcmp(ps_alg, 'gnat')
                inds = gnat(nl_bas, ndeim_pts, ndeim_pts + n_os_points);
            else
                throw(MException('Unknown point selection algorithm %s', ps_alg));
            end;

            if false;
               [P, inds] = calc_qdeim_proj_mat(nl_bas); % Should select points based on the basis
            else;
            end;
            % For testing, use all rows
            %inds = [1:size(nl_snaps,1)];
            %divider = size(nl_u_snaps,1);
            %nl_u_inds = inds(inds <= divider);
            %nl_v_inds = inds(inds > divider);
            %nl_u_inds
            %exit
            %inds = [nl_u_inds, nl_v_inds];
            disp("ENDING NL POD");
        end;
        inds = inds(:,1:ndeim_pts);

       
        % Kento's ROM approach. Calculate the gradients of the POD modes
        ux_pods = [];
        uy_pods = [];
        vx_pods = [];
        vy_pods = [];
        for i = 1:nb;
            [ux_pod, uy_pod] = my_lgrad(reshape(pod_u(:,i),size(x)),0);
            [vx_pod, vy_pod] = my_lgrad(reshape(pod_v(:,i),size(x)),0);
            ux_pods = [ux_pods, reshape(ux_pod, nL,1)];
            uy_pods = [uy_pods, reshape(uy_pod, nL,1)];
            vx_pods = [vx_pods, reshape(vx_pod, nL,1)];
            vy_pods = [vy_pods, reshape(vy_pod, nL,1)];
        end;

      
        % Is this problematic? What does it mean to integrate on the DEIM points?
        u_deim = pod_u;
        v_deim = pod_v;
%       u_deim = pod_u;
%       v_deim = pod_v;
    
        if separate
        % Delete this
        exit;
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

        Me_pod = [Me.*pod_u(:,2:end); Me.*pod_v(:,2:end)]; 
        proj_mat = Me_pod'*nl_bas;
        
        % This is equivalent to c1
        %c1 = Me_pod'*([pod_u(:,1);pod_u(:,1)].*[ux_pods(:,1);vx_pods(:,1)] + ...
        %                                    [pod_v(:,1);pod_v(:,1)].*[uy_pods(:,1);vy_pods(:,1)]);

        c2 = Me_pod'*([pod_u;pod_u].*[ux_pods(:,1);vx_pods(:,1)] + ...
                      [pod_v;pod_v].*[uy_pods(:,1);vy_pods(:,1)]);

        c3 = Me_pod'*([pod_u(:,1);pod_u(:,1)].*[ux_pods;vx_pods] + ...
                      [pod_v(:,1);pod_v(:,1)].*[uy_pods;vy_pods]);


        zeroth_mode_contribution = c2 + c3;
        % The first row was counted twice.
        zeroth_mode_contribution(:,1) = zeroth_mode_contribution(:,1)/2;
        
        
        %zeroth_mode_contribution
        nl_bas_inds = nl_bas(inds,:);
        if size(nl_bas,2) == size(inds,1);   
            interp_mat = inv(nl_bas(inds,:));
        else;
            interp_mat = pinv(nl_bas(inds,:));
        end;
       
        proj_and_interp_mat = proj_mat*interp_mat;
  
%       proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*nl_bas*inv(nl_bas);
%       proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*([Me;Me].*nl_bas)*inv(nl_bas(inds,:)); 
        % For testing
        %proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]';

        % Matrices for CLSDEIM
        % Doesn't seem like the 2 should be necessary
        %Ainv = inv(2*nl_bas(inds,:)'*nl_bas(inds,:));
        if clsdeim
            Ainv = inv(nl_bas(inds,:)'*nl_bas(inds,:));

            % Matrices for MCLSDEIM
            nl_snapshot_proj = nl_bas'*[nl_snaps_u; nl_snaps_v];
            tau = inv(cov(nl_snapshot_proj'));
            size(nl_snapshot_proj)
            size(tau)
            size(nl_bas_inds)
            mu = mean(nl_snapshot_proj,2);
            alpha = 1e-14;
            A_tau_inv = inv(nl_bas(inds,:)'*nl_bas(inds,:) + alpha*tau); 
        end;
         
    end;

    mclsdeim = false;
    separate=false;

    if separate
    % Delete this
    exit;
    else
        % Note that this does not seem to work. Need to explicitly include the zeroth
        % mode interactions apparently
        %conv_deim = ((u_deim_stack(:,1:end)*ucoef(1:end)).*(ux_deim_stack(:,1:end)*ucoef(1:end)) + ... 
        %             (v_deim_stack(:,1:end)*ucoef(1:end)).*(uy_deim_stack(:,1:end)*ucoef(1:end)));

        conv_deim = ((u_deim_stack(:,2:end)*ucoef(2:end)).*(ux_deim_stack(:,2:end)*ucoef(2:end)) + ...
                     (v_deim_stack(:,2:end)*ucoef(2:end)).*(uy_deim_stack(:,2:end)*ucoef(2:end)));

        %out_coef = interp_mat*conv_deim;

        if clsdeim;
            out_coef = interp_mat*conv_deim;
            b = proj_mat'*ucoef(2:end);
            out_coef = out_coef - ((b'*out_coef)/(b'*Ainv*b))*(Ainv*b);
            out_coef = proj_mat*out_coef;
        elseif mclsdeim;
            %disp('here');

            out_coef = interp_mat*conv_deim;
            b = proj_mat'*ucoef(2:end);
            out_coef = A_tau_inv*(nl_bas_inds'*conv_deim + alpha*tau*mu); 
            out_coef = out_coef - ((b'*out_coef)/(b'*A_tau_inv*b))*(A_tau_inv*b);
            out_coef = proj_mat*out_coef;
    
        %    This does not work. It keeps the problem from blowing up for longer, but it still eventually blows up.
        %    disp('HERE');
        %    options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
        %    cmin_func = @(x) inv_p_nl@x;
        %    [out_coef,fval,exitflag,output] = fmincon(cmin_func,conv_deim,[],[],[],[],nl_max_coef,nl_min_coef,[],options);
        else 
            % Standard DEIM
            out_coef = proj_and_interp_mat*conv_deim;
        end;

        %zeroth_mode_contribution = conv_tensor_reduced(ucoef, pod_u, pod_v, x, y, [1, nb, nb-1])
        out_coef = out_coef + zeroth_mode_contribution*ucoef;%c2*ucoef(1:end) + c3*ucoef(1:end) - c1;
    
        % This does the same thing as the above, but the above is likely faster
        %out_coef = out_coef + conv_tensor_sparse(ucoef, pod_u, pod_v, x, y, validator);

    end;
end

