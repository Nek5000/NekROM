function [out_coef] = conv_tensor(ucoef, pod_u, pod_v, x, y)

    % Compute the full convection tensor
    % Mostly replaced by conv_tensor_dense

    persistent Me rx ry sx sy jaci d lgrad nL nb tensor
    %persistent tensor

    if isempty(tensor)
        nx1 = size(x,1);
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
        nL = prod(size(x));
        nb = size(pod_u,2)
        Me = reshape(jac.*(w*w'),nL,1);

    if false; % Normal way of calculating the gradient
        [ux_fom, uy_fom] = lgrad(u_fom, 0);
        [vx_fom, vy_fom] = lgrad(v_fom, 0);
    else
        % Kento's ROM approach. Calculate the gradients of the POD modes
        ux_pods = zeros(size(pod_u));
        uy_pods = zeros(size(pod_u));
        vx_pods = zeros(size(pod_u));
        vy_pods = zeros(size(pod_u));
        
        % For Pseudo-FOM version
        %{
        ux_fom = zeros(size(x));
        uy_fom = zeros(size(x));
        vx_fom = zeros(size(x));
        vy_fom = zeros(size(x));
        %}
        for i = 1:size(pod_u,2);
            [ux_pod, uy_pod] = lgrad(reshape(pod_u(:,i),size(x)),0);
            [vx_pod, vy_pod] = lgrad(reshape(pod_v(:,i),size(x)),0);
           
            %{
            % For Pseudo-FOM version 
            ux_fom = ux_fom + ux_pod*ucoef(i);
            uy_fom = uy_fom + uy_pod*ucoef(i);
            vx_fom = vx_fom + vx_pod*ucoef(i);
            vy_fom = vy_fom + vy_pod*ucoef(i);
            %}

            ux_pods(:,i) = reshape(ux_pod, nL,1);
            uy_pods(:,i) = reshape(uy_pod, nL,1);
            vx_pods(:,i) = reshape(vx_pod, nL,1);
            vy_pods(:,i) = reshape(vy_pod, nL,1);
        end;
        %ux_fom = reshape(ux_pods*ucoef, size(x));
        %uy_fom = reshape(uy_pods*ucoef, size(x));
        %vx_fom = reshape(vx_pods*ucoef, size(x));
        %vy_fom = reshape(vy_pods*ucoef, size(x));

        %pod_u_weak = Me.*pod_u;
        %pod_v_weak = Me.*pod_v;

        % Could also apply Me to pod_u and pod_v inside the loop instead
        % but this seems more efficient.
        pod_weak = [Me.*pod_u(:,2:nb);Me.*pod_v(:,2:nb)];
        tensor = zeros(nb,nb,nb-1);
        for i=1:nb;
            %pod_u_weak = Me.*pod_u(:,i);
            %pod_v_weak = Me.*pod_v(:,i); 
            %tensor(i,:,:) = [pod_u_weak.*ux_pods + pod_v_weak.*uy_pods;
            %                 pod_u_weak.*vx_pods + pod_v_weak.*vy_pods]'*pod(:,2:nb);
            tensor(i,:,:) = [pod_u(:,i).*ux_pods + pod_v(:,i).*uy_pods;
                             pod_u(:,i).*vx_pods + pod_v(:,i).*vy_pods]'*pod_weak;

            %{
            % More readable version
            for j=1:nb;
                % Calculate pencil
                [i,j]
                tensor(i,j,:) = pod_weak'*[pod_u(:,i).*ux_pods(:,j) + pod_v(:,i).*uy_pods(:,j); 
                                              pod_u(:,i).*vx_pods(:,j) + pod_v(:,i).*vy_pods(:,j)];
                % Reduce pencil
                %ijk = pod(:,2:nb)'*ij;
                %tensor(i,j,:) = ijk;
            end;
            %}
        end
    end
        
    end; 
    % End of if block to persist the tensor. Note, if comparing against the conv_fom approach
    % this will need to be moved back up

    %{
    % For Pseudo-FOM version
    u_fom = reshape(Me.*(pod_u*ucoef), size(x));
    v_fom = reshape(Me.*(pod_v*ucoef), size(x));
    conv_u_fom = reshape(u_fom.*ux_fom + v_fom.*uy_fom, nL,1);
    conv_v_fom = reshape(u_fom.*vx_fom + v_fom.*vy_fom, nL,1);    
    out_coef = [pod_u(:,2:end); pod_v(:,2:end)]'*[conv_u_fom; conv_v_fom]
    %}     

    outprod = tensorprod(tensor, ucoef, 1,1);
    out_coef = tensorprod(outprod,ucoef,1,1);
    %out_coef
    %exit;
end
