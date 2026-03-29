% Pseudo-ROM convection operator
% (pseudo because the work still scales
% with the size of the original problem)
% This needs to do the same thing as
% reshape(cu*utmp(:,1),nb,nb+1)*u(:,1);
%
% Note: Dealiasing is not currently implemented. Is it needed?
function [out_coef] = conv_fom(ucoef, pod_u, pod_v, x, y)

    %persistent Me rx ry sx sy jaci d lgrad nL
    persistent ux_pods uy_pods vx_pods vy_pods pod_uv

    if isempty(pod_uv)
        nx1 = size(x,1);
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
        Me = reshape(jac.*(w*w'),nL,1);

        nL = prod(size(x));
        % Assume the geometry is not moving
        % and pre-calcuate the pod derivatives
        ux_pods = zeros(size(pod_u));
        uy_pods = zeros(size(pod_u));
        vx_pods = zeros(size(pod_u));
        vy_pods = zeros(size(pod_u));
        for i = 1:size(pod_u,2);
            [ux_pod, uy_pod] = lgrad(reshape(pod_u(:,i),size(x)),0);
            [vx_pod, vy_pod] = lgrad(reshape(pod_v(:,i),size(x)),0);
            ux_pods(:,i) = reshape(ux_pod, nL,1);
            uy_pods(:,i) = reshape(uy_pod, nL,1);
            vx_pods(:,i) = reshape(vx_pod, nL,1);
            vy_pods(:,i) = reshape(vy_pod, nL,1);
        end;
        pod_uv = [Me.*pod_u(:,2:end);Me.*pod_v(:,2:end)];
    end;

    u_fom = pod_u*ucoef;
    v_fom = pod_v*ucoef;

    if false; % Normal way of calculating the gradient
        [ux_fom, uy_fom] = lgrad(u_fom, 0);
        [vx_fom, vy_fom] = lgrad(v_fom, 0);
    else
        ux_fom = ux_pods*ucoef;
        uy_fom = uy_pods*ucoef;
        vx_fom = vx_pods*ucoef;
        vy_fom = vy_pods*ucoef;
    end

    conv_u_fom = u_fom.*ux_fom + v_fom.*uy_fom;
    conv_v_fom = u_fom.*vx_fom + v_fom.*vy_fom;
    
    out_coef = pod_uv'*[conv_u_fom; conv_v_fom];
end
