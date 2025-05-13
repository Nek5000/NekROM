% This needs to do the same thing as
% reshape(cu*utmp(:,1),nb,nb+1)*u(:,1);
function [out_coef] = conv_fom(ucoef, pod_u, pod_v, snaps)
    x=snaps.flds{1}.x;
    y=snaps.flds{1}.y;
    nx1 = size(x,1);
    [zi, w] = zwgll(nx1-1);
    d = deriv_mat(zi);
    [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
    lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
    nL = prod(size(x));
    Me = reshape(jac.*(w*w'),nL,1);

    if false; % Normal way of calculating the gradient
        [ux_fom, uy_fom] = lgrad(u_fom, 0);
        [vx_fom, vy_fom] = lgrad(v_fom, 0);
    else
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
        ux_fom = reshape(ux_pods*ucoef, size(x));
        uy_fom = reshape(uy_pods*ucoef, size(x));
        vx_fom = reshape(vx_pods*ucoef, size(x));
        vy_fom = reshape(vy_pods*ucoef, size(x));
    end

    u_fom = reshape(Me.*pod_u*ucoef, size(x));
    v_fom = reshape(Me.*pod_v*ucoef, size(x));

    conv_u_fom = reshape(u_fom.*ux_fom + v_fom.*uy_fom, nL,1);
    conv_v_fom = reshape(u_fom.*vx_fom + v_fom.*vy_fom, nL,1);
    
    out_coef = [pod_u(:,2:end); pod_v(:,2:end)]'*[conv_u_fom; conv_v_fom];
end

