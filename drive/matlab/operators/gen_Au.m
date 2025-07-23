% ROM diffusion operator.
% This should have the same action as the Au operator from the Fortran code
% Should probably add convection tensor calculation as well.
function[Au, Bu, u0] = gen_Au(pod_u, pod_v, snaps)
    x=snaps.flds{1}.x;
    y=snaps.flds{1}.y;
    nx1 = size(x,1);
    [zi, w] = zwgll(nx1-1);
    Dh = deriv_mat(zi);
    Dht = Dh';
    [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,Dh);
    nL = prod(size(x));
    % The rx etc. arrays aren't multiplied by jaci. Assigning Me = massmatrix / Jac does the same thing.
    Me = jaci.*(w*w');

    % Calculate geometric factors
    Grr = reshape(Me.*(rx.*rx + ry.*ry),nL,1);
    Grs = reshape(Me.*(rx.*sx + ry.*sy),nL,1);
    Gss = reshape(Me.*(sx.*sx + sy.*sy),nL,1);

    % Kento's ROM approach. Calculate the derivatives of the POD modes
    ur_pods = [];
    us_pods = [];
    vr_pods = [];
    vs_pods = [];
    for i = 1:size(pod_u,2);
        pod_u_vec = reshape(pod_u(:,i),size(x));
        pod_v_vec = reshape(pod_v(:,i),size(x));
        ur_pods = [ur_pods, reshape(pagemtimes(Dh,pod_u_vec), nL,1)];
        us_pods = [us_pods, reshape(pagemtimes(pod_u_vec,Dht), nL,1)];
        vr_pods = [vr_pods, reshape(pagemtimes(Dh,pod_v_vec), nL,1)];
        vs_pods = [vs_pods, reshape(pagemtimes(pod_v_vec,Dht), nL,1)];
        %ur_pods = [ur_pods, reshape(tensorprod(Dh, pod_u_vec, [2],[1]), nL,1)];
        %us_pods = [us_pods, reshape(tensorprod(pod_u_vec, Dht,[2],[1]), nL,1)];
        %vr_pods = [vr_pods, reshape(tensorprod(Dh, pod_v_vec, [2],[1]), nL,1)];
        %vs_pods = [vs_pods, reshape(tensorprod(pod_v_vec, Dht,[2],[1]), nL,1)];

    end;
    Au = [ur_pods;us_pods]'* [diag(Grr)*ur_pods + diag(Grs)*us_pods; diag(Grs)*ur_pods + diag(Gss)*us_pods] + [vr_pods;vs_pods]'* [diag(Grr)*vr_pods + diag(Grs)*vs_pods; diag(Grs)*vr_pods + diag(Gss)*vs_pods]; 

    Me = reshape(jac.*(w*w'),nL,1);
    bas = [pod_u;pod_v];
    Bu = bas'*sparse(diag([Me;Me]))*bas; 
end
