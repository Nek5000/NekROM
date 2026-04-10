function [Au, Bu] = gen_Au(pod_u, pod_v, x, y)
% GEN_AU  SEM-consistent ROM diffusion and mass operators
% Compatible with MATLAB and GNU Octave

%% Grid sizes
nx1 = size(x,1);
ne  = size(x,3);
nPOD = size(pod_u,2);
nL = nx1*nx1*ne;

use_pagemtimes = (exist('pagemtimes','builtin') == 5);

%% SEM operators
[zi,w] = zwgll(nx1-1);
Dh  = deriv_mat(zi);
Dht = Dh';

W = w*w';

%% Geometry
[~,~,~,~,rx,ry,sx,sy,jac,jaci,~] = deriv_geo(x,y,Dh);

%% Geometric coefficients
Me_inv = jaci .* W;

Grr = reshape(Me_inv .* (rx.^2 + ry.^2), nL, 1);
Grs = reshape(Me_inv .* (rx.*sx + ry.*sy), nL, 1);
Gss = reshape(Me_inv .* (sx.^2 + sy.^2), nL, 1);

%% Reshape POD modes to SEM layout
u = reshape(pod_u, nx1, nx1, ne, nPOD);
v = reshape(pod_v, nx1, nx1, ne, nPOD);

%% Allocate derivatives
ur = zeros(nx1,nx1,ne,nPOD);
us = ur;
vr = ur;
vs = ur;

%% Compute derivatives

if use_pagemtimes

    % MATLAB fast path
    for e = 1:ne
        ur(:,:,e,:) = pagemtimes(Dh,  u(:,:,e,:));
        us(:,:,e,:) = pagemtimes(u(:,:,e,:),  Dht);

        vr(:,:,e,:) = pagemtimes(Dh,  v(:,:,e,:));
        vs(:,:,e,:) = pagemtimes(v(:,:,e,:),  Dht);
    end

else

    % Octave fallback
    for e = 1:ne
        for k = 1:nPOD

            u_mode = u(:,:,e,k);
            v_mode = v(:,:,e,k);

            ur(:,:,e,k) = Dh * u_mode;
            us(:,:,e,k) = u_mode * Dht;

            vr(:,:,e,k) = Dh * v_mode;
            vs(:,:,e,k) = v_mode * Dht;

        end
    end

end

%% Flatten derivatives
ur = reshape(ur,nL,nPOD);
us = reshape(us,nL,nPOD);
vr = reshape(vr,nL,nPOD);
vs = reshape(vs,nL,nPOD);

%% Diffusion operator
Au = ...
    (ur' * (Grr .* ur + Grs .* us) + ...
     us' * (Grs .* ur + Gss .* us)) ...
  + ...
    (vr' * (Grr .* vr + Grs .* vs) + ...
     vs' * (Grs .* vr + Gss .* vs));

%% Mass matrix
Me_vec = reshape(jac .* W, nL, 1);

Bu = pod_u' * (Me_vec .* pod_u) ...
   + pod_v' * (Me_vec .* pod_v);

end
