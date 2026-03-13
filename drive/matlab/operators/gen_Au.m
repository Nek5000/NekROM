function [Au, Bu] = gen_Au(pod_u, pod_v, x, y)
% GEN_AU  Build ROM diffusion (Au) and mass (Bu) operators
%
% Drop-in replacement compatible with MATLAB and Octave.
%
% Inputs
%   pod_u, pod_v : POD velocity modes (nGrid × nPOD)
%   x, y         : grid coordinates
%
% Outputs
%   Au : diffusion operator
%   Bu : mass matrix
%
%
    %% Grid sizes
    [nx, ny] = size(x);
    nL   = nx * ny;
    nPOD = size(pod_u, 2);

    %% 1. Spectral element operators
    [zi, w] = zwgll(nx-1);
    Dh  = deriv_mat(zi);
    Dht = Dh';

    W = w * w';

    %% 2. Geometry factors
    [~,~,~,~,rx,ry,sx,sy,jac,jaci,~] = deriv_geo(x, y, Dh);

    %% 3. Geometric coefficients
    Me_inv = jaci .* W;

    Grr = reshape(Me_inv .* (rx.^2 + ry.^2), nL, 1);
    Grs = reshape(Me_inv .* (rx.*sx + ry.*sy), nL, 1);
    Gss = reshape(Me_inv .* (sx.^2 + sy.^2), nL, 1);

    %% 4. Compute derivatives of POD modes

    % ---- d/dr derivatives (fully vectorized) ----
    ur = reshape(Dh * reshape(pod_u, nx, ny * nPOD), nL, nPOD);
    vr = reshape(Dh * reshape(pod_v, nx, ny * nPOD), nL, nPOD);

    % ---- d/ds derivatives (requires columnwise multiply) ----
    us = zeros(nL, nPOD);
    vs = zeros(nL, nPOD);

    for k = 1:nPOD

        u_mode = reshape(pod_u(:,k), nx, ny);
        v_mode = reshape(pod_v(:,k), nx, ny);

        us(:,k) = reshape(u_mode * Dht, nL, 1);
        vs(:,k) = reshape(v_mode * Dht, nL, 1);

    end

    %% 5. Assemble diffusion operator (symmetric form)

    Au = ...
        (ur' * (Grr .* ur + Grs .* us) + ...
         us' * (Grs .* ur + Gss .* us)) ...
      + ...
        (vr' * (Grr .* vr + Grs .* vs) + ...
         vs' * (Grs .* vr + Gss .* vs));

    %% 6. Mass matrix

    Me_vec = reshape(jac .* W, nL, 1);

    Bu = pod_u' * (Me_vec .* pod_u) ...
       + pod_v' * (Me_vec .* pod_v);

end
