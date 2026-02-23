function [Au, Bu] = gen_Au(pod_u, pod_v, x, y)
    % ROM diffusion operator - Hybrid MATLAB/Octave Implementation
    
    [nx, ny] = size(x);
    nL = nx * ny;
    nPOD = size(pod_u, 2);
    
    % 1. Setup Operators & Geometry
    [zi, w] = zwgll(nx-1);
    Dh = deriv_mat(zi);
    Dht = Dh';
    W = w * w'; 
    [~,~,~,~,rx,ry,sx,sy,jac,jaci,~] = deriv_geo(x, y, Dh);
    
    % 2. Calculate Geometric Factors
    Me_inv = jaci .* W;
    Grr = reshape(Me_inv .* (rx.^2 + ry.^2), nL, 1);
    Grs = reshape(Me_inv .* (rx.*sx + ry.*sy), nL, 1);
    Gss = reshape(Me_inv .* (sx.*sx + sy.*sy), nL, 1);

    % 3. Branching Logic for Differentiation
    if exist('pagemtimes', 'builtin') || exist('pagemtimes', 'file')
        % --- MATLAB Optimized Path ---
        u_3d = reshape(pod_u, nx, ny, nPOD);
        v_3d = reshape(pod_v, nx, ny, nPOD);
        
        ur = reshape(pagemtimes(Dh, u_3d), nL, nPOD);
        us = reshape(pagemtimes(u_3d, Dht), nL, nPOD);
        vr = reshape(pagemtimes(Dh, v_3d), nL, nPOD);
        vs = reshape(pagemtimes(v_3d, Dht), nL, nPOD);
    else
        % --- Octave / Legacy Path ---
        % d/dr (Left derivative) is vectorized via matrix stacking
        ur = reshape(Dh * reshape(pod_u, nx, ny * nPOD), nL, nPOD);
        vr = reshape(Dh * reshape(pod_v, nx, ny * nPOD), nL, nPOD);
        
        % d/ds (Right derivative) requires a loop for the column-wise mult
        us = zeros(nL, nPOD); 
        vs = zeros(nL, nPOD);
        for i = 1:nPOD
            u_m = reshape(pod_u(:,i), nx, ny);
            v_m = reshape(pod_v(:,i), nx, ny);
            us(:,i) = reshape(u_m * Dht, nL, 1);
            vs(:,i) = reshape(v_m * Dht, nL, 1);
        end
    end

    % 4. Build Au (Symmetric Diffusion Operator)
    % Grouping u and v terms for numerical stability and clarity
    Au = (ur' * (Grr .* ur + Grs .* us) + us' * (Grs .* ur + Gss .* us)) + ...
         (vr' * (Grr .* vr + Grs .* vs) + vs' * (Grs .* vr + Gss .* vs));

    % 5. Build Bu (Mass Matrix)
    Me_vec = reshape(jac .* W, nL, 1);
    Bu = pod_u' * (Me_vec .* pod_u) + pod_v' * (Me_vec .* pod_v);
end
