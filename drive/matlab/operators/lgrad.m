function [ux, uy] = lgrad(u, x, y, mode, reset)
    % LGRAD Computes the local gradient on a spectral element.
    % Optional 'reset' argument forces recalculation of geometry.

    persistent geo
    
    % Initialize or Reset geometry data
    if isempty(geo) || (nargin > 4 && reset)
        nx1 = size(x, 1);
        [zi, ~] = zwgll(nx1-1);
        D = deriv_mat(zi);
        
        % Pack geometry into a struct for cleaner handling
        [~,~,~,~, rx, ry, sx, sy, ~, jaci, ~] = deriv_geo(x, y, D);
        
        geo.rx = rx; geo.ry = ry;
        geo.sx = sx; geo.sy = sy;
        geo.jaci = jaci;
        geo.D = D;
    end

    % Compute Gradients
    % Note: Direct calling is often faster than an anonymous function handle in a loop
    [ux, uy] = grad(u, geo.rx, geo.ry, geo.sx, geo.sy, geo.jaci, geo.D, mode);
end
