function [inds, weights, info] = compressed_quadrature(phi, mass_weights, target_count, candidate_inds)
    % COMPRESSED_QUADRATURE Choose a sparse positive quadrature rule.
    %
    % This builds a positive weight vector `weights` on a subset of row
    % indices `inds` so that the weighted Gram matrix of `phi` matches the
    % full Gram matrix induced by `mass_weights` in a least-squares sense:
    %
    %   phi' * diag(mass_weights) * phi  ~=  phi(inds,:)' * diag(weights) * phi(inds,:)
    %
    % The intended use is a dealiased (overintegrated) DEIM-family runtime:
    % - Evaluate nonlinear terms at `inds`
    % - Multiply by sqrt(weights)
    % - Use weighted gappy-POD / CLSDEIM / MCLSDEIM machinery unchanged.
    %
    % Inputs:
    %   phi          [n x r] basis values on the quadrature grid (unweighted)
    %   mass_weights [n x 1] positive quadrature/mass weights (e.g. Me)
    %   target_count number of points to keep
    %   candidate_inds (optional) candidate point set to select from
    %
    % Outputs:
    %   inds    selected point indices into rows of phi
    %   weights positive weights aligned with inds
    %   info    diagnostics struct

    if nargin < 3 || isempty(target_count)
        target_count = min(size(phi, 1), 3 * size(phi, 2));
    end
    if nargin < 4
        candidate_inds = [];
    end

    n = size(phi, 1);
    r = size(phi, 2);
    mass_weights = mass_weights(:);
    assert(numel(mass_weights) == n, 'compressed_quadrature:WeightMismatch', ...
        'mass_weights must have length %d.', n);
    assert(all(mass_weights > 0), 'compressed_quadrature:NonPositiveWeights', ...
        'mass_weights must be strictly positive.');

    target_count = max(target_count, r);
    target_count = min(target_count, n);

    if isempty(candidate_inds)
        % Use QDEIM on the weighted basis to bias selection toward points with
        % large contribution in the mass inner product.
        weighted_phi = bsxfun(@times, sqrt(mass_weights), phi);
        inds = qdeim(weighted_phi, target_count, []);
    else
        candidate_inds = candidate_inds(:);
        candidate_inds = candidate_inds(candidate_inds >= 1 & candidate_inds <= n);
        candidate_inds = unique(candidate_inds, 'stable');
        assert(~isempty(candidate_inds), 'compressed_quadrature:EmptyCandidates', ...
            'candidate_inds is empty after filtering.');
        weighted_phi = bsxfun(@times, sqrt(mass_weights(candidate_inds)), phi(candidate_inds, :));
        local = qdeim(weighted_phi, min(target_count, numel(candidate_inds)), []);
        inds = candidate_inds(local);
    end
    inds = inds(:);

    % Target Gram matrix vector (upper triangle).
    gram_full = phi' * bsxfun(@times, mass_weights, phi);
    tri_mask = triu(true(r));
    g = gram_full(tri_mask);

    % Build linear system A * w ~= g, where each point contributes
    % w_i * vec(triu(phi_i * phi_i')) to the Gram matrix.
    k = numel(inds);
    A = zeros(numel(g), k);
    for j = 1:k
        v = phi(inds(j), :).';
        outer = v * v.';
        A(:, j) = outer(tri_mask);
    end

    % Nonnegative least squares for positive quadrature weights.
    [weights, nnls_info] = nnls_lawson_hanson(A, g);

    % Prune tiny weights (keeps runtime cheaper and avoids ill-conditioning).
    wmax = max(weights);
    if wmax <= 0
        error('compressed_quadrature:ZeroWeights', 'Failed to compute a positive quadrature rule.');
    end
    prune_tol = 1e-14 * wmax;
    keep = weights > prune_tol;
    if sum(keep) < r
        % Ensure the reduced rule remains usable for r-dimensional least squares.
        [~, order] = sort(weights, 'descend');
        keep = false(size(weights));
        keep(order(1:r)) = true;
    end
    inds = inds(keep);
    weights = weights(keep);

    info = struct();
    info.requested_count = target_count;
    info.selected_count = numel(inds);
    info.nnls = nnls_info;
    info.residual_norm = norm(A(:, keep) * weights - g);
    info.target_norm = norm(g);
end

function [x, info] = nnls_lawson_hanson(A, b)
    % NNLS_LAWSON_HANSON Simple active-set NNLS solver (no toolboxes).
    %
    % Solves min ||A x - b||_2 s.t. x >= 0.

    [m, n] = size(A);
    b = b(:);
    assert(numel(b) == m, 'nnls_lawson_hanson:DimMismatch', 'b must have length %d.', m);

    % Tolerances loosely following Lawson-Hanson recommendations.
    tol = 10 * eps(max(norm(A, 1), 1));
    max_iter = max(30 * n, 200);

    x = zeros(n, 1);
    passive = false(n, 1);
    iter = 0;

    % Dual variables.
    w = A' * (b - A * x);

    while true
        iter = iter + 1;
        if iter > max_iter
            break;
        end

        % Find the most positive w among active (non-passive) indices.
        w_active = w;
        w_active(passive) = -inf;
        [wmax, t] = max(w_active);
        if ~(wmax > tol)
            break;
        end

        passive(t) = true;

        while true
            p_idx = find(passive);
            if isempty(p_idx)
                break;
            end

            s = zeros(n, 1);
            % Least squares over passive set.
            s(p_idx) = A(:, p_idx) \ b;

            if all(s(p_idx) > tol)
                x = s;
                break;
            end

            % Move toward s until some component hits zero.
            neg = (s <= tol) & passive;
            alpha = inf;
            neg_idx = find(neg);
            for ii = 1:numel(neg_idx)
                i = neg_idx(ii);
                denom = x(i) - s(i);
                if denom > 0
                    alpha = min(alpha, x(i) / denom);
                end
            end
            if ~isfinite(alpha)
                alpha = 0;
            end
            x = x + alpha * (s - x);
            x(x < tol) = 0;
            passive = passive & (x > 0);
        end

        w = A' * (b - A * x);
    end

    info = struct();
    info.iterations = iter;
    info.tolerance = tol;
    info.max_iter = max_iter;
    info.primal_residual = norm(A * x - b);
end
