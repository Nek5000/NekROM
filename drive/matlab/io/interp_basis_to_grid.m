function [u_tgt, v_tgt] = interp_basis_to_grid(u_src, v_src, x_src, y_src, x_tgt, y_tgt)
    % INTERP_BASIS_TO_GRID Interpolate stacked basis vectors between tensor-product grids.
    %
    % The source and target grids must have the same number of elements and use
    % square 2D spectral element layouts. If the grids already match, the
    % inputs are returned unchanged.

    if isequal(size(x_src), size(x_tgt)) && isequal(size(y_src), size(y_tgt))
        u_tgt = u_src;
        v_tgt = v_src;
        return;
    end

    assert(ndims(x_src) == 3 && ndims(x_tgt) == 3, ...
        'interp_basis_to_grid:GridRank', 'Only 2D tensor-product grids are supported.');

    nx_src = size(x_src, 1);
    ny_src = size(x_src, 2);
    nx_tgt = size(x_tgt, 1);
    ny_tgt = size(x_tgt, 2);
    nelt_src = size(x_src, 3);
    nelt_tgt = size(x_tgt, 3);

    assert(nx_src == ny_src && nx_tgt == ny_tgt, ...
        'interp_basis_to_grid:SquareGrid', 'Only square 2D elements are supported.');
    assert(nelt_src == nelt_tgt, ...
        'interp_basis_to_grid:ElementCount', 'Source and target grids must have the same number of elements.');
    assert(size(u_src, 1) == nx_src * ny_src * nelt_src && size(v_src, 1) == nx_src * ny_src * nelt_src, ...
        'interp_basis_to_grid:SourceSize', 'Basis rows must match the source grid.');

    [zi_src, ~] = zwgll(nx_src - 1);
    [zi_tgt, ~] = zwgll(nx_tgt - 1);
    interp_1d = interp_mat(zi_tgt, zi_src);

    n_modes = size(u_src, 2);
    u_tgt = zeros(nx_tgt * ny_tgt * nelt_tgt, n_modes);
    v_tgt = zeros(nx_tgt * ny_tgt * nelt_tgt, n_modes);

    src_block = nx_src * ny_src;
    tgt_block = nx_tgt * ny_tgt;

    for ie = 1:nelt_src
        src_idx = (ie - 1) * src_block + (1:src_block);
        tgt_idx = (ie - 1) * tgt_block + (1:tgt_block);

        for im = 1:n_modes
            u_elem = reshape(u_src(src_idx, im), nx_src, ny_src);
            v_elem = reshape(v_src(src_idx, im), nx_src, ny_src);

            u_tgt(tgt_idx, im) = reshape(interp_1d * u_elem * interp_1d', [], 1);
            v_tgt(tgt_idx, im) = reshape(interp_1d * v_elem * interp_1d', [], 1);
        end
    end
end
