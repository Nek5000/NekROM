function [bas, u0, uk] = get_pod_basis(snaps_obj, nb, reorder, subtract_mean, conserve_momentum)
        % Assemble the snapshots
        [u_snaps, v_snaps] = get_snaps(snaps_obj, reorder);
        snaps = [u_snaps; v_snaps];

        [x,y] = get_grid(snaps_obj, reorder);

        [bas, u0, uk] = get_pod_basis_from_arrays(snaps, x, y, nb, subtract_mean, conserve_momentum);
end
