function [bas, u0, uk] = get_pod_basis(snaps_obj, nb, reorder, subtract_mean, conserve_momentum, method, inner_product, varargin)
        snaps_obj = apply_snapshot_file_list_subset(snaps_obj);
        % Assemble the snapshots
        [u_snaps, v_snaps] = get_snaps(snaps_obj, reorder);

        [x,y] = get_grid(snaps_obj, reorder);

        mode0_ref = [];
        if ~subtract_mean
            mode0_ref = load_pod_mode0_reference(snaps_obj, reorder);
        end

        if isempty(mode0_ref)
            [bas, u0, uk] = get_pod_basis_from_arrays(u_snaps, v_snaps, x, y, nb, ...
                subtract_mean, conserve_momentum, method, inner_product, varargin{:});
        else
            [bas, u0, uk] = get_pod_basis_from_arrays(u_snaps, v_snaps, x, y, nb, ...
                subtract_mean, conserve_momentum, method, inner_product, mode0_ref, varargin{:});
        end
end
