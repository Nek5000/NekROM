function [u, v] = get_snaps(snaps, reorder)
    % Default reorder to 0 if not provided
    if nargin < 2
        reorder = false;
    end

    % Pre-calculate dimensions
    num_fields = numel(snaps.flds);
    nL = numel(snaps.flds{1}.u);
    
    % Initialize output matrices
    u = zeros(nL, num_fields);
    v = zeros(nL, num_fields);

    if reorder
        % Pre-calculate sort order once
        Ie = get_sort_order(snaps.flds{1}.x, snaps.flds{1}.y);
        
        for i = 1:num_fields
            % Use linear indexing with Ie for faster reordering and flattening
            temp_u = snaps.flds{i}.u(:,:,Ie);
            temp_v = snaps.flds{i}.v(:,:,Ie);
            u(:, i) = temp_u(:);
            v(:, i) = temp_v(:);
        end
    else
        for i = 1:num_fields
            % Direct flattening using the colon operator
            u(:, i) = snaps.flds{i}.u(:);
            v(:, i) = snaps.flds{i}.v(:);
        end
    end
end
