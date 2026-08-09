function [u, v, w, t] = get_snaps(snaps, reorder)
    % Default reorder to 0 if not provided
    if nargin < 2
        reorder = false;
    end

    % Pre-calculate dimensions
    num_fields = numel(snaps.flds);
    nL = numel(snaps.flds{1}.u);

    has_w = isfield(snaps.flds{1}, 'w') && ~isempty(snaps.flds{1}.w);
    has_t = isfield(snaps.flds{1}, 't') && ~isempty(snaps.flds{1}.t);
    
    % Initialize output matrices
    u = zeros(nL, num_fields);
    v = zeros(nL, num_fields);
    w = [];
    t = [];
    if has_w
        w = zeros(nL, num_fields);
    end
    if has_t
        t = zeros(nL, num_fields);
    end

    if reorder
        % Pre-calculate sort order once
        Ie = get_sort_order(snaps.flds{1}.x, snaps.flds{1}.y);
        
        for i = 1:num_fields
            % Use linear indexing with Ie for faster reordering and flattening
            temp_u = snaps.flds{i}.u(:,:,Ie);
            temp_v = snaps.flds{i}.v(:,:,Ie);
            u(:, i) = temp_u(:);
            v(:, i) = temp_v(:);
            if has_w
                temp_w = snaps.flds{i}.w(:,:,Ie);
                w(:, i) = temp_w(:);
            end
            if has_t
                temp_t = snaps.flds{i}.t(:,:,Ie);
                t(:, i) = temp_t(:);
            end
        end
    else
        for i = 1:num_fields
            % Direct flattening using the colon operator
            u(:, i) = snaps.flds{i}.u(:);
            v(:, i) = snaps.flds{i}.v(:);
            if has_w
                w(:, i) = snaps.flds{i}.w(:);
            end
            if has_t
                t(:, i) = snaps.flds{i}.t(:);
            end
        end
    end
end
