function [indices] = gappy_pod_fast(U_nl, n, method)
    % Gappy POD with Optimized QR Path (lambda = 0)
    % U_nl:   Basis matrix (N x p)
    % n:      Total number of indices to select
    % method: (Optional) 'synced' [default] or 'spread'

    % See https://epubs.siam.org/doi/pdf/10.1137/22M1484018
    % This is essentially DEIM with oversampling.
    % Code optimized using Gemini
    % 'Spread' is a roundrobin distribution of the points
    % over the basis vectors. 'Sync
    
    if nargin < 3 || isempty(method)
        method = 'synced'; 
    end
    
    [N, p] = size(U_nl);
    indices = zeros(1, n);
    
    % Initialize with the first index based on the first mode
    [~, indices(1)] = max(abs(U_nl(:, 1)));
    
    if strcmp(method, 'synced')
        % Paper logic: nested loops to sync samples to specific modes
        % For the oversampled case, this spends more of the budget when
        % fewer basis vectors are present.
        n_iter = ceil((n - 1) / (p - 1));
        curr = 1;
        for i = 2:p
            for k = 1:n_iter
                if curr >= n, break; end
                
                % Standard DEIM/Gappy: Approx mode 'i' using modes '1:i-1'
                U_sub = U_nl(:, 1:i-1);
                idx_prev = indices(1:curr);
                
                % High-speed QR solve for coefficients
                c = U_sub(idx_prev, :) \ U_nl(idx_prev, i);
                
                % Residual calculation and selection
                r = U_nl(:, i) - U_sub * c;
                [~, next_idx] = max(abs(r));
                
                curr = curr + 1;
                indices(curr) = next_idx;
            end
            if curr >= n, break; end
        end
    else
        % Spread logic: single loop, distributing modes gradually
        % For the oversampled case, this spends more of the budget
        % after all basis vectors are present.
        for j = 2:n
            col_idx = min(p, ceil(j * p / n));
            U_sub = U_nl(:, 1:col_idx-1);
            idx_prev = indices(1:j-1);
            
            c = U_sub(idx_prev, :) \ U_nl(idx_prev, col_idx);
            
            r = U_nl(:, col_idx) - U_sub * c;
            [~, next_idx] = max(abs(r));
            indices(j) = next_idx;
        end
    end
end
