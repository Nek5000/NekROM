function [indices] = gappy_pod_new(U_nl, n, method, tol)
    % Gappy POD with Optimized QR Path and Adaptive Early Exit
    % U_nl:   Basis matrix (N x p)
    % n:      Maximum number of indices to select
    % method: (Optional) 'synced' [default] or 'spread'
    % tol:    (Optional) Early exit threshold. If empty, uses 100*eps*norm(U_nl, 'inf')
    %
    % See https://epubs.siam.org/doi/pdf/10.1137/22M1484018
    % This is essentially DEIM with oversampling.
    % Code optimized using Gemini.
    % The difference between "spread" and "synced" is that in the
    % oversampling case,  synced
    % spends more of its budget when only a few basis vectors
    % are present while spread spends more of its budget when
    % all of them are present.
    
    [N, p] = size(U_nl);
    
    % --- Handle Optional Arguments ---
    if nargin < 3 || isempty(method), method = 'synced'; end
    if nargin < 4 || isempty(tol)
        % Adaptive tolerance based on matrix scale and machine precision
        tol = 100 * eps(class(U_nl)) * max(N, p) * max(abs(U_nl(:)));
    end
    
    indices = [];
    [max_val, first_idx] = max(abs(U_nl(:, 1)));
    
    if max_val < tol
        warning('Basis values are below tolerance. Returning empty indices.');
        return; 
    end
    
    indices = [first_idx];
    curr = 1;

    % --- Selection Logic ---
    if strcmp(method, 'synced')
        n_iter = ceil((n - 1) / (p - 1));
        for i = 2:p
            for k = 1:n_iter
                if curr >= n, break; end
                
                U_sub = U_nl(:, 1:i-1);
                % High-speed QR solve
                c = U_sub(indices, :) \ U_nl(indices, i);
                r = U_nl(:, i) - U_sub * c;
                
                [max_gain, next_idx] = max(abs(r));
                
                if max_gain < tol
                    return; % Early exit: Basis is effectively spanned
                end
                
                curr = curr + 1;
                indices(end+1) = next_idx;
            end
            if curr >= n, break; end
        end
    else
        % Spread logic: Distributes points evenly across modes
        for j = 2:n
            col_idx = min(p, ceil(j * p / n));
            U_sub = U_nl(:, 1:col_idx-1);
            
            c = U_sub(indices, :) \ U_nl(indices, col_idx);
            r = U_nl(:, col_idx) - U_sub * c;
            
            [max_gain, next_idx] = max(abs(r));
            
            if max_gain < tol
                return; 
            end
            
            curr = curr + 1;
            indices(end+1) = next_idx;
        end
    end
end
