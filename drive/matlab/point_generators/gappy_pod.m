function [indices] = gappy_pod(U_nl, n, method, tol)
% GAPPY_POD_QR Optimized Gappy POD / DEIM with oversampling.
%
% USAGE:
%    [indices] = gappy_pod(U_nl, n, 'synced', 1e-10)
%
% INPUTS:
%    U_nl   - Basis matrix (N x p), usually POD modes.
%    n      - Maximum number of indices to select.
%    method - 'synced' (default) or 'spread'. 
%    tol    - Early exit threshold for residual norm.

    [N, p] = size(U_nl);
    
    % --- Handle Optional Arguments ---
    if nargin < 3 || isempty(method), method = 'synced'; end
    if nargin < 4 || isempty(tol)
        tol = 100 * eps(class(U_nl)) * max(N, p) * max(abs(U_nl(:)));
    end
    
    % --- Initialization ---
    indices = zeros(n, 1);
    [max_val, first_idx] = max(abs(U_nl(:, 1)));
    
    if max_val < tol
        warning('Basis values are below tolerance. Returning empty.');
        indices = []; return; 
    end
    
    indices(1) = first_idx;
    curr_n = 1;

    if strcmpi(method, 'synced')
        %% --- Synced Logic ---
        % Distribute oversampling iterations across the available basis modes
        n_iter = ceil((n - 1) / (p - 1));
        
        for i = 2:p
            % U_curr contains the modes we are approximating AGAINST
            U_curr = U_nl(:, 1:i-1);
            
            for k = 1:n_iter
                if curr_n >= n, break; end
                
                % Solve the small system: U_curr(indices, :) * c = U_nl(indices, i)
                % We use the backslash operator on the small sampled subset.
                % For very large n, you could maintain a QR here, but for 
                % typical DEIM/Gappy POD, n is small enough that this is fast.
                U_samp = U_curr(indices(1:curr_n), :);
                rhs_samp = U_nl(indices(1:curr_n), i);
                
                % Solve for coefficients
                c = U_samp \ rhs_samp;
                
                % Calculate residual across the full domain
                % Dimension Check: (N x i-1) * (i-1 x 1) = (N x 1)
                r = U_nl(:, i) - U_curr * c;
                [max_gain, next_idx] = max(abs(r));
                
                if max_gain < tol
                    indices = indices(1:curr_n); 
                    return; 
                end
                
                curr_n = curr_n + 1;
                indices(curr_n) = next_idx;
            end
            if curr_n >= n, break; end
        end
        
    else
        %% --- Spread Logic ---
        % 'spread' dynamically increases the basis size as we add points
        for j = 2:n
            % Determine which mode to use as the "target"
            col_limit = min(p, ceil(j * p / n));
            U_curr = U_nl(:, 1:col_limit-1);
            target_col = U_nl(:, col_limit);
            
            U_samp = U_curr(indices(1:curr_n), :);
            rhs_samp = target_col(indices(1:curr_n));
            
            % Solve and find the point of maximum error
            c = U_samp \ rhs_samp;
            r = target_col - U_curr * c;
            
            [max_gain, next_idx] = max(abs(r));
            
            if max_gain < tol
                indices = indices(1:curr_n); 
                return; 
            end
            
            curr_n = curr_n + 1;
            indices(curr_n) = next_idx;
        end
    end
    
    % Trim pre-allocated vector
    indices = indices(1:curr_n);
end
