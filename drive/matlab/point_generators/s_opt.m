function [index] = s_opt(Vo, M, index, outfile)
% S_OPT Generates S-optimal indices of a matrix
%   Optimized implementation focusing on vectorization and memory efficiency.
%   Obtained by querying Gemini for improvements to the original algorithm.
%   Original code in the `extra` folder.  
%
%   Parameters:
%       Vo      : Candidate Matrix (N_Bm x N)
%       M       : Number of S-optimal points to compute
%       index   : Initial index set (use [] if none)
%       outfile : Optional filename string for saving results

    [N_Bm, N] = size(Vo);
    nVo = Vo.^2; % Pre-compute squared elements

    % Pre-allocate the full index vector
    final_index = zeros(M, 1);
    inum = length(index);
    if inum > 0
        final_index(1:inum) = index;
    else
        % Initial step: Maximize the magnitude of the first basis vector
        [~, i_idx] = max(abs(Vo(:, 1)));
        final_index(1) = i_idx;
        inum = 1;
    end

    % Main optimization loop
    for i = inum + 1 : M
        curr_idx = final_index(1:i-1);
        
        if i <= N
            % --- Under-determined or Exact Case ---
            % Current basis subset
            V_prev = Vo(curr_idx, 1:i-1);
            V_curr_col = Vo(curr_idx, i);
            
            % Compute intermediate terms
            atA0 = V_curr_col' * V_prev;
            ata = sum(V_curr_col.^2);
            
            % bbb results in (i-1) x (1 + N_Bm)
            bbb = (V_prev' * V_prev) \ [atA0; Vo(:, 1:i-1)]';
            
            c = bbb(:, 2:end); % (i-1) x N_Bm
            g2 = bbb(:, 1);    % (i-1) x 1
            
            % Vectorized calculation of b (N_Bm x 1)
            % Vo(:, 1:i-1) is (N_Bm x i-1), c' is (N_Bm x i-1)
            b = 1 + sum(Vo(:, 1:i-1) .* c', 2);
            
            % Vectorized calculation of g1 ((i-1) x N_Bm)
            % Implicit expansion: (1 x i-1) + (N_Bm x i-1) -> (N_Bm x i-1)
            g1 = (atA0 + (Vo(:, 1:i-1) .* Vo(:, i)))';
            
            % Vectorized calculation of g3 and GG
            g3 = sum(c' .* g1', 2) ./ b;
            GG = g2 + (c .* (Vo(:, i) - g3)'); 
            
            % Compute selection criterion A
            A_val = ata + Vo(:, i).^2 - sum(g1' .* GG', 2);
            A_val(A_val < 0) = 0;
            
            % Implicit expansion for log-determinant terms
            nV = sum(nVo(curr_idx, 1:i), 1);
            noM = sum(log(nV + nVo(:, 1:i)), 2);
            
            A = log(abs(A_val) + eps) + log(b) - noM;
            
        else
            % --- Over-determined Case ---
            V_sub = Vo(curr_idx, :);
            b_coeff = (V_sub' * V_sub) \ Vo';
            
            nV = sum(nVo(curr_idx, :), 1);
            noM = sum(log(nV + nVo), 2);
            
            A = log(1 + sum(Vo .* b_coeff', 2)) - noM;
        end
        
        % Avoid selecting existing points
        A(curr_idx) = -inf;
        
        % Select the next optimal index
        [~, best_idx] = max(A);
        final_index(i) = best_idx(1);

        % Progress reporting
        if ismember(i, floor(M * [0.25, 0.5, 0.75]))
            fprintf('Progress: %d/%d points computed...\n', i, M);
        end
    end
    
    index = final_index;
    if nargin > 3 && ~isempty(outfile)
        save([outfile, '.mat'], 'index');
    end
end
