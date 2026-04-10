function [index] = s_opt_qr_fast(Vo, M, index, outfile)
    [N_Bm, N] = size(Vo);
    nVo = Vo.^2; 
    final_index = zeros(M, 1);
    
    % --- Initialization ---
    inum = length(index);
    if inum > 0
        final_index(1:inum) = index;
        [~, R] = qr(Vo(final_index(1:inum), :), 0);
    else
        [~, i_idx] = max(sum(nVo, 2));
        final_index(1) = i_idx;
        inum = 1;
        R = zeros(N, N); % Initialize as square for consistent updates
        R(1, :) = Vo(i_idx, :);
    end

    for i = inum + 1 : M
        % 1. Leverage Scores: Score = ||v / R||^2
        % Using the triangular property of R for speed
        y = Vo / R; 
        scores = sum(y.^2, 2);
        
        % 2. S-Optimization Criteria
        nV = sum(nVo(final_index(1:i-1), :), 1);
        noM = sum(log(nV + nVo), 2);
        A = log(1 + scores + eps) - noM;
        
        % Selection
        A(final_index(1:i-1)) = -inf;
        [~, best_idx] = max(A);
        best_idx = best_idx(1);
        final_index(i) = best_idx;

        % 3. Manual O(N^2) Givens Update
        % We "push" the new row into the R matrix
        new_row = Vo(best_idx, :);
        for j = 1:N
            % Generate rotation to zero out new_row(j) using R(j,j)
            [G, y_rot] = planerot([R(j,j); new_row(j)]);
            
            % Update the diagonal element
            R(j,j) = y_rot(1);
            new_row(j) = y_rot(2); % This will be 0, but helps track the vector
            
            % Apply rotation to the rest of the rows (Vectorized)
            if j < N
                combined = G * [R(j, j+1:end); new_row(j+1:end)];
                R(j, j+1:end) = combined(1, :);
                new_row(j+1:end) = combined(2, :);
            end
        end
        
        if mod(i, floor(M/4)) == 0
            fprintf('Progress: %d/%d (Fast QR Update)\n', i, M);
        end
    end
    
    index = final_index;
    if nargin > 3 && ~isempty(outfile), save([outfile, '.mat'], 'index'); end
end
