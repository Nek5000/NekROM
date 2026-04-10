function [index] = s_opt_final(Vo, M, index, outfile)
% S_OPT_FINAL High-performance S-optimal index selection.
%   Uses QR-decomposition updates and incremental log-sum tracking.

    [N_Bm, N] = size(Vo);
    nVo = Vo.^2; % Pre-compute squared elements
    
    % --- Initialization ---
    final_index = zeros(M, 1);
    if isempty(index)
        % Initial step: Maximize row norm as a strong starting point
        [~, best_idx] = max(sum(nVo, 2));
        final_index(1) = best_idx;
        inum = 1;
    else
        inum = length(index);
        final_index(1:inum) = index;
    end

    % Pre-compute running sum of squares (S) for the denominator term
    % This avoids O(N*N_Bm*i) re-calculations
    S = sum(nVo(final_index(1:inum), :), 1);
    
    % Initialize R (Upper triangular factor)
    % Building initial R from provided indices
    [~, R] = qr(Vo(final_index(1:inum), :), 0);

    % --- Main Optimization Loop ---
    for i = inum + 1 : M
        % 1. Compute Leverage Scores via QR substitution
        % This is the stable equivalent of diag(Vo * inv(V'*V) * Vo')
        % complexity: O(N_Bm * N^2)
        tmp = R' \ Vo';
        leverage = sum(Vo .* (R \ tmp)', 2);
        
        % 2. Vectorized Incremental Log-Determinant Calculation
        % Complexity: O(N_Bm * N)
        noM = sum(log(S + nVo + eps), 2);
        
        % 3. Calculate Selection Criterion
        % Score = log(1 + leverage) - sum(log(Selected_Squares + Candidate_Squares))
        A = log(1 + leverage + eps) - noM;
        
        % Avoid re-selecting existing indices
        A(final_index(1:i-1)) = -inf;
        
        % 4. Select Best Index
        [~, best_idx] = max(A);
        best_idx = best_idx(1);
        final_index(i) = best_idx;
        
        % 5. Rank-1 Update of S and R
        % Update square sums
        S = S + nVo(best_idx, :);
        
        % Update QR factor (Appending a row and re-triangularizing)
        % For N > 500, consider using 'qrupdate' if available in your toolbox
        [~, R] = qr([R; Vo(best_idx, :)], 0);

        % Progress Reporting
        if mod(i, max(1, floor(M/4))) == 0
            fprintf('S-Opt Progress: %d/%d points...\n', i, M);
        end
    end
    
    index = final_index;
    if nargin > 3 && ~isempty(outfile)
        save([outfile, '.mat'], 'index');
    end
end
