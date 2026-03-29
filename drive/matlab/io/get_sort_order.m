function [sort_order] = get_sort_order(x, y)
    % 1. Calculate means by nesting (compatible with Octave & MATLAB)
    % We mean across dim 1, then dim 2
    mean_x = squeeze(mean(mean(x, 1), 2));
    mean_y = squeeze(mean(mean(y, 1), 2));

    % 2. Create a matrix of [index, mean_x, mean_y]
    phys_co = [(1:size(x, 3))', mean_x, mean_y];

    % 3. Sort by Y (col 3) then X (col 2)
    phys_co = sortrows(phys_co, [3, 2]);

    % 4. Extract the sorted indices
    sort_order = phys_co(:, 1);
end
