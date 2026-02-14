function [x_fom, y_fom] = get_grid(snaps, reorder)
  % Handle default value for reorder
  if nargin < 2 || isempty(reorder)
    reorder = false;
  end

  % Extract grid data
  x_fom = snaps.flds{1}.x;
  y_fom = snaps.flds{1}.y;

  % Apply reordering if requested
  if reorder
    % Ie is typically a vector of indices for the 3rd dimension
    Ie = get_sort_order(x_fom, y_fom);
    
    x_fom = x_fom(:, :, Ie);
    y_fom = y_fom(:, :, Ie);
  end
end
