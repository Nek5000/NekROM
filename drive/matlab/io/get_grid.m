function [x_fom, y_fom] = get_grid(snaps, reorder)
  if nargin < 2
    reorder = 0;
  end;

  x_fom = snaps.flds{1}.x;
  y_fom = snaps.flds{1}.y;

  if reorder;
    Ie = get_sort_order(x_fom, y_fom);
    x_fom = x_fom(:,:,Ie);
    y_fom = y_fom(:,:,Ie);
  end
end

