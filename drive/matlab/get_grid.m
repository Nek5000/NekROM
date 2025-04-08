function [x_fom, y_fom] = get_grid_and_pod(snaps)
  x_fom = snaps.flds{1}.x;
  y_fom = snaps.flds{1}.y;
end

