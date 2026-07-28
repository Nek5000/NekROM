function[u, v] = get_snaps(snaps, reorder)

  if nargin < 2;
    reorder = 0;
  end

  if reorder;
    Ie = get_sort_order(snaps.flds{1}.x, snaps.flds{1}.y);
  end

  u = snaps.flds{1}.u;
  nL = prod(size(u),"all");%nr*ns*nE;
  disp('Field size');
  size(snaps.flds)
  [nbasis, nbasis1] = size(snaps.flds);
  u = zeros(nL,nbasis);
  v = zeros(nL,nbasis);

  for i=1:nbasis;
    if reorder
        u(:,i) = reshape(snaps.flds{i}.u(:,:,Ie),nL,1);
        v(:,i) = reshape(snaps.flds{i}.v(:,:,Ie),nL,1);
    else
        u(:,i) = reshape(snaps.flds{i}.u, nL,1);
        v(:,i) = reshape(snaps.flds{i}.v,nL,1);
    end;
  end;
end


