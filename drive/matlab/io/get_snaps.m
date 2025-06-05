function[u, v] = get_snaps(snaps)
  u = snaps.flds{1}.u;
  %y_fom = snaps.flds{1}.y;

  %[nr, ns, nE] = size(u);
  nL = prod(size(u),"all");%nr*ns*nE;
  [nbasis, nbasis1] = size(snaps.flds)
  u = [];
  v = [];
  for i=1:nbasis;
    %u_snap = snaps.flds{i}.u;
    u = [u, reshape(snaps.flds{i}.u, nL,1)];
    v = [v, reshape(snaps.flds{i}.v,nL,1)];
  end;


  % The bases are not 
  %[u;v]'*[u;v]
  %exit;

  %avg_snaps = NekSnaps(avg_cname);
  %u_avg = reshape(avg_snaps.flds{1}.u,nL,1);
  %v_avg = reshape(avg_snaps.flds{1}.v,nL,1);
  %% Technically not the POD anymore
  %exit;
  %u(:,1) = u_avg(:,1);
  %v(:,1) = v_avg(:,1);
end


