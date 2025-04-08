function[pod_u, pod_v, x_fom, y_fom] = get_pod(snaps)
  u = snaps.flds{1}.u;
  %y_fom = snaps.flds{1}.y;

  [nr, ns, nE] = size(u);
  nL = prod(size(u),"all");%nr*ns*nE;
  [nbasis, nbasis1] = size(snaps.flds)
  pod_u = [];
  pod_v = [];
  for i=1:nbasis;
    %u_snap = snaps.flds{i}.u;
    pod_u = [pod_u, reshape(snaps.flds{i}.u, nL,1)];
    pod_v = [pod_v, reshape(snaps.flds{i}.v,nL,1)];
  end;


  % The bases are not 
  %[pod_u;pod_v]'*[pod_u;pod_v]
  %exit;

  %avg_snaps = NekSnaps(avg_cname);
  %u_avg = reshape(avg_snaps.flds{1}.u,nL,1);
  %v_avg = reshape(avg_snaps.flds{1}.v,nL,1);
  %% Technically not the POD anymore
  %exit;
  %pod_u(:,1) = u_avg(:,1);
  %pod_v(:,1) = v_avg(:,1);
end


