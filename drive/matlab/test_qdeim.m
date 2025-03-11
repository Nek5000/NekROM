cname='../csncyl';
snaps=NekSnaps(cname);
snaps.pfun = 'u_abs';
snaps.plot_type='contour';
snaps

%exit;
%size(snaps.isnaps)
%snaps.isnaps
%class(snaps.flds{100})
%snaps.flds{100}.u

u_snaps = [];
v_snaps = [];
class(snaps.flds{1}.u);
[nr, ns, nE] = size(snaps.flds{1}.x);
nL = nr*ns*nE;
for i=1:500;
    %u_snap = snaps.flds{i}.u;
    u_snaps = [u_snaps, reshape(snaps.flds{i}.u,nL,1)];
    v_snaps = [v_snaps, reshape(snaps.flds{i}.v,nL,1)];
end;

%size(u_snaps)
[P_u, u_indices] = calc_qdeim_proj_mat(u_snaps);
[P_v, v_indices] = calc_qdeim_proj_mat(v_snaps);
[P_uv, uv_indices] = calc_qdeim_proj_mat([u_snaps;v_snaps]);


xvals = reshape(snaps.flds{1}.x,nL,1);
yvals = reshape(snaps.flds{1}.y,nL,1);

xxvals = [xvals;xvals];
yyvals = [yvals;yvals];

snaps.first();
%exit;

pause(1);
for i=1:100;
snaps.next();
hold on;
plot(xvals(u_indices), yvals(u_indices), '.', 'MarkerSize', 22); hold on;
plot(xvals(v_indices), yvals(v_indices), '.', 'MarkerSize', 22); hold on;
plot(xxvals(uv_indices), yyvals(uv_indices), '.', 'MarkerSize', 22); hold off;

pause(.01);
end;
pause(30);

