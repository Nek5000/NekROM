% Code to generate the basis functions in Matlab for testing purposes
%Maybe define some unit test and move these tests to them?

if 0
    subtract_mean = 1;
    conserve_momentum = 0;
    snaps_obj = NekSnaps(strcat(snaps_path,casename)); % Should the snaps object have reordering capability?
    [pod_ml, u0_full_ml, uk_full_ml] = get_pod_basis(snaps_obj,nb,reorder,subtract_mean,conserve_momentum);

    pod_u_ml = pod_ml(1:size(pod_ml,1)/2,1:nb+1);
    pod_v_ml = pod_ml(size(pod_ml,1)/2 + 1:end,1:nb+1);

    [au_full_ml, bu_full_ml] = gen_Au(pod_u_ml, pod_v_ml,x_fom, y_fom);
end;


%% Test that the basis vectors are the same, modulo sign differences
if 0;
    pod = [pod_u; pod_v];

    % Make sure momentum conservation is off
    assert(norm(abs(pod) - abs(pod_ml))/norm(abs(pod)) < 1e-5);

    %figure(1);
    %patch_plot(x_fom,y_fom, reshape(pod_v(:,1), size(x_fom)), [], 'PlotType', 'surface');
    %title("NekROM U-POD avg");
    %figure(2);
    %patch_plot(x_fom,y_fom, reshape(pod_v_ml(:,1), size(x_fom)), [], 'PlotType', 'surface');
    %title("Matlab U-POD avg");
end;


%% Test that MATLAB and Fortran Au and Bu operators are the same
if 0;
    %bu_full
    %bu_full_ml
    %abs(bu_full - bu_full_ml
    assert(norm(abs(bu_full) - abs(bu_full_ml))/norm(abs(bu_full)) < 1e-5);
    assert(norm(abs(au_full) - abs(au_full_ml))/norm(abs(au_full)) < 1e-5)
end;

% Use all of the matlab defined basis functions
if 0;
    % Note that the non-linear evaluations are dumped on the same
    % grid as the NekROM basis coordinates. Not necessarily the
    % coordinates from the snapshots (these two coordinates are not necessarily the same
    % apparently)

    %[x_fom_ml, y_fom_ml] = get_grid(snaps_obj, 1);
    %change_basis = [pod_u_ml; pod_v_ml]'*[pod_u;pod_v]; 
    pod_orig = [pod_u;pod_v];
    pod_u = pod_u_ml(:,1:nb+1);
    pod_v = pod_v_ml(:,1:nb+1);
    uk_full = uk_full_ml;
    u0_full = u0_full_ml;
    %u0_full = change_basis*u0_full;
    %get_sort_order(x_fom, y_fom)
    %get_sort_order(x_fom_ml, y_fom_ml)

    %exit;
    %x_fom = x_fom_ml;
    %y_fom = y_fom_ml;
    size(pod_u)
    size(pod_u_ml)

    au_full_orig = au_full;
    bu_full_orig = bu_full;
    au_full = au_full_ml;
    bu_full = bu_full_ml;
end;


%u0_full_orig = u0_full
%u0_full = change_basis*u0_full;
%u0_full
%pod_ml*u0_full
%pod_orig*u0_full_orig
%norm([pod_u_new; pod_v_new]*u0_full - pod_orig*u0_full_orig)/norm(pod_orig*u0_full_orig)
%exit;

%u0_full_roundtrip = change_basis'*u0_full

%% Test out reconstructing the snapshots
if 0; 
[u_snaps, v_snaps] = get_snaps(snaps_obj, 1);
for i=1:size(uk_full,2);
    i
    u = uk_full(:,i);
    u_proj = u_snaps(:,i);
    v_proj = v_snaps(:,i);
    %u_proj = pod_u(:,1:nb+1)*u(1:end,1);
    %v_proj = pod_v(:,1:nb+1)*u(1:end,1);
    %u_proj = nl_snaps_u(:,i);
    %v_proj = nl_snaps_v(:,i);

    if plot_vel_mag      
        u_abs = sqrt(u_proj.^2 + v_proj.^2);
        plot_field = reshape(u_abs, size(x_fom));
    elseif plot_vort;
        plot_field = lcurl(reshape(u_proj,size(x_fom)), reshape(v_proj,size(x_fom)), x_fom, y_fom);%vx-uy;
    end;
    %norm(u_abs)
    hold off;
    % Can plot surface or contour, contour not currently available in main branch though
    patch_plot(x_fom,y_fom, reshape(plot_field,size(x_fom)), [], 'PlotType', 'contour');
    pause(0.01);
end;
exit;
end;
