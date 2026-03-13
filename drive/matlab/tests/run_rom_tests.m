function run_rom_tests(snaps_path,casename,nb,x_fom,y_fom,pod_u,pod_v,au_full,bu_full)

disp("Running ROM validation tests")

reorder=1;
subtract_mean=1;
conserve_momentum=0;

snaps_obj=NekSnaps([snaps_path,casename]);

[pod_ml,u0_ml,uk_ml] = ...
    get_pod_basis(snaps_obj,nb,reorder,subtract_mean,conserve_momentum);

pod_u_ml = pod_ml(1:size(pod_ml,1)/2,1:nb+1);
pod_v_ml = pod_ml(size(pod_ml,1)/2+1:end,1:nb+1);

[au_ml,bu_ml] = gen_Au(pod_u_ml,pod_v_ml,x_fom,y_fom);

%% ----------------------------------------------------
%% POD BASIS COMPARISON
%% ----------------------------------------------------

pod=[pod_u;pod_v];

err = norm(abs(pod)-abs(pod_ml))/norm(abs(pod));

disp(["POD basis relative error: ",num2str(err)])

assert(err < 1e-5)

%% ----------------------------------------------------
%% OPERATOR COMPARISON
%% ----------------------------------------------------

errA = norm(abs(au_full)-abs(au_ml))/norm(abs(au_full));
errB = norm(abs(bu_full)-abs(bu_ml))/norm(abs(bu_full));

disp(["Au error: ",num2str(errA)])
disp(["Bu error: ",num2str(errB)])

assert(errA < 1e-5)
assert(errB < 1e-5)

disp("All ROM tests passed")

end
