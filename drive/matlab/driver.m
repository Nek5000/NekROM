%#######################################################
%
%# Matlab Driver for Galerkin-based reduced order model 
%# v0.0.1 - Refactored
%
%# Ping-Hsuan Tsai / Refactored
%# 2024-09-05
%
%#######################################################

% Clear variables and command window (avoids performance hit of 'clear all')
clear variables; close all; clc;

% Add any important scripts to path
addpath('./point_generators');
addpath('./io');
addpath('./operators');

%% Specify the case path and case name
cases = ["ldc", "cyl", "shear", "t2d"];
thiscase = cases(3); % Select "shear"

switch thiscase
    case 'ldc'
        path = '../../examples/ldc/';
        nsteps = 10 * 1e5; 
        dt     = 1.0e-03;
        iostep = 1000;
        nu     = 1./15000;
        nb     = 30;
    case 'cyl'
        path = '../../examples/cyl/';
        nsteps = 10 * 1.25e05; 
        dt     = 4.0e-03;
        iostep = 500;
        nu     = 0.01;
        nb     = 20;
    case 'shear'
        path = '../../examples/shear/';
        nsteps = 10 * 4000;
        dt     = 1e-3;
        iostep = 100;
        nu     = 1/1000;
        nb     = 20; 
    case 't2d'
        path = '../../examples/t2d/';
        nsteps = 800000;
        dt     = 0.002;
        iostep = 100;
        nu     = 0.0001;
        nb     = 3;  
    otherwise
        error("Unhandled case name: %s", thiscase);
end

casename = char(thiscase);
snaps_path = strcat(path, 'snaps/');

%% IO parameters
ifvort  = true;  
ifwrite = true;  
ifvis   = false; 

%% ROM stabilization strategies
ifcopt  = false;
ifleray = false;
ifefr   = false;
iftr    = false;

if ifcopt
    reg_str = 'CROM';
elseif ifleray
    reg_str = 'Leray';
elseif ifefr
    reg_str = 'EFR';
elseif iftr
    reg_str = 'TR';
else
    reg_str = 'GROM';
end

if ifleray
   radius = 0.01; 
   fprintf('The filter radius used in the L-ROM is: %f\n', radius);
elseif ifefr || iftr
   radius = 0.01; 
   relax = dt; 
   fprintf('The filter radius and relaxation used in the ROM are: %f, %f\n', radius, relax);
end

if_run_tests=1;

%% Point selection algorithm for DEIM
ps_algs = ["sopt", "gpode", "gappy_pod", "gnat"];
ps_alg = ps_algs(1);

conv_approaches = ["fom", "ftensor", "rtensor", "deim", "clsdeim", "mclsdeim"];
conv_approach = conv_approaches(6);

switch conv_approach
    case 'rtensor'
        ts1 = idivide(int32(nb), int32(2));
        tensor_size = [ts1, ts1, ts1]; 
    case {'deim', 'clsdeim', 'mclsdeim'}
        clsdeim = false;
        ndeim_pts = 200;
        assert(ndeim_pts > 0, 'ndeim_pts must be greater than 0');
        os_multiplier = 2;
        n_os_points = ceil(os_multiplier * ndeim_pts);
end

%% Get the grid and POD bases
reorder = 1; 

cname = fullfile(snaps_path, strcat('bas', casename));
bas_snaps = NekSnaps(cname); 
[pod_u, pod_v] = get_snaps(bas_snaps, reorder);
[x_fom, y_fom] = get_grid(bas_snaps, reorder);
inde = bas_snaps.flds{1}.inde;

%% Define path to dump output in
casedir = sprintf('output/%s_%s/', casename, datestr(now, 'yyyy-mm-dd-HH-MM-SS'));
mkdir(casedir);
basepath = fullfile(casedir, 'fields', casename);
if ~exist(fullfile(casedir, 'fields'), 'dir'), mkdir(fullfile(casedir, 'fields')); end

%% Get the non-linear snapshots and calculate the DEIM points
if ismember(conv_approach, ["deim", "clsdeim", "mclsdeim"])
    matlab_pod_basis = 0;

    if matlab_pod_basis || conv_approach == "mclsdeim"
        nl_cname = fullfile(snaps_path, strcat('csn', casename));
        nl_snaps_obj = NekSnaps(nl_cname);        
        [nl_snaps_u, nl_snaps_v] = get_snaps(nl_snaps_obj, reorder);
    else
        nl_snaps_u = [];
        nl_snaps_v = [];
    end
    
    if matlab_pod_basis
        [nl_bas, ~, ~] = get_pod_basis_from_arrays(nl_snaps_u, nl_snaps_v, x_fom, y_fom, ndeim_pts, 0, 0);
    else
        nl_bas_cname = fullfile(snaps_path, strcat('cba', casename));
        nl_bas_obj_nr = NekSnaps(nl_bas_cname);
        [nl_bas_u_nr, nl_bas_v_nr] = get_snaps(nl_bas_obj_nr, reorder);
        nl_bas = [nl_bas_u_nr; nl_bas_v_nr];
    end

    deim_data = setup_conv_deim(pod_u, pod_v, nl_bas, nl_snaps_u, nl_snaps_v, x_fom, y_fom, ndeim_pts, n_os_points, ps_alg);
end

[au_full, bu_full, cu_full, u0_full, uk_full, mb, ns] = load_full_ops(fullfile(path, 'ops'));

%% Generate the FOM mass matrix
Me = get_Me(x_fom, y_fom);

%% Call validation tests
% Refactored to pass explicit variables rather than relying on global workspace
if if_run_tests
    run_tests(snaps_path, casename, nb, reorder, pod_u, pod_v, au_full, bu_full);
end;

% Note that these are in the original ordering from the Nek5000 simulation
[au, a0, bu, cu, c0, c1, c2, c3, u0, uk, ukmin, ukmax] = get_r_dim_ops(au_full, bu_full, cu_full, u0_full, uk_full, nb);

%% Initialization
time   = 0;
rhs    = zeros(nb, 1);
ext    = zeros(nb, 3);
hufac  = [];

% Preallocate arrays for performance
num_outputs = floor(nsteps / iostep);
ucoef = zeros(num_outputs, nb+1);
kes = zeros(num_outputs, 1);
momentums = zeros(num_outputs, 2);
io_idx = 1;

if ifleray || ifefr || iftr
   dfHfac = [];
   dfHfac = set_df(au, bu, radius, 1, dfHfac);
end

%% Begin integrate ROM with BDFk/EXTk
u = zeros(nb+1, 3); 
u(:,1) = u0;
[alphas, betas] = setcoef();

u_proj = pod_u(:, 1:nb+1) * u(:, 1);
v_proj = pod_v(:, 1:nb+1) * u(:, 1);

field_data = struct('u', u_proj, 'v', v_proj, 'x', x_fom, 'y', y_fom, 'inde', inde, 'size', size(x_fom), 'time', 0.0, 'iostep', 0);
output_fields(basepath, field_data, ifvort, ifwrite, ifvis); 

for istep = int32(1:nsteps)
    time = double(istep) * dt;
    ito = min(istep, 3);
    
    if istep <= 3
        hufac = [];
    end

    rhs = zeros(nb, 1);
    ext(:,3) = ext(:,2);
    ext(:,2) = ext(:,1);
    ext(:,1) = zeros(nb, 1);

    if ifleray
        utmp = [1; (dfHfac \ (dfHfac' \ u(2:end, 1)))];
    else 
        utmp = u;
    end

    % Convection approaches
    switch conv_approach 
        case 'fom'
            c_coef = conv_fom(u(:,1), pod_u, pod_v, x_fom, y_fom, true);
        case 'ftensor'
            c_coef = (reshape(c0*utmp(:,1), nb, nb+1) * u(:,1));
        case 'rtensor'
            c_coef = conv_tensor_dense(u(:,1), pod_u, pod_v, x_fom, y_fom, tensor_size);
        case {'deim', 'clsdeim', 'mclsdeim'}
            c_coef = conv_deim(u(:,1), deim_data, conv_approach);
        otherwise
            error('Unrecognized conv_approach');
    end

    ext(:,1) = ext(:,1) - c_coef;
    ext(:,1) = ext(:,1) - nu * a0;

    if iftr
        utmp = [1; (dfHfac \ (dfHfac' \ u(2:end, 1)))];
        ext(:,1) = ext(:,1) - relax * (u(2:end, 1) - utmp(2:end));
    end

    rhs = rhs + ext * alphas(:, ito);
    rhs = rhs - bu * (u(2:end, :) * betas(2:end, ito)) / dt;

    if ifcopt
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');
        [x, fval, exitflag, output] = fmincon(@(x)rom_residual(x, au, bu, nu, betas, dt, ito, rhs), ...
            u(2:end, 1), [], [], [], [], ukmin(2:end), ukmax(2:end), [], options);
        u_new = [1; x];
    else
        if isempty(hufac)
            h = bu * betas(1, ito) / dt + au * nu;
            hfac = chol(h);
        end
        u_new = [1; (hfac \ (hfac' \ rhs))];
    end

    if ifefr
        utmp = [1; (dfHfac \ (dfHfac' \ u_new(2:end)))];
        u_new = (1 - relax) * u_new + relax * utmp;
    end
        
    u = shift(u, u_new, 3);

    if any(isnan(u(:,1)))
        fprintf('NaN detected at step %d. Aborting.\n', istep);
        break;
    end

    if mod(istep, iostep) == 0
        fprintf('IOSTEP = %d\n', istep); 
        ucoef(io_idx, :) = u(:, 1)';

        u_proj = pod_u(:, 1:nb+1) * u(:, 1);
        v_proj = pod_v(:, 1:nb+1) * u(:, 1);

        ke = 0.5 * (u_proj' * (Me .* u_proj) + v_proj' * (Me .* v_proj));
        kes(io_idx) = ke;

        momentum = [sum(Me .* u_proj), sum(Me .* v_proj)];
        momentums(io_idx, :) = momentum;

        field_data.u = u_proj;
        field_data.v = v_proj;
        field_data.time = time;
        field_data.iostep = idivide(istep, iostep);

        output_fields(basepath, field_data, ifvort, ifwrite, ifvis);
        io_idx = io_idx + 1;
    end
end

%% Write outputs and generate plots
fileID = fopen(fullfile(casedir, "ucoef"), 'w');
fprintf(fileID, "%24.15e\n", ucoef'); 
fclose(fileID);

figure(2);
plot(kes, 'LineWidth', 1.5)
xlabel("Time")
ylabel("Kinetic energy")
title("Energy Conservation")
exportgraphics(gca, fullfile(casedir, "ke.pdf"), 'ContentType', 'vector');

figure(3)
plot(momentums(:,1), 'LineWidth', 1.5); hold on;
plot(momentums(:,2), 'LineWidth', 1.5);
xlabel("Time")
ylabel("Momentum components")
title("Momentum Conservation");
legend('U Momentum', 'V Momentum');
exportgraphics(gca, fullfile(casedir, "momentum.pdf"), 'ContentType', 'vector');

disp("Simulation complete.");

%#####################################
%# Auxiliary functions
%#####################################

function hfac = set_df(a, b, dfRadius, dfOrder, hfac)
    if isempty(hfac)
        bfac = chol(b);
        h = (dfRadius^2) * (bfac \ (bfac' \ a));
        for i = 2:dfOrder
            h = h * ((dfRadius^2) * (bfac \ (bfac' \ a)));
        end
        h = h + eye(size(a));
        hfac = chol(h);
    end
end

function [alphas, betas] = setcoef()
    alphas = zeros(3, 3);
    betas  = zeros(4, 3);

    alphas(1,1) =  1.0;
    alphas(1,2) =  2.0; alphas(2,2) = -1.0;
    alphas(1,3) =  3.0; alphas(2,3) = -3.0; alphas(3,3) =  1.0;

    betas(1,1) =  1.0; betas(2,1) = -1.0;
    betas(1,2) =  1.5; betas(2,2) = -2.0; betas(3,2) =  0.5;
    betas(1,3) =  11.0/6; betas(2,3) = -3.0; betas(3,3) =  1.5; betas(4,3) = -1.0/3;
end

function a = shift(a, b, n)
    for i = n:-1:2 
        a(:,i) = a(:,i-1); 
    end
    a(:,1) = b;
end

function F = rom_residual(x, a, b, diff, betas, dt, ito, rhs)                                                                                                                                                                                                                                                                                             
    h = b * betas(1, ito) / dt + a * diff;
    F1 = h * x - rhs;
    F = norm(F1);
end
