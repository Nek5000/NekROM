%#######################################################
%
% Matlab/Octave Driver for Galerkin-based Reduced Order Model
%
% Ping-Hsuan Tsai
% 2024-09-05
%
% Cleaned version with tests separated
%
%#######################################################

clear all; close all;

addpath('./point_generators');
addpath('./io');
addpath('./operators');
addpath('./tests');

%% ----------------------------------------------------
%% USER PARAMETERS
%% ----------------------------------------------------

cases = ["ldc","cyl","shear","t2d"];
thiscase = cases(3);

run_tests = false;

switch thiscase
    case "ldc"
        path='../../examples/ldc/';
        snaps_path=[path,'snaps/'];
        casename='ldc';

        nsteps = 10*1e5;
        dt     = 1e-3;
        iostep = 1000;
        nu     = 1./15000;
        nb     = 30;

    case "cyl"
        path='../../examples/cyl/';
        snaps_path=[path,'snaps/'];
        casename='cyl';

        nsteps = 10*1.25e5;
        dt     = 4e-3;
        iostep = 500;
        nu     = 0.01;
        nb     = 20;

    case "shear"
        path='../../examples/shear/';
        snaps_path=[path,'snaps/'];
        casename='shear';

        nsteps = 2*4000;
        dt     = 1e-3;
        iostep = 100;
        nu     = 1/1000;
        nb     = 20;

    case "td2"
        path='../../examples/t2d/';
        casename='t2d';

        nsteps = 800000;
        dt=0.002;
        iostep=100;
        nu=0.0001;
        nb=3;

    otherwise
        error("Unhandled case name")
end

%% ----------------------------------------------------
%% IO OPTIONS
%% ----------------------------------------------------

ifvort  = true;
ifwrite = true;
ifvis   = false;

%% ----------------------------------------------------
%% ROM STABILIZATION
%% ----------------------------------------------------

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

%% ----------------------------------------------------
%% CONVECTION APPROACH
%% ----------------------------------------------------

conv_approaches = ["fom","ftensor","rtensor","deim","clsdeim","mclsdeim"];
conv_approach   = conv_approaches(1);

ps_alg = "sopt";

switch conv_approach
    case 'rtensor'
        ts1 = idivide(int32(nb),int32(2));
        tensor_size = [ts1,ts1,ts1];

    case {'deim','clsdeim','mclsdeim'}
        ndeim_pts = 200;
        os_multiplier = 2;
        n_os_points = ceil(os_multiplier*ndeim_pts);
end

%% ----------------------------------------------------
%% LOAD SNAPSHOTS / GRID
%% ----------------------------------------------------

reorder = 1;

cname=[snaps_path,'bas',casename];
bas_snaps = NekSnaps(cname);

[pod_u,pod_v] = get_snaps(bas_snaps,0);
[x_fom,y_fom] = get_grid(bas_snaps,0);

inde = bas_snaps.flds{1}.inde;

%% ----------------------------------------------------
%% OUTPUT DIRECTORY
%% ----------------------------------------------------

mkdir('output');

casedir = sprintf('output/%s_%s',casename,datestr(now,'yyyy-mm-dd-HH-MM-SS/'));
mkdir(casedir);

basepath=[casedir,'fields/',casename];

%% ----------------------------------------------------
%% SETUP DEIM (if used)
%% ----------------------------------------------------

if any(strcmp(conv_approach,["deim","clsdeim","mclsdeim"]))

    nl_bas_cname=[snaps_path,'cba',casename];
    nl_bas_obj_nr=NekSnaps(nl_bas_cname);

    [nl_bas_u_nr,nl_bas_v_nr]=get_snaps(nl_bas_obj_nr,reorder);
    nl_bas=[nl_bas_u_nr;nl_bas_v_nr];

    deim_data = setup_conv_deim( ...
        pod_u,pod_v,nl_bas,[],[],x_fom,y_fom, ...
        ndeim_pts,n_os_points,ps_alg);

end

%% ----------------------------------------------------
%% LOAD OPERATORS
%% ----------------------------------------------------

[au_full,bu_full,cu_full,u0_full,uk_full,mb,ns] = ...
    load_full_ops([path,'ops']);

Me = get_Me(x_fom,y_fom);

%% ----------------------------------------------------
%% OPTIONAL TESTS
%% ----------------------------------------------------

if run_tests
    run_rom_tests(snaps_path,casename,nb,x_fom,y_fom,pod_u,pod_v,au_full,bu_full)
end
% ----------------------------------------------------
%% REDUCED OPERATORS
%% ----------------------------------------------------

[au,a0,bu,cu,c0,c1,c2,c3,u0,uk,ukmin,ukmax] = ...
get_r_dim_ops(au_full,bu_full,cu_full,u0_full,uk_full,nb);

%% ----------------------------------------------------
%% INITIALIZE TIME INTEGRATION
%% ----------------------------------------------------

time=0;
rhs=zeros(nb,1);
ext=zeros(nb,3);

ucoef=zeros((nsteps/iostep),nb+1);

u=zeros(nb+1,3);
u(:,1)=u0;

[alphas,betas]=setcoef();

kes=[];
momentums=[];

%% ----------------------------------------------------
%% INITIAL FIELD OUTPUT
%% ----------------------------------------------------

u_proj = pod_u(:,1:nb+1)*u(:,1);
v_proj = pod_v(:,1:nb+1)*u(:,1);

field_data = struct( ...
'u',u_proj,'v',v_proj, ...
'x',x_fom,'y',y_fom, ...
'inde',inde, ...
'size',size(x_fom), ...
'time',0.0, ...
'iostep',0);

output_fields(basepath,field_data,ifvort,ifwrite,ifvis);

%% ----------------------------------------------------
%% MAIN TIME LOOP
%% ----------------------------------------------------

for istep=int32(1:nsteps)

    time=istep*dt;
    ito=min(istep,3);

    rhs=zeros(nb,1);

    ext(:,3)=ext(:,2);
    ext(:,2)=ext(:,1);
    ext(:,1)=0;

    switch conv_approach

        case 'fom'
            c_coef = conv_fom(u(:,1),pod_u,pod_v,x_fom,y_fom,true);

        case 'rtensor'
            c_coef = conv_tensor_dense( ...
                u(:,1),pod_u,pod_v,x_fom,y_fom,tensor_size);

        case {'deim','clsdeim','mclsdeim'}
            c_coef = conv_deim(u(:,1),deim_data,conv_approach);

        otherwise
            error("Unknown convection approach")

    end

    ext(:,1)=ext(:,1)-c_coef;
    ext(:,1)=ext(:,1)-nu*a0;

    rhs = rhs + ext*alphas(:,ito);
    rhs = rhs - bu*(u(2:end,:)*betas(2:end,ito))/dt;

    if istep==1
        h = bu*betas(1,ito)/dt + au*nu;
        hfac = chol(h);
    end

    u_new=[1;(hfac\(hfac'\rhs))]';

    u=shift(u,u_new,3);

    if mod(istep,iostep)==0

        disp(["IOSTEP=",num2str(istep)])

        ucoef(istep/iostep,:)=u(:,1);

        u_proj = pod_u(:,1:nb+1)*u(:,1);
        v_proj = pod_v(:,1:nb+1)*u(:,1);

        ke = 0.5*(u_proj'*(Me.*u_proj)+v_proj'*(Me.*v_proj));
        kes=[kes;ke];

        momentum=[sum(Me.*u_proj),sum(Me.*v_proj)];
        momentums=[momentums;momentum];

        field_data.u=u_proj;
        field_data.v=v_proj;
        field_data.time=time;
        field_data.iostep=idivide(istep,iostep);

        output_fields(basepath,field_data,ifvort,ifwrite,ifvis);

    end

end

%% ----------------------------------------------------
%% SAVE RESULTS
%% ----------------------------------------------------

fileID=fopen([casedir,"/ucoef"],'w');
fprintf(fileID,"%24.15e\n",ucoef);
fclose(fileID);

figure(1)
plot(kes)
xlabel("Time")
ylabel("Kinetic Energy")

figure(2)
plot(momentums(:,1)); hold on
plot(momentums(:,2))
xlabel("Time")
ylabel("Momentum")

disp("Finished ROM run")

%% ====================================================
%% AUXILIARY FUNCTIONS
%% ====================================================

function [alphas,betas]=setcoef()

alphas=zeros(3,3);
betas=zeros(4,3);

alphas(:,1)=[1;0;0];
alphas(:,2)=[2;-1;0];
alphas(:,3)=[3;-3;1];

betas(:,1)=[1;-1;0;0];
betas(:,2)=[1.5;-2;0.5;0];
betas(:,3)=[11/6;-3;1.5;-1/3];

end

function a=shift(a,b,n)

for i=n:-1:2
    a(:,i)=a(:,i-1);
end

a(:,1)=b;

end
