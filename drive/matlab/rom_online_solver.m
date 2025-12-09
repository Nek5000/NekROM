%#######################################################
%
%# Matlab Driver for Galerkin-based reduced order model 
%# v0.0.0
%
%# Ping-Hsuan Tsai
%# 2024-09-05
%
%#######################################################

%TODO: Move user parameters to separate file.

%% User Parameters
% Problem parmeters:
%   nsteps : int
%     The total number of time steps to integrate the ROM ODEs
%   dt     : float
%     Time step size used for integrating the ROM ODEs, could be larger than the FOM dt
%   iostep : int
%     Every iostep you want to store the ROM coefficients
%   nu     : float
%     Viscosity
%   nb     : int 
%     Size of your ROM
%
% Stabilization methods:
%   ifcopt : bool 
%     Set to true if you want to use constrained ROM (C-ROM)
%   ifleray : bool 
%     Set to true if you want to use Leray-ROM (L-ROM)
%   ifefr : bool 
%     Set to true if you want to use evolve-filter-relax ROM (EFR-ROM)
%   iftr : bool 
%     Set to true if you want to use time-relaxataion ROM (TR-ROM)

% Clear memory
clear all; close all;

% Add any important scripts to path
addpath('./point_generators');
addpath('./io');
addpath('./operators')

%% Specify the case path and case name

% Use one of the pre-existing cases or add your own
cases = ['ldc', 'cyl', 'shear', 't2d'];
thiscase = cases(2);

% TODO: Should just use the values from the .rea or MOR file by default
switch thiscase
    case 'ldc'
        path='../../examples/ldc_v2/';
        snaps_path=strcat(path,'snaps/');
        casename='ldc';

        nsteps = 10*1e5;%80000;%1.25000E+05;%20000; 
        dt     = 1.000000E-03;%0.001;
        iostep = 1000;%5*1000;%500;%250;%500;%10;
        nu     = 1./15000;%0.01;
        nb     = 30;
    case 'cyl'
        path='../../examples/cyl/';
        snaps_path=strcat(path,'snaps/');
        casename='cyl';

        nsteps = 10*1.25000E+05;%20000; 
        dt     = 4.000000E-03;%0.001;
        iostep = 500;%250;%500;%10;
        nu     = 0.01;
        nb     = 20;
    case 'shear'
        path='../../examples/shear4/';
        snaps_path=strcat(path,'snaps/');
        casename='nick'%'shear4';%'thin';

        %nsteps = 4000; %Reconstruction
        %nsteps = 8000; % Extrapolation
        nsteps = 10*4000;
        dt     = 1e-3;
        iostep = 100;
        nu     = 1/1000;
        nb     = 30; 
    case 'td2'
        path='../../examples/t2d/';
        casename='t2d';

        nsteps=800000;
        dt=0.002;
        iostep=100;
        nu=0.0001;
        nb=3;  
    otherwise
        error("unhandled case name");
end

% Whether or not to plot on an iostep
bool_plot = true;
% Plot the vorticity, otherwise only calculates the velocity magnitude
plot_vort = true;

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
elseif iftr;
    reg_str = 'TR';
else
    reg_str = 'GROM';
end;   

if (ifleray)
   radius = 0.01; % Feel free to change it
   if exist('radius', 'var')
      fprintf('The filter radius used in the L-ROM is: %f\n', radius);
   end
elseif (ifefr) || (iftr)
   radius = 0.01; % Feel free to change it
   relax = dt; % Feel free to change it
   if exist('radius', 'var') && exist('relax', 'var')
      fprintf('The filter radius and relaxation used in the ROM are: %f %f\n', radius, relax);
   end
end

%% Point selection algorithm
ps_algs = ['sopt', 'gpode', 'gappy_pod', 'gnat'];
ps_alg = ps_algs(1);

%% Hyperreduction algorithms
%TODO support multiple hyperreduction algorithms using a string
%hr_alg="clsdeim";
clsdeim = false;

% number of deim points
ndeim_pts = 128;%100;%100;%400;%10;%800;%400;
% number of oversample points
os_multiplier = 2;
n_os_points=ceil(os_multiplier*ndeim_pts);


%% Get the grid and POD bases for plotting purposes
% NekROM may re-arrange the elements of the basis vectors 
% so if comparing basis vectors from MATLAB and Fortran
% this may be necessary. It sorts puts the elements in
% lexicographical order.
% TODO: Support re-arranging elements based on an inde array
reorder = 1; 

% Load the grid and the snapshots 
cname=strcat(snaps_path,strcat('bas',casename));
bas_snaps = NekSnaps(cname);
[pod_u, pod_v] = get_snaps(bas_snaps,0);
[x_fom, y_fom] = get_grid(bas_snaps,0);
inde = bas_snaps.flds{1}.inde;

%% Define path to dump output in
%basepath = sprintf('%s_rom_snaps_reduced_tensor/%s',casename,casename);
basepath = sprintf('%s_rom_snaps_copt_deim_%i/%s',casename,ndeim_pts,casename);

%% Test writing field
%write_field(sprintf('%s_rom_/%s',casename,casename), inde, x_fom, y_fom, pod_u(:,1), pod_v(:,1), 0.0, 0)
%exit; 

%% Get the non-linear snapshots and calculate the DEIM points
if ndeim_pts > 0;

    % Needed for CLSDEIM and computing POD basis in MATLAB
    % Seems like a flaw in CLSDEIM to require loading all of
    % the snapshots. Is there another way to do this?
    nl_cname = strcat(snaps_path,strcat('csn',casename));
    nl_snaps_obj = NekSnaps(nl_cname);        
    [nl_snaps_u, nl_snaps_v] = get_snaps(nl_snaps_obj,reorder);
    nl_snaps = [nl_snaps_u; nl_snaps_v];
    
    if 0
        % Set to 1 to generate non-linear POD basis in MATLAB
        [nl_bas, ~, ~] = get_pod_basis_from_arrays(nl_snaps_u, nl_snaps_v, x_fom, y_fom, ndeim_pts, 0, 0);
    else
        % Get the NekROM non-linear POD basis
        nl_bas_cname = strcat(snaps_path, strcat('cba',casename));
        nl_bas_obj_nr = NekSnaps(nl_bas_cname);
        [nl_bas_u_nr, nl_bas_v_nr] = get_snaps(nl_bas_obj_nr,reorder);
        nl_bas = [nl_bas_u_nr; nl_bas_v_nr];
    end;
end;

[au_full, bu_full, cu_full, u0_full, uk_full, mb, ns] = load_full_ops(strcat(path,'ops'));

%% Create POD in MATLAB if desired
if 0
    subtract_mean = 1;
    conserve_momentum = 0;
    snaps_obj = NekSnaps(strcat(snaps_path,casename)); % Should the snaps object have reordering capability?
    [pod_ml, u0_full_ml, uk_full_ml] = get_pod_basis(snaps_obj,nb,reorder,subtract_mean,conserve_momentum);

    pod_u_ml = pod_ml(1:size(pod_ml,1)/2,1:nb+1);
    pod_v_ml = pod_ml(size(pod_ml,1)/2 + 1:end,1:nb+1);
end;

% Isn't this the same as bu_full? We can probably avoid calling this.
Me = get_Me(x_fom, y_fom);
npf = size(Me, 1);
Me_vec = spdiags([Me;Me], 0, 2*npf,2*npf);

% Maybe define some unit test and move these tests to them?

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

[au_full_ml, bu_full_ml] = gen_Au(pod_u_ml, pod_v_ml,x_fom, y_fom);

%% Check that MATLAB and Fortran Au and Bu operators are the same
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

if 0; % Test out reconstructing the snapshots
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


% Note that these are in the original ordering
[au, a0, bu, cu, c0, c1, c2, c3, u0, uk, ukmin, ukmax] = get_r_dim_ops(au_full, bu_full, cu_full, u0_full, uk_full, nb);

% Initialize variables
time   = 0.;
rhs    = zeros(nb,1);
ext    = zeros(nb,3);
hufac  = [];
ucoef  = zeros((nsteps/iostep),nb+1);
if (ifleray) || (ifefr) || (iftr)
   dfHfac = [];
   [dfHfac] = set_df(au, bu, radius, 1, dfHfac);
end

%h=bu*betas(1,ito)/dt+au*nu;


% Begin integrate ROM with BDFk/EXTk
u     = zeros(nb+1,3); % vectors for BDF3/EXT3
u(:,1)=u0;
[alphas, betas] = setcoef();

kes = [];
momentums = [];

%Au_ml
%au_full


%Au_ml*u
%au_full*u

%(au_full_ml*u)./(au_full*u)
%(bu_full_ml*u)./(bu_full*u)
%exit;

u_proj = pod_u(:,1:nb+1)*u(1:end,1);
v_proj = pod_v(:,1:nb+1)*u(1:end,1);
vort   = lcurl(reshape(u_proj,size(x_fom)), reshape(v_proj,size(x_fom)), x_fom, y_fom);
write_field(basepath, inde, struct('x', x_fom, 'y', y_fom, 'u', u_proj, 'v', v_proj, 't', vort), size(x_fom), 0.0, 0);

for istep=1:nsteps
   istep
   ito=min(istep,3);
   if istep<= 3
      hufac = [];
   end

   % Compute the right-handed side of the fully discretized system associated to BDFk/EXTk
   rhs = zeros(nb,1);

   ext(:,3)=ext(:,2);
   ext(:,2)=ext(:,1);
   ext(:,1)=ext(:,1)*0;

   if ifleray
      % Compute filtered u
      utmp = [1;(dfHfac\(dfHfac'\u(2:end,1)))];
   else 
      utmp = u;
   end

   %% Compare the NL evaluation results.

   if ndeim_pts == 0;
    if 0;
        % NekROM convection tensor version
        c_coef = (reshape(c0*utmp(:,1),nb,nb+1)*u(:,1));
        %c_coef' %ext(:,1)=ext(:,1)-reshape(cu*utmp(:,1),nb,nb+1)*u(:,1);
        %c_coef
        %exit;
    else
        % MATLAB approach
        % Pseudo ROM version
        % This should be identical to the above.

        %c_coef = conv_tensor(u(:,1), pod_u, pod_v, x_fom, y_fom);


        %c_coef = conv_tensor_dense(u(:,1), pod_u, pod_v, x_fom, y_fom);
        c_coef = conv_tensor_dense(u(:,1), pod_u, pod_v, x_fom, y_fom, [nb/2,nb/2,nb]);
        %c_coef = conv_tensor_reduced(u(:,1), pod_u, pod_v, x_fom, y_fom);

        %c_coef = (conv_fom(u(:,1), pod_u, pod_v, x_fom, y_fom));

        %norm(c_coef - c_coef1)/norm(c_coef)
        %exit;
        %norm(pod_u(:,1:nb+1)*c_coef)
    end;
   else;
    % DEIM version
    % Should be close, but not identical to, the above

    c_coef = conv_deim(u(:,1), pod_u, pod_v, nl_bas, nl_snaps_u, nl_snaps_v, x_fom, y_fom,  ndeim_pts,istep,clsdeim,n_os_points,ps_alg, Me);
    
    % Note that c1, c2, c3 need to be reordered or need to use original order for everything.
    %c_coef = c_coef - c1+c2*utmp(:,1)+c3*utmp(:,1); % Remind me why this is needed? 
    %c1
    %c2
    %c3
    %c2*utmp(:,1)
    %c3*utmp(:,1)
    %exit
    %Can this be incorporated into conv_deim (NJC):
    % (Ping-Hsuan added that, seems to do something with the zeroth modes.)
    %c_coef
   end;

   %norm(pod_u(:,1:nb+1)*c_coef)
   ext(:,1)=ext(:,1)-c_coef;

   %ext(:,1)=ext(:,1)-0*c_coef; % Try turning off convection

   %% End result comparison

   ext(:,1)=ext(:,1)-nu*a0;

   if iftr
      % Compute filtered u
      utmp = [1;(dfHfac\(dfHfac'\u(2:end,1)))];
      ext(:,1) = ext(:,1)-relax*(u(2:end,1)-utmp(2:end));
   end

   rhs=rhs+ext*alphas(:,ito);
   rhs=rhs-bu*(u(2:end,:)*betas(2:end,ito))/dt;

   % Solve the linear system to get next step solution
   if (ifcopt)
%     if any(u(2:end,1) > ukmax(2:end)) || any(u(2:end,1) < ukmin(2:end))
         %fprintf('in constrained %d \n',istep);
         options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
         [x,fval,exitflag,output] = fmincon(@(x)rom_residual(x,au,bu,nu,betas,dt,ito,rhs),u(2:end,1),[],[],[],[],ukmin(2:end),ukmax(2:end),[],options);
         u_new = [1,x'];
%     else
%  nsteps = 80000;%1.25000E+05;%20000; 
%dt     = 1.000000E-03;%0.001;
%iostep = 100;%500;%250;%500;%10;
      % if constraints are satisfied, do normal solve
%        if isempty(hufac)
%           h=bu*betas(1,ito)/dt+au*nu;
%           hfac=chol(h);
%        end
%        u_new = [1,(hfac\(hfac'\rhs))'];
%     end
   else
      if isempty(hufac)
         h=bu*betas(1,ito)/dt+au*nu;
         hfac=chol(h);
      end
      u_new = [1,(hfac\(hfac'\rhs))'];
   end

   % EFR-ROM additional step
   if (ifefr)
      utmp = u_new;
      utmp = [1,(dfHfac\(dfHfac'\u_new(1,2:end)'))'];
      u_new = (1-relax)*u_new + relax*utmp;
   end
        
   u = shift(u,u_new,3);
   time = istep*dt;

   if any(isnan(u(:,1)));
      break;
   end;

   if (mod(istep,iostep) == 0);
      ucoef(istep/iostep,:)=u(:,1);
      u(:,1)
      u(2,1)
      u(3,1)

      % Calculate quantities of interest
      u_proj = pod_u(:,1:nb+1)*u(1:end,1);
      v_proj = pod_v(:,1:nb+1)*u(1:end,1);
      ke = 0.5*[u_proj; v_proj]'*Me_vec*[u_proj; v_proj]
      kes = [kes; ke];
      momentum = [sum(Me.*u_proj), sum(Me.*v_proj)]
      momentums = [momentums;momentum];
      
      if bool_plot;
 
        data = struct('u', u_proj, 'v', v_proj);
        if plot_vort;
            data.t = lcurl(reshape(u_proj,size(x_fom)), reshape(v_proj,size(x_fom)), x_fom, y_fom);%vx-uy;
        end;

        if 0; % Set to 1 to enable plotting in MATLAB. This slows the code considerably though.
            hold off;
            % Surface or contour. Note that contours are not supported on the main branch of NekToolkit
            if plot_vort
                patch_plot(x_fom,y_fom, reshape(data.t,size(x_fom)), [], 'PlotType', 'surface');
            else
                patch_plot(x_fom,y_fom, reshape(u_proj.^2 + v_proj.^2),size(x_fom)), [], 'PlotType', 'surface');
            end;
        end;

        disp(sprintf('Writing output %i', istep));
        write_field(basepath, inde, data, size(x_fom), time, istep/iostep);
     end;
   end
end

%% Output results
if ndeim_pts > 0;
    if clsdeim
        clsdeimstr='clsdeim';
    else
        clsdeimstr='';
    end;
    casedir= sprintf('%s_nb%d_results_ndeim_pts%d_%s_%s',casename,nb,ndeim_pts,clsdeimstr,reg_str)
else
    casedir= sprintf('%s_nb%d_results_%s',casename,nb,reg_str)
end;
mkdir(casedir);

% Dump out ucoef in casedir
fileID = fopen(casedir+"/ucoef",'w');
fprintf(fileID,"%24.15e\n",ucoef);
fclose(fileID);

% Plot quantities of interest (here momentum and KE)
% Does the FOM even conserve these? It would depend
% on the time-stepper right?
figure(2);
plot(kes)
xlabel("Time")
ylabel("Kinetic energy")
title("Momentum and Energy Conserving")
ax = gca;
exportgraphics(ax, "ke.pdf", 'ContentType', 'vector');

figure(3)
plot(momentums(:,1)); hold on;
plot(momentums(:,2))
xlabel("Time")
ylabel("Momentum components")
title("Momentum and Energy Conserving");
ax = gca;
exportgraphics(ax, "momentum.pdf", 'ContentType', 'vector');

disp("Paused")
pause()
pause()
pause()

%end

%#####################################
%
%# Auxiliary functions
%# v0.0.0
%
%# Ping-Hsuan Tsai
%# 2024-09-05
%
%#####################################

function [hfac] = set_df(a,b,dfRadius,dfOrder,hfac)
% Construct the mth order differential filter (df) 
% I + (\delta^2 B^{-1}A)^m
% with radius delta
   if isempty(hfac)
      bfac = chol(b);
      h = (dfRadius^2)*(bfac\(bfac'\a));
      for i=2:dfOrder
         h = h*((dfRadius^2)*(bfac\(bfac'\a)));
      end
      h = h + eye(size(a));
      hfac=chol(h);
   end
end

function [alphas, betas] = setcoef()
% Setup BFDk/EXTk coefficients
   alphas=zeros(3,3);
   betas=zeros(4,3);

   alphas(1,1)=  1.0;

   alphas(1,2)=  2.0;
   alphas(2,2)= -1.0;

   alphas(1,3)=  3.0;
   alphas(2,3)= -3.0;
   alphas(3,3)=  1.0;

   betas(1,1)=  1.0;
   betas(2,1)= -1.0;

   betas(1,2)=  1.5;
   betas(2,2)= -2.0;
   betas(3,2)=  0.5;

   betas(1,3)=  11.0/6;
   betas(2,3)= -3.0;
   betas(3,3)=  1.5;
   betas(4,3)= -1.0/3;
end

function a = shift(a,b,n)
   for i=n:-1:2; 
      a(:,i)=a(:,i-1); 
   end
   a(:,1)=b;
end

function F = rom_residual(x,a,b,diff,betas,dt,ito,rhs)                                                                                                                                                                                            
   h=b*betas(1,ito)/dt+a*diff;
   F1 = h*x-rhs;
   F = norm(F1);
end

