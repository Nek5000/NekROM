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
% Needs to be string, not a character array.
cases = ["ldc", "cyl", "shear", "t2d"];
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

%% IO parameters
ifvort = true; % Whether or not to calculate the vorticity on an IO step and store as temperature
ifwrite = true; % Whether or not to write field files on an IO step
ifvis = false; % Whether or not to visualize a 2D field in Matlab on an IO step

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

%% Point selection algorithm for DEIM
ps_algs = ["sopt", "gpode", "gappy_pod", "gnat"];
ps_alg = ps_algs(1);

conv_approaches = ["fom", "ftensor", "rtensor", "deim", "clsdeim"];
conv_approach = conv_approaches(5);

switch conv_approach
    case 'fom'

    case 'ftensor'

    case 'rtensor'
        ts1 = idivide(nb,int32(2));
        tensor_size = [ts1, ts1, ts1]; 
    case {'deim', 'clsdeim'}
        clsdeim = false;
        if conv_approach == 'clsdeim';
            % TODO: Make this a string rather than a boolean
            % to handle multiple methods
            clsdeim = true;
        end;
        ndeim_pts = 200;
        assert(ndeim_pts > 0);
        os_multiplier = 2;
        % Number of oversample points
        n_os_points=ceil(os_multiplier*ndeim_pts);
    otherwise
        error("Unrecogized convection operator approach");
end

%% Get the grid and POD bases for plotting purposes
% Nek5000 doesn't guarantee the elements are in lexicographical
% order, 
% so if comparing basis vectors from MATLAB and Fortran
% this may be necessary. It sorts puts the elements in
% lexicographical order.
% Note: Perhaps the inde array the NekToolKit produces could
% help with this.
% TODO: Support re-arranging elements based on an inde array
reorder = 1; 

% Load the grid and the snapshots 
cname=strcat(snaps_path,strcat('bas',casename));
bas_snaps = NekSnaps(cname);
[pod_u, pod_v] = get_snaps(bas_snaps,0);
[x_fom, y_fom] = get_grid(bas_snaps,0);
inde = bas_snaps.flds{1}.inde;

%% Define path to dump output in
%% Output results
%if ndeim_pts > 0;
%    if clsdeim
%        clsdeimstr='clsdeim';
%    else
%        clsdeimstr='';
%    end;
%    casedir= sprintf('%s_nb%d_results_ndeim_pts%d_%s_%s',casename,nb,ndeim_pts,clsdeimstr,reg_str)
%else
%    casedir= sprintf('%s_nb%d_results_%s',casename,nb,reg_str)

%end;


casedir = sprintf('%s_%s', casename, datestr(now, 'yyyy-mm-dd-HH-MM-SS/'));
mkdir(casedir);
basepath = strcat(casedir, 'fields/', casename);
%logfile = fopen(strcat(casedir,'logfile', 'wt'));

%fprintf(logfile, 'nb = %d', nb);


%% Test writing field
%write_field(sprintf('%s_rom_/%s',casename,casename), inde, x_fom, y_fom, pod_u(:,1), pod_v(:,1), 0.0, 0)
%exit; 

%% Get the non-linear snapshots and calculate the DEIM points
if conv_approach == "deim" || conv_approach == "clsdeim";

    matlab_pod_basis = 0;
    % Needed for CLSDEIM and computing POD basis in MATLAB
    % Seems like a flaw in CLSDEIM to require loading all of
    % the snapshots. Is there another way to do this?

    if matlab_pod_basis || conv_approach == "clsdeim"
        nl_cname = strcat(snaps_path,strcat('csn',casename));
        nl_snaps_obj = NekSnaps(nl_cname);        
        [nl_snaps_u, nl_snaps_v] = get_snaps(nl_snaps_obj,reorder);
        %nl_snaps = [nl_snaps_u; nl_snaps_v];
    else
        nl_snaps_u = [];
        nl_snaps_v = [];
    end;
    
    if matlab_pod_basis
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

%% Generate the FOM mass matrix.
% Used to calculate the energy
Me = get_Me(x_fom, y_fom);

%% Can call tests here if desired
% tests

% Note that these are in the original ordering from the Nek5000 simulation
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


%% Begin integrate ROM with BDFk/EXTk
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
%vort   = lcurl(reshape(u_proj,size(x_fom)), reshape(v_proj,size(x_fom)), x_fom, y_fom);

field_data = struct('u', u_proj, 'v', v_proj, 'x', x_fom, 'y', y_fom, 'inde', inde, 'size', size(x_fom), 'time', 0.0, 'iostep', 0);

output_fields(basepath, field_data, ifvort, ifwrite, ifvis); 

%write_field(basepath, inde, struct('x', x_fom, 'y', y_fom, 'u', u_proj, 'v', v_proj, 't', vort), size(x_fom), 0.0, 0);

for istep=int32(1:nsteps);
   %disp(istep)
   time = istep*dt;

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

    switch conv_approach 
        case 'fom'
            c_coef = conv_fom(u(:,1), pod_u, pod_v, x_fom, y_fom);
        case 'ftensor'
            c_coef = (reshape(c0*utmp(:,1),nb,nb+1)*u(:,1));
            %c_coef' %ext(:,1)=ext(:,1)-reshape(cu*utmp(:,1),nb,nb+1)*u(:,1);
        case 'rtensor'
            c_coef = conv_tensor_dense(u(:,1), pod_u, pod_v, x_fom, y_fom, tensor_size);
        case {'deim', 'clsdeim'}
            c_coef = conv_deim(u(:,1), pod_u, pod_v, nl_bas, nl_snaps_u, nl_snaps_v, x_fom, y_fom,  ndeim_pts,istep,clsdeim,n_os_points,ps_alg);
        otherwise
            error('Unrecognized conv_approach');
    end;

    %{
   if ndeim_pts == 0;
    if 0;
        % NekROM convection tensor version
        %c_coef
        %exit;
    else
        % MATLAB approach
        % Pseudo ROM version
        % This should be identical to the above.

        %c_coef = conv_tensor(u(:,1), pod_u, pod_v, x_fom, y_fom);


        %c_coef = conv_tensor_dense(u(:,1), pod_u, pod_v, x_fom, y_fom);
        %c_coef = conv_tensor_reduced(u(:,1), pod_u, pod_v, x_fom, y_fom);

        %c_coef = (conv_fom(u(:,1), pod_u, pod_v, x_fom, y_fom));

        %norm(c_coef - c_coef1)/norm(c_coef)
        %exit;
        %norm(pod_u(:,1:nb+1)*c_coef)
    end;
   else;
    % DEIM version
    % Should be close, but not identical to, the above

    
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
    %}

   %norm(pod_u(:,1:nb+1)*c_coef)
   % Isn't ext(:,1) zero at this point?
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

   if any(isnan(u(:,1)));
      break;
   end;

   if (mod(istep,iostep) == 0);
      disp(sprintf('IOSTEP = %d\n', istep)); 
      ucoef(istep/iostep,:)=u(:,1);
      %u(:,1)
      %u(2,1)
      %u(3,1)

      % Calculate quantities of interest
      u_proj = pod_u(:,1:nb+1)*u(1:end,1);
      v_proj = pod_v(:,1:nb+1)*u(1:end,1);
      data = struct('u', u_proj, 'v', v_proj);

      %ke = 0.5*[u_proj; v_proj]'*([Me;Me].*[u_proj; v_proj]) 
      ke = 0.5*(u_proj'*(Me.*u_proj) + v_proj'*(Me.*v_proj));
      kes = [kes; ke];

      momentum = [sum(Me.*u_proj), sum(Me.*v_proj)]
      momentums = [momentums;momentum];

      % Update the field struct and write it out
      field_data.u = u_proj;
      field_data.v = v_proj;
      field_data.time = time;
      field_data.iostep = idivide(istep,iostep);

      output_fields(basepath, field_data, ifvort, ifwrite, ifvis);

   end
end



% Dump out ucoef in casedir
fileID = fopen(casedir+"/ucoef",'w');
fprintf(fileID,"%24.15e\n",ucoef);
fclose(fileID);
%fclose(logfile);

% Plot quantities of interest (here momentum and KE)
% Does the FOM even conserve these? It would depend
% on the time-stepper right?
figure(2);
plot(kes)
xlabel("Time")
ylabel("Kinetic energy")
title("Momentum and Energy Conserving")
ax = gca;
exportgraphics(ax, strcat(casedir, "ke.pdf"), 'ContentType', 'vector');

figure(3)
plot(momentums(:,1)); hold on;
plot(momentums(:,2))
xlabel("Time")
ylabel("Momentum components")
title("Momentum and Energy Conservation");
ax = gca;
exportgraphics(ax, strcat(casedir, "momentum.pdf"), 'ContentType', 'vector');

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

