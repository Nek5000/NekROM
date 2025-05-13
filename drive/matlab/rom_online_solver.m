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

clear all; close all;

% Add any important scripts to path
addpath('./point_generators');
addpath('./io');

%% Specify the case path and case name

path='../../examples/ldc/';
casename='ldc';

%path='../../examples/conv/';
%casename='cyl';

%path='../../examples/shear4/';
%casename='shear4';%'thin';

% Should just use the values from the .rea or MOR file by default
% allowing for overrides
if contains(path, 'ldc')
    nsteps = 10*1e5;%80000;%1.25000E+05;%20000; 
    dt     = 1.000000E-03;%0.001;
    iostep = 5*1000;%500;%250;%500;%10;
    nu     = 1./15000;%0.01;
    nb     = 20;
elseif contains(path,'conv')
    nsteps = 10*1.25000E+05;%20000; 
    dt     = 4.000000E-03;%0.001;
    iostep = 5000;%250;%500;%10;
    nu     = 0.01;
    nb     = 20;
elseif contains(path,'shear4')
    nsteps = 2*40000;
    dt     = 1e-4;
    iostep = 100;
    nu     = 1/40000;
    nb     = 30;    
else
  disp('Error: unrecognized case');
  exit;
end;

% Whether or not to plot on an iostep
bool_plot = true;

% ROM stabilization strategies
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

%deims= [50];%[170];%[2,4,8,16,20 40 80 100 200 400 800 1000]
%% Point selection algorithm
%ps_alg='sopt';
ps_alg ='gpode';
%ps_alg = 'gappy_pod';
%ps_alg = 'gnat';

%% Hyperreduction algorithms
hr_alg="clsdeim";
clsdeim = false;

% number of deim points
%for j=1:length(deims)
ndeim_pts = 256;%deims(j);
os_multiplier = 2;
n_os_points=os_multiplier*ndeim_pts;

[au_full, bu_full, cu_full, u0_full, uk_full, mb, ns] = load_full_ops(strcat(path,'ops'));
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

%% Get the grid and POD bases for plotting purposes
cname=strcat(path,strcat('bas',casename));
%avg_cname='../avgcyl';
bas_snaps = NekSnaps(cname);
[pod_u, pod_v] = get_snaps(bas_snaps);
[x_fom, y_fom] = get_grid(bas_snaps);

%% Get the non-linear snapshots and calculate the DEIM points
if ndeim_pts > 0;
    nl_cname = strcat(path,strcat('csn',casename));
    nl_snaps = NekSnaps(nl_cname);
end;
%h=bu*betas(1,ito)/dt+au*nu;


% Begin integrate ROM with BDFk/EXTk
u     = zeros(nb+1,3); % vectors for BDF3/EXT3
u(:,1)=u0;
[alphas, betas] = setcoef();


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

   % ROM convection tensor version
   if ndeim_pts == 0;
    c_coef = (reshape(c0*utmp(:,1),nb,nb+1)*u(:,1));
    %c_coef' %ext(:,1)=ext(:,1)-reshape(cu*utmp(:,1),nb,nb+1)*u(:,1);

    % Pseudo ROM version
    % This should be identical to the above.
    %  c_coef = (conv_fom(u(:,1), pod_u, pod_v, bas_snaps))
    %norm(pod_u(:,1:nb+1)*c_coef)
   else;
    % DEIM version
    % Should be close, but not identical to, the above

    c_coef = conv_deim(u(:,1), pod_u, pod_v, x_fom, y_fom, nl_snaps,ndeim_pts,istep,clsdeim,n_os_points,ps_alg);
    c_coef = c_coef - c1+c2*utmp(:,1)+c3*utmp(:,1);
   end;

   %norm(pod_u(:,1:nb+1)*c_coef)
   ext(:,1)=ext(:,1)-c_coef;

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
   time = time+dt;

   if any(isnan(u(:,1)));
      break;
   end;

   if (mod(istep,iostep) == 0);
      ucoef(istep/iostep,:)=u(:,1);
      u(:,1)

      if bool_plot;      
      u_proj = pod_u(:,1:nb+1)*u(1:end,1);
      v_proj = pod_v(:,1:nb+1)*u(1:end,1);
     
      u_abs = sqrt(u_proj.^2 + v_proj.^2);
      %norm(u_abs)
      hold off;
      patch_plot(x_fom,y_fom, reshape(u_abs,size(x_fom)), [], 'PlotType', 'surface');
      pause(0.01);
     end;
   end

end

% Output results
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

