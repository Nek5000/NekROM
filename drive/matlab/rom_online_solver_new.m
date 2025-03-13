%#######################################################
%
%# Matlab Driver for Galerkin-based reduced order model 
%# v0.0.0
%
%# Ping-Hsuan Tsai
%# 2024-09-05
%
%#######################################################


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

nsteps = 10*1e5;%80000;%1.25000E+05;%20000; 
dt     = 1.000000E-03;%0.001;
iostep = 5*1000;%500;%250;%500;%10;
nu     = 1./15000;%0.01;
nb     = 20;

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

%deims= [1 2 3 4 5 6 7 8 9 10]
%deims= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]
deims= [200];%[2,4,8,16,20 40 80 100 200 400 800 1000]
n_os_points=0;
% number of deim points
for j=1:length(deims)
ndeim_pts = deims(j);
clsdeim = false;

[au_full, bu_full, cu_full, u0_full, uk_full, mb, ns] = load_full_ops('../../examples/ldc/ops');

%size(au_full)

[au, a0, bu, cu, c0, c1, c2, c3, u0, uk, ukmin, ukmax] = get_r_dim_ops(au_full, bu_full, cu_full, u0_full, uk_full, nb);

%size(uk_full)
%size(uk)
%exit;

%size(au)
%exit;

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
cname='../../examples/ldc/basldc'; % Modify this for your case
%avg_cname='../avgcyl';
bas_snaps = NekSnaps(cname);
[pod_u, pod_v, x_fom, y_fom] = get_grid_and_pod(bas_snaps);

%% Get the non-linear snapshots and calculate the DEIM points
nl_cname = '../../examples/ldc/csnldc'; % Modify this for your case
nl_snaps = NekSnaps(nl_cname);

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
%  c_coef' %ext(:,1)=ext(:,1)-reshape(cu*utmp(:,1),nb,nb+1)*u(:,1);

   % Pseudo ROM version
   % This should be identical to the above.
%  c_coef = (conv_fom(u(:,1), pod_u, pod_v, bas_snaps))
   %norm(pod_u(:,1:nb+1)*c_coef)
   else;
   % DEIM version
   % Should be close, but not identical to, the above
    c_coef = conv_deim(u(:,1), pod_u, pod_v, nl_snaps,ndeim_pts,istep,clsdeim,n_os_points)-c1+c2*utmp(:,1)+c3*utmp(:,1);
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
dt     = 1.000000E-03;%0.001;
iostep = 100;%500;%250;%500;%10;
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

   bool_plot = true;
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

%size(ucoef)
%ucoef
%exit;

%plot(ucoef(:,2)); hold on;
%plot(flip(uk(2,:)));
%pause(60);

if ndeim_pts > 0;
    if clsdeim
        clsdeimstr='clsdeim';
    else
        clsdeimstr='';
    end;
    casedir= sprintf('nb%d_results_ndeim_pts%d_%s_%s',nb,ndeim_pts,clsdeimstr,reg_str)
else
    casedir= sprintf('nb%d_results_%s',nb,reg_str)
end;
mkdir(casedir);

% Dump out ucoef in casedir
fileID = fopen(casedir+"/ucoef",'w');
fprintf(fileID,"%24.15e\n",ucoef);
fclose(fileID);
end

%#####################################
%
%# Auxiliary functions
%# v0.0.0
%
%# Ping-Hsuan Tsai
%# 2024-09-05
%
%#####################################

function [a0_full, b0_full, c0_full, u0_full, uk_full, mb, ms] =  load_full_ops(path)
% load_full_ops function loads the full ROM operators stored in the specified path
% specified by input variable path.
%
% Output:
% a0_full : Full stiffness matrix of size mb+1 x mb+1 (The +1 comes from the zeroth mode)
% b-1_full : Full   mass    matrix of size mb+1 x mb+1
% cu_full : Full advection tensor of size mb x mb+1 x mb+1
% u0_full : Vector of size mb+1, which contains the ROM coefficients of the
%           projection of initial conditons onto mb-reduced space
% uk_full : Matrix of size mb+1 x ns. Each column contains the ROM coefficients
%           of the projection of the snapshot onto mb-reduced space
% mb : total number of modes
% ns : nubmer of snapshots used to create those operators
% 
% This function can take times to load operators with nb >= 300.

   fprintf('Loading ROM operators and vectors... \n');
   fprintf('Currently only support velocity... \n');

   mb=dlmread(fullfile(path,"nb"));

   % load stiffness matrix
   a0_full = dlmread(fullfile(path,"au"));
   a0_full = reshape(a0_full,mb+1,mb+1);

   % load mass matrix
   b0_full = dlmread(fullfile(path,"bu"));
   b0_full = reshape(b0_full,mb+1,mb+1);

   % load advection tensor
   c0_full = dlmread(fullfile(path,"cu"));
   c0_full = reshape(c0_full,mb,mb+1,mb+1);

   u0_full = dlmread(fullfile(path,"u0"));

   ms = dlmread(fullfile(path,"ns"));
   uk_full = dlmread(fullfile(path,"uk"));
   uk_full = reshape(uk_full,mb+1,ms);

   fprintf("done loading ... \n");
   
end

function [a, a0, b, c, c0, c1, c2, c3, u0, uk, ukmin, ukmax] =  get_r_dim_ops(au_full, bu_full, cu_full, u0_full, uk_full, nb)
   index  = [1:nb+1];
   index1 = [1:nb];
   index2 = [2:nb+1];
   a      = au_full(index2,index2);
   a0     = au_full(index2,1);
   b      = bu_full(index2,index2);

   cutmp  = cu_full(index1,index,index);
   c0      = reshape(cutmp,nb*(nb+1),nb+1);
   c1 = cutmp(:,1,1);                                                       
   c2 = reshape(cutmp(:,1,:),nb,nb+1);                                      
   c3 = reshape(cutmp(:,:,1),nb,nb+1); 

   c      = cu_full(index1,index2,index2);
   c      = reshape(c,nb*(nb),nb);

   u0      = u0_full(index);

   uk = uk_full(index,:);
   ukmin = min(uk,[],2);
   ukmax = max(uk,[],2);
end

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


% This needs to do the same thing as
% reshape(cu*utmp(:,1),nb,nb+1)*u(:,1);
function [out_coef] = conv_fom(ucoef, pod_u, pod_v, snaps)
    x=snaps.flds{1}.x;
    y=snaps.flds{1}.y;
    nx1 = size(x,1);
    [zi, w] = zwgll(nx1-1);
    d = deriv_mat(zi);
    [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
    lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
    nL = prod(size(x));
    Me = reshape(jac.*(w*w'),nL,1);

    if false; % Normal way of calculating the gradient
        [ux_fom, uy_fom] = lgrad(u_fom, 0);
        [vx_fom, vy_fom] = lgrad(v_fom, 0);
    else
        % Kento's ROM approach. Calculate the gradients of the POD modes
        ux_pods = [];
        uy_pods = [];
        vx_pods = [];
        vy_pods = [];
        for i = 1:size(pod_u,2);
            [ux_pod, uy_pod] = lgrad(reshape(pod_u(:,i),size(x)),0);
            [vx_pod, vy_pod] = lgrad(reshape(pod_v(:,i),size(x)),0);
            ux_pods = [ux_pods, reshape(ux_pod, nL,1)];
            uy_pods = [uy_pods, reshape(uy_pod, nL,1)];
            vx_pods = [vx_pods, reshape(vx_pod, nL,1)];
            vy_pods = [vy_pods, reshape(vy_pod, nL,1)];
        end;
        ux_fom = reshape(ux_pods*ucoef, size(x));
        uy_fom = reshape(uy_pods*ucoef, size(x));
        vx_fom = reshape(vx_pods*ucoef, size(x));
        vy_fom = reshape(vy_pods*ucoef, size(x));
    end

    u_fom = reshape(Me.*pod_u*ucoef, size(x));
    v_fom = reshape(Me.*pod_v*ucoef, size(x));

    conv_u_fom = reshape(u_fom.*ux_fom + v_fom.*uy_fom, nL,1);
    conv_v_fom = reshape(u_fom.*vx_fom + v_fom.*vy_fom, nL,1);
    
    out_coef = [pod_u(:,2:end); pod_v(:,2:end)]'*[conv_u_fom; conv_v_fom];
end

function[P,indices] = calc_qdeim_proj_mat(U_nl)
    [A,B,P] = qr(U_nl');
    P = P(:,1:size(U_nl,2));
    indices = [];
    for col=1:size(P,2);
        [val, ind] = max(P(:,col));
        indices = [indices, ind];
    end;
end


function[P,indices] = calc_deim_proj_mat(U_nl)
    [maxval, ind] = max(abs(U_nl(:,1)), [], 1);
    [n,m] = size(U_nl);
    P = zeros(n,m);
    indices = zeros(1,m);
    P(ind,1) = 1;
    indices(1,1) = ind;
    for i=2:m;
    % These two ways of calculating c should give the same result
        %c = (P(:,1:i-1)'*U_nl(:,1:i-1)) \ P(:,1:i-1)'*U_nl(:,i)
        c = U_nl(indices(1:i-1),1:i-1) \ U_nl(indices(1:i-1),i);
    r = U_nl(:,i) - U_nl(:,1:i-1)*c;
    [maxval, ind] = max(abs(r), [], 1);
    indices(1,i) = ind;
    P(ind,i) = 1;
    end;
end


