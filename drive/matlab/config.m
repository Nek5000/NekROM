% config.m - Configuration for Galerkin-based reduced order model 
% Refactored: 2024-09-05

%% Case Selection
cases = ["ldc", "cyl", "shear", "t2d"];
thiscase = cases(3); % Select "shear"
casename = char(thiscase);

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

snaps_path = strcat(path, 'snaps/');

%% IO & Physics Flags
ifvort  = true;  
ifwrite = true;  
ifvis   = false; 
if_run_tests = 1;

%% ROM Stabilization Strategies
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

% Filter parameters
radius = 0.01; 
relax = dt; 

%% Convection & DEIM Settings
ps_algs = ["sopt", "gpode", "gappy_pod", "gnat"];
ps_alg = ps_algs(1);

conv_approaches = ["fom", "ftensor", "rtensor", "deim", "clsdeim", "mclsdeim"];
conv_approach = conv_approaches(6);

switch conv_approach
    case 'rtensor'
        ts1 = idivide(int32(nb), int32(2));
        tensor_size = [ts1, ts1, ts1]; 
    case {'deim', 'clsdeim', 'mclsdeim'}
        ndeim_pts = 200;
        os_multiplier = 2;
        n_os_points = ceil(os_multiplier * ndeim_pts);
end
