% config.m - Configuration for Galerkin-based reduced order model
% Refactored for Octave Compatibility: 2026
%
% This file defines all simulation parameters. Values can be overridden
% using environment variables (see README.md).
%
% Key Parameters:
%   casename       - Test case: 'ldc', 'cyl', 'shear', 't2d', 'os7000'
%   nsteps         - Total timesteps
%   dt             - Timestep size
%   iostep         - Output frequency (every N steps)
%   nu             - Kinematic viscosity
%   nb             - Number of POD modes (excluding mean)
%   conv_approach  - Convection method: 'fom', 'ftensor', 'rtensor',
%                    'deim', 'clsdeim', 'mclsdeim'
%   ps_alg         - DEIM point selection: 'gpode', 'gappy_pod', 'sopt', 'gnat'
%   ndeim_pts      - Number of DEIM sample points (for DEIM methods)
%   ifvort         - Compute vorticity field
%   ifwrite        - Write field files
%   ifvis          - Real-time visualization (slow)
%   ifplot         - Save energy/momentum plots
%
% Stabilization Options (advanced):
%   ifleray        - Leray regularization (spectral viscosity)
%   ifefr          - Evolve-then-filter
%   iftr           - Tikhonov regularization
%   ifcopt         - Constrained optimization
%
% See README.md for detailed documentation.

%% Case Selections 
% Use Cell Arrays {} instead of String Arrays [] for Octave compatibility
cases = {'ldc', 'cyl', 'shear', 't2d', 'os7000'};
thiscase = cases{3}; % Use curly braces {} to extract string from cell
env_case = getenv('NEKROM_CASE');
if ~isempty(env_case)
    thiscase = lower(strtrim(env_case));
end
casename = thiscase; 

switch thiscase
    case 'ldc'
        case_path = '../../examples/ldc/';
        nsteps = 10 * 1e5; 
        dt     = 1.0e-03;
        iostep = 1000;
        nu     = 1/15000;
        nb     = 30;
    case 'cyl'
        case_path = '../../examples/cyl/';
        nsteps = 10 * 1.25e05; 
        dt     = 4.0e-03;
        iostep = 500;
        nu     = 0.01;
        nb     = 20;
    case 'shear'
        case_path = '../../examples/shear/';
        nsteps = 10 * 4000;
        dt     = 1e-3;
        iostep = 100;
        nu     = 1/1000;
        nb     = 20; 
    case 't2d'
        case_path = '../../examples/t2d/';
        nsteps = 800000;
        dt     = 0.002;
        iostep = 100;
        nu     = 0.0001;
        nb     = 3;  
    case {'os7000', 'u3_t020_n13'}
        case_path = '../../examples/os7000/';
        casename  = 'u3_t020_n13';
        nsteps = 10000;
        dt     = 2.0e-02;
        iostep = 50;
        nu     = 1/7500;
        nb     = 1;
    otherwise
        % Octave error() uses standard formatting
        error(['Unhandled case name: ', thiscase]);
end

case_meta = load_case_metadata(case_path, casename);
if ~isempty(case_meta.nb)
    nb = case_meta.nb;
end
pod_mode0 = case_meta.pod_mode0;
subtract_mean = case_meta.subtract_mean;

env_nsteps = getenv('NEKROM_NSTEPS');
if ~isempty(env_nsteps)
    nsteps = str2double(env_nsteps);
end
env_iostep = getenv('NEKROM_IOSTEP');
if ~isempty(env_iostep)
    iostep = str2double(env_iostep);
end

path = case_path;
snaps_path = [case_path, 'snaps/']; % Standard concatenation

%% IO & Physics Flags
ifvort  = read_env_bool('NEKROM_IFVORT', true);
ifwrite = read_env_bool('NEKROM_IFWRITE', true);
ifvis   = read_env_bool('NEKROM_IFVIS', false);
if_run_tests = read_env_bool('NEKROM_IF_RUN_TESTS', true);
ifplot  = read_env_bool('NEKROM_IFPLOT', true);

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
% Converted to cell arrays
ps_algs = {'sopt', 'gpode', 'gappy_pod', 'gnat'};
ps_alg = ps_algs{3};

conv_approaches = {'fom', 'ftensor', 'rtensor', 'deim', 'clsdeim', 'mclsdeim'};
conv_approach = conv_approaches{6};
deim_finegrid = false;
deim_alpha = 1.0e-12;
if ~isempty(case_meta.conv_approach)
    conv_approach = case_meta.conv_approach;
end
if ~isempty(case_meta.deim_alpha)
    deim_alpha = case_meta.deim_alpha;
end
env_conv_approach = getenv('NEKROM_CONV_APPROACH');
if ~isempty(env_conv_approach)
    conv_approach = lower(strtrim(env_conv_approach));
end
deim_finegrid = read_env_bool('NEKROM_DEIM_FINEGRID', deim_finegrid);
deim_dealias = read_env_bool('NEKROM_DEIM_DEALIAS', false);
deim_dealias_cquad = read_env_bool('NEKROM_DEIM_DEALIAS_CQUAD', false);
deim_dealias_quad = read_env_bool('NEKROM_DEIM_DEALIAS_QUAD', false);
if deim_dealias_quad
    deim_dealias = false;
    deim_dealias_cquad = false;
end
if deim_dealias_cquad
    deim_dealias = false;
end

if deim_dealias_quad
    deim_dealias_mode = 'quad';
elseif deim_dealias_cquad
    deim_dealias_mode = 'cquad';
elseif deim_dealias
    deim_dealias_mode = 'sample';
else
    deim_dealias_mode = 'none';
end

% Validate conv_approach early to catch typos
valid_conv_approaches = {'fom', 'ftensor', 'rtensor', 'deim', 'clsdeim', 'mclsdeim'};
if ~ismember(conv_approach, valid_conv_approaches)
    error('Invalid conv_approach: "%s". Valid options: %s', ...
        conv_approach, strjoin(valid_conv_approaches, ', '));
end

switch conv_approach
    case 'rtensor'
        % Replaced idivide with floor for generic compatibility
        ts1 = floor(nb / 2);
        tensor_size = [ts1, ts1, ts1];
    case {'deim', 'clsdeim', 'mclsdeim'}
        ndeim_pts = 20;
        os_multiplier = 0;
        n_os_points = ceil(os_multiplier * ndeim_pts);

        % Validate ps_alg
        valid_ps_algs = {'sopt', 'gpode', 'gappy_pod', 'gnat'};
        if ~ismember(ps_alg, valid_ps_algs)
            error('Invalid ps_alg: "%s". Valid options: %s', ...
                ps_alg, strjoin(valid_ps_algs, ', '));
        end
end
