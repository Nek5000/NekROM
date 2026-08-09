# NekROM MATLAB/Octave Driver

Fast, interactive reduced-order model (ROM) driver for testing and prototyping. Supports multiple convection approximation methods including FOM, tensor, and DEIM variants.

## Features

- **Multiple ROM approaches**: FOM pseudo-ROM, precomputed tensor, runtime tensor, DEIM, CLS-DEIM, MCLS-DEIM
- **Octave compatible**: Works with free/open-source GNU Octave
- **Thermo-fluid support**: Advances coupled velocity + temperature ROMs when `ops/{at,bt,ct,t0,tk}` is present; optional TDEIM via `ops/tdeim_*`
- **Comprehensive diagnostics**: Energy/momentum conservation tracking, NaN detection, stability metrics
- **Automated testing**: Unit tests for DEIM operators, stability comparison tools
- **Example cases**: Lid-driven cavity, cylinder flow, shear layer, Taylor-Green vortex, Orr-Sommerfeld/Poiseuille (plus thermo cases like `ann`, `rb_axi`, `cylbig_abm`)

## Quick Start

### Prerequisites

**Required:**
- **MATLAB** (R2020a or later) or **Octave** (6.0 or later)
- **NekToolKit** - Spectral element utilities library **(MANDATORY)**
  - Repository: https://github.com/kent0/NekToolKit
  - Used for: Snapshot reading, quadrature nodes/weights, derivative operators, gradient computation
  - Driver will fail immediately if NekToolKit functions are not available

**Installation:**
```bash
# Clone NekToolKit as sibling to NekROM (driver auto-detects)
cd /path/to/workspace
git clone https://github.com/kent0/NekToolKit.git

# Directory structure:
#   workspace/
#   ├── NekROM/
#   └── NekToolKit/

# OR manually add to MATLAB path:
addpath('/path/to/NekToolKit/matlab')
```

**Octave users:**
```bash
# Install optional packages (for constrained optimization)
pkg install -forge optim
```

**Verify dependencies:**
```matlab
cd drive/matlab
check_dependencies  % Checks for NekToolKit and all requirements
```

### Running Your First Case

```bash
cd drive/matlab

# Option 1: Run in MATLAB GUI
matlab
>> driver

# Option 2: Run from command line
matlab -batch "driver"

# Option 3: Run with Octave
octave --eval "driver"
```

**Default case**: Shear layer (20 modes, 40,000 timesteps, mclsdeim convection)

Results appear in `output/shear_YYYY-MM-DD-HH-MM-SS/`

## Selecting a Case

Edit `config.m` line 7:
```matlab
cases = {'ldc', 'cyl', 'shear', 't2d', 'os7000', 'ann', 'rb_axi', 'cylbig_abm'};
thiscase = cases{3};  % Change index
```

Or use environment variables:
```bash
export NEKROM_CASE=cyl
matlab -batch "driver"
```

### Case Descriptions

| Case | Description | Re/Ra | ROM modes | Challenge |
|------|-------------|-------|-----------|-----------|
| `ldc` | Lid-driven cavity | Re=15,000 | 30 | Recirculation zones |
| `cyl` | Flow past cylinder | Re=100 | 20 | Vortex shedding (periodic) |
| `shear` | Free shear layer | Re=1,000 | 20 | Kelvin-Helmholtz instability |
| `t2d` | Taylor-Green 2D | Re=10,000 | 3 | Decaying turbulence |
| `os7000` | Orr-Sommerfeld/Poiseuille TS wave | Re=7,500 | 1 | Linear instability growth |

## Convection Approximation Methods

Set in `config.m` line 95:
```matlab
conv_approaches = {'fom', 'ftensor', 'rtensor', 'deim', 'clsdeim', 'mclsdeim'};
conv_approach = conv_approaches{6};  % Choose by index
```

Or via environment variable:
```bash
export NEKROM_CONV_APPROACH=deim
matlab -batch "driver"
```

### Method Comparison

| Method | Offline Cost | Online Cost | Accuracy | Stability | Use Case |
|--------|-------------|-------------|----------|-----------|----------|
| `fom` | Low | High O(N×M) | Exact | Excellent | Debugging, baseline |
| `ftensor` | Very High | Low O(N³) | Good | Excellent | Production (precomputed) |
| `rtensor` | Medium | Medium O(N³) | Good | Excellent | Moderate basis size |
| `deim` | Low | Very Low O(pN) | Fair | Poor | Fast prototyping |
| `clsdeim` | Low | Very Low O(pN) | Good | Good | Recommended default |
| `mclsdeim` | Medium | Very Low O(pN) | Best | Best | High-fidelity required |

**Legend:**
- N = ROM dimension (nb)
- M = FOM spatial resolution
- p = DEIM sample points (typically 2-4N)

### DEIM Variants Explained

**DEIM** (Discrete Empirical Interpolation Method):
- Sparse sampling of nonlinear term
- Fastest but least stable
- Can go unstable for long integrations

**CLS-DEIM** (Constrained Least-Squares DEIM):
- Adds linear constraint to enforce ROM subspace consistency
- Better stability than DEIM
- Minimal additional cost

**MCLS-DEIM** (Modified CLS-DEIM):
- Adds statistical regularization from training snapshots
- Best accuracy and stability
- Requires training data (automatically generated)

**Recommendation**: Start with `clsdeim`, upgrade to `mclsdeim` if accuracy critical.

## Configuration Parameters

Key settings in `config.m`:

### Time Integration
```matlab
nsteps = 40000;      % Total timesteps
dt = 1e-3;           % Timestep size
iostep = 100;        % Output frequency (every N steps)
```

### ROM Settings
```matlab
nb = 20;             % Number of POD modes (excluding mean)
nu = 1/1000;         % Kinematic viscosity
```

### Thermo-Fluid Settings
```matlab
kappa = 1.0;  % thermal diffusion coefficient (temperature Laplacian prefactor)
```
Environment overrides:
- `NEKROM_KAPPA` (thermal diffusion coefficient)
- `NEKROM_TDEIM_FROM_OPS=1` (use `ops/tdeim_*` when present)
- `MOR_DISABLE_TDEIM=1` (disable TDEIM even if `ops/tdeim_*` exists)
- `NEKROM_GX`, `NEKROM_GY`, `NEKROM_GZ` (gravity vector for buoyancy when `ops/buxt` etc exist)

### POD Inner Product
The MATLAB validation helpers read the offline `ops/ips` tag and support the same POD inner products as the Fortran path: `L2`, `H10`, and `HLM`.

`HLM` uses the Helmholtz metric from Fortran, approximated in MATLAB as `1/Re * H10 + (11/6)/dt * L2`. It therefore needs the Reynolds number and timestep when recomputing POD bases for validation.

The driver also reads stable case metadata from the case directory when it is available:
- `ops/nb`, `ops/ns`, `ops/ips`
- `.mor` fields `pod:type`, `pod:mode0`, `deim:mode`, and `deim:alpha`

If `<casename>.mor` is not present, the loader will use the single `.mor`
file in the case directory when that is unambiguous. This covers cases such
as `examples/rb_axi/rb.mor` without requiring the directory name and file
name to match.

`pod:mode0` controls the zeroth-mode reference used during validation (`avg` or `state`).

Runtime knobs such as `nsteps`, `dt`, `iostep`, and the visualization flags still come from `config.m` and environment-variable overrides.

### DEIM Settings (for deim/clsdeim/mclsdeim)
```matlab
ndeim_pts = 20;              % Number of DEIM sample points
os_multiplier = 2;           % Oversampling ratio
n_os_points = ceil(os_multiplier * ndeim_pts);

% Point selection algorithms
ps_algs = {'sopt', 'gpode', 'gappy_pod', 'gnat'};
ps_alg = ps_algs{1};         % 1=sopt, 2=gpode, 3=gappy_pod, 4=gnat

% Stability options
deim_finegrid = false;       % Interpolate POD to fine grid
deim_dealias = false;        % Use oversampled DEIM points
deim_dealias_cquad = false;  % EXPERIMENTAL: compressed quadrature on a 3/2 grid (not persisted unless opted in)
deim_dealias_quad = false;   % Use full 3/2-rule quadrature (MATLAB-only; not saved to ops/)
```

**Future research**: structure-preserving hyper-reduction that enforces a skew-adjoint convection operator in the physical $L^2$ energy inner product (so the online nonlinearity does near-zero kinetic-energy work), rather than relying on sampling stability heuristics alone.

**Note on `deim_dealias_cquad`**: this mode is experimental. By default, the driver keeps its operators in memory and does not overwrite `ops/deim_*`. To persist anyway, set `NEKROM_DEIM_DEALIAS_CQUAD_PERSIST=1`.

**Point selection algorithm recommendations**:
- `gpode` (QDEIM): Best default - fast, stable, adaptive
- `gappy_pod`: Simpler alternative, good for debugging
- `sopt`: Theoretically optimal but 100x slower (not recommended)
- `gnat`: Only for GNAT hyper-reduction method

### Output Flags
```matlab
ifvort = true;       % Compute vorticity field
ifwrite = true;      % Write field files
ifvis = false;       % Real-time visualization (slows execution)
ifplot = true;       % Save energy/momentum plots
if_run_tests = true; % Run validation tests on startup
```

## Environment Variable Override

All config parameters can be overridden via environment variables:

```bash
export NEKROM_CASE=cyl
export NEKROM_CONV_APPROACH=clsdeim
export NEKROM_NSTEPS=10000
export NEKROM_IOSTEP=200
export NEKROM_DEIM_FINEGRID=1
export NEKROM_DEIM_DEALIAS_CQUAD=1
export NEKROM_DEIM_CQUAD_MULT=3
export NEKROM_DEIM_DEALIAS_CQUAD_PERSIST=0
export NEKROM_IFVORT=0
export NEKROM_IFWRITE=1
export NEKROM_IFVIS=0
export NEKROM_IFPLOT=1
export NEKROM_IF_RUN_TESTS=0
export NEKROM_CONV_ENFORCE_SKEW_ADJOINT=1
export NEKROM_CONV_SKEW_INNER=l2
export NEKROM_CONV_SKEW_APPLY=deim

matlab -batch "driver"
```

**Useful for**:
- CI/CD pipelines
- Parameter sweeps (loop over env vars)
- Batch job submission

## Output Files

Results saved to timestamped directory: `output/{case}_{timestamp}/`

```
output/shear_2026-06-30-14-32-15/
├── fields/             # Nek5000-format field files
│   ├── shear0.f00001
│   ├── shear0.f00002
│   └── ...
├── ucoef               # ROM coefficients (ASCII, nb columns)
├── stability.mat       # Diagnostic metrics (MATLAB struct)
├── ke.pdf              # Kinetic energy plot
└── momentum.pdf        # Momentum conservation plot
```

### Reading Results

**MATLAB/Octave:**
```matlab
load('output/shear_2026-06-30-14-32-15/stability.mat');
results

% Access fields:
results.completed       % true if finished without NaN
results.nan_detected    % true if simulation went unstable
results.ke_rel_drift    % (KE_final - KE_initial) / KE_initial
results.max_ucoef_norm  % max ||u|| (detecting blow-up)
results.ucoef           % [num_outputs × (nb+1)] coefficient history
results.kes             % [num_outputs × 1] kinetic energy
results.momentums       % [num_outputs × 2] [u_momentum, v_momentum]
```

**Visualize:**
```matlab
plot(results.kes);
xlabel('Output step'); ylabel('Kinetic Energy');
```

## Testing

### Unit Tests
Test DEIM operator setup and basic functionality:
```matlab
run_deim_tests(snaps_path, casename, reorder, ndeim_pts, n_os_points, ps_alg)

% Or with defaults from config.m:
run_deim_tests
```

Validates:
- DEIM point count and uniqueness
- Matrix dimensions
- No NaN/Inf values
- Oversampling/quadrature paths

### Stability Comparison
Compare coarse-grid vs fine-grid DEIM stability:
```matlab
comparison = run_deim_stability_compare();

% Run for 10,000 steps with output every 100:
comparison = run_deim_stability_compare(10000, 100);
```

Returns struct with:
- `comparison.coarse` - results using standard grid
- `comparison.fine` - results with `deim_finegrid=true`
- Completion status, NaN detection, energy drift, wall time

### MATLAB vs Fortran Comparison
Compare the latest MATLAB run against the current Fortran case snapshots:
```matlab
comparison = run_driver_compare();

% Or limit to a subset of cases:
comparison = run_driver_compare({'ldc', 'cyl'});
```

The harness looks for:
- MATLAB output under `output/{case}_*/fields/`
- Fortran output under `../../examples/{case}/snaps/` or a direct `../../examples/{case}/` driver output directory

It reports snapshot counts plus relative errors for the field history, kinetic
energy, and momentum derived from the loaded snapshots.

## Troubleshooting

### "Operators not found" Error
**Cause**: Offline phase hasn't been run.

**Solution**:
```bash
cd ../../examples/shear  # Or your case directory
makerom shear            # Runs Nek5000 offline phase
```

This generates files in `ops/`:
- `au`, `bu`, `cu` - reduced operators
- `u0` - initial condition
- `uk` - constraints

### NaN Detected During Simulation
**Cause**: ROM went unstable.

**Possible fixes**:
1. Reduce timestep: `dt = 5e-4` (half current value)
2. Increase ROM modes: `nb = 30` (more accuracy)
3. Switch to more stable convection: `conv_approach = 'mclsdeim'`
4. Enable dealiasing: `deim_dealias_cquad = true` (recommended), `deim_dealias = true`, or `deim_dealias_quad = true`
5. Add stabilization: `ifleray = true` or `ifefr = true` in config

### "optim package not found" (Octave)
**Cause**: Missing optional Octave package.

**Solution**:
```bash
octave
> pkg install -forge optim
```

Only needed if `ifcopt = true` (constrained optimization for ROM solve).

### Slow Execution
**Causes**:
- Real-time visualization enabled: Set `ifvis = false`
- Writing too frequently: Increase `iostep` (e.g., 500 instead of 100)
- FOM convection: Switch to `clsdeim` or `ftensor`

**Expected runtimes** (nb=20, nsteps=10000, typical workstation):
- `fom`: ~30-60 seconds
- `ftensor`: ~5-10 seconds
- `clsdeim`/`mclsdeim`: ~3-5 seconds

### Memory Issues
**Cause**: Large FOM mesh or many snapshots.

**Solutions**:
1. Reduce ROM dimension: `nb = 10`
2. Use tensor instead of FOM: `conv_approach = 'ftensor'`
3. Increase swap space (OS-level)

## Advanced Usage

### Custom Initial Condition
Modify `driver.m` line 132:
```matlab
u(:,1) = u0;  % Default: load from ops/u0

% Custom IC examples:
u(:,1) = [1; randn(nb,1)*0.1];           % Small random perturbation
u(:,1) = [1; sin((1:nb)'*pi/nb)*0.5];    % Smooth modes
u(:,1) = load('my_restart.mat').ucoef;  % Restart from file
```

### ROM Stabilization Techniques
In `config.m`, enable one of:

```matlab
% Leray regularization (spectral viscosity)
ifleray = true;
radius = 0.01;  % Filter radius

% Evolve-then-filter (post-step smoothing)
ifefr = true;
relax = dt;     % Relaxation parameter

% Tikhonov regularization (energy penalty)
iftr = true;
relax = dt;

% Constrained optimization (box constraints on modes)
ifcopt = true;
```

**Recommendation**: Try `ifleray` first (most common in literature).

### Parameter Sweeps
Bash script example:
```bash
#!/bin/bash
for nb in 10 15 20 25 30; do
    for approach in clsdeim mclsdeim; do
        export NEKROM_CONV_APPROACH=$approach
        export NEKROM_NSTEPS=5000
        matlab -batch "nb=$nb; driver"
    done
done
```

### Parallel Execution
MATLAB Parallel Computing Toolbox (if available):

Edit `setup_conv_deim.m` line 29-36:
```matlab
% Change from:
for i = 1:nb

% To:
parfor i = 1:nb  % Parallel loop
```

Speedup: 2-4× for offline phase (one-time cost).

## File Organization

```
drive/matlab/
├── driver.m                 # Main integration loop
├── config.m                 # User-editable parameters
├── run_deim_tests.m         # Unit test suite
├── run_deim_stability_compare.m  # Stability benchmarks
├── operators/               # Convection operator implementations
│   ├── conv_fom.m           # Pseudo-ROM (full reconstruction)
│   ├── conv_tensor_*.m      # Tensor-based methods
│   ├── conv_deim.m          # DEIM family (deim/clsdeim/mclsdeim)
│   └── setup_conv_deim.m    # Offline DEIM setup
├── point_generators/        # DEIM point selection algorithms
│   ├── qdeim.m              # Q-DEIM (recommended)
│   ├── gappy_pod.m          # Gappy POD / classical DEIM
│   ├── s_opt.m              # S-optimality (slow, not recommended)
│   └── gnat.m               # GNAT point distribution
├── io/                      # File I/O and data structures
│   ├── load_full_ops.m      # Load offline operators
│   ├── output_fields.m      # Write solution fields
│   └── interp_basis_to_grid.m  # Grid interpolation
│   (Note: NekSnaps.m now provided by NekToolKit)
└── plotting/                # Post-processing scripts
```

## Dependencies

**Required** (built-in):
- Linear algebra: `chol`, `svd`, `qr`, `\` (backslash)
- I/O: `fopen`, `fread`, `fwrite`, `save`, `load`
- Plotting: `plot`, `patch`, `figure`

**Required** (external):
- **NekToolKit** (https://github.com/kent0/NekToolKit) **MANDATORY**
  - Functions used: `NekSnaps`, `zwgll`, `deriv_mat`, `deriv_geo`, `grad`, `interp_mat`
  - Purpose: Snapshot reading and spectral element utilities
  - Driver checks for these functions at startup and errors if missing

**Optional**:
- `pagemtimes` (MATLAB R2020b+): 2-3× faster tensor assembly
- `exportgraphics` (MATLAB R2020a+): Better PDF output (fallback: `print`)
- Parallel Computing Toolbox: Multi-core offline phase (marginal speedup)
- Octave `optim` package: Required only if `ifcopt=true` (constrained optimization)

## Performance Tips

1. **First run is slower**: MATLAB JIT compiler optimizes on second run (~30% speedup)
2. **Disable visualization**: `ifvis = false` for production runs (10× faster)
3. **Batch mode**: `matlab -batch` avoids GUI overhead
4. **Use tensor for production**: Precompute with `ftensor`, reuse for many runs
5. **Profile bottlenecks**:
   ```matlab
   profile on
   driver
   profile viewer
   ```

## Citing

If you use this driver for research, please cite:

```bibtex
@software{nekrom_matlab,
  title = {NekROM MATLAB Driver},
  author = {{NekROM Development Team}},
  url = {https://github.com/Nek5000/NekROM},
  year = {2024}
}
```

## Getting Help

1. Check this README for common issues
2. Review example cases in `../../examples/*/README.md`
3. Open issue: https://github.com/Nek5000/NekROM/issues
4. Read NekROM documentation: https://nekrom.readthedocs.io

## Contributing

See `../../doc/conventions.md` for coding style.

When adding a new convection operator:
1. Create `operators/conv_mymethod.m` with signature `out_coef = conv_mymethod(ucoef, ...)`
2. Add case to switch in `driver.m` line 160
3. Add to `conv_approaches` list in `config.m` line 94
4. Add tests in `run_tests.m`
5. Update this README with method description

## License

See `../../LICENSE.txt`
