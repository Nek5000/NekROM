# NekROM MATLAB Driver - Quick Reference

One-page cheat sheet for common operations.

## Run Simulation

```bash
cd drive/matlab
matlab -batch "driver"              # MATLAB
octave --eval "driver"              # Octave
```

## Change Case

**Method 1: Edit config.m**
```matlab
thiscase = cases{2};  % 1=ldc, 2=cyl, 3=shear, 4=t2d, 5=os7000
```

**Method 2: Environment variable**
```bash
export NEKROM_CASE=cyl
matlab -batch "driver"
```

## Change Convection Method

**Method 1: Edit config.m**
```matlab
conv_approach = conv_approaches{4};  % 1=fom, 2=ftensor, 3=rtensor,
                                     % 4=deim, 5=clsdeim, 6=mclsdeim
```

**Method 2: Environment variable**
```bash
export NEKROM_CONV_APPROACH=clsdeim
matlab -batch "driver"
```

## Common Parameters

| Parameter | config.m | Environment Variable | Default |
|-----------|----------|---------------------|---------|
| Case | `thiscase = cases{3}` | `NEKROM_CASE=shear` | shear |
| Timesteps | `nsteps = 10000` | `NEKROM_NSTEPS=10000` | 40000 |
| Output freq | `iostep = 200` | `NEKROM_IOSTEP=200` | 100 |
| ROM modes | `nb = 15` | - | 20 |
| Convection | See above | `NEKROM_CONV_APPROACH` | mclsdeim |
| Visualization | `ifvis = true` | `NEKROM_IFVIS=1` | false |

## Convection Methods (Speed vs Accuracy)

```
Fastest        →        Most Accurate
deim < clsdeim < mclsdeim < rtensor < ftensor < fom
│                        │                       │
Unstable              Recommended            Exact
```

**Recommendation**: Start with `clsdeim`, upgrade to `mclsdeim` if needed.

## Troubleshooting

### "Case directory not found"
```bash
# Check case name spelling
export NEKROM_CASE=cyl  # not 'cylinder'
```

### "Operators not found"
```bash
cd ../../examples/shear
makerom shear           # Run offline phase
cd ../../drive/matlab
matlab -batch "driver"
```

### NaN detected during run
```matlab
% Option 1: Reduce timestep
dt = 5e-4;  % Half the current value

% Option 2: Switch to more stable method
conv_approach = conv_approaches{6};  % mclsdeim

% Option 3: Enable dealiasing
deim_dealias = true;
```

### Too slow
```matlab
ifvis = false;          % Disable real-time plotting
iostep = 500;           % Write less frequently
conv_approach = 'clsdeim';  % Faster than fom/ftensor
```

## File Locations

```
Results:     output/{case}_{timestamp}/
             ├── fields/               (field files)
             ├── ucoef                 (ROM coefficients)
             ├── stability.mat         (diagnostics)
             ├── ke.pdf                (energy plot)
             └── momentum.pdf          (momentum plot)

Operators:   ../../examples/{case}/ops/
             ├── au, bu, cu            (reduced operators)
             ├── ips                  (offline POD inner product)
             ├── u0                    (initial condition)
             └── uk                    (constraints)

Snapshots:   ../../examples/{case}/snaps/
             ├── bas{case}0.f*         (POD basis)
             ├── cba{case}0.f*         (nonlinear basis)
             └── csn{case}0.f*         (training snapshots)
```

## Quick Checks

```matlab
% Test installation
test_installation

% Verify operators exist
ls('../../examples/shear/ops/')

% Check ROM dimension
config; fprintf('ROM dimension: %d\n', nb);

% View last results
load('output/latest_run/stability.mat');
plot(results.kes); title('Kinetic Energy');
```

## Batch Runs (Parameter Sweep)

```bash
#!/bin/bash
for nb in 10 15 20 25; do
    for method in clsdeim mclsdeim; do
        export NEKROM_CONV_APPROACH=$method
        matlab -batch "nb=$nb; driver"
    done
done
```

## One-Liners

```bash
# Quick test (100 steps)
matlab -batch "config; nsteps=100; iostep=10; driver"

# No visualization, fast
matlab -batch "config; ifvis=false; ifwrite=false; driver"

# Change case and method
export NEKROM_CASE=cyl NEKROM_CONV_APPROACH=clsdeim
matlab -batch "driver"

# Run tests only
matlab -batch "run_deim_tests"
```

## Performance Expectations

| Case | Modes | Method | Steps | Time |
|------|-------|--------|-------|------|
| shear | 20 | clsdeim | 10k | ~5s |
| cyl | 20 | mclsdeim | 10k | ~8s |
| ldc | 30 | ftensor | 10k | ~12s |
| t2d | 3 | fom | 10k | ~3s |

*Typical workstation, MATLAB R2022a*

## Getting Help

1. **README.md** - Full documentation
2. **test_installation** - Verify setup
3. **Issue tracker** - https://github.com/Nek5000/NekROM/issues

## Pro Tips

- First run is slower (JIT compilation)
- Use `fom` for debugging (exact, easy to verify)
- Use `clsdeim` for production (fast + stable)
- Keep `ifvis=false` except when debugging
- Check `results.nan_detected` after runs
- Monitor `results.ke_rel_drift` for energy conservation

---

*For complete documentation, see README.md*
