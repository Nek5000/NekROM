# Breaking Changes - NekToolKit Now Mandatory

**Date:** 2026-06-30

## Summary

NekToolKit is now a **mandatory dependency** for NekROM MATLAB driver. The bundled `NekSnaps.m` has been removed in favor of the NekToolKit version.

## What Changed

### 1. NekToolKit is Now Required
Previously, NekToolKit was treated as optional with bundled fallbacks. Now it is **strictly required** and checked at driver startup.

### 2. Bundled NekSnaps.m Removed
- **Removed:** `drive/matlab/io/NekSnaps.m`
- **Reason:** NekToolKit provides a more feature-complete version
- **Archived as:** `NekSnaps.m.deprecated` (for reference only)

### 3. Mandatory Dependency Check
The driver now validates all NekToolKit functions at startup:
- `NekSnaps` - Snapshot reader
- `zwgll` - GLL quadrature nodes/weights
- `deriv_mat` - Spectral derivative matrix
- `deriv_geo` - Geometric derivatives
- `grad` - Gradient operator
- `interp_mat` - Interpolation matrix

**If any function is missing, the driver will fail immediately** with installation instructions.

## Migration Guide

### Before (Old Behavior)
```matlab
cd NekROM/drive/matlab
driver  % Might work without NekToolKit, but fail later
```

### After (New Behavior)
```matlab
# Install NekToolKit first
cd /path/to/workspace
git clone https://github.com/kent0/NekToolKit.git

# Verify installation
cd NekROM/drive/matlab
check_dependencies()  # Must pass before running driver

# Run driver
driver  # Will fail immediately if NekToolKit not available
```

## Error Messages

### Old Error (Confusing)
```
Undefined function 'zwgll' for input arguments of type 'double'.
Error in setup_conv_deim (line 24)
```

### New Error (Clear)
```
Error: NekToolKit dependency not satisfied.
Missing functions: NekSnaps, zwgll, deriv_mat

NekToolKit is required for NekROM MATLAB driver.
Install from: https://github.com/kent0/NekToolKit

Quick install:
  cd /path/to/workspace
  git clone https://github.com/kent0/NekToolKit.git

For detailed diagnostics, run: check_dependencies()
```

## Why This Change?

### Problems with Old Approach
1. **Hidden dependencies** - Users got cryptic errors deep in execution
2. **Code duplication** - Maintaining bundled copies of NekToolKit code
3. **Feature gaps** - Bundled NekSnaps.m lacked visualization features
4. **Version drift** - Risk of bundled code becoming stale

### Benefits of New Approach
1. **Clear requirements** - Users know upfront what's needed
2. **Better errors** - Immediate failure with installation guidance
3. **Single source** - NekToolKit maintained separately
4. **Full features** - Access to all NekToolKit capabilities

## Affected Workflows

### ✅ No Change Needed
If you already have NekToolKit installed (common setup), nothing changes.

### ⚠️ Action Required
If you were relying on bundled fallbacks:
1. Install NekToolKit from https://github.com/kent0/NekToolKit
2. Add to MATLAB path or place as sibling to NekROM
3. Verify with `check_dependencies()`

## Installation Options

### Option 1: Auto-Detection (Recommended)
```bash
cd /path/to/workspace
git clone https://github.com/kent0/NekToolKit.git

# Directory structure (driver auto-detects):
workspace/
├── NekROM/
└── NekToolKit/
```

### Option 2: Manual Path
```matlab
% Add to startup.m or before running driver
addpath('/path/to/NekToolKit/matlab')
```

### Option 3: Permanent Path
```matlab
% In MATLAB, add permanently
addpath('/path/to/NekToolKit/matlab')
savepath
```

## Verification

Run the dependency checker:
```matlab
cd NekROM/drive/matlab
check_dependencies()
```

Expected output:
```
=== NekROM MATLAB Dependency Check ===

Environment:
  ✓ Running on MATLAB R2023b

NekToolKit Functions (required):
  ✓ NekSnaps (from NekToolKit)
  ✓ zwgll (from NekToolKit)
  ✓ deriv_mat (from NekToolKit)
  ✓ deriv_geo (from NekToolKit)
  ✓ grad (from NekToolKit)
  ✓ interp_mat (from NekToolKit)

=== Summary ===
✓ All dependencies satisfied!
```

## Rollback (Not Recommended)

If you need the old behavior temporarily:
1. Restore bundled NekSnaps: `mv NekSnaps.m.deprecated NekSnaps.m`
2. Remove dependency check from `driver.m`

**Warning:** This is unsupported and may break in future updates.

## Questions?

- Check README.md for full documentation
- Run `check_dependencies()` for diagnostics
- Open issue on NekROM repository

## Related Changes

See also:
- `NEKTOOLKIT_DEPENDENCY_CHANGES.md` - Complete change log
- `io/NekSnaps_README.txt` - Notes on removed bundled version
- `check_dependencies.m` - New dependency validation tool
