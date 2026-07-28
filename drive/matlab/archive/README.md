# Archive Directory

This directory contains deprecated or superseded code kept for reference.

**These files are NOT maintained and may not work with current codebase.**

## Contents

### `old/` - Deprecated main driver code
- `rom_online_solver.m` - Original monolithic driver (superseded by `driver.m`)
- `driver.m` - Early version of modular driver
- Various test scripts

**Replaced by**: Current `driver.m` with modular structure

### `io/old/` - Deprecated I/O functions
- `write_field.m` - Old field writer with TODOs for 3D support
- `output_fields.m` - Early version of output routine
- `get_grid.m` - Superseded version

**Replaced by**: `io/` module with improved implementations

### `operators/old/` - Deprecated operator implementations
- `conv_deim.m` - Early DEIM with mixed concerns (should separate DEIM/CLS-DEIM)
- `conv_tensor_*.m` - Various tensor assembly approaches tried and abandoned
- `gen_Au.m` - Early linear operator assembly

**Replaced by**: `operators/` module with clean separation (conv_deim.m, conv_tensor_dealiased.m, etc.)

### `point_generators/old/` - Deprecated point selection algorithms
- `s_opt/` - Multiple iterations of S-optimality (original, improved, performance)
- `gpode/` - Original GPODE implementation before the current refactoring
- `gappy_pod/` - Multiple Gappy POD variants (original, new, QR, performance)
- `gnat/` - GNAT performance tests

**Replaced by**: `point_generators/` with cleaned-up implementations (qdeim.m, gappy_pod.m, s_opt.m, gnat.m)

## Why Keep This?

1. **Historical reference**: Shows evolution of algorithms and design decisions
2. **Performance comparisons**: Benchmark scripts show speedup of current versions
3. **Alternative approaches**: Some abandoned methods may inspire future work
4. **Code archaeology**: Helps understand why certain choices were made

## Should I Use Code From Here?

**No.** Use the current implementations in the parent directory:
- `../driver.m` - Current driver
- `../io/*.m` - Current I/O functions
- `../operators/*.m` - Current operator implementations
- `../point_generators/*.m` - Current point selection algorithms

## Cleaning Up

If you're sure you don't need this history, you can delete the entire `archive/` directory:
```bash
cd drive/matlab
rm -rf archive/
```

## Last Updated

Archive created: 2026-06-30

Original "old/" directories moved here during codebase cleanup to:
- Clarify what's maintained vs deprecated
- Remove clutter from main path
- Preserve history without confusion
