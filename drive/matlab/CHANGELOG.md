# Changelog - NekROM MATLAB Driver

All notable changes to the MATLAB/Octave driver are documented here.

## [Unreleased] - 2026-06-30

### Added
- **README.md** - Comprehensive user guide (682 lines)
  - Quick start instructions for MATLAB and Octave
  - Case descriptions and selection guide
  - Convection method comparison table with online costs
  - DEIM variants explained (DEIM/CLS-DEIM/MCLS-DEIM)
  - Complete parameter reference
  - Environment variable override documentation
  - Troubleshooting section with solutions
  - Advanced usage examples (custom IC, stabilization, parameter sweeps)
  - Performance tips and benchmarks
  - File organization reference
  
- **QUICKREF.md** - One-page quick reference card
  - Common commands and parameters
  - Speed vs accuracy tradeoff chart
  - Troubleshooting quick fixes
  - Performance expectations table
  
- **test_installation.m** - Installation verification script
  - Tests configuration loading
  - Verifies case directories exist
  - Checks required functions on path
  - Tests Octave compatibility
  - Validates point selection algorithms available
  - Summary report with actionable feedback
  
- **archive/README.md** - Documentation for archived deprecated code
  - Explains what's in archive and why
  - Guidance on whether to use archived code (no)
  - Historical context for code evolution

### Modified
- **config.m**
  - Added comprehensive parameter documentation at file header
  - Added early validation for `conv_approach` with helpful error messages
  - Added validation for `ps_alg` when using DEIM methods
  - Validates parameters immediately after loading from environment
  
- **driver.m**
  - Added file existence checks with actionable error messages:
    - Case directory validation
    - Snapshots directory validation  
    - POD basis file validation
    - Operators directory and file validation
  - Improved progress reporting:
    - Replaced verbose per-iostep print with 5% increments
    - Shows percentage, step numbers, and current time
    - Cleaner output for long runs
  - Added NekToolKit missing warning (informative, not error)
    - Explains NekToolKit is optional
    - Provides installation link
    - Only shows in MATLAB (silent in Octave)

### Changed
- **Code organization**
  - Moved `old/` → `archive/old/`
  - Moved `io/old/` → `archive/io_old/`
  - Moved `operators/old/` → `archive/operators_old/`
  - Moved `point_generators/old/` → `archive/point_generators_old/`
  - Clarifies which code is current vs deprecated
  - Reduces clutter in main directories
  - Preserves history for reference

### Fixed
- Cryptic "file not found" errors now show:
  - Exactly what file is missing
  - Where it should be located
  - How to generate it (run makerom command)
- Parameter typos caught early instead of failing mid-simulation
- Progress output less verbose (was printing every iostep)

### Deprecated
- Code in `archive/` directories is no longer maintained
- See `archive/README.md` for details on what was superseded

---

## Previous Versions

### [Unversioned] - Pre-2026-06-30

Original driver implementation. Notable features:
- Clean modular architecture (`operators/`, `io/`, `point_generators/`)
- Multiple convection methods (FOM, tensor, DEIM variants)
- Octave compatibility layer
- Comprehensive testing infrastructure
- Multiple example cases (ldc, cyl, shear, t2d)

**Issues addressed in 2026-06-30 update**:
- ❌ No user documentation
- ❌ Cryptic error messages
- ❌ Late parameter validation
- ❌ Verbose output
- ❌ Deprecated code mixed with current
- ❌ No installation test

---

## Migration Guide

### From Pre-2026-06-30 to 2026-06-30

**No breaking changes.** All existing workflows continue to work.

**Optional improvements you can adopt:**

1. **Better error messages** - automatic, no action needed
2. **Cleaner progress output** - automatic, no action needed
3. **Installation test** - run `test_installation` to verify setup
4. **Documentation** - read `README.md` for tips and best practices
5. **Quick reference** - bookmark `QUICKREF.md` for common tasks

**If you used `old/` code:**
- Check `archive/` for moved files
- Consider migrating to current implementations
- See `archive/README.md` for guidance

**If you have automation scripts:**
- Environment variables still work the same way
- Return values/output format unchanged
- Progress messages changed but exit codes same

---

## Semantic Versioning

This project uses informal versioning. Dates mark significant updates.

**Version scheme**: `YYYY-MM-DD` for major updates, unreleased for minor

**Breaking change policy**: Avoid whenever possible, document clearly when necessary

**Backward compatibility**: Maintained unless security/correctness requires break

---

## Contributing

When making changes:
1. Update this CHANGELOG with your changes
2. Document new parameters in `config.m` header and `README.md`
3. Add tests to `test_installation.m` if adding dependencies
4. Update `QUICKREF.md` if adding common operations
5. Follow existing code style (see `../../doc/conventions.md`)

---

## Acknowledgments

**2026-06-30 improvements**: Documentation, error handling, code organization
- Identified issues through code review and user feedback simulation
- Addressed all critical usability concerns
- Maintained backward compatibility and Octave support

**Original implementation**: Modular architecture, multiple methods, testing
- Clean separation of concerns
- Octave compatibility from start
- Comprehensive operator implementations

---

## See Also

- **README.md** - Complete user documentation
- **QUICKREF.md** - Quick reference for common tasks
- **archive/README.md** - Information about deprecated code
- **test_installation.m** - Verify your setup
- **../../doc/** - NekROM project-wide documentation

---

*This CHANGELOG follows [Keep a Changelog](https://keepachangelog.com/) principles.*
