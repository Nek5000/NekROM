NekSnaps.m has been moved to NekToolKit (required dependency).

The bundled NekSnaps.m has been removed as of 2026-06-30.
NekROM now requires NekToolKit to be installed.

Install NekToolKit:
  https://github.com/kent0/NekToolKit

The NekToolKit version provides:
- All functionality of the bundled version
- Additional visualization features (show, next, prev)
- Better maintained and feature-complete

If you encounter "Undefined function 'NekSnaps'":
  1. Ensure NekToolKit is installed
  2. Add NekToolKit/matlab to your MATLAB path
  3. Or place NekToolKit as sibling to NekROM (auto-detected)
  4. Run check_dependencies() for diagnostics

For reference, the deprecated bundled version is available as:
  NekSnaps.m.deprecated
