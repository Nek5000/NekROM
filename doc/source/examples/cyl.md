# Cylinder Example

`examples/cyl` is the shortest end-to-end workflow in NekROM. It generates snapshots for 2D flow past a cylinder, builds the
ROM, and post-processes drag and lift data.

## FOM

The cylinder case runs a 2D incompressible flow problem with a unit-diameter cylinder in a rectangular domain. The `run_fom`
script generates the mesh, runs the Nek5000 FOM, and copies the output snapshots into `snaps/`.

## ROM

The `run_rom` script builds the POD bases from the snapshot set and then runs the ROM in Nek5000. It also extracts the drag
history into `rom.dragx.dat` and `rom.dragy.dat`.

## Reference

The original case notes are in `examples/cyl/README.md`, and the same workflow is summarized in the main quickstart guide.
