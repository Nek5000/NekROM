## FOM: Flow past a cylinder in 2D.
A unit-diameter cylinder is centered at (0,0) in a box on [-2.5,17] x [-5,5].

With unit inflow velocity at x=-2.5, the Reynolds number is given by Re =
1/viscosity

Run (with NekROM/bin and Nek5000/bin in $PATH and $MOR_DIR set to NekROM) `./run_fom` to generate the FOM fields and drag/lift data.

## ROM: Flow past a cylinder in 2D

To obtain the ROM solution, run 

```
mkdir snaps && cp cyl0.f* snaps && ls snaps/*0.f* > file.list
makerom cyl_rom
nekbmpi cyl 1
```

```
grep dragx logfile | sed 's/dragx//g' > rom.dragx.dat
grep dragy logfile | sed 's/dragy//g' > rom.dragy.dat
```

The usr file for the ROM, cyl_rom.usr, mimics the FOM user file. In particular the rom_userchk subroutine is only a slight modification of the FOM userchk subroutine to accomodate variable name conflicts.

To see the drag vs. time:

grep drag logfile

The detail computation of drag in the ROM is in Kento Kaneko An Augmented Basis Method for Reduced Order Models of Turbulent Flow. Ph.D. thesis (2022).
