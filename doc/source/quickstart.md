(quickstart-section-tag)=

# Quickstart

## Running the flow past a cylinder example

In its most straightforward workflow, using NekROM requires running a full-order model simulation to generate solution snapshots
and then executing NekROM targeting these snapshots to generate POD bases and run the reduced-order simulation. Here we illustrate
the process with the 2D flow past a cylinder example.

1. Clone NekROM and Nek5000

```shell
git clone https://github.com/Nek5000/NekROM.git
git clone https://github.com/Nek5000/Nek5000.git
```

2. Set environment variables

```shell
export PATH=$PATH:$(pwd)/NekROM/bin:$(pwd)/Nek5000/bin
export MOR_DIR=$(pwd)/NekROM/bin
```

3. Build Nek5000 tools

```shell
cd Nek5000/tools
./maketools
cd ../../
```

4. Change to cylinder example directory

```shell
cd NekROM/examples/cyl
```

5. Generate the FOM snapshots

```shell
./run_fom
```

This run script executes the Nek5000 cylinder case and extracts drag and lift data.
See below for more details.

6. Run the ROM using the FOM snapshots

```shell
./run_rom
```

This run script creates the POD bases and uses them to run the reduced-order simulation.
See below for more details.

7. Visualize the ROM

```shell
visnek romcyl
```

Open the resulting `romcyl.nek5000` file with Visit or Paraview to visualize the ROM output

### Understanding the run scripts

The `run_fom` and `run_rom` scripts automate several tasks in the NekROM workflow. NekROM users
are encouraged to write similar run scripts for their own cases.

The `run_fom` script,
shown below, builds the Nek5000 case for the FOM, gnerates the mesh, and runs the Nek5000 simulation.
After the simulation completes, the script copies the output snapshots to the `snaps` directory
and creates `file.list` to tell NekROM which snapshots to use. Finally, the script
extracts some drag and lift data from the log file.

```{literalinclude} ../../examples/cyl/run_fom
:language: shell
```

The `run_rom` script builds the `cyl_rom` case, executes the ROM with Nek5000, and extracts the
data lift and drag data from the log file of the ROM simulation.

```{literalinclude} ../../examples/cyl/run_rom
:language: shell
```

### NekROM output files

Executing NekROM generates several different kinds of files.

`romcyl0.*`: Snapshots of the ROM solution over time, similar to the snapshots of the FOM.

`bascyl0.*`: Snapshots of the POD modes. The mode number is indicated in the file extension.

`avgcyl0.*`: The average modes. These files include the average of the snapshots and the reconstructed average.

`uiccyl0.*`: Initial conditions

`tkecyl0.*`: Turbulent kinetic energy

`lapcyl0.*`:

`tmncyl0.*`:

### NekROM input files

`LMOR`: Compile-time parameters for NekROM. This is generated when `makerom` is executed.

```fortran
! MOR Compile-Time Allocation Parameters

parameter (lmu=1)   ! 0 -> disable velocity allocation
parameter (lmp=1)   ! 0 -> disable pressure allocation
parameter (lmt=1)   ! 0 -> disable temperature allocation

parameter (ls=500)  ! max number of snapshots
parameter (lcs=ls)  ! max number of coefficient set

parameter (lei=0)   ! 0 -> one residual, 1 -> affine decomp
parameter (lb=20)       ! max number of basis
parameter (lelm=lelt)   ! number of local elements
parameter (ltr=1+0*399) ! max number of tensor rank
parameter (ledvis=0)   ! 0 -> disable eddy viscosity allocation

parameter (lk=1)   ! largest wave number for pdrag calculation
parameter (lmsk=1) ! number of partitions

parameter (lintp=1) ! max number of interpolation points
parameter (lbat=1024) ! max size of batch vector (inner-product iter.)

parameter (lbavg=1+0*(lb-1)) ! size of average field allocation

! Auxiliary

parameter (lsu=(lcs-1)*lmu+1) ! size of velocity snapshots allocation
parameter (lsp=(lcs-1)*lmp+1) ! size of pressure snapshots allocation
parameter (lst=(lcs-1)*lmt+1) ! size of temperature snapshots allocation

parameter (lub=(lb-1)*lmu+1) ! size of velocity basis allocation
parameter (lpb=(lb-1)*lmp+1) ! size of pressure basis allocation
parameter (ltb=(lb-1)*lmt+1) ! size of temperature basis allocation

parameter (lres=1) ! size of residual storage
parameter (lres_u=((3*lb+lb**2)-1)*lei+1) ! size of residual storage for vel
parameter (lres_t=((2*lb+lb**2)-1)*lei+1) ! size of residual storage for temp
```

`cyl.mor`: Run-time parameters for NekROM

```{literalinclude} ../../examples/cyl/cyl.mor
:language: text
```

`cyl_rom.usr`: User specified functions for NekROM case. This is similar to the Nek5000 `cyl_fom.usr` file, but
also has several NekROM specific functions. Additionally, `param(170) = -1` is added to `userchk` in this
file to tell NekROM to read from `cyl.mor` rather than the FOM `cyl.rea` file. NekROM specific functions
include the following: `rom_userchk`, `rom_userbases`, `rom_userfop`, and `rom_userrhs`.

### Running parametrically

TODO
