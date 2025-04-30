.. _quickstart_section_tag:
  
Quickstart
==========

## Running the flow past a cylinder example

In its most straightforward workflow, using NekROM requires running a full-order model simulation to generate solution snapshots
and then executing NekROM targeting these snapshots to generate POD bases and run the reduced-order simulation. Here we illustrate
the process with the 2D flow past a cylinder example.

### 1. Clone NekROM and Nek5000
```sh
git clone https://github.com/Nek5000/NekROM.git
git clone https://github.com/Nek5000/Nek5000.git
```
Maybe build Nek5000 tools?

### 2. Set environment variables
```sh
export PATH=$PATH:NekROM/bin:Nek5000/bin
export MOR_DIR=NekROM/bin
```

### 3. Change to cylinder example directory
```sh
cd NekROM/examples/cyl
```

### 4. Generate the FOM snapshots
```sh
./run_fom
```

Also generate drag and lift data
Explain the run script here. Maybe can include the script text as well?

### 5. Run the ROM using the FOM snapshots
```sh
./run_rom
```
Creates the POD bases and uses them to run the simulation as a reduced-order model

Explain the run script here. Maybe can include the script text as well.

### 6. Visualize the ROM
```sh
./visnek romcyl
```
Open the resulting `romcyl.nek5000` file with Visit or Paraview to visualize the ROM output

The `run_fom` and `run_rom` scripts automate several tasks in the NekROM workflow. We explain these scripts in more detail below. TODO

## Running the lid-driven cavity example 

(Might be the same process, worth adding?)
ADD ME

## What are all of these files NekROM produces? (Rename)
ADD ME

## How do I run parametric ROM?
