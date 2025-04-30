.. _quickstart_section_tag:
  
Quickstart
==========

----------------------------------------
Running the flow past a cylinder example
----------------------------------------

In its most straightforward workflow, using NekROM requires running a full-order model simulation to generate solution snapshots
and then executing NekROM targeting these snapshots to generate POD bases and run the reduced-order simulation. Here we illustrate
the process with the 2D flow past a cylinder example.

1. Clone NekROM and Nek5000

.. code-block:: shell

    git clone https://github.com/Nek5000/NekROM.git
    git clone https://github.com/Nek5000/Nek5000.git

2. Set environment variables

.. code-block:: shell
   
    export PATH=$PATH:$(pwd)/NekROM/bin:$(pwd)/Nek5000/bin
    export MOR_DIR=$(pwd)/NekROM/bin

3. Build Nek5000 tools

.. code-block:: shell

    cd Nek5000/tools
    ./maketools
    cd ../../

4. Change to cylinder example directory

.. code-block:: shell

    cd NekROM/examples/cyl

5. Generate the FOM snapshots

.. code-block:: shell

    ./run_fom

This run script executes the Nek5000 cylinder case and extracts drag and lift data.
See below for more details.

6. Run the ROM using the FOM snapshots

.. code-block:: shell

    ./run_rom

This run script creates the POD bases and uses them to run the reduced-order simulation.
See below for more details.


7. Visualize the ROM

.. code-block:: shell

    visnek romcyl

Open the resulting `romcyl.nek5000` file with Visit or Paraview to visualize the ROM output

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Understanding the run scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `run_fom` and `run_rom` scripts automate several tasks in the NekROM workflow. NekROM users
are encouraged to write similar run scripts for their own cases. 

The `run_fom` script,
shown below, builds the Nek5000 case for the FOM, gnerates the mesh, and runs the Nek5000 simulation.
After the simulation completes, the script copies the output snapshots to the `snaps` directory
and creates `file.list` to tell NekROM which snapshots to use. Finally, the script
extracts some drag and lift data from the log file.

.. literalinclude:: ../../examples/cyl/run_fom
    :language: shell

The `run_rom` script builds the `cyl_rom` case, executes the ROM with Nek5000, and extracts the
data lift and drag data from the log file of the ROM simulation.

.. literalinclude:: ../../examples/cyl/run_rom
    :language: shell

^^^^^^^^^^^^^^^^^^^
NekROM output files
^^^^^^^^^^^^^^^^^^^
Executing NekROM generates several different kinds of files.

`romcyl0.*`: Snapshots of the ROM solution over time, similar to the snapshots of the FOM.

`bascyl0.*`: Snapshots of the POD modes. The mode number is indicated in the file extension.

`avgcyl0.*`: The average mode I assume. Why are there multiple ones of these.

`lapcyl0.*`: ????

`tkecyl0.*`: Turbulent kinetic energy

`tmncyl0.*`: ????

`uiccyl0.*`: ????

^^^^^^^^^^^^^^^^^^
NekROM input files
^^^^^^^^^^^^^^^^^^

`LMOR`: Compile-time parameters for NekROM

`cyl.mor`: Run-time parameters for NekROM

`cyl_rom.usr`: Similar to `.usr` files in Nek5000 with several NekROM specific functions.

^^^^^^^^^^^^^^^^^^^^^^
Running parametrically
^^^^^^^^^^^^^^^^^^^^^^
TODO
