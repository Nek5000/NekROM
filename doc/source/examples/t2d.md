# Taylor-Green Vortex / Kovasznay

`examples/t2d` is a 2D Taylor-Green vortex / Kovasznay convergence test rather than a full ROM walkthrough.

## Purpose

The case demonstrates spectral convergence for the PN-PN discretization and can be switched to PN-PN-2 by editing `SIZE`.

## Running

The case README tells users to run `run_all`, which executes a small collection of cases and writes the resulting errors to
`err.tot`.

## Post-Processing

The output can be plotted with the provided `gnu.in` script in gnuplot.
