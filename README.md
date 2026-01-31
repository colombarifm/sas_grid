# sas_grid

sas_grid generates grids based on the solvent accessible surface of a given structure. 
It reads a .xyz file and writes out a new file containing the grid point coordinates (labeled as XX). It can generate
a SAS grid around the molecular structure or considering the free space in a cavity or other regions of interest.

sas_grid is a free software written in Fortran 2003 language, being available at https://github.com/colombarifm/sas_grid under the GPLv3+ License. 
It runs under Linux environment with gfortran/gcc 5.4+ compilers.

## Build with make
Download the github repository and build it with make

```bash
   git clone https://github.com/colombarifm/themis.git
   cd sas_grid
   make
```


## usage

For SAS grid generation, one should use the following command:

`sas_grid --input [FILE] --radius [RADIUS] --factor [FACTOR] --type sas`

`[FILE]` is a .xyz coordinate file.

`[RADIUS]` is the solvent probe radius, in Angstrom.

`[FACTOR]` is an integer factor for the tessellation sphere. For each atom-centered sphere, N_points = 2 + factor^2 * 10.

For a cavity grid generation, one should use the following command:

`sas_grid --input [FILE] --radius [RADIUS] --type box --min [Xmin] [Ymin] [Zmin] --max [Xmax] [Ymax] [Zmax]`

`[Xmin] [Ymin] [Zmin]` are the minimum values for xyz coordinates

`[Xmax] [Ymax] [Zmax]` are the maximum values for xyz coordinates

## TODO

vdW radii for transition elements are still missing. The code will be updated soon.
