# Speedy Physics
This is a directory that exclusively contains the SPEEDY fortran routines version 41 and related files used in 
the physical parameterizations. The routine calling all parameterization routines is 'phy_phypar_grid.f'

For a documentation of the SPEEDY physical parameterization scheme, please consult the web-page:
http://users.ictp.it/~kucharsk/speedy-net.html 
If you use the physical parameterization schemes in any publication, please cite the 3 papers listed there. 

The licence for these routines is in the text file licence.txt. 

## Requirements

- Linux/macOS/WSL with GNU Fortran
- Anaconda/Miniconda

## Set-up

```sh
conda env create -f environment.yml
conda activate speedy_physics
```

## How to run the test

```sh
pushd src
make
popd
python test.py
```
