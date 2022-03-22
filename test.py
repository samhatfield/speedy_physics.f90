from physics import physics
import numpy as np
import xarray as xr

# Open NetCDF file containing an example of prognostic variables (taken from a
# run of speedy.f90)
progs = xr.open_dataset("prognostics.nc").squeeze().stack(ngp=("lat", "lon"))

# Infer grid dimensions from u-wind shape
nlev, ngp = progs["u"].shape

# Define sigma levels
σ_half = np.array([0.000, 0.050, 0.140, 0.260, 0.420, 0.600, 0.770, 0.900, 1.000])
σ = 0.5*(σ_half[1:] + σ_half[:-1])

# Extract prognostics from input data file, converting units to the ones used internally by Speedy,
# where appropriate
prog_u = np.transpose(progs["u"].data)
prog_v = np.transpose(progs["v"].data)
prog_t = np.transpose(progs["t"].data)
prog_q = np.transpose(progs["q"].data)*1000.0 # Convert from kg/kg to g/kg
prog_phi = np.transpose(progs["phi"].data)*9.81 # Convert to geopotential
prog_sp = np.transpose(progs["ps"].data/1.0e5) # Convert from Pa to 10^5 Pa

tend_u, tend_v, tend_t, tend_q = physics.get_physical_tendencies(prog_u, prog_v, prog_t, prog_q,
                                                                 prog_phi, prog_sp, σ)


