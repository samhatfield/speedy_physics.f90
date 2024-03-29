from src.physics import physics
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Open NetCDF file containing an example of prognostic variables (taken from a
# run of speedy.f90)
progs_ds = xr.open_dataset("prognostics.nc")

# Extract Gaussian latitudes
lats = progs_ds.coords["lat"].values

# Squeeze lat-lon into a single dimension "ngp"
progs = progs_ds.squeeze().stack(ngp=("lat", "lon"))

# Infer grid dimensions from u-wind shape
nlev, ngp = progs["u"].shape

# Define sigma levels
σ_half = np.array([0.000, 0.050, 0.140, 0.260, 0.420, 0.600, 0.770, 0.900, 1.000])

# Open land-sea mask
lsm = xr.open_dataset("lsm.nc").squeeze().stack(ngp=("lat", "lon"))["lsm"].data

# Extract prognostics from input data file, converting units to the ones used internally by Speedy,
# where appropriate
prog_u = np.transpose(progs["u"].data)
prog_v = np.transpose(progs["v"].data)
prog_t = np.transpose(progs["t"].data)
prog_q = np.transpose(progs["q"].data)*1000.0 # Convert from kg/kg to g/kg
prog_phi = np.transpose(progs["phi"].data)*9.81 # Convert to geopotential
prog_sp = np.transpose(progs["ps"].data/1.0e5) # Convert from Pa to 10^5 Pa

# The time in the year as a fraction (0.0 - 1.0)
# Needed to compute the solar insolation
tyear = 0.0

tend_u, tend_v, tend_t, tend_q = physics.get_physical_tendencies(prog_u, prog_v, prog_t, prog_q,
                                                                 prog_phi, prog_sp,
                                                                 σ_half, lats, tyear, lsm)

# Model levels
levels = range(1, nlev+1)

# Only plot 10 grid points
npoints = 10
np.random.seed(1)
grid_points_to_plot = np.random.choice(tend_t.shape[0], npoints)

fig, ax = plt.subplots(ncols=2, figsize=(12, 6), sharey=True)
plt.suptitle("Parametrised physical tendencies for 10 selected grid points")

for i in grid_points_to_plot:
    ax[0].plot(tend_t[i,:], levels)
ax[0].invert_yaxis()
ax[0].set_title("Temperature tendency")
ax[0].set_ylabel("Model level")

for i in grid_points_to_plot:
    ax[1].plot(tend_q[i,:], levels)
ax[1].set_title("Humidity tendency")

plt.savefig("tendencies.png", bbox_inches="tight")
plt.show()

