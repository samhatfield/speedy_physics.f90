from physics import physics
import numpy as np

nlat, nlon, nlev = 96, 48, 8
ngp = nlat*nlon

σ_half = np.array([0.000, 0.050, 0.140, 0.260, 0.420, 0.600, 0.770, 0.900, 1.000])
σ = 0.5*(σ_half[1:] + σ_half[:-1])

prog_u = np.random.rand(ngp,nlev)
prog_v = np.random.rand(ngp,nlev)
prog_t = np.random.rand(ngp,nlev)
prog_q = np.random.rand(ngp,nlev)
prog_phi = np.random.rand(ngp,nlev)
prog_sp = np.random.rand(ngp)

tend_u, tend_v, tend_t, tend_q = physics.get_physical_tendencies(prog_u, prog_v, prog_t, prog_q, prog_phi, prog_sp, σ)

assert np.isclose(10.0*tend_u, prog_u).all()
assert np.isclose(10.0*tend_v, prog_v).all()
assert np.isclose(10.0*tend_t, prog_t).all()
assert np.isclose(10.0*tend_q, prog_q).all()

print("All tests passed")

