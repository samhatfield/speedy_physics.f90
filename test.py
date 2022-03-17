from params import params
from physics import physics
import numpy as np

prog_u = np.random.rand(params.ngp,params.nlev)
prog_v = np.random.rand(params.ngp,params.nlev)
prog_t = np.random.rand(params.ngp,params.nlev)
prog_q = np.random.rand(params.ngp,params.nlev)
prog_phi = np.random.rand(params.ngp,params.nlev)
prog_sp = np.random.rand(params.ngp)

tend_u, tend_v, tend_t, tend_q = physics.get_physical_tendencies(prog_u, prog_v, prog_t, prog_q, prog_phi, prog_sp)

assert np.isclose(10.0*tend_u, prog_u).all()
assert np.isclose(10.0*tend_v, prog_v).all()
assert np.isclose(10.0*tend_t, prog_t).all()
assert np.isclose(10.0*tend_q, prog_q).all()

print("All tests passed")

