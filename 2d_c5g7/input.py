import numpy as np
import mcdc
from c5g7_2d_model import model, ZERO
import shutil

mcdc = model()
ZERO = 1e-3

# =============================================================================
# iQMC Parameters
# =============================================================================
pitch = 1.26
N = 15_000
Nx = 17 * 3 * 2
Ny = 17 * 3 * 2
G = 7
x_grid = np.linspace(ZERO, pitch * 17 * 3, Nx + 1)
y_grid = np.linspace(-pitch * 17 * 3, -ZERO, Ny + 1)

phi0 = np.ones((G, Nx, Ny))

mcdc.iQMC(
    x=x_grid,
    y=y_grid,
    maxit=150,
    g=np.ones(G),
    phi0=phi0,
    mode='batched',
    sample_method='halton',
    scores=['source-x', 'source-y'],
)

# =============================================================================
# run mcdc
# =============================================================================

# Setting
mcdc.setting(N_particle=N, rng_seed=123456)
mcdc.eigenmode(N_inactive=100, N_active=50)
# Run
mcdc.run()
# shutil.move('output.h5', f'1_results/constant/N={N}.h5')


