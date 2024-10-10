#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 09:42:43 2023
@author: sampasmann

The Takeda 1 Benchmark is a 2G 3D PWR quarter core problem. See:
    https://www.tandfonline.com/doi/pdf/10.1080/18811248.1991.9731408
"""
import numpy as np
import mcdc
import argparse
import h5py
# Two cases: (1) control rod is withdrawn, (2) control inserted
case = 1

# =============================================================================
# Calculate XS
# =============================================================================

chi = np.array([[1.0, 1.0], [0.0, 0.0]])


#### core ####
core_sigmaT= np.array([2.23775e-1, 1.03864])
core_absorption = np.array([8.52709e-3, 1.58196e-1])
core_production = np.array([9.09319e-3, 2.90183e-1])
core_sigmaS = np.array([[1.92423e-1, 0.0], 
                        [2.28253e-2, 8.80439e-1]])

core_nu = np.array([2.5, 2.5])
core_sigmaF = core_production / core_nu
core_sigmaC = (core_sigmaT - core_sigmaF - core_sigmaS.sum(axis=0))


#### control rod ####
rod_sigmaT= np.array([8.52325e-2, 2.1746e-1])
rod_absorption = np.array([1.74439e-2, 1.82224e-1])
rod_production = np.array([0.0, 0.0])
rod_sigmaS = np.array([[6.77241e-2, 0.0], 
                       [6.45461e-5, 3.52358e-2]])

rod_nu = np.array([2.5, 2.5])
rod_sigmaF = rod_production / rod_nu
rod_sigmaC = (rod_sigmaT - rod_sigmaF - rod_sigmaS.sum(axis=0))


#### reflector ####
reflector_sigmaT= np.array([2.50367e-1, 1.64482])
reflector_absorption = np.array([4.16392e-4, 2.02999e-2])
reflector_production = np.array([0.0, 0.0])
reflector_sigmaS = np.array([[1.93446e-1, 0.0], 
                             [5.65042e-2, 1.62452]])

reflector_nu = np.array([2.5, 2.5])
reflector_sigmaF = reflector_production / reflector_nu
reflector_sigmaC = (reflector_sigmaT - reflector_sigmaF - reflector_sigmaS.sum(axis=0))


#### void ####
void_sigmaT= np.array([1.28407e-2, 2.40997e-5])
void_absorption = np.array([4.65132e-5, 1.3289e-3])
void_production = np.array([0.0, 0.0])
void_sigmaS = np.array([[1.277e-2, 0.0], 
                        [2.40997e-5, 1.07387e-2]])

void_nu = np.array([2.5, 2.5])
void_sigmaF = void_production / void_nu
void_sigmaC = (void_sigmaT - void_sigmaF - void_sigmaS.sum(axis=0))




# =============================================================================
# Set Matierals
# =============================================================================

m_core = mcdc.material(
    capture=core_sigmaC,
    scatter=core_sigmaS,
    fission=core_sigmaF,
    nu_p=core_nu,
    chi_p=chi,
)

m_reflector = mcdc.material(
    capture=reflector_sigmaC,
    scatter=reflector_sigmaS,
    fission=reflector_sigmaF,
    nu_p=reflector_nu,
    chi_p=chi,
)

m_rod = mcdc.material(
    capture=rod_sigmaC,
    scatter=rod_sigmaS,
    fission=rod_sigmaF,
    nu_p=rod_nu,
    chi_p=chi,
)

m_void = mcdc.material(
    capture=void_sigmaC,
    scatter=void_sigmaS,
    fission=void_sigmaF,
    nu_p=void_nu,
    chi_p=chi,
)


# =============================================================================
# Set Model
# =============================================================================

# x-surfaces
sx0 = mcdc.surface("plane-x", x=0.0, bc="reflective")
sx1 = mcdc.surface("plane-x", x=15.0)
sx2 = mcdc.surface("plane-x", x=20.0)
sx3 = mcdc.surface("plane-x", x=25.0, bc="vacuum")

# y-surfaces
sy0 = mcdc.surface("plane-y", y=0.0, bc="reflective")
sy1 = mcdc.surface("plane-y", y=5.0)
sy2 = mcdc.surface("plane-y", y=15.0)
sy3 = mcdc.surface("plane-y", y=25.0, bc="vacuum")

# z-surfaces
sz0 = mcdc.surface("plane-z", z=0.0, bc="reflective")
sz1 = mcdc.surface("plane-z", z=15.0)
sz2 = mcdc.surface("plane-z", z=25.0, bc="vacuum")

# set cells
if case==1:
    # case 1 : control rod removed
    mcdc.cell([+sx1, -sx2, +sy0, -sy1, +sz0, -sz2], m_void)
elif case==2:
    # case 2 : control rod inserted
    mcdc.cell([+sx1, -sx2, +sy0, -sy1, +sz0, -sz2], m_rod)
# core
mcdc.cell([+sx0, -sx1, +sy0, -sy2, +sz0, -sz1], m_core)

# reflecting regions
mcdc.cell([+sx0, -sx1, +sy0, -sy2, +sz1, -sz2], m_reflector)
mcdc.cell([+sx2, -sx3, +sy0, -sy1, +sz0, -sz2], m_reflector)
mcdc.cell([+sx1, -sx3, +sy1, -sy3, +sz0, -sz2], m_reflector)
mcdc.cell([+sx0, -sx3, +sy2, -sy3, +sz0, -sz2], m_reflector)

# =============================================================================
# iQMC Settings
# =============================================================================
# Parse arguments
parser = argparse.ArgumentParser(description="MC/DC: Monte Carlo Dynamic Code")
parser.add_argument(
    "--N", type=int, help="Particle Histories"
)
# parser.add_argument(
#     "--Nx", type=int, help="Spatial Cells in x-direction"
# )
args, unargs = parser.parse_known_args()
# N = args.N
N = 1000

Nx = Ny = Nz = 25
maxit = 10
tol = 1e-3
G = 2
x = np.linspace(0, 25, num= Nx + 1)
y = np.linspace(0, 25, num= Ny + 1)
z = np.linspace(0, 25, num= Nz + 1)
generator = "halton"
eig_solver = "power_iteration"
fixed_source_solver = "gmres"
restart = 100


with h5py.File("iqmc_N_per_Nc=800.h5", "r") as f:
    phi0 = f["iqmc/flux"][:]
#     source0 = f["iqmc/source/constant"][:]
    fission_source0=f["iqmc/fission_source"][:]
    k_eff = f["k_eff"][()]
# phi0 = np.ones((G, Nx, Ny, Nz))

mcdc.iQMC(
    x=x,
    y=y,
    z=z,
    g=np.ones(G),
    phi0=phi0,
    fixed_source=np.zeros_like(phi0),
    # source0=source0,
    fission_source0=fission_source0,
    maxitt=maxit,
    tol=tol,
    generator=generator,
    eigenmode_solver=eig_solver,
    fixed_source_solver=fixed_source_solver,
    krylov_restart=restart,
    score=['tilt-x', 'tilt-y', 'tilt-z']
    # source_tilt=1,
)

# Setting
mcdc.setting(N_particle=N, output_name='test')#, output_name='test')
mcdc.eigenmode(k_init=k_eff)
# Run
mcdc.run()

# =============================================================================
# MC Settings
# =============================================================================
# Uniform isotropic source throughout the domain

# mcdc.source(x=[0.0, 25.0], y=[0.0, 25.0], z=[0.0, 25.0], energy=np.ones(2),isotropic=True)

# # # # Tally: cell-average and cell-edge angular fluxes and currents
# mcdc.tally(scores=["flux"], 
#             x=np.linspace(0.0, 25.0, 25), 
#             y=np.linspace(0.0, 25.0, 25),
#             z=np.linspace(0.0, 25.0, 25))

# # # # Setting
# N = 5e5
# mcdc.setting(N_particle=N)
# mcdc.eigenmode(N_inactive=10, N_active=20)
# mcdc.population_control()
# # # mcdc.implicit_capture()


# # # # Run
# mcdc.run()








