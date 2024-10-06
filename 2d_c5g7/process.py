#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 10:41:23 2023

@author: sampasmann

This file compares ST=0, ST=1, ST=2, & MC results

"""

import numpy as np
import h5py
from tabulate import tabulate


# =============================================================================
# Assebmley Error Script
# =============================================================================
# =============================================================================
def Get_Assembly_Flux(phi, idx, idy):
    global first, last
    return phi[first:last, idx[0] : idx[-1]+1, idy[0] : idy[-1]+1]



def Get_Assembly_Power(phi, xo, xf, x_mid, yo, yf, y_mid):
    """
    Assembly Powers
    """
    idx = np.where((x_mid >= xo) & (x_mid <= xf))[0]
    idy = np.where((-y_mid >= -yo) & (-y_mid <= -yf))[0]
    flux = Get_Assembly_Flux(phi, idx, idy)
    return flux.sum()

    
# =============================================================================
# Load Results
# =============================================================================
norm_factor = 1056

k_ref = 1.186550
Assembly_I_ref = 492.8
Assembly_II_ref = 211.7
Assembly_III_ref = 211.7
Assembly_IV_ref = 139.8

# # MC Reference
with h5py.File("reference.h5", "r") as f:
    phi_avg = f["tally/flux/mean"][:]
    phi_sd = f["tally/flux/sdev"][:]
    x = f["tally/grid/x"][:]
    y = f["tally/grid/y"][:]
    x_mid = 0.5 * (x[1:] + x[:-1])
    y_mid = 0.5 * (y[1:] + y[:-1])
    N = f['input_deck/setting/N_particle'][()]
    f.close()
    

# # Monte Carlo results
# with h5py.File("results/mc.h5", "r") as f:
#     phi0 = f["tally/flux/mean"][:]
#     phi0_sd = f["tally/flux/sdev"][:]
#     x = f["tally/grid/x"][:]
#     y = f["tally/grid/y"][:]
#     x_mid0 = 0.5 * (x[1:] + x[:-1])
#     y_mid0 = 0.5 * (y[1:] + y[:-1])
#     N = f['input_deck/setting/N_particle'][()]
#     f.close()

    
# iQMC results
with h5py.File("results/51x51/Nr=200_st=0.h5", "r") as f:
    # phi1 = f["iqmc/flux"][...]
    phi1 = f["iqmc/fission_power"][:,0,:,:,0]
    phi1 *= norm_factor / phi1.sum()
    x = f["iqmc/grid/x"][:]
    y = f["iqmc/grid/y"][:]
    x_mid1 = 0.5 * (x[1:] + x[:-1])
    y_mid1 = 0.5 * (y[1:] + y[:-1])
    N1 = f['input_deck/setting/N_particle'][()]
    k1 = f['k_eff'][()]
    f.close()

with h5py.File("results/51x51/Nr=200_st=1.h5", "r") as f:
    # phi1 = f["iqmc/flux"][...]
    phi2 = f["iqmc/fission_power"][:,0,:,:,0]
    phi2 *= norm_factor / phi2.sum()
    x = f["iqmc/grid/x"][:]
    y = f["iqmc/grid/y"][:]
    x_mid2 = 0.5 * (x[1:] + x[:-1])
    y_mid2 = 0.5 * (y[1:] + y[:-1])
    N1 = f['input_deck/setting/N_particle'][()]
    k2 = f['k_eff'][()]
    f.close()


with h5py.File("results/102x102/Nr=200_st=0.h5", "r") as f:
    # phi1 = f["iqmc/flux"][...]
    phi3 = f["iqmc/fission_power"][:,0,:,:,0]
    phi3 *= norm_factor / phi3.sum()
    x = f["iqmc/grid/x"][:]
    y = f["iqmc/grid/y"][:]
    x_mid3 = 0.5 * (x[1:] + x[:-1])
    y_mid3 = 0.5 * (y[1:] + y[:-1])
    N1 = f['input_deck/setting/N_particle'][()]
    k3 = f['k_eff'][()]
    f.close()

with h5py.File("results/102x102/Nr=200_st=1.h5", "r") as f:
    # phi1 = f["iqmc/flux"][...]
    phi4 = f["iqmc/fission_power"][:,0,:,:,0]
    phi4 *= norm_factor / phi4.sum()
    x = f["iqmc/grid/x"][:]
    y = f["iqmc/grid/y"][:]
    x_mid4 = 0.5 * (x[1:] + x[:-1])
    y_mid4 = 0.5 * (y[1:] + y[:-1])
    N1 = f['input_deck/setting/N_particle'][()]
    k4 = f['k_eff'][()]
    f.close()

with h5py.File("results/204x204/Nr=200_st=0.h5", "r") as f:
    # phi1 = f["iqmc/flux"][...]
    phi5 = f["iqmc/fission_power"][:,0,:,:,0]
    phi5 *= norm_factor / phi5.sum()
    x = f["iqmc/grid/x"][:]
    y = f["iqmc/grid/y"][:]
    x_mid5 = 0.5 * (x[1:] + x[:-1])
    y_mid5 = 0.5 * (y[1:] + y[:-1])
    N1 = f['input_deck/setting/N_particle'][()]
    k5 = f['k_eff'][()]
    f.close()

with h5py.File("results/204x204/Nr=400_st=2.h5", "r") as f:
    # phi1 = f["iqmc/flux"][...]
    phi6 = f["iqmc/fission_power"][:,0,:,:,0]
    phi6 *= norm_factor / phi6.sum()
    x = f["iqmc/grid/x"][:]
    y = f["iqmc/grid/y"][:]
    x_mid6 = 0.5 * (x[1:] + x[:-1])
    y_mid6 = 0.5 * (y[1:] + y[:-1])
    N1 = f['input_deck/setting/N_particle'][()]
    k6 = f['k_eff'][()]
    f.close()
    
# with h5py.File("results/408x408/Nr=200_st=0.h5", "r") as f:
#     # phi1 = f["iqmc/flux"][...]
#     phi7 = f["iqmc/fission_power"][:,0,:,:,0]
#     phi7 *= norm_factor / phi7.sum()
#     x = f["iqmc/grid/x"][:]
#     y = f["iqmc/grid/y"][:]
#     x_mid7 = 0.5 * (x[1:] + x[:-1])
#     y_mid7 = 0.5 * (y[1:] + y[:-1])
#     N1 = f['input_deck/setting/N_particle'][()]
    # k7 = f['k_eff'][()]
#     f.close()

# with h5py.File("results/408x408/Nr=200_st=1.h5", "r") as f:
#     # phi1 = f["iqmc/flux"][...]
#     phi8 = f["iqmc/fission_power"][:,0,:,:,0]
#     phi8 *= norm_factor / phi8.sum()
#     x = f["iqmc/grid/x"][:]
#     y = f["iqmc/grid/y"][:]
#     x_mid8 = 0.5 * (x[1:] + x[:-1])
#     y_mid8 = 0.5 * (y[1:] + y[:-1])
#     N1 = f['input_deck/setting/N_particle'][()]
    # k8 = f['k_eff'][()]
#     f.close()

# =============================================================================
# Geometry
# =============================================================================
pitch = 1.26
radius = 0.54

# surfaces
x0 = 0.0
x1 = pitch * 17
x2 = pitch * 17 * 2
x3 = pitch * 17 * 3

y0 = 0.0
y1 = -pitch * 17
y2 = -pitch * 17 * 2
y3 = -pitch * 17 * 3

first = 0
last = 7

# =============================================================================
# Plot to make sure I'm grabbing the right indices of the flux
# =============================================================================

import matplotlib.pyplot as plt
xo = x0
xf = x3
yo = y0
yf = y3
idx = np.where((x_mid1 >= xo) & (x_mid1 <= xf))[0]
idy = np.where((-y_mid1 >= -yo) & (-y_mid1 <= -yf))[0]
X,Y = np.meshgrid(x_mid1[idx], y_mid1[idy])
X = X.T
Y = Y.T
flux = Get_Assembly_Flux(phi1, idx, idy)[0:-1,...].sum(axis=0)

plt.figure(dpi=300, figsize=(8, 4))
plt.pcolormesh(X, Y, flux, shading="nearest")
# plt.colorbar().set_label(r"Scalar Flux", rotation=270, labelpad=15)
ax = plt.gca()
ax.set_aspect("equal")
plt.xlabel(r"$x$ [cm]")
plt.ylabel(r"$y$ [cm]")
plt.title(r"Thermal Neutron Flux")

# =============================================================================
# Assembly I
# =============================================================================
# x-dimensions
xo = x0
xf = x1
# y-dimensions
yo = y0
yf = y1

# Assembly_I_ref = Get_Assembly_Power(phi_avg, xo, xf, x_mid, yo, yf, y_mid)
# Assembly_I_ref_sd = Get_Assembly_Power(phi_sd, xo, xf, x_mid, yo, yf, y_mid)

# Assembly_I_mc = Get_Assembly_Power(phi0, xo, xf, x_mid0, yo, yf, y_mid0)
# Assembly_I_mc_sd = Get_Assembly_Power(phi0_sd, xo, xf, x_mid0, yo, yf, y_mid0)

Assembly_I_1 = Get_Assembly_Power(phi1, xo, xf, x_mid1, yo, yf, y_mid1)
Assembly_I_2 = Get_Assembly_Power(phi2, xo, xf, x_mid2, yo, yf, y_mid2)
Assembly_I_3 = Get_Assembly_Power(phi3, xo, xf, x_mid3, yo, yf, y_mid3)
Assembly_I_4 = Get_Assembly_Power(phi4, xo, xf, x_mid4, yo, yf, y_mid4)
Assembly_I_5 = Get_Assembly_Power(phi5, xo, xf, x_mid5, yo, yf, y_mid5)
Assembly_I_6 = Get_Assembly_Power(phi6, xo, xf, x_mid6, yo, yf, y_mid6)
# Assembly_I_7 = Get_Assembly_Power(phi7, xo, xf, x_mid7, yo, yf, y_mid7)
# Assembly_I_8 = Get_Assembly_Power(phi8, xo, xf, x_mid8, yo, yf, y_mid8)


# =============================================================================
# Assembly II
# =============================================================================
# x-dimensions
xo = x1
xf = x2
# y-dimensions
yo = y0
yf = y1

# Assembly_II_ref = Get_Assembly_Power(phi_avg, xo, xf, x_mid, yo, yf, y_mid)
# Assembly_II_ref_sd = Get_Assembly_Power(phi_sd, xo, xf, x_mid, yo, yf, y_mid)

# Assembly_II_mc = Get_Assembly_Power(phi0, xo, xf, x_mid0, yo, yf, y_mid0)
# Assembly_II_mc_sd = Get_Assembly_Power(phi0_sd, xo, xf, x_mid0, yo, yf, y_mid0)

Assembly_II_1 = Get_Assembly_Power(phi1, xo, xf, x_mid1, yo, yf, y_mid1)
Assembly_II_2 = Get_Assembly_Power(phi2, xo, xf, x_mid2, yo, yf, y_mid2)
Assembly_II_3 = Get_Assembly_Power(phi3, xo, xf, x_mid3, yo, yf, y_mid3)
Assembly_II_4 = Get_Assembly_Power(phi4, xo, xf, x_mid4, yo, yf, y_mid4)
Assembly_II_5 = Get_Assembly_Power(phi5, xo, xf, x_mid5, yo, yf, y_mid5)
Assembly_II_6 = Get_Assembly_Power(phi6, xo, xf, x_mid6, yo, yf, y_mid6)
# Assembly_II_7 = Get_Assembly_Power(phi7, xo, xf, x_mid7, yo, yf, y_mid7)
# Assembly_II_8 = Get_Assembly_Power(phi8, xo, xf, x_mid8, yo, yf, y_mid8)


# =============================================================================
# Assembly III
# =============================================================================
# x-dimensions
xo = x0
xf = x1
# y-dimensions
yo = y1
yf = y2

# Assembly_III_ref = Get_Assembly_Power(phi_avg, xo, xf, x_mid, yo, yf, y_mid)
# Assembly_III_ref_sd = Get_Assembly_Power(phi_sd, xo, xf, x_mid, yo, yf, y_mid)

# Assembly_III_mc = Get_Assembly_Power(phi0, xo, xf, x_mid0, yo, yf, y_mid0)
# Assembly_III_mc_sd = Get_Assembly_Power(phi0_sd, xo, xf, x_mid0, yo, yf, y_mid0)

Assembly_III_1 = Get_Assembly_Power(phi1, xo, xf, x_mid1, yo, yf, y_mid1)
Assembly_III_2 = Get_Assembly_Power(phi2, xo, xf, x_mid2, yo, yf, y_mid2)
Assembly_III_3 = Get_Assembly_Power(phi3, xo, xf, x_mid3, yo, yf, y_mid3)
Assembly_III_4 = Get_Assembly_Power(phi4, xo, xf, x_mid4, yo, yf, y_mid4)
Assembly_III_5 = Get_Assembly_Power(phi5, xo, xf, x_mid5, yo, yf, y_mid5)
Assembly_III_6 = Get_Assembly_Power(phi6, xo, xf, x_mid6, yo, yf, y_mid6)
# Assembly_III_7 = Get_Assembly_Power(phi7, xo, xf, x_mid7, yo, yf, y_mid7)
# Assembly_III_8 = Get_Assembly_Power(phi8, xo, xf, x_mid8, yo, yf, y_mid8)


# =============================================================================
# Assembly IV
# =============================================================================
# x-dimensions
xo = x1
xf = x2
# y-dimensions
yo = y1
yf = y2

# Assembly_IV_ref = Get_Assembly_Power(phi_avg, xo, xf, x_mid, yo, yf, y_mid)
# Assembly_IV_ref_sd = Get_Assembly_Power(phi_sd, xo, xf, x_mid, yo, yf, y_mid)

# Assembly_IV_mc = Get_Assembly_Power(phi0, xo, xf, x_mid0, yo, yf, y_mid0)
# Assembly_IV_mc_sd = Get_Assembly_Power(phi0_sd, xo, xf, x_mid0, yo, yf, y_mid0)

Assembly_IV_1 = Get_Assembly_Power(phi1, xo, xf, x_mid1, yo, yf, y_mid1)
Assembly_IV_2 = Get_Assembly_Power(phi2, xo, xf, x_mid2, yo, yf, y_mid2)
Assembly_IV_3 = Get_Assembly_Power(phi3, xo, xf, x_mid3, yo, yf, y_mid3)
Assembly_IV_4 = Get_Assembly_Power(phi4, xo, xf, x_mid4, yo, yf, y_mid4)
Assembly_IV_5 = Get_Assembly_Power(phi5, xo, xf, x_mid5, yo, yf, y_mid5)
Assembly_IV_6 = Get_Assembly_Power(phi6, xo, xf, x_mid6, yo, yf, y_mid6)
# Assembly_IV_7 = Get_Assembly_Power(phi7, xo, xf, x_mid7, yo, yf, y_mid7)
# Assembly_IV_8 = Get_Assembly_Power(phi8, xo, xf, x_mid8, yo, yf, y_mid8)


# =============================================================================
# k-Eigenvalue
# =============================================================================
D = 3
table = []

def error(k, k_ref):
    return abs(k - k_ref) / k_ref * 100

err1 = error(k1, k_ref)
err2 = error(k2, k_ref)
err3 = error(k3, k_ref)
err4 = error(k4, k_ref)
err5 = error(k5, k_ref)
err6 = error(k6, k_ref)
# err7 = error(k7, k_ref)
# err8 = error(k8, k_ref)

row = ["MCNP",
       "---",
       "---",
       k_ref,
       ""]
table.append(row)

row = ["iQMC",
       "51x51",
       "Constant",
        k1,
        np.round(err1,  decimals=D)]
table.append(row)

row = ["iQMC",
       "51x51",
       "Linear",
        k2,
        np.round(err2,  decimals=D)]
table.append(row)

row = ["iQMC",
       "102x102",
       "Constant",
        k3,
        np.round(err3,  decimals=D)]
table.append(row)

row = ["iQMC",
       "102x102",
       "Linear",
        k4,
        np.round(err4,  decimals=D)]
table.append(row)

row = ["iQMC",
       "204x204",
       "Constant",
        k5,
        np.round(err5,  decimals=D)]
table.append(row)

row = ["iQMC",
       "204x204",
       "Linear",
        k6,
        np.round(err6,  decimals=D)]
table.append(row)

# row = ["iQMC",
#        "408x408",
#        "Constant",
#         k7,
#         np.round(err7,  decimals=D)]
# table.append(row)

# row = ["iQMC",
#        "408x408",
#        "Linear",
#         k8,
#         np.round(err8,  decimals=D)]
# table.append(row)

print(
    tabulate(
        table,
        headers=[
            "Code",
            "Mesh",
            "Source",
            "Eigenvalue",
            "Percent Error",
        ], tablefmt="simple"
    )
)


# =============================================================================
# Assembly Power Table
# =============================================================================
D = 5

table = []


row = ["MCNP",
       "---",
       "---",
        np.round(Assembly_I_ref, decimals=D),
        np.round(Assembly_II_ref, decimals=D),
        np.round(Assembly_III_ref, decimals=D),
        np.round(Assembly_IV_ref, decimals=D)]
table.append(row)

# row = ["MC Ref", 
#         u"{} \u00B1 {}".format(np.round(Assembly_I_ref,decimals=D),np.round(Assembly_I_ref_sd,decimals=D)),
#         u"{} \u00B1 {}".format(np.round(Assembly_II_ref,decimals=D),np.round(Assembly_II_ref_sd,decimals=5)),
#         u"{} \u00B1 {}".format(np.round(Assembly_III_ref,decimals=D),np.round(Assembly_III_ref_sd,decimals=D)),
#         u"{} \u00B1 {}".format(np.round(Assembly_IV_ref,decimals=D),np.round(Assembly_IV_ref_sd,decimals=D))]
# table.append(row)


# row = ["MC", 
#         u"{} \u00B1 {}".format(np.round(Assembly_I_mc,decimals=D),np.round(Assembly_I_mc_sd,decimals=D)),
#         u"{} \u00B1 {}".format(np.round(Assembly_II_mc,decimals=D),np.round(Assembly_II_mc_sd,decimals=5)),
#         u"{} \u00B1 {}".format(np.round(Assembly_III_mc,decimals=D),np.round(Assembly_III_mc_sd,decimals=D)),
#         u"{} \u00B1 {}".format(np.round(Assembly_IV_mc,decimals=D),np.round(Assembly_IV_mc_sd,decimals=D))]
# table.append(row)


row = ["iQMC",
       "51x51",
       "Constant",
        np.round(Assembly_I_1, decimals=D),
        np.round(Assembly_II_1, decimals=D),
        np.round(Assembly_III_1, decimals=D),
        np.round(Assembly_IV_1, decimals=D)]
table.append(row)


row = ["iQMC",
       "51x51",
       "Linear",
        np.round(Assembly_I_2, decimals=D),
        np.round(Assembly_II_2, decimals=D),
        np.round(Assembly_III_2, decimals=D),
        np.round(Assembly_IV_2, decimals=D)]
table.append(row)

row = ["iQMC",
       "102x102",
       "Constant",
        np.round(Assembly_I_3, decimals=D),
        np.round(Assembly_II_3, decimals=D),
        np.round(Assembly_III_3, decimals=D),
        np.round(Assembly_IV_3, decimals=D)]
table.append(row)

row = ["iQMC",
       "102x102",
       "Linear",
        np.round(Assembly_I_4, decimals=D),
        np.round(Assembly_II_4, decimals=D),
        np.round(Assembly_III_4, decimals=D),
        np.round(Assembly_IV_4, decimals=D)]
table.append(row)

row = ["iQMC",
       "204x204",
       "Constant",
        np.round(Assembly_I_5, decimals=D),
        np.round(Assembly_II_5, decimals=D),
        np.round(Assembly_III_5, decimals=D),
        np.round(Assembly_IV_5, decimals=D)]
table.append(row)

row = ["iQMC",
       "204x204",
       "Linear",
        np.round(Assembly_I_6, decimals=D),
        np.round(Assembly_II_6, decimals=D),
        np.round(Assembly_III_6, decimals=D),
        np.round(Assembly_IV_6, decimals=D)]
table.append(row)

# row = [
        # "iQMC",
        # "408x408",
#        "Constant",
#         np.round(Assembly_I_7, decimals=D),
#         np.round(Assembly_II_7, decimals=D),
#         np.round(Assembly_III_7, decimals=D),
#         np.round(Assembly_IV_7, decimals=D)]
# table.append(row)

# row = [
        # "iQMC",
        # "408x408",
#        "Linear",
#         np.round(Assembly_I_8, decimals=D),
#         np.round(Assembly_II_8, decimals=D),
#         np.round(Assembly_III_8, decimals=D),
#         np.round(Assembly_IV_8, decimals=D)]
# table.append(row)

print(
    tabulate(
        table,
        headers=[
            "Code",
            "Mesh",
            "Source",
            "Inner UO2",
            "MOX I",
            "MOX II",
            "Outter UO2",
        ], tablefmt="simple"
    )
)
