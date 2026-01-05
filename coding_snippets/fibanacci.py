#!/usr/bin/env python3
"""
Fibonacci lattice point distribution on a sphere (or hemisphere / quarter‑hemisphere).

Author:  ChatGPT
Promnt:
    make a python numpy-based script that distributes points as evenly as possible on a sphere using Fibonacci Lattice method. write the core math directly in phi and theta. Assume hemisphere is that of z>0 and quarter of hemisphere is z>0, y>0 and x>0. for a demo set default variables of N=20 and r=1. provide list of phi and theta values.  add a 3D plot using matplotlib showing the resulting grid with equal axes representation. Add an estimate of uniformity of points and the delta parameter
Date:    2025-11-23
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401  (needed for 3D projection)
from itertools import combinations

# --------------------------------------------------------------------------- #
# Core routine – Fibonacci lattice on a sphere/hemisphere/quarter‑hemisphere
# --------------------------------------------------------------------------- #

def fibonacci_lattice(N, r=1.0, hemisphere=False, quarter=False):
    """
    Generate N points on the surface of a sphere (radius r) using the Fibonacci lattice.
    
    Parameters
    ----------
    N : int
        Number of points to generate.
    r : float, optional
        Radius of the sphere. Default is 1.0.
    hemisphere : bool, optional
        If True, only keep points with z > 0 (upper hemisphere).
    quarter : bool, optional
        If True, further restrict to x>0, y>0, z>0 (one quarter of the upper hemisphere).
    
    Returns
    -------
    phi   : ndarray shape (N,)
        Azimuthal angles in [0, 2π].
    theta : ndarray shape (N,)
        Polar angles in [0, π] (measured from +z axis).
    xyz   : ndarray shape (N,3)
        Cartesian coordinates of the points.
    """
    # Golden angle (in radians) – the key to even spacing
    golden_angle = np.pi * (3. - np.sqrt(5.))  # ≈ 2.399963229
    epsilon = 0.5
    # --------------------------------------------------------------------- #
    # 1. Compute y (cos(theta)) values
    # --------------------------------------------------------------------- #
    if hemisphere:
        # Upper hemisphere: y ∈ [0, 1]
        y = (np.arange(N) + epsilon) / (N+2*epsilon)          # uniform in [0,1)
    else:
        # Full sphere: y ∈ [-1, 1]
        y = 1 - (np.arange(N) + epsilon) * 2 / (N+2*epsilon)   # uniform in (-1,1)

    theta = np.arccos(y)                      # polar angle

    # --------------------------------------------------------------------- #
    # 2. Compute phi values
    # --------------------------------------------------------------------- #
    phi = (golden_angle * np.arange(N)) % (2*np.pi)

    # --------------------------------------------------------------------- #
    # 3. Convert to Cartesian coordinates
    # --------------------------------------------------------------------- #
    x = r * np.sin(theta) * np.cos(phi)
    y_cart = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    xyz = np.column_stack((x, y_cart, z))

    # --------------------------------------------------------------------- #
    # 4. Apply hemisphere / quarter restrictions
    # --------------------------------------------------------------------- #
    if hemisphere:
        mask = z > 0
        phi   = phi[mask]
        theta = theta[mask]
        xyz   = xyz[mask]

    if quarter:
        mask = (xyz[:,0] > 0) & (xyz[:,1] > 0) & (xyz[:,2] > 0)
        phi   = phi[mask]
        theta = theta[mask]
        xyz   = xyz[mask]

    return phi, theta, xyz

# --------------------------------------------------------------------------- #
# Utility – uniformity estimate
# --------------------------------------------------------------------------- #

def uniformity_metrics(xyz, N):
    """
    Compute simple metrics that quantify how uniformly the points are spaced.
    
    Parameters
    ----------
    xyz : ndarray shape (M,3)
        Cartesian coordinates of M points on a unit sphere.
    
    Returns
    -------
    mean_dist   : float
        Mean great‑circle distance between all distinct pairs.
    std_dist    : float
        Standard deviation of those distances.
    min_dist    : float
        Minimum pairwise distance.
    max_dist    : float
        Maximum pairwise distance.
    delta       : float
        (max - min) / mean  – a dimensionless measure of spread.
    """
    # Normalise to unit sphere for great‑circle distance calculation
    xyz_unit = xyz / np.linalg.norm(xyz, axis=1)[:,None]

    # Compute all pairwise dot products
    dots = np.dot(xyz_unit, xyz_unit.T)
    # Clip due to numerical errors
    dots = np.clip(dots, -1.0, 1.0)

    # Great‑circle distances (in radians)
    dists = np.arccos(dots)

    # We only need the upper triangle (excluding diagonal)
    iu = np.triu_indices_from(dists, k=1)
    pair_dists = dists[iu]

    mean_dist = pair_dists.mean()
    std_dist  = pair_dists.std(ddof=0)
    min_dist  = np.sqrt(N)*pair_dists.min()
    max_dist  = pair_dists.max()
    delta     = (max_dist - min_dist) / mean_dist

    return mean_dist, std_dist, min_dist, max_dist, delta

# --------------------------------------------------------------------------- #
# Demo – default parameters
# --------------------------------------------------------------------------- #

def main():
    # Default demo values
    N   = 150          # number of points to generate
    r   = 1.0         # radius of the sphere
    hemi = True       # use upper hemisphere (z > 0)
    quarter = False   # set to True for x>0, y>0, z>0

    phi, theta, xyz = fibonacci_lattice(N=N, r=r,
                                        hemisphere=hemi,
                                        quarter=quarter)

    print(f"Generated {len(phi)} points on the {'hemisphere' if hemi else 'sphere'}")
    print("\nPhi (azimuthal) angles [rad]:")
    print(np.round(phi, 4))
    print("\nTheta (polar) angles [rad]:")
    print(np.round(theta, 4))

    # --------------------------------------------------------------------- #
    # Uniformity metrics
    # --------------------------------------------------------------------- #
    mean_dist, std_dist, min_dist, max_dist, delta = uniformity_metrics(xyz, N)

    print("\nUniformity metrics (great‑circle distances in radians):")
    print(f"  Mean distance : {mean_dist:.4f}")
    print(f"  Std. dev.     : {std_dist:.4f}")
    print(f"  Min distance  : {min_dist:.4f}")
    print(f"  Max distance  : {max_dist:.4f}")
    print(f"  Delta (spread): {delta:.4f}")

    # --------------------------------------------------------------------- #
    # Plotting
    # --------------------------------------------------------------------- #
    fig = plt.figure(figsize=(8,8))
    ax  = fig.add_subplot(111, projection='3d')
    ax.scatter(xyz[:,0], xyz[:,1], xyz[:,2],
               c=xyz[:,2], cmap='viridis', s=50, edgecolor='k')

    # Equal aspect ratio
    max_range = np.array([xyz[:,0].max()-xyz[:,0].min(),
                          xyz[:,1].max()-xyz[:,1].min(),
                          xyz[:,2].max()-xyz[:,2].min()]).max() / 2.0

    mid_x = (xyz[:,0].max()+xyz[:,0].min()) * 0.5
    mid_y = (xyz[:,1].max()+xyz[:,1].min()) * 0.5
    mid_z = (xyz[:,2].max()+xyz[:,2].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Fibonacci lattice on {"hemisphere" if hemi else "sphere"} (N={len(phi)})')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
