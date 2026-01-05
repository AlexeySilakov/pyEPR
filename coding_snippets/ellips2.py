import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Ellipsoid and sphere parameters
a = 3.0  # x-axis radius of ellipsoid
b = 2.0  # y-axis radius of ellipsoid
c = 1.5  # z-axis radius of ellipsoid
R = 2.5  # sphere radius

# Precompute constants
alpha = 1/a**2 - 1/c**2
beta  = 1/b**2 - 1/c**2
gamma = 1 - R**2/c**2

# Parametric equations for intersection
def x_t(t): return np.sqrt(gamma/alpha) * np.cos(t)
def y_t(t): return np.sqrt(gamma/beta) * np.sin(t)
def z_t(t): return np.sqrt(R**2 - x_t(t)**2 - y_t(t)**2)

# Derivative of z
def dz_dt_abs(t):
    numerator = gamma * (1/alpha - 1/beta) * np.cos(t) * np.sin(t)
    denominator = np.sqrt(R**2 - gamma/alpha * np.cos(t)**2 - gamma/beta * np.sin(t)**2)
    return np.abs(numerator / denominator)

# Arc length integrand
def integrand(t):
    dx_dt2 = np.abs((-np.sqrt(gamma/alpha) * np.sin(t))**2)
    dy_dt2 = np.abs(( np.sqrt(gamma/beta) * np.cos(t))**2)
    dz_dt2 = dz_dt_abs(t)**2
    return np.sqrt(dx_dt2 + dy_dt2 + dz_dt2)

# Integrate over [0, pi/2] and multiply by 4
L, _ = quad(integrand, 0, np.pi/2)
L_total = 4 * L
print(f"Estimated length of intersection curve: {L_total:.6f}")

# Create t values for plotting
t_vals = np.linspace(0, 2*np.pi, 400)
x_vals = x_t(t_vals)
y_vals = y_t(t_vals)
z_vals = z_t(t_vals)

# 3D plotting
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')

# Sphere mesh
u = np.linspace(0, np.pi, 50)
v = np.linspace(0, 2*np.pi, 50)
xs = R * np.sin(u)[:,None] * np.cos(v)
ys = R * np.sin(u)[:,None] * np.sin(v)
zs = R * np.cos(u)[:,None] * np.ones_like(v)
ax.plot_surface(xs, ys, zs, color='cyan', alpha=0.3, edgecolor='k')

# Ellipsoid mesh
xe = a * np.sin(u)[:,None] * np.cos(v)
ye = b * np.sin(u)[:,None] * np.sin(v)
ze = c * np.cos(u)[:,None] * np.ones_like(v)
ax.plot_surface(xe, ye, ze, color='magenta', alpha=0.3, edgecolor='k')

# Intersection curve
ax.plot(x_vals, y_vals, z_vals, color='red', linewidth=2, label='Intersection')

# Labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Intersection of Sphere and Ellipsoid')
ax.legend()
ax.set_box_aspect([1,1,1])

plt.show()
