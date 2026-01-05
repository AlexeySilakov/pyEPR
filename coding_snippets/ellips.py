import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ================================================================
# Parameters
# ================================================================
# User-provided, but reordered so c > b > a
params = np.array([3.0, 2.0, 1.5])
a, b, c = np.sort(params)
R = 1.9
nPoints = 100

print(f"Using parameters ordered as: a={a}, b={b}, c={c}")

# ================================================================
# Analytical intersection curve
# ================================================================
alpha = 1/a**2 - 1/c**2
beta  = 1/b**2 - 1/c**2
gamma = 1 - (R**2)/c**2

if alpha <= 0 or beta <= 0 or gamma <= 0:
    raise RuntimeError("Parameters produce no real intersection curve.")

X = np.sqrt(gamma / alpha)
Y = np.sqrt(gamma / beta)

# Parametric curve r(t)
t = np.linspace(0, np.pi/2, nPoints)

x = X * np.cos(t)
y = Y * np.sin(t)
z = np.sqrt(R**2 - x**2 - y**2)
#z[t > np.pi] *= -1   # continuous closed curve

dt = np.pi/2 / (nPoints-1)
x1 = X * np.cos(t+dt)
y1 = Y * np.sin(t+dt)
z1 = np.sqrt(R**2 - x1**2 - y1**2)

dx = x1-x
dy = y1-y
dz = z1-z

print(np.sum( np.sqrt((dx**2+dy**2+dz**2))))

st, ct = np.sin(t), np.cos(t)
numerator = np.abs(R**2*(X**2*st**2 + Y**2*ct**2) - X**2*Y**2)
denominator = np.abs(R**2 - X**2*ct**2 - Y**2*st**2)


f = np.sqrt(numerator / denominator)
dt = np.diff(t)

# CUMSUM version of trapezoid rule:
#   ∫ f(t) dt ≈ Σ 0.5 * (f_i + f_{i+1}) * Δt
# integrand_trap = 0.5 * (f[:-1] + f[1:]) * dt
#L_analytical = np.cumsum(integrand_trap)[-1]

L_analytical = np.sum(f[1:]*dt)

# ================================================================
# Numerical intersection (grid sampling)
# ================================================================
gridN = 200
xs = np.linspace(-R, R, gridN)
ys = np.linspace(-R, R, gridN)

pts = []

for xi in xs:
    for yi in ys:
        inside = R**2 - xi**2 - yi**2
        if inside < 0:
            continue
        for zi in (np.sqrt(inside), -np.sqrt(inside)):
            if abs((xi**2)/a**2 + (yi**2)/b**2 + (zi**2)/c**2 - 1) < 0.002:
                pts.append([xi, yi, zi])

pts = np.array(pts)

# Sort and compute numerical length
angles = np.arctan2(pts[:,1], pts[:,0])
order = np.argsort(angles)
pts_sorted = pts[order]

d = np.sqrt(np.sum(np.diff(pts_sorted, axis=0)**2, axis=1))
L_numeric = d.sum()

# ================================================================
# Print comparison
# ================================================================
print("\n===== CURVE LENGTH COMPARISON =====")
print(f"Analytical (cumsum): {L_analytical:.6f}")
print(f"Point-cloud numerical: {L_numeric:.6f}")
print("Difference:", abs(L_analytical - L_numeric))

# ================================================================
# Plot sphere, ellipsoid, and intersections
# ================================================================
fig = plt.figure(figsize=(10, 9))
ax = fig.add_subplot(111, projection='3d')

# Sphere
u = np.linspace(0, 2*np.pi, 30)
v = np.linspace(0, np.pi, 15)
xs = R*np.cos(u)[:,None]*np.sin(v)[None,:]
ys = R*np.sin(u)[:,None]*np.sin(v)[None,:]
zs = R*np.ones_like(xs)*np.cos(v)[None,:]
ax.plot_wireframe(xs, ys, zs, color='gray', alpha=0.4)

# Ellipsoid
xe = a*np.cos(u)[:,None]*np.sin(v)[None,:]
ye = b*np.sin(u)[:,None]*np.sin(v)[None,:]
ze = c*np.ones_like(xe)*np.cos(v)[None,:]
ax.plot_wireframe(xe, ye, ze, color='blue', alpha=0.35)

# Intersection curves
ax.plot(x, y, z, 'r', linewidth=2, label="Analytical")
ax.scatter(pts[:,0], pts[:,1], pts[:,2], s=4, c='k', label="Numerical")

ax.legend()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Sphere / Ellipsoid Intersection")
ax.set_box_aspect([1,1,1])
plt.show()
