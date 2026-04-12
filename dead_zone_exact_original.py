import numpy as np
import matplotlib.pyplot as plt

# Critical Thiele modulus (from paper)
def mu_star(n):
    return np.sqrt(2.0 * (n + 1.0)) / (1.0 - n)

# Parameters
n = 0.5          # reaction exponent (0 < n < 1)
phi = 6.0        # Thiele modulus

phi_star = mu_star(n)

phi=phi_star

print(f"n={n}, phi={phi}, phi*={phi_star:.6f}")

# Dead-zone length
xdz = 1 - (phi_star / phi)

# Domain
x = np.linspace(0, 1, 500)

# Exact solution
u = np.zeros_like(x)

# Compute solution
mask = x > xdz
u[mask] = ((x[mask] - xdz) / (1 - xdz)) ** (2 / (1 - n))

# Plot
plt.figure()
plt.plot(x, u, label="Exact solution")
plt.axvline(xdz, linestyle="--", label=f"Dead-zone boundary x_dz={xdz:.6f}")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Exact Solution with Dead Zone")
plt.legend()
plt.grid()

plt.show()