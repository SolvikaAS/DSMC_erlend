import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Settings
# ============================================================

filename = "DSMC2D_Field.plt"

x_target = -5.0

kB = 1.380662e-23  # Boltzmann constant [J/K]


# ============================================================
# Read DSMC data
# ============================================================

data = np.loadtxt(
    filename,
    skiprows=2
)

X = data[:, 0]
Y = data[:, 1]
N = data[:, 2]
T = data[:, 7]


# ============================================================
# Pressure
# ============================================================

P = N * kB * T


# ============================================================
# Find grid column closest to x = -5
# ============================================================

x_values = np.unique(X)

x_actual = x_values[
    np.argmin(np.abs(x_values - x_target))
]

print(f"Requested x = {x_target:.3f} m")
print(f"Using x = {x_actual:.3f} m")


# ============================================================
# Extract data at this x position
# ============================================================

mask = np.isclose(X, x_actual)

y_line = Y[mask]
T_line = T[mask]
P_line = P[mask]


# Sort by y
order = np.argsort(y_line)

y_line = y_line[order]
T_line = T_line[order]
P_line = P_line[order]


# ============================================================
# Temperature plot
# ============================================================

plt.figure(figsize=(7, 5))

plt.plot(
    y_line,
    T_line
)

plt.xlabel("y [m]")
plt.ylabel("Temperature [K]")

plt.title(
    f"Temperature at x = {x_actual:.2f} m"
)

plt.grid()

plt.tight_layout()
plt.show()


# ============================================================
# Pressure plot
# ============================================================

plt.figure(figsize=(7, 5))

plt.plot(
    y_line,
    P_line
)

plt.xlabel("y [m]")
plt.ylabel("Pressure [Pa]")

plt.title(
    f"Pressure at x = {x_actual:.2f} m"
)

plt.grid()

plt.tight_layout()
plt.show()
