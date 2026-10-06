
import numpy as np
import matplotlib.pyplot as plt

# Constants
kB = 1.380649e-23       # Boltzmann constant [J/K]
T = 300.0               # Temperature [K]
m = 6.63e-26            # Argon molecular mass [kg]

# Load simulated speeds
speed = np.loadtxt("speed.dat")

# Theoretical Maxwell-Boltzmann speed distribution
v = np.linspace(0, speed.max(), 1000)

f = 4 * np.pi * (m / (2 * np.pi * kB * T))**1.5 \
    * v**2 * np.exp(-m * v**2 / (2 * kB * T))

# Plot simulation
plt.hist(
    speed,
    bins=100,
    density=True,
    alpha=0.6,
    label="Simulation"
)

# Plot theory
plt.plot(
    v,
    f,
    linewidth=2,
    label="Maxwell-Boltzmann"
)

plt.xlabel("Speed [m/s]")
plt.ylabel("Probability density")
plt.title("Argon speed distribution at 300 K")
plt.legend()
plt.grid(alpha=0.2)

plt.show()
