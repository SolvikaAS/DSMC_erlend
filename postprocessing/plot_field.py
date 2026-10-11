from plot_cli import parse_args, finish_figure, read_run_metadata

args = parse_args("Plot DSMC density, temperature, pressure, and velocity fields.")

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


# ============================================================
# Settings
# ============================================================

filename = args.filename

kB = 1.380662e-23  # Boltzmann constant [J/K]

# Mesh visualization
mesh_step = 1


# ============================================================
# Geometry saved with the DSMC field
# ============================================================

metadata = read_run_metadata(filename)
obstacle = metadata["obstacle"]
if obstacle["type"] != "box":
    raise ValueError("This plotting script expects a box obstacle.")
BoxX, BoxY = obstacle["center"]
BoxWidth = obstacle["width"]
BoxHeight = obstacle["height"]
if "angle_of_attack_degrees" in metadata:
    print(f"Angle of attack: {metadata['angle_of_attack_degrees']:g} degrees")


# ============================================================
# Read data
# ============================================================

data = np.loadtxt(
    filename,
    skiprows=2
)

X  = data[:, 0]
Y  = data[:, 1]
N  = data[:, 2]
Ux = data[:, 3]
Uy = data[:, 4]
Tx = data[:, 5]
Ty = data[:, 6]
T  = data[:, 7]

print(f"Number of points: {len(X)}")


# ============================================================
# Remove invalid values
# ============================================================

valid = (
    np.isfinite(X)
    & np.isfinite(Y)
    & np.isfinite(N)
    & np.isfinite(Ux)
    & np.isfinite(Uy)
    & np.isfinite(Tx)
    & np.isfinite(Ty)
    & np.isfinite(T)
)

X  = X[valid]
Y  = Y[valid]
N  = N[valid]
Ux = Ux[valid]
Uy = Uy[valid]
Tx = Tx[valid]
Ty = Ty[valid]
T  = T[valid]


# ============================================================
# Grid
# ============================================================

x_values = np.unique(X)
y_values = np.unique(Y)

print(f"Unique x values: {len(x_values)}")
print(f"Unique y values: {len(y_values)}")

xi = np.linspace(
    X.min(),
    X.max(),
    len(x_values)
)

yi = np.linspace(
    Y.min(),
    Y.max(),
    len(y_values)
)

Xgrid, Ygrid = np.meshgrid(xi, yi)


# ============================================================
# Interpolation helper
# ============================================================

def make_grid(values):

    return griddata(
        (X, Y),
        values,
        (Xgrid, Ygrid),
        method="nearest"
    )


# ============================================================
# Generate field grids
# ============================================================

Ngrid  = make_grid(N)

Uxgrid = make_grid(Ux)
Uygrid = make_grid(Uy)

Txgrid = make_grid(Tx)
Tygrid = make_grid(Ty)

Tgrid  = make_grid(T)


# ============================================================
# Solid-body mask
# ============================================================

body_mask = (
    (np.abs(Xgrid - BoxX) <= BoxWidth / 2)
    & (np.abs(Ygrid - BoxY) <= BoxHeight / 2)
)


# ============================================================
# Derived quantities
# ============================================================

Ugrid = np.sqrt(
    Uxgrid**2 +
    Uygrid**2
)

Pgrid = (
    Ngrid *
    kB *
    Tgrid
)

Pxx_grid = (
    Ngrid *
    kB *
    Txgrid
)

Pyy_grid = (
    Ngrid *
    kB *
    Tygrid
)


# ============================================================
# Field ranges
# ============================================================

def outside_body_values(field):

    mask = (
        body_mask |
        ~np.isfinite(field)
    )

    result = np.ma.array(
        field,
        mask=mask
    )

    return result.compressed()


print()
print("Field ranges")
print("--------------------------------")

N_valid = outside_body_values(Ngrid)
T_valid = outside_body_values(Tgrid)
U_valid = outside_body_values(Ugrid)
P_valid = outside_body_values(Pgrid)

print(
    f"N: {np.min(N_valid):.3e} - "
    f"{np.max(N_valid):.3e} 1/m^3"
)

print(
    f"T: {np.min(T_valid):.3f} - "
    f"{np.max(T_valid):.3f} K"
)

print(
    f"U: {np.min(U_valid):.3f} - "
    f"{np.max(U_valid):.3f} m/s"
)

print(
    f"P: {np.min(P_valid):.3e} - "
    f"{np.max(P_valid):.3e} Pa"
)


# ============================================================
# Mesh
# ============================================================

def add_mesh(ax):

    """
    Draw the DSMC computational mesh.

    mesh_step = 1  -> every grid line
    mesh_step = 5  -> every fifth grid line
    """

    # Vertical grid lines
    for x in xi[::mesh_step]:

        ax.axvline(
            x,
            linewidth=0.25,
            alpha=1,
            zorder=10
        )

    # Horizontal grid lines
    for y in yi[::mesh_step]:

        ax.axhline(
            y,
            linewidth=0.25,
            alpha=0.35,
            zorder=10
        )


# ============================================================
# Generic field plotting function
# ============================================================

def field_plot(
    field,
    title,
    cbar_label,
    cmap="viridis",
    nlevels=80,
    percentiles=(1, 99)
):

    fig, ax = plt.subplots(
        figsize=(10, 7)
    )

    # --------------------------------------------------------
    # Mask solid body + invalid cells
    # --------------------------------------------------------

    plot_field = np.ma.array(
        field,
        mask=(
            body_mask |
            ~np.isfinite(field)
        )
    )

    # --------------------------------------------------------
    # Color limits
    # --------------------------------------------------------

    finite = plot_field.compressed()

    vmin = np.percentile(
        finite,
        percentiles[0]
    )

    vmax = np.percentile(
        finite,
        percentiles[1]
    )

    levels = np.linspace(
        vmin,
        vmax,
        nlevels
    )

    # --------------------------------------------------------
    # Field
    # --------------------------------------------------------

    contour = ax.contourf(
        Xgrid,
        Ygrid,
        plot_field,
        levels=levels,
        cmap=cmap,
        extend="both",
        corner_mask=False  # Preserve square corners at masked grid cells.
    )

    # --------------------------------------------------------
    # Mesh
    # --------------------------------------------------------

    add_mesh(ax)

    # --------------------------------------------------------
    # Color bar
    # --------------------------------------------------------

    cbar = fig.colorbar(
        contour,
        ax=ax,
        pad=0.02
    )

    cbar.set_label(
        cbar_label,
        fontsize=12
    )

    # --------------------------------------------------------
    # Axes
    # --------------------------------------------------------

    ax.set_xlabel(
        "x [m]",
        fontsize=12
    )

    ax.set_ylabel(
        "y [m]",
        fontsize=12
    )

    ax.set_title(
        title,
        fontsize=14
    )

    ax.set_aspect("equal")

    ax.tick_params(
        labelsize=10
    )

    fig.tight_layout()

    finish_figure(fig, args, title)


# ============================================================
# Number density
# ============================================================

field_plot(
    Ngrid,
    "DSMC Number Density",
    r"Number density $n$ [m$^{-3}$]",
    cmap="viridis"
)


# ============================================================
# Translational temperature
# ============================================================

field_plot(
    Tgrid,
    "DSMC Translational Temperature",
    "Temperature [K]",
    cmap="inferno"
)


# ============================================================
# Pressure
# ============================================================

field_plot(
    Pgrid,
    "DSMC Static Pressure",
    "Pressure [Pa]",
    cmap="viridis"
)


# ============================================================
# Tx
# ============================================================

field_plot(
    Txgrid,
    r"Translational Temperature $T_x$",
    r"$T_x$ [K]",
    cmap="inferno"
)


# ============================================================
# Ty
# ============================================================

field_plot(
    Tygrid,
    r"Translational Temperature $T_y$",
    r"$T_y$ [K]",
    cmap="inferno"
)


# ============================================================
# Velocity vector plot
# ============================================================

fig, ax = plt.subplots(
    figsize=(10, 7)
)


# ------------------------------------------------------------
# Mask solid body
# ------------------------------------------------------------

plot_U = np.ma.array(
    Ugrid,
    mask=(
        body_mask |
        ~np.isfinite(Ugrid)
    )
)


# ------------------------------------------------------------
# Color range
# ------------------------------------------------------------

finite = plot_U.compressed()

vmin = np.percentile(
    finite,
    1
)

vmax = np.percentile(
    finite,
    99
)

levels = np.linspace(
    vmin,
    vmax,
    80
)


# ------------------------------------------------------------
# Velocity magnitude
# ------------------------------------------------------------

contour = ax.contourf(
    Xgrid,
    Ygrid,
    plot_U,
    levels=levels,
    cmap="plasma",
    extend="both",
    corner_mask=False
)


# ------------------------------------------------------------
# Color bar
# ------------------------------------------------------------

cbar = fig.colorbar(
    contour,
    ax=ax,
    pad=0.02
)

cbar.set_label(
    "Velocity magnitude [m/s]"
)


# ============================================================
# Computational mesh
# ============================================================

add_mesh(ax)


# ============================================================
# Velocity arrows
# ============================================================

step = 6

Xq = Xgrid[::step, ::step]
Yq = Ygrid[::step, ::step]

Uqx = Uxgrid[::step, ::step]
Uqy = Uygrid[::step, ::step]

mask_q = body_mask[::step, ::step]

Uqx = np.ma.array(
    Uqx,
    mask=mask_q
)

Uqy = np.ma.array(
    Uqy,
    mask=mask_q
)


ax.quiver(
    Xq,
    Yq,
    Uqx,
    Uqy,
    color="white",
    pivot="mid",
    scale=None,
    width=0.002
)


# ============================================================
# Axis formatting
# ============================================================

ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")

ax.set_title(
    "DSMC Velocity Field"
)

ax.set_aspect("equal")

fig.tight_layout()

finish_figure(fig, args, "DSMC Velocity Field")
