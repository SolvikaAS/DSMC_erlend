"""Overlay saved DSMC profiles and an optional exported SPARTA profile."""

import json
from pathlib import Path

import numpy as np

from plot_cli import save_figure


# Change these settings here; no command-line options are needed.
PROFILE_AXIS = "y"       # "x": along x at fixed y; "y": along y at fixed x.
PROFILE_POSITION = -5   # Fixed coordinate in metres.
RESULTS_DIRECTORY = Path(__file__).resolve().parents[1] / "results"
OUTPUT_DIRECTORY = RESULTS_DIRECTORY / "plots" / "comparison"
SPARTA_FILES = {
    "specular": RESULTS_DIRECTORY.parent / "sparta_y_x_-5_specular.plt",
    "diffuse": RESULTS_DIRECTORY.parent / "sparta_y_x_-5_diffuse.plt",
    "cll": RESULTS_DIRECTORY.parent / "sparta_y_x_-5_cll.plt",
}  # Set a path to None to disable its curve.


def resolve_field(path):
    path = Path(path)
    if path.is_dir():
        path = path / "DSMC2D_Field.plt"
    if not path.is_file():
        raise ValueError(f"Field file does not exist: {path}")
    return path


def read_field(field):
    data = np.loadtxt(field, skiprows=2, ndmin=2)
    if data.shape[1] != 8 or not np.isfinite(data[:, :2]).all():
        raise ValueError(f"Invalid field coordinates or columns in {field}")
    return data


def load_run(path):
    field = resolve_field(path)
    metadata_file = field.with_suffix(".json")
    if not metadata_file.is_file():
        raise ValueError(f"Missing run metadata: {metadata_file}")
    with metadata_file.open() as stream:
        metadata = json.load(stream)
    if metadata.get("obstacle", {}).get("type") != "box":
        raise ValueError(f"Expected box geometry in {metadata_file}")
    return {"path": field, "metadata": metadata, "data": read_field(field)}


def run_label(run):
    if "label" in run:
        return run["label"]
    model = run["metadata"].get("wall_model", run["path"].parent.name)
    if model == "cll":
        accommodation = run["metadata"].get("cll_accommodation")
        if accommodation is not None:
            return (f"CLL (an={accommodation['normal']:g}, "
                    f"at={accommodation['tangential']:g})")
        return "CLL"
    return model.capitalize()


def load_sparta(path, wall_model="diffuse"):
    field = resolve_field(path)
    return {
        "path": field, "metadata": {}, "data": read_field(field),
        "label": f"SPARTA ({'CLL' if wall_model == 'cll' else wall_model})", "reference": True,
    }


def check_settings(first, second):
    for key in ("obstacle", "angle_of_attack_degrees", "freestream", "grid",
                "sampling", "collision_model", "wall_temperature"):
        a = first["metadata"].get(key)
        b = second["metadata"].get(key)
        if a != b:
            raise ValueError(f"Runs have different {key} settings; use matching cases for comparison")


def extract_profile(run, axis, target):
    data = run["data"]
    varying_column, fixed_column = (0, 1) if axis == "x" else (1, 0)
    positions = np.unique(data[:, fixed_column])
    actual = positions[np.argmin(np.abs(positions - target))]
    if run.get("reference"):
        if len(np.unique(data[:, varying_column])) < 2:
            raise ValueError(f"SPARTA file has no profile along {axis}; choose the axis exported in {run['path']}")
        if len(positions) == 1 and not np.isclose(actual, target, rtol=0, atol=1e-10):
            raise ValueError(f"SPARTA profile is at {actual:g} m, not {target:g} m; update PROFILE_POSITION")
    line = data[data[:, fixed_column] == actual]
    line = line[np.argsort(line[:, varying_column])]

    obstacle = run["metadata"].get("obstacle")
    solid = np.zeros(len(line), dtype=bool)
    if obstacle is not None:
        cx, cy = obstacle["center"]
        solid = ((np.abs(line[:, 0] - cx) <= obstacle["width"] / 2)
                 & (np.abs(line[:, 1] - cy) <= obstacle["height"] / 2))
    valid = np.isfinite(line).all(axis=1) & (line[:, 2] > 0) & ~solid
    if not valid.any():
        raise ValueError(f"No sampled gas cells on the requested line in {run['path']}")
    values = {
        "density": line[:, 2].copy(),
        "temperature": line[:, 7].copy(),
        "pressure": line[:, 2] * 1.380662e-23 * line[:, 7],
        "speed": np.hypot(line[:, 3], line[:, 4]),
    }
    # Retain gaps so plotted lines never connect through the solid obstacle.
    for value in values.values():
        value[~valid] = np.nan
    return line[:, varying_column], values, actual


def main():
    paths = [RESULTS_DIRECTORY / "specular", RESULTS_DIRECTORY / "diffuse"]
    if (RESULTS_DIRECTORY / "cll" / "DSMC2D_Field.plt").is_file():
        paths.append(RESULTS_DIRECTORY / "cll")
    axis = PROFILE_AXIS
    target = PROFILE_POSITION
    if axis not in ("x", "y") or not np.isfinite(target):
        raise SystemExit("Set PROFILE_AXIS to x or y and PROFILE_POSITION to a finite coordinate.")
    try:
        runs = [load_run(path) for path in paths]
        for run in runs[1:]:
            check_settings(runs[0], run)
        for wall_model, path in SPARTA_FILES.items():
            if path is not None:
                runs.append(load_sparta(path, wall_model))
        profiles = [extract_profile(run, axis, target) for run in runs]
        if not all(np.isclose(profiles[0][2], profile[2], rtol=0, atol=1e-10)
                   for profile in profiles[1:]):
            raise ValueError("The runs select different grid lines; choose runs with matching grids")
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(str(error)) from error

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    fixed_axis = "y" if axis == "x" else "x"
    actual = profiles[0][2]
    print(f"Requested {fixed_axis} = {target:g} m; using {fixed_axis} = {actual:g} m")
    for run in runs:
        label = run_label(run)
        print(f"{label}: {run['path']}")
    quantities = {
        "temperature": ("Temperature", "K"),
        "pressure": ("Pressure", "Pa"),
        "density": ("Number density", "m⁻³"),
        "speed": ("In-plane speed", "m/s"),
    }
    for quantity, (title, unit) in quantities.items():
        fig, ax = plt.subplots(figsize=(8, 5))
        for index, (run, (coordinates, values, _)) in enumerate(zip(runs, profiles)):
            label = run_label(run)
            ax.plot(coordinates, values[quantity], label=label,
                    linestyle=("-", "--", "-.", ":")[index % 4])
        ax.set_xlabel(f"{axis} [m]")
        ax.set_ylabel(f"{title} [{unit}]")
        ax.set_title(f"{title} comparison at {fixed_axis} = {actual:.3f} m")
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        save_figure(fig, OUTPUT_DIRECTORY, f"compare_{quantity}_along_{axis}_{fixed_axis}_{actual:.3f}")


if __name__ == "__main__":
    main()
