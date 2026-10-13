"""Shared command-line options and figure output for DSMC plotting scripts."""

import argparse
import json
from pathlib import Path
import re

# Select the model to plot here when using the single-field plotting scripts.
WALL_MODEL = "cll"  # "specular", "diffuse", or "cll"


def read_run_metadata(filename):
    """Read geometry saved beside the field; use the default box for older files."""
    metadata_file = filename.with_suffix(".json")
    if metadata_file.is_file():
        with metadata_file.open() as stream:
            return json.load(stream)
    return {
        "obstacle": {"type": "box", "center": [0.0, 0.0], "width": 6.0, "height": 6.0},
    }


def parse_args(description, *, x_slice=False):
    results_dir = Path(__file__).resolve().parents[1] / "results"
    default_file = results_dir / WALL_MODEL / "DSMC2D_Field.plt"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "filename", nargs="?", type=Path, default=default_file,
        help=f"DSMC field file (default: project results/{WALL_MODEL}/DSMC2D_Field.plt)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=results_dir / "plots",
        help="Save PNG plots in this directory (default: project results/plots)",
    )
    if x_slice:
        parser.add_argument("--x", type=float, default=-5.0, help="Requested x position in metres")
    args = parser.parse_args()
    if not args.filename.is_file():
        parser.error(f"field file does not exist: {args.filename}")
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        # Scripts parse options before importing pyplot. File output should
        # not attempt to connect to an X11/Wayland desktop display.
        import matplotlib
        matplotlib.use("Agg")
    return args


def finish_figure(fig, args, name):
    import matplotlib.pyplot as plt

    if args.output_dir is None:
        plt.show()
    else:
        save_figure(fig, args.output_dir, name)


def save_figure(fig, output_dir, name):
    import matplotlib.pyplot as plt

    filename = re.sub(r"[^a-z0-9-]+", "_", name.lower()).strip("_") + ".png"
    destination = output_dir / filename
    fig.savefig(destination, dpi=150)
    print(f"Saved {destination}")
    plt.close(fig)
