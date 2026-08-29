"""Command-line gear generator for FGPG2.

Usage:
    python FGPG2_CLI.py <path/to/Results.csv>

Reads the gear parameters from the given CSV (two columns: parameter/value,
as written by the GUI's Inputs.csv) and writes Result.csv, Result.dxf,
Result1.png and Result2.png into the same directory as the input file.
"""

import argparse
import os
import sys

import pandas as pd

from fgpg2.gear import GearParams
from fgpg2.plotter import generate

DEFAULT_PARAMS = GearParams(
    m=1.0, z=18, alpha=20, x=0.0, b=0.05, a=1.0, d=1.25, c=0.2, e=0.1,
    x_0=0.0, y_0=0.0, seg_circle=360, seg_involute=15, seg_edge_r=5,
    seg_root_r=5, seg_outer=5, seg_root=5, scale=0.7,
)

_INT_FIELDS = {
    "z", "seg_circle", "seg_involute", "seg_edge_r",
    "seg_root_r", "seg_outer", "seg_root",
}


def load_params(csv_path: str) -> GearParams:
    """Read the parameter/value CSV into GearParams, defaulting any missing key."""
    df = pd.read_csv(csv_path)
    if df.shape[1] < 2:
        raise ValueError(f"Expected a parameter/value CSV, got {csv_path}")
    values = df.set_index(df.columns[0])[df.columns[1]]

    kwargs = {}
    for key in DEFAULT_PARAMS.__dataclass_fields__:
        if key in values.index:
            raw = values.loc[key]
            kwargs[key] = int(raw) if key in _INT_FIELDS else float(raw)
        else:
            kwargs[key] = getattr(DEFAULT_PARAMS, key)
    return GearParams(**kwargs)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate gear result files from a parameters CSV.",
    )
    parser.add_argument(
        "csv",
        help="Path to the input parameters CSV (e.g. Results.csv).",
    )
    args = parser.parse_args()

    csv_path = os.path.abspath(args.csv)
    if not os.path.exists(csv_path):
        print(f"Error: file not found: {csv_path}", file=sys.stderr)
        return 1

    try:
        p = load_params(csv_path)
    except Exception as exc:
        print(f"Error: failed to read parameters: {exc}", file=sys.stderr)
        return 1

    work_dir = os.path.dirname(csv_path)
    os.makedirs(work_dir, exist_ok=True)
    generate(p, work_dir)

    for name in ("Result.csv", "Result.dxf", "Result1.png", "Result2.png"):
        path = os.path.join(work_dir, name)
        status = "OK" if os.path.exists(path) else "MISSING"
        print(f"[{status}] {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())