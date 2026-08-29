"""DXF drawing and CSV spec output for generated gears."""

import os

import ezdxf
import numpy as np
import pandas as pd

from .gear import GearParams, GearGeometry, combine_tooth_dxf, rotation, transform


def _deg(rad: float) -> float:
    return rad * 360.0 / (2.0 * np.pi)


def save_dxf(
    inv_l, edge_l, root_l,
    p: GearParams,
    outer_dia: float,
    root_dia: float,
    offset_dia: float,
    work_dir: str,
) -> str:
    """Write a DXF file containing the half-tooth splines, arcs and centring lines."""
    doc = ezdxf.new("R2000")
    msp = doc.modelspace()

    x, y = combine_tooth_dxf(inv_l, edge_l, root_l)
    align_angle = np.pi / 2 + np.pi / p.z
    x3, y3 = rotation(x, y, align_angle, 1)  # align tooth to top
    x4, y4 = -x3, y3  # mirror half for the second spline
    x5, y5 = transform(x3, y3, p.x_0, p.y_0)
    x6, y6 = transform(x4, y4, p.x_0, p.y_0)

    spline5 = [(x5[i], y5[i]) for i in range(len(x5))]
    spline6 = [(x6[i], y6[i]) for i in range(len(x6))]
    msp.add_spline(spline5)
    msp.add_spline(spline6)

    outer_angle = _deg(np.arctan(x3[0] / y3[0]))
    msp.add_arc(
        center=(p.x_0, p.y_0),
        radius=outer_dia / 2,
        start_angle=90 + outer_angle,
        end_angle=90 - outer_angle,
    )

    root_start_angle = 180 / p.z
    root_end_angle = _deg(np.arctan(x3[-1] / y3[-1]))
    msp.add_arc(
        center=(p.x_0, p.y_0), radius=root_dia / 2,
        start_angle=90 - root_start_angle, end_angle=90 + root_end_angle,
    )
    msp.add_arc(
        center=(p.x_0, p.y_0), radius=root_dia / 2,
        start_angle=90 - root_end_angle, end_angle=90 + root_start_angle,
    )

    px1 = -np.sin(np.pi / p.z) * root_dia / 2 + p.x_0
    py1 = np.cos(np.pi / p.z) * root_dia / 2 + p.y_0
    px2 = -np.sin(-np.pi / p.z) * root_dia / 2 + p.x_0
    py2 = np.cos(-np.pi / p.z) * root_dia / 2 + p.y_0
    msp.add_line([p.x_0, p.y_0], [px1, py1])
    msp.add_line([p.x_0, p.y_0], [px2, py2])

    msp.add_circle((p.x_0, p.y_0), radius=offset_dia / 2)

    path = os.path.join(work_dir, "Result.dxf")
    doc.saveas(path)
    return path


def save_spec(
    p: GearParams,
    base_dia: float | None,
    pitch_dia: float,
    offset_dia: float,
    root_dia: float,
    outer_dia: float,
    work_dir: str,
) -> str:
    """Write the gear specification spreadsheet as CSV."""
    is_cycloid = p.profile == "cycloid"
    rows = []
    if is_cycloid:
        rows.append(("Type", "Cycloid"))
        rows.append(("Rolling Radius Factor (Pd/m)", p.gen_ratio))
    else:
        rows.append(("Type", "Standard" if p.alpha == 20 else "Non-Standard"))
        rows.append(("Pressure Angle", p.alpha, "deg"))
    rows += [
        ("Module", p.m, "mm"),
        ("Teeth Number", p.z, "ea"),
        ("Offset Factor", p.x),
        ("Offset", p.x * p.m, "mm"),
        ("Backlash Factor", p.b),
        ("Backlash", p.b * p.m, "mm"),
        ("Addendum Factor", p.a),
        ("Addendum", p.a * p.m, "mm"),
        ("Dedendum Factor", p.d),
        ("Dedendum", p.d * p.m, "mm"),
        ("Total Tooth Height", p.a * p.m + p.d * p.m, "mm"),
    ]
    if base_dia is not None:
        rows.append(("Base Circle Dia", base_dia, "mm"))
    rows += [
        ("Pitch Circle Dia", pitch_dia, "mm"),
        ("Offset Circle Dia", offset_dia, "mm"),
        ("Root Circle Dia", root_dia, "mm"),
        ("Outer Circle Dia", outer_dia, "mm"),
    ]
    path = os.path.join(work_dir, "Result.csv")
    pd.DataFrame(rows).to_csv(path, header=False, index=False)
    return path