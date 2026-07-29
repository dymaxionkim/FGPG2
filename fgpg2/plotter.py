"""Render gear geometry with matplotlib and persist preview images + DXF + CSV."""

import os

import matplotlib.pyplot as plt
import numpy as np

from .exporters import save_dxf, save_spec
from .gear import (
    GearParams,
    as_internal_gear,
    circle,
    combine_tooth,
    compute_geometry,
    edge_round_curve,
    involute_curve,
    outer_arc,
    root_arc,
    root_round_curve,
    rotation,
    symmetry_y,
    transform,
)


def generate(p: GearParams, work_dir: str) -> None:
    """Plot the whole gear and one tooth, saving Result1/2.png, Result.dxf, Result.csv."""
    p = as_internal_gear(p)
    geo = compute_geometry(p)

    inv_r = involute_curve(p, geo)
    inv_l = symmetry_y(inv_r[0], inv_r[1])
    edge_r = edge_round_curve(p, geo, *inv_r[:2])
    edge_l = symmetry_y(edge_r[0], edge_r[1])
    root_r = root_round_curve(p, geo)
    root_l = symmetry_y(root_r[0], root_r[1])
    outer_r = outer_arc(p, geo)
    outer_l = symmetry_y(outer_r[0], outer_r[1])
    arc_r = root_arc(p, geo)
    arc_l = symmetry_y(arc_r[0], arc_r[1])

    x1, y1 = combine_tooth(
        inv_r[:2], inv_l, edge_r[:2], edge_l, root_r[:2], root_l,
        outer_r[:2], outer_l, arc_r[:2], arc_l,
    )

    x2, y2 = rotation(x1, y1, geo.align_angle, 1)

    fig = plt.figure(figsize=(5, 5))
    ax = plt.axes()
    ax.set_aspect("equal")
    ax.set_title("Fine Gear Profile Generator 3")
    ax.grid(True)

    for i in range(p.z):
        xt, yt = rotation(x2, y2, geo.p_angle, i)
        xt, yt = transform(xt, yt, p.x_0, p.y_0)
        ax.plot(xt, yt, "-", linewidth=1.5, color="black")

    base_dia = p.m * p.z * np.cos(geo.alpha_0)
    pitch_dia = p.m * p.z
    offset_dia = 2 * p.m * (p.z / 2 + p.x)
    outer_dia = 2 * p.m * (p.z / 2 + p.x + p.a)
    root_dia = 2 * p.m * (p.z / 2 + p.x - p.d)

    for dia, style, color, label in (
        (base_dia, ":", "cyan", "Base Circle"),
        (pitch_dia, ":", "magenta", "Pitch Circle"),
        (offset_dia, ":", "red", "Offset Circle"),
        (outer_dia, ":", "brown", "Outer Circle"),
        (root_dia, ":", "grey", "Root Circle"),
    ):
        xc, yc = transform(*circle(dia, p.seg_circle), p.x_0, p.y_0)
        ax.plot(xc, yc, style, linewidth=1.0, color=color, label=label)

    _draw_annotations(ax, p, base_dia, pitch_dia, offset_dia, outer_dia, root_dia)

    whole_path = os.path.join(work_dir, "Result1.png")
    fig.savefig(whole_path, dpi=100)

    y_height = ((outer_dia - root_dia) / 2.0) / p.scale
    y_cen = (outer_dia + root_dia) / 4.0 + p.y_0
    ax.set_xlim(-(y_height / 2.0) + p.x_0, (y_height / 2.0) + p.x_0)
    ax.set_ylim(y_cen - y_height / 2.0, y_cen + y_height / 2.0)
    tooth_path = os.path.join(work_dir, "Result2.png")
    fig.savefig(tooth_path, dpi=100)
    plt.close(fig)

    save_dxf(inv_l, edge_l, root_l, p, outer_dia, root_dia, offset_dia, work_dir)
    save_spec(p, base_dia, pitch_dia, offset_dia, root_dia, outer_dia, work_dir)


def _draw_annotations(ax, p, base_dia, pitch_dia, offset_dia, outer_dia, root_dia):
    height = root_dia / 22.5
    nrow = 15 * height
    green = []
    colored = []

    for i, (fmt, val) in enumerate(
        (
            ("Module m=%s[mm]", p.m),
            ("Teeth Number z=%s[ea]", p.z),
            ("Pressure angle alpha=%s[deg]", p.alpha),
            ("Offset factor x=%s", p.x),
            ("Backlash factor b=%s", p.b),
            ("Addendum factor a=%s", p.a),
            ("Dedendum factor d=%s", p.d),
            ("Radius Factor of Edge Round of Hob c=%s", p.c),
            ("Radius Factor of Edge Round of Tooth e=%s", p.e),
            ("Center Position = %s[mm]", [p.x_0, p.y_0]),
        )
    ):
        green.append((fmt % (val,), "green"))
    for i, (fmt, val, color) in enumerate(
        (
            ("Base Circle Dia = %s[mm]", base_dia, "cyan"),
            ("Pitch Circle Dia = %s[mm]", pitch_dia, "magenta"),
            ("Offset Circle Dia = %s[mm]", offset_dia, "red"),
            ("Outer Circle Dia = %s[mm]", outer_dia, "brown"),
            ("Root Circle Dia = %s[mm]", root_dia, "grey"),
        ),
        start=10,
    ):
        colored.append((fmt % (val,), color))

    for i, (text, color) in enumerate(green + colored):
        ax.text(
            p.x_0, p.y_0 + nrow / 2 - height * i,
            text,
            verticalalignment="center", horizontalalignment="center",
            color=color, fontsize="x-small",
        )