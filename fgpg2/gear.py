"""Gear tooth geometry: input parameters, computed geometry, and curve construction."""

from dataclasses import dataclass, replace

import numpy as np


@dataclass
class GearParams:
    """User-controllable gear input parameters."""

    m: float
    z: int
    alpha: float
    x: float
    b: float
    a: float
    d: float
    c: float
    e: float
    x_0: float
    y_0: float
    seg_circle: int
    seg_involute: int
    seg_edge_r: int
    seg_root_r: int
    seg_outer: int
    seg_root: int
    scale: float


@dataclass
class GearGeometry:
    """Calculated intermediate values describing one tooth."""

    alpha_0: float
    alpha_m: float
    alpha_is: float
    theta_is: float
    theta_ie: float
    alpha_e: float
    x_e: float
    y_e: float
    x_e0: float
    y_e0: float
    alpha_ts: float
    theta_te: float
    e: float
    p_angle: float
    align_angle: float


def as_internal_gear(p: GearParams) -> GearParams:
    """Convert external-gear parameters to the equivalent internal-gear set."""
    if p.z >= 0:
        return p
    return replace(
        p,
        z=-p.z,
        x=-p.x,
        b=-p.b,
        a=p.d,
        d=p.a,
        c=p.e,
        e=p.c,
    )


def compute_geometry(p: GearParams) -> GearGeometry:
    """Compute the angular/positional constants of a single tooth."""
    M, Z, ALPHA, X, B, A, D, C, E = (p.m, p.z, p.alpha, p.x, p.b, p.a, p.d, p.c, p.e)

    alpha_0 = ALPHA * (2 * np.pi / 360)
    alpha_m = np.pi / Z
    alpha_is = (
        alpha_0
        + np.pi / (2 * Z)
        + B / (Z * np.cos(alpha_0))
        - (1 + 2 * X / Z) * np.sin(alpha_0) / np.cos(alpha_0)
    )
    theta_is = np.sin(alpha_0) / np.cos(alpha_0) + 2 * (
        C * (1 - np.sin(alpha_0)) + X - D
    ) / (Z * np.cos(alpha_0) * np.sin(alpha_0))
    theta_ie = (
        2 * E / (Z * np.cos(alpha_0))
        + np.sqrt(((Z + 2 * (X + A - E)) / (Z * np.cos(alpha_0))) ** 2 - 1)
    )
    alpha_e = alpha_is + theta_ie - np.arctan(
        np.sqrt(((Z + 2 * (X + A - E)) / (Z * np.cos(alpha_0))) ** 2 - 1)
    )
    x_e = M * (Z / 2 + X + A) * np.cos(alpha_e)
    y_e = M * (Z / 2 + X + A) * np.sin(alpha_e)
    x_e0 = M * (Z / 2 + X + A - E) * np.cos(alpha_e)
    y_e0 = M * (Z / 2 + X + A - E) * np.sin(alpha_e)
    alpha_ts = (
        (2 * (C * (1 - np.sin(alpha_0)) - D) * np.sin(alpha_0) + B) / (Z * np.cos(alpha_0))
        - 2 * C * np.cos(alpha_0) / Z
        + np.pi / (2 * Z)
    )
    theta_te = (
        2 * C * np.cos(alpha_0) / Z
        - 2 * (D - X - C * (1 - np.sin(alpha_0))) * np.cos(alpha_0) / (Z * np.sin(alpha_0))
    )

    # Modify E when the involute end overlaps the tooth-centre symmetry line.
    if (
        alpha_e > alpha_m
        and alpha_m > alpha_is + theta_ie - np.arctan(theta_ie)
    ):
        E = (E / 2) * np.cos(alpha_0) * (
            theta_ie - np.sqrt((1 / np.cos(alpha_is + theta_ie - alpha_m)) ** 2 - 1)
        )

    p_angle = 2 * np.pi / Z
    align_angle = np.pi / 2 - np.pi / Z

    return GearGeometry(
        alpha_0=alpha_0,
        alpha_m=alpha_m,
        alpha_is=alpha_is,
        theta_is=theta_is,
        theta_ie=theta_ie,
        alpha_e=alpha_e,
        x_e=x_e,
        y_e=y_e,
        x_e0=x_e0,
        y_e0=y_e0,
        alpha_ts=alpha_ts,
        theta_te=theta_te,
        e=E,
        p_angle=p_angle,
        align_angle=align_angle,
    )


def symmetry_y(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Reflect a curve about the X axis, reversing point order."""
    return x[::-1], -y[::-1]


def involute_curve(p: GearParams, geo: GearGeometry) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta1 = np.linspace(geo.theta_is, geo.theta_ie, p.seg_involute)
    radius = (1 / 2) * p.m * p.z * np.cos(geo.alpha_0) * np.sqrt(1 + theta1 ** 2)
    angle = geo.alpha_is + theta1 - np.arctan(theta1)
    x11 = radius * np.cos(angle)
    y11 = radius * np.sin(angle)
    return x11, y11, theta1


def edge_round_curve(
    p: GearParams,
    geo: GearGeometry,
    x_inv: np.ndarray,
    y_inv: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    E = geo.e
    theta3_min = np.arctan((y_inv[-1] - geo.y_e0) / (x_inv[-1] - geo.x_e0))
    theta3_max = np.arctan((geo.y_e - geo.y_e0) / (geo.x_e - geo.x_e0))
    theta3 = np.linspace(theta3_min, theta3_max, p.seg_edge_r)
    x21 = p.m * E * np.cos(theta3) + geo.x_e0
    y21 = p.m * E * np.sin(theta3) + geo.y_e0
    return x21, y21, theta3, theta3_min, theta3_max


def root_round_curve(
    p: GearParams, geo: GearGeometry
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    M, Z, X, D, C = p.m, p.z, p.x, p.d, p.c
    theta_t = np.linspace(0, geo.theta_te, p.seg_root_r)

    if C != 0 and (D - X - C) == 0:
        # Trochoid traces a full half-circle on the arc centreline.
        theta_s = (np.pi / 2) * np.ones(len(theta_t))
    elif (D - X - C) != 0:
        theta_s = np.arctan((M * Z * theta_t / 2) / (M * D - M * X - M * C))
    else:
        theta_s = np.zeros(len(theta_t))

    ang = theta_t + geo.alpha_ts
    x31 = M * (
        (Z / 2 + X - D + C) * np.cos(ang)
        + (Z / 2) * theta_t * np.sin(ang)
        - C * np.cos(theta_s + ang)
    )
    y31 = M * (
        (Z / 2 + X - D + C) * np.sin(ang)
        - (Z / 2) * theta_t * np.cos(ang)
        - C * np.sin(theta_s + ang)
    )
    return x31, y31, theta_t, theta_s


def outer_arc(p: GearParams, geo: GearGeometry) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta6 = np.linspace(geo.alpha_e, geo.alpha_m, p.seg_outer)
    r = p.m * (p.z / 2 + p.a + p.x)
    return r * np.cos(theta6), r * np.sin(theta6), theta6


def root_arc(p: GearParams, geo: GearGeometry) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta7 = np.linspace(0, geo.alpha_ts, p.seg_root)
    r = p.m * (p.z / 2 - p.d + p.x)
    return r * np.cos(theta7), r * np.sin(theta7), theta7


def combine_tooth(
    inv_r, inv_l, edge_r, edge_l, root_r, root_l,
    outer_r, outer_l, arc_r, arc_l,
):
    """Stitch the left/right tooth halves into a single closed outline.

    Each argument is an (x, y) array pair. The *_l segments are the
    X-axis mirror of the *_r segments; their first point is dropped to
    avoid duplicating the endpoints shared with outer/root arcs.
    """
    edge_lx, edge_ly = np.delete(edge_l[0], 0), np.delete(edge_l[1], 0)
    inv_lx, inv_ly = np.delete(inv_l[0], 0), np.delete(inv_l[1], 0)
    root_lx, root_ly = np.delete(root_l[0], 0), np.delete(root_l[1], 0)
    root_rx, root_ry = np.delete(root_r[0], 0), np.delete(root_r[1], 0)
    inv_rx, inv_ry = np.delete(inv_r[0], 0), np.delete(inv_r[1], 0)
    edge_rx, edge_ry = np.delete(edge_r[0], 0), np.delete(edge_r[1], 0)

    x = np.concatenate(
        (outer_l[0], edge_lx, inv_lx, root_lx,
         arc_l[0], arc_r[0],
         root_rx, inv_rx, edge_rx, outer_r[0])
    )
    y = np.concatenate(
        (outer_l[1], edge_ly, inv_ly, root_ly,
         arc_l[1], arc_r[1],
         root_ry, inv_ry, edge_ry, outer_r[1])
    )
    return x, y


def transform(x: np.ndarray, y: np.ndarray, dx: float, dy: float):
    return x + dx, y + dy


def rotation(x: np.ndarray, y: np.ndarray, angle: float, i: int):
    c, s = np.cos(angle * i), np.sin(angle * i)
    return c * x - s * y, s * x + c * y


def circle(dia: float, seg: int):
    theta0 = np.linspace(0.0, 2 * np.pi, seg)
    return dia / 2 * np.sin(theta0), dia / 2 * np.cos(theta0)


def combine_tooth_dxf(inv_l, edge_l, root_l):
    """Half-tooth spline used as input for the DXF export."""
    inv_x, inv_y = np.delete(inv_l[0], 0), np.delete(inv_l[1], 0)
    root_x, root_y = np.delete(root_l[0], 0), np.delete(root_l[1], 0)
    x = np.concatenate((edge_l[0], inv_x, root_x))
    y = np.concatenate((edge_l[1], inv_y, root_y))
    return x, y