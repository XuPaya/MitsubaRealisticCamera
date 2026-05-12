from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path

import drjit as dr
import mitsuba as mi
import numpy as np

from glass_dictionary import Eta_lookup


RAY_T_EPSILON = 1e-7

@dataclass
class LensElementLegacy:
    curvature_radius: float
    thickness: float
    eta: float
    aperture_radius: float

@dataclass
class LensElement:
    curvature_radius: float
    thickness: float
    eta: Eta_lookup
    aperture_radius: float


@dataclass
class Bounds2:
    min_x: float = math.inf
    min_y: float = math.inf
    max_x: float = -math.inf
    max_y: float = -math.inf

    @property
    def area(self) -> float:
        return max(0.0, self.max_x - self.min_x) * max(0.0, self.max_y - self.min_y)

    def is_degenerate(self) -> bool:
        return self.max_x <= self.min_x or self.max_y <= self.min_y

    def expand(self, delta: float) -> "Bounds2":
        if self.is_degenerate():
            return self
        return Bounds2(
            self.min_x - delta,
            self.min_y - delta,
            self.max_x + delta,
            self.max_y + delta,
        )


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def normalize(v):
    inv_len = 1.0 / math.sqrt(dot(v, v))
    return (v[0] * inv_len, v[1] * inv_len, v[2] * inv_len)


def get_property(props, names, default):
    for name in names:
        if name in props:
            return props.get(name)
    return default


def resolve_file(filename: str, directory: str = "") -> Path:
    path = Path(filename)
    candidates = []
    if directory:
        candidates.append(Path(directory) / path)
    candidates.append(path)

    try:
        candidates.append(Path(mi.file_resolver().resolve(str(path))))
    except Exception:
        pass

    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    if not path.is_absolute():
        for scene_dir in Path.cwd().glob("scenes/*"):
            candidate = scene_dir / path
            if candidate.is_file():
                return candidate.resolve()
    return candidates[0].resolve()


def load_lens_file(filename: str | os.PathLike) -> list[float]:
    values: list[float] = []
    with open(filename, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.replace(",", " ").split()
            try:
                values.extend(float(p) for p in parts)
            except ValueError as exc:
                raise ValueError(f"{filename}:{line_no}: invalid lens value") from exc
    if not values or len(values) % 4:
        raise ValueError(f"{filename}: lens files must contain groups of four floats")
    return values


def intersect_spherical_element(radius, z_center, origin, direction):
    ox, oy, oz = origin[0], origin[1], origin[2] - z_center
    dx, dy, dz = direction
    a = dot(direction, direction)
    b = 2.0 * (dx * ox + dy * oy + dz * oz)
    c = ox * ox + oy * oy + oz * oz - radius * radius
    discrim = b * b - 4.0 * a * c
    if discrim < 0:
        print("no intersection with spherical element: discriminant =", discrim)
        return None
    root = math.sqrt(discrim)
    t0 = (-b - root) / (2.0 * a)
    t1 = (-b + root) / (2.0 * a)
    use_closer_t = (dz > 0.0) ^ (radius < 0.0)
    t = min(t0, t1) if use_closer_t else max(t0, t1)
    if t < -RAY_T_EPSILON:
        print("intersection with spherical element is behind the ray origin: t =", t)
        return None
    n = normalize((ox + t * dx, oy + t * dy, oz + t * dz))
    if dot(n, (-dx, -dy, -dz)) < 0:
        n = (-n[0], -n[1], -n[2])
    return t, n


def refract(wi, n, eta):
    _eta = eta
    if type(eta) is not float:
        _eta = eta.numpy()
    if _eta == 0.0:
        _eta = 1.0
    cos_theta_i = dot(n, wi)
    if cos_theta_i < 0:
        _eta = 1.0 / _eta
        cos_theta_i = -cos_theta_i
        n = (-n[0], -n[1], -n[2])
    sin2_theta_i = max(0.0, 1.0 - cos_theta_i * cos_theta_i)
    sin2_theta_t = sin2_theta_i / (_eta * _eta)
    if sin2_theta_t >= 1.0:
        return None
    cos_theta_t = math.sqrt(1.0 - sin2_theta_t)
    return (
        -wi[0] / _eta + (cos_theta_i / _eta - cos_theta_t) * n[0],
        -wi[1] / _eta + (cos_theta_i / _eta - cos_theta_t) * n[1],
        -wi[2] / _eta + (cos_theta_i / _eta - cos_theta_t) * n[2],
    )


def compute_cardinal_points(r_in_o, r_in_d, r_out_o, r_out_d):
    tf = -r_out_o[0] / r_out_d[0]
    fz = -(r_out_o[2] + tf * r_out_d[2])
    tp = (r_in_o[0] - r_out_o[0]) / r_out_d[0]
    pz = -(r_out_o[2] + tp * r_out_d[2])
    return pz, fz


def intersect_spherical_element_vec(radius, z_center, ox, oy, oz, dx, dy, dz):
    oz_rel = oz - z_center
    a = dx * dx + dy * dy + dz * dz
    b = 2.0 * (dx * ox + dy * oy + dz * oz_rel)
    c = ox * ox + oy * oy + oz_rel * oz_rel - radius * radius
    discrim = b * b - 4.0 * a * c
    root = dr.sqrt(dr.maximum(discrim, 0.0))
    t0 = (-b - root) / (2.0 * a)
    t1 = (-b + root) / (2.0 * a)
    use_closer = (dz > 0.0) ^ (radius < 0.0)
    t = dr.select(use_closer, dr.minimum(t0, t1), dr.maximum(t0, t1))
    hit = (discrim >= 0.0) & (t >= 0.0)
    n = dr.normalize(mi.Vector3f(ox + t * dx, oy + t * dy, oz_rel + t * dz))
    n = dr.select(dr.dot(n, mi.Vector3f(-dx, -dy, -dz)) < 0.0, -n, n)
    return t, n.x, n.y, n.z, hit


def refract_vec(wi, n, eta):
    cos_i = dr.dot(n, wi)
    flip = cos_i < 0.0
    _eta = eta.array
    # if eta is 0, set it to 1;
    # _eta = dr.select(_eta == 0.0, 1.0, _eta)
    eta_eff = dr.select(flip, 1.0 / _eta, _eta)
    cos_i = dr.select(flip, -cos_i, cos_i)
    n = dr.select(flip, -n, n)
    sin2_i = dr.maximum(0.0, 1.0 - cos_i * cos_i)
    sin2_t = sin2_i / (eta_eff * eta_eff)
    ok = sin2_t < 1.0
    cos_t = dr.sqrt(dr.maximum(0.0, 1.0 - sin2_t))
    return -wi / eta_eff + (cos_i / eta_eff - cos_t) * n, ok


def intersect_spherical_element_np(radius, z_center, ox, oy, oz, dx, dy, dz):
    oz_rel = oz - z_center
    a = dx * dx + dy * dy + dz * dz
    b = 2.0 * (dx * ox + dy * oy + dz * oz_rel)
    c = ox * ox + oy * oy + oz_rel * oz_rel - radius * radius
    discrim = b * b - 4.0 * a * c
    root = np.sqrt(np.maximum(discrim, 0.0))
    t0 = (-b - root) / (2.0 * a)
    t1 = (-b + root) / (2.0 * a)
    use_closer = (dz > 0.0) ^ (radius < 0.0)
    t = np.where(use_closer, np.minimum(t0, t1), np.maximum(t0, t1))
    hit = (discrim >= 0.0) & (t >= 0.0)
    nx, ny, nz = ox + t * dx, oy + t * dy, oz_rel + t * dz
    inv_len = 1.0 / np.sqrt(nx * nx + ny * ny + nz * nz)
    nx, ny, nz = nx * inv_len, ny * inv_len, nz * inv_len
    flip = nx * -dx + ny * -dy + nz * -dz < 0.0
    return t, (
        np.where(flip, -nx, nx),
        np.where(flip, -ny, ny),
        np.where(flip, -nz, nz),
    ), hit


def refract_np(wi, n, eta):
    if eta == 0.0:
        eta = 1.0
    cos_i = n[0] * wi[0] + n[1] * wi[1] + n[2] * wi[2]
    flip = cos_i < 0.0
    eta_eff = np.where(flip, 1.0 / eta, eta)
    cos_i = np.where(flip, -cos_i, cos_i)
    nx = np.where(flip, -n[0], n[0])
    ny = np.where(flip, -n[1], n[1])
    nz = np.where(flip, -n[2], n[2])
    sin2_i = np.maximum(0.0, 1.0 - cos_i * cos_i)
    sin2_t = sin2_i / (eta_eff * eta_eff)
    ok = sin2_t < 1.0
    cos_t = np.sqrt(np.maximum(0.0, 1.0 - sin2_t))
    scale = cos_i / eta_eff - cos_t
    return (
        -wi[0] / eta_eff + scale * nx,
        -wi[1] / eta_eff + scale * ny,
        -wi[2] / eta_eff + scale * nz,
    ), ok


def aperture_lookup_np(image, u, v):
    if image is None:
        return np.zeros_like(u)
    h, w = image.shape
    x = u * w - 0.5
    y = v * h - 0.5
    x0 = np.floor(x).astype(np.int64)
    y0 = np.floor(y).astype(np.int64)
    tx = x - x0
    ty = y - y0

    def sample(ix, iy):
        ok = (ix >= 0) & (ix < w) & (iy >= 0) & (iy < h)
        return np.where(ok, image[np.clip(iy, 0, h - 1), np.clip(ix, 0, w - 1)], 0.0)

    return (
        (1 - tx) * (1 - ty) * sample(x0, y0)
        + tx * (1 - ty) * sample(x0 + 1, y0)
        + (1 - tx) * ty * sample(x0, y0 + 1)
        + tx * ty * sample(x0 + 1, y0 + 1)
    )


def aperture_lookup_vec(image, flat_image, u, v):
    h, w = image.shape
    x = u * w - 0.5
    y = v * h - 0.5
    x0 = mi.Int32(dr.floor(x))
    y0 = mi.Int32(dr.floor(y))
    tx = x - mi.Float(x0)
    ty = y - mi.Float(y0)

    def sample(ix, iy):
        ok = (ix >= 0) & (ix < w) & (iy >= 0) & (iy < h)
        idx = mi.UInt32(
            dr.maximum(0, dr.minimum(w - 1, ix))
            + dr.maximum(0, dr.minimum(h - 1, iy)) * w
        )
        return dr.select(ok, dr.gather(mi.Float, flat_image, idx), 0.0)

    return (
        (1 - tx) * (1 - ty) * sample(x0, y0)
        + tx * (1 - ty) * sample(x0 + 1, y0)
        + (1 - tx) * ty * sample(x0, y0 + 1)
        + tx * ty * sample(x0 + 1, y0 + 1)
    )


def radical_inverse(base: int, values: np.ndarray) -> np.ndarray:
    n = values.copy()
    inv_base = 1.0 / base
    inv = inv_base
    reversed_digits = np.zeros(n.shape, dtype=np.float64)
    while np.any(n):
        next_n = n // base
        digit = n - next_n * base
        reversed_digits += digit.astype(np.float64) * inv
        inv *= inv_base
        n = next_n
    return reversed_digits


def make_aperture_image(name: str, base_dir: Path, res: int = 256):
    if not name:
        return None
    if name == "gaussian":
        y, x = np.mgrid[0:res, 0:res]
        uvx = -1.0 + 2.0 * (x + 0.5) / res
        uvy = -1.0 + 2.0 * (y + 0.5) / res
        image = np.maximum(0.0, np.exp(-(uvx * uvx + uvy * uvy)) - math.exp(-1.0))
    elif name == "square":
        image = np.zeros((res, res), dtype=np.float32)
        image[res // 4 : 3 * res // 4, res // 4 : 3 * res // 4] = 4.0
    elif name in ("pentagon", "star"):
        image = rasterize_aperture(name, res)
    else:
        bitmap = mi.Bitmap(str(resolve_file(name, str(base_dir))))
        image = np.array(bitmap, dtype=np.float32)
        if image.ndim == 3:
            image = np.mean(image[..., :3], axis=2)
    image = np.flipud(image.astype(np.float32))
    avg = float(np.mean(image))
    return image * ((math.pi / 4.0) / avg) if avg > 0.0 else image


def rasterize_aperture(name: str, res: int):
    if name == "pentagon":
        c1 = (math.sqrt(5.0) - 1.0) / 4.0
        c2 = (math.sqrt(5.0) + 1.0) / 4.0
        s1 = math.sqrt(10.0 + 2.0 * math.sqrt(5.0)) / 4.0
        s2 = math.sqrt(10.0 - 2.0 * math.sqrt(5.0)) / 4.0
        vertices = [
            (0.0, 0.8),
            (0.8 * s1, 0.8 * c1),
            (0.8 * s2, -0.8 * c2),
            (-0.8 * s2, -0.8 * c2),
            (-0.8 * s1, 0.8 * c1),
        ]
    else:
        vertices = []
        for i in reversed(range(10)):
            radius = 1.0 if i & 1 else math.cos(math.radians(72.0)) / math.cos(math.radians(36.0))
            vertices.append((radius * math.cos(math.pi * i / 5.0), radius * math.sin(math.pi * i / 5.0)))

    image = np.zeros((res, res), dtype=np.float32)
    for y in range(res):
        for x in range(res):
            px = -1.0 + 2.0 * (x + 0.5) / res
            py = -1.0 + 2.0 * (y + 0.5) / res
            winding = 0
            for i, v0 in enumerate(vertices):
                v1 = vertices[(i + 1) % len(vertices)]
                edge = (px - v0[0]) * (v1[1] - v0[1]) - (py - v0[1]) * (v1[0] - v0[0])
                if v0[1] <= py:
                    winding += 1 if v1[1] > py and edge > 0 else 0
                elif v1[1] <= py and edge < 0:
                    winding -= 1
            image[y, x] = 1.0 if winding else 0.0
    return image
