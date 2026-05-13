from __future__ import annotations

import math
import os
import warnings
from dataclasses import dataclass
from pathlib import Path

import drjit as dr
import mitsuba as mi

from glass_dictionary import Eta_lookup


RAY_T_EPSILON = 1e-7

@dataclass
class LensSurface:
    def intersect(self, vertex_z, origin, direction):
        raise NotImplementedError

    @property
    def curvature_radius(self) -> float:
        return 0.0

    @property
    def is_stop(self) -> bool:
        return False


@dataclass
class SphericalSurface(LensSurface):
    radius: float

    @property
    def curvature_radius(self) -> float:
        return self.radius

    def intersect(self, vertex_z, origin, direction):
        return intersect_spherical_element(self.radius, vertex_z + self.radius, origin, direction)


@dataclass
class ApertureStop(LensSurface):
    @property
    def is_stop(self) -> bool:
        return True

    def intersect(self, vertex_z, origin, direction):
        origin = as_vector3(origin)
        direction = as_vector3(direction)
        denom_ok = direction.z != 0.0
        t = (vertex_z - origin.z) / dr.select(denom_ok, direction.z, 1.0)
        hit = denom_ok & (t >= -RAY_T_EPSILON)
        if dr.width(hit) == 1 and not bool(hit):
            return None
        return t, mi.Vector3f(0.0), hit


@dataclass
class LensElement:
    surface: LensSurface
    thickness: float
    eta: Eta_lookup
    aperture_radius: float

    @property
    def curvature_radius(self) -> float:
        return self.surface.curvature_radius

    @property
    def is_stop(self) -> bool:
        return self.surface.is_stop


@dataclass
class Bounds2:
    min_x: float = math.inf
    min_y: float = math.inf
    max_x: float = -math.inf
    max_y: float = -math.inf

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


@dataclass
class ApertureImage:
    width: int
    height: int
    pixels: list[float]

    def flat(self) -> list[float]:
        return self.pixels

    def at(self, x: int, y: int) -> float:
        return self.pixels[y * self.width + x]


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


def load_lens_elements(
    filename: str | os.PathLike,
    aperture_diameter_mm: float,
    mm_to_world: float,
) -> list[LensElement]:
    values = []
    with open(filename, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.replace(",", " ").split()
            if len(parts) != 4:
                raise ValueError(f"{filename}:{line_no}: expected 4 values per line, got {len(parts)}")
            values.extend((float(parts[0]), float(parts[1]), parse_eta(parts[2]), float(parts[3])))

    if not values or len(values) % 4:
        raise ValueError(f"{filename}: lens files must contain groups of four values")

    elements = []
    requested_aperture = aperture_diameter_mm * mm_to_world
    for i in range(0, len(values), 4):
        curvature_radius = values[i] * mm_to_world
        thickness = values[i + 1] * mm_to_world
        eta = values[i + 2]
        element_aperture = values[i + 3] * mm_to_world

        if curvature_radius == 0:
            if requested_aperture > element_aperture:
                warnings.warn(
                    f"aperture_diameter {aperture_diameter_mm:g} mm exceeds "
                    f"lens stop {element_aperture / mm_to_world:g} mm; clamping"
                )
            else:
                element_aperture = requested_aperture
            surface = ApertureStop()
        else:
            surface = SphericalSurface(curvature_radius)

        elements.append(LensElement(surface, thickness, eta, 0.5 * element_aperture))
    return elements


def parse_eta(value: str) -> Eta_lookup:
    try:
        return Eta_lookup({"type": "constant", "value": float(value)})
    except ValueError:
        return Eta_lookup({"type": "glass", "glass_name": value})


def eta_at(eta: Eta_lookup, wavelength):
    value = eta(wavelength)
    if wavelength is None:
        return scalar_value(value)
    return value


def eta_is_zero(eta: Eta_lookup) -> bool:
    return eta.type == "constant" and eta.value == 0.0


def scalar_value(value) -> float:
    if hasattr(value, "array"):
        return float(value.array[0])
    return float(value)


def is_scalar_false(mask) -> bool:
    return dr.width(mask) == 1 and not bool(mask)


def as_vector3(value):
    return value if hasattr(value, "x") else mi.Vector3f(*value)


def as_tuple3(value):
    return scalar_value(value.x), scalar_value(value.y), scalar_value(value.z)


def refract(wi, n, eta):
    wt, ok = refract_direction(wi, n, eta)
    if not bool(ok):
        return None
    return as_tuple3(wt)


def refract_direction(wi, n, eta):
    _eta = eta.array if hasattr(eta, "array") else eta
    _eta = dr.select(_eta == 0.0, 1.0, _eta)
    wi = as_vector3(wi)
    n = mi.Normal3f(as_vector3(n))
    f, cos_theta_t, _, eta_ti = mi.fresnel(dr.dot(n, wi), _eta)
    return mi.refract(wi, n, cos_theta_t, eta_ti), f < 1.0


def compute_cardinal_points(r_in_o, r_in_d, r_out_o, r_out_d):
    tf = -r_out_o[0] / r_out_d[0]
    fz = -(r_out_o[2] + tf * r_out_d[2])
    tp = (r_in_o[0] - r_out_o[0]) / r_out_d[0]
    pz = -(r_out_o[2] + tp * r_out_d[2])
    return pz, fz


def intersect_spherical_element(radius, z_center, origin, direction):
    origin = as_vector3(origin)
    direction = as_vector3(direction)
    oz_rel = origin.z - z_center
    a = dr.dot(direction, direction)
    b = 2.0 * (direction.x * origin.x + direction.y * origin.y + direction.z * oz_rel)
    c = origin.x * origin.x + origin.y * origin.y + oz_rel * oz_rel - radius * radius
    discrim = b * b - 4.0 * a * c
    root = dr.sqrt(dr.maximum(discrim, 0.0))
    t0 = (-b - root) / (2.0 * a)
    t1 = (-b + root) / (2.0 * a)
    use_closer = (direction.z > 0.0) ^ (radius < 0.0)
    t = dr.select(use_closer, dr.minimum(t0, t1), dr.maximum(t0, t1))
    hit = (discrim >= 0.0) & (t >= -RAY_T_EPSILON)
    if dr.width(hit) == 1 and not bool(hit):
        return None
    n = dr.normalize(mi.Vector3f(
        origin.x + t * direction.x,
        origin.y + t * direction.y,
        oz_rel + t * direction.z,
    ))
    n = dr.select(dr.dot(n, -direction) < 0.0, -n, n)
    return t, n, hit


def aperture_lookup(image: ApertureImage | None, u: float, v: float) -> float:
    if image is None or u < 0.0 or u > 1.0 or v < 0.0 or v > 1.0:
        return 0.0
    x = u * image.width - 0.5
    y = v * image.height - 0.5
    x0 = math.floor(x)
    y0 = math.floor(y)
    tx = x - x0
    ty = y - y0

    def sample(ix, iy):
        return image.at(ix, iy) if 0 <= ix < image.width and 0 <= iy < image.height else 0.0

    return (
        (1 - tx) * (1 - ty) * sample(x0, y0)
        + tx * (1 - ty) * sample(x0 + 1, y0)
        + (1 - tx) * ty * sample(x0, y0 + 1)
        + tx * ty * sample(x0 + 1, y0 + 1)
    )


def aperture_lookup_vec(image: ApertureImage, flat_image, u, v):
    h, w = image.height, image.width
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


def make_aperture_image(name: str, base_dir: Path, res: int = 256):
    if not name:
        return None
    width = height = res
    if name == "gaussian":
        pixels = []
        for y in range(res):
            for x in range(res):
                uvx = -1.0 + 2.0 * (x + 0.5) / res
                uvy = -1.0 + 2.0 * (y + 0.5) / res
                pixels.append(max(0.0, math.exp(-(uvx * uvx + uvy * uvy)) - math.exp(-1.0)))
    elif name == "square":
        pixels = []
        for y in range(res):
            for x in range(res):
                pixels.append(4.0 if res // 4 <= x < 3 * res // 4 and res // 4 <= y < 3 * res // 4 else 0.0)
    elif name in ("pentagon", "star"):
        pixels = rasterize_aperture(name, res)
    else:
        bitmap = mi.Bitmap(str(resolve_file(name, str(base_dir))))
        bitmap = bitmap.convert(mi.Bitmap.PixelFormat.Y, mi.Struct.Type.Float32, False)
        width, height = bitmap.width(), bitmap.height()
        pixels = [float(v) for v in mi.TensorXf(bitmap).array]

    pixels = flip_y(pixels, width, height)
    avg = sum(pixels) / len(pixels)
    if avg > 0.0:
        scale = (math.pi / 4.0) / avg
        pixels = [p * scale for p in pixels]
    return ApertureImage(width, height, pixels)


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

    pixels = []
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
            pixels.append(1.0 if winding else 0.0)
    return pixels


def flip_y(pixels: list[float], width: int, height: int) -> list[float]:
    flipped = []
    for y in range(height):
        row = height - 1 - y
        flipped.extend(pixels[row * width : (row + 1) * width])
    return flipped
