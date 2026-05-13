from __future__ import annotations

import math

import drjit as dr
import mitsuba as mi

from utils import *


class RealisticCamera(mi.Sensor):
    def __init__(self, props):
        super().__init__(props)
        self.m_needs_sample_2 = True
        self.m_needs_sample_3 = True

        lens_file = str(get_property(props, ("lens_file", "lensfile"), ""))
        lens_directory = str(get_property(props, ("lens_directory",), ""))
        if not lens_file:
            raise ValueError("realisticcamera requires a lens_file string")

        self.lens_file = resolve_file(lens_file, lens_directory)
        self.aperture_diameter = float(
            get_property(props, ("aperture_diameter", "aperturediameter"), 1.0)
        )
        self.focus_distance = float(get_property(props, ("focus_distance", "focusdistance"), 10.0))
        self.film_diagonal = float(get_property(props, ("film_diagonal",), 35.0))
        self.mm_to_world = float(get_property(props, ("mm_to_world",), 0.001))
        self.use_exit_pupil = bool(get_property(props, ("use_exit_pupil",), True))
        self.debug_lens_trace = bool(get_property(props, ("debug_lens_trace",), False))
        self.exit_pupil_sample_count = int(get_property(props, ("exit_pupil_sample_count",), 65536))
        self.aperture = str(get_property(props, ("aperture", "aperture_image"), ""))

        self.aperture_image = make_aperture_image(self.aperture, self.lens_file.parent)
        self._radical_inverse = mi.RadicalInverse(3)
        self.elements = load_lens_elements(
            self.lens_file,
            self.aperture_diameter,
            self.mm_to_world,
        )
        self._aperture_flat = (
            mi.Float(self.aperture_image.flat())
            if dr.is_jit_v(mi.Float) and self.aperture_image is not None
            else None
        )

        film_size = self.film().size()
        self.full_resolution = (int(film_size.x), int(film_size.y))
        crop_size = self.film().crop_size()
        self.crop_size = (int(crop_size.x), int(crop_size.y))
        crop_offset = self.film().crop_offset()
        self.crop_offset = (int(crop_offset.x), int(crop_offset.y))
        aspect = self.full_resolution[1] / self.full_resolution[0]
        diagonal = self.film_diagonal * self.mm_to_world
        x = math.sqrt(diagonal * diagonal / (1.0 + aspect * aspect))
        y = aspect * x
        self.physical_extent = Bounds2(-0.5 * x, -0.5 * y, 0.5 * x, 0.5 * y)
        self._film_diag_world = diagonal

        self.elements[-1].thickness = self.focus_thick_lens(self.focus_distance)
        self.exit_pupil_bounds = self._compute_exit_pupil_bounds()
        self._prepare_jit_tables()

        if self.debug_lens_trace:
            print(self.to_string())

    def _prepare_jit_tables(self):
        if not dr.is_jit_v(mi.Float):
            self._exit_min_x = self._exit_min_y = None
            self._exit_max_x = self._exit_max_y = None
            self._aperture_flat = None
            return
        self._exit_min_x = mi.Float([b.min_x for b in self.exit_pupil_bounds])
        self._exit_min_y = mi.Float([b.min_y for b in self.exit_pupil_bounds])
        self._exit_max_x = mi.Float([b.max_x for b in self.exit_pupil_bounds])
        self._exit_max_y = mi.Float([b.max_y for b in self.exit_pupil_bounds])
        self._aperture_flat = (
            mi.Float(self.aperture_image.flat())
            if self.aperture_image is not None
            else None
        )

    def _radical_inverse_samples(self, start, end, n_samples):
        if dr.is_jit_v(mi.Float):
            idx = dr.arange(mi.UInt32, start, end)
        else:
            idx = mi.UInt32(start)
        return (
            self._radical_inverse.eval(0, idx),
            self._radical_inverse.eval(1, idx),
            (mi.Float(idx) + 0.5) / n_samples,
        )

    def _exit_bound_values(self, r_film):
        n = len(self.exit_pupil_bounds)
        if self._exit_min_x is not None:
            idx = mi.UInt32(dr.minimum(mi.Float(n - 1), r_film / (self._film_diag_world * 0.5) * n))
            return (
                dr.gather(mi.Float, self._exit_min_x, idx),
                dr.gather(mi.Float, self._exit_min_y, idx),
                dr.gather(mi.Float, self._exit_max_x, idx),
                dr.gather(mi.Float, self._exit_max_y, idx),
                mi.Mask(True),
            )

        r_index = min(n - 1, int(scalar_value(r_film) / (self._film_diag_world * 0.5) * n))
        bound = self.exit_pupil_bounds[r_index]
        return bound.min_x, bound.min_y, bound.max_x, bound.max_y, not bound.is_degenerate()

    def _merge_exit_bound_sample(self, bound, active, rx, ry):
        if dr.width(active) > 1:
            if bool(dr.any(active)):
                bound.min_x = min(bound.min_x, float(dr.min(dr.select(active, rx, math.inf))[0]))
                bound.min_y = min(bound.min_y, float(dr.min(dr.select(active, ry, math.inf))[0]))
                bound.max_x = max(bound.max_x, float(dr.max(dr.select(active, rx, -math.inf))[0]))
                bound.max_y = max(bound.max_y, float(dr.max(dr.select(active, ry, -math.inf))[0]))
        elif bool(active):
            bound.min_x = min(bound.min_x, scalar_value(rx))
            bound.min_y = min(bound.min_y, scalar_value(ry))
            bound.max_x = max(bound.max_x, scalar_value(rx))
            bound.max_y = max(bound.max_y, scalar_value(ry))

    def lens_rear_z(self) -> float:
        return self.elements[-1].thickness

    def lens_front_z(self) -> float:
        return sum(e.thickness for e in self.elements)

    def rear_element_radius(self) -> float:
        return self.elements[-1].aperture_radius

    def trace_lenses_from_film(self, origin, direction):
        return self._trace_lenses(
            origin[0], origin[1], origin[2],
            direction[0], direction[1], direction[2],
            from_film=True,
        )[:3]

    def trace_lenses_from_scene(self, origin, direction):
        return self._trace_lenses(
            origin[0], origin[1], origin[2],
            direction[0], direction[1], direction[2],
            from_film=False,
        )[:3]

    def _element_eta(self, index, wavelength):
        if index < 0 or eta_is_zero(self.elements[index].eta):
            return 1.0
        return eta_at(self.elements[index].eta, wavelength)

    def _eta_pair(self, index, wavelength, from_film):
        if from_film:
            return self._element_eta(index, wavelength), self._element_eta(index - 1, wavelength)
        return self._element_eta(index - 1, wavelength), self._element_eta(index, wavelength)

    def _trace_lenses(self, ox, oy, oz, dx, dy, dz, from_film, wavelength=None, active=True):
        origin = mi.Vector3f(ox, oy, -oz)
        direction = mi.Vector3f(dx, dy, -dz)
        element_z = 0.0 if from_film else -self.lens_front_z()
        weight = mi.Float(1.0)
        active = mi.Mask(active)
        indices = range(len(self.elements) - 1, -1, -1) if from_film else range(len(self.elements))
        if is_scalar_false(active):
            return 0.0, None, None, False

        for i in indices:
            element = self.elements[i]
            if from_film:
                element_z -= element.thickness
            hit_info = element.surface.intersect(element_z, origin, direction)
            if hit_info is None:
                return 0.0, None, None, False
            t, normal, hit = hit_info

            hit_point = origin + t * direction
            if element.is_stop and from_film and self.aperture_image is not None:
                u = (hit_point.x / element.aperture_radius + 1.0) * 0.5
                v = (hit_point.y / element.aperture_radius + 1.0) * 0.5
                if self._aperture_flat is not None:
                    weight = aperture_lookup_vec(self.aperture_image, self._aperture_flat, u, v)
                else:
                    weight = aperture_lookup(self.aperture_image, u, v)
                aperture_ok = weight != 0.0
            else:
                aperture_ok = (
                    hit_point.x * hit_point.x + hit_point.y * hit_point.y
                    <= element.aperture_radius * element.aperture_radius
                )
            active &= hit & aperture_ok
            if is_scalar_false(active):
                return 0.0, None, None, False

            origin = hit_point
            if not element.is_stop:
                eta_i, eta_t = self._eta_pair(i, wavelength, from_film)
                direction, refr_ok = refract_direction(dr.normalize(-direction), normal, eta_t / eta_i)
                active &= refr_ok
                if is_scalar_false(active):
                    return 0.0, None, None, False
            if not from_film:
                element_z += element.thickness

        return (
            weight,
            mi.Point3f(origin.x, origin.y, -origin.z),
            mi.Vector3f(direction.x, direction.y, -direction.z),
            active,
        )

    def compute_thick_lens_approximation(self):
        x = 0.001 * self._film_diag_world
        r_scene_o = (x, 0.0, self.lens_front_z() + 1.0)
        r_scene_d = (0.0, 0.0, -1.0)
        _, r_film_o, r_film_d = self.trace_lenses_from_scene(r_scene_o, r_scene_d)
        if r_film_o is None:
            raise RuntimeError("unable to trace scene ray for thick lens approximation: error code %.1f" % _)
        pz0, fz0 = compute_cardinal_points(r_scene_o, r_scene_d, r_film_o, r_film_d)

        r_film_o = (x, 0.0, self.lens_rear_z() - 1.0)
        r_film_d = (0.0, 0.0, 1.0)
        _, r_scene_o, r_scene_d = self.trace_lenses_from_film(r_film_o, r_film_d)
        if r_scene_o is None:
            raise RuntimeError("unable to trace film ray for thick lens approximation")
        pz1, fz1 = compute_cardinal_points(r_film_o, r_film_d, r_scene_o, r_scene_d)
        return (pz0, pz1), (fz0, fz1)

    def focus_thick_lens(self, focus_distance):
        pz, fz = self.compute_thick_lens_approximation()
        pz = tuple(scalar_value(v) for v in pz)
        fz = tuple(scalar_value(v) for v in fz)
        f = fz[0] - pz[0]
        z = -focus_distance
        c = (pz[1] - z - pz[0]) * (pz[1] - z - 4.0 * f - pz[0])
        if c <= 0:
            raise ValueError(f"focus_distance {focus_distance:g} is too short")
        delta = (pz[1] - z + pz[0] - math.sqrt(c)) * 0.5
        return self.elements[-1].thickness + delta

    def _compute_exit_pupil_bounds(self):
        n = 64
        if not self.use_exit_pupil:
            rear = self.rear_element_radius()
            return [Bounds2(-rear, -rear, rear, rear) for _ in range(n)]
        bounds = []
        for i in range(n):
            r0 = i / n * self._film_diag_world * 0.5
            r1 = (i + 1) / n * self._film_diag_world * 0.5
            bounds.append(self.bound_exit_pupil(r0, r1))
        return bounds

    def bound_exit_pupil(self, film_x0, film_x1):
        rear_radius = self.rear_element_radius()
        proj = Bounds2(-1.5 * rear_radius, -1.5 * rear_radius, 1.5 * rear_radius, 1.5 * rear_radius)
        bound = Bounds2()
        vectorized = dr.is_jit_v(mi.Float)
        n_samples = max(1, self.exit_pupil_sample_count)
        chunk = min(65536, n_samples)
        if not vectorized:
            chunk = 1
        for start in range(0, n_samples, chunk):
            end = min(start + chunk, n_samples)
            u0, u1, s = self._radical_inverse_samples(start, end, n_samples)
            pfx = film_x0 + s * (film_x1 - film_x0)
            rx = proj.min_x + u0 * (proj.max_x - proj.min_x)
            ry = proj.min_y + u1 * (proj.max_y - proj.min_y)

            _, _, _, active = self._trace_lenses(
                pfx, mi.Float(0.0), 0.0,
                rx - pfx, ry, mi.Float(self.lens_rear_z()),
                from_film=True,
                active=mi.Mask(True),
            )

            self._merge_exit_bound_sample(bound, active, rx, ry)
        return self._expand_exit_pupil_bound(bound, proj, n_samples)

    def _expand_exit_pupil_bound(self, bound, proj, n_samples):
        if bound.is_degenerate():
            return bound
        diag = math.hypot(proj.max_x - proj.min_x, proj.max_y - proj.min_y)
        return bound.expand(2.0 * diag / math.sqrt(n_samples))

    def _sample_exit_pupil(self, pfx, pfy, u_lens):
        r_film = dr.sqrt(pfx * pfx + pfy * pfy)
        min_x, min_y, max_x, max_y, bound_ok = self._exit_bound_values(r_film)
        inv_r = dr.select(r_film != 0, 1.0 / r_film, 0.0)
        sin_theta = pfy * inv_r
        cos_theta = dr.select(r_film != 0, pfx * inv_r, 1.0)
        area = (max_x - min_x) * (max_y - min_y)
        px = min_x + u_lens.x * (max_x - min_x)
        py = min_y + u_lens.y * (max_y - min_y)
        pplx = cos_theta * px - sin_theta * py
        pply = sin_theta * px + cos_theta * py
        pupil = mi.Point3f(pplx, pply, mi.Float(self.lens_rear_z()))
        ok = bound_ok & (area > 0.0)
        return pupil, dr.select(ok, 1.0 / area, 0.0), ok

    def sample_ray(self, time, wavelength_sample, position_sample, aperture_sample, active=True):
        empty_si = mi.SurfaceInteraction3f()
        ray_wavelengths, ray_wav_weights = self.sample_wavelengths(empty_si, wavelength_sample, active)

        pixel_sample = (position_sample.x > 1.0) | (position_sample.y > 1.0)
        sx = dr.select(
            pixel_sample,
            (position_sample.x + self.crop_offset[0]) / self.full_resolution[0],
            (position_sample.x * self.crop_size[0] + self.crop_offset[0]) / self.full_resolution[0],
        )
        sy = dr.select(
            pixel_sample,
            (position_sample.y + self.crop_offset[1]) / self.full_resolution[1],
            (position_sample.y * self.crop_size[1] + self.crop_offset[1]) / self.full_resolution[1],
        )

        pfx = self.physical_extent.min_x + sx * (self.physical_extent.max_x - self.physical_extent.min_x)
        pfy = self.physical_extent.min_y + sy * (self.physical_extent.max_y - self.physical_extent.min_y)
        pupil, pdf, pupil_ok = self._sample_exit_pupil(pfx, pfy, aperture_sample)
        d_film = mi.Vector3f(pupil.x - pfx, pupil.y - pfy, pupil.z)
        wvl = ray_wavelengths[0] if ray_wavelengths.shape[0] > 0 else None

        weight, ro, rd, lens_ok = self._trace_lenses(
            pfx, pfy, 0.0,
            d_film.x, d_film.y, d_film.z,
            from_film=True,
            wavelength=wvl,
            active=active & pupil_ok,
        )
        if is_scalar_false(lens_ok):
            return mi.Ray3f(), mi.Spectrum(0.0)

        d_norm = dr.normalize(d_film)
        weight *= dr.power(d_norm.z, 4.0) / (pdf * self.lens_rear_z() * self.lens_rear_z())
        weight = dr.select(lens_ok, weight, 0.0)
        ray = self.world_transform() @ mi.Ray3f(ro, dr.normalize(rd), time, ray_wavelengths)
        if ray_wavelengths.shape[0] > 0:
            ray_wav_weights[1:] = 0.0
            ray_wav_weights[0] *= 4.0
        return ray, ray_wav_weights * weight

    def sample_ray_differential(self, time, sample1, sample2, sample3, active=True):
        ray, weight = self.sample_ray(time, sample1, sample2, sample3, active)
        rx, _ = self.sample_ray(
            time, sample1, mi.Point2f(sample2.x + 1.0 / self.crop_size[0], sample2.y), sample3, active
        )
        ry, _ = self.sample_ray(
            time, sample1, mi.Point2f(sample2.x, sample2.y + 1.0 / self.crop_size[1]), sample3, active
        )
        result = mi.RayDifferential3f(ray)
        result.o_x, result.d_x = rx.o, rx.d
        result.o_y, result.d_y = ry.o, ry.d
        result.has_differentials = True
        return result, weight

    def to_string(self):
        return (
            "RealisticCamera["
            f"lens_file={self.lens_file}, elements={len(self.elements)}, "
            f"focus_distance={self.focus_distance:g}, film_diagonal={self.film_diagonal:g}]"
        )


_registered_variants = set()


def register_realistic_camera():
    variant = mi.variant()
    if variant is None or variant in _registered_variants:
        return
    try:
        mi.register_sensor("realisticcamera", lambda props: RealisticCamera(props))
    except RuntimeError as exc:
        if "already" not in str(exc).lower() and "registered" not in str(exc).lower():
            raise
    _registered_variants.add(variant)


register_realistic_camera()
