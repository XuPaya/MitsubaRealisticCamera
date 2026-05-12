from __future__ import annotations

import math
import warnings

import drjit as dr
import mitsuba as mi
import numpy as np

import os

from utils import *

from glass_dictionary import sellmeier_ior


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
        # self.legacy_elements = self._make_elements(load_lens_file(self.lens_file))
        self.elements = self._load_elements(self.lens_file)

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

    def _make_elements(self, lens_parameters: list[float]) -> list[LensElementLegacy]:
        elements = []
        set_aperture = self.aperture_diameter * self.mm_to_world
        for i in range(0, len(lens_parameters), 4):
            curvature_radius = lens_parameters[i] * self.mm_to_world
            thickness = lens_parameters[i + 1] * self.mm_to_world
            eta = lens_parameters[i + 2]
            aperture_diameter = lens_parameters[i + 3] * self.mm_to_world
            if curvature_radius == 0:
                if set_aperture > aperture_diameter:
                    warnings.warn(
                        f"aperture_diameter {self.aperture_diameter:g} mm exceeds "
                        f"lens stop {aperture_diameter / self.mm_to_world:g} mm; clamping"
                    )
                else:
                    aperture_diameter = set_aperture
            elements.append(LensElementLegacy(curvature_radius, thickness, eta, 0.5 * aperture_diameter))
        return elements

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
            mi.Float(self.aperture_image.ravel().astype(np.float32))
            if self.aperture_image is not None
            else None
        )
    
    def _load_elements(self, filename: str | os.PathLike) -> list[LensElement]:
        
        lens_parameters = []
        with open(filename, "r", encoding="utf-8") as f:
            for line_no, line in enumerate(f, 1):
                line = line.split("#", 1)[0].strip()
                if not line:
                    continue
                parts = line.replace(",", " ").split()
                if len(parts) != 4:
                    raise ValueError(f"{filename}:{line_no}: expected 4 values per line, got {len(parts)}")
                lens_parameters.append(float(parts[0])) # curvature radius
                lens_parameters.append(float(parts[1])) # thickness
                try: 
                    eta_value = float(parts[2])
                    lens_parameters.append(Eta_lookup({"type": "constant", "value": eta_value})) # eta value
                except ValueError:
                    lens_parameters.append(Eta_lookup({"type": "glass", "glass_name": parts[2]})) # eta name
                lens_parameters.append(float(parts[3])) # aperture diameter
        if not lens_parameters or len(lens_parameters) % 4:
            raise ValueError(f"{filename}: lens files must contain groups of four floats")
    
        elements = []
        set_aperture = self.aperture_diameter * self.mm_to_world
        for i in range(0, len(lens_parameters), 4):
            curvature_radius = lens_parameters[i] * self.mm_to_world
            thickness = lens_parameters[i + 1] * self.mm_to_world
            eta = lens_parameters[i + 2] # Eta_lookup object
            aperture_diameter = lens_parameters[i + 3] * self.mm_to_world
            if curvature_radius == 0:
                if set_aperture > aperture_diameter:
                    warnings.warn(
                        f"aperture_diameter {self.aperture_diameter:g} mm exceeds "
                        f"lens stop {aperture_diameter / self.mm_to_world:g} mm; clamping"
                    )
                else:
                    aperture_diameter = set_aperture
                
            elements.append(LensElement(curvature_radius, thickness, eta, 0.5 * aperture_diameter))
        return elements

    def lens_rear_z(self) -> float:
        return self.elements[-1].thickness

    def lens_front_z(self) -> float:
        return sum(e.thickness for e in self.elements)

    def rear_element_radius(self) -> float:
        return self.elements[-1].aperture_radius

    def trace_lenses_from_film(self, origin, direction):
        ox, oy, oz = origin[0], origin[1], -origin[2]
        dx, dy, dz = direction[0], direction[1], -direction[2]
        element_z = 0.0
        weight = 1.0

        for i in range(len(self.elements) - 1, -1, -1):
            element = self.elements[i]
            element_z -= element.thickness
            is_stop = element.curvature_radius == 0
            if is_stop:
                if dz == 0:
                    return 0.0, None, None
                t = (element_z - oz) / dz
                if t < -RAY_T_EPSILON:
                    return 0.0, None, None
            else:
                hit = intersect_spherical_element(
                    element.curvature_radius,
                    element_z + element.curvature_radius,
                    (ox, oy, oz),
                    (dx, dy, dz),
                )
                if hit is None:
                    return 0.0, None, None
                t, n = hit

            hx, hy, hz = ox + t * dx, oy + t * dy, oz + t * dz
            if is_stop and self.aperture_image is not None:
                weight = float(aperture_lookup_np(
                    self.aperture_image,
                    np.array([(hx / element.aperture_radius + 1.0) * 0.5]),
                    np.array([(hy / element.aperture_radius + 1.0) * 0.5]),
                )[0])
                if weight == 0:
                    return 0.0, None, None
            elif hx * hx + hy * hy > element.aperture_radius * element.aperture_radius:
                return 0.0, None, None

            ox, oy, oz = hx, hy, hz
            if not is_stop:
                eta_i = element.eta(None)
                eta_t = self.elements[i - 1].eta(None) if i > 0 and self.elements[i - 1].eta(None)>0.0 else 1.0
                wt = refract(normalize((-dx, -dy, -dz)), n, eta_t / eta_i)
                if wt is None:
                    return 0.0, None, None
                dx, dy, dz = wt

        return weight, (ox, oy, -oz), (dx, dy, -dz)

    def trace_lenses_from_scene(self, origin, direction):
        ox, oy, oz = origin[0], origin[1], -origin[2]
        dx, dy, dz = direction[0], direction[1], -direction[2]
        element_z = -self.lens_front_z()
        for i, element in enumerate(self.elements):
            is_stop = element.curvature_radius == 0
            if is_stop:
                if dz == 0:
                    return 0.0, None, None
                t = (element_z - oz) / dz
                if t < -RAY_T_EPSILON:
                    return 1.0, None, None
            else:
                hit = intersect_spherical_element(
                    element.curvature_radius,
                    element_z + element.curvature_radius,
                    (ox, oy, oz),
                    (dx, dy, dz),
                )
                if hit is None:
                    return 2.0, None, None
                t, n = hit
            hx, hy, hz = ox + t * dx, oy + t * dy, oz + t * dz
            if hx * hx + hy * hy > element.aperture_radius * element.aperture_radius:
                return 0.0, None, None
            ox, oy, oz = hx, hy, hz
            if not is_stop:
                eta_i = self.elements[i - 1].eta(None) if i > 0 and self.elements[i - 1].eta(None)>0.0 else 1.0
                eta_t = element.eta(None) if element.eta(None)>0.0 else 1.0
                wt = refract(normalize((-dx, -dy, -dz)), n, eta_t / eta_i)
                if wt is None:
                    return 3.0, None, None
                dx, dy, dz = wt
            element_z += element.thickness
        return 1.0, (ox, oy, -oz), (dx, dy, -dz)

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
        n_samples = max(1, self.exit_pupil_sample_count)
        chunk = min(65536, n_samples)
        for start in range(0, n_samples, chunk):
            end = min(start + chunk, n_samples)
            idx = np.arange(start, end, dtype=np.uint64)
            u0 = radical_inverse(2, idx)
            u1 = radical_inverse(3, idx)
            s = (idx.astype(np.float64) + 0.5) / n_samples
            pfx = film_x0 + s * (film_x1 - film_x0)
            rx = proj.min_x + u0 * (proj.max_x - proj.min_x)
            ry = proj.min_y + u1 * (proj.max_y - proj.min_y)
            active = self._trace_lenses_from_film_np(
                pfx,
                np.zeros_like(pfx),
                np.zeros_like(pfx),
                rx - pfx,
                ry,
                np.full_like(pfx, self.lens_rear_z()),
            )
            
            if np.any(active):
                bound.min_x = min(bound.min_x, float(np.min(rx[active])))
                bound.min_y = min(bound.min_y, float(np.min(ry[active])))
                bound.max_x = max(bound.max_x, float(np.max(rx[active])))
                bound.max_y = max(bound.max_y, float(np.max(ry[active])))
        if bound.is_degenerate():
            return bound
        diag = math.hypot(proj.max_x - proj.min_x, proj.max_y - proj.min_y)
        return bound.expand(2.0 * diag / math.sqrt(n_samples))

    def _sample_exit_pupil_scalar(self, p_film, u_lens):
        r_film = math.hypot(p_film[0], p_film[1])
        r_index = min(len(self.exit_pupil_bounds) - 1, int(r_film / (self._film_diag_world * 0.5) * len(self.exit_pupil_bounds)))
        bound = self.exit_pupil_bounds[r_index]
        if bound.is_degenerate():
            return None, 0.0
        px = bound.min_x + float(u_lens[0]) * (bound.max_x - bound.min_x)
        py = bound.min_y + float(u_lens[1]) * (bound.max_y - bound.min_y)
        sin_theta = p_film[1] / r_film if r_film != 0 else 0.0
        cos_theta = p_film[0] / r_film if r_film != 0 else 1.0
        pupil = (
            cos_theta * px - sin_theta * py,
            sin_theta * px + cos_theta * py,
            self.lens_rear_z(),
        )
        return pupil, 1.0 / bound.area

    def _sample_exit_pupil_vec(self, pfx, pfy, u_lens):
        r_film = dr.sqrt(pfx * pfx + pfy * pfy)
        n = len(self.exit_pupil_bounds)
        idx = mi.UInt32(dr.minimum(mi.Float(n - 1), r_film / (self._film_diag_world * 0.5) * n))
        min_x = dr.gather(mi.Float, self._exit_min_x, idx)
        min_y = dr.gather(mi.Float, self._exit_min_y, idx)
        max_x = dr.gather(mi.Float, self._exit_max_x, idx)
        max_y = dr.gather(mi.Float, self._exit_max_y, idx)
        area = (max_x - min_x) * (max_y - min_y)
        px = min_x + u_lens.x * (max_x - min_x)
        py = min_y + u_lens.y * (max_y - min_y)
        inv_r = dr.select(r_film != 0, 1.0 / r_film, 0.0)
        sin_theta = pfy * inv_r
        cos_theta = dr.select(r_film != 0, pfx * inv_r, 1.0)
        pplx = cos_theta * px - sin_theta * py
        pply = sin_theta * px + cos_theta * py
        pplz = dr.auto.Float(self.lens_rear_z())
        pupil = mi.Point3f(
            pplx,
            pply,
            pplz,
        )
        return pupil, dr.select(area > 0.0, 1.0 / area, 0.0), area > 0.0

    def sample_ray(self, time, wavelength_sample, position_sample, aperture_sample, active=True):
        if not dr.is_jit_v(mi.Float):
            return self._sample_ray_scalar(time, wavelength_sample, position_sample, aperture_sample, active)
        return self._sample_ray_vec(time, wavelength_sample, position_sample, aperture_sample, active)

    def _sample_ray_scalar(self, time, wavelength_sample, position_sample, aperture_sample, active=True):
        if float(position_sample.x) > 1.0 or float(position_sample.y) > 1.0:
            sx = (float(position_sample.x) + self.crop_offset[0]) / self.full_resolution[0]
            sy = (float(position_sample.y) + self.crop_offset[1]) / self.full_resolution[1]
        else:
            sx = (float(position_sample.x) * self.crop_size[0] + self.crop_offset[0]) / self.full_resolution[0]
            sy = (float(position_sample.y) * self.crop_size[1] + self.crop_offset[1]) / self.full_resolution[1]
        p2x = self.physical_extent.min_x + sx * (self.physical_extent.max_x - self.physical_extent.min_x)
        p2y = self.physical_extent.min_y + sy * (self.physical_extent.max_y - self.physical_extent.min_y)
        p_film = (p2x, p2y, 0.0)
        p_pupil, pdf = self._sample_exit_pupil_scalar((p_film[0], p_film[1]), aperture_sample)
        if p_pupil is None:
            return mi.Ray3f(), mi.Spectrum(0.0)
        d_film = (p_pupil[0] - p_film[0], p_pupil[1] - p_film[1], p_pupil[2])
        weight, ro, rd = self.trace_lenses_from_film(p_film, d_film)
        if weight == 0.0:
            return mi.Ray3f(), mi.Spectrum(0.0)
        cos_theta = normalize(d_film)[2]
        weight *= cos_theta**4 / (pdf * self.lens_rear_z() * self.lens_rear_z())
        rd = normalize(rd)
        
        empty_si = mi.SurfaceInteraction3f()
        ray_wavelengths, ray_wav_weights = self.sample_wavelengths(empty_si, wavelength_sample, active)
        
        ray = mi.Ray3f(mi.Point3f(*ro), mi.Vector3f(*rd), float(time), ray_wavelengths)
        ray = self.world_transform() @ ray
        return ray, ray_wav_weights * weight

    def _sample_ray_vec(self, time, wavelength_sample, position_sample, aperture_sample, active=True):
        empty_si = mi.SurfaceInteraction3f()
        ray_wavelengths, ray_wav_weights = self.sample_wavelengths(empty_si, wavelength_sample, active)
        
        sx = (position_sample.x * self.crop_size[0] + self.crop_offset[0]) / self.full_resolution[0]
        sy = (position_sample.y * self.crop_size[1] + self.crop_offset[1]) / self.full_resolution[1]
        p2x = self.physical_extent.min_x + sx * (self.physical_extent.max_x - self.physical_extent.min_x)
        p2y = self.physical_extent.min_y + sy * (self.physical_extent.max_y - self.physical_extent.min_y)
        pfx, pfy = p2x, p2y
        pupil, pdf, pupil_ok = self._sample_exit_pupil_vec(pfx, pfy, aperture_sample)
        dfx, dfy, dfz = pupil.x - pfx, pupil.y - pfy, pupil.z
        weight, ro, rd, lens_ok = self._trace_lenses_from_film_vec(pfx, pfy, 0.0, dfx, dfy, dfz, ray_wavelengths, active & pupil_ok)
        d_norm = dr.normalize(mi.Vector3f(dfx, dfy, dfz))
        weight *= dr.power(d_norm.z, 4.0) / (pdf * self.lens_rear_z() * self.lens_rear_z())
        weight = dr.select(lens_ok, weight, 0.0)
        
        ray = mi.Ray3f(ro, dr.normalize(rd), time, ray_wavelengths)
        ray = self.world_transform() @ ray
        if ray_wavelengths.shape[0] > 0: # spectral rendering enabled: ray is warped according to first sampled wavelength
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

    def _trace_lenses_from_film_vec(self, ox, oy, oz, dx, dy, dz, ray_wavelengths, active):
        ox, oy, oz = ox, oy, -oz
        dx, dy, dz = dx, dy, -dz
        active = mi.Mask(active)
        weight = mi.Float(1.0)
        element_z = 0.0

        for i in range(len(self.elements) - 1, -1, -1):
            element = self.elements[i]
            element_z -= element.thickness
            is_stop = element.curvature_radius == 0
            if is_stop:
                denom_ok = dz != 0.0
                t = (element_z - oz) / dr.select(denom_ok, dz, 1.0)
                hit = denom_ok & (t >= -RAY_T_EPSILON)
                nx = ny = nz = mi.Float(0.0)
            else:
                t, nx, ny, nz, hit = intersect_spherical_element_vec(
                    element.curvature_radius, element_z + element.curvature_radius, ox, oy, oz, dx, dy, dz
                )
            hx, hy, hz = ox + t * dx, oy + t * dy, oz + t * dz
            aperture_ok = hx * hx + hy * hy <= element.aperture_radius * element.aperture_radius
            if is_stop and self._aperture_flat is not None:
                aw = aperture_lookup_vec(
                    self.aperture_image,
                    self._aperture_flat,
                    (hx / element.aperture_radius + 1.0) * 0.5,
                    (hy / element.aperture_radius + 1.0) * 0.5,
                )
                weight = aw
                aperture_ok = aw != 0.0
            active &= hit & aperture_ok
            ox, oy, oz = hx, hy, hz
            if not is_stop:
                wi = dr.normalize(mi.Vector3f(-dx, -dy, -dz))
                if ray_wavelengths.shape[0] > 0:
                    wvl = ray_wavelengths[0]
                else:
                    wvl = None                
                eta_i = element.eta(wvl)
                eta_t = self.elements[i - 1].eta(wvl) if i > 0 and self.elements[i - 1].eta(None)>0.0 else 1.0
                wt, refr_ok = refract_vec(wi, mi.Vector3f(nx, ny, nz), eta_t / eta_i)
                active &= refr_ok
                dx, dy, dz = wt.x, wt.y, wt.z
        return (
            weight,
            mi.Point3f(ox, oy, -oz),
            mi.Vector3f(dx, dy, -dz),
            active,
        )

    def _trace_lenses_from_film_np(self, ox, oy, oz, dx, dy, dz):
        ox, oy, oz = ox.copy(), oy.copy(), -oz.copy()
        dx, dy, dz = dx.copy(), dy.copy(), -dz.copy()
        active = np.ones_like(ox, dtype=bool)
        element_z = 0.0
        for i in range(len(self.elements) - 1, -1, -1):
            element = self.elements[i]
            element_z -= element.thickness
            if element.curvature_radius == 0:
                denom_ok = dz != 0.0
                t = (element_z - oz) / np.where(denom_ok, dz, 1.0)
                hit = denom_ok & (t >= -RAY_T_EPSILON)
                n = None
            else:
                t, n, hit = intersect_spherical_element_np(
                    element.curvature_radius, element_z + element.curvature_radius, ox, oy, oz, dx, dy, dz
                )
            hx, hy, hz = ox + t * dx, oy + t * dy, oz + t * dz
            aperture_ok = hx * hx + hy * hy <= element.aperture_radius * element.aperture_radius
            if element.curvature_radius == 0 and self.aperture_image is not None:
                aw = aperture_lookup_np(
                    self.aperture_image,
                    (hx / element.aperture_radius + 1.0) * 0.5,
                    (hy / element.aperture_radius + 1.0) * 0.5,
                )
                aperture_ok = aw != 0.0
            active &= hit & aperture_ok
            ox, oy, oz = hx, hy, hz
            if element.curvature_radius != 0:
                inv_len = 1.0 / np.sqrt(dx * dx + dy * dy + dz * dz)
                wi = (-dx * inv_len, -dy * inv_len, -dz * inv_len)
                eta_i = element.eta(None)
                eta_t = self.elements[i - 1].eta(None) if i > 0 and self.elements[i - 1].eta(None)>0.0 else 1.0
                wt, ok = refract_np(wi, n, eta_t / eta_i)
                
                active &= ok
                dx, dy, dz = wt
        return active


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
