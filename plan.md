# Plan: Mitsuba RealisticCamera + Spectral Lens + CMOS/CCD Film

## Project goal

Implement a Mitsuba camera/sensor pipeline in three stages:

1. **Phase 1:** Port a PBRT-v4-style `RealisticCamera` to Mitsuba, including PBRT `.dat` lens loader and constant-index lens tracing.
2. **Phase 2:** Extend the lens model with wavelength-dependent glass IOR, chromatic aberration, and optional 2011-style spectral lens effects such as internal reflection ghosting / flare.
3. **Phase 3:** Add a CMOS/CCD sensor-response/noise pipeline as a film/development module, based on the CCD/CMOS simulation paper.

The implementation should be incremental. Do **not** start Phase 2 until Phase 1 has a working constant-eta camera with reproducible tests. Do **not** start Phase 3 until Phase 2A chromatic aberration works.

---

## Important roadmap adjustment

The original roadmap is valid, but Phase 2 is too broad if interpreted as “implement every effect from the 2011 lens simulation paper.” Treat Phase 2 as two sub-phases:

- **Phase 2A, required:** wavelength-dependent IOR + chromatic aberration.
- **Phase 2B, optional:** flare/ghosting/coatings/aperture diffraction-like effects.

Do **not** block the main camera implementation on flare or diffraction. Those are much less stable and much easier to turn into a noisy science-fair volcano.

---

## Target architecture

### Plugin split

Use this separation of responsibility:

```text
Mitsuba Sensor plugin: realisticcamera
    - film/sample to lens ray generation
    - lens prescription loading
    - spherical element intersection
    - refraction through lens stack
    - aperture stop / housing clipping
    - exit pupil sampling
    - focus adjustment
    - wavelength-dependent IOR in Phase 2
    - rolling shutter only if needed later, but not in Phase 1

Mitsuba Film or development module: cmosfilm / cmos_develop
    - spectral sensor response / RGB SRFs
    - photons -> electrons
    - shot noise
    - read noise
    - dark current
    - fixed pattern noise
    - full-well saturation
    - ADC quantization
    - optional blooming
```

For Phase 3, begin with an offline Python or C++ postprocess/development step if a full Mitsuba `Film` plugin is too invasive. A postprocess path is acceptable as long as the file format and units are documented.

---

## Baseline references to inspect before coding

Codex should inspect these files/docs locally or online before implementation:

```text
PBRT v4:
  src/pbrt/cameras.h
  src/pbrt/cameras.cpp
  src/pbrt/util/sampling.h
  src/pbrt/util/geometry.h
  scenes/lenses/*.dat from pbrt-v4-scenes or pbrt-v3-scenes

Mitsuba 3:
  src/sensors/perspective.cpp
  src/sensors/thinlens.cpp, if present
  include/mitsuba/render/sensor.h
  src/render/sensor.cpp
  src/films/hdrfilm.cpp
  src/films/specfilm.cpp, if present
  plugin registration macros / CMake plugin conventions
```

If Mitsuba source layout differs, search for existing sensor plugins and mirror their registration pattern.

---

## Global implementation rules

1. Keep units explicit.
   - Lens files use **millimeters**.
   - Mitsuba world space often uses scene units; do not silently mix world meters and lens millimeters.
   - Store a `world_scale` or `mm_to_world` parameter if necessary.

2. Keep coordinate conventions isolated.
   - PBRT and Mitsuba may use different camera-space conventions.
   - Implement lens tracing in an internal `LensCameraSpace` that matches PBRT as closely as possible.
   - Convert to Mitsuba camera/world space only at the sensor boundary.

3. Preserve PBRT `.dat` compatibility in Phase 1.
   - Format:
     ```text
     radius_mm thickness_mm eta aperture_radius_mm
     ```
   - `radius_mm == 0` and/or `eta == 0` represents aperture stop, depending on PBRT convention.
   - Comments and blank lines must be ignored.

4. Phase 1 must support constant-eta lenses only.
   - Do not add glass dictionaries in Phase 1.
   - Do not add wavelength-dependent refraction in Phase 1.

5. Phase 2 must preserve Phase 1 behavior.
   - Constant-eta `.dat` files must still render exactly or near-exactly as before.

6. Write tests before large feature additions.
   - At minimum, add parser tests, spherical intersection tests, refraction tests, and rendered sanity scenes.

---

# Phase 1: PBRT-style RealisticCamera in Mitsuba

## Phase 1 objective

Implement a Mitsuba `realisticcamera` sensor that reproduces the PBRT-v4 `RealisticCamera` model closely enough for:

- compound spherical lens tracing;
- PBRT `.dat` lens loading;
- aperture stop clipping;
- lens housing clipping;
- focus adjustment;
- exit pupil sampling;
- bokeh, cat’s-eye bokeh, vignetting, distortion;
- comparison against Mitsuba `thinlens` and, if possible, PBRT output.

Do not implement chromatic aberration yet.

---

## Phase 1 deliverables

### Source files

Create or adapt file names to match the repo layout:

```text
src/sensors/realisticcamera.cpp
include/mitsuba/render/realistic_lens.h       # optional helper header
src/libcore/realistic_lens.cpp                # optional helper implementation
resources/lenses/wide.22mm.dat
resources/lenses/dgauss.dat                   # if available
tests/test_realistic_camera.py
tools/visualize_lens.py                       # optional but strongly recommended
```

### Scene examples

```text
scenes/realistic_camera/grid_distortion.xml
scenes/realistic_camera/bokeh_lights.xml
scenes/realistic_camera/vignetting_flatfield.xml
scenes/realistic_camera/focus_sweep.xml
```

### Documentation

```text
docs/realistic_camera.md
```

This document should describe:

- lens file format;
- coordinate convention;
- supported Mitsuba variants;
- known unsupported sensor methods;
- examples.

---

## Phase 1 plugin XML interface

Target XML shape:

```xml
<sensor type="realisticcamera">
    <string name="lens_file" value="resources/lenses/wide.22mm.dat"/>
    <float name="aperture_diameter" value="8.0"/>
    <float name="focus_distance" value="3.0"/>
    <float name="film_diagonal" value="35.0"/>
    <float name="mm_to_world" value="0.001"/>
    <boolean name="use_exit_pupil" value="true"/>
    <boolean name="debug_lens_trace" value="false"/>

    <sampler type="independent">
        <integer name="sample_count" value="64"/>
    </sampler>

    <film type="hdrfilm">
        <integer name="width" value="1024"/>
        <integer name="height" value="768"/>
    </film>
</sensor>
```

Parameter meanings:

| Parameter | Type | Phase | Meaning |
|---|---:|---:|---|
| `lens_file` | string | 1 | PBRT `.dat` lens prescription path |
| `aperture_diameter` | float | 1 | Effective aperture diameter in mm |
| `focus_distance` | float | 1 | Focus distance in scene units or meters; document choice clearly |
| `film_diagonal` | float | 1 | Sensor/film diagonal in mm |
| `mm_to_world` | float | 1 | Conversion from lens mm to Mitsuba world units |
| `use_exit_pupil` | bool | 1 | Use PBRT-style exit pupil bounds sampling |
| `debug_lens_trace` | bool | 1 | Print/emit debug info for failed rays |
| `aperture_image` | string | 1 optional | Non-circular aperture image mask, if easy to port |

---

## Phase 1 data structures

Implement these concepts first in a renderer-independent helper module.

```cpp
struct LensElement {
    float radius_mm;           // spherical curvature radius; 0 for stop
    float thickness_mm;        // distance to next surface
    float eta;                 // constant IOR of medium after surface
    float aperture_radius_mm;  // physical aperture/lens housing radius

    bool is_stop() const;
};

struct LensSystem {
    std::vector<LensElement> elements;
    float film_diagonal_mm;
    float aperture_diameter_mm;
    float focus_distance_world;
    float mm_to_world;
};
```

Do not add glass names or Sellmeier coefficients yet.

---

## Phase 1 implementation steps

### Step 1. Build and plugin skeleton

1. Add a new Mitsuba sensor plugin named `realisticcamera`.
2. Start by copying the structure of an existing sensor plugin, preferably `thinlens` or `perspective`.
3. Implement constructor parameter parsing.
4. Implement `to_string()` / debug print.
5. Return a simple pinhole ray first to verify plugin registration.

Acceptance:

- Scene loads with `<sensor type="realisticcamera">`.
- A pinhole fallback render matches `perspective` qualitatively.

---

### Step 2. PBRT `.dat` loader

Implement a robust parser:

- skip empty lines;
- skip `#` comments;
- parse four floats per lens element;
- preserve lens order exactly as stored;
- report line number on error;
- validate nonnegative aperture radius;
- allow zero radius for stop;
- allow zero eta only for stop;
- allow negative curvature radius.

Add tests:

```text
test_parse_wide_22mm_dat
test_parse_comments_and_blank_lines
test_reject_invalid_column_count
test_reject_negative_aperture_radius
test_identify_stop_element
```

Acceptance:

- Parser can load PBRT lens files.
- `LensSystem::to_string()` prints the same values in order.

---

### Step 3. Lens-space coordinate convention

Implement an internal coordinate convention and document it in code.

Recommended:

- Use PBRT’s lens tracing convention internally.
- Implement helper functions:

```cpp
Point3f film_to_lens_space(Point2f film_sample_mm);
Ray3f lens_space_to_mitsuba_camera_ray(Ray3f lens_ray_mm);
Ray3f mitsuba_camera_to_lens_space(Ray3f camera_ray_world);
```

Do not scatter sign flips throughout the code. One sign flip hidden in the wrong place will waste a day and then have the audacity to pass half the tests.

Acceptance:

- Central film point with central aperture sample produces a ray along the optical axis.
- Symmetric film samples produce symmetric outgoing rays.

---

### Step 4. Spherical surface intersection

Port PBRT’s spherical element intersection logic into a helper:

```cpp
bool intersect_spherical_element(
    float radius_mm,
    float z_center_mm,
    const Ray3f &ray_mm,
    float *t,
    Normal3f *normal);
```

Requirements:

- support positive and negative radius;
- choose correct root based on ray direction and curvature sign;
- reject negative `t`;
- return a normal facing against incoming ray or document convention;
- test with analytic rays.

Tests:

```text
test_intersect_positive_radius_on_axis
test_intersect_negative_radius_on_axis
test_intersect_off_axis_inside_aperture
test_intersect_miss
test_normal_orientation
```

---

### Step 5. Refraction through spherical interfaces

Implement Snell refraction:

```cpp
bool refract_dielectric(
    const Vector3f &wi,
    const Normal3f &n,
    float eta_i,
    float eta_t,
    Vector3f *wt);
```

Requirements:

- handle air-to-glass and glass-to-air;
- handle total internal reflection by returning false;
- use normalized directions;
- match Mitsuba/PBRT direction sign conventions;
- do not apply Fresnel weighting in Phase 1 main imaging path unless PBRT does so in the target path.

Tests:

```text
test_refraction_normal_incidence
test_refraction_oblique_air_to_glass
test_refraction_oblique_glass_to_air
test_total_internal_reflection
```

---

### Step 6. Trace lenses from film

Implement:

```cpp
std::optional<LensTraceResult> trace_lenses_from_film(
    const Ray3f &ray_from_film_mm,
    float wavelength_nm_or_dummy);
```

Phase 1 ignores wavelength.

Algorithm:

1. Start with a ray from film plane toward rear lens element.
2. For each element from rear to front:
   - advance to the element plane/surface;
   - if stop: intersect stop plane, reject if radius exceeds stop aperture;
   - if spherical element:
     - intersect sphere;
     - reject if intersection radius exceeds aperture radius;
     - refract from current medium to next medium;
   - update ray origin/direction.
3. If ray exits front lens, convert to Mitsuba camera/world ray.
4. Return throughput/weight according to PBRT-style camera sampling.

Debug counters:

```text
num_lens_rays
num_stop_rejected
num_housing_rejected
num_tir_rejected
num_success
```

Acceptance:

- Central rays pass through simple biconvex lens.
- Off-axis rays are clipped at high field angles.
- A flat white wall render shows vignetting for wide lenses.

---

### Step 7. Sampling strategy

Implement two modes.

#### Mode A: simple rear-element sampling

- Sample a disk on the rear element or nominal aperture.
- Trace through lens.
- Rejection sampling handles invalid rays.

This is easier and should be implemented first.

#### Mode B: PBRT-style exit pupil sampling

Port PBRT’s approach:

- precompute exit pupil bounds as a function of film x coordinate;
- sample bounds for current film point;
- trace ray through lens;
- return PDF/weight correction.

Acceptance:

- Mode A works but may be noisy.
- Mode B has fewer rejected rays and stable edge illumination.

---

### Step 8. Focus adjustment

Port or reimplement PBRT’s thick-lens focus approximation:

```cpp
void compute_thick_lens_approximation(...);
float focus_thick_lens(float focus_distance);
```

Acceptance scenes:

- Focus target at `focus_distance` appears sharp.
- Near/far targets blur correctly.
- Focus sweep changes the sharp plane monotonically.

If focus adjustment is too time-consuming, provide manual `film_distance_offset_mm` as temporary override and mark thick-lens focus as TODO. Do not let focus perfection block the whole plugin.

---

### Step 9. Rendered Phase 1 sanity tests

Create scenes:

1. **Grid distortion**
   - Render a planar checker/grid.
   - Compare `perspective`, `thinlens`, `realisticcamera`.
   - Expected: RealisticCamera may show barrel/pincushion/fisheye-like distortion depending on lens.

2. **Bokeh lights**
   - Small bright spheres or disks at different depths.
   - Expected: defocus blur and aperture-shaped bokeh.

3. **Cat’s-eye bokeh**
   - Bright out-of-focus lights near frame edge.
   - Expected: edge bokeh clipped into cat’s-eye shapes.

4. **Flat field / vignetting**
   - Uniform white plane or environment.
   - Expected: edge falloff for appropriate lens.

5. **Focus sweep**
   - Slanted plane or row of objects.
   - Expected: sharp region moves with focus distance.

---

## Phase 1 Definition of Done

Phase 1 is done when all of the following are true:

- `realisticcamera` plugin compiles and loads.
- PBRT `.dat` lens files load correctly.
- Simple lens tracing works for at least one toy biconvex lens and one PBRT lens.
- Rendered images show at least three of:
  - distortion;
  - vignetting;
  - depth of field;
  - cat’s-eye bokeh;
  - focus sweep.
- Unit tests pass.
- Documentation explains coordinate convention and known limitations.

Known acceptable limitations at end of Phase 1:

- no chromatic aberration;
- no flare;
- no diffraction;
- no sensor noise;
- no aspheric elements;
- no polarization;
- `sample_direction` / sensor importance methods may be unimplemented if not needed by target integrators.

---

# Phase 2: Spectral glass, chromatic aberration, optional flare

## Phase 2 objective

Extend the Phase 1 camera from constant-index lens tracing to wavelength-dependent glass tracing. Required output is visible chromatic aberration under spectral rendering.

---

## Phase 2A required deliverables

```text
resources/glass/schott_basic.yaml
resources/lenses/wide_22mm_spectral.yaml
tests/test_glass_catalog.py
tests/test_chromatic_aberration.py
scenes/realistic_camera/chromatic_grid.xml
scenes/realistic_camera/longitudinal_ca_focus.xml
```

---

## Phase 2A data model

Keep PBRT `.dat` loader, but add a richer YAML/JSON lens format.

Example:

```yaml
units: mm
elements:
  - radius_mm: 35.98738
    thickness_mm: 1.21638
    medium_after: N-BK7
    aperture_radius_mm: 23.716

  - radius_mm: 11.69718
    thickness_mm: 9.99570
    medium_after: AIR
    aperture_radius_mm: 17.996

  - stop: true
    thickness_mm: 2.27766
    aperture_radius_mm: 8.756
```

Glass dictionary:

```yaml
AIR:
  model: constant
  n: 1.0

N-BK7:
  model: sellmeier
  wavelength_unit: um
  range_um: [0.3, 2.5]
  B: [1.03961212, 0.231792344, 1.01046945]
  C: [0.00600069867, 0.0200179144, 103.560653]
```

Implement only these models first:

```text
constant
sellmeier_3term
```

Optional later:

```text
cauchy
sellmeier_extended
tabulated_nk
```

---

## Phase 2A glass IOR implementation

Implement:

```cpp
class GlassModel {
public:
    float ior(float wavelength_nm) const;
};
```

Sellmeier formula:

```text
n^2(lambda) = 1 + sum_i B_i * lambda^2 / (lambda^2 - C_i)
```

Use micrometers if the catalog says micrometers.

Tests:

```text
test_air_constant_ior
test_bk7_ior_587nm_approx
test_sellmeier_rejects_out_of_range_wavelength
test_unknown_glass_name_error
test_constant_eta_dat_compatibility
```

---

## Phase 2A Mitsuba spectral integration

Implementation approach:

1. Start with `scalar_spectral` if possible.
2. Retrieve sampled wavelength(s) from Mitsuba’s spectral path state / sensor sampling API.
3. For each primary camera ray, trace the lens using the active wavelength.
4. If Mitsuba carries multiple wavelengths per sample, either:
   - trace one ray per wavelength and return spectral contribution consistently, or
   - initially restrict plugin to scalar spectral mode and document limitation.

Do not fake CA by shifting RGB channels in screen space. That is a fallback demo only, not the actual Phase 2 feature.

---

## Phase 2A chromatic aberration tests

1. **Lateral CA test**
   - White grid or black-white edges near frame edge.
   - Expected: colored fringes near high-contrast edges.
   - Compare constant-eta vs spectral glass.

2. **Longitudinal CA test**
   - White point lights or slanted high-contrast plane.
   - Expected: different wavelengths focus at different depths.

3. **Aperture dependence test**
   - Render same scene at different aperture diameters.
   - Expected: axial CA / defocus-related fringing changes with aperture.

4. **Spot diagram tool**
   - Trace monochromatic bundles for 450nm, 550nm, 650nm.
   - Save 2D points at film plane.
   - Expected: different centroids or blur sizes.

---

## Phase 2B optional: flare / ghosting

Treat flare as optional and separate from the main imaging path.

Recommended first implementation: two-reflection ghosting only.

### Feature interface

```xml
<boolean name="enable_ghosting" value="true"/>
<integer name="max_ghost_reflections" value="2"/>
<float name="coating_reflectance_scale" value="0.05"/>
<float name="ghost_exposure_scale" value="1.0"/>
<boolean name="output_ghost_aov" value="true"/>
```

### Algorithm

For each selected pair of lens surfaces `(i, j)`:

1. Trace a path where all surfaces transmit except `i` and `j`, which reflect.
2. Weight by Fresnel reflectance at reflected surfaces and transmission at other surfaces.
3. Clip by apertures/housing.
4. Accumulate to a separate ghost contribution.
5. Composite:

```text
final = main_image + ghost_scale * ghost_image
```

Start with strong point/directional lights only. Do not attempt full-scene unbiased ghosting first.

### Acceptance

- Rendering with a bright light source shows faint ghost spots/images.
- Disabling ghosting returns exactly the Phase 2A image.
- Ghost contribution can be output separately.

---

## Phase 2B optional: aperture diffraction / glare

Do not implement full wave optics unless time remains.

Acceptable approximations:

1. Postprocess diffraction PSF for very bright pixels.
2. Aperture-shape Fourier approximation as a separate image-space effect.
3. Clearly label it as hybrid approximation, not full lens wave simulation.

Do not mix this into the geometric lens tracer in a way that breaks Phase 2A.

---

## Phase 2 Definition of Done

Phase 2A is done when:

- YAML/JSON spectral lens files load.
- Glass dictionary loads.
- Constant `.dat` compatibility remains intact.
- Spectral lens tracing changes ray directions as a function of wavelength.
- Rendered images show visible CA compared to constant-eta baseline.
- Spot diagram or numeric test demonstrates wavelength-dependent focus/projection.

Phase 2B is optional and done only if:

- ghosting can be toggled;
- ghost AOV can be separated;
- the main image path remains stable;
- report clearly states approximations.

---

# Phase 3: CMOS/CCD sensor response and noise

## Phase 3 objective

Implement a physically motivated CMOS/CCD sensor response pipeline as either:

1. a Mitsuba `Film` plugin named `cmosfilm`, or
2. an offline development tool named `cmos_develop` that consumes EXR/specfilm output.

Prefer the offline development path first if Film plugin integration becomes risky.

---

## Phase 3A recommended first target: offline develop tool

Create:

```text
tools/cmos_develop.py
resources/sensors/example_cmos.yaml
resources/sensors/example_rgb_srf.csv
tests/test_cmos_develop.py
```

Input:

```text
linear EXR from Mitsuba, preferably high-spp and linear radiance/exposure-like units
```

Output:

```text
noisy_linear.exr
noisy_raw.tiff or .npy
noisy_srgb.png
noise_aovs.exr, optional
```

Reason: sensor noise should be applied to final exposure/electron counts, not confused with Monte Carlo path tracing noise.

---

## Phase 3B later target: Mitsuba Film plugin

Only after Phase 3A works, port to a `Film` plugin:

```xml
<film type="cmosfilm">
    <integer name="width" value="1920"/>
    <integer name="height" value="1080"/>
    <string name="sensor_profile" value="resources/sensors/example_cmos.yaml"/>
    <string name="response_curves" value="resources/sensors/example_rgb_srf.csv"/>
    <boolean name="shot_noise" value="true"/>
    <boolean name="read_noise" value="true"/>
    <boolean name="fixed_pattern_noise" value="true"/>
    <boolean name="adc_quantization" value="true"/>
</film>
```

---

## Phase 3 sensor profile schema

Example:

```yaml
sensor:
  width: 1920
  height: 1080
  pixel_pitch_um: 4.0
  cfa_pattern: none        # none, RGGB, BGGR, GRBG, GBRG
  bit_depth: 12

exposure:
  exposure_time_s: 0.01
  analog_gain: 1.0
  digital_gain: 1.0
  black_level_dn: 64
  white_level_dn: 4095

charge:
  full_well_e: 30000
  conversion_gain_uV_per_e: 20.0
  dark_current_e_per_s: 0.1
  dark_current_temp_ref_C: 20.0
  dark_current_temp_doubling_C: 6.0

noise:
  enable_shot_noise: true
  enable_read_noise: true
  read_noise_e: 2.0
  enable_prnu: true
  prnu_sigma: 0.005
  enable_dsnu: true
  dsnu_sigma_e: 1.0
  enable_column_fpn: true
  column_fpn_sigma_dn: 0.5
  hot_pixel_fraction: 0.0001
  hot_pixel_dark_current_multiplier: 100.0
  seed: 42

adc:
  bits: 12
  enable_quantization: true
  adc_offset_dn: 0.0
  adc_gain_dn_per_e: auto
  enable_adc_nonlinearity: false

blooming:
  enabled: false
  threshold_e: 30000
  spill_fraction: 0.5
  mode: vertical       # vertical, isotropic
```

---

## Phase 3 effects to implement

Implement in this order.

### 1. Exposure scaling and electrons

```text
linear_sensor_value -> photons/electrons
```

Use a simple exposure scale first:

```text
electrons = exposure_scale * linear_value
```

Later, if spectral radiometry is calibrated:

```text
photons(lambda) = irradiance(lambda) * pixel_area * exposure_time * lambda / (h c)
electrons = photons * QE(lambda)
```

### 2. Full well and saturation

```text
electrons_clipped = min(electrons, full_well_e)
```

### 3. Photon shot noise

```text
electrons_noisy ~ Poisson(electrons)
```

### 4. Dark current

```text
dark_e = dark_current_e_per_s * exposure_time_s
```

Add dark current shot noise:

```text
dark_noisy ~ Poisson(dark_e)
```

### 5. Read noise

```text
electrons += Normal(0, read_noise_e)
```

### 6. PRNU

Fixed multiplicative per-pixel gain:

```text
electrons *= 1 + prnu_map
```

The PRNU map must be fixed for a given seed.

### 7. DSNU / dark FPN

Fixed additive per-pixel offset in electrons:

```text
electrons += dsnu_map
```

### 8. Column/row FPN

Add fixed column/row offsets:

```text
digital += column_offset[x]
digital += row_offset[y]
```

### 9. Hot pixels

Generate a fixed mask. Hot pixels have elevated dark current or fixed offsets.

### 10. ADC quantization

```text
dn = round(gain * electrons + black_level)
dn = clamp(dn, 0, 2^bits - 1)
```

### 11. Optional blooming

Not in the CCD/CMOS paper’s minimal model; implement only as optional extension.

Simple algorithm:

```text
excess = max(electrons - full_well_e, 0)
electrons = min(electrons, full_well_e)
spill excess to neighbors according to blooming kernel
```

CCD-like mode can spill mostly along columns.

---

## Phase 3 tests

### Numeric tests

```text
test_poisson_mean_variance
test_read_noise_variance
test_full_well_clipping
test_adc_quantization_range
test_seed_reproducibility
test_prnu_is_fixed_across_frames
test_hot_pixel_mask_is_fixed
test_column_fpn_shape
```

### Visual tests

1. **Flat field**
   - Expect PRNU texture and column banding if enabled.

2. **Dark frame**
   - Expect dark current, read noise, hot pixels.

3. **Exposure ramp**
   - Expect saturation/clipping at high exposure.

4. **Photon transfer curve**
   - Render flat fields at increasing exposure.
   - Plot variance vs mean.
   - Shot-noise-dominated region should have variance proportional to mean.

5. **High ISO vs low ISO**
   - Same image, different gain/noise settings.

---

## Phase 3 Definition of Done

Phase 3 is done when:

- The CMOS development path can consume a Mitsuba EXR/specfilm output.
- At least these effects work:
  - shot noise;
  - read noise;
  - full-well saturation;
  - quantization;
  - PRNU;
  - dark current;
  - hot pixels or column FPN.
- Noise is deterministic for a fixed seed.
- Sensor noise is visually separable from path-tracing noise.
- Documentation explains units and limitations.

---

# Validation plan across all phases

## Render comparison grid

For every major milestone, render the same scenes with:

```text
Mitsuba perspective
Mitsuba thinlens
RealisticCamera Phase 1 constant eta
RealisticCamera Phase 2 spectral glass
RealisticCamera Phase 2 spectral + optional ghosting
RealisticCamera + CMOS development
```

## Required figures for final report

1. Lens layout diagram.
2. Exit pupil visualization.
3. Thin lens vs realistic lens grid distortion.
4. Vignetting flat-field image.
5. Bokeh/cat’s-eye bokeh image.
6. Constant eta vs spectral CA image.
7. Spot diagram for R/G/B wavelengths.
8. Optional ghost-only AOV.
9. Clean render vs CMOS noisy render.
10. Dark frame / flat field noise examples.
11. Photon transfer curve, if Phase 3 is completed.

---

# Known limitations to document

## Phase 1 limitations

- Spherical elements only.
- Constant refractive index only.
- No coatings.
- No internal reflections.
- No diffraction.
- No polarization.
- No sensor noise.
- May not support all Mitsuba sensor importance-sampling methods.

## Phase 2 limitations

- Spectral glass depends on catalog quality.
- If using sampled wavelengths, CA noise may require high spp.
- Flare/ghosting, if implemented, is physically inspired but not a full unbiased lens light transport simulation.
- Aperture diffraction, if implemented as postprocess, is hybrid and approximate.

## Phase 3 limitations

- Sensor model is high-level, not transistor-level.
- Real camera ISP is not implemented.
- Bayer/demosaic may be omitted or approximate.
- Noise parameters are synthetic unless calibrated from real sensor data.
- Blooming is optional and not part of the referenced CCD/CMOS base model.

---

# Codex execution checklist

Use this checklist literally.

## First task batch

1. Inspect Mitsuba sensor plugins.
2. Create `realisticcamera` plugin skeleton.
3. Add XML parsing for Phase 1 parameters.
4. Render pinhole fallback.
5. Commit.

## Second task batch

1. Implement `.dat` parser.
2. Add parser tests.
3. Add toy lens data.
4. Commit.

## Third task batch

1. Implement spherical intersection.
2. Implement Snell refraction.
3. Add unit tests.
4. Commit.

## Fourth task batch

1. Implement trace from film through lens.
2. Add debug counters.
3. Render simple biconvex lens scene.
4. Commit.

## Fifth task batch

1. Implement sampling and weighting.
2. Add exit pupil sampling if simple sampling works.
3. Render vignetting and bokeh scenes.
4. Commit.

## Sixth task batch

1. Implement focus adjustment.
2. Render focus sweep.
3. Write Phase 1 docs.
4. Tag Phase 1 complete.

## Seventh task batch

1. Implement glass catalog.
2. Implement spectral YAML lens file.
3. Add Sellmeier IOR.
4. Add glass tests.
5. Commit.

## Eighth task batch

1. Connect sampled wavelength to lens tracing.
2. Render CA tests.
3. Add spot diagram tool.
4. Tag Phase 2A complete.

## Ninth task batch, optional

1. Implement two-reflection ghosting in separate code path.
2. Add ghost AOV.
3. Render bright-light ghost scene.
4. Document approximation.
5. Tag Phase 2B optional complete.

## Tenth task batch

1. Implement offline `cmos_develop.py`.
2. Add sensor profile YAML.
3. Add shot/read/saturation/quantization.
4. Add deterministic tests.
5. Commit.

## Eleventh task batch

1. Add PRNU/DSNU/hot pixels/column FPN.
2. Add dark frame and flat field examples.
3. Add photon transfer curve script.
4. Tag Phase 3A complete.

## Twelfth task batch, optional

1. Port CMOS development to Mitsuba `Film` plugin.
2. Keep offline tool as reference implementation.
3. Verify both produce matching output for same seed/profile.
4. Tag Phase 3B complete.

---

# Do-not-do list

Do not:

- implement chromatic aberration as RGB screen-space channel shifts;
- rewrite the whole renderer;
- add flare before basic lens tracing works;
- add CMOS noise per path sample and call it sensor noise;
- mix sensor blooming with optical lens flare;
- silently change units between mm and scene units;
- require aspheric lenses for Phase 1 or Phase 2;
- make Phase 2 depend on all 2011 paper effects;
- make Phase 3 depend on full Film plugin integration if a postprocess tool works.

---

# Minimal final project fallback

If time runs out, the minimum successful project is:

```text
Phase 1 complete
+ Phase 2A chromatic aberration complete
+ Phase 3A simple CMOS postprocess with shot/read/saturation/quantization
```

Skip:

```text
flare
diffraction
Bayer/demosaic
full Mitsuba Film plugin
blooming
advanced ADC nonlinearity
```

This fallback is still a coherent project:

> A Mitsuba compound-lens realistic camera with spectral glass dispersion and a physically motivated CMOS development/noise pipeline.
