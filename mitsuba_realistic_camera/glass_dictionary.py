import drjit as dr
import mitsuba as mi


# -----------------------------------------------------------------------------
# Glass database
# -----------------------------------------------------------------------------
# All wavelength-dependent entries assume lambda_um is in micrometers.
#
# Supported models:
#   constant:
#       n(lambda) = value
#
#   sellmeier:
#       n^2(lambda) = 1 + sum_i B_i * lambda^2 / (lambda^2 - C_i)
#
#   poly_n2:
#       n^2(lambda) =
#           a0
#         + a2   * lambda^2
#         + a4   * lambda^4
#         + am2  * lambda^-2
#         + am4  * lambda^-4
#         + am6  * lambda^-6
#         + am8  * lambda^-8
#
#   sampled:
#       piecewise-linear interpolation over sampled wavelength / ior pairs.
#       This is useful if you later want to paste PBRT's raw sampled glass arrays.
# -----------------------------------------------------------------------------

glass_dictionary = {
    # -------------------------------------------------------------------------
    # Constant media
    # -------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------
    # SCHOTT / common optical glasses, Sellmeier form
    # -------------------------------------------------------------------------
    "N-BK7": {
        "model": "sellmeier",
        "B": [1.03961212, 0.231792344, 1.01046945],
        "C": [0.00600069867, 0.0200179144, 103.560653],
    },
    "BK7": {
        "alias": "N-BK7",
    },
    "glass-BK7": {
        "alias": "N-BK7",
    },
    "glassbk7": {
        "alias": "N-BK7",
    },

    "N-BAF10": {
        "model": "sellmeier",
        "B": [1.5851495, 0.143559385, 1.08521269],
        "C": [0.00926681282, 0.0424489805, 105.613573],
    },
    "BAF10": {
        "alias": "N-BAF10",
    },
    "glass-BAF10": {
        "alias": "N-BAF10",
    },
    "glassbaf10": {
        "alias": "N-BAF10",
    },

    "N-FK51A": {
        "model": "sellmeier",
        "B": [0.971247817, 0.216901417, 0.904651666],
        "C": [0.00472301995, 0.0153575612, 168.68133],
    },
    "N-FK51": {
        "alias": "N-FK51A",
    },
    "FK51A": {
        "alias": "N-FK51A",
    },
    "glass-FK51A": {
        "alias": "N-FK51A",
    },
    "glassfk51a": {
        "alias": "N-FK51A",
    },

    "N-LASF9": {
        "model": "sellmeier",
        "B": [2.00029547, 0.298926886, 1.80691843],
        "C": [0.0121426017, 0.0538736236, 156.530829],
    },
    "LASF9": {
        "alias": "N-LASF9",
    },
    "glass-LASF9": {
        "alias": "N-LASF9",
    },
    "glasslasf9": {
        "alias": "N-LASF9",
    },

    "N-SF5": {
        "model": "sellmeier",
        "B": [1.52481889, 0.187085527, 1.42729015],
        "C": [0.011254756, 0.0588995392, 129.141675],
    },
    "SF5": {
        "alias": "N-SF5",
    },
    "glass-SF5": {
        "alias": "N-SF5",
    },
    "glasssf5": {
        "alias": "N-SF5",
    },

    "N-SF10": {
        "model": "sellmeier",
        "B": [1.62153902, 0.256287842, 1.64447552],
        "C": [0.0122241457, 0.0595736775, 147.468793],
    },
    "SF10": {
        # PBRT's "glass-SF10" is commonly mapped to SCHOTT N-SF10.
        "alias": "N-SF10",
    },
    "glass-SF10": {
        "alias": "N-SF10",
    },
    "glasssf10": {
        "alias": "N-SF10",
    },

    "N-SF11": {
        "model": "sellmeier",
        "B": [1.73759695, 0.313747346, 1.89878101],
        "C": [0.013188707, 0.0623068142, 155.23629],
    },
    "SF11": {
        "alias": "N-SF11",
    },
    "glass-SF11": {
        "alias": "N-SF11",
    },
    "glasssf11": {
        "alias": "N-SF11",
    },

    # -------------------------------------------------------------------------
    # Extra entries useful for the Brendel/Tessar example
    # -------------------------------------------------------------------------
    "LLF1": {
        "model": "sellmeier",
        "B": [1.21640125, 0.13366454, 0.883399468],
        "C": [0.00857807248, 0.0420143003, 107.59306],
    },
    "LLF7": {
        # Historical LLF7 is not always easy to find as a modern SCHOTT entry.
        # LLF1 has very close nd/Vd and is a useful substitute.
        "alias": "LLF1",
    },

    # This is HIKARI/Nikon SF7, not SCHOTT three-term Sellmeier.
    # nd ~= 1.639800, Vd ~= 34.56.
    "SF7": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 2.602245,
            "a2": -0.002282728,
            "a4": 0.0,
            "am2": 0.0313239,
            "am4": -0.001109815,
            "am6": 0.0002844724,
            "am8": -1.205447e-05,
        },
    },
    "E-SF7": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 2.60224509,
            "a2": -0.00228272809,
            "a4": 0.0,
            "am2": 0.0313239043,
            "am4": -0.00110981516,
            "am6": 0.00028447243,
            "am8": -1.2054467e-05,
        },
    },

    # -------------------------------------------------------------------------
    # HOYA/OHARA glasses used by the JP2022-033487 ZMX example
    # -------------------------------------------------------------------------
    "FC5": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 2.1894054,
            "a2": -0.0099044908,
            "a4": 0.0,
            "am2": 0.008640337,
            "am4": 0.00022263067,
            "am6": -1.2291942e-05,
            "am8": 5.9386349e-07,
        },
    },
    "FCD100": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 2.0482157,
            "a2": -0.0043211222,
            "a4": 0.0,
            "am2": 0.0061826755,
            "am4": 3.141148e-05,
            "am6": 3.5370793e-06,
            "am8": -1.6694497e-07,
        },
    },
    "TAFD33": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 3.4316955,
            "a2": -0.013156273,
            "a4": 0.0,
            "am2": 0.035631862,
            "am4": 0.00091378109,
            "am6": -4.0366886e-06,
            "am8": 3.3867688e-06,
        },
    },
    "MC-TAF105": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 3.0705638,
            "a2": -0.013310234,
            "a4": 0.0,
            "am2": 0.025280355,
            "am4": 0.00021829234,
            "am6": 3.2555293e-05,
            "am8": -8.1300421e-07,
        },
    },
    "FF5": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 2.4743324,
            "a2": -0.010955338,
            "a4": 0.0,
            "am2": 0.019293801,
            "am4": 0.0014497732,
            "am6": -0.00011038744,
            "am8": 1.1136008e-05,
        },
    },
    "FCD10A": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 2.1085888,
            "a2": -0.0046553722,
            "a4": 0.0,
            "am2": 0.006954748,
            "am4": 3.2184141e-05,
            "am6": 5.2670946e-06,
            "am8": -2.7164844e-07,
        },
    },
    "FCD705": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 2.3742832,
            "a2": -0.0062014418,
            "a4": 0.0,
            "am2": 0.010582967,
            "am4": 6.3100703e-05,
            "am6": 8.196715e-06,
            "am8": -3.4424446e-07,
        },
    },
    "TAFD45": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 3.672962,
            "a2": -0.0151805,
            "a4": 0.0,
            "am2": 0.0457187,
            "am4": 0.002186221,
            "am6": -9.93742e-05,
            "am8": 1.490466e-05,
        },
    },
    "NBFD25": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 3.2879389,
            "a2": -0.015856356,
            "a4": 0.0,
            "am2": 0.045638807,
            "am4": 0.0033027661,
            "am6": -0.00021690686,
            "am8": 2.9625863e-05,
        },
    },
    "NBFD29": {
        "model": "poly_n2",
        "coefficients": {
            "a0": 3.023484,
            "a2": -0.014446333,
            "a4": 0.0,
            "am2": 0.034673303,
            "am4": 0.0020637437,
            "am6": -0.00011126222,
            "am8": 1.4468973e-05,
        },
    },
    "S-NBM51": {
        "model": "sellmeier",
        "B": [1.37023101, 0.177665568, 1.30515471],
        "C": [0.00871920342, 0.0405725552, 112.703058],
    },

    # -------------------------------------------------------------------------
    # Existing / useful aliases from your previous dictionary
    # -------------------------------------------------------------------------
    "N-LAK9": {
        "model": "sellmeier",
        "B": [1.46231905, 0.344399589, 1.15508372],
        "C": [0.00724270181, 0.0243353131, 85.4686868],
    },
    "LAK9": {
        # For older lens prescriptions that simply say LAK9.
        # This is a practical modern mapping, not necessarily the exact
        # historical glass used in the original patent.
        "alias": "N-LAK9",
    },

    "N-SF6": {
        "model": "sellmeier",
        "B": [1.77931763, 0.338149866, 2.08734474],
        "C": [0.0133714182, 0.0617533621, 174.01759],
    },
    "SF6": {
        "alias": "N-SF6",
    },

    "N-F2": {
        "model": "sellmeier",
        "B": [1.39757037, 0.159201403, 1.2686543],
        "C": [0.00995906143, 0.0546931752, 119.248346],
    },
    "F2": {
        "model": "sellmeier",
        "B": [1.34533359, 0.209073176, 0.937357162],
        "C": [0.00997743871, 0.0470450767, 111.886764],
    },
    "N-BAK2": {
        "model": "sellmeier",
        "B": [1.01662154, 0.319903051, 0.937232995],
        "C": [0.00592383763, 0.0203828415, 113.118417],
    },

    "N-LAF34": {
        "model": "sellmeier",
        "B": [1.75836958, 0.313537785, 1.18925231],
        "C": [0.00872810026, 0.0293020832, 85.1780644],
    },

    "N-LAK8": {
        "model": "sellmeier",
        "B": [1.33183167, 0.546623206, 1.19084015],
        "C": [0.00620023871, 0.0216465439, 82.5827736],
    },
}


def _resolve_glass_entry(name, max_depth=16):
    """Resolve aliases in glass_dictionary."""
    if name not in glass_dictionary:
        raise KeyError(f"Unknown glass name: {name}")

    entry = glass_dictionary[name]
    depth = 0

    while "alias" in entry:
        name = entry["alias"]
        if name not in glass_dictionary:
            raise KeyError(f"Glass alias points to unknown entry: {name}")

        entry = glass_dictionary[name]
        depth += 1

        if depth > max_depth:
            raise RuntimeError(f"Glass alias chain seems cyclic or too deep: {name}")

    return entry


def _as_tensor(x):
    return dr.auto.TensorXf(x)


def sellmeier_ior(lambda_um, B, C):
    """
    Three-term Sellmeier IOR.

    lambda_um:
        Wavelength in micrometers.
    B, C:
        Sellmeier coefficients. C is in um^2.
    """
    if len(B) != len(C):
        raise ValueError(f"B and C must have the same length, got {len(B)} and {len(C)}")

    # to micrometers
    lam = _as_tensor(lambda_um) / 1000.0
    lam2 = dr.sqr(lam)

    eta2 = dr.ones_like(lam2)

    for b, c in zip(B, C):
        eta2 += b * lam2 / (lam2 - c)

    return dr.sqrt(eta2)


def poly_n2_ior(lambda_um, coefficients):
    """
    Polynomial-in-lambda model for n^2(lambda).

    Supports terms:
        a0
        a2  * lambda^2
        a4  * lambda^4
        am2 * lambda^-2
        am4 * lambda^-4
        am6 * lambda^-6
        am8 * lambda^-8
    """
    # to micrometers
    lam = _as_tensor(lambda_um) / 1000.0
    lam2 = dr.sqr(lam)

    inv_lam2 = 1.0 / lam2

    a0 = coefficients.get("a0", 0.0)
    a2 = coefficients.get("a2", 0.0)
    a4 = coefficients.get("a4", 0.0)
    am2 = coefficients.get("am2", 0.0)
    am4 = coefficients.get("am4", 0.0)
    am6 = coefficients.get("am6", 0.0)
    am8 = coefficients.get("am8", 0.0)

    eta2 = (
        a0
        + a2 * lam2
        + a4 * dr.sqr(lam2)
        + am2 * inv_lam2
        + am4 * dr.sqr(inv_lam2)
        + am6 * dr.sqr(inv_lam2) * inv_lam2
        + am8 * dr.sqr(dr.sqr(inv_lam2))
    )

    return dr.sqrt(eta2)


def sampled_ior(lambda_um, wavelengths_nm, values):
    """
    Piecewise-linear sampled IOR.

    This is meant to mimic PBRT's PiecewiseLinearSpectrum-style sampled data.
    Input lambda is still lambda_um; table wavelengths are in nm by default.
    Values outside the tabulated domain are clamped to the endpoint values.
    """
    if len(wavelengths_nm) != len(values):
        raise ValueError(
            f"wavelengths_nm and values must have same length, got "
            f"{len(wavelengths_nm)} and {len(values)}"
        )
    if len(wavelengths_nm) < 2:
        raise ValueError("sampled_ior requires at least two samples.")

    lam_nm = _as_tensor(lambda_um) * 1000.0

    # Clamp by default. This is often better than returning zero for IOR.
    result = dr.full_like(lam_nm, float(values[0]))
    result = dr.select(lam_nm >= wavelengths_nm[-1], float(values[-1]), result)

    for i in range(len(wavelengths_nm) - 1):
        x0 = float(wavelengths_nm[i])
        x1 = float(wavelengths_nm[i + 1])
        y0 = float(values[i])
        y1 = float(values[i + 1])

        t = (lam_nm - x0) / (x1 - x0)
        y = y0 + t * (y1 - y0)

        mask = (lam_nm >= x0) & (lam_nm < x1)
        result = dr.select(mask, y, result)

    return result


def eval_glass_ior(glass_name, lambda_um):
    entry = _resolve_glass_entry(glass_name)
    model = entry["model"]

    if model == "sellmeier":
        return sellmeier_ior(lambda_um, entry["B"], entry["C"])

    if model == "poly_n2":
        return poly_n2_ior(lambda_um, entry["coefficients"])

    if model == "sampled":
        return sampled_ior(lambda_um, entry["wavelengths_nm"], entry["values"])

    raise ValueError(f"Unsupported glass model: {model}")


class Eta_lookup:
    def __init__(self, prop):
        self.type = prop["type"]

        if self.type == "glass":
            self.glass_name = prop["glass_name"]
            # Resolve once here, mostly to fail early if the name is invalid.
            _resolve_glass_entry(self.glass_name)

        elif self.type == "constant":
            self.value = float(prop["value"])

        else:
            raise ValueError(f"Unsupported eta type: {self.type}")

    def __call__(self, lambda_nm):
        if lambda_nm is None:
            lambda_nm = 550  # Default to green if no wavelength provided.
        if self.type == "glass":
            return eval_glass_ior(self.glass_name, lambda_nm)

        if self.type == "constant":
            lam = _as_tensor(lambda_nm)
            # return dr.full_like(lam, self.value)
            return dr.full(dr.auto.TensorXf, self.value, lam.shape)

        raise RuntimeError("Unreachable eta type.")

if __name__ == "__main__":
    mi.set_variant('cuda_ad_spectral')
    # Quick test: print IOR at 550nm for all entries.
    import numpy as np
    spec_np = np.ones((4, 50), dtype=np.float32) * 550
    spec = mi.Spectrum(spec_np)
            
    # glass_types = glass_dictionary.keys()
    glass_types = ["SF7"]
    for glass_name in glass_types:
        eta_i = Eta_lookup({"type": "glass", "glass_name": glass_name})
        ior = eta_i(spec)
        print(f"{glass_name}: IOR at 550nm = {ior[0]}")
