"""Helper layer for the Near Wellbore Stresses Streamlit app.

Everything in this module is input preparation, book-keeping and unit
handling. **No stress equation lives here.** The near wellbore stresses and
the principal stresses are produced exclusively by

    geomechpy.near_wellbore_stresses.NearWellboreStressesCalculation
        .calculate_kirsch_borehole_wall_stresses(...)
        .calculate_principal_stresses_analytical(...)

which in turn use the rotations in `geomechpy.toolbox`. The results returned
by those functions are passed through untouched.
"""

from __future__ import annotations

import io

import numpy as np
import pandas as pd

from geomechpy.near_wellbore_stresses import NearWellboreStressesCalculation

# ---------------------------------------------------------------------------
# Unit constants
# ---------------------------------------------------------------------------

MPA_TO_PSI = 145.0377377             # MPa   -> psi
BAR_TO_PSI = 14.50377377             # bar   -> psi
KPA_TO_PSI = 0.1450377377            # kPa   -> psi
M_PER_FT = 0.3048                    # metres per foot
FT_PER_M = 1.0 / M_PER_FT

# Mud weight -> pressure gradient. Plain unit conversions applied to the input
# before it reaches geomechpy; they are not part of any stress solution.
PPG_TO_PSI_PER_FT = 0.0519479        # 1 lb/gal  = 0.0519479 psi/ft
SG_TO_PSI_PER_FT = 0.4335275         # 1 s.g.    = 0.4335275 psi/ft

# Input pressure unit -> psi
PRESSURE_INPUT_UNITS = {"psi": 1.0, "MPa": MPA_TO_PSI, "bar": BAR_TO_PSI, "kPa": KPA_TO_PSI}
# psi -> output display unit
PRESSURE_OUTPUT_UNITS = {"psi": 1.0, "MPa": 1.0 / MPA_TO_PSI, "bar": 1.0 / BAR_TO_PSI, "kPa": 1.0 / KPA_TO_PSI}

DEPTH_UNITS = ["ft", "m"]
MUD_WEIGHT_UNITS = {"ppg": PPG_TO_PSI_PER_FT, "s.g.": SG_TO_PSI_PER_FT}

NULL_SENTINELS = [-999.0, -999.25, -9999.0, -9999.25, -99999.0, 9999.0]

# ---------------------------------------------------------------------------
# Curves the app understands
# ---------------------------------------------------------------------------

REQUIRED_CURVES = ["DEPTH"]
OPTIONAL_CURVES = ["SV", "SHMAX", "SHMIN", "PP", "PW", "SHMAX_AZI", "PR", "INC", "AZI"]
ALL_CURVES = REQUIRED_CURVES + OPTIONAL_CURVES

CURVE_LABELS = {
    "DEPTH": "Depth (TVD)",
    "SV": "Vertical stress Sv",
    "SHMAX": "Max. horizontal stress SHmax",
    "SHMIN": "Min. horizontal stress Shmin",
    "PP": "Pore pressure Pp",
    "PW": "Wellbore / mud pressure Pw",
    "SHMAX_AZI": "SHmax azimuth (deg)",
    "PR": "Static Poisson's ratio",
    "INC": "Borehole inclination (deg)",
    "AZI": "Borehole azimuth (deg)",
}

# Curves carrying a pressure; everything else is an angle or unitless.
PRESSURE_CURVES = ["SV", "SHMAX", "SHMIN", "PP", "PW"]

# Parameter table columns handed to geomechpy (pressures already in psi).
PARAM_COLUMNS = [
    "DEPTH", "TVD_FT", "SV", "SHMAX", "SHMIN", "PP", "PW",
    "SHMAX_AZI", "PR", "INC", "AZI",
]

# ---------------------------------------------------------------------------
# Stress component book-keeping
# ---------------------------------------------------------------------------

# Fixed categorical slots: a component keeps its colour in every chart.
COMPONENT_COLORS_LIGHT = {
    "SIGMA_RR": "#2a78d6",  # blue
    "SIGMA_TT": "#eb6834",  # orange
    "SIGMA_ZZ": "#1baf7a",  # aqua
    "SIGMA_TZ": "#eda100",  # yellow
    "SIGMA_RT": "#e87ba4",  # magenta
    "SIGMA_RZ": "#008300",  # green
    "SIGMA_1": "#4a3aa7",   # violet
    "SIGMA_2": "#e34948",   # red
}
COMPONENT_COLORS_DARK = {
    "SIGMA_RR": "#3987e5",
    "SIGMA_TT": "#d95926",
    "SIGMA_ZZ": "#199e70",
    "SIGMA_TZ": "#c98500",
    "SIGMA_RT": "#d55181",
    "SIGMA_RZ": "#008300",
    "SIGMA_1": "#9085e9",
    "SIGMA_2": "#e66767",
}

COMPONENT_LABELS = {
    "SIGMA_RR": "σrr — radial",
    "SIGMA_TT": "σθθ — tangential",
    "SIGMA_ZZ": "σzz — axial",
    "SIGMA_TZ": "στz — tangential-axial shear",
    "SIGMA_RT": "σrθ — radial-tangential shear",
    "SIGMA_RZ": "σrz — radial-axial shear",
    "SIGMA_1": "σ1 — max. principal",
    "SIGMA_2": "σ2 — min. principal",
}

WALL_COMPONENTS = ["SIGMA_RR", "SIGMA_TT", "SIGMA_ZZ", "SIGMA_TZ", "SIGMA_RT", "SIGMA_RZ"]
PRINCIPAL_COMPONENTS = ["SIGMA_1", "SIGMA_2"]
# σrθ and σrz are identically zero on the borehole wall (free-surface boundary
# condition), so they are reported in the table rather than drawn as flat lines.
PLOTTED_COMPONENTS = ["SIGMA_RR", "SIGMA_TT", "SIGMA_ZZ", "SIGMA_TZ"]

RESULT_COLUMNS = ["THETA"] + WALL_COMPONENTS + PRINCIPAL_COMPONENTS + ["THETA_TORTUOSITY"]
# Columns holding a pressure -> converted together when the display unit changes.
RESULT_PRESSURE_COLUMNS = WALL_COMPONENTS + PRINCIPAL_COMPONENTS

DISPLAY_NAMES = {
    "THETA": "θ (deg)",
    "SIGMA_RR": "σrr",
    "SIGMA_TT": "σθθ",
    "SIGMA_ZZ": "σzz",
    "SIGMA_TZ": "στz",
    "SIGMA_RT": "σrθ",
    "SIGMA_RZ": "σrz",
    "SIGMA_1": "σ1",
    "SIGMA_2": "σ2",
    "THETA_TORTUOSITY": "θ tortuosity (deg)",
}


def component_colors(theme_base: str = "light") -> dict[str, str]:
    """Fixed component -> colour map, stepped for the active Streamlit theme."""
    return COMPONENT_COLORS_DARK if str(theme_base).lower() == "dark" else COMPONENT_COLORS_LIGHT


# ---------------------------------------------------------------------------
# Data loading & cleaning
# ---------------------------------------------------------------------------

def _drop_unit_rows(df: pd.DataFrame, max_rows: int = 3) -> tuple[pd.DataFrame, int]:
    """Drop leading non-numeric rows (unit headers such as 'ft', 'psi')."""
    dropped = 0
    while dropped < max_rows and len(df) > 0:
        first = df.iloc[0]
        numeric = pd.to_numeric(first, errors="coerce")
        if numeric.notna().sum() >= max(1, int(0.5 * len(first))):
            break
        df = df.iloc[1:]
        dropped += 1
    return df.reset_index(drop=True), dropped


def clean_dataframe(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Coerce to numeric, drop unit rows and null sentinels, report what happened."""
    messages: list[str] = []

    df.columns = [str(c).strip() for c in df.columns]
    df, dropped = _drop_unit_rows(df)
    if dropped:
        messages.append(f"Dropped {dropped} unit/header row(s) below the column names.")

    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(axis=1, how="all")
    replaced = int(df.isin(NULL_SENTINELS).sum().sum())
    if replaced:
        df = df.replace(NULL_SENTINELS, np.nan)
        messages.append(f"Replaced {replaced} null sentinel value(s) (-999.25 and similar) with NaN.")

    before = len(df)
    df = df.dropna(axis=0, how="all").reset_index(drop=True)
    if before != len(df):
        messages.append(f"Dropped {before - len(df)} empty row(s).")

    return df, messages


def load_data(uploaded_file) -> tuple[pd.DataFrame, list[str]]:
    """Read an uploaded CSV, Excel or LAS file and clean it up.

    Returns (dataframe, informational messages).
    Raises ValueError with a user-friendly message on failure.
    """
    name = uploaded_file.name.lower()
    try:
        if name.endswith(".las"):
            import lasio  # optional dependency; only needed for LAS input
            try:
                uploaded_file.seek(0)
            except Exception:
                pass
            raw = uploaded_file.getvalue() if hasattr(uploaded_file, "getvalue") else uploaded_file.read()
            text = raw.decode("utf-8", errors="replace") if isinstance(raw, (bytes, bytearray)) else raw
            las = lasio.read(io.StringIO(text))
            df = las.df().reset_index()  # depth curve is the index -> expose as a column
        elif name.endswith((".csv", ".txt")):
            df = pd.read_csv(uploaded_file, skip_blank_lines=True)
            if df.shape[1] == 1:  # non-comma delimiter: re-sniff
                uploaded_file.seek(0)
                df = pd.read_csv(uploaded_file, sep=None, engine="python", skip_blank_lines=True)
        elif name.endswith((".xls", ".xlsx")):
            df = pd.read_excel(uploaded_file)
        else:
            raise ValueError("Unsupported file type. Please upload a .las, .csv, .xls or .xlsx file.")
    except ValueError:
        raise
    except Exception as exc:
        raise ValueError(f"Could not parse file '{uploaded_file.name}': {exc}") from exc

    if df.empty:
        raise ValueError("The uploaded file contains no rows.")
    if df.shape[1] < 2:
        raise ValueError("The uploaded file needs at least a depth column and one more curve.")

    df, messages = clean_dataframe(df)
    if df.empty:
        raise ValueError("No data rows remained after cleaning the file.")
    return df, messages


def guess_column(curve: str, columns: list[str]) -> int:
    """Best-effort index of the column matching a curve mnemonic (selectbox default).

    Returns 0 ('-- not mapped --') when nothing matches.
    """
    aliases = {
        "DEPTH": ["depth", "dept", "tvd", "md"],
        "SV": ["sv", "svert", "ovb", "obg", "overburden", "sigv", "sigmav", "vertical"],
        "SHMAX": ["shmax", "sh_max", "shmaxi", "sigh", "sigmah", "smax"],
        "SHMIN": ["shmin", "sh_min", "sigh_min", "sigmah_min", "smin", "fg"],
        "PP": ["pp", "pore", "porepressure", "pore_pressure", "ppg", "pnorm"],
        "PW": ["pw", "pmud", "mud", "mw", "ecd", "wellbore"],
        "SHMAX_AZI": ["shmax_azi", "shmaxazi", "azi_shmax", "shazim", "shmax_azimuth", "shaz"],
        "PR": ["pr", "poisson", "nu", "poissons", "pois"],
        "INC": ["inc", "incl", "dev", "deviation", "inclination", "drift"],
        "AZI": ["azi", "azim", "azimuth", "hazi", "bearing"],
    }
    lowered = [c.lower().strip() for c in columns]
    for alias in aliases[curve]:
        for i, col in enumerate(lowered):
            if col == alias or col.startswith(alias):
                return i + 1  # +1 for the '-- not mapped --' placeholder at index 0
    return 0


# ---------------------------------------------------------------------------
# Default configuration & sample data
# ---------------------------------------------------------------------------

def default_config() -> dict:
    """Configuration used before the sidebar has been touched."""
    return dict(
        depth_unit="ft",
        pressure_unit="psi",
        output_unit="psi",
        sv_source="gradient", sv_gradient=1.00,
        shmax_source="gradient", shmax_gradient=0.90,
        shmin_source="gradient", shmin_gradient=0.75,
        pp_source="gradient", pp_gradient=0.465,
        pw_source="mud_weight", pw_gradient=0.50,
        mud_weight=9.6, mud_weight_unit="ppg",
        shmax_azimuth=0.0,
        poisson_ratio=0.25,
        inclination=0.0,
        azimuth=0.0,
        air_gap=0.0,
        theta_step=2.0,
    )


def generate_sample_data(top: float = 6000.0, base: float = 12000.0, step: float = 50.0) -> pd.DataFrame:
    """Synthetic deviated-well data set (TVD in ft, pressures in psi).

    Gradients drift mildly with depth so the near wellbore profile visibly
    changes as the depth player is stepped through the well.
    """
    tvd = np.arange(top, base + step, step, dtype=float)
    n = len(tvd)
    frac = (tvd - tvd[0]) / max(tvd[-1] - tvd[0], 1.0)

    sv_grad = 0.98 + 0.06 * frac
    pp_grad = 0.452 + 0.075 * frac ** 2          # mild overpressure ramp
    shmin_grad = 0.74 + 0.05 * frac
    shmax_grad = 0.86 + 0.13 * frac              # stress regime rotates with depth
    mud_grad = pp_grad + 0.045

    return pd.DataFrame({
        "DEPTH": tvd,
        "SV": sv_grad * tvd,
        "SHMAX": shmax_grad * tvd,
        "SHMIN": shmin_grad * tvd,
        "PP": pp_grad * tvd,
        "PW": mud_grad * tvd,
        "SHMAX_AZI": np.full(n, 135.0),
        "PR": 0.22 + 0.06 * frac,
        "INC": np.linspace(5.0, 65.0, n),   # builds angle with depth
        "AZI": np.full(n, 45.0),
    }).round(4)


def sample_csv_bytes() -> bytes:
    """Correctly formatted example CSV for the sidebar template download."""
    return generate_sample_data().to_csv(index=False).encode("utf-8")


# ---------------------------------------------------------------------------
# Parameter table
# ---------------------------------------------------------------------------

def _column_values(data: pd.DataFrame, column: str) -> np.ndarray | None:
    if not column or column not in data.columns:
        return None
    return pd.to_numeric(data[column], errors="coerce").to_numpy(dtype=float)


def _resolve_pressure(
    data: pd.DataFrame,
    column: str,
    source: str,
    unit_factor: float,
    gradient_psift: float,
    tvd_ft: np.ndarray,
    label: str,
) -> np.ndarray:
    """Pressure curve in psi, taken from a mapped column or from a gradient."""
    if source == "column":
        values = _column_values(data, column)
        if values is None:
            raise ValueError(f"{label}: source is 'column' but no column is mapped.")
        return values * unit_factor
    return gradient_psift * tvd_ft


def _resolve_scalar(
    data: pd.DataFrame,
    column: str,
    source: str,
    constant: float,
    n: int,
    label: str,
) -> np.ndarray:
    """Angle / Poisson's ratio curve, taken from a mapped column or a constant."""
    if source == "column":
        values = _column_values(data, column)
        if values is None:
            raise ValueError(f"{label}: source is 'column' but no column is mapped.")
        return values
    return np.full(n, float(constant))


def build_parameter_table(data: pd.DataFrame, column_map: dict[str, str], config: dict) -> pd.DataFrame:
    """Assemble the per-depth inputs required by `calculate_kirsch_borehole_wall_stresses`.

    Returns one row per depth with all pressures already converted to psi and
    all angles in degrees. Rows with a missing value are dropped, since the
    geomechpy call needs every argument.
    """
    if data is None or data.empty:
        raise ValueError("No input data available.")

    depth_col = column_map.get("DEPTH", "")
    depth = _column_values(data, depth_col)
    if depth is None:
        raise ValueError("A depth column must be mapped before the stresses can be computed.")

    depth_unit = config.get("depth_unit", "ft")
    tvd_ft = depth * FT_PER_M if depth_unit == "m" else depth.copy()
    tvd_ft = tvd_ft + (config.get("air_gap", 0.0) * (FT_PER_M if depth_unit == "m" else 1.0))

    unit_factor = PRESSURE_INPUT_UNITS[config.get("pressure_unit", "psi")]
    n = len(depth)

    # Wellbore pressure: a mapped column, a gradient, or a mud weight.
    pw_source = config.get("pw_source", "mud_weight")
    if pw_source == "mud_weight":
        mw_factor = MUD_WEIGHT_UNITS[config.get("mud_weight_unit", "ppg")]
        pw = float(config.get("mud_weight", 9.6)) * mw_factor * tvd_ft
    else:
        pw = _resolve_pressure(data, column_map.get("PW", ""), pw_source, unit_factor,
                               config.get("pw_gradient", 0.5), tvd_ft, "Wellbore pressure Pw")

    params = pd.DataFrame({
        "DEPTH": depth,
        "TVD_FT": tvd_ft,
        "SV": _resolve_pressure(data, column_map.get("SV", ""), config.get("sv_source", "gradient"),
                                unit_factor, config.get("sv_gradient", 1.0), tvd_ft, "Vertical stress Sv"),
        "SHMAX": _resolve_pressure(data, column_map.get("SHMAX", ""), config.get("shmax_source", "gradient"),
                                   unit_factor, config.get("shmax_gradient", 0.9), tvd_ft, "SHmax"),
        "SHMIN": _resolve_pressure(data, column_map.get("SHMIN", ""), config.get("shmin_source", "gradient"),
                                   unit_factor, config.get("shmin_gradient", 0.75), tvd_ft, "Shmin"),
        "PP": _resolve_pressure(data, column_map.get("PP", ""), config.get("pp_source", "gradient"),
                                unit_factor, config.get("pp_gradient", 0.465), tvd_ft, "Pore pressure Pp"),
        "PW": pw,
        "SHMAX_AZI": _resolve_scalar(data, column_map.get("SHMAX_AZI", ""), config.get("shmax_azi_source", "constant"),
                                     config.get("shmax_azimuth", 0.0), n, "SHmax azimuth"),
        "PR": _resolve_scalar(data, column_map.get("PR", ""), config.get("pr_source", "constant"),
                              config.get("poisson_ratio", 0.25), n, "Poisson's ratio"),
        "INC": _resolve_scalar(data, column_map.get("INC", ""), config.get("inc_source", "constant"),
                               config.get("inclination", 0.0), n, "Borehole inclination"),
        "AZI": _resolve_scalar(data, column_map.get("AZI", ""), config.get("azi_source", "constant"),
                               config.get("azimuth", 0.0), n, "Borehole azimuth"),
    })

    params = params.replace([np.inf, -np.inf], np.nan).dropna(subset=PARAM_COLUMNS)
    params = params.sort_values("DEPTH").reset_index(drop=True)
    if params.empty:
        raise ValueError("No depth sample has a complete set of inputs. Check the column mapping and units.")
    return params


def manual_depth_frame(top: float, base: float, step: float) -> pd.DataFrame:
    """Depth grid for the manual-input mode (no file uploaded)."""
    if step <= 0:
        raise ValueError("The depth step must be greater than zero.")
    if base <= top:
        raise ValueError("The base depth must be deeper than the top depth.")
    return pd.DataFrame({"DEPTH": np.arange(top, base + step * 0.5, step, dtype=float)})


def parameter_warnings(params: pd.DataFrame) -> list[str]:
    """Physical sanity checks on the assembled inputs (advisory only)."""
    notes: list[str] = []
    if (params["SHMAX"] < params["SHMIN"]).any():
        notes.append("SHmax is below Shmin at one or more depths — check the stress inputs.")
    if (params["PP"] > params["SHMIN"]).any():
        notes.append("Pore pressure exceeds Shmin at one or more depths.")
    if (params["PW"] < params["PP"]).any():
        notes.append("Wellbore pressure is below pore pressure at one or more depths (underbalanced).")
    if ((params["PR"] <= 0) | (params["PR"] >= 0.5)).any():
        notes.append("Poisson's ratio falls outside 0–0.5 at one or more depths.")
    return notes


# ---------------------------------------------------------------------------
# geomechpy calls
# ---------------------------------------------------------------------------

def theta_grid(step_deg: float = 2.0) -> np.ndarray:
    """Azimuthal angles around the borehole wall, 0°–360° inclusive."""
    step_deg = float(step_deg)
    if step_deg <= 0:
        raise ValueError("The θ step must be greater than zero.")
    n = int(round(360.0 / step_deg)) + 1
    return np.linspace(0.0, 360.0, n)


def wall_stresses_at_depth(row, theta: np.ndarray) -> pd.DataFrame:
    """Near wellbore stresses at a single depth, straight from geomechpy.

    Calls `calculate_kirsch_borehole_wall_stresses` for the six wall stress
    components and `calculate_principal_stresses_analytical` for σ1, σ2 and the
    tortuosity angle. No value is modified after the call.
    """
    wall = NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses(
        shmin=float(row["SHMIN"]),
        shmax=float(row["SHMAX"]),
        svert=float(row["SV"]),
        pore_pressure=float(row["PP"]),
        shmax_azimuth=float(row["SHMAX_AZI"]),
        mud_pressure=float(row["PW"]),
        theta=theta,
        poisson_ratio_static=float(row["PR"]),
        borehole_deviation=float(row["INC"]),
        borehole_azimuth=float(row["AZI"]),
    )

    # σ1/σ2 use arctan(2·στz / (σθθ − σzz)); the denominator can vanish, which
    # numpy reports as a warning. The library's result is kept as returned.
    with np.errstate(divide="ignore", invalid="ignore"):
        principal = NearWellboreStressesCalculation.calculate_principal_stresses_analytical(
            sigma_tt=wall.sigma_tt,
            sigma_zz=wall.sigma_zz,
            sigma_tz=wall.sigma_tz,
        )

    return pd.DataFrame({
        "THETA": theta,
        "SIGMA_RR": wall.sigma_rr,
        "SIGMA_TT": wall.sigma_tt,
        "SIGMA_ZZ": wall.sigma_zz,
        "SIGMA_TZ": wall.sigma_tz,
        "SIGMA_RT": wall.sigma_rt,
        "SIGMA_RZ": wall.sigma_rz,
        "SIGMA_1": principal.sigma_1,
        "SIGMA_2": principal.sigma_2,
        "THETA_TORTUOSITY": principal.theta_tortuosity,
    })


def depth_profile(params: pd.DataFrame, theta: np.ndarray) -> pd.DataFrame:
    """Per-depth summary of the geomechpy output over the full 0°–360° sweep."""
    records = []
    for _, row in params.iterrows():
        res = wall_stresses_at_depth(row, theta)
        tt = res["SIGMA_TT"].to_numpy()
        zz = res["SIGMA_ZZ"].to_numpy()
        s1 = res["SIGMA_1"].to_numpy()
        s2 = res["SIGMA_2"].to_numpy()
        records.append({
            "DEPTH": row["DEPTH"],
            "SIGMA_RR": res["SIGMA_RR"].iloc[0],           # constant around the wall
            "SIGMA_TT_MAX": np.nanmax(tt),
            "SIGMA_TT_MIN": np.nanmin(tt),
            "THETA_AT_TT_MAX": theta[int(np.nanargmax(tt))],
            "THETA_AT_TT_MIN": theta[int(np.nanargmin(tt))],
            "SIGMA_ZZ_MAX": np.nanmax(zz),
            "SIGMA_ZZ_MIN": np.nanmin(zz),
            "SIGMA_1_MAX": np.nanmax(s1),
            "SIGMA_2_MIN": np.nanmin(s2),
        })
    return pd.DataFrame(records)


def all_depth_results(params: pd.DataFrame, theta: np.ndarray) -> pd.DataFrame:
    """Long-format θ sweep for every depth (used by the full-well download)."""
    frames = []
    for _, row in params.iterrows():
        res = wall_stresses_at_depth(row, theta)
        res.insert(0, "DEPTH", row["DEPTH"])
        frames.append(res)
    return pd.concat(frames, ignore_index=True)


def wall_summary(result: pd.DataFrame) -> dict:
    """Headline numbers at the selected depth, read off the geomechpy output."""
    theta = result["THETA"].to_numpy()
    tt = result["SIGMA_TT"].to_numpy()
    s1 = result["SIGMA_1"].to_numpy()
    s2 = result["SIGMA_2"].to_numpy()
    return {
        "sigma_rr": float(result["SIGMA_RR"].iloc[0]),
        "sigma_tt_max": float(np.nanmax(tt)),
        "sigma_tt_min": float(np.nanmin(tt)),
        "theta_tt_max": float(theta[int(np.nanargmax(tt))]),
        "theta_tt_min": float(theta[int(np.nanargmin(tt))]),
        "sigma_1_max": float(np.nanmax(s1)),
        "theta_s1_max": float(theta[int(np.nanargmax(s1))]),
        "sigma_2_min": float(np.nanmin(s2)),
        "theta_s2_min": float(theta[int(np.nanargmin(s2))]),
    }


# ---------------------------------------------------------------------------
# Presentation helpers
# ---------------------------------------------------------------------------

def convert_pressures(df: pd.DataFrame, columns: list[str], output_unit: str) -> pd.DataFrame:
    """Copy of `df` with the listed psi columns expressed in `output_unit`."""
    factor = PRESSURE_OUTPUT_UNITS.get(output_unit, 1.0)
    out = df.copy()
    if factor != 1.0:
        for col in columns:
            if col in out.columns:
                out[col] = out[col] * factor
    return out


def display_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Rename internal columns to their symbols for on-screen tables."""
    return df.rename(columns=DISPLAY_NAMES)


def to_csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")
