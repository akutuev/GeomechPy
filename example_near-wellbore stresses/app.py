"""
Near Wellbore Stresses — Calculation & Visualisation (Streamlit front end).

A one-page tool built directly on top of `geomechpy.near_wellbore_stresses`.
Upload well data (CSV / Excel / LAS) or define the inputs by hand, and the app
computes the stresses on the borehole wall around the full 0°–360° sweep for
any borehole orientation, then lets you step through the well with a depth
player and watch the profile change.

Every stress number shown comes from the library:

  * `NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses`
    -> σrr, σθθ, σzz, στz, σrθ, σrz
  * `NearWellboreStressesCalculation.calculate_principal_stresses_analytical`
    -> σ1, σ2 and the wellbore tortuosity angle

No stress equation is reimplemented here; the app only prepares the inputs,
converts units and draws the results.

Run locally:   streamlit run app.py
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from utils import near_wellbore as nw

# ---------------------------------------------------------------------------
# Page setup & session state
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="Near Wellbore Stresses",
    page_icon="🕳️",
    layout="wide",
    initial_sidebar_state="expanded",
)

for key, default in {
    "raw_df": None,          # uploaded / sample / manual input data
    "params": None,          # per-depth inputs handed to geomechpy (psi, deg)
    "profile": None,         # per-depth summary of the 0-360 sweep
    "config": None,          # configuration used for the last run
    "data_source": None,     # label shown in the sidebar
    "load_messages": [],     # info messages from the data cleaner
    "input_mode": "manual",  # 'manual' or 'file'
    "playing": False,        # depth player auto-advance
}.items():
    st.session_state.setdefault(key, default)


def _theme_base() -> str:
    """'light' or 'dark' for the active Streamlit theme (best effort)."""
    try:
        return str(st.context.theme.type or "light")
    except Exception:
        try:
            return str(st.get_option("theme.base") or "light")
        except Exception:
            return "light"


COLORS = nw.component_colors(_theme_base())
GRID_COLOR = "rgba(128,128,128,0.25)"
REF_COLOR = "rgba(128,128,128,0.55)"


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _annotation_offset(theta_at: float) -> int:
    """Nudge a θ annotation inward so it clears the axis at 0° and 360°."""
    if theta_at < 45:
        return 45
    if theta_at > 315:
        return -45
    return 0


def _radial_range(values: np.ndarray) -> list[float]:
    """Radial axis bounds with a little head room; keeps negative stresses visible."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return [0.0, 1.0]
    lo, hi = float(np.min(finite)), float(np.max(finite))
    pad = max((hi - lo) * 0.08, abs(hi) * 0.02, 1.0)
    return [lo - pad, hi + pad]


def polar_figure(result: pd.DataFrame, components: list[str], unit: str, depth_label: str) -> go.Figure:
    """Stress distribution around the borehole wall, θ = 0° at top of hole."""
    fig = go.Figure()
    stacked = []
    for comp in components:
        values = result[comp].to_numpy(dtype=float)
        stacked.append(values)
        fig.add_trace(go.Scatterpolar(
            r=values,
            theta=result["THETA"],
            mode="lines",
            name=nw.COMPONENT_LABELS[comp],
            line=dict(color=COLORS[comp], width=2),
            hovertemplate=f"{nw.DISPLAY_NAMES[comp]} %{{r:,.0f}} {unit}<br>θ %{{theta:.0f}}°<extra></extra>",
        ))

    rng = _radial_range(np.concatenate(stacked)) if stacked else [0.0, 1.0]
    fig.update_layout(
        title=f"Stress distribution around the borehole wall — {depth_label}",
        polar=dict(
            # angle=0 lays the radial axis out horizontally, the one orientation in
            # which Plotly draws its tick labels upright.
            radialaxis=dict(range=rng, gridcolor=GRID_COLOR, angle=0, tickangle=0, nticks=5,
                            tickformat=",.0f", ticksuffix=f" {unit}", tickfont=dict(size=10)),
            angularaxis=dict(direction="clockwise", rotation=90, dtick=30,
                             gridcolor=GRID_COLOR, ticksuffix="°"),
        ),
        legend=dict(orientation="h", yanchor="top", y=-0.05, x=0),
        margin=dict(l=40, r=40, t=60, b=40),
        height=560,
    )
    return fig


def _add_direction_markers(fig: go.Figure, directions: dict) -> None:
    """Dotted guides at the SHmax- and Shmin-facing wall points."""
    for theta_at, label in ((directions["theta_shmax"], "SHmax"), (directions["theta_shmin"], "Shmin")):
        for x in (theta_at, theta_at + 180.0):
            if 0 <= x <= 360:
                fig.add_vline(x=x, line=dict(color=REF_COLOR, width=1, dash="dot"),
                              annotation_text=label, annotation_position="top",
                              annotation_font=dict(size=10))


def components_figure(result: pd.DataFrame, components: list[str], unit: str, summary: dict,
                      directions: dict) -> go.Figure:
    """Stress components against the angle around the borehole wall."""
    fig = go.Figure()
    for comp in components:
        fig.add_trace(go.Scatter(
            x=result["THETA"], y=result[comp], mode="lines",
            name=nw.COMPONENT_LABELS[comp],
            line=dict(color=COLORS[comp], width=2),
            hovertemplate=f"{nw.DISPLAY_NAMES[comp]} %{{y:,.0f}} {unit}<extra></extra>",
        ))

    if "SIGMA_TT" in components:
        for theta_at, value, text in (
            (summary["theta_tt_max"], summary["sigma_tt_max"], "σθθ max"),
            (summary["theta_tt_min"], summary["sigma_tt_min"], "σθθ min"),
        ):
            fig.add_annotation(x=theta_at, y=value, text=f"{text} @ {theta_at:.0f}°",
                               showarrow=True, arrowhead=2, arrowsize=0.8,
                               arrowcolor=REF_COLOR, ax=_annotation_offset(theta_at), ay=-28,
                               font=dict(size=11))

    _add_direction_markers(fig, directions)
    fig.update_layout(
        title="Stress components vs. angle around the borehole wall",
        xaxis=dict(title="θ from top of hole (deg)", dtick=45, range=[0, 360], gridcolor=GRID_COLOR),
        yaxis=dict(title=f"Stress ({unit})", gridcolor=GRID_COLOR, tickformat=",.0f"),
        hovermode="x unified",
        legend=dict(orientation="h", yanchor="top", y=-0.18, x=0),
        margin=dict(l=60, r=30, t=60, b=60),
        height=460,
    )
    return fig


def principal_figure(result: pd.DataFrame, unit: str, summary: dict, directions: dict) -> go.Figure:
    """Maximum and minimum principal stresses on the borehole wall."""
    fig = go.Figure()
    for comp in nw.PRINCIPAL_COMPONENTS:
        fig.add_trace(go.Scatter(
            x=result["THETA"], y=result[comp], mode="lines",
            name=nw.COMPONENT_LABELS[comp],
            line=dict(color=COLORS[comp], width=2),
            hovertemplate=f"{nw.DISPLAY_NAMES[comp]} %{{y:,.0f}} {unit}<extra></extra>",
        ))
    fig.add_hline(y=0, line=dict(color=REF_COLOR, width=1, dash="dot"))
    fig.add_annotation(x=summary["theta_s1_max"], y=summary["sigma_1_max"],
                       text=f"σ1 max @ {summary['theta_s1_max']:.0f}°", showarrow=True,
                       arrowhead=2, arrowsize=0.8, arrowcolor=REF_COLOR,
                       ax=_annotation_offset(summary["theta_s1_max"]), ay=-28, font=dict(size=11))
    fig.add_annotation(x=summary["theta_s2_min"], y=summary["sigma_2_min"],
                       text=f"σ2 min @ {summary['theta_s2_min']:.0f}°", showarrow=True,
                       arrowhead=2, arrowsize=0.8, arrowcolor=REF_COLOR,
                       ax=_annotation_offset(summary["theta_s2_min"]), ay=28, font=dict(size=11))
    _add_direction_markers(fig, directions)
    fig.update_layout(
        title="Principal stresses at the borehole wall",
        xaxis=dict(title="θ from top of hole (deg)", dtick=45, range=[0, 360], gridcolor=GRID_COLOR),
        yaxis=dict(title=f"Stress ({unit})", gridcolor=GRID_COLOR, tickformat=",.0f"),
        hovermode="x unified",
        legend=dict(orientation="h", yanchor="top", y=-0.18, x=0),
        margin=dict(l=60, r=30, t=60, b=60),
        height=420,
    )
    return fig


def tortuosity_figure(result: pd.DataFrame) -> go.Figure:
    """Wellbore tortuosity angle — a single series, so the title names it."""
    fig = go.Figure(go.Scatter(
        x=result["THETA"], y=result["THETA_TORTUOSITY"], mode="lines",
        line=dict(color=COLORS["SIGMA_1"], width=2), showlegend=False,
        hovertemplate="θ tortuosity %{y:.1f}°<extra></extra>",
    ))
    fig.add_hline(y=0, line=dict(color=REF_COLOR, width=1, dash="dot"))
    fig.update_layout(
        title="Wellbore tortuosity angle (plane of zero shear stress)",
        xaxis=dict(title="θ from top of hole (deg)", dtick=45, range=[0, 360], gridcolor=GRID_COLOR),
        yaxis=dict(title="θ tortuosity (deg)", gridcolor=GRID_COLOR),
        hovermode="x unified",
        margin=dict(l=60, r=30, t=60, b=50),
        height=320,
    )
    return fig


def profile_figure(profile: pd.DataFrame, unit: str, depth_unit: str, selected_depth: float) -> go.Figure:
    """Borehole-wall extremes against depth, with the player position marked."""
    fig = go.Figure()
    # σ1 ≥ max(σθθ, σzz) everywhere, so σ1 max would draw straight over σθθ max;
    # only the four curves that carry distinct information are drawn.
    series = [
        ("SIGMA_1_MAX", "σ1 max", COLORS["SIGMA_1"], "solid"),
        ("SIGMA_TT_MIN", "σθθ min", COLORS["SIGMA_TT"], "dash"),
        ("SIGMA_RR", "σrr (Pw − Pp)", COLORS["SIGMA_RR"], "solid"),
        ("SIGMA_2_MIN", "σ2 min", COLORS["SIGMA_2"], "solid"),
    ]
    for col, label, color, dash in series:
        fig.add_trace(go.Scatter(
            x=profile[col], y=profile["DEPTH"], mode="lines", name=label,
            line=dict(color=color, width=2, dash=dash),
            hovertemplate=f"{label} %{{x:,.0f}} {unit}<extra></extra>",
        ))
    fig.add_hline(y=selected_depth, line=dict(color=REF_COLOR, width=1.5, dash="dot"),
                  annotation_text=f"selected {selected_depth:,.1f} {depth_unit}",
                  annotation_position="top right")
    fig.update_layout(
        title="Borehole-wall stress extremes vs. depth",
        xaxis=dict(title=f"Stress ({unit})", gridcolor=GRID_COLOR, tickformat=",.0f"),
        yaxis=dict(title=f"Depth ({depth_unit})", autorange="reversed", gridcolor=GRID_COLOR),
        hovermode="y unified",
        legend=dict(orientation="h", yanchor="top", y=-0.12, x=0),
        margin=dict(l=70, r=30, t=60, b=60),
        height=620,
    )
    return fig


# ---------------------------------------------------------------------------
# Cached geomechpy sweeps
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner=False)
def cached_profile(params: pd.DataFrame, theta_step: float) -> pd.DataFrame:
    return nw.depth_profile(params, nw.theta_grid(theta_step))


@st.cache_data(show_spinner=False)
def cached_all_depths(params: pd.DataFrame, theta_step: float) -> pd.DataFrame:
    return nw.all_depth_results(params, nw.theta_grid(theta_step))


# ---------------------------------------------------------------------------
# Sidebar — inputs & parameters
# ---------------------------------------------------------------------------

def source_selector(curve: str, mapped: bool, kind: str) -> str:
    """Radio choosing between a mapped log column and a manual value.

    `mapped` is part of the widget key so the default flips to the log column
    as soon as one is mapped.
    """
    manual_label = "Gradient" if kind == "gradient" else "Constant"
    labels = {"Log column": "column", manual_label: kind if kind == "gradient" else "constant"}
    choice = st.radio(
        f"{nw.CURVE_LABELS[curve]} from",
        list(labels.keys()),
        index=0 if mapped else 1,
        horizontal=True,
        key=f"src_{curve}_{int(mapped)}",
    )
    return labels[choice]


with st.sidebar:
    st.title("🕳️ Near Wellbore Stresses")
    st.caption(
        "Kirsch borehole-wall stresses powered by "
        "[geomechpy](https://github.com/sohwaisheng1/GeomechPy_WS) — "
        "`geomechpy.near_wellbore_stresses`."
    )

    with st.expander("1. Data input", expanded=True):
        uploaded = st.file_uploader(
            "Upload well data (CSV / Excel / LAS)",
            type=["csv", "txt", "xls", "xlsx", "las"],
            help="One row per depth sample. Any of DEPTH, SV, SHMAX, SHMIN, PP, PW, "
                 "SHMAX_AZI, PR, INC, AZI can be supplied; anything missing is filled "
                 "in from a gradient or a constant below.",
        )
        if uploaded is not None:
            try:
                st.session_state.raw_df, st.session_state.load_messages = nw.load_data(uploaded)
                st.session_state.data_source = f"📄 {uploaded.name}"
                st.session_state.input_mode = "file"
                st.session_state.params = None
                st.session_state.profile = None
            except ValueError as exc:
                st.error(str(exc))

        c_sample, c_dl = st.columns(2)
        if c_sample.button("🧪 Load Sample", width="stretch"):
            st.session_state.raw_df = nw.generate_sample_data()
            st.session_state.data_source = "🧪 Synthetic deviated well (6000–12000 ft TVD)"
            st.session_state.load_messages = []
            st.session_state.input_mode = "file"
            st.session_state.params = None
            st.session_state.profile = None
        c_dl.download_button(
            "⬇️ Template", data=nw.sample_csv_bytes(),
            file_name="near_wellbore_example.csv", mime="text/csv",
            width="stretch", help="Correctly formatted example CSV.",
        )

        st.markdown("**…or define a depth range manually**")
        m1, m2, m3 = st.columns(3)
        manual_top = m1.number_input("Top", value=6000.0, step=100.0, format="%.1f")
        manual_base = m2.number_input("Base", value=12000.0, step=100.0, format="%.1f")
        manual_step = m3.number_input("Step", value=50.0, min_value=0.1, step=10.0, format="%.1f")
        if st.button("📐 Use manual depth range", width="stretch"):
            try:
                st.session_state.raw_df = nw.manual_depth_frame(manual_top, manual_base, manual_step)
                st.session_state.data_source = (
                    f"📐 Manual depth range {manual_top:,.0f}–{manual_base:,.0f} (step {manual_step:g})"
                )
                st.session_state.load_messages = []
                st.session_state.input_mode = "manual"
                st.session_state.params = None
                st.session_state.profile = None
            except ValueError as exc:
                st.error(str(exc))

        if st.session_state.data_source:
            st.success(f"Loaded: {st.session_state.data_source}")
        for msg in st.session_state.load_messages:
            st.info(msg)

    raw_df = st.session_state.raw_df
    has_data = raw_df is not None
    columns = list(raw_df.columns) if has_data else []

    with st.expander("2. Column mapping", expanded=has_data and st.session_state.input_mode == "file"):
        column_map: dict[str, str] = {}
        if has_data:
            options = ["-- not mapped --"] + columns
            for curve in nw.ALL_CURVES:
                required = curve in nw.REQUIRED_CURVES
                label = f"{nw.CURVE_LABELS[curve]} {'(required)' if required else '(optional)'}"
                choice = st.selectbox(label, options, index=nw.guess_column(curve, columns), key=f"map_{curve}")
                column_map[curve] = "" if choice == "-- not mapped --" else choice
        else:
            st.info("Load data or define a depth range first.")

    mapped = {curve: bool(column_map.get(curve)) for curve in nw.ALL_CURVES}

    with st.expander("3. Units", expanded=False):
        depth_unit = st.radio("Depth (TVD) unit", nw.DEPTH_UNITS, index=0, horizontal=True,
                              help="Gradients are applied in psi/ft, so metric depths are converted first.")
        pressure_unit = st.selectbox("Input pressure unit", list(nw.PRESSURE_INPUT_UNITS.keys()), index=0,
                                     help="Applies to the mapped Sv / SHmax / Shmin / Pp / Pw columns.")
        output_unit = st.selectbox("Output pressure unit", list(nw.PRESSURE_OUTPUT_UNITS.keys()), index=0,
                                   help="Results are computed in psi; this only changes how they are displayed.")
        air_gap = st.number_input(f"Air gap / KB ({depth_unit})", value=0.0, min_value=0.0, format="%.1f",
                                  help="Added to the depth before a gradient is applied.")

    with st.expander("4. Far-field stresses & pore pressure", expanded=True):
        st.caption("Take each magnitude from a mapped log column, or build it from a gradient × TVD.")
        sv_source = source_selector("SV", mapped["SV"], "gradient")
        sv_gradient = st.number_input("Sv gradient (psi/ft)", value=1.00, min_value=0.30, max_value=2.00,
                                      step=0.005, format="%.3f", disabled=(sv_source == "column"))
        shmax_source = source_selector("SHMAX", mapped["SHMAX"], "gradient")
        shmax_gradient = st.number_input("SHmax gradient (psi/ft)", value=0.90, min_value=0.20, max_value=2.00,
                                         step=0.005, format="%.3f", disabled=(shmax_source == "column"))
        shmin_source = source_selector("SHMIN", mapped["SHMIN"], "gradient")
        shmin_gradient = st.number_input("Shmin gradient (psi/ft)", value=0.75, min_value=0.20, max_value=2.00,
                                         step=0.005, format="%.3f", disabled=(shmin_source == "column"))
        pp_source = source_selector("PP", mapped["PP"], "gradient")
        pp_gradient = st.number_input("Pp gradient (psi/ft)", value=0.465, min_value=0.20, max_value=1.20,
                                      step=0.005, format="%.3f", disabled=(pp_source == "column"))

    with st.expander("5. Wellbore pressure", expanded=True):
        pw_mode = st.radio(
            "Wellbore pressure Pw from",
            ["Mud weight", "Gradient", "Log column"],
            index=2 if mapped["PW"] else 0, horizontal=True,
            key=f"pw_mode_{int(mapped['PW'])}",
            help="Mud weight and gradient are converted to a pressure at TVD before the geomechpy call.",
        )
        pw_source = {"Mud weight": "mud_weight", "Gradient": "gradient", "Log column": "column"}[pw_mode]
        mw_col, mwu_col = st.columns([2, 1])
        mud_weight = mw_col.number_input("Mud weight", value=9.6, min_value=1.0, max_value=25.0, step=0.1,
                                         format="%.2f", disabled=(pw_source != "mud_weight"))
        mud_weight_unit = mwu_col.selectbox("Unit", list(nw.MUD_WEIGHT_UNITS.keys()), index=0,
                                            disabled=(pw_source != "mud_weight"))
        pw_gradient = st.number_input("Pw gradient (psi/ft)", value=0.50, min_value=0.10, max_value=2.00,
                                      step=0.005, format="%.3f", disabled=(pw_source != "gradient"))
        if pw_source == "column" and not mapped["PW"]:
            st.warning("No Pw column is mapped — map one in section 2 or switch to mud weight / gradient.")

    with st.expander("6. Wellbore geometry & rock", expanded=True):
        shmax_azi_source = source_selector("SHMAX_AZI", mapped["SHMAX_AZI"], "constant")
        shmax_azimuth = st.slider("SHmax azimuth (deg from North)", 0.0, 360.0, 135.0, 1.0,
                                  disabled=(shmax_azi_source == "column"))
        inc_source = source_selector("INC", mapped["INC"], "constant")
        inclination = st.slider("Borehole inclination (deg)", 0.0, 90.0, 30.0, 1.0,
                                disabled=(inc_source == "column"),
                                help="0° = vertical, 90° = horizontal.")
        azi_source = source_selector("AZI", mapped["AZI"], "constant")
        azimuth = st.slider("Borehole azimuth (deg from North)", 0.0, 360.0, 45.0, 1.0,
                            disabled=(azi_source == "column"))
        pr_source = source_selector("PR", mapped["PR"], "constant")
        poisson_ratio = st.slider("Static Poisson's ratio", 0.05, 0.49, 0.25, 0.01,
                                  disabled=(pr_source == "column"))

    with st.expander("7. Angular resolution", expanded=False):
        theta_step = st.select_slider("θ step around the borehole wall (deg)",
                                      options=[0.5, 1.0, 2.0, 5.0, 10.0, 15.0], value=2.0,
                                      help="The 0°–360° sweep passed to geomechpy as the `theta` array.")

    st.divider()
    run_clicked = st.button("🚀 Run Calculation", type="primary", width="stretch", disabled=not has_data)

# Bundle the configuration once so the depth player can reuse it without a re-run.
config = dict(
    depth_unit=depth_unit, pressure_unit=pressure_unit, output_unit=output_unit, air_gap=air_gap,
    sv_source=sv_source, sv_gradient=sv_gradient,
    shmax_source=shmax_source, shmax_gradient=shmax_gradient,
    shmin_source=shmin_source, shmin_gradient=shmin_gradient,
    pp_source=pp_source, pp_gradient=pp_gradient,
    pw_source=pw_source, pw_gradient=pw_gradient,
    mud_weight=mud_weight, mud_weight_unit=mud_weight_unit,
    shmax_azi_source=shmax_azi_source, shmax_azimuth=shmax_azimuth,
    inc_source=inc_source, inclination=inclination,
    azi_source=azi_source, azimuth=azimuth,
    pr_source=pr_source, poisson_ratio=poisson_ratio,
    theta_step=theta_step,
) if has_data else nw.default_config()

# ---------------------------------------------------------------------------
# Run the calculation
# ---------------------------------------------------------------------------

if run_clicked:
    try:
        with st.spinner("Computing near wellbore stresses with geomechpy..."):
            params = nw.build_parameter_table(raw_df, column_map, config)
            profile = cached_profile(params, float(theta_step))
        st.session_state.params = params
        st.session_state.profile = profile
        st.session_state.config = config
        st.session_state.playing = False
        st.toast("Calculation complete ✅")
    except ValueError as exc:
        st.session_state.params = None
        st.session_state.profile = None
        st.error(f"⚠️ {exc}")
    except Exception as exc:  # keep the app alive on unexpected input
        st.session_state.params = None
        st.session_state.profile = None
        st.error(f"Unexpected error during calculation: {exc}")

# ---------------------------------------------------------------------------
# Main area
# ---------------------------------------------------------------------------

st.title("Near Wellbore Stresses — Calculation & Visualisation")
st.markdown(
    "Borehole-wall stresses for **any borehole orientation**, computed with the Kirsch "
    "solution in `geomechpy.near_wellbore_stresses`. Configure the inputs in the sidebar, "
    "press **🚀 Run Calculation**, then use the **depth player** to sweep the well."
)

params = st.session_state.params
profile = st.session_state.profile
run_config = st.session_state.config or config
out_unit = run_config.get("output_unit", "psi")
d_unit = run_config.get("depth_unit", "ft")

if params is None:
    st.info(
        "No results yet. Load the sample data (or upload your own / define a depth range) "
        "in the sidebar, then press **🚀 Run Calculation**."
    )
    ref = pd.DataFrame({
        "geomechpy argument": ["shmin", "shmax", "svert", "pore_pressure", "mud_pressure",
                               "shmax_azimuth", "borehole_deviation", "borehole_azimuth",
                               "poisson_ratio_static", "theta"],
        "App input": ["Shmin", "SHmax", "Sv", "Pore pressure Pp", "Wellbore pressure Pw",
                      "SHmax azimuth", "Borehole inclination", "Borehole azimuth",
                      "Static Poisson's ratio", "0°–360° sweep at the chosen θ step"],
        "Unit": ["psi", "psi", "psi", "psi", "psi", "deg", "deg", "deg", "–", "deg"],
    })
    st.subheader("Inputs required by `calculate_kirsch_borehole_wall_stresses`")
    st.dataframe(ref, width="stretch", hide_index=True)
    if st.session_state.raw_df is not None:
        st.subheader("Input data preview")
        st.dataframe(st.session_state.raw_df, width="stretch", height=320)
    st.stop()

# --- Depth player -----------------------------------------------------------

depths = params["DEPTH"].to_numpy(dtype=float)
n_depths = len(depths)
depth_options = [float(d) for d in depths]

# `depth_idx` is the authoritative player position. The slider is a mirror of
# it: Streamlit forbids writing a widget's key after the widget exists, so the
# index is pushed into the slider *before* it is created, and the slider's
# on_change pushes the user's drag back into the index.
st.session_state.setdefault("depth_idx", 0)
st.session_state.depth_idx = int(np.clip(st.session_state.depth_idx, 0, n_depths - 1))


def _step_to(index: int) -> None:
    st.session_state.depth_idx = int(np.clip(index, 0, n_depths - 1))


def _sync_from_slider() -> None:
    value = float(st.session_state.depth_slider)
    st.session_state.depth_idx = int(np.argmin(np.abs(depths - value)))


st.subheader("🎚️ Depth player")
st.caption(f"{n_depths} depth samples — move the slider or press play to watch the profile evolve.")

ctrl = st.columns([1, 1, 1, 1, 1, 3])
ctrl[0].button("⏮", width="stretch", help="First depth",
               on_click=_step_to, args=(0,), disabled=st.session_state.playing)
ctrl[1].button("◀", width="stretch", help="Previous depth",
               on_click=_step_to, args=(st.session_state.depth_idx - 1,),
               disabled=st.session_state.playing)
if st.session_state.playing:
    ctrl[2].button("⏸", width="stretch", type="primary", help="Pause",
                   on_click=lambda: st.session_state.update(playing=False))
else:
    ctrl[2].button("▶", width="stretch", type="primary", help="Play through the well",
                   on_click=lambda: st.session_state.update(playing=True))
ctrl[3].button("▶|", width="stretch", help="Next depth",
               on_click=_step_to, args=(st.session_state.depth_idx + 1,),
               disabled=st.session_state.playing)
ctrl[4].button("⏭", width="stretch", help="Last depth",
               on_click=_step_to, args=(n_depths - 1,), disabled=st.session_state.playing)
play_speed = ctrl[5].select_slider("Play speed (s per depth)", options=[0.1, 0.25, 0.5, 1.0, 2.0],
                                   value=0.25, key="play_speed")

st.session_state.depth_slider = depth_options[st.session_state.depth_idx]
st.select_slider(
    f"Depth ({d_unit})", options=depth_options, key="depth_slider",
    on_change=_sync_from_slider, format_func=lambda d: f"{d:,.1f} {d_unit}",
)
st.checkbox("Loop playback", value=True, key="loop_play")

idx = st.session_state.depth_idx
row = params.iloc[idx]
depth_label = f"{row['DEPTH']:,.1f} {d_unit}"
st.progress((idx + 1) / n_depths, text=f"Depth {idx + 1} of {n_depths} — {depth_label}")

# --- Stresses at the selected depth (geomechpy) -----------------------------

theta = nw.theta_grid(run_config.get("theta_step", 2.0))
result_psi = nw.wall_stresses_at_depth(row, theta)
summary_psi = nw.wall_summary(result_psi)

result = nw.convert_pressures(result_psi, nw.RESULT_PRESSURE_COLUMNS, out_unit)
summary = nw.wall_summary(result)
directions = nw.stress_direction_thetas(row)

st.divider()
m = st.columns(5)
m[0].metric(f"σrr (Pw − Pp) [{out_unit}]", f"{summary['sigma_rr']:,.0f}")
m[0].caption("constant around the wall")
for col, label, value, theta_at in (
    (m[1], "σθθ max", summary["sigma_tt_max"], summary["theta_tt_max"]),
    (m[2], "σθθ min", summary["sigma_tt_min"], summary["theta_tt_min"]),
    (m[3], "σ1 max", summary["sigma_1_max"], summary["theta_s1_max"]),
    (m[4], "σ2 min", summary["sigma_2_min"], summary["theta_s2_min"]),
):
    col.metric(f"{label} [{out_unit}]", f"{value:,.0f}")
    col.caption(f"at θ = {theta_at:.0f}°")

st.caption(
    f"θ = 0° is the top of hole, in the vertical plane bearing {row['AZI']:.0f}° from North. "
    f"SHmax bears {row['SHMAX_AZI']:.0f}°, so the SHmax-facing wall sits at "
    f"θ ≈ {directions['theta_shmax']:.0f}° / {directions['theta_shmax'] + 180:.0f}° and the "
    f"Shmin-facing wall at θ ≈ {directions['theta_shmin']:.0f}° / "
    f"{directions['theta_shmin'] + 180:.0f}°. The Kirsch solution puts σθθ at its **minimum** "
    f"facing SHmax and its **maximum** facing Shmin"
    + ("." if directions["exact"] else
       f" — exact for a vertical hole; at {row['INC']:.0f}° inclination the wall circle is "
       f"tilted, so treat these θ values as a guide.")
)

with st.expander("Inputs handed to geomechpy at this depth"):
    inputs = pd.DataFrame({
        "Argument": ["shmin", "shmax", "svert", "pore_pressure", "mud_pressure",
                     "shmax_azimuth", "borehole_deviation", "borehole_azimuth", "poisson_ratio_static"],
        "Value (psi / deg)": [row["SHMIN"], row["SHMAX"], row["SV"], row["PP"], row["PW"],
                              row["SHMAX_AZI"], row["INC"], row["AZI"], row["PR"]],
    })
    st.dataframe(inputs.style.format({"Value (psi / deg)": "{:,.4g}"}),
                 width="stretch", hide_index=True)

warnings = nw.parameter_warnings(params)
for note in warnings:
    st.warning(note)

# --- Tabs -------------------------------------------------------------------

tab_polar, tab_cart, tab_principal, tab_depth, tab_data = st.tabs(
    ["🧭 Polar view", "📈 Components vs θ", "🎯 Principal stresses", "🪜 Depth profile", "📋 Results & download"]
)

with tab_polar:
    chosen = st.multiselect(
        "Components to draw", nw.PLOTTED_COMPONENTS, default=nw.PLOTTED_COMPONENTS,
        format_func=lambda c: nw.COMPONENT_LABELS[c], key="polar_components",
    ) or ["SIGMA_TT"]
    st.plotly_chart(polar_figure(result, chosen, out_unit, depth_label), width="stretch")
    st.caption(
        "θ is measured from the top of hole (TOH) and increases clockwise looking down the "
        "borehole axis. The radial axis starts at the minimum plotted stress, not at zero, so "
        "negative (tensile) values stay visible — read the tick labels. σθθ usually dwarfs the "
        "other components; de-select it above to inspect σrr and στz on their own scale. "
        "σθθ is smallest on the SHmax-facing wall and largest on the Shmin-facing wall — see "
        "the orientation note above the tabs for where those fall on the θ axis."
    )

with tab_cart:
    st.plotly_chart(components_figure(result, nw.PLOTTED_COMPONENTS, out_unit, summary, directions),
                    width="stretch")
    st.caption(
        "σrθ and σrz are identically zero on the borehole wall (traction-free surface), so they "
        "are reported in the results table rather than drawn as flat lines."
    )

with tab_principal:
    st.plotly_chart(principal_figure(result, out_unit, summary, directions), width="stretch")
    st.plotly_chart(tortuosity_figure(result), width="stretch")
    st.caption(
        "σ1, σ2 and the tortuosity angle come from "
        "`calculate_principal_stresses_analytical(σθθ, σzz, στz)`. The tortuosity angle is "
        "undefined where σθθ = σzz."
    )

with tab_depth:
    profile_disp = nw.convert_pressures(
        profile,
        ["SIGMA_RR", "SIGMA_TT_MAX", "SIGMA_TT_MIN", "SIGMA_ZZ_MAX", "SIGMA_ZZ_MIN",
         "SIGMA_1_MAX", "SIGMA_2_MIN"],
        out_unit,
    )
    st.plotly_chart(profile_figure(profile_disp, out_unit, d_unit, float(row["DEPTH"])),
                    width="stretch")
    st.caption(
        "Each point is the extreme of a full 0°–360° geomechpy sweep at that depth. σ1 max is "
        "the largest principal stress anywhere on the wall; it equals σθθ max wherever the "
        "tangential-axial shear vanishes, so σθθ max is left out to avoid two identical curves."
    )

with tab_data:
    st.subheader(f"Borehole-wall stresses at {depth_label}")
    st.dataframe(nw.display_frame(result[nw.RESULT_COLUMNS]).style.format(precision=2),
                 width="stretch", height=380)

    st.subheader("Downloads")
    d1, d2, d3 = st.columns(3)
    d1.download_button(
        "⬇️ Selected depth (CSV)", data=nw.to_csv_bytes(nw.display_frame(result[nw.RESULT_COLUMNS])),
        file_name=f"near_wellbore_stresses_{row['DEPTH']:.0f}{d_unit}.csv",
        mime="text/csv", width="stretch",
    )
    d2.download_button(
        "⬇️ Depth profile (CSV)", data=nw.to_csv_bytes(profile_disp),
        file_name="near_wellbore_depth_profile.csv", mime="text/csv", width="stretch",
    )
    if d3.button("📦 Build all-depth export", width="stretch",
                 help="Full θ sweep at every depth — may take a moment for fine θ steps."):
        st.session_state.all_depth_csv = nw.to_csv_bytes(
            nw.convert_pressures(
                cached_all_depths(params, float(run_config.get("theta_step", 2.0))),
                nw.RESULT_PRESSURE_COLUMNS, out_unit,
            ).rename(columns=nw.DISPLAY_NAMES)
        )
    if st.session_state.get("all_depth_csv"):
        d3.download_button(
            "⬇️ All depths (CSV)", data=st.session_state.all_depth_csv,
            file_name="near_wellbore_stresses_all_depths.csv", mime="text/csv",
            width="stretch",
        )

    with st.expander("Per-depth inputs used (psi / deg)"):
        st.dataframe(params.style.format(precision=2), width="stretch", height=320)

st.caption(
    "All stress values are produced by `geomechpy.near_wellbore_stresses."
    "NearWellboreStressesCalculation`; the app only prepares inputs, converts units and plots."
)

# ---------------------------------------------------------------------------
# Depth player auto-advance (kept last so the plots render before the re-run)
# ---------------------------------------------------------------------------

if st.session_state.playing:
    next_index = idx + 1
    if next_index >= n_depths:
        if st.session_state.get("loop_play", True):
            next_index = 0
        else:
            st.session_state.playing = False
            next_index = idx
    if st.session_state.playing:
        time.sleep(float(st.session_state.get("play_speed", 0.25)))
        st.session_state.depth_idx = next_index
    st.rerun()
