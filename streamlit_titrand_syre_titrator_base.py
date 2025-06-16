import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from io import BytesIO

st.set_page_config(layout="centered")
st.markdown("<h3 style='text-align: center;'> 🔬 JKT's titrerkurvegenerator</h3>", unsafe_allow_html=True)

# Inputs
with st.sidebar:
    st.markdown("<p style='font-size:13px;'>Volumen af analyseopløsningen (mL)</p>", unsafe_allow_html=True)
    vanalytic = st.number_input("", value=25.0, step=1.0, min_value=10.0, max_value=30.0)

    st.markdown("<p style='font-size:13px;'>Koncentration af syren (M)</p>", unsafe_allow_html=True)
    cacid = st.number_input("", value=0.06, step=0.01, min_value=0.04, max_value=0.10)

    st.markdown("<p style='font-size:13px;'>Koncentration af basen (M)</p>", unsafe_allow_html=True)
    ctitr = st.number_input("", value=0.10, step=0.01, min_value=0.05, max_value=0.12)

    st.markdown("<p style='font-size:13px;'>pKa af syren</p>", unsafe_allow_html=True)
    pKa = st.number_input("", value=5.00, step=0.1, min_value=2.0, max_value=6.5)

#vanalytic = st.sidebar.number_input("Volumen af analyseopløsningen (mL)", value=25.0, step=1.0, min_value=10.0, max_value=30.0)
#cacid = st.sidebar.number_input("Koncentration af syren (M)", value=0.06, step=0.01, min_value=0.05, max_value=0.10)
#ctitr = st.sidebar.number_input("Koncentration af basen (M)", value=0.1, step=0.01, min_value=0.05, max_value=0.12)
#pKa = st.sidebar.number_input("pKa af syren", value=5.00, step=0.1, min_value=2.0, max_value=6.50)

# --- LICENSINFO I SIDEBAR ---
st.sidebar.markdown("---", unsafe_allow_html=True)
st.sidebar.markdown(
    "<p style='font-size:12px; color:gray;'>© 2025 Jens Kaalby Thomsen.<br>"
    "Dette værktøj er udgivet under <a href='https://opensource.org/licenses/MIT' target='_blank'>MIT-licens</a>.</p>",
    unsafe_allow_html=True
)
st.sidebar.markdown(
    "<p style='font-size:12px; color:gray;'>Hvis figurer fra denne applikation anvendes i opgaver eller lignende, "
    "bedes du kreditere forfatteren.</p>",
    unsafe_allow_html=True
)

Ka = 10 ** (-pKa)
Ve = cacid * vanalytic / ctitr
max_volume = 30  # Fikseret beregnings- og plotområde til 30 mL

def pH_buffer_min(pKa):
    if pKa <= 3:
        return pKa
    elif pKa >= 6:
        return pKa - 2.2
    elif 5.5 <= pKa < 6:
        return pKa - 1.8
    elif 5 <= pKa < 5.5:
        return pKa - 1.6
    elif 4.5 <= pKa < 5:
        return pKa - 1.3 
    elif 4 <= pKa < 5:
        return pKa - 1.1
    else:
        return pKa - (pKa - 3) * (1.5 / 3)

# Pufferregion
pH_buf_start = pH_buffer_min(pKa)
pH_buffer = np.linspace(pH_buf_start, pKa + 2.5, 200)
v_buffer = (cacid * 10**pH_buffer * vanalytic) / (ctitr * (10**pKa + 10**pH_buffer))

# Begræns buffer til max_volume
v_buffer = v_buffer[v_buffer <= max_volume]
pH_buffer = pH_buffer[:len(v_buffer)]

# Initial pH (løses numerisk)
def solve_initial_pH():
    def f(H):
        return H**2 / (cacid - H) - Ka
    H_init = brentq(f, 1e-15, cacid - 1e-15)
    return -np.log10(H_init)

pH_initial = solve_initial_pH()

# Blød transition fra 0 til starten af pufferområdet
v_transition = np.linspace(0, v_buffer[0], 20)
pH_transition = np.linspace(pH_initial, pH_buffer[0], 20)

# Post-ækvivalens område
v_post_start = max(Ve + 0.01, v_buffer[-1])
v_post = np.linspace(v_post_start, max_volume, 100)
OH_conc = (ctitr * v_post - cacid * vanalytic) / (v_post + vanalytic)
pH_post = 14 + np.log10(OH_conc)

# Samler alle datapunkter
v_combined = np.concatenate([v_transition, v_buffer, v_post])
pH_combined = np.concatenate([pH_transition, pH_buffer, pH_post])

# Plot
fig, ax = plt.subplots(figsize=(10, 5.5))
ax.plot(v_combined, pH_combined, color='darkred')

ax.set_xlabel("Volumen af NaOH (mL)", fontsize=14)
ax.set_ylabel("pH", fontsize=14)
ax.set_xlim(0, max_volume)
ax.set_ylim(0, 14)
ax.set_xticks(np.arange(0, max_volume + 1, 1))
ax.set_yticks(np.arange(0, 15, 1))
ax.tick_params(labelsize=10)
ax.set_title("Titrand: svag syre. Titrator: stærk base (NaOH)", fontsize=14)
ax.grid(True, linestyle='-', alpha=0.5)

st.pyplot(fig)

# Save figure to in-memory buffer
def save_figure(fig, fmt):
    buf = BytesIO()
    fig.savefig(buf, format=fmt, bbox_inches='tight')
    buf.seek(0)
    return buf

# Create download buttons side by side
col1, col2 = st.columns([1, 1])
with col1:
    st.download_button(
        label="Download PNG",
        data=save_figure(fig, "png"),
        file_name="titrerkurve.png",
        mime="image/png",
        key="download_png"
    )

with col2:
    st.download_button(
        label="Download SVG",
        data=save_figure(fig, "svg"),
        file_name="titrerkurve.svg",
        mime="image/svg+xml",
        key="download_svg"
    )
