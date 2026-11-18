"""
createMockCatalogue_sigma.py

Build a GW mock catalogue from a given model rate file (.npz with z, R, DL),
deriving the distance-error relation sigma(D_L) along the way, with no prior
knowledge of sigma.

Stage 1 - calibration
    Draw N_CAL events from the rate, give each masses and angles from the same
    priors as the original generator, and compute its distance error with the
    bilby Fisher matrix. Fit one shared relation  sigma(D) = a * D^b  to these.

Stage 2 - mock
    Draw the catalogue events from the same rate (independent random draws) and
    add noise with the fitted relation evaluated at each event's TRUE distance.
    The likelihood uses the same relation at each trial distance, so mock and
    likelihood are consistent.

The per-event Fisher errors are used ONLY to fit (a, b); they are never attached
to the mock events, because errors computed from true values carry information
about the true distance that the likelihood does not model.
"""

import pickle
import numpy as np
import bilby
from numpy.linalg import inv
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline, PchipInterpolator

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
RATE_FILE   = "./gwmodel_fiducial_rate.npz"   # model rate at the truth: z, R, DL
T_OBS       = 4.0                                # years
SEED_CAL    = 1                                  # calibration draws
SEED_MOCK   = 65                                 # mock draws
N_CAL       = 2000                               # events used to fit sigma(D)
ERROR_SCALE = 1.0                                # e.g. 0.25 for errors four times smaller
USE_POISSON = True
N_FINE      = 200001
OUT         = "./CE_mock_sigmafit_Tobs%.1f_x%.2f_seed%d.pkl" % (T_OBS, ERROR_SCALE, SEED_MOCK)


# ---------------------------------------------------------------------------
# Rate, distances, and a proper sampler
# ---------------------------------------------------------------------------
d       = np.load(RATE_FILE)
z_grid  = np.asarray(d["z"])
R_grid  = np.asarray(d["R"])
DL_grid = np.asarray(d["DL"])

pos    = R_grid > 0
lnR    = PchipInterpolator(z_grid[pos], np.log(R_grid[pos]), extrapolate=False)
z_fine = np.linspace(z_grid[pos][0], z_grid[pos][-1], N_FINE)
R_fine = np.nan_to_num(np.exp(lnR(z_fine)), nan=0.0)
cdf    = cumulative_trapezoid(R_fine, z_fine, initial=0.0)
DL_of_z = CubicSpline(z_grid, DL_grid)


def draw_redshifts(rng, n):
    return np.interp(rng.random(n), cdf / cdf[-1], z_fine)


# ---------------------------------------------------------------------------
# bilby setup (as in the original generator)
# ---------------------------------------------------------------------------
minimum_frequency  = 5.
sampling_frequency = 4096
duration           = 16.

waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    waveform_arguments=dict(waveform_approximant="IMRPhenomD",
                            reference_frequency=20.,
                            minimum_frequency=minimum_frequency),
)
ifos = bilby.gw.detector.InterferometerList(["CE"])
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency, duration=duration, start_time=0)


def sigma_ln_DL_fisher(p):
    """Fisher-matrix error on ln D_L for parameters (ln Mc, ln D_L, iota)."""
    ifo  = ifos[0]
    df   = waveform_generator.frequency_array[1] - waveform_generator.frequency_array[0]
    mask = ifo.strain_data.frequency_mask
    psd  = ifo.power_spectral_density_array[mask]
    fin  = np.isfinite(psd) & (psd > 0)
    psd_f = psd[fin]

    fp = ifo.antenna_response(p["ra"], p["dec"], p["geocent_time"], p["psi"], "plus")
    fc = ifo.antenna_response(p["ra"], p["dec"], p["geocent_time"], p["psi"], "cross")

    def strain(q):
        h = waveform_generator.frequency_domain_strain(q)
        return (fp * h["plus"] + fc * h["cross"])[mask][fin]

    h0  = strain(p)
    eps = 1e-2

    dmc = eps * p["chirp_mass"]
    pp, pm = dict(p), dict(p)
    pp["chirp_mass"] += dmc
    pm["chirp_mass"] -= dmc
    dh_dlnmc = (strain(pp) - strain(pm)) / (2 * dmc) * p["chirp_mass"]

    dh_dlnDL = -h0

    dio = eps * max(p["theta_jn"], 0.01)
    pp, pm = dict(p), dict(p)
    pp["theta_jn"] += dio
    pm["theta_jn"] -= dio
    dh_dio = (strain(pp) - strain(pm)) / (2 * dio)

    def inner(a, b):
        return 4 * np.real(np.sum(np.conjugate(a) * b / psd_f) * df)

    derivs = [dh_dlnmc, dh_dlnDL, dh_dio]
    F = np.array([[inner(a, b) for b in derivs] for a in derivs])
    s = np.sin(p["theta_jn"])
    if abs(s) > 1e-6:
        F[2, 2] += (np.cos(p["theta_jn"]) / s)**2 + 1.0     # prior term on inclination
    return np.sqrt(inv(F)[1, 1])


def random_source(rng, z):
    """Masses and angles from the same priors as the original generator."""
    return dict(
        chirp_mass=rng.uniform(20, 40) * (1 + z),             # detector frame
        mass_ratio=rng.uniform(0.5, 1.0),
        luminosity_distance=float(DL_of_z(z)),
        theta_jn=np.arccos(rng.uniform(-1, 1)),
        psi=rng.uniform(0, np.pi),
        phase=rng.uniform(0, 2 * np.pi),
        geocent_time=0.0,
        ra=rng.uniform(0, 2 * np.pi),
        dec=np.arcsin(rng.uniform(-1, 1)),
        a_1=0.0, a_2=0.0, tilt_1=0.0, tilt_2=0.0, phi_12=0.0, phi_jl=0.0,
    )


# ---------------------------------------------------------------------------
# Stage 1: calibrate sigma(D) = a * D^b
# ---------------------------------------------------------------------------
rng_cal = np.random.default_rng(SEED_CAL)
z_cal   = draw_redshifts(rng_cal, N_CAL)
DL_cal  = DL_of_z(z_cal)
sig_cal = np.empty(N_CAL)
for i, zc in enumerate(z_cal):
    sig_cal[i] = DL_cal[i] * sigma_ln_DL_fisher(random_source(rng_cal, zc))
    if (i + 1) % 200 == 0:
        print(f"calibration: {i + 1}/{N_CAL}")

ok      = np.isfinite(sig_cal) & (sig_cal > 0)
b, ln_a = np.polyfit(np.log(DL_cal[ok]), np.log(sig_cal[ok]), 1)
a       = ERROR_SCALE * np.exp(ln_a)
scatter = np.std(np.log(sig_cal[ok]) - (ln_a + b * np.log(DL_cal[ok])))
print(f"fitted sigma(D) = {np.exp(ln_a):.4g} * D^{b:.3f};  used a = {a:.4g} (scale {ERROR_SCALE})")
print(f"scatter of per-event errors about the fit: {scatter:.3f} in ln(sigma)")

sigma_of_D = lambda D: a * np.abs(D)**b

# --- plot of the fit, to judge its quality -----------------------------------
import matplotlib
matplotlib.use("Agg")                          # works without a display
import matplotlib.pyplot as plt

D_plot  = np.geomspace(DL_cal[ok].min(), DL_cal[ok].max(), 200)
fit_raw = np.exp(ln_a) * D_plot**b             # fit before ERROR_SCALE
resid   = np.log(sig_cal[ok]) - (ln_a + b * np.log(DL_cal[ok]))

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 6.5), sharex=True,
                               gridspec_kw={"height_ratios": [3, 1]})
ax1.loglog(DL_cal[ok] / 1e3, sig_cal[ok] / 1e3, ".", ms=3, alpha=0.4, label="Fisher errors (calibration)")
ax1.loglog(D_plot / 1e3, fit_raw / 1e3, "k-", lw=2,
           label=rf"fit: $\sigma = {np.exp(ln_a):.3g}\,D_L^{{{b:.3f}}}$")
ax1.fill_between(D_plot / 1e3, fit_raw * np.exp(-scatter) / 1e3, fit_raw * np.exp(scatter) / 1e3,
                 color="k", alpha=0.12, label=rf"$\pm 1$ scatter ({scatter:.2f} in $\ln\sigma$)")
if ERROR_SCALE != 1.0:
    ax1.loglog(D_plot / 1e3, ERROR_SCALE * fit_raw / 1e3, "r--", lw=1.5,
               label=f"used in mock (x{ERROR_SCALE})")
ax1.set_ylabel(r"$\sigma_{D_L}$ [Gpc]")
ax1.legend(fontsize=8)

ax2.semilogx(DL_cal[ok] / 1e3, resid, ".", ms=3, alpha=0.4)
ax2.axhline(0, color="k", lw=1)
edges = np.geomspace(DL_cal[ok].min(), DL_cal[ok].max(), 11)
cent  = np.sqrt(edges[1:] * edges[:-1])
med   = [np.median(resid[(DL_cal[ok] >= lo) & (DL_cal[ok] < hi)])
         if np.any((DL_cal[ok] >= lo) & (DL_cal[ok] < hi)) else np.nan
         for lo, hi in zip(edges[:-1], edges[1:])]
ax2.plot(cent / 1e3, med, "ro-", ms=4, label="binned median")
ax2.set_xlabel(r"$D_L$ [Gpc]")
ax2.set_ylabel(r"$\ln(\sigma / \sigma_{\rm fit})$")
ax2.legend(fontsize=8)

fig.tight_layout()
plot_file = OUT.replace(".pkl", "_sigma_fit.png")
fig.savefig(plot_file, dpi=150)
print("fit plot saved:", plot_file)


# ---------------------------------------------------------------------------
# Stage 2: the mock catalogue
# ---------------------------------------------------------------------------
rng  = np.random.default_rng(SEED_MOCK)
zeta = T_OBS * cdf[-1]
N    = rng.poisson(zeta) if USE_POISSON else int(round(zeta))
z    = draw_redshifts(rng, N)

DL_true = DL_of_z(z)
DL_obs  = DL_true + rng.normal(0.0, sigma_of_D(DL_true))

cat = {
    "z_true":         z,
    "R_true":         np.nan_to_num(np.exp(lnR(z)), nan=0.0),
    "DL_true":        DL_true,
    "DL_obs":         DL_obs,
    "T_obs":          T_OBS,
    "expected_event": zeta,
    "observed_event": N,
    "sigma_model":    {"a": a, "b": b},
    "error_scale":    ERROR_SCALE,
    "noise_seed":     SEED_MOCK,
    "calibration":    {"z": z_cal, "DL": DL_cal, "sigma_DL": sig_cal,
                       "a_unscaled": np.exp(ln_a), "b": b, "scatter": scatter},
    "rate_grid":      {"z": z_grid, "R": R_grid, "DL": DL_grid},
}
with open(OUT, "wb") as f:
    pickle.dump(cat, f)

r = (DL_obs - DL_true) / sigma_of_D(DL_true)
print(f"zeta = {zeta:.1f}   N = {N}   (N - zeta)/sqrt(zeta) = {(N - zeta) / np.sqrt(zeta):+.2f}")
print(f"noise residuals: mean {r.mean():+.4f}, std {r.std():.4f}")
print("saved:", OUT)
