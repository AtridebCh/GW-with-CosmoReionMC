import numpy as np
import bilby
from scipy.integrate import simpson
from numpy.linalg import inv
from scipy.interpolate import interp1d
import pickle
import sys
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Reionization model — Fortran f2py compiled module
# change the path
# ---------------------------------------------------------------------------
sys.path.append('../ReionizationWithF2Py/python_dir')
import reion_f as f
from reion_f import run_model


np.random.seed(42)

# ---------------------------------------------------------------------------
# Redshift grid — must match Fortran exactly
# ---------------------------------------------------------------------------
zstart, zend, dz = 30.0, 0.0, 0.2
Z_arraySize = int(round(abs((zend - zstart) / dz)))
Z_reion = np.array([ik * (-abs(dz)) + zstart for ik in range(Z_arraySize + 1)])


# ---------------------------------------------------------------------------
# Run reionization model once at fiducial parameters
# ---------------------------------------------------------------------------
e_SF_II =0.01
(Z_reion, QH_Q, dNLLdz, gamma_PI, 
         sfr_popII, sfr_popIII,
         dvc_dz, D_L, age_Gyr,
         tau_factor, ierr
            ) = run_model(
                h0=6.711000e+01,
                ombh2=2.260000e-02,
                omch2=1.179535e-01,
                ns=9.600000e-01,
                sigma_8=0.8022147706311814,
                esc_popii=0.0036,
                esc_popiii= 0.0,
                lambda0= 2.44,
                zstart_in=zstart,
                zend_in=zend,
                dz_in=dz,
                z_arraysize = Z_arraySize
            )

print(f"sfr_popII  range: {sfr_popII.min():.3e} – {sfr_popII.max():.3e}")
print(f"lumdist    range: {D_L.min():.3e} – {D_L.max():.3e}")
print(f"dvc_dz     range: {dvc_dz.min():.3e} – {dvc_dz.max():.3e}")


t_to_z_interp = interp1d(age_Gyr, Z_reion, fill_value   = 'extrapolate')
z_to_t_interp = interp1d(Z_reion, age_Gyr, fill_value   = 'extrapolate')
DL_interp   = interp1d(Z_reion, D_L, fill_value   = 'extrapolate')
_dVc_interp = interp1d(Z_reion, dvc_dz, bounds_error=False, fill_value='extrapolate')
# ---------------------------------------------------------------------------
# Log-space SFRD interpolator
# ---------------------------------------------------------------------------


_sfr_log     = np.log10(np.maximum(sfr_popII, 1e-30))
_sfrd_interp = interp1d(
    Z_reion, _sfr_log,
    kind         = 'linear',
    bounds_error = False,
    fill_value   = 'extrapolate'   # exp(-inf) = 0 outside range
)


# ---------------------------------------------------------------------------
# Cosmological functions
# ---------------------------------------------------------------------------
def hubble_gyr(z, H0, ombh2, omch2):
    """H(z) in Gyr^-1 for flat LCDM."""
    h   = H0 / 100.0
    Om0 = (ombh2 + omch2) / h**2
    OmL = 1.0 - Om0
    E_z = np.sqrt(Om0 * (1 + z)**3 + OmL)
    return H0 * E_z * 1.022e-3   # Gyr^-1


# ---------------------------------------------------------------------------
# Volumetric merger rate script_R(z_m) — Eq 11
# With delta-function delay time the integral collapses to:
#   script_R(z_m) = R_f(z_f*) * |dt_f/dz_f|^{-1}
#                 = sfrd(z_f*) * H(z_f*) / (1 + z_f*)
# where dt/dz = (1+z) / H(z)  =>  dz/dt = H(z) / (1+z)
# ---------------------------------------------------------------------------
N_events      = 100
z_min         = 6.0
z_max         = 25.0
snr_threshold = 8.0

t_delay = 0.1   # Gyr


z_grid = np.linspace(z_min, z_max, 500)
t_merger_grid = z_to_t_interp(z_grid) 
t_form_grid   = t_merger_grid - t_delay
valid_grid    = t_form_grid > 0

zf_grid = np.where(
    valid_grid,
    t_to_z_interp(np.maximum(t_form_grid, 1e-3)),
    0.0
)

dzf_dtf = hubble_gyr(zf_grid, 6.711000e+01, 2.260000e-02, 1.179535e-01) / (1.0 + zf_grid) / 1e9

script_R_unnorm = np.where(
    valid_grid,
    10**_sfrd_interp(zf_grid) * dzf_dtf,
    0.0
)

# --- Normalise: script_R(z_m=4.822) = 20 Gpc^-3 yr^-1 ---
script_R_at_4822 = interp1d(z_grid, script_R_unnorm,
                                        bounds_error=False,
                                        fill_value='extrapolate')(4.822)
norm             = 20.0 / script_R_at_4822

# R(z_m) = dVc/dz(z_m) / (1+z_m) * script_R(z_m) — fully vectorised
R_z = (_dVc_interp(z_grid) / 1e9) / (1.0 + z_grid) * norm * script_R_unnorm


p_z = np.maximum(R_z, 0.0)
p_z /= simpson(p_z, z_grid)

cdf = np.cumsum(p_z) * (z_grid[1] - z_grid[0])
cdf /= cdf[-1]


def sample_redshift():
    u = np.random.rand()
    return np.interp(u, cdf, z_grid)

# ---------------------------------------------------------------------------
# Redshift sampling from physical merger rate R(z_m)
# ---------------------------------------------------------------------------

minimum_frequency  = 5.
sampling_frequency = 4096
duration           = 16.

waveform_arguments = dict(
    waveform_approximant="IMRPhenomD",
    reference_frequency=20.,
    minimum_frequency=minimum_frequency,
)

waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    waveform_arguments=waveform_arguments,
)

ifos = bilby.gw.detector.InterferometerList(["CE"])


def compute_snr(params):
    ifos.set_strain_data_from_power_spectral_densities(
        sampling_frequency=sampling_frequency,
        duration=duration,
        start_time=0,
    )
    h_f_dict    = waveform_generator.frequency_domain_strain(params)
    h_projected = ifos[0].get_detector_response(h_f_dict, params)
    snr         = np.sqrt(np.real(ifos[0].optimal_snr_squared(signal=h_projected)))
    return snr


def fisher_gaussian_posterior(injection_parameters, nsamples=3000):
    ifo = ifos[0]
    df  = waveform_generator.frequency_array[1] - waveform_generator.frequency_array[0]
    mask        = ifo.strain_data.frequency_mask
    psd         = ifo.power_spectral_density_array[mask]
    finite_mask = np.isfinite(psd) & (psd > 0)
    psd_f       = psd[finite_mask]

    fp = ifo.antenna_response(
        injection_parameters["ra"], injection_parameters["dec"],
        injection_parameters["geocent_time"], injection_parameters["psi"], "plus"
    )
    fc = ifo.antenna_response(
        injection_parameters["ra"], injection_parameters["dec"],
        injection_parameters["geocent_time"], injection_parameters["psi"], "cross"
    )

    def projected_strain(params):
        h = waveform_generator.frequency_domain_strain(params)
        return (fp * h["plus"] + fc * h["cross"])[mask][finite_mask]

    h0 = projected_strain(injection_parameters)

    eps      = 1e-2
    delta_mc = eps * injection_parameters["chirp_mass"]
    p_plus   = injection_parameters.copy(); p_plus["chirp_mass"]  += delta_mc
    p_minus  = injection_parameters.copy(); p_minus["chirp_mass"] -= delta_mc
    dh_dln_mc = (projected_strain(p_plus) - projected_strain(p_minus)) / (2 * delta_mc)

    dh_dln_DL = -h0

    delta_iota = eps * max(injection_parameters["theta_jn"], 0.01)
    p_plus  = injection_parameters.copy(); p_plus["theta_jn"]  += delta_iota
    p_minus = injection_parameters.copy(); p_minus["theta_jn"] -= delta_iota
    dh_diota = (projected_strain(p_plus) - projected_strain(p_minus)) / (2 * delta_iota)

    def inner(a, b):
        return 4 * np.real(np.sum(np.conjugate(a) * b / psd_f) * df)

    derivs = [dh_dln_mc, dh_dln_DL, dh_diota]
    n      = len(derivs)
    F      = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            F[i, j] = inner(derivs[i], derivs[j])

    iota     = injection_parameters["theta_jn"]
    sin_iota = np.sin(iota)
    if abs(sin_iota) > 1e-6:
        F[2, 2] += (np.cos(iota) / sin_iota)**2 + 1.0

    cov         = inv(F)
    
    '''
    rho_DL_iota = cov[1, 2] / np.sqrt(cov[1, 1] * cov[2, 2])

    print(f"  Fisher cond: {np.linalg.cond(F):.2e}  "
          f"sigma_ln_DL: {np.sqrt(cov[1,1]):.4f}  "
          f"1/SNR: {1/inner(h0,h0)**0.5:.4f}  "
          f"rho(DL,iota): {rho_DL_iota:.3f}")
    '''

    return cov


catalogue = []

while len(catalogue) < N_events:

    z  = sample_redshift()
    DL = DL_interp(z)

    Mc_source = np.random.uniform(20, 40)
    q         = np.random.uniform(0.5, 1.0)
    Mc_det    = Mc_source * (1 + z)

    theta_jn = np.arccos(np.random.uniform(-1, 1))
    psi      = np.random.uniform(0, np.pi)
    ra       = np.random.uniform(0, 2 * np.pi)
    dec      = np.arcsin(np.random.uniform(-1, 1))
    phase    = np.random.uniform(0, 2 * np.pi)

    injection_parameters = dict(
        chirp_mass=Mc_det,
        mass_ratio=q,
        luminosity_distance=DL,
        theta_jn=theta_jn,
        psi=psi,
        phase=phase,
        geocent_time=0.0,
        ra=ra,
        dec=dec,
        a_1=0.0,
        a_2=0.0,
        tilt_1=0.0,
        tilt_2=0.0,
        phi_12=0.0,
        phi_jl=0.0
    )

    snr = compute_snr(injection_parameters)

    if snr < snr_threshold:
        continue

    cov_samples = fisher_gaussian_posterior(injection_parameters)

    sigma_ln_DL = np.sqrt(cov_samples[1, 1])
    sigma_DL    = DL * sigma_ln_DL

    catalogue.append({
        "z_true":            z,
        "injection":         injection_parameters, 
        "snr":               snr,
        "sigma_DL":          sigma_DL,
        "sigma_ln_DL":       sigma_ln_DL,
    })

with open("CE_mock_catalog.pkl", "wb") as f:
    pickle.dump(catalogue, f)

print(f"\nCatalogue saved: {len(catalogue)} events")
