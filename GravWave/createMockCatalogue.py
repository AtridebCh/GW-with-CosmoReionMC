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
sys.path.append('/Users/users/achatterjee/Data/GW_with_CosmoReionMC/ReionizationWithF2Py/python_dir')
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
fzero    = 10**(-0.43567402)
alpha_lo = 0.0 
alpha_hi = 0.0

h0 = 6.7380389E+01
ombh2 = 2.2358184E-02
omch2 = 1.1992849E-01

(Z_reion, QH_Q, dNLLdz, gamma_PI, 
         sfr_popII, sfr_popIII,
         dvc_dz, D_L, age_Gyr,
         tau_factor, ndot_H, ierr
            ) = run_model(
                h0     = h0,
                ombh2  = ombh2,
                omch2  = omch2,
                ns     = 9.6480108E-01,
                sigma_8=  0.815188053487577, # 0.814733 was waht Barun checked with best-fit from CMB+Reion
                fzero    = fzero, 
                alpha_lo = alpha_lo, 
                alpha_hi = alpha_hi,
                alpha_z  = 0.0,
                esc_popii = 10**(-1.7893972),
                lambda0= 5.5645855E+00,
                zstart_in=zstart,
                zend_in=zend,
                dz_in=dz,
                z_arraysize = Z_arraySize
            )




fzero    = 10**(-0.43567402)
alpha_lo = 0.0
alpha_hi = 0.0



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
#N_events      = [6662] #638 was for 70 days, 6662 for 2 yr
N_events_list = [13274] #[100] #[6636]  #np.logspace(3, 4, 10)

#print(N_events_list)
z_min         = 6.0
z_max         = 25.0
snr_threshold = 8.0

t_delay = 0.1   # Gyr


z_grid = np.linspace(z_min, z_max, 280)
t_merger_grid = z_to_t_interp(z_grid) 
t_form_grid   = t_merger_grid - t_delay
valid_grid    = t_form_grid > 0

zf_grid = np.where(
    valid_grid,
    t_to_z_interp(np.maximum(t_form_grid, 1e-3)),
    0.0
)

dzf_dtf = hubble_gyr(zf_grid, h0, ombh2, omch2) / (1.0 + zf_grid) / 1e9


script_R_unnorm = np.where(
    valid_grid,
    10**_sfrd_interp(zf_grid) * dzf_dtf,
    0.0
)

# --- Normalise: script_R(z_m=4.822) = 20 Gpc^-3 yr^-1 ---
script_R_at_4822 = interp1d(z_grid, script_R_unnorm,
                                        bounds_error=False,
                                        fill_value='extrapolate')(4.822)
norm             = 20.0 / script_R_at_4822 #20.0 was from Tiwari 26, with 3.5 I have checked that 100 events occur within T_obs (as shown in the output if you run this code)
#print(norm)
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

    return cov



def run_injection_campaign(n_inj=200_000, z_min=z_min, z_max=z_max, seed=1234):
    """
    Broad injection campaign to estimate p_det(z), independent of Lambda.
    IMPORTANT: z is drawn UNIFORMLY here (the injection/proposal distribution)
    -- deliberately NOT via sample_redshift()/cdf (which encodes the
    FIDUCIAL R(z|Lambda_fid)). Using the fiducial-weighted draw here would
    entangle the selection function with one particular Lambda and give poor
    coverage anywhere your fit explores away from it.
    Mc_source, q, and angles ARE drawn from the same fixed priors as the
    real catalogue loop below, since those don't depend on Lambda.
    """
    rng = np.random.default_rng(seed)
 
    z         = rng.uniform(z_min, z_max, size=n_inj)
    Mc_source = rng.uniform(20, 40, size=n_inj)
    q         = rng.uniform(0.5, 1.0, size=n_inj)
    Mc_det    = Mc_source * (1 + z)
 
    theta_jn = np.arccos(rng.uniform(-1, 1, size=n_inj))
    psi      = rng.uniform(0, np.pi, size=n_inj)
    ra       = rng.uniform(0, 2 * np.pi, size=n_inj)
    dec      = np.arcsin(rng.uniform(-1, 1, size=n_inj))
    phase    = rng.uniform(0, 2 * np.pi, size=n_inj)
 
    DL = DL_interp(z)
 
    snr = np.empty(n_inj)
    for i in range(n_inj):
        params = dict(
            chirp_mass=Mc_det[i], mass_ratio=q[i], luminosity_distance=DL[i],
            theta_jn=theta_jn[i], psi=psi[i], phase=phase[i],
            geocent_time=0.0, ra=ra[i], dec=dec[i],
            a_1=0.0, a_2=0.0, tilt_1=0.0, tilt_2=0.0, phi_12=0.0, phi_jl=0.0,
        )
        snr[i] = compute_snr(params)
 
    found = snr >= snr_threshold
    return dict(z=z, Mc_source=Mc_source, q=q, snr=snr, found=found)

def estimate_p_det(injections, n_bins=60, z_min=z_min, z_max=z_max):
    """
    Binned found/injected ratio + monotone smoothing (isotonic regression;
    falls back to a moving average if sklearn isn't installed). Detection
    probability is expected to be non-increasing in z.
    """
    z           = injections["z"]
    bin_edges   = np.linspace(z_min, z_max, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
 
    counts_tot,   _ = np.histogram(z, bins=bin_edges)
    counts_found, _ = np.histogram(z[injections["found"]], bins=bin_edges)
 
    with np.errstate(invalid="ignore", divide="ignore"):
        p_hat = np.where(counts_tot > 0, counts_found / counts_tot, 0.0)
        se    = np.where(counts_tot > 0,
                          np.sqrt(p_hat * (1 - p_hat) / np.maximum(counts_tot, 1)),
                          np.nan)
 
    try:
        from sklearn.isotonic import IsotonicRegression
        valid = counts_tot > 0
        iso = IsotonicRegression(increasing=False, out_of_bounds="clip")
        iso.fit(bin_centers[valid], p_hat[valid])
        p_smooth = iso.predict(bin_centers)
    except ImportError:
        kernel = np.ones(5) / 5
        p_smooth = np.clip(np.convolve(p_hat, kernel, mode="same"), 0, 1)
 
    p_det_interp = interp1d(bin_centers, p_smooth, bounds_error=False,
                             fill_value=(p_smooth[0], p_smooth[-1]))
    diagnostics = dict(bin_centers=bin_centers, p_hat=p_hat, se=se,
                        counts_tot=counts_tot, counts_found=counts_found,
                        p_smooth=p_smooth)
    return p_det_interp, diagnostics
 
 
def mu_of_lambda(R_zgrid, p_det_interp, T_obs_yr, z_g):
    """
    mu(Lambda) as a cheap 1D integral, reusing the SAME p_det(z) for every
    trial Lambda (n_s, A_s) -- no new injections needed per Lambda.
    R_of_z_func(z) should be built the SAME way as R_z above (script_R,
    dVc/dz, norm, etc.) but recomputed at the trial n_s / sigma_8 via
    run_model(...), then wrapped in interp1d, e.g.:
        R_trial_interp = interp1d(z_grid, R_z_trial, bounds_error=False, fill_value=0.0)
    """
    return T_obs_yr * simpson(R_zgrid, z_g)
 
'''
 
# --- Run it once, here, before generating the catalogue -------------------
injections = run_injection_campaign(n_inj=2000)
print(f"Found fraction overall: {injections['found'].mean():.4f}")

p_det_interp, p_det_diag = estimate_p_det(injections, n_bins=60)
'''

# sanity check: mu(fiducial Lambda) should come out ~ N_events if T_obs is
# set/calibrated correctly -- if this disagrees a lot with your catalogue
# acceptance rate, something (units/normalization/z-range) is inconsistent
# between the generation loop and this selection-function estimate.
#R_fid_interp = interp1d(z_grid, R_z, bounds_error=False, fill_value=0.0)
T_obs_yr = 4.0 #0.03 #2.0 #70.0/365.25  # <-- set this to your actual CE observing time
mu_fid = mu_of_lambda(R_z, 1.0, T_obs_yr, z_grid)
print(f"mu(fiducial Lambda) = {mu_fid:.3f}  (compare to N_events = {N_events_list[0]})")

'''
with open("../GravWave/selection_function_p_det.pkl", "wb") as fpick:
    pickle.dump({"p_det_interp": p_det_interp, "diagnostics": p_det_diag}, fpick)
'''

def make_catalogue(event_No):

    catalogue = {
    "z_true":      [],
    "R_true":      [],
    "injection":   [],
    "snr":         [],
    "sigma_DL":    [],
    "sigma_ln_DL": [],
}

    #while len(catalogue) < event_No:
    for i in range(event_No):

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

        if snr < 0.0: #:snr_threshold:
            continue

        cov_samples = fisher_gaussian_posterior(injection_parameters)

        sigma_ln_DL = np.sqrt(cov_samples[1, 1])
        sigma_DL    = DL * sigma_ln_DL
        '''
        catalogue.append({
            "z_true":            z,
            "R_true":      float(np.interp(z, z_grid, R_z)),   # fiducial R(z_true)
            "injection":         injection_parameters,
            "snr":               snr,
            "sigma_DL":          sigma_DL,
            "sigma_ln_DL":       sigma_ln_DL,
            })
        '''
        catalogue["z_true"].append(z)
        catalogue["R_true"].append(float(np.interp(z, z_grid, R_z)))
        catalogue["injection"].append(injection_parameters)
        catalogue["snr"].append(snr)
        catalogue["sigma_DL"].append(sigma_DL)
        catalogue["sigma_ln_DL"].append(sigma_ln_DL)
    catalogue["expected_event"] = event_No
    catalogue["observed_event"] = len(catalogue["z_true"])
    catalogue["T_obs"] = T_obs_yr
    
    with open("../GravWave/CE_mock_catalog_%d_events_without_snr.pkl"%event_No, "wb") as f:
        pickle.dump(catalogue, f)
    print(f"\nCatalogue saved: {catalogue['expected_event']} events")
    return 1.0

for N_event in N_events_list:
    make_catalogue(N_event)

