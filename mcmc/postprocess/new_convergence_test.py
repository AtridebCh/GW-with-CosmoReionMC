import os
import numpy as np
import matplotlib.pyplot as plt
from emcee.autocorr import function_1d
# ============================================================
# USER INPUTS
# ============================================================


ndim = 14


file_root = 'test_MCMC'

data_dir = '/Users/users/achatterjee/Data/GW_with_CosmoReionMC/mcmc/MCMC_sampler/chain_storage/run_with_merger_rate/'

nwalkers = 52



param_name_arr = [
    "hubble",
    "ombh2",
    "omch2",
    "As",
    "ns",
    "f_zero",
    "alpha_lo",
    "alpha_hi",
    "alpha_z",
    "esc_II",
    "lambda0",
    "f_X",
    "f_alpha",
    "tau"
]

burnin_fraction = 0.10

# ============================================================
# READ CHAINS
# ============================================================

for k in range(nwalkers):

    filename = os.path.join(
        data_dir,
        f"{file_root}_{k+1}.txt"
    )

    all_data = np.loadtxt(filename, unpack=True)

    if k == 0:

        mcmc_steps = all_data.shape[1]

        lnprob = np.zeros(
            (nwalkers, mcmc_steps)
        )

        chains = np.zeros(
            (nwalkers, mcmc_steps, ndim)
        )

    lnprob[k, :] = -all_data.T[:mcmc_steps, 1]

    chains[k, :, :] = all_data.T[
        :mcmc_steps,
        2:ndim+2
    ]

# ============================================================
# BURN-IN
# ============================================================

burnin = int(burnin_fraction * mcmc_steps)

post_chain = chains[:, burnin:, :]

nsteps_post = post_chain.shape[1]

print()
print("====================================================")
print(f"Number of walkers      = {nwalkers}")
print(f"Number of dimensions   = {ndim}")
print(f"Total MCMC steps       = {mcmc_steps}")
print(f"Burn-in steps          = {burnin}")
print(f"Post burn-in steps     = {nsteps_post}")
print("====================================================")

# ============================================================
# AUTOCORRELATION FUNCTIONS
# (from emcee tutorial)
# ============================================================

def next_pow_two(n):

    i = 1

    while i < n:
        i <<= 1

    return i



def autocorr_func_1d(x, norm=True):

    x = np.atleast_1d(x)

    if len(x.shape) != 1:
        raise ValueError(
            "invalid dimensions for 1D autocorrelation"
        )

    n = next_pow_two(len(x))

    f = np.fft.fft(
        x - np.mean(x),
        n=2*n
    )

    acf = np.fft.ifft(
        f * np.conjugate(f)
    )[:len(x)].real

    acf /= 4*n

    if norm:
        acf /= acf[0]

    return acf


def gelman_rubin_param(chains):
    """
    chains shape = (nchains, nsteps)
    """

    m, n = chains.shape

    chain_means = np.mean(chains, axis=1)
    chain_vars  = np.var(chains, axis=1, ddof=1)

    W = np.mean(chain_vars)

    B = n * np.var(chain_means, ddof=1)

    var_hat = ((n - 1)/n) * W + B/n

    Rhat = np.sqrt(var_hat / W)

    return Rhat - 1.0

def auto_window(taus, c):

    m = np.arange(len(taus)) < c * taus

    if np.any(m):
        return np.argmin(m)

    return len(taus) - 1


def autocorr_new(y, c=5.0):

    # y shape = (nwalkers, nsteps)

    f = np.zeros(y.shape[1])

    for yy in y:

        f += autocorr_func_1d(yy)

    f /= len(y)

    taus = 2.0 * np.cumsum(f) - 1.0

    window = auto_window(taus, c)

    return taus[window]

# ============================================================
# ESTIMATE TAU FOR EACH PARAMETER
# ============================================================

taus = []

print()
print("====================================================")
print("AUTOCORRELATION ANALYSIS")
print("====================================================")

for j in range(ndim):

    y = post_chain[:, :, j]

    rminus1 = gelman_rubin_param(
        post_chain[:, :, j]
    )

    tau = autocorr_new(y)

    taus.append(tau)

    print(
        f"{param_name_arr[j]:15s} "
        f"tau = {tau:.2f}  "
        f"R-1 = {rminus1:.5f}"
    )

# ============================================================
# CONVERGENCE CHECK
# ============================================================

print()
print("====================================================")
print("CONVERGENCE CHECK")
print("====================================================")

for name, tau in zip(param_name_arr, taus):

    ratio = nsteps_post / tau

    status = (
        "GOOD"
        if ratio > 50
        else "MORE STEPS NEEDED"
    )

    print(
        f"{name:15s} "
        f"N/tau = {ratio:.1f}   "
        f"--> {status}"
    )

# ============================================================
# EFFECTIVE SAMPLE SIZE
# ============================================================

print()
print("====================================================")
print("EFFECTIVE SAMPLE SIZE")
print("====================================================")

Ntot = nwalkers * nsteps_post

for name, tau in zip(param_name_arr, taus):

    neff = Ntot / tau

    print(
        f"{name:15s} "
        f"N_eff = {neff:.0f}"
    )

# ============================================================
# TUTORIAL-STYLE TAU VS CHAIN LENGTH PLOTS
# ============================================================

N = np.exp(
    np.linspace(
        np.log(100),
        np.log(nsteps_post),
        15
    )
).astype(int)

for j in range(ndim):

    y = post_chain[:, :, j]

    tau_est = np.empty(len(N))

    for i, n in enumerate(N):

        tau_est[i] = autocorr_new(
            y[:, :n]
        )

    plt.figure(figsize=(6, 4))

    plt.loglog(
        N,
        tau_est,
        "o-",
        label=param_name_arr[j]
    )

    plt.loglog(
        N,
        N / 50.0,
        "--",
        label=r"$N/50$"
    )

    plt.xlabel("Chain length")
    plt.ylabel(r"$\tau$")
    plt.title(param_name_arr[j])

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        f"tau_convergence_{param_name_arr[j]}.pdf"
    )

    plt.close()

print()
print("Saved convergence plots.")
