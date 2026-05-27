import sys
import csv
import numpy as np
import emcee
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 10


# ============================================================
# SETTINGS
# ============================================================

file_root = 'test_MCMC'

data_dir = '/home/atri/GW_with_CosmoReionMC/mcmc/MCMC_sampler/chain_storage/test_run25_05_2026/'

burnin = 0.3        # fraction or integer
nwalkers = 60


# ============================================================
# LOAD PARAMETER NAMES
# ============================================================

filename = data_dir + file_root + '.paramnames'

with open(filename) as f:

    reader = csv.reader(f, delimiter=" ")

    param_name_arr = np.asarray(list(zip(*reader))[0])


ndim = len(param_name_arr)

print("Parameters:")
print(param_name_arr)
print()


# ============================================================
# LOAD CHAINS
# ============================================================

for k in range(nwalkers):

    filename = data_dir + '/' + file_root + "_" + str(k+1) + ".txt"

    all_data = np.loadtxt(filename, unpack=True)

    if k == 0:

        mcmc_steps = np.shape(all_data)[1]

        lnprob = np.zeros([nwalkers, mcmc_steps])

        chains = np.zeros([nwalkers, mcmc_steps, ndim])

    lnprob[k, :] = -all_data.T[:mcmc_steps, 1]

    chains[k, :, :] = all_data.T[:mcmc_steps, 2:ndim+2]


print("Number of walkers =", nwalkers)
print("Number of dimensions =", ndim)
print("Total MCMC steps =", mcmc_steps)


# ============================================================
# BURN-IN
# ============================================================

if burnin > 1:

    burnin = int(burnin)

    if burnin >= mcmc_steps:
        burnin = 0

else:

    burnin = int(mcmc_steps * burnin)

print("Burn-in steps =", burnin)

effective_steps = mcmc_steps - burnin

print("Post burn-in steps =", effective_steps)
print()


# ============================================================
# AUTOCORRELATION ANALYSIS
# ============================================================

print("====================================================")
print("AUTOCORRELATION ANALYSIS")
print("====================================================")

taus = []

for j in range(ndim):

    try:

        tau = emcee.autocorr.integrated_time(
            chains[:, burnin:, j],
            quiet=True
        )

        tau_mean = np.mean(tau)

        taus.append(tau_mean)

        print(
            f"{param_name_arr[j]:15s} "
            f"tau = {tau_mean:.2f}"
        )

    except Exception as e:

        taus.append(np.nan)

        print(
            f"{param_name_arr[j]:15s} "
            f"FAILED ({e})"
        )

print()


# ============================================================
# CONVERGENCE CHECK
# ============================================================

print("====================================================")
print("CONVERGENCE CHECK")
print("====================================================")

for j in range(ndim):

    tau = taus[j]

    if np.isnan(tau):
        continue

    ratio = effective_steps / tau

    if ratio < 10:
        status = "NOT converged"

    elif ratio < 30:
        status = "Poor"

    elif ratio < 50:
        status = "Probably OK"

    else:
        status = "GOOD"

    print(
        f"{param_name_arr[j]:15s} "
        f"N/tau = {ratio:.1f}   --> {status}"
    )

print()


# ============================================================
# EFFECTIVE SAMPLE SIZE
# ============================================================

print("====================================================")
print("EFFECTIVE SAMPLE SIZE")
print("====================================================")

for j in range(ndim):

    tau = taus[j]

    if np.isnan(tau):
        continue

    n_eff = nwalkers * effective_steps / tau

    print(
        f"{param_name_arr[j]:15s} "
        f"N_eff = {n_eff:.0f}"
    )

print()


# ============================================================
# FLATTEN SAMPLES
# ============================================================

samples = chains[:, burnin:, :]

s = samples.shape

samples = samples.reshape(s[0] * s[1], s[2])

print("Flattened samples shape =", samples.shape)
print()


# ============================================================
# TRACE PLOTS
# ============================================================

print("Making trace plots...")

fig, axes = plt.subplots(
    ndim,
    figsize=(10, 2*ndim),
    sharex=True
)

if ndim == 1:
    axes = [axes]

for j in range(ndim):

    ax = axes[j]

    for i in range(nwalkers):

        ax.plot(
            chains[i, :, j],
            alpha=0.3,
            lw=0.5
        )

    ax.axvline(
        burnin,
        color='red',
        ls='--',
        lw=1
    )

    ax.set_ylabel(param_name_arr[j])

axes[-1].set_xlabel("Step")

plt.tight_layout()

plt.savefig(
    data_dir + 'trace_plots.png',
    dpi=300
)

plt.show()


# ============================================================
# LOG PROBABILITY PLOT
# ============================================================

print("Making log-probability plot...")

plt.figure(figsize=(8, 4))

for i in range(nwalkers):

    plt.plot(
        lnprob[i],
        alpha=0.3,
        lw=0.5
    )

plt.axvline(
    burnin,
    color='red',
    ls='--'
)

plt.xlabel("Step")
plt.ylabel("-log posterior")

plt.tight_layout()

plt.savefig(
    data_dir + 'logprob_plot.png',
    dpi=300
)

plt.show()


# ============================================================
# FINAL SUMMARY
# ============================================================

print()
print("====================================================")
print("FINAL SUMMARY")
print("====================================================")

max_tau = np.nanmax(taus)

print(f"Maximum autocorrelation time = {max_tau:.2f}")

criterion = effective_steps / max_tau

print(f"N/tau_max = {criterion:.1f}")

if criterion > 50:

    print()
    print("Chain appears WELL CONVERGED")

elif criterion > 30:

    print()
    print("Chain is PROBABLY converged")

else:

    print()
    print("Chain is likely NOT converged")

print()
print("Done.")
