'''
import pickle
import numpy as np

src = "/Users/users/achatterjee/Data/GW_with_CosmoReionMC/GravWave/CE_mock_catalog_6636_events_without_snr.pkl"
dst = src.replace(".pkl", "_noisy_sigmaD_seed65.pkl")

with open(src, "rb") as f:
    cat = pickle.load(f)

DL_true  = np.array([inj["luminosity_distance"] for inj in cat["injection"]])   # Mpc
sigma_DL = np.array(cat["sigma_DL"])                                            # per-event errors (from true SNR)

# --- fit one error relation sigma(D) = a * D^b, shared by all events ---
ok      = (DL_true > 0) & (sigma_DL > 0)
b, ln_a = np.polyfit(np.log(DL_true[ok]), np.log(sigma_DL[ok]), 1)
a       = np.exp(ln_a)
sigma_of_D = lambda D: a * np.abs(D)**b

resid = np.log(sigma_DL[ok]) - np.log(sigma_of_D(DL_true[ok]))
print(f"sigma(D) = {a:.4g} * D^{b:.3f}   (D in Mpc)")
print(f"scatter of per-event errors about the fit: {resid.std():.3f} in ln(sigma)")

# --- noise drawn from the shared relation, evaluated at the TRUE distance ---
rng          = np.random.default_rng(65)
sigma_used   = sigma_of_D(DL_true)
cat["DL_obs"]      = DL_true + rng.normal(0.0, sigma_used)   # observed = true + Gaussian error
cat["sigma_model"] = {"a": a, "b": b}                       # the relation the likelihood must use
cat["noise_seed"]  = 65

with open(dst, "wb") as f:
    pickle.dump(cat, f)

print("saved:", dst)
'''

import pickle
import numpy as np

src         = "/Users/users/achatterjee/Data/GW_with_CosmoReionMC/GravWave/CE_mock_catalog_13274_events_without_snr.pkl"  # new catalogue
error_scale = 0.25
seed        = 65
dst         = src.replace(".pkl", f"_noisy_sigmaD_x{error_scale}_seed{seed}.pkl")

with open(src, "rb") as f:
    cat = pickle.load(f)

DL_true  = np.array([inj["luminosity_distance"] for inj in cat["injection"]])   # Mpc
sigma_DL = np.array(cat["sigma_DL"])                                            # per-event errors from the generator

# fit the shared relation sigma(D) = a * D^b, then scale it
ok      = (DL_true > 0) & (sigma_DL > 0)
b, ln_a = np.polyfit(np.log(DL_true[ok]), np.log(sigma_DL[ok]), 1)
a       = error_scale * np.exp(ln_a)
sigma_of_D = lambda D: a * np.abs(D)**b
print(f"sigma(D) = {a:.4g} * D^{b:.3f}   (scale {error_scale}, D in Mpc)")

# noise drawn with the scaled error at the true distance
rng = np.random.default_rng(seed)
cat["DL_obs"]      = DL_true + rng.normal(0.0, sigma_of_D(DL_true))
cat["sigma_model"] = {"a": a, "b": b}
cat["error_scale"] = error_scale
cat["noise_seed"]  = seed

with open(dst, "wb") as f:
    pickle.dump(cat, f)

# quick checks
r = (cat["DL_obs"] - DL_true) / sigma_of_D(DL_true)
print("events:", len(DL_true), "  T_obs:", cat["T_obs"])
print("noise residuals: mean %.3f, std %.3f" % (r.mean(), r.std()))
print("saved:", dst)
