import numpy as np

import emcee

from getdist import MCSamples
from getdist import plots

backend = emcee.backends.HDFBackend("test_run_CMB_Reion_GW.h5", read_only=True)

# Remove burn-in (example: first 1000 iterations)
burnin = 1
names  = ['hubble', 'ombh2', 'omch2', 'As', 'ns', 'f_zero', 'alpha_lo', 'alpha_hi', 'alpha_z', 'log_esc_II', 'lambda0', 'f_X', 'f_alpha', 'tau']
labels = ['H_{0}','\omega_{b} h^{2}', '\omega_{c} h^{2}', 'A_s', 'n_s', 'f_{0}', '\\alpha_{lo}', '\\alpha_{hi}', '\\alpha_{z}', 'f_{esc, II}', '\lambda_0', 'f_{X}', 'f_{\\alpha}', '\\tau']

# Shape: (nsteps-burnin, nwalkers, ndim)
samples = backend.get_chain(
    discard=burnin,
    flat=True
)

loglikes = -backend.get_log_prob(
    discard=burnin,
    flat=True
)
blobs = backend.get_blobs(
        discard=burnin, 
        flat=True)

all_samples = np.column_stack([samples, blobs[:, 0]])

mc = MCSamples(
    samples=samples,
    loglikes=loglikes,
    names=names,
    labels=labels,
)

g = plots.get_subplot_plotter()

g.triangle_plot(
    [mc],
    filled=True,
)
g.export('test_run_CMB_Reion_GW_getdist_from_backend.pdf')
