from getdist import loadMCSamples
from getdist import plots

import matplotlib.pyplot as plt

samples = loadMCSamples("/Users/users/achatterjee/Data/GW_with_CosmoReionMC/mcmc/MCMC_sampler/chains_with_varying_log_f_esc/run_with_CMB_Reion_GW")

#samples = loadMCSamples("/Users/users/achatterjee/Data/GW_with_CosmoReionMC/mcmc/MCMC_sampler/chains_with_varying_f_esc/run_with_CMB_Reion_GW")
print('total number of sample', samples.numrows)

likestats = samples.getLikeStats()
print (likestats)

margestats = samples.getMargeStats()
print (margestats)


g = plots.get_subplot_plotter()

g.triangle_plot(
    [samples], ['hubble','ombh2','omch2', 'As', 'ns', 'f_zero', 'log_esc_II', 'lambda0', 'tau'],
    filled=True
)
plt.savefig("cobaya_triangle_varying_log_f_esc.pdf")
