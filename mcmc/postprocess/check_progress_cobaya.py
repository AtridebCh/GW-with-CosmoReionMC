import matplotlib.pyplot as plt
from cobaya.samplers.mcmc import plot_progress

# Assuming chain saved at `chains/gaussian`

plot_progress("/Users/users/achatterjee/Data/GW_with_CosmoReionMC/mcmc/MCMC_sampler/chains/run_with_CMB_Reion_GW")
plt.tight_layout()
plt.show()
plt.savefig('check_progress.pdf')
