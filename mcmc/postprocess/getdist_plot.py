from getdist import loadMCSamples
from getdist import plots
import matplotlib.pyplot as plt
#import pylab as plt

import numpy as np


file_root = 'Infclipped'
#file_root = 'test_MCMC'
burnin = 0.0
#burnin = 0.3
plt.rcParams['font.size']=8


#params = ['epsilon_a0','epsilon_a1', 'epsilon_a2','epsilon_a3', 'epsilon_a4','epsilon_a5', 'lambda_0',  'tau']

data_path='/Users/users/achatterjee/Data/GW_with_CosmoReionMC/mcmc/postprocess/run_with_merger_rate_Infclipped/'
#data_path='/Users/users/achatterjee/Data/GW_with_CosmoReionMC/mcmc/MCMC_sampler/chain_storage/run_with_merger_rate/'

samples = loadMCSamples(data_path + file_root, settings={'ignore_rows':burnin})

likestats = samples.getLikeStats()
print (likestats)

'''
worst_value = samples.getGelmanRubin()
print(worst_value)
'''
samples.addDerived(np.log10(samples.getParams().esc_II), 'esc_II_log', 'log(f_{\\rm esc})')
samples.addDerived(np.log10(samples.getParams().f_zero), 'f_zero_log', 'log(f_{0})')
#print (samples.getNumSampleSummaryText())


margestats = samples.getMargeStats()
print (margestats)

g = plots.get_subplot_plotter(width_inch=6)
g.axes_fontsize=25

'''
hubble H_{0}
ombh2 \Omega_{b} h^{2}
omch2 \Omega_{c} h^{2}
As A_s
ns n_s
f_zero f_{0}
alpha_lo \alpha_{lo}
alpha_hi \alpha_{hi}
alpha_z \alpha_{z}
esc_II f_{esc, II}
lambda0 \lambda_0
f_X f_{X}
f_alpha f_{\alpha}
tau \tau
["hubble", "ombh2", "omch2", "As", "ns", "f_zero", "alpha_lo", "alpha_hi", "alpha_z", "esc_II", "tau"]       
'''
g.triangle_plot(samples, ["hubble", "ombh2", "omch2", "As", "ns", "f_zero", "alpha_lo", "alpha_hi", "alpha_z", "esc_II", "lambda0", "tau"], filled=True, contour_colors=['blue'], legend_labels=[r'$\mathbf{with \, \chi^2_R+\chi^2_{D_L}}$']) ##, ['log(M1)', 'log(M2)'])

     
g.export(file_root + '_getdist_run_with_merger_rate.pdf')



