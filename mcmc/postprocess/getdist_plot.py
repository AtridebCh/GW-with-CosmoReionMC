from getdist import loadMCSamples
from getdist import plots
import matplotlib.pyplot as plt
#import pylab as plt

import numpy as np


file_root = 'Infclipped'
burnin = 0.0
plt.rcParams['font.size']=8




params = ['epsilon_a0','epsilon_a1', 'epsilon_a2','epsilon_a3', 'lambda_0', 'tau']

#params = ['epsilon_a0','epsilon_a1', 'epsilon_a2','epsilon_a3', 'epsilon_a4','epsilon_a5', 'lambda_0',  'tau']

data_path='/home/atridebchatterjee/reion_GPR/chain_storage/4_param_run25_09_2022_Infclipped/'


samples = loadMCSamples(data_path + file_root, settings={'ignore_rows':burnin})


#print (samples.getNumSampleSummaryText())

likestats = samples.getLikeStats()
print (likestats)
margestats = samples.getMargeStats()
print (margestats)

g = plots.get_subplot_plotter(width_inch=6)
g.axes_fontsize=25


g.triangle_plot(samples, params, filled=True, contour_colors=['blue'],legend_labels=[r'$\mathbf{GPR-4params}$']) ##, ['log(M1)', 'log(M2)'])


     
g.export(file_root + '_getdist_triangle_4params.pdf')



