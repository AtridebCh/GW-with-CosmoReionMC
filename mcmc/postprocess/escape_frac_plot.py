import sys

from getdist import loadMCSamples
from fgivenx import plot_contours, samples_from_getdist_chains, plot_lines, plot_dkl
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager


from pylab import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.patches import Rectangle

import numpy as np

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF

from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning


Planck={'H0':67.70, 'ombh2': 0.0223, 'omch2': 0.12, 'ns':0.96, 'sigma8':0.81} 
Astro_dict={'esc_pop_III':0.0}

sys.path.append('../reion_21cm_code/Reionization_code')
import running_reion
from running_reion import Redshift_Evolution



plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['figure.titlesize'] = 10


plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 20
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']



#lymanRedshift, lymanLimitData, lyman_error=np.loadtxt('./Lyman_limit.dat',usecols=(0,1,2),unpack=True)
#redshiftGamma, Gamma, Gamma_max, Gamma_min=np.loadtxt('./gamma_data_all_combined.dat',skiprows=1,usecols=(0,1,2,3),unpack=True)

'''
asymmetric_error = [Gamma-Gamma_min,Gamma_max-Gamma]

f_esc_upperlim_data=np.loadtxt('f_esc_upper_limit.dat')
z, z_low, z_up, f_esc_upper = f_esc_upperlim_data[:,0], f_esc_upperlim_data[:,1], f_esc_upperlim_data[:,2], f_esc_upperlim_data[:,3]

zerr=[z-z_low, z_up-z]
'''
def reion_out(Z, theta, count=3):
	eps_param = theta[0:np.size(theta)-1]
	lambda_0 = theta[-1]
	#print(lambda_0)
	esc_PopII_redshift=np.array([3, 8, 13, 18]).reshape(-1, 1)
	#esc_PopII_redshift=np.array([3, 6, 9, 12, 15, 18]).reshape(-1, 1)
	esc_Pop_II_val=np.array(eps_param).reshape(-1, 1)
	#print(f'esc_param: {eps_param}')
	GPR_interpolator=gaussian_regress_process(esc_Pop_II_val, esc_PopII_redshift)
	mean_prediction, std_prediction = GPR_interpolator.predict(Z.reshape(-1,1), return_std=True)
	return mean_prediction, GPR_interpolator, lambda_0



@ignore_warnings(category=ConvergenceWarning)
def gaussian_regress_process(esc_Pop_II_array, esc_PopII_redshift, n_restarts_optimizer=9 ):
	kernel = 1 * RBF(8.0)
	gaussian_process = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=n_restarts_optimizer)
	gaussian_process.fit(esc_PopII_redshift, esc_Pop_II_array)
	return gaussian_process

zstart = 2.0
zend = 20.0
dz_step=0.2
n=int(abs(((zend-zstart)/dz_step)))+1

Z=np.linspace(zend,zstart,n)

Z_plot=np.linspace(20.0, 2.0, int((20.0-2.0)/0.2+1))

file_root = 'Infclipped'
burnin = 0.0



#params = ['epsilon_a0','epsilon_a1', 'epsilon_a2','epsilon_a3', 'epsilon_a4','epsilon_a5', 'lambda_0']
params = ['epsilon_a0','epsilon_a1', 'epsilon_a2','epsilon_a3', 'lambda_0']
samples_fgx, weights_fgx = samples_from_getdist_chains(params, '../chain_storage/4_param_run25_09_2022_Infclipped/' + file_root, settings={'ignore_rows':burnin})



name_string='_4_param'

def epsilon_out(Z, theta):
	epsilon_out, GPR_interpolator, lambda_0 = reion_out(Z, theta)
	return epsilon_out / 0.04 #for plotting with \eps* = 0.01

np.random.seed(7)	
random_index=np.random.randint(len(samples_fgx), size=20) # change "size" to make plots more densed

fig=plt.figure(figsize=(7,5.7))

gs = gridspec.GridSpec(nrows=1,
                       ncols=1,
                       figure=fig,
                       #width_ratios= [1, 1, 1],
                       #height_ratios=[1, 1, 1],
                       wspace=0.0,
                       hspace=0.0)


ax_escfrac = fig.add_subplot(gs[0, 0])
ax_escfrac.set_ylim([-0.01, 0.5])


esc=[]
for count, params in enumerate(samples_fgx[random_index]):
	esc.append(np.array(epsilon_out(Z_plot, params)))
esc_mean = np.mean(np.asarray(esc), axis=0)
esc_std = np.std(np.asarray(esc), axis=0)


esc_upper = esc_mean+esc_std
esc_upper[esc_upper>1.0]=1.0
esc_upper2 = esc_mean+2*esc_std
esc_upper2[esc_upper2>1.0]=1.0

plt.fill_between(Z, esc_upper, esc_mean-esc_std, color='b')
plt.fill_between(Z, esc_upper2, esc_upper, alpha=0.5, color='b')
plt.fill_between(Z, esc_mean-2*esc_std, esc_mean-esc_std, alpha=0.5, color= 'b')


### esc_frac_plot ######
#cbar = plot_contours(epsilon_out, Z_plot, samples_fgx[random_index], weights=weights_fgx[random_index], ax=ax_escfrac, ntrim=50, colors='Blues_r',lines=False,contour_line_levels=[1,2],contour_color_levels=[0,1])

z, z_low, z_up, f_esc_upper = 5.8, 5.8, 5.8, 0.08
zerr=z-z_low#[z-z_low, z_up-z]


# Cristiani et al. (2016)
# Arxiv:1603.09351
zdata = 3.8
taudata = 0.053
zerr = 0.2
errorlow = [0.012]
errorup = [0.027]
asymmetric_error = np.array(list(zip(errorlow, errorup))).T
ax_escfrac.errorbar(zdata, taudata, yerr=asymmetric_error, xerr=zerr, label='Cristiani+16', marker='^', ms=10, color='r', ls='', capsize=3, mec='w')

# Matthee et al. (2017)
#https://academic.oup.com/mnras/article/465/3/3637/2544379
zdata=[2.2, 3.8, 4.9]
taudata=[0.059, 0.027, 0.06]
errorlow=[0.042, 0.023, 0.052]
errorupp=[0.145, 0.072, 0.139]

ax_escfrac.errorbar(zdata, taudata, yerr=[errorlow, errorupp], label='Matthee+17', marker='o', ms=8, color='r', ls='', capsize=3, mec='w')

#Rutkowski+17
# Arxiv: 1705.06355
zdata = 2.64
taudata = 0.056
zerr = 0.26
tauerr = 0.02
ax_escfrac.errorbar(zdata, taudata, yerr=tauerr, xerr=zerr, uplims=True, label='Rutkowski+17', marker='d', ms=10, color='r', ls='', capsize=3, mec='w')

# Kakiichi et al (2018)
# Arxiv: 1803.02981
zdata = 5.85
taudata = 0.11
zerr = 0.55
tauerr = 0.05
ax_escfrac.errorbar(zdata, taudata, yerr=tauerr, xerr=zerr, label='Kakiichi+18', marker='h', ms=10, color='r', ls='', capsize=3, mec='w')

# Mestric (2021)
# https://doi.org/10.1093/mnras/stab2615
zdata = 3.28
taudata = 0.06
zerr = 0.28
tauerr = 0.022
ax_escfrac.errorbar(zdata, taudata, yerr=tauerr, xerr=zerr, uplims=True, label='Mestric+21', marker='X', ms=10, color='r', ls='', capsize=3, mec='w')

# Pahl et al. (2022)
zdata = 3.075
taudata = 0.06
zerr = 0.325
tauerr = 0.01
ax_escfrac.errorbar(zdata, taudata, yerr=tauerr, xerr=zerr, label='Pahl+22', marker='s', ms=8, color='r', ls='', capsize=3, mec='w')

# Saldana-Lopez et al. (2022b)
# Arxiv: 2201.11800
zdata = 5.0
taudata = 0.1
zerr = 1.0
tauerr = 0.02
ax_escfrac.errorbar(zdata, taudata, yerr=tauerr, xerr=zerr, uplims=True, label='Saldana-Lopez+22', marker='p', ms=10, color='r', ls='', capsize=3, mec='w')

# Finkelstein et al. (2019)
zdata=[4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
errorlow=[0.0025, 0.004, 0.0055, 0.007, 0.008, 0.0125, 0.016, 0.0185, 0.02]
errorupp=[0.0025, 0.004, 0.0055, 0.007, 0.008, 0.0125, 0.016, 0.0185, 0.02]
taudata=[0.0075, 0.010, 0.0155, 0.022, 0.031, 0.0445, 0.058, 0.0655, 0.07]
ax_escfrac.errorbar(zdata, taudata, yerr=[errorlow, errorupp], label='Finkelstein+19', marker='s', ms=8, color='k', ls='', capsize=3, mec='k', mfc='none')


# Faisst A. L. (2016)
# Arxiv: 1605.06507
zdata = np.linspace(2,12,100)
taudata = 0.023*((1+zdata)/3.0)**1.17
ax_escfrac.plot(zdata,taudata,'k', ls='dotted',label='Faisst+16')

# Puchwein et al. (2018) - fitting model eq. 14
# Arxiv: 1801.04931
zdata = np.linspace(2,12,100)
taudata = np.minimum(0.18, 6.9*1e-5 * (1+zdata)**3.97)
ax_escfrac.plot(zdata,taudata,'k--',label='Puchwein+19')

# Faucher-Giguere (2020) - fitting model eq. 12
# https://doi.org/10.1093/mnras/staa302

zdata = np.linspace(3,12,100)
taudata = np.minimum(0.2, 0.01 * ((1+zdata)/4.0)**2.5)
ax_escfrac.plot(zdata,taudata,'k-.',label='Faucher-Giguere+20')

# Mitra+15
#f_esc_data=np.loadtxt('fesc.dat')
#z_esc, f_esc, low_err, upper_error = f_esc_data[:,0], f_esc_data[:,3], f_esc_data[:,4], #f_esc_data[:,5]
#ax_escfrac.errorbar(z_esc, f_esc, yerr=[low_err, upper_error], label='Mitra+15', marker='o', ms=8, color='k', ls='', capsize=3, mec='k', mfc='none')
ax_escfrac.plot(Z, esc_mean, color='darkorange', ls='-', label='mean', lw=5)

'''
ax_escfrac.errorbar(z_esc, f_esc, yerr=[low_err, upper_error], fmt='ro', ms=12, mec='k', capsize=3)


upperlim_f_esc=ax_escfrac.errorbar(z, f_esc_upper,  xerr=zerr, yerr=[ 0.02, 0.01, 0.02, 0.02], uplims=np.array([True]), label=['Kakiichi+18', 'Tanvir+19', 'Pahl+21', 'Mestric+21'], barsabove=True,fmt=['ro','rs','r^','rp'], ms=12, mec='k')


#Ma 2015 FIRE simulation
simulation_f_esc=ax_escfrac.errorbar([10., 8.7, 6.9], np.array([0.10, 0.15, 0.04]), yerr=[0.02, 0.03, 0.009], uplims=np.array([True]), barsabove=True,fmt='gs', ms=12, mec='k')
#Sharma 2016, Eagle simulation
lowerlim_f_esc=ax_escfrac.errorbar([3.0, 6.0, 8.0], np.array([0.09, 0.12, 0.15]),  yerr=[0.035, 0.04, 0.04], lolims=np.array([True]), barsabove=True,fmt='gs', ms=12, mec='g', mfc='none')
#error=ax_escfrac.errorbar(z_esc, f_esc, error=[low_err, upper_error], barsabove=True,fmt='.r')
'''

#ax_escfrac.set_ylim(1, 5)
ax_escfrac.set_xlabel('${\\rm redshift}~ (z)$', fontsize=20)
ax_escfrac.set_ylabel('$f_{\\rm esc}$', fontsize=20, labelpad=15)
#ax_escfrac.set_yscale('log')

### for major minor ticks location ###
majorLocator   = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
ax_escfrac.xaxis.set_major_locator(majorLocator)
minorLocator   = MultipleLocator(0.05)
#ax_escfrac.yaxis.set_minor_locator(minorLocator)
minorLocator   = MultipleLocator(1)
ax_escfrac.xaxis.set_minor_locator(minorLocator)


#minorLocator   = MultipleLocator(0.05)
#ax.yaxis.set_minor_locator(minorLocator)

ax_escfrac.set_xticks([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
ax_escfrac.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax_escfrac.tick_params(axis='both', which='both', labelsize=20)


ax_escfrac.set_yticklabels(ax_escfrac.get_yticks(), weight='bold')
ax_escfrac.set_xticklabels(ax_escfrac.get_xticks(), weight='bold')

'''
plt.subplots_adjust(left=0.15, bottom=0.13, right=0.9, top=0.97)
cax = plt.axes([0.95, 0.13, 0.01, 0.84])
plt.colorbar(cbar, cax=cax)
'''
font = font_manager.FontProperties(weight='black',
                                   style='normal', size=12)
ax_escfrac.legend(loc='upper left', prop=font)

plt.savefig('esc_frac.pdf', format='pdf', dpi=80)
plt.show()
