import sys

from getdist import loadMCSamples
from fgivenx import plot_contours, samples_from_getdist_chains, plot_lines, plot_dkl
import matplotlib.pyplot as plt

from pylab import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.patches import Rectangle

import matplotlib.font_manager as font_manager

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

#asymmetric_error = [Gamma-Gamma_min,Gamma_max-Gamma]

f_esc_data=np.loadtxt('fesc.dat')
z_esc, f_esc, low_err, upper_error = f_esc_data[:,0], f_esc_data[:,3]*100.0, f_esc_data[:,4]*100.0, f_esc_data[:,5]*100.0

f_esc_upperlim_data=np.loadtxt('f_esc_upper_limit.dat')
z, z_low, z_up, f_esc_upper = f_esc_upperlim_data[:,0], f_esc_upperlim_data[:,1], f_esc_upperlim_data[:,2], f_esc_upperlim_data[:,3]*100.0

zerr=[z-z_low, z_up-z]

def reion_out(Z, theta, count=3):
	eps_param = theta[0:np.size(theta)-1]
	lambda_0 = theta[-1]
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
	return epsilon_out / 0.04 #for plotting
	
def LL_out(Z, theta):
	mean_prediction, GPR_interpolator, lambda_0 = reion_out(Z, theta)	
	if np.all(mean_prediction>0.0):
		Redshift, QHII, tau, dnlldz, gamma_PI, QHII5point8, x_HI=Redshift_Evolution(Planck['H0'], Planck['ombh2'], Planck['omch2'], Planck['sigma8'], Planck['ns'], 
Astro_dict['esc_pop_III'], GPR_interpolator, lambda_0 ).quntity_for_MCMC()
	#print(dnlldz[Redshift<8.2].shape)
	return dnlldz[Redshift<12.2]
	
	
def Gamma_out(Z, theta):
	mean_prediction, GPR_interpolator, lambda_0 = reion_out(Z, theta)	
	if np.all(mean_prediction>0.0):
		Redshift, QHII, tau, dnlldz, gamma_PI, QHII5point8, x_HI=Redshift_Evolution(Planck['H0'], Planck['ombh2'], Planck['omch2'], Planck['sigma8'], Planck['ns'], 
Astro_dict['esc_pop_III'], GPR_interpolator, lambda_0 ).quntity_for_MCMC()
	return gamma_PI[Redshift<12.2]
	
def xHI_out(Z, theta):
	mean_prediction, GPR_interpolator, lambda_0 = reion_out(Z, theta)	
	if np.all(mean_prediction>0.0):
		Redshift, QHII, tau, dnlldz, gamma_PI, QHII5point8, x_HI=Redshift_Evolution(Planck['H0'], Planck['ombh2'], Planck['omch2'], Planck['sigma8'], Planck['ns'], 
Astro_dict['esc_pop_III'], GPR_interpolator, lambda_0 ).quntity_for_MCMC()
	return x_HI[Redshift<12.2]

def QHII_out(Z, theta):
	mean_prediction, GPR_interpolator, lambda_0 = reion_out(Z, theta)	
	if np.all(mean_prediction>0.0):
		Redshift, QHII, tau, dnlldz, gamma_PI, QHII5point8, x_HI=Redshift_Evolution(Planck['H0'], Planck['ombh2'], Planck['omch2'], Planck['sigma8'], Planck['ns'], 
Astro_dict['esc_pop_III'], GPR_interpolator, lambda_0 ).quntity_for_MCMC()
	QHII[QHII>1.0] = 1.0
	return QHII[Redshift<20.2]

np.random.seed(7)	
random_index=np.random.randint(len(samples_fgx), size=20) # change "size" to make plots more densed
random_index1 = np.random.randint(len(samples_fgx), size=20000)

HII=[]












fig=plt.figure(figsize=(7,5.7))

gs = gridspec.GridSpec(nrows=1,
                       ncols=1,
                       figure=fig,
                       #width_ratios= [1, 1, 1],
                       #height_ratios=[1, 1, 1],
                       wspace=0.0,
                       hspace=0.0)


ax_QHII = fig.add_subplot(gs[0, 0])
ax_QHII.set_ylim([0.0, 1.1])

for param in samples_fgx[random_index]:
	#plt.plot(Z_plot, QHII_out(Z, param))
	HII.append(QHII_out(Z, param))
QHII_mean=np.mean(np.asarray(HII), axis=0)

QHII_std = np.std(np.asarray(HII), axis=0)
QHII_upper = QHII_mean+QHII_std
QHII_upper[QHII_upper>1.0]=1.0
QHII_upper2 = QHII_mean+2*QHII_std
QHII_upper2[QHII_upper2>1.0]=1.0

plt.fill_between(Z, QHII_upper, QHII_mean-QHII_std, color='b')
plt.fill_between(Z, QHII_upper2, QHII_upper, alpha=0.5, color='b')
plt.fill_between(Z, QHII_mean-2*QHII_std, QHII_mean-QHII_std, alpha=0.5, color= 'b')

################## QHII plot ###############
#cbar = plot_contours(QHII_out, Z_plot, samples_fgx[random_index1], weights=weights_fgx[random_index1], ax=ax_QHII, ntrim=50, colors='Blues_r',lines=False,contour_line_levels=[1,2],contour_color_levels=[0,1,2]) #,contour_line_levels=[1,2],contour_color_levels=[0,1,2])
ax_QHII.set_ylabel('$Q_{\\rm HII}$', fontsize=20, labelpad=15)
ax_QHII.set_xlabel('${\\rm redshift}~ (z)$', fontsize=20)
#ax_xHI.set_yscale('log')


xmin = 2.0
xmax = 20.0


# DATA POINTS ############

#McGreer+15
zdata=[5.58,5.87,6.1]
taudata=[0.96,0.94,0.62]
errorlow=[0.05,0.05,0.2]
errorupp=[0.0,0.0,0.0]
ax_QHII.errorbar(zdata,taudata,yerr=[errorlow, errorupp], label='McGreer+15', marker='^', ms=11, color='gray', ls='', mec='k')
ax_QHII.errorbar(zdata, taudata, marker='', ms=4, yerr=[0.07,0.07,0.07], lolims=True, ls='none', capsize=3, color='gray', mec='k')

#GRB
zdata=[5.95,6.3]
taudata=[0.9+0.02,0.5+0.01]
a=[[5.95,0.9],[6.3,0.5]]

ax_QHII.plot(*zip(*a), label='Chornock+13 ($z=5.9$); \n Totani+06 ($z=6.3$)', marker='^', ms=12, color='r', ls='', mec='r', mfc='none')
ax_QHII.errorbar(zdata, taudata, marker='', ms=4, yerr=[0.07, 0.07], lolims=True, ls='none', capsize=3, color='r')

#Quasar near zone (Schroeder+13 & Bolton+11)
zdata=[6.3,7.1]
taudata=[0.9,0.9]
a=[[6.3,0.9],[7.1,0.9]]
ax_QHII.plot(*zip(*a), label='Schroeder+13 ($z=6.3$); \n Bolton+11 ($z=7.1$)', marker='d', ms=12, color='darkorange', ls='', mec='k', mfc='darkorange')
ax_QHII.errorbar(zdata, taudata, marker='', ms=1, yerr=0.1, uplims=True, ls='none', capsize=3, color='darkorange', mec='k')

#Lyalpha emision fraction (Schenker+14)
zdata=[8.0]
taudata=[1-0.65]
a=[[8.0,1-0.65]]
ax_QHII.plot(*zip(*a), label='Schenker+14', marker='s', ms=10, color='g', ls='', mec='w', mfc='g')
ax_QHII.errorbar(zdata, taudata, marker='', ms=1, yerr=0.1, uplims=True, ls='none', capsize=3, color='g', mec='w')
zdata=[7.0]
taudata=[1-0.34]
errorlow=[0.12]
errorupp=[0.09]
ax_QHII.errorbar(zdata, taudata, marker='s', ms=10, yerr=[errorlow,errorupp], ls='', mfc='g', capsize=3, color='g', mec='w')

#Mason+18, LBG
zdata=[6.95]
taudata=[0.41]
errorlow=[0.15]
errorupp=[0.11]
a=[[6.95,0.41]]
ax_QHII.plot(*zip(*a), label='Mason+18', marker='p', ms=12, color='darkmagenta', ls='', mec='w', mfc='darkmagenta')
ax_QHII.errorbar(zdata, taudata, marker='', ms=1, yerr=[errorlow,errorupp], ls='', mfc='darkmagenta', capsize=3, color='darkmagenta', mec='w')

# Davies+18
# https://arxiv.org/abs/1802.06066
zdata=[7.09, 7.54]
taudata=[0.52,0.4]
errorlow=[0.26,0.23]
errorupp=[0.26,0.20]
a=[[7.09,0.52],[7.54,0.4]]
ax_QHII.plot(*zip(*a), label='Davies+18', marker='h', ms=12, color='k', ls='', mec='w', mfc='k')
ax_QHII.errorbar(zdata, taudata, marker='', ms=1, yerr=[errorlow,errorupp], ls='', mfc='k', capsize=3, color='k', mec='w')

# Greig+18
# https://arxiv.org/abs/1807.01593
zdata=[7.54]
taudata=[0.79]
errorlow=[0.19]
errorupp=[0.17]
a=[[7.54,0.79]]
ax_QHII.plot(*zip(*a), label='Greig+18', marker='o', ms=12, color='deeppink', ls='', mec='w', mfc='deeppink')
ax_QHII.errorbar(zdata, taudata, marker='', ms=1, yerr=[errorlow,errorupp], ls='', mfc='deeppink', capsize=3, color='deeppink', mec='w')

# Ďurovčíková+20
# https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.4256D/abstract
zdata=[7.0851,  7.5413]
taudata=[0.75,0.4]
errorlow=[0.05,0.11]
errorupp=[0.05,0.11]
a=[[7.0851,0.75],[7.5413,0.4]]
ax_QHII.plot(*zip(*a), label='Ďurovčíková+20', marker='o', ms=12, color='red', ls='', mec='red', mfc='none')
ax_QHII.errorbar(zdata, taudata, marker='', ms=1, yerr=[errorlow,errorupp], ls='', mfc='red', capsize=3, color='red')


# Jin+23
# https://ui.adsabs.harvard.edu/abs/2023ApJ...942...59J/abstract
zdata=[6.3, 6.5, 6.7]
taudata=[ 1-0.79+0.017, 1-0.87+0.017, 1-0.94+0.017]

a=[[6.3,1-0.79],[6.5,1-0.87],[6.7,1-0.94]]
ax_QHII.plot(*zip(*a), label='Jin+23', marker='s', ms=10, color='darkcyan', ls='', mec='darkcyan', mfc='none')
ax_QHII.errorbar(zdata, taudata, marker='', ms=1, yerr=[0.07, 0.07, 0.07], lolims=True, ls='none', capsize=3, color='darkcyan', mec='darkcyan', mfc='none')

##########################
ax_QHII.plot(Z_plot, QHII_mean, ls='-', color='darkorange', lw=3, label='mean')


### for major minor ticks location ###
majorLocator   = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
ax_QHII.xaxis.set_major_locator(majorLocator)
minorLocator   = MultipleLocator(0.1)
ax_QHII.yaxis.set_minor_locator(minorLocator)
minorLocator   = MultipleLocator(1)
ax_QHII.xaxis.set_minor_locator(minorLocator)


#minorLocator   = MultipleLocator(0.05)
#ax.yaxis.set_minor_locator(minorLocator)

ax_QHII.set_xticks([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
ax_QHII.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax_QHII.tick_params(axis='both', which='both', labelsize=20)


ax_QHII.set_yticklabels(ax_QHII.get_yticks(), weight='bold')
ax_QHII.set_xticklabels(ax_QHII.get_xticks(), weight='bold')

'''
plt.subplots_adjust(left=0.15, bottom=0.13, right=0.9, top=0.97)
cax = plt.axes([0.95, 0.13, 0.01, 0.84])
plt.colorbar(cbar, cax=cax)
'''
x1 = np.array([2, 5.8])
y1 = np.array([1, 1])

ax_QHII.plot(x1, y1, color='b', ls='solid')

font = font_manager.FontProperties(weight='black',
                                   style='normal', size=12)
ax_QHII.legend(loc='upper right', prop=font)

plt.savefig('QHII.pdf', format='pdf', dpi=80)
plt.show()


