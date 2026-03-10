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

import numpy as np

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF

from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning

Planck={'H0':67.70, 'ombh2': 0.0223, 'omch2': 0.12, 'ns':0.96, 'sigma8':0.81} 
Astro_dict={'esc_pop_III':0.0}

sys.path.append('/home/atridebchatterjee/reion_GPR/reion_21cm_code/Reionization_code')
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



lymanRedshift, lymanLimitData, lyman_error=np.loadtxt('/home/atridebchatterjee/parameter_estimation_reionization/ObsData/Lyman_limit.dat',usecols=(0,1,2),unpack=True)
redshiftGamma, Gamma, Gamma_max, Gamma_min=np.loadtxt('/home/atridebchatterjee/parameter_estimation_reionization/ObsData/gamma_data_all_combined.dat',skiprows=1,usecols=(0,1,2,3),unpack=True)

asymmetric_error = [Gamma-Gamma_min,Gamma_max-Gamma]

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

Z_plot=np.linspace(10.0, 2.0, int((10.0-2.0)/0.2+1))

file_root = 'Infclipped'
burnin = 0.0



#params = ['epsilon_a0','epsilon_a1', 'epsilon_a2','epsilon_a3', 'epsilon_a4','epsilon_a5', 'lambda_0']
params = ['epsilon_a0','epsilon_a1', 'epsilon_a2','epsilon_a3', 'lambda_0']
samples_fgx, weights_fgx = samples_from_getdist_chains(params, '/home/atridebchatterjee/reion_GPR/chain_storage/4_param_run25_09_2022_Infclipped/' + file_root, settings={'ignore_rows':burnin})


name_string='_4_param'

def epsilon_out(Z, theta):
	epsilon_out, GPR_interpolator, lambda_0 = reion_out(Z, theta)
	return epsilon_out*10**3 #for plotting
	
def LL_out(Z, theta):
	mean_prediction, GPR_interpolator, lambda_0 = reion_out(Z, theta)	
	if np.all(mean_prediction>0.0):
		Redshift, QHII, tau, dnlldz, gamma_PI, QHII5point8, x_HI=Redshift_Evolution(Planck['H0'], Planck['ombh2'], Planck['omch2'], Planck['sigma8'], Planck['ns'], 
Astro_dict['esc_pop_III'], GPR_interpolator, lambda_0 ).quntity_for_MCMC()
	#print(dnlldz[Redshift<8.2].shape)
	return dnlldz[Redshift<10.2]
	
	
def Gamma_out(Z, theta):
	mean_prediction, GPR_interpolator, lambda_0 = reion_out(Z, theta)	
	if np.all(mean_prediction>0.0):
		Redshift, QHII, tau, dnlldz, gamma_PI, QHII5point8, x_HI=Redshift_Evolution(Planck['H0'], Planck['ombh2'], Planck['omch2'], Planck['sigma8'], Planck['ns'], 
Astro_dict['esc_pop_III'], GPR_interpolator, lambda_0 ).quntity_for_MCMC()
	return gamma_PI[Redshift<10.2]
	
def xHI_out(Z, theta):
	mean_prediction, GPR_interpolator, lambda_0 = reion_out(Z, theta)	
	if np.all(mean_prediction>0.0):
		Redshift, QHII, tau, dnlldz, gamma_PI, QHII5point8, x_HI=Redshift_Evolution(Planck['H0'], Planck['ombh2'], Planck['omch2'], Planck['sigma8'], Planck['ns'], 
Astro_dict['esc_pop_III'], GPR_interpolator, lambda_0 ).quntity_for_MCMC()
	return x_HI[Redshift<10.2]

def QHII_out(Z, theta):
	mean_prediction, GPR_interpolator, lambda_0 = reion_out(Z, theta)	
	if np.all(mean_prediction>0.0):
		Redshift, QHII, tau, dnlldz, gamma_PI, QHII5point8, x_HI=Redshift_Evolution(Planck['H0'], Planck['ombh2'], Planck['omch2'], Planck['sigma8'], Planck['ns'], 
Astro_dict['esc_pop_III'], GPR_interpolator, lambda_0 ).quntity_for_MCMC()
	return Q_HII[Redshift<10.2]

np.random.seed(7)	
random_index=np.random.randint(len(samples_fgx), size=20) # change "size" to make plots more densed


fig = plt.figure(figsize=(24,7))
gs = gridspec.GridSpec(nrows=1, 
                       ncols=3, 
                       figure=fig, 
                       wspace=0.6,
                       hspace=0.0)
fig.subplots_adjust(left=0.1, right=0.9)

#ax_escfrac = fig.add_subplot(gs[0, 0])
ax_xHI=fig.add_subplot(gs[0, 2])
ax_LL=fig.add_subplot(gs[0, 0])
ax_gamma_PI=fig.add_subplot(gs[0, 1])
#ax_QHII=fig.add_subplot(gs[0, 2])



'''
### esc_frac_plot ######
cbar = plot_contours(epsilon_out, Z_plot, samples_fgx[random_index], weights=weights_fgx[random_index], ax=ax_escfrac, ntrim=100, contour_line_levels=[1,2], colors='Blues_r',lines=False) #,contour_line_levels=[1,2],contour_color_levels=[0,1,2])
ax_escfrac.errorbar(z_esc, f_esc, yerr=[low_err, upper_error], fmt='.k', barsabove=True)

#ax_escfrac.set_ylim(1, 5)
ax_escfrac.set_xlabel(r'redshift $(\mathbf{Z})$')
ax_escfrac.set_ylabel(r'$\mathbf{f_{esc} \times 10^2}$')
ax_escfrac.set_yscale('log')
ax_escfrac.text(1.5, 30, '(a)',color='k', fontsize=20, fontweight="bold")
ax_escfrac.axes.get_xaxis().set_visible(False)
'''
'''
################## QHII plot ###############
cbar_QHII = plot_contours(QHII_out, Z_plot, samples_fgx[random_index], weights=weights_fgx[random_index], ax=ax_xHI, ntrim=100, colors='Blues_r',lines=False) #,contour_line_levels=[1,2],contour_color_levels=[0,1,2])
ax_QHII.set_ylabel(r'$\mathbf{Q_{HII}}$')
ax_QHII.set_xlabel(r'redshift $(\mathbf{Z}$)')
#ax_xHI.set_yscale('log')
'''


### LL plot #########
cbar_LL = plot_contours(LL_out, Z_plot, samples_fgx[random_index], weights=weights_fgx[random_index], ax=ax_LL, ntrim=100,  contour_line_levels=[1,2], colors='Blues_r',lines=False) #,contour_line_levels=[1,2],contour_color_levels=[0,1,2])
#plt.yscale('log')
ax_LL.set_ylabel(r'$\mathbf{dN_{LL}/dz}$')
ax_LL.set_xlabel(r'redshift ($\mathbf{Z}$)')
ax_LL.text(2.5, 10.5, '(a)',color='k', fontsize=20, fontweight="bold")


################xHI plot#################
cbar_xHI = plot_contours(xHI_out, Z_plot, samples_fgx[random_index], weights=weights_fgx[random_index], ax=ax_xHI, ntrim=100, colors='Blues_r',lines=False) #,contour_line_levels=[1,2],contour_color_levels=[0,1,2])
#ax_xHI.set_ylabel(r'$\mathbf{x_{HI}}$')
ax_xHI.set_xlabel(r'redshift $(\mathbf{Z}$)')
#ax_xHI.set_yscale('log')



###############################################


####### Gamma_PI plot #########
cbar_gamma = plot_contours(Gamma_out, Z_plot, samples_fgx[random_index], weights=weights_fgx[random_index], ax=ax_gamma_PI, ntrim=100, colors='Blues_r',lines=False) #,contour_line_levels=[1,2],contour_color_levels=[0,1,2])

ax_gamma_PI.set_ylabel(r'$\mathbf{\Gamma_{PI}/10^{-12} sec^{-1}}$')
ax_gamma_PI.set_xlabel(r'redshift ($\mathbf{Z})$')
ax_gamma_PI.set_yticks([0.5, 1.0, 1.5])
ax_gamma_PI.text(8, 1.7, '(b)',color='k', fontsize=20, fontweight="bold")


#ax_gamma_PI.axes.get_xaxis().set_visible(False)
#########################################################


#upperlim_x_HI=ax_x_HI.errorbar(np.array([5.6,5.8,6.0]),np.array([0.06,0.1,0.20]),xerr=0.06,yerr=0.06,uplims=np.array([True]),marker='o',barsabove=True,fmt='.k',zorder=10)
#plt.errorbar(lymanRedshift,lymanLimitData,yerr=lyman_error,fmt='.k',  barsabove=True,zorder=10)
lineObsGamma=ax_gamma_PI.errorbar(redshiftGamma, Gamma, yerr=asymmetric_error, fmt='.k', barsabove=True)
lineObsLL=ax_LL.errorbar(lymanRedshift, lymanLimitData, yerr=lyman_error, fmt='.k', barsabove=True,zorder=10)
#upperlim_f_esc=ax_escfrac.errorbar(z, f_esc_upper,  xerr=zerr, yerr=0.5,  uplims=np.array([True]), barsabove=True,fmt='.k')





ymin = 0.00002
ymax = 1.05

xmin = 2.0
xmax = 10.0


#ax_xHI.text(0.1, 0.5, '$\mathbf{x_{HI}}$', va='center', rotation='vertical')
ax_xHI.set_ylabel(r'$\mathbf{x_{HI}}$')
ax_xHI.set_yscale('log')
ax_xHI.set_ylim((ymin, 0.09))
ax_xHI.spines['top'].set_visible(False)

divider = make_axes_locatable(ax_xHI)
axLin = divider.append_axes("top", size=2.75, pad=0.0, sharex=ax_xHI)


axLin.set_xscale('linear')
axLin.set_ylim((0.0001, ymax))
ax_xHI.axis([xmin, xmax, ymin, 0.0001])

cbar_xHI = plot_contours(xHI_out, Z_plot, samples_fgx[random_index], weights=weights_fgx[random_index], ax=axLin, ntrim=100, colors='Blues_r',lines=False)

ax_xHI.set_xlabel(r'redshift ($\mathbf{Z}$)')
plt.setp(axLin.get_xticklabels())

#McGreer+15
zdata=[5.58,5.87,6.1]
taudata=[0.09-0.05,0.11-0.055,0.58-0.055]
a=[[5.58,0.09],[5.87,0.11],[6.1,0.58]]

axLin.errorbar(zdata, taudata, marker='', ms=4, yerr=[0.001, 0.001, 0.001], uplims=True, ls='none', capsize=4, color='k', mec='w')
axLin.text(8, 0.9, '(c)',color='k', fontsize=20, fontweight="bold")

#Jin+22
zdata=[6.3, 6.5, 6.7]
taudata=[ 0.79, 0.87, 0.94]

axLin.errorbar(zdata, taudata, marker='', ms=4, yerr=[0.07, 0.07, 0.07], uplims=True, ls='none', capsize=4, color='r', mec='w')

'''
#Fan+06b
zdata=[5.02,5.25,5.45,5.65,5.85,6.1]
taudata=[pow(10.0,-4.2596),pow(10.0,-4.1871),pow(10.0,-4.2076),pow(10.0,-4.0706),pow(10.0,-3.9586),pow(10.0,-3.3468)]
errorlow=[pow(10.0,-4.2596)-pow(10.0,-4.3979),pow(10.0,-4.1871)-pow(10.0,-4.3468),pow(10.0,-4.2076)-pow(10.0,-4.3979),pow(10.0,-4.0706)-pow(10.0,-4.3010),pow(10.0,-3.9586)-pow(10.0,-4.0969),pow(10.0,-3.3468)-pow(10.0,-3.9208)]
errorupp=[pow(10.0,-4.1549)-pow(10.0,-4.2596),pow(10.0,-4.0457)-pow(10.0,-4.1871),pow(10.0,-4.0223)-pow(10.0,-4.2076),pow(10.0,-3.8861)-pow(10.0,-4.0706),pow(10.0,-3.9586)-pow(10.0,-3.9586),pow(10.0,-3.3468)-pow(10.0,-3.3468)]
ax_xHI.errorbar(zdata,taudata,yerr=[errorlow, errorupp], marker='.', ms=8, color='r', ls='', capsize=1, mec='w')
'''
zdata=[5.85,6.1]
taudata=[pow(10.0,-3.9586)+0.0002,pow(10.0,-3.3468)+0.001]
ax_xHI.errorbar(zdata, taudata, marker='', ms=4, yerr=[0.0002,0.001], uplims=True, ls='none', capsize=4, color='r', mec='w')

#GRB
zdata=[5.95,6.3]
taudata=[0.1-0.055,0.5-0.055]
a=[[5.95,0.1],[6.3,0.5]]

axLin.errorbar(zdata, taudata, marker='', ms=4, yerr=0.035, lolims=True, ls='none', capsize=4, color='r')

#Lyalpha emitter (at z=6.6 - Ouchi et al. 2010; at z=7 - Ota et al. 2008)
zdata=[6.6]
taudata=[0.4-0.06]
a=[[6.6,0.4]]

axLin.errorbar(zdata, taudata, marker='', ms=4, yerr=0.04, lolims=True, ls='none', capsize=4, color='r', mec='w')
zdata=[7.0]
taudata=[0.48]
axLin.errorbar(zdata, taudata, marker='s', ms=5, yerr=[0.48-0.32], ls='', mec='r', mfc='none', capsize=8, color='r')

#Lyalpha emision fraction (Schenker+14)
zdata=[8.0]
taudata=[0.65+0.07]
a=[[8.0,0.65]]

axLin.errorbar(zdata, taudata, marker='', ms=4, yerr=0.05, uplims=True, ls='none', capsize=4, color='r', mec='w')
zdata=[7.0]
taudata=[0.34]
errorlow=[0.12]
errorupp=[0.09]
axLin.errorbar(zdata, taudata, marker='s', ms=5, yerr=[errorlow,errorupp], ls='', mfc='k', capsize=4, color='r', mec='w')

#Quasar near zone (Schroeder+13)
zdata=[6.3,7.1]
taudata=[0.1+0.06,0.1+0.06]
a=[[6.3,0.1],[7.1,0.1]]
axLin.errorbar(zdata, taudata, marker='', ms=4, yerr=0.06, uplims=True, ls='none', capsize=4, color='r', mec='w')

axLin.axes.get_xaxis().set_visible(False)

axLin.set_ylabel(r'')





plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = plt.axes([0.85, 0.1, 0.01, 0.8])
plt.colorbar(cbar_LL, cax=cax)


plt.savefig('all_together_xHI.pdf')
plt.show()




