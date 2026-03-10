from joblib import Parallel, delayed
import multiprocessing
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from scipy import interpolate
import sys
sys.path.append('/home/atridebchatterjee/reion_GPR/reion_21cm_code/Reionization_code') 
import running_reion
from running_reion import Redshift_Evolution

from warnings import filterwarnings
filterwarnings('ignore')

from pylab import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

Lymanlimitdatafile='/home/atridebchatterjee/reion_GPR/reion_21cm_code/ObsData/Lyman_limit.dat'
gammadatafile='/home/atridebchatterjee/reion_GPR/reion_21cm_code/ObsData/gamma_data_all_combined.dat'

sample=np.loadtxt('/home/atridebchatterjee/reion_GPR/chain_storage/chains_Infclipped_Aug18/Infclipped.txt')
#print(f'num of rows: {sample.shape[0]}\n num of columns: {sample.shape[1]}')
Hubble, ombh2, omch2, sigma8, ns, esc_pop_III= 67.70, 0.0223, 0.12, 0.81, 0.96, 0.0


sample_plot=sample[np.random.randint(len(sample), size=5000)][:,1:-1]

def likelihood(LL, gamma_PI, tau, Z):
	gamma_PI_interpolate=interpolate.interp1d(Z, gamma_PI, fill_value='extrapolate')
	gamma_model=gamma_PI_interpolate(redshiftGamma)
		
	LL_interpolate=interpolate.interp1d(Z, LL, fill_value='extrapolate')
	LL_model=LL_interpolate(lymanRedshift)
	loglikeLymanSystem=-0.5 * np.sum((lymanLimitData - LL_model) ** 2 /lymanError**2)
	loglikeGamma=-0.5 * np.sum((Gamma_log - np.log10(gamma_model)) ** 2 /error_log**2)
	logliketau=-0.5 * (tau - 0.054)**2 / (0.007**2)
	loglike_tot=loglikeLymanSystem+loglikeGamma+logliketau
	return loglike_tot


def ObsData(gammadatafile, Lymanlimitdatafile):
	redshiftGamma, Gamma, Gamma_max, Gamma_min=np.loadtxt(gammadatafile, usecols=(0,1,2,3), unpack=True)        
	Gamma_log=np.log10(Gamma) 
	error_up=Gamma_max-Gamma
	error_log=error_up/Gamma
	lymanRedshift, lymanLimitData, lymanError=np.loadtxt(Lymanlimitdatafile, usecols=(0,1,2), unpack=True)
	return lymanRedshift, redshiftGamma, lymanLimitData, lymanError, Gamma_log, error_log, error_up



def data_for_plotting(a0, a1, a2, a3, a4, a5, lambda_0, tau_chain):
	eps_param= np.array([a0, a1, a2, a3, a4, a5])
	GPR_interpolator, esc_popII=reion_out(Z, eps_param)
	Redshift, QHII, tau, dnlldz, gamma_PI, QHII5point8=Redshift_Evolution(Hubble, ombh2, omch2, sigma8, ns, esc_pop_III, GPR_interpolator, lambda_0).quntity_for_MCMC()
	#print(QHII5point8)
	
	if np.allclose(tau, tau_chain) and tau<0.06:
		#count=count+1
		#samples[count]=np.array([a0, a1, a2, a3, a4, a5])
		#weghts_array[count]=likelihood
		Z_arr.append(Redshift)
		QHII_arr.append(QHII)
		dnLLdz_arr.append(dnlldz)
		gamma_PI_arr.append(gamma_PI)
		esc_PopII_array.append(esc_popII)
		return Z_arr, QHII_arr, dnLLdz_arr, gamma_PI_arr, esc_PopII_array #tau_array
	else:
		pass
	#return Redshift, QHII, tau, dnlldz, gamma_PI, likelihood

def reion_out(Z, eps_param):
	#esc_PopII_redshift=np.array([3, 8, 13, 18]).reshape(-1, 1)
	esc_PopII_redshift=np.array([3, 6, 9, 12, 15, 18]).reshape(-1, 1)
	esc_Pop_II_val=np.array(eps_param).reshape(-1, 1)
	GPR_interpolator=gaussian_regress_process(esc_Pop_II_val, esc_PopII_redshift)
	mean_prediction, std_prediction = GPR_interpolator.predict(Z.reshape(-1,1), return_std=True)
	print(GPR_interpolator.score(esc_PopII_redshift, esc_Pop_II_val))
	#print(eps_param)
	#print(mean_prediction)
	return GPR_interpolator, mean_prediction



def gaussian_regress_process(esc_Pop_II_array, esc_PopII_redshift, n_restarts_optimizer=10 ):
	kernel = 1 * RBF(8.0)
	gaussian_process = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=n_restarts_optimizer)
	gaussian_process.fit(esc_PopII_redshift, esc_Pop_II_array)
	return gaussian_process

lymanRedshift, redshiftGamma, lymanLimitData, lymanError, Gamma_log, error_log, error_up=ObsData(gammadatafile, Lymanlimitdatafile)


zstart = 2.0
zend = 25.0 #To change this, you have to change also change this in CoredataGenerator.py, start.py, and running_reion.py (the last two files are inside the folder Reioniation_code)
dz_step=0.2
n=int(abs(((zend-zstart)/dz_step)))+1
Z=np.linspace(zend,zstart,n)


num_cores = multiprocessing.cpu_count()	


Z_arr, QHII_arr, tau_arr, dnLLdz_arr, gamma_PI_arr, esc_PopII_array=[], [], [], [], [], []



best_fit=sample[np.argmin(sample[:,1]), 1:-1]
print(f'best fit params:{best_fit}')	
lambda_0_best=best_fit[-1]
GPR_best, PopII_esc_best=reion_out(Z, best_fit[1:7])

plt.plot(Z,PopII_esc_best)
plt.show()












fig = plt.figure(figsize=(15,7))
gs = gridspec.GridSpec(nrows=2, 
                       ncols=2, 
                       figure=fig, 
                       wspace=0.6,
                       hspace=0.0)
fig.subplots_adjust(left=0.1, right=0.9)


ax_QHII = fig.add_subplot(gs[0, 0])


ax_dnLLdz = fig.add_subplot(gs[0, 1],sharex=ax_QHII)


ax_Gamma_PI_H=fig.add_subplot(gs[1, 0],sharex=ax_QHII)


ax_PopII_esc=fig.add_subplot(gs[1, 1],sharex=ax_QHII)



'''
ax_text=fig.add_subplot(gs[0, 2],frameon=False)
ax_text.text(0.1, 0.6, r'$\mathbf{-}$ \textbf{best fit} ($\mathbf{CMB+QSO}$)',color='red', fontsize=20, fontweight="bold",va="center", ha="center")
ax_text.text(0.3, 0.4, r'$\mathbf{-}$ \textbf{best fit} ($\mathbf{CMB+QSO+EDGES}$)',color='blue', fontsize=20, fontweight="bold",va="center", ha="center")
ax_text.text(0.4, 0.2, r'$\mathbf{-}$  \textbf{ random samples} ($\mathbf{CMB+QSO}$)',color='salmon', fontsize=20, fontweight="bold",va="center", ha="center")
ax_text.text(0.5, 0.0, r'$\mathbf{-}$  \textbf{ random samples} ($\mathbf{CMB+QSO+EDGES}$) ',color='cyan', fontsize=20, fontweight="bold",va="center", ha="center")
'''




ax_QHII.set_ylabel(r'$\mathbf{Q_{HII}}$')
ax_QHII.set_xlim(2.0, 9.0)
ax_QHII.set_xlabel(r'$\mathbf{redshift (z)}')
ax_QHII.axes.get_xaxis().set_visible(False)
ax_QHII.text(5, 0.8, '(a)',color='k', fontsize=20, fontweight="bold")



ax_Gamma_PI_H.set_ylabel(r'$\mathbf{\Gamma_{PI}^{HII}/10^{-12} \,sec^{-1}}$')
ax_Gamma_PI_H.set_ylim(10**-2,8) 
ax_Gamma_PI_H.set_yticks([10**-2,1.0])
ax_Gamma_PI_H.set_xlim(2.0, 8.0)
ax_Gamma_PI_H.set_xlabel(r'$\mathbf{redshift (z)}$')
ax_Gamma_PI_H.set_yscale('log') #nonposy='clip'
ax_Gamma_PI_H.text(5, 5, '(b)',color='k', fontsize=20, fontweight="bold")

ax_dnLLdz.set_ylabel(r'$\mathbf{dN_{LL}/dz}$')
ax_dnLLdz.set_xlim(2.0, 8.0)
ax_dnLLdz.set_ylim(0.0,15.0)
ax_dnLLdz.set_yticks([0,4.0,8.0,8.0])
ax_dnLLdz.axes.get_xaxis().set_visible(False)
ax_dnLLdz.set_xlabel(r'$\mathbf{redshift (z)}$')
ax_dnLLdz.text(5, 13, '(c)',color='k', fontsize=20, fontweight="bold")


#ax_PopII_esc

ax_PopII_esc.set_ylabel(r'$\mathbf{f_{esc, II}}$')
ax_PopII_esc.set_xlim(2.0, 8.0)
ax_PopII_esc.set_ylim(0.0, 0.015)
#ax_PopII_esc.set_yticks([0,4.0,8.0,8.0])
ax_PopII_esc.set_xlabel(r'$\mathbf{redshift (z)}$')
ax_PopII_esc.text(5, 0.01, '(d)',color='k', fontsize=20, fontweight="bold")




lineObsGamma=ax_Gamma_PI_H.errorbar(redshiftGamma,10**Gamma_log,yerr=error_up,fmt='.k',barsabove=True,zorder=10)
lineObsLL=ax_dnLLdz.errorbar(lymanRedshift,lymanLimitData,yerr=lymanError,fmt='.k', barsabove=True,zorder=10)

################################################################################################################################


Z_best, QH_best, tau_best, dnLLdz_best, gamma_PI_H_best, QHII5point8=Redshift_Evolution(Hubble, ombh2, omch2, sigma8, ns, esc_pop_III, GPR_best, best_fit[-2]).quntity_for_MCMC()

best_likelihood=likelihood(dnLLdz_best, gamma_PI_H_best, tau_best, Z_best)
print(f'best fit now:{np.abs(best_likelihood)}\nbest_fit MCMC: {best_fit[0]}')
if np.allclose(np.abs(best_likelihood), best_fit[0], atol=1e-2):
	print('best_fit likelihood is matching')
else:
	print('best fit likelihood is not matching so stopping the run')
	exit()
	
linebestQH,=ax_QHII.plot(Z_best, QH_best, color='red', lw=2, zorder=10)
linesbestGamma_H,=ax_Gamma_PI_H.plot(Z_best, gamma_PI_H_best, color='red',lw=2, zorder=10)
linebestdnLLdz,=ax_dnLLdz.plot(Z_best, dnLLdz_best,color='red',lw=2, zorder=10)
linbest_PopII,=ax_PopII_esc.plot(Z, PopII_esc_best, color='red', lw=2, zorder=10)


Observables = Parallel(n_jobs=num_cores)(delayed(data_for_plotting)(*params[1:]) for params in sample_plot)


Sample_plottingdata = np.squeeze(np.array(list(filter(lambda item: item is not None, Observables))))
if Sample_plottingdata.ndim==2:
	print('Too few sample, so exiting')
	exit()
#print(f'Length of Sample_plottingdata {len(Sample_plottingdata)}')
#print(f'shape of None removed array: {Sample_plottingdata.shape}')

for observed in Sample_plottingdata:
	lineQH,=ax_QHII.plot(np.array(observed[0]),np.array(observed[1]),c='salmon',alpha=0.1, zorder=5)
	linednLLdz,=ax_dnLLdz.plot(np.array(observed[0]),np.array(observed[2]),c='salmon',alpha=0.1,zorder=5)
		
	lineGamma_H,=ax_Gamma_PI_H.plot(np.array(observed[0]),np.array(observed[3]),c='salmon',alpha=0.1,zorder=5)
	lienePopII_esc,=ax_PopII_esc.plot(np.array(observed[0]),np.array(observed[4]),c='salmon',alpha=0.1,zorder=5)


'''
fig.savefig('./redshift_evolution_together_Aug18.pdf')
plt.draw()
plt.show()
'''
