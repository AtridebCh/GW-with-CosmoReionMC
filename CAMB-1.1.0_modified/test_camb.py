from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import camb
import time
from scipy.optimize import curve_fit

A=time.time()
#data_exter=np.loadtxt('QHe.dat',usecols=(0,1),unpack=True)
Z_array,x_HI=np.loadtxt('QHe.dat',usecols=(0,1),unpack=True)
#x_HI=data_exter[:,1]
#camb_path='/home/atrideb/CAMB-1.1.0_modified/'

#print(len(x_HI))
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency

#plt.plot(,np.flip(extra_x_HI))
#plt.show()
#pars=camb.read_ini('/home/atrideb/CAMB-1.1.0_modified/inifiles/planck_2018.ini')
def tanh_func(z,z_re,delta_z):
	#delta_z=3.0
	y=(1.0+z)**1.5	
	y_re=(1.0+z_re)**1.5
	delta_y=1.5*delta_z*(1+z)**0.5
	f=1.08
	result=f/2.0*(1+np.tanh((y-y_re)/delta_y))
	index=result>1.0
	result[index]=1.0
	return 1-result
popt, pcov = curve_fit(tanh_func, Z_array, x_HI,bounds=(0, [12.0,10.0]))
print(popt[0])
pars.redshift_array_external=np.flip(Z_array)
pars.reionization_history=np.flip(x_HI)

pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0)
pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0)
pars.Reion.redshift=popt[0]
pars.Reion.delta_redshift=0.6
#print('here')
data= camb.get_background(pars)
redshift = 10**(np.linspace(0, 2,1000))
back_ev = data.get_background_redshift_evolution(redshift, ['x_e', 'visibility'])
print(time.time()-A)
plt.plot(redshift,back_ev['x_e'],'ro')
plt.plot(Z_array,x_HI,'bo')
plt.show()

results=camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print(totCL.shape)
#Python CL arrays are all zero based (starting at L=0), Note L=0,1 entries will be zero by default.
#The different CL are always in the order TT, EE, BB, TE (with BB=0 for unlensed scalar results).

ls = np.arange(totCL.shape[0])
fig, ax = plt.subplots(2,2, figsize = (12,12))
ax[0,0].plot(ls,totCL[:,0], color='k')
ax[0,0].set_ylim(0.0,6000.0)
ax[0,0].plot(ls,unlensedCL[:,0], color='r')
ax[0,0].set_title('TT')
ax[0,1].plot(ls[2:], 1-unlensedCL[2:,0]/totCL[2:,0]);
ax[0,1].set_title(r'$\Delta TT$')
ax[1,0].plot(ls,totCL[:,1], color='k')
ax[1,0].plot(ls,unlensedCL[:,1], color='r')
ax[1,0].set_title(r'$EE$')
ax[1,1].plot(ls,totCL[:,3], color='k')
ax[1,1].plot(ls,unlensedCL[:,3], color='r')
for ax in ax.reshape(-1): ax.set_xlim([2,2500]);
plt.show()
#There are functions get plot evolution of variables, e.g. for the background as a function of conformal time:

