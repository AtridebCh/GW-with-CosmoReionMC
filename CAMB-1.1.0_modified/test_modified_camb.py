import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import curve_fit
import camb

start_time=time.time()

#Currently use the precomputed reionization history from QHII.dat file for the input and the optical depth is calculated from the reionization history written in the QHII.dat file

Redshift, QHII=np.loadtxt('QHII.dat', unpack=True)

optical_depth=0.054

pars = camb.CAMBparams()
pars.ReionExternal=True
pars.set_cosmology(H0=70.5, ombh2=0.0226, omch2=0.1131, mnu=0.06, omk=0, tau=optical_depth)
pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0)
pars.Reion.UsePCReion=True
pars.length_array = len(Redshift) #length of redshift array, best if you don't change
pars.redshift_array_external=Redshift
pars.reionization_history=QHII
#print(pars)

data= camb.get_background(pars)
redshift_1= np.linspace(0.0, 10.0, 101)
back_ev = data.get_background_redshift_evolution(redshift_1, ['x_e', 'visibility'])

#print('time to calculate cls',time.time()-start_time)


plt.plot(redshift_1, back_ev['x_e'],'ro', label ='x_e')
plt.plot(Redshift, QHII, label='QHII')
plt.legend()
plt.savefig('plot_xe_QHII.pdf')
plt.show()

results=camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']


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
plt.savefig('cls_with_external_reion.pdf')
plt.show()

plt.show()
