import numpy as np
import csv
import os

file_root = 'GPR_MCMC'
burnin = 0.3
#data_dir = '/home/atridebchatterjee/reion_GPR/chain_storage/chains_Aug18/'
data_dir='/home/atridebchatterjee/reion_GPR/chain_storage/4_param_run25_09_2022/'

infclip_dir=os.makedirs('/home/atridebchatterjee/reion_GPR/chain_storage/4_param_run25_09_2022_Infclipped', exist_ok=True)

paramname_file = data_dir + file_root+ '.paramnames'
with open(paramname_file) as f:
   reader = csv.reader(f, delimiter=" ")
   param_name_arr=np.asarray(list(zip(*reader))[0])

#param_latex_name_arr=['$H_{0}$', '$\\Omega_b\\,h^2$','$\\Omega_c\\,h^2$','$A_s$','$n_s$','$f^{II}_{esc}$', '$\\lambda_{0}$', '$m_{x}$', '$f_{X}$', '$f_{\\alpha}$', '$\\tau$', '$\sigma_8$']
#param_latex_name_arr=['$H_{0}$', '$\\Omega_b\\,h^2$','$\\Omega_c\\,h^2$','$A_s$','$n_s$','$f^{II}_{esc}$', '$\\lambda_{0}$', '$m_{x}$', '$\\tau$', '$\sigma_8$']

ndim=int(len(param_name_arr))
nwalkers=20


for k in range(nwalkers):
    filename = data_dir + '/' + file_root + "_"+str(k+1)+".txt"
    all_data =  np.loadtxt(filename,unpack=True)

    if k==0:
        mcmc_steps = all_data.shape[1]
        burnin = int(mcmc_steps * burnin)
        lnprob = np.zeros([nwalkers, mcmc_steps-burnin])
        chains = np.zeros([nwalkers, mcmc_steps-burnin, ndim])
    

    lnprob[k] = all_data[1,burnin:]
    chains[k,:,:] = all_data[2:ndim+2, burnin:].T

#chains[:,:,7]=1.0/chains[:,:,7]
#burnin = int(mcmc_steps * burnin)

lnprob= lnprob[:,:]
samples = chains[:, :, :]

chains_2d=samples.reshape(-1,ndim)
lnprob_flat=lnprob.flatten()

print(len(lnprob_flat))

print ('best-fit likelihood = ', np.min(lnprob_flat))
param_best = chains_2d[np.argmax(lnprob_flat),:]

print(param_best)

index=np.logical_not(np.isinf(lnprob_flat))

lnprobInfClipped=lnprob_flat[index]

#print(chains_2d.shape, len(lnprob_flat))
chains_2d=chains_2d[index,:]

np.savetxt('/home/atridebchatterjee/reion_GPR/chain_storage/4_param_run25_09_2022_Infclipped/Infclipped.txt', np.c_[np.ones(len(lnprobInfClipped)), lnprobInfClipped, chains_2d[:,] ], fmt='%1.4e')

'''
import corner
levels = 1.0 - np.exp(-0.5 * np.array([1.0]) ** 2)
#levels = np.exp(-0.5 * np.arange(1.0, 2.0, 3.0) ** 2)
#print levels


all_latex_labels_arr = param_latex_name_arr
truth_arr = param_best



fig = corner.corner(chains_2d, labels=param_latex_name_arr, 
                    truths=truth_arr, 
                    plot_contours=True, quantiles=[0.16, 0.5, 0.84], show_titles=True)
fig.savefig(file_root + "_triangle_10thFeb.pdf")
print ('done plotting triangles')
'''
