import sys
import numpy as np
import emcee
import matplotlib.pyplot as plt

plt.rcParams['font.size']=8



file_root = 'GPR_MCMC'
burnin = 0.8
data_dir = '/home/atridebchatterjee/reion_GPR/chain_storage/chains_Aug18/'


import csv
filename = data_dir + file_root+ '.paramnames'
with open(filename) as f:
   reader = csv.reader(f, delimiter=" ")
   param_name_arr=np.asarray(list(zip(*reader))[0])

param_latex_name_arr=['a_{0}', 'a_{1}', 'a_{2}', 'a_{3}', 'a_{4}', 'a_{5}', '\\lambda_{0}', '\\tau', 'Q_{HII}']
#print(param_latex_name_arr)
ndim=int(len(param_name_arr))
nwalkers=28

nblobs=0

derived_name_arr=None
    

for k in range(nwalkers):
    filename = data_dir + '/' + file_root + "_"+str(k+1)+".txt"
    all =  np.loadtxt(filename, skiprows=0, unpack=True)

    if k==0:
        mcmc_steps = np.shape(all)[1]
        lnprob = np.zeros([nwalkers, mcmc_steps])
        chains = np.zeros([nwalkers, mcmc_steps, ndim])
        #if derived_name_arr is not None: blob_chains = np.zeros([nwalkers, mcmc_steps, nblobs])
    
    #print(np.shape(all))
    lnprob[k,:] = -all.T[:mcmc_steps,1]
    chains[k,:,:] = all.T[:mcmc_steps,2:ndim+2]
    #blob_chains[k,:,:] = all.T[:mcmc_steps,ndim+2:ndim+nblobs+2]
    

print ('steps = ', mcmc_steps)
#print(np.shape(blobs), np.shape(chains))

samples = chains[:,:,:]
s = samples.shape



if burnin > 1.0:
    burnin = int(burnin)
    if burnin >= mcmc_steps: burnin = 0
else:
    burnin = int(mcmc_steps * burnin)
print ('burnin steps = ', burnin)


auto_corr = np.zeros([nwalkers, ndim])
'''
for j in range(ndim):
    for i in range(nwalkers):
        auto_corr[i,j] = emcee.autocorr.integrated_time(chains[i, burnin:, j])
        #chains=chains.reshape(-1,j)
    print('autocorrelation time for ', param_name_arr[j], ': ', "{0:.2f}".format(np.mean(auto_corr[:,j])))
'''


samples = chains[:, burnin:, :]
s = samples.shape

samples = samples.reshape(s[0] * s[1], s[2])  # Flatten the sample list.
s = samples.shape
print ('number of samples used = ', s[0])


for j in range(ndim-2):
    auto_corr[j] = emcee.autocorr.integrated_time(samples[:, j])
    #print(auto_corr[j])
    print('autocorrelation time for ', param_name_arr[j], ': ', "{0:.2f}".format(np.mean(auto_corr[:, j])))

if derived_name_arr is not None:
    blobs = blob_chains[:,burnin:,:]
    b = blobs.shape
    blobs = blobs.reshape(b[0] * b[1], b[2])
    b = blobs.shape
    #print(np.shape(blobs), np.shape(samples))

lnprob = lnprob[:,burnin:]
l = lnprob.shape
lnprob = lnprob.reshape(l[0] * l[1])




print ('best-fit likelihood = ', np.max(lnprob))
param_best = samples[np.argmax(lnprob),:]
if derived_name_arr is not None: blobs_best = blobs[np.argmax(lnprob),:]
for i in range(ndim): print ("best-fit ", param_name_arr[i], ": ", param_best[i])
if derived_name_arr is not None:
    for i in range(nblobs): print ("best-fit ", derived_name_arr[i], ": ", blobs_best[i])

if derived_name_arr is not None:
    samples_all = np.zeros([s[0], s[1]+b[1]])
    samples_all[:,0:s[1]] = samples[:,:]
    samples_all[:,s[1]:s[1]+b[1]+1] = blobs[:,:]
else:
    samples_all = np.zeros([s[0], s[1]])
    samples_all[:,0:s[1]] = samples[:,:]



median=np.zeros(ndim)

for i in range(ndim):
    param_mcmc = np.percentile(samples_all[:, i], [16, 50, 84])
    median[i]=param_mcmc[1]
    q = np.diff(param_mcmc)
    print ("median +- 1-sigma ", param_name_arr[i], ": ", "{0:.6f}".format(param_mcmc[1]), "{0:.6f}".format(q[0]), "{0:.6f}".format(q[1]))
    
if derived_name_arr is not None:
    for i in range(nblobs):
        param_mcmc = np.percentile(samples_all[:, i+ndim], [16, 50, 84])
        q = np.diff(param_mcmc)
        print ("median +- 1-sigma ", derived_name_arr[i], ": ", "{0:.6f}".format(param_mcmc[1]), "{0:.6f}".format(q[0]), "{0:.6f}".format(q[1]))




for i in range(ndim): 
    param_latex_name_arr[i] = "$"+ param_latex_name_arr[i] +"$"

if derived_name_arr is not None:
    for i in range(nblobs): derived_latex_name_arr[i] = r"$" + derived_latex_name_arr[i] + "$"

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#         CHAIN PLOTTING             #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

plt.clf()
scale = 1.0
fig, axes = plt.subplots(ndim-2, 1, sharex=True, figsize=(scale * 8, scale * (ndim+nblobs)))

for i in range(ndim-2):
    axes[i].plot(np.log10(chains[:, :, i]).T, color="r", alpha=0.2)
    #axes[i].yaxis.set_major_locator(MaxNLocator(5))
    #axes[i].axhline(param_mcmc_arr[i,0], color="#888888", lw=2)
    axes[i].set_ylabel(param_latex_name_arr[i])


for i in range(nblobs):
    axes[ndim+i].plot(blob_chains[:, :, i].T, color="k", alpha=0.2)
    #axes[ndim+i].yaxis.set_major_locator(MaxNLocator(5))
    #axes[ndim+i].axhline(blob_mcmc_arr[i,0], color="#888888", lw=2)
    axes[ndim+i].set_ylabel(derived_latex_name_arr[i])
    

axes[ndim-2+nblobs-1].set_xlabel("step number")
    


fig.tight_layout(h_pad=0.0)
fig.savefig(file_root + "_chain_18thAug.pdf")
print ('done plotting chains')

