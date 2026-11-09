import sys, os

import numpy as np

import emcee
from utils import Params
from Likelihood import Likelihood
from CoredataGenerator import CoreModule
from sampler import MCMCsampler


Lymanlimitdatafile = '/Users/users/achatterjee/Data/GW_with_CosmoReionMC/ReionizationWithF2Py/ObsData/Lyman_limit.dat'
gammadatafile      = '/Users/users/achatterjee/Data/GW_with_CosmoReionMC/ReionizationWithF2Py/ObsData/gamma_data_all_combined.dat'
T_b_Datafile       = '/Users/users/achatterjee/Data/GW_with_CosmoReionMC/global_21cm_signal/brightness_temp.txt'

'''
params = Params(("hubble",                [6.71554e+01, 65, 80, 0.01*(80-65.0)]),
                ("ombh2",                 [2.23010e-02, 0.01, 0.03, 0.01*(0.03-0.01)]),
                ("omch2",                 [1.21812e-01, 0.09, 0.2, 0.01*(0.2-0.09)]),
                ("As",                    [2.1176e-09, 1.8e-9, 2.35e-9, 0.01*(2.35e-09-1.8e-09)]),
                ("ns",                    [9.6419e-01, 0.8, 1.2, 0.01*(1.2-0.8)]),
                ("f_zero",                [0.17, 0.0, 1.0, 0.01]),
                ("alpha_lo",              [0.54, -2.0, 1.0, 0.1]),
                ("alpha_hi",              [0.14, -2.0, 5.0, 0.1]),
                ("alpha_z",               [-0.18, -1.0, 1.0, 0.01]),
                ("esc_II",                [0.039, 0.0, 1.0, 0.1]),
                ("lambda0",               [4.564590e+00, 1.0,10.0,0.01*(6.0-1.0)]),
                ('f_X',                   [1.977400e-01, 0.01, 1.0, 0.001]),
                ('f_alpha',               [1.000560e+00, 0.01, 1.0, 0.001]))
'''
'''
params = Params(("hubble", [6.7266e+01, 65, 80, 0.01*(80-65.0)]), #6.726630e+01
                ("ombh2", [2.23740e-02, 0.01, 0.03, 0.01*(0.03-0.01)]), # 2.2341e-02
                ("omch2", [1.206109e-01, 0.09, 0.2, 0.01*(0.2-0.09)]), #1.206109e-01
                ("As", [2.112747e-09, 1.8e-9, 2.35e-9, 0.01*(2.35e-09-1.8e-09)]), #2.112747e-09
                ("ns", [9.655997e-01, 0.8, 1.2, 0.01*(1.2-0.8)]), #9.655997e-01
                ("f_zero", [0.19, 0.0, 1.0, 0.01]),#0.19 GW depends
                ("alpha_lo", [0.31, -2.0, 2.0, 0.1]),#0.31 GW depends
                ("alpha_hi", [0.11, -2.0, 2.0, 0.1]),#0.11 GW depdends
                ("alpha_z", [0.0, -1.0, 1.0, 0.01]),#0.0
                ("esc_II", [0.035, 0.0, 1.0, 0.001]),#0.035 GW depends
                ("lambda0", [4.98e+00, 1.0, 12.0,0.01*(6.0-1.0)]),#4.98 GW does not depend
                ('f_X', [2.00e-01, 0.01, 1.0, 0.001]),#2.00e-1
                ('f_alpha', [1.0e+00, 0.01, 1.0, 0.001]))#1.0
'''

params = Params(("hubble",                [6.721184e+01, 65, 80, 0.01*(80-65.0)]),
                ("ombh2",                 [2.261010e-02, 0.01, 0.03, 0.01*(0.03-0.01)]),
                ("omch2",                 [1.171812e-01, 0.09, 0.2, 0.01*(0.2-0.09)]),
                ("As",                    [2.121276e-09, 1.8e-9, 2.35e-9, 0.01*(2.35e-09-1.8e-09)]),
                ("ns",                    [9.605819e-01, 0.8, 1.2, 0.01*(1.2-0.8)]),
                ("f_zero",                [0.26,0.0, 1.0, 0.01]),
                ("alpha_lo",              [0.2, 0, 2.0, 0.1]),
                ("alpha_hi",              [0.2, 0, 2.0, 0.1]),
                ("alpha_z",               [0, -4.0, 4.0, 0.1]),
                ("esc_II",                [0.035, 0.0, 1.0, 0.01]),
                ("lambda0",               [5.364590, 1.0, 12.0, 0.1*(12.0-1.0)]))
'''
                ('f_X',                   [1.977400e-01, 0.01, 1.0, 0.001]),
                ('f_alpha',               [1.000560e+00, 0.01, 2.0, 0.001]))
'''


from datetime import date

today = date.today()

# dd/mm/YY
d1 = today.strftime("%m_%Y")

param_latex_name_arr=['H_{0}','\Omega_{b} h^{2}', '\Omega_{c} h^{2}', 'A_s', 'n_s', 'f_{0}', '\\alpha_{lo}', '\\alpha_{hi}', '\\alpha_{z}', 'f_{esc, II}', '\lambda_0', 'f_{X}', 'f_{\\alpha}']



#change the path if needed
parent_dir="/Users/users/achatterjee/Data/GW_with_CosmoReionMC/mcmc/MCMC_sampler/chain_storage_CMB_Reion_GW/" 
MCMC_dir="run_with_merger_rate_6th_July"

chain_storage_path = os.path.join(parent_dir, MCMC_dir)

os.makedirs(chain_storage_path, exist_ok=True)

file_root='test_MCMC'

write_file=Params.write_paramnames_ranges_file(params[:,1], params[:,2], params.keys, param_latex_name_arr, chain_storage_path, file_root)



like = Likelihood(CoreModule, min_param = params[:,1], max_param= params[:,2])
like.setup(gammadatafile, Lymanlimitdatafile, T_b_Datafile)

filename = "/Users/users/achatterjee/Data/GW_with_CosmoReionMC/mcmc/MCMC_sampler/test_run_CMB_Reion_GW.h5"
backend = emcee.backends.HDFBackend(filename)

sampler = MCMCsampler(
                params = params, 
                likelihoodComputation = like, 
                filePrefix = '/Users/users/achatterjee/Data/GW_with_CosmoReionMC/mcmc/MCMC_sampler/%s'%str(d1), 
                chain_storage_path = chain_storage_path, 
                fileroot = file_root,
                walkersRatio = 4,  
                sampleIterations = 1000000,
                n_cores = 1,
                backend=backend
                )


print("start sampling")
sampler.startSampling()
print("done!")


