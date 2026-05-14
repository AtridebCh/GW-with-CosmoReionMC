import sys, os

import numpy as np

from utils import Params
from Likelihood import Likelihood
from CoredataGenerator import CoreModule
from sampler import MCMCsampler


Lymanlimitdatafile = '../../ReionizationWithF2Py/ObsData/Lyman_limit.dat'
gammadatafile      = '../../ReionizationWithF2Py/ObsData/gamma_data_all_combined.dat'
T_b_Datafile       = '../../global_21cm_signal/brightness_temp.txt'


params = Params(("hubble",                [6.821184e+01, 65, 80, 0.01*(80-65.0)]),
                ("ombh2",                 [2.261010e-02, 0.01, 0.03, 0.01*(0.03-0.01)]),
                ("omch2",                 [1.171812e-01, 0.09, 0.2, 0.01*(0.2-0.09)]),
                ("As",                    [2.121276e-09, 1.8e-9, 2.35e-9, 0.01*(2.35e-09-1.8e-09)]),
                ("ns",                    [9.605819e-01, 0.8, 1.2, 0.01*(1.2-0.8)]),
                ("omega_zero",            [-1.0, -2.0, 2.0, 0.1]),
                ("omega_a",               [0.0, -2.0, 2.0, 0.1]),
                ("f_zero",                [0.17, 0.0, 1.0, 0.01]),
                ("alpha_lo",              [0.61, 0.0, 1.0, 0.1]),
                ("alpha_hi",              [1.1, 0.0, 5.0, 0.1]),
                ("alpha_z",               [1e-05, 0.0, 1.0, 1e-06]),
                ("esc_II",                [0.36, 0.0, 1.0, 0.1]),
                ("lambda0",               [12.364590e+00, 1.0,10.0,0.01*(6.0-1.0)]),
                ('f_X',                   [1.977400e-01, 0.01, 1.0, 0.001]),
                ('f_alpha',               [1.000560e+00, 0.01, 1.0, 0.001]))

from datetime import date

today = date.today()

# dd/mm/YY
d1 = today.strftime("%d_%m_%Y")

param_latex_name_arr=['H_{0}','\omega_{b} h^{2}', '\omega_{c} h^{2}', 'A_s', 'n_s', 'f_{0}', 'alpha_{lo}', 'alpha_{hi}', 'alpha_{z}', 'f_{esc, II}', '\lambda_0', 'f_{X}', 'f_{alpha}']



#change the path if needed
parent_dir="./chain_storage/" 
MCMC_dir="test_run"+d1

chain_storage_path = os.path.join(parent_dir, MCMC_dir)

os.makedirs(chain_storage_path, exist_ok=True)

file_root='test_MCMC'

write_file=Params.write_paramnames_ranges_file(params[:,1], params[:,2], params.keys, param_latex_name_arr, chain_storage_path, file_root)



like = Likelihood(CoreModule)
like.setup(gammadatafile, Lymanlimitdatafile, T_b_Datafile)

sampler = MCMCsampler(
                params= params, 
                likelihoodComputation=like, 
                filePrefix=str(d1), 
                chain_storage_path=chain_storage_path, 
                fileroot=file_root,
                walkersRatio=4,  
                sampleIterations=1000000,
                n_cores = 1
                )


print("start sampling")
sampler.startSampling()
print("done!")


