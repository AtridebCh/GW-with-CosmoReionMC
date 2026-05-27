import numpy as np

from utils import Params
from Likelihood import Likelihood
from CoredataGenerator import CoreModule

 
Lymanlimitdatafile = '../../ReionizationWithF2Py/ObsData/Lyman_limit.dat'
gammadatafile      = '../../ReionizationWithF2Py/ObsData/gamma_data_all_combined.dat'
T_b_Datafile       = '../../global_21cm_signal/brightness_temp.txt'

params = Params(("hubble",                [6.821184e+01, 65, 80, 0.01*(80-65.0)]),
                ("ombh2",                 [2.261010e-02, 0.01, 0.03, 0.01*(0.03-0.01)]),
                ("omch2",                 [1.171812e-01, 0.09, 0.2, 0.01*(0.2-0.09)]),
                ("As",                    [2.121276e-09, 1.8e-9, 2.35e-9, 0.01*(2.35e-09-1.8e-09)]),
                ("ns",                    [9.605819e-01, 0.8, 1.2, 0.01*(1.2-0.8)]),
                ("omega_zero",            [-1.0, -2.0, 2.0, 0.1]),
                ("omega_a",               [0.0, -1.0, 2.0, 0.1]),
                ("f_zero",                [0.006, 0.0, 1.0, 0.01]),
                ("alpha_lo",              [0.0, 0.0, 1.0, 0.1]),
                ("alpha_hi",              [0.0, 0.0, 5.0, 0.1]),
                ("alpha_z",               [1e-05, 0.0, 1.0, 0.01]),
                ("esc_II",                [0.25, 0.0, 1.0, 0.1]),
                ("lambda0",               [2.64590e+00, 1.0,10.0,0.01*(6.0-1.0)]),
                ('f_X',                   [1.977400e-01, 0.01, 1.0, 0.001]),
                ('f_alpha',               [1.000560e+00, 0.01, 1.0, 0.001]))


print('params_value', params[:,0])


like = Likelihood(CoreModule)
like.setup(gammadatafile, Lymanlimitdatafile, T_b_Datafile)


likelihood, blobs=like(params[:,0], single_run=True)
print('likelihhod=%s, tau=%s, QHII_5point8=%s' %(likelihood, blobs[0], blobs[1]))
