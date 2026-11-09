import numpy as np

from utils import Params
from Likelihood import Likelihood
from CoredataGenerator import CoreModule

 
Lymanlimitdatafile = '/Users/users/achatterjee/Data/GW_with_CosmoReionMC/ReionizationWithF2Py/ObsData/Lyman_limit.dat'
gammadatafile      = '/Users/users/achatterjee/Data/GW_with_CosmoReionMC/ReionizationWithF2Py/ObsData/gamma_data_all_combined.dat'
T_b_Datafile       = '/Users/users/achatterjee/Data/GW_with_CosmoReionMC/global_21cm_signal/brightness_temp.txt'

'''

parameter bestfit
hubble 6.6389000E+01
ombh2 2.2172000E-02
omch2 1.2189000E-01
As 2.3213000E-09
ns 9.7157000E-01
f_zero 1.0773000E-02
alpha_lo 6.7158000E-01
alpha_hi -6.4627000E-01
alpha_z 9.8518000E-01
esc_II 1.2568000E-03
'''
'''
params =  Params(("hubble", [6.5685000E+01,  6.5000000E+01,  6.7481000E+01, 1.0]),
                 ("ombh2",  [2.2100000E-02,  2.1045000E-02,  2.2905000E-02, 0.001]),  
                 ("omch2",     [1.2349000E-01,  1.2025000E-01,  1.2612000E-01,  0.001]),
                 ("As",        [2.0704000E-09,  2.0488000E-09,  2.1365000E-09,  1.0488000E-09]),
                 ("ns",        [9.6368000E-01,  9.3983000E-01,  9.8098000E-01,  0.01]),
                 ("f_zero",    [7.4385000E-01,  2.1280000E-01,  8.5916000E-01,  0.01]),
                 ("alpha_lo",  [7.0989000E-01,  1.6211000E-01,  9.9986000E-01,  0.01]),  
                 ("alpha_hi",  [6.9776000E-01,  2.6815000E-01,  7.6084000E-01,  0.01]),  
                 ("alpha_z",   [1.9512000E-04, -1.9437000E-01,  3.8487000E-01, 0.01]),  
                 ("esc_II",    [3.9177000E-02,  2.9526000E-02,  5.9033000E-02,  0.01]),  
                 ("lambda0",   [4.0635000E+00,  3.7502000E+00,  5.5415000E+00,  0.1]),  
                 ("f_X",       [1.8414000E-01,  1.6663000E-01,  2.2074000E-01,  0.01]),  
                 ("f_alpha",   [9.9623000E-01,  9.8317000E-01,  9.9999000E-01,  0.1]))


parameter   bestfit        lower1         upper1         lower2         upper2
hubble      6.7838000E+01  6.5991000E+01  7.0212000E+01  6.5991000E+01  7.0484000E+01   H_{0}
ombh2       2.2374000E-02  2.1670000E-02  2.2928000E-02  2.1670000E-02  2.3099000E-02   \Omega_{b} h^{2}
omch2       1.1879000E-01  1.1354000E-01  1.2228000E-01  1.1186000E-01  1.2288000E-01   \Omega_{c} h^{2}
As          2.2884000E-09  2.2074000E-09  2.3479000E-09  2.2066000E-09  2.3493000E-09   A_s
ns          9.6895000E-01  9.5228000E-01  9.8662000E-01  9.5228000E-01  9.8980000E-01   n_s
f_zero      4.2830000E-02  2.9707000E-02  1.1905000E-01  2.7073000E-02  1.9523000E-01   f_{0}
alpha_lo    5.0793000E-01  2.1389000E-01  8.4751000E-01  1.7353000E-01  8.4751000E-01   \alpha_{lo}
alpha_hi   -4.8872000E-01 -9.7282000E-01 -2.6587000E-01 -9.7282000E-01 -1.2709000E-01   \alpha_{hi}
alpha_z     3.2596000E-01 -9.9961000E-01  9.8228000E-01 -9.9961000E-01  9.8871000E-01   \alpha_{z}
esc_II      6.1687000E-03  2.5241000E-03  1.6143000E-02  2.5241000E-03  2.7548000E-02   f_{esc, II}
lambda0     9.7823000E+00  7.9077000E+00  1.0000000E+01  5.7017000E+00  1.0000000E+01   \lambda_0
f_X         5.6418000E-01  1.0865000E-02  9.9989000E-01  1.0865000E-02  9.9989000E-01   f_{X}
f_alpha     9.8455000E-01  3.7886000E-01  9.9997000E-01  3.6143000E-01  9.9997000E-01   f_{\alpha}
tau         1.0066000E-01  8.2009000E-02  1.1482000E-01  8.2009000E-02  1.2109000E-01   \tau


'''
params = Params(("hubble", [6.7380389E+01, 65, 80, 0.01*(80-65.0)]), #6.726630e+01
                ("ombh2", [2.2358184E-02, 0.01, 0.03, 0.01*(0.03-0.01)]), # 2.2341e-02
                ("omch2", [1.1992849E-01, 0.09, 0.2, 0.01*(0.2-0.09)]), #1.206109e-01
                ("As", [2.1208786E-09, 1.8e-9, 2.35e-9, 0.01*(2.35e-09-1.8e-09)]), #2.112747e-09
                ("ns", [9.6480108E-01, 0.8, 1.2, 0.01*(1.2-0.8)]), #9.655997e-01
                ("log_f_zero", [-4.3567402E-01, -4.0, 0.0, 0.1]),#0.19 GW depends 0.63408256
                ("alpha_lo", [0.0, 0.0, 1.0, 0.1]),#0.31 GW depends
                ("alpha_hi", [0.0, 0.0, 1.0, 0.1]),#0.11 GW depdends
                ("alpha_z", [0.0, -1.0, 1.0, 0.01]),#0.0
                ("log_esc_II", [-1.7893972E+00, -4.0, 0.0, 0.1]),#0.035; -2.077428 GW depends
                ("lambda0", [5.5645855E+00, 1.0, 12.0, 0.01*(6.0-1.0)]), #4.98 GW does not depend
                ("log_f_X", [np.log10(0.2), -4.0, 0.0, 0.1]),
                ("log_f_alpha", [0.0, -1.0, 1.0, 0.01]))
#5       490.02998       66.954859     0.022302521      0.12140578   2.1155203e-09      0.96220522      0.63408256       -2.077428       7.0204522      0.91200322     0.057762674       -21.77744        -21.77744       1023.6148                              1023.6148                
#('f_X', [2.00e-01, 0.01, 1.0, 0.001]),#2.00e-1
#('f_alpha', [1.0e+00, 0.01, 1.0, 0.001]))#1.0

print('params_value', params[:,0])


like = Likelihood(CoreModule, min_param = params[:,1], max_param= params[:,2])
like.setup(gammadatafile, Lymanlimitdatafile, T_b_Datafile)


likelihood, blobs=like(params[:,0], single_run=True)
print('likelihood=%s, tau=%s, QHII_5point8=%s' %(likelihood, blobs[0], blobs[1]))
