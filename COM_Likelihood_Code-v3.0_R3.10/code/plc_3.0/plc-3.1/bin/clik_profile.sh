# this code cannot be run directly
# do 'source /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/bin/clik_profile.sh' from your sh shell or put it in your profile

function addvar () {
local tmp="${!1}" ;
tmp="${tmp//:${2}:/:}" ; tmp="${tmp/#${2}:/}" ; tmp="${tmp/%:${2}/}" ;
export $1="${2}:${tmp}" ;
} 

if [ -z "${PATH}" ]; then 
PATH=/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/bin
export PATH
else
addvar PATH /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/bin
fi
if [ -z "${PYTHONPATH}" ]; then 
PYTHONPATH=/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib/python3.1/site-packages
export PYTHONPATH
else
addvar PYTHONPATH /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib/python3.1/site-packages
fi
if [ -z "${LD_LIBRARY_PATH}" ]; then 
LD_LIBRARY_PATH=/usr/lib
export LD_LIBRARY_PATH
else
addvar LD_LIBRARY_PATH /usr/lib
fi
if [ -z "${LD_LIBRARY_PATH}" ]; then 
LD_LIBRARY_PATH=/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib
export LD_LIBRARY_PATH
else
addvar LD_LIBRARY_PATH /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib
fi
if [ -z "${LD_LIBRARY_PATH}" ]; then 
LD_LIBRARY_PATH=/home/atri/anaconda3/envs/camb_gw_cosmo/lib
export LD_LIBRARY_PATH
else
addvar LD_LIBRARY_PATH /home/atri/anaconda3/envs/camb_gw_cosmo/lib
fi
if [ -z "${LD_LIBRARY_PATH}" ]; then 
LD_LIBRARY_PATH=/home/atri/anaconda3/envs/camb_gw_cosmo/lib
export LD_LIBRARY_PATH
else
addvar LD_LIBRARY_PATH /home/atri/anaconda3/envs/camb_gw_cosmo/lib
fi

CLIK_PATH=/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1
export CLIK_PATH

CLIK_DATA=/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/share/clik
export CLIK_DATA

CLIK_PLUGIN=rel2015
export CLIK_PLUGIN

