# this code cannot be run directly
# do 'source /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/bin/clik_profile.csh' from your csh shell or put it in your profile

 

if !($?PATH) then
setenv PATH /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/bin
else
set newvar=$PATH
set newvar=`echo ${newvar} | sed s@:/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/bin:@:@g`
set newvar=`echo ${newvar} | sed s@:/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/bin\$@@` 
set newvar=`echo ${newvar} | sed s@^/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/bin:@@`  
set newvar=/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/bin:${newvar}                     
setenv PATH /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/bin:${newvar} 
endif
if !($?PYTHONPATH) then
setenv PYTHONPATH /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib/python3.1/site-packages
else
set newvar=$PYTHONPATH
set newvar=`echo ${newvar} | sed s@:/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib/python3.1/site-packages:@:@g`
set newvar=`echo ${newvar} | sed s@:/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib/python3.1/site-packages\$@@` 
set newvar=`echo ${newvar} | sed s@^/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib/python3.1/site-packages:@@`  
set newvar=/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib/python2.7/site-packages:${newvar}                     
setenv PYTHONPATH /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib/python2.7/site-packages:${newvar} 
endif
if !($?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH /usr/lib
else
set newvar=$LD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@:/usr/lib:@:@g`
set newvar=`echo ${newvar} | sed s@:/usr/lib\$@@` 
set newvar=`echo ${newvar} | sed s@^/usr/lib:@@`  
set newvar=/usr/lib:${newvar}                     
setenv LD_LIBRARY_PATH /usr/lib:${newvar} 
endif
if !($?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib
else
set newvar=$LD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@:/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib:@:@g`
set newvar=`echo ${newvar} | sed s@:/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib\$@@` 
set newvar=`echo ${newvar} | sed s@^/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib:@@`  
set newvar=/home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib:${newvar}                     
setenv LD_LIBRARY_PATH /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/lib:${newvar} 
endif
if !($?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH /home/atri/anaconda3/envs/camb_gw_cosmo/lib
else
set newvar=$LD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@:/home/atri/anaconda3/envs/camb_gw_cosmo/lib:@:@g`
set newvar=`echo ${newvar} | sed s@:/home/atri/anaconda3/envs/camb_gw_cosmo/lib\$@@` 
set newvar=`echo ${newvar} | sed s@^/home/atri/anaconda3/envs/camb_gw_cosmo/lib:@@`  
set newvar=/home/atri/anaconda3/envs/camb_gw_cosmo/lib:${newvar}                     
setenv LD_LIBRARY_PATH /home/atri/anaconda3/envs/camb_gw_cosmo/lib:${newvar} 
endif
if !($?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH /home/atri/anaconda3/envs/camb_gw_cosmo/lib
else
set newvar=$LD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@:/home/atri/anaconda3/envs/camb_gw_cosmo/lib:@:@g`
set newvar=`echo ${newvar} | sed s@:/home/atri/anaconda3/envs/camb_gw_cosmo/lib\$@@` 
set newvar=`echo ${newvar} | sed s@^/home/atri/anaconda3/envs/camb_gw_cosmo/lib:@@`  
set newvar=/home/atri/anaconda3/envs/camb_gw_cosmo/lib:${newvar}                     
setenv LD_LIBRARY_PATH /home/atri/anaconda3/envs/camb_gw_cosmo/lib:${newvar} 
endif

setenv CLIK_PATH /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1
setenv CLIK_DATA /home/atri/Downloads/COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1/share/clik

setenv CLIK_PLUGIN rel2015

