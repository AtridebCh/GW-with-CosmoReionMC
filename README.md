############# Installation Guide  ####################
You need to install the following package in order to run the code succesfully.

The required packages are:


1. CAMB-1.1.0_modified
2. Planck "clik" code
3. Bilby (GW)



1.Create a new conda environment for hassle free installation

conda create -n camb_gw_cosmo python=3.10 -y
conda activate camb_gw_cosmo

conda install -c conda-forge numpy scipy cython matplotlib astropy gfortran emcee corner -y

2. Install GW packages

pip install bilby (for GW part)
pip install lalsuite. (for GW part)




########################
Go inside the directory "GW_with_CosmoReionMC". From now on, we will call this directory as our parent directory.

Here We will describe how to install and run different components of the package.
###############################

##################################
3. CAMB-1.1.0_modified:

(i) Inside parent directory you will see a directory called CAMB-1.1.0_modified.

(ii) Open the terminal inside the CAMB-1.1.0_modified directory

(iii)Please note that in order to run the camb succesfully you have to have gfortran version> 6.3 or higher (check with gfortran --version). This is usually the case. In case, it is not please go to 
   "https://camb.readthedocs.io/en/latest/fortran_compilers.html#fortran-compilers"  to sort out the gfortran issue. If you sort out the gfortran issue succesfully you will see gfortran version >6.3 when you check the gfortran version. I will strongly recommend you to make this gfortran (i.e. gfortran >6.3) as your default gfortran version.

(iv)Give the following command (inside conda environment)
         python setup.py build
         python setup.py install --user  #This --user is very important please don't miss this. If you use --user here then use --user in all python setup.py install --user
         
(v) If you edit or change anything in camb, you need to do
rm -rf build/
python setup.py build
python setup.py install --user
         

(vi) Inside this directory, there is a code called "test_modified_camb.py". As you will see, we put an external reionization history from outside (I kept a precomputed reionization history in QHII.dat) camb in this code. Run this, if it runs successfully it will generate a plot which shows the redshift evolution of Q_{HII} and x_e (i.e., including the effect of Helium reionization). It will also produce plots showing cl_TT, cl_EE, cl_BB, cl_TE.
#######################################################

##################################################
4.Planck "clik" code

# Undobtedly, the most tricky installation in this package.
# One way could be to install cobaya, which automatically install planck likelihood, but I have never installed it that way.

(i) Download the data and likelihood from here
https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/CMB_spectrum_%26_Likelihood_Code#2018_Likelihood

(ii) cd /path/to/planck_code/plc_3.0/plc-3.1 (wherever you kept the planck code)


# activate your environment first
conda activate camb_gw_cosmo


#if you want to install with waf (recommended) then (a detailed but a bit old guideline is given in https://cosmologist.info/cosmomc/readme_planck.html)

#(i) CFITSIO — required by clik
sudo apt install libcfitsio-dev

# (ii) Optional but recommended
sudo apt install libblas-dev liblapack-dev
 
#Configure — waf is the build system bundled with plc
(iv) python waf configure --lapack_mkl=/path/to/mkl   # if using MKL
# or
(v)python waf configure --lapack_blas               # if using system BLAS

(vi) python waf install

##If configuration fails, the most common fix is:
python waf configure --cfitsio_prefix=/usr --lapack_blas

#Step 4 — Set environment variables
#Add to your .bashrc or .bash_profile
export CLIK_PATH=/path/to/plc_3.0
export CLIK_DATA=/path/to/planck_data   # where you extracted the data tarball
export LD_LIBRARY_PATH=$CLIK_PATH/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$CLIK_PATH/lib/pkgconfig:$PKG_CONFIG_PATH

# Then
source ~/.bashrc

#Step 5 — Install Python wrapper

source clik_profile.sh   # sets up environment

pip install ./
# or
python setup.py install

#Step 6 — Test
python
import clik

# Load one of the likelihood files from the data tarball
clik_file = "/path/to/planck_data/plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik"
like = clik.clik(clik_file)

print(like.get_lmax())        # should print [2508, 2508, 2508, -1, -1, -1]
print(like.get_extra_parameter_names())   # nuisance parameter names

#If the waf fails then using make
 
#### Makefile installation ###########
conda install -c conda-forge cfitsio
conda install -c conda-forge lapack blas
conda install cython

make clean

make install \
     CFITSIOPATH=$CONDA_PREFIX \
     FC=gfortran \
     LAPACK="-L$CONDA_PREFIX/lib -llapack -lblas" \
     LAPACKLIBPATH=$CONDA_PREFIX/lib \
     FRUNTIME="-lgfortran -lgomp -lm" \
     "SIMALLLKL=\$(addprefix \$(ODIR)/,clik_simall.o)" \
     CFLAGS="-m64 -fopenmp -fPIC -Isrc -Isrc/minipmc -Isrc/cldf -Isrc/simall -Isrc/plik -Isrc/lenslike/plenslike -D HAS_LAPACK -D LAPACK_CLIK -D NOHEALPIX -D CLIK_LENSING -D CAMSPEC_V1 -I$CONDA_PREFIX/include"

make install_python \
     CFITSIOPATH=$CONDA_PREFIX \
     FC=gfortran \
     LAPACK="-L$CONDA_PREFIX/lib -llapack -lblas" \
     LAPACKLIBPATH=$CONDA_PREFIX/lib \
     FRUNTIME="-lgfortran -lgomp -lm" \
     "SIMALLLKL=\$(addprefix \$(ODIR)/,clik_simall.o)"


#if it runs successfully then try
python -c "import clik; print('clik imported successfully')"
python -c " import clik 
l = clik.clik('/home/atri/GW_with_CosmoReionMC/COM_Likelihood_Data-baseline_R3.00/baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik')
print('simall loaded, lmax:', l.lmax)
print('nuisance:', l.get_extra_parameter_names())
"
#if not then
cat svnversion 
#this should show you plc_3.1, if not then
sed -i 's/version = open("svnversion").read()+" MAKEFILE"/version = "3.1"/' setup.py

#then
make install_python      CFITSIOPATH=$CONDA_PREFIX \
     FC=gfortran \
     LAPACK="-L$CONDA_PREFIX/lib -llapack -lblas"      
     LAPACKLIBPATH=$CONDA_PREFIX/lib      
     FRUNTIME="-lgfortran -lm"

If this works then in the terminal (othetwise clik won't be imported next time you activate conda environment)

cat > $CONDA_PREFIX/etc/conda/activate.d/clik_env.sh << EOF
export CLIK_PATH=/path_to_COM_Likelihood_Code-v3.0_R3.10/code/plc_3.0/plc-3.1
export CLIK_PLUGIN=rel2015
export LINK_CLIK="-L\$CLIK_PATH/lib -lclik"
export LD_LIBRARY_PATH=\$CLIK_PATH/lib:\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH
export PYTHONPATH=\$CLIK_PATH/lib/python3.1/site-packages:\$PYTHONPATH
EOF

####### End of Makefile installation ###########################################

5. Compiling Reionization with F2py

(i) Inside the parent directory, there is a folder called ReionizationWithF2py
(ii)  cd to that directory, you will see a folder called src.
(iii) I kept all the fortran code in src, open the SEDreader_mod.f90, you have to change the path to the SED files over there (line no. 27 and 28)
(v) make distclean (pwd should show you are inside ReionizationWithF2py, but not inside src)
(vi) make
(vii) make python
If you make any change in fortran files, then repeat (steps v-vii)
(viii) once they run successfully, run python ./python_dir/test_model.py, it will produce a few plots. Check if those plots make sense.

##############################
### After this step, one has to generate mock data for GW and 21cm. In this case, I have already created (for fiducial cosmo and astrophysical parameters) and put them in respective folders, i.e, GW catalogue in GravWave and Global signal in global_21cm_signal ###
(i) If you want to change/create a new GW catalogue, go to GravWave folder inside parent directory. Change path in createMockSignal.py and run it
(ii) If you want to change/create a new 21cm signal, go to global_21cm_signal inside parent directory. Change path in crearMockCatalogu.py and run it
###################################################       
Running the MCMC code:

Inside parent directory there is a folder called "mcmc". Entering the folder you will get another directory called MCMC_sampler. In that folder,

(i) open the Liklihood.py, then you have to change the path of few data files and likleihood, as mentioned inside the code. find the word "change" with ctrl+F, then change the path.
(ii) open the CoredataGenerator.py, change the paths as mentioned in the file. Similar to previous step
(iii) Once you are done, open the python file called single_run.py
(iv) run that (after changing the path of the data file), it should produce a few plots for the fiducial cosmology, showing different observables, take a look at those plots, check if those plots make sense.
(v) run the mcmc code, python mcmc_run.py (again you have to change the path of the data files)
(vi) to run in multiple cores, change the variable n_cores inside the code.

################################################


