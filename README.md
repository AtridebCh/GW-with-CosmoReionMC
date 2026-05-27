# GW_with_CosmoReionMC — Installation Guide

## Prerequisites

Before starting, create and activate a dedicated conda environment:

```bash
conda create -n camb_gw_cosmo python=3.10 -y (python 3.11 gives completely wrong answer, python3.9 is also safe)
conda activate camb_gw_cosmo

conda install -c conda-forge numpy scipy cython matplotlib astropy gfortran emcee corner -y
```

Install gravitational wave packages:

```bash
pip install bilby
pip install lalsuite
```

> **Note:** All steps below assume you are working from inside the `GW_with_CosmoReionMC/` directory (referred to as the **parent directory**).

---

## 1. CAMB-1.1.0_modified

**Location:** `<parent>/CAMB-1.1.0_modified/`

### Requirements

- `gfortran >= 6.3` — check with `gfortran --version`
- If your version is too old, follow the [CAMB Fortran compiler guide](https://camb.readthedocs.io/en/latest/fortran_compilers.html#fortran-compilers) and set the upgraded version as your default.

### Build & Install

```bash
cd CAMB-1.1.0_modified/
python setup.py build
python setup.py install --user   # The --user flag is required
```

> **Tip:** If you modify any CAMB source files, rebuild with:
> ```bash
> rm -rf build/
> python setup.py build
> python setup.py install --user
> ```

### Verify Installation

Run the test script:

```bash
python test_modified_camb.py
```

If successful, this will:
- Read an external reionization history from `QHII.dat`
- Generate plots of Q_HII and x_e (including Helium reionization)
- Generate C_ℓ plots: TT, EE, BB, TE

---

### Installing in kapteyn cluster

1. module load anaconda3/2022.10; otherwise the default anaconda environment has python 3.11 which will give incorrect result.

2. use full path for every file i.e., change file path in CoredataGenerator.py and Likleihood.py, MCMC_run.py etc

3. cd /path/to/CAMB-1.1.0_modified
   First delete all the files that was generated in earlier installation
   find . -name "*.o" -delete
   find . -name "*.so" -delete
   find . -name "*.mod" -delete
   rm -rf build/

   Now rebuild with HARD override:

   export FFLAGS="-O0 -march=x86-64"
   export FCFLAGS="-O0 -march=x86-64"
   export CFLAGS="-O0 -march=x86-64"

   Then:

   python setup.py build 
   python setup.py install --user

Although installing clipy-like was easy, while running (in python3.9), it gives an error with .astype, so had to open /Users/users/achatterjee/.local/lib/python3.9/site-packages/clipy/simall.py and change

dl = jnp.astype(jnp.floor(cls[1,self.lmin:]*self.llp1/self.stepEE),jnp.int32)

to:

dl = jnp.floor(
    cls[1,self.lmin:] * self.llp1 / self.stepEE
).astype(jnp.int32)

so, change all jnp.astype in that code in the above mentioned way


## 2. Planck "clik" Likelihood Code

We would strongly recommend to install using clipy from https://github.com/benabed/clipy

If this installation is succesfull, then skip the "clik" installation completely. If using clipy two changes are required

(i) replace "import clik" with "import clipy as clik" in path/to/GW_with_CosmoReionMC/mcmc/MCMC_sampler/Likelihood.py
(ii) replace loglike += l(input)[0] with loglike += l(input) in path/to/GW_with_CosmoReionMC/mcmc/MCMC_sampler/Likelihood.py


> **Warning:** This is the most complex installation in the package.

### Step 1 — Download Data & Likelihood

Download from the [Planck Legacy Archive](https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/CMB_spectrum_%26_Likelihood_Code#2018_Likelihood).

### Step 2 — Navigate to the plc Directory

```bash
cd /path/to/planck_code/plc_3.0/plc-3.1
conda activate camb_gw_cosmo
```

### Step 3 — Install System Dependencies

```bash
sudo apt install libcfitsio-dev          # Required
sudo apt install libblas-dev liblapack-dev  # Optional but recommended
```

### Step 4 — Configure & Build (waf, recommended)

```bash
# With MKL:
python waf configure --lapack_mkl=/path/to/mkl

# With system BLAS:
python waf configure --lapack_blas

# If configuration fails, try:
python waf configure --cfitsio_prefix=/usr --lapack_blas

python waf install
```

For a detailed (if slightly dated) guide, see [cosmologist.info/cosmomc/readme_planck.html](https://cosmologist.info/cosmomc/readme_planck.html).

### Step 5 — Set Environment Variables

Add the following to your `~/.bashrc` or `~/.bash_profile`:

```bash
export CLIK_PATH=/path/to/plc_3.0
export CLIK_DATA=/path/to/planck_data
export LD_LIBRARY_PATH=$CLIK_PATH/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$CLIK_PATH/lib/pkgconfig:$PKG_CONFIG_PATH
```

Then reload:

```bash
source ~/.bashrc
```

### Step 6 — Install Python Wrapper

```bash
source clik_profile.sh
pip install ./
# or
python setup.py install
```

### Step 7 — Verify

```python
import clik

clik_file = "/path/to/planck_data/plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik"
like = clik.clik(clik_file)

print(like.get_lmax())                    # Expected: [2508, 2508, 2508, -1, -1, -1]
print(like.get_extra_parameter_names())   # Nuisance parameter names
```

---

### Alternative: Makefile Installation

If `waf` fails, use `make` instead:

```bash
conda install -c conda-forge cfitsio lapack blas
conda install cython

make clean

make install \
  CFITSIOPATH=$CONDA_PREFIX \
  FC=gfortran \
  LAPACK="-L$CONDA_PREFIX/lib -llapack -lblas" \
  LAPACKLIBPATH=$CONDA_PREFIX/lib \
  FRUNTIME="-lgfortran -lgomp -lm" \
  "SIMALLLKL=\$(addprefix \$(ODIR)/,clik_simall.o)" \
  CFLAGS="-m64 -fopenmp -fPIC -Isrc -Isrc/minipmc -Isrc/cldf -Isrc/simall \
          -Isrc/plik -Isrc/lenslike/plenslike \
          -D HAS_LAPACK -D LAPACK_CLIK -D NOHEALPIX -D CLIK_LENSING -D CAMSPEC_V1 \
          -I$CONDA_PREFIX/include"

make install_python \
  CFITSIOPATH=$CONDA_PREFIX \
  FC=gfortran \
  LAPACK="-L$CONDA_PREFIX/lib -llapack -lblas" \
  LAPACKLIBPATH=$CONDA_PREFIX/lib \
  FRUNTIME="-lgfortran -lgomp -lm" \
  "SIMALLLKL=\$(addprefix \$(ODIR)/,clik_simall.o)"
```

Verify:

```bash
python -c "import clik; print('clik imported successfully')"
python -c "
import clik
l = clik.clik('/path/to/COM_Likelihood_Data/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik')
print('simall loaded, lmax:', l.lmax)
print('nuisance:', l.get_extra_parameter_names())
"
```

If the import fails, check the version file:

```bash
cat svnversion   # Should show: plc_3.1
```

If it doesn't, patch it:

```bash
sed -i 's/version = open("svnversion").read()+" MAKEFILE"/version = "3.1"/' setup.py
```

Then reinstall:

```bash
make install_python \
  CFITSIOPATH=$CONDA_PREFIX \
  FC=gfortran \
  LAPACK="-L$CONDA_PREFIX/lib -llapack -lblas" \
  LAPACKLIBPATH=$CONDA_PREFIX/lib \
  FRUNTIME="-lgfortran -lm"
```

### Persist Environment Variables Across conda Activations

So that `clik` is importable every time you activate the environment:

```bash
cat > $CONDA_PREFIX/etc/conda/activate.d/clik_env.sh << EOF
export CLIK_PATH=/path/to/plc_3.0/plc-3.1
export CLIK_PLUGIN=rel2015
export LINK_CLIK="-L\$CLIK_PATH/lib -lclik"
export LD_LIBRARY_PATH=\$CLIK_PATH/lib:\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH
export PYTHONPATH=\$CLIK_PATH/lib/python3.1/site-packages:\$PYTHONPATH
EOF
```

---

## 3. Reionization Module (F2py)

**Location:** `<parent>/ReionizationWithF2py/`

### Setup

Before compiling, open `src/SEDreader_mod.f90` and update the SED file paths at **lines 27–28**.

### Build

```bash
cd ReionizationWithF2py/   # Must be here, NOT inside src/
make distclean
make
make python                # May need setuptools version < 70
```

> If you modify any Fortran source files, repeat the three `make` steps above.

### Verify

```bash
python ./python_dir/test_model.py
```

This produces several diagnostic plots — inspect them to confirm the output looks physically reasonable.

---

## 4. Generating Mock Data

Pre-generated mock data for fiducial cosmological and astrophysical parameters is already included:
- **GW catalogue** — `<parent>/GravWave/`
- **21cm global signal** — `<parent>/global_21cm_signal/`

To regenerate or modify run the following code from ReionizationF2py directory:

| Dataset | Directory | Script |
|---|---|---|
| GW catalogue | `GravWave/` | `createMockCatalogue.py`|
| 21cm signal | `global_21cm_signal/` | `createMockSignal.py`  |

Update the file paths inside each script before running.

---

## 5. Running the MCMC

**Location:** `<parent>/mcmc/MCMC_sampler/`

### Setup

1. Open `Likelihood.py` — search for the keyword `change` (`Ctrl+F`) and update all data file and likelihood paths.
2. Open `CoredataGenerator.py` — do the same.

### Test Run

```bash
python single_run.py
```

This runs the pipeline for the fiducial cosmology and generates diagnostic plots for all observables. Verify the plots look sensible before proceeding.

### MCMC Run

```bash
python mcmc_run.py
```

## 6. Removing or adding new parameters in MCMC
**Location:** `<parent>/mcmc/MCMC_sampler/`
1. open CoredataGenerator.py; add/remove the name of the parameters in FREE_PARAM_MAPPING_DEFAULT 
2. If the parameter is required in CAMB/Reionizatiom mode then make the change in the argument while calling CAMB/run_model inside _call__ method 
3. In single_run.py add/remove the name of the parameters along with its initial values, range, step size in params = Params(( ...))
4. Remember to run the code for generating the mock catalogues. To do this one needs to open the code that create the catalogues and then change the parameter values inside

## 7. Updates on 27/05/2026
1. Added dynamical dark energy with "ppf" to even consider phantom models, the fluid model often gives error that "w crossing -1".
2. Added isValid function in /path/to/GW_with_CosmoreionMC/mcmc/MCMC_sampler/Likelihood.py which ensures that the free parameters are inside the prior range.
3. "new" method for calculating df_coll/dt has been implemented correctly, the calculations are also updated in the SFRD notes in overleaf.
4. Removed data points at z<3.5 in gamma_data and Lyman Limit data to ensure numerical convergence
5. Added condition to ensure the reionization history to be monotonic in CoredataGenerator.py 
6. Updated the convergence test and postprocessing codes in the postprocess directory.
7. Replaced clik with clipy (https://github.com/benabed/clipy) 
