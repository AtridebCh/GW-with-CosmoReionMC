# ===========================================================================
# Compiler and flags
# ===========================================================================
FC      = gfortran
FFLAGS = -O2 -Wall -Wextra -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fPIC -ffree-line-length-none
FFLAGS_PY = -O2 -fPIC -ffree-line-length-none -J$(MODDIR) -I$(MODDIR)
# ===========================================================================
# Directories
# ===========================================================================
SRCDIR  = src
OBJDIR  = obj
MODDIR  = mod
PYDIR   = python_dir
# ===========================================================================
# Source files — ORDER MATTERS for module dependencies
# ===========================================================================
SRCS =  $(SRCDIR)/kinds_mod.f90                 \
        $(SRCDIR)/constants_mod.f90             \
        $(SRCDIR)/parameters_mod.f90            \
        $(SRCDIR)/state_mod.f90                 \
        $(SRCDIR)/variables_mod.f90             \
        $(SRCDIR)/adaptint_mod.f90              \
        $(SRCDIR)/findroot_mod.f90              \
        $(SRCDIR)/onedspline_mod.f90            \
        $(SRCDIR)/backgroundCosmology_mod.f90   \
        $(SRCDIR)/crossSection_mod.f90          \
        $(SRCDIR)/recrates_mod.f90              \
        $(SRCDIR)/SEDreader_mod.f90             \
        $(SRCDIR)/thermal_mod.f90               \
        $(SRCDIR)/stellar_mod.f90               \
        $(SRCDIR)/inhomoReion_mod.f90           \
        $(SRCDIR)/utils_mod.f90                 \
        $(SRCDIR)/reionization_mod.f90          \
        $(SRCDIR)/driver.f90

OBJS = $(patsubst $(SRCDIR)/%.f90, $(OBJDIR)/%.o, $(SRCS))
# all objects except driver.o — f2py provides its own entry point
OBJS_PY = $(filter-out $(OBJDIR)/driver.o, $(OBJS))
# ===========================================================================
# Target
# ===========================================================================
EXE = reion

.PHONY: all clean distclean

all: $(OBJDIR) $(MODDIR) $(EXE)

$(EXE): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(MODDIR):
	mkdir -p $(MODDIR)

# ===========================================================================
# Python extension via f2py
# ===========================================================================
python: $(OBJS_PY) $(PYDIR)
	f2py -c $(SRCDIR)/python_wrapper.f90 $(OBJS_PY) \
	     -m reion_f \
	     --f90flags='$(FFLAGS_PY)'
	mv reion_f*.so $(PYDIR)/

$(PYDIR):
	mkdir -p $(PYDIR)

clean:
	rm -rf $(OBJDIR) $(MODDIR)

distclean: clean
	rm -f $(EXE)
