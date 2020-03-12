# fortran compiler
# FC90 = gfortran
FC90 = ifort

# OpenMP-compiled library (intel)
LIB = -L/Users/Michael/phreeqcRM-3.4-OpenMP-intel/lib -lphreeqcrm

# optimization level
# FLAGS = -O3

# for debugging with ifort
FLAGS = -O0 -g -traceback -check all
# generate warnings and direct to compile.log
WARN = -warn all > compile.log 2>&1

# for debugging with GNU
# FLAGS = -g -fcheck=all -fbacktrace -Og -Wall > compile.log 2>&1
# generate warnings and direct to compile.log
# WARN = -Wall > compile.log 2>&1

# add LIB to FLAGS
FLAGS += $(LIB)
# add WARN to FLAGS
FLAGS += $(WARN)

miRPT_hMetal_files = miRPT_hMetal.f90 mod_miRPT_hMetal.f90 RM_interface.F90 kdtree2.f90

# executable name
executable = exe_hMetal_miRPT
# directory in which to place the binary
binDir = ./
# executable with path appended
dirExecutable = $(addprefix $(binDir),$(executable))
all: $(dirExecutable)

$(dirExecutable): $(miRPT_hMetal_files)
	$(FC90) RM_interface.F90 -c $(FLAGS)
	$(FC90) kdtree2.f90 -c $(FLAGS)
	$(FC90) mod_miRPT_hMetal.f90 -c $(FLAGS)
	$(FC90) $(miRPT_hMetal_files) -o $@ $(FLAGS)

clean:
	rm -r exe_* *.mod *.o *.out *.err myrun* *.dSYM *.inp compile.log
