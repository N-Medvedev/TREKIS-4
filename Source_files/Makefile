# 0000000000000000000000000000000000000000000000000000000000000
# This file is part of TREKIS-4
# available at: https://github.com/N-Medvedev/TREKIS-4
# The makefile compiles (type make) and makes an executable called TREKIS.x
# This file was written by N.Medvedev 
# in 2018-2026
#----------------------------------------------------- 

EXEPATH = .

ifeq ($(db),y)
	EXE = $(EXEPATH)/TREKIS_DEBUG.x
	exe_name = /TREKIS_DEBUG.x
else
	EXE = $(EXEPATH)/TREKIS.x
	exe_name = /TREKIS.x
endif
LD = ld

# choose the compiles among gfortran, gfortran6.1.0, ifort2011, ifort2013
ifeq ($(c),gf)
	F90 = gfortran
	F90FLAGS = -ffree-line-length-0 -fdec -fdec-format-defaults -cpp -fdefault-real-8 -fdefault-double-8 -O2
	LDFLAGS += -fopenmp
	ifeq (${OMP},y)
		F90FLAGS += -fopenmp
	endif
else ifeq ($(c),gf6)
	F90 = gfortran6.1.0
	F90FLAGS = -ffree-line-length-0 -fdec -fdec-format-defaults -cpp -fdefault-real-8 -fdefault-double-8 -O2
	LDFLAGS += -fopenmp
	ifeq (${OMP},y)
		F90FLAGS += -fopenmp
	endif
else ifeq ($(c),if13)
	F90 = ifort2013

	ifneq (${OMP},) # OMP=no means no omp
#		F90FLAGS = -mkl
	else  # there is omp by default
		F90FLAGS = -qopenmp
#		F90FLAGS += -mkl=parallel
	endif

	ifeq ($(db),y)
		# flags for debugging
		#F90FLAGS += -g -fbounds-check
		F90FLAGS +=  -real-size 64 -debug all -check all -fpe0 -fp-stack-check -O0 -g -fp-model precise -traceback -gen-interfaces -warn interfaces -fpp
	else
		# Flags for maximum performance
		F90FLAGS += -O5 -fpp -ipo -real-size 64
	endif

	F90FLAGS += -standard-semantics -assume nofpe_summary
else
	F90 = ifort

	F90FLAGS = -fpp -cpp

	ifneq (${OMP},) # OMP=no means no omp
#		F90FLAGS += -qmkl
	else  # there is omp by default
		F90FLAGS += -qopenmp
#		F90FLAGS += -qopenmp -mkl=parallel
	endif

	ifeq ($(db),y)
		# flags for debugging
		#F90FLAGS += -g -fbounds-check
		F90FLAGS += -debug all -check all -fpe0 -fp-stack-check -O0 -g -fp-model precise -traceback -gen-interfaces -warn interfaces
	else
		# Flags for maximum performance
		F90FLAGS += -O5 -ipo -real-size 64
	endif

	F90FLAGS += -standard-semantics -assume nofpe_summary
endif

# list the files necessary for the compilation in the order of their dependencies
OBJS = Universal_constants.o Geometries.o Little_subroutines.o Objects.o Variables.o Gnuplotting.o Periodic_table.o Dealing_with_files.o Dealing_with_EADL.o Dealing_with_cdf_files.o Read_numerical_parameters.o Relativity.o SHI_charge_state.o CS_general_tools.o CS_integration_limits.o CDF_Ritchi.o CDF_Mermin.o CDF_delta.o Dealing_with_XYZ_files.o MD_general_tools.o Dealing_with_LAMMPS.o Read_MD_parameters.o Read_input_data.o Dealing_with_DOS.o Initial_conditions.o CDF_get_from_data.o CS_photons_pair_creation.o CS_photons_Compton.o CS_photons_Rayleigh.o CS_electrons_inelastic.o CS_electrons_elastic.o CS_electrons_Bremsstrahlung.o CS_ion_inelastic.o CS_positrons_inelastic.o CS_positrons_elastic.o CS_positrons_Bremsstrahlung.o CS_positrons_annihilation.o CS_holes_elastic.o CS_holes_inelastic.o CS_muon_inelastic.o CS_muon_elastic.o CS_muon_Bremsstrahlung.o Output.o MC_general_tools.o MC_data_analysis.o MC_photon.o MC_electron.o MC_positron.o MC_hole.o MC_SHI.o MC_muon.o MC.o MD_data_analysis.o MD_Pot_Simple.o MD_Pot_ZBL.o MD_Pot_Buck.o MD_Pot_SW.o MD_Pot_Coulomb.o MD.o TREKIS_main_file.o


# print explanations
default:
	@echo "**********************************************"
ifeq ($(c),gf)
	@echo "Compiling TREKIS-4 with gfortran"
	@echo "----------------------------------------------"
	@echo "Run the code as: ." $(exe_name)
	@echo "----------------------------------------------"
else ifeq ($(c),gf6)
	@echo "Compiling TREKIS-4 with gfortran 6.1.0"
	@echo "----------------------------------------------"
	@echo "Run the code as: ." $(exe_name)
	@echo "----------------------------------------------"
else ifeq ($(c),if13)
	@echo "Compiling TREKIS-4 with ifort 2013"
	@echo "----------------------------------------------"
	@echo "Run the code as: ." $(exe_name)
	@echo "----------------------------------------------"
else
	@echo "Compiling TREKIS-4 with ifort"
	@echo "----------------------------------------------"
	@echo "Run the code as: ." $(exe_name)
	@echo "----------------------------------------------"
endif


# autoclean before the compilation
#	rm -f *.o
#	rm -f *.mod
	$(MAKE) $(EXE)
# autoclean after the compilation
#	rm *.o
#	rm *.mod

# compile and load
$(EXE):	$(OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $(EXE)  $(OBJS)

.SUFFIXES: .f90 .o .mod
.f90.o:
	$(F90) $(F90FLAGS) -c $*.f90

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.yaml
	rm -f *.optrpt
	rm -f *.x

