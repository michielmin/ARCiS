# makefile for OpacityTool (with comments!)
# MacOSX 10.10.1 gfortran gcc version 5.0.0 20141005
#

# compiler= FC, flags = FFlags
# linker= LINKER, flags= LDFLAGS, libraries=LIBS

FC	  = ifort
LINKER  = ifort

ifeq ($(gfort),true)
	FC	  = gfortran
	LINKER	  = gfortran
endif

# array boundary check
ifeq ($(debug),true)
  ifeq ($(gfort),true)
    DEBUGGING = -Og -Wall -Wno-unused-variable -fbounds-check -fbacktrace -fcheck=all -ffpe-trap=zero,overflow -g3
  else	
    DEBUGGING = -check all -g -traceback -check bounds -check uninit -O0
  endif
endif

# enforce multi core compilation with:
# cl> make multi=true
ifeq ($(multi),true)
	ifeq ($(gfort),true)
		MULTICORE = -fopenmp -DUSE_OPENMP
	else
		MULTICORE = -qopenmp -parallel -fp-model strict -DUSE_OPENMP
		ifeq ($(shell uname),Linux)
			MULTICORE = -qopenmp -fp-model strict -DUSE_OPENMP
		endif
	endif
endif

ifeq ($(multinest),false)
  LIBS_MN      = -DNO_MULTINEST
else
  LIBS_MN      = -lmultinest -DUSE_MULTINEST
endif

ifeq ($(mcmc),true)
  LIBS_MCMC      = -lmcmcrun -framework Accelerate -DUSE_MCMCF90
else
  LIBS_MCMC      = -DNO_MCMCF90
endif

# Platform specific compilation options
ifeq ($(gfort),true)
  FLAG_ALL      = -O5 -finit-local-zero $(MULTICORE) -I$(HOME)/include -I/usr/local/modules $(LIBS_MN) $(LIBS_MCMC) -frecursive -finit-derived -Wuninitialized
  FLAG_LINUX    = -ffixed-line-length-none -cpp -malign-double
  FLAG_MAC      = -m64 -ffixed-line-length-none -cpp -malign-double
else
  FLAG_ALL      = -O3 -g -extend-source -zero -prec-div $(MULTICORE) -assume buffered_io -I/usr/local/modules -fp-model strict -heap-arrays $(LIBS_MN) $(LIBS_MCMC)
  FLAG_LINUX    = -xHOST -fpp
  FLAG_MAC      = -xHOST -qopt-prefetch -static-intel -fpp -heap-arrays
endif

LIBS_FITS		= -lcfitsio

ifeq ($(shell uname),Linux)
  FFLAGS   = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS) -diag-disable vec $(DEBUGGING) 
  LDFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS) -I$(HOME)/include $(DEBUGGING) 
  LIBS     = -L$(HOME)/lib -lm $(LIBS_FITS) $(LIBS_MN) -llapack Version.f 
else
  FFLAGS  = $(FLAG_ALL) $(FLAG_MAC) $(FLAG_FITS) $(DEBUGGING)  
  LDFLAGS = $(FLAG_ALL) $(FLAG_MAC) $(FLAG_FITS) $(DEBUGGING) 
  LIBS    =  -L/usr/local/lib $(LIBS_FITS) $(LIBS_MN) -llapack Version.f
endif


# files to make
OBJS	= Modules.o \
		InputOutput.o \
		Main.o \
		easy_chem_extra.o \
		easy_chem.o \
		ComputeT.o \
		DiffuseCloud.o \
		diseq_diffusion.o \
		diseq_rate.o \
		diseq_timescale.o \
		diseq_calc.o \
		OpacityFITS.o \
		SetupStructure.o \
		SetupOpacities.o \
		Raytrace.o \
		WriteOutput.o \
		Subroutines.o \
		CIA.o \
		ComputePart.o \
		SetupCloud.o \
		ReadParticleFits.o \
		RefIndData.o \
		EnstatiteData.o \
		ForsteriteData.o \
		BrookiteData.o \
		WaterData.o \
		CorrundumData.o \
		SiO2Data.o \
		SiOData.o \
		FeOData.o \
		Mg0.6Fe0.4OData.o \
		MgOData.o \
		IronData.o \
		SiCData.o \
		OrganicsHenning.o \
		Soot.o \
		TholinData.o \
		LabradoriteData.o \
		H2SO4.o \
		KuruczData.o \
		MCRad.o \
		Genetic.o \
		MakeAI.o \
		dlsei.o \
		Lapack.o \
		writeFITS.o \
		params_multinest.o \
		MultiNestARCiS.o \
		MCMC_ARCiS.o \
		mrqmin.o \
		truncated_normal.o \
		amoeba.o \
		ggchem/datamod.o \
		ggchem/database.o \
		ggchem/supersat.o \
		ggchem/stindex.o \
		ggchem/smchem16.o \
		ggchem/smchem8.o \
		ggchem/nucleation.o \
		ggchem/ggchem.o \
		ggchem/init_dustchem.o \
		ggchem/init_chemistry.o \
		ggchem/gauss_nm.o \
		ggchem/gauss8.o \
		ggchem/gauss16.o \
		ggchem/equil_cond.o \
		ggchem/upper.o \
		GGchemARCiS.o \
		nasa_polynomial.o \
		GGchem_linpack_q.o \
		GGchem_is_nan.o \
		Init.o \
		Run3D.o \
		ReadFull3D.o \
		Retrieval.o \
		PostEqualWeights.o \
		LightCurve.o \
		ConvertColors.o \
		Formation.o \
		ComputePAH.o

# program name and install location
PROGRAM       = ARCiS
DEST	      = ${HOME}/bin

# make actions 
all:		$(PROGRAM)
version:;	echo "#define gitversion \"$(shell git rev-parse HEAD)\"" > gitversion.h
clean:;		rm -f $(OBJS) $(PROGRAM) *.mod *.i *.i90
install:	$(PROGRAM)
		mv $(PROGRAM) $(DEST)

# how to compile program 
.SUFFIXES : .o .f .f90 .F

.f.o:
	$(FC) $(LDFLAGS) -c $< -o $@ 

.f90.o:
	$(FC) $(LDFLAGS) -c $< -o $@ 

.F.o:
	$(FC) $(LDFLAGS) -c $< -o $@ 

$(PROGRAM):     version  $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

# recompile everything if Modules.f has changed 
$(OBJS):	Modules.f

# recompile everything if InputOutput.f has changed 
$(OBJS):	InputOutput.f


