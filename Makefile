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
    DEBUGGING = -fbounds-check -fbacktrace -fcheck=all
  else	
    DEBUGGING = -check all -g -traceback -check bounds -check uninit -O0
  endif
endif

# enforce multi core compilation with:
# cl> make multi=true
ifeq ($(multi),true)
	ifeq ($(gfort),true)
		MULTICORE = -openmp -DUSE_OPENMP
	else
		MULTICORE = -openmp -fp-model strict -DUSE_OPENMP
		ifeq ($(shell uname),Linux)
			MULTICORE = -qopenmp -fp-model strict -DUSE_OPENMP
		endif
	endif
	ifeq ($(debug),true)
		MULTICORE = -openmp -DUSE_OPENMP
		ifeq ($(shell uname),Linux)
			MULTICORE = -qopenmp -DUSE_OPENMP
		endif
	endif
endif


# Platform specific compilation options
ifeq ($(gfort),true)
  FLAG_ALL      = -O5 -g $(MULTICORE) -lgfortran -I$(HOME)/include
  FLAG_LINUX    = -ffixed-line-length-132 -cpp
  FLAG_MAC      = -m64 -ffixed-line-length-132 -cpp
else
  FLAG_ALL      = -O3 -g -extend-source -zero -prec-div $(MULTICORE) -assume buffered_io -I/usr/local/modules -fp-model strict -heap-arrays
  FLAG_LINUX    = -xHOST -fpp
  FLAG_MAC      = -xHOST -opt-prefetch -static-intel -fpp -heap-arrays 
endif

LIBS_FITS		= -lcfitsio

ifeq ($(shell uname),Linux)
  FFLAGS   = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS) -diag-disable vec $(DEBUGGING) 
  LDFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS) -I$(HOME)/include $(DEBUGGING) 
  LIBS     = -L$(HOME)/lib -lm $(LIBS_FITS) -llapack -lmultinest Version.f 
else
  FFLAGS  = $(FLAG_ALL) $(FLAG_MAC) $(FLAG_FITS) $(DEBUGGING) 
  LDFLAGS = $(FLAG_ALL) $(FLAG_MAC) $(FLAG_FITS) $(DEBUGGING) 
  LIBS    =  -L/usr/local/lib $(LIBS_FITS) -lmultinest Version.f
endif


# files to make
OBJS	= Modules.o \
		InputOutput.o \
		Main.o \
		easy_chem_extra.o \
		easy_chem.o \
		ComputeT.o \
		DiffuseCloud.o \
		Init.o \
		diseq_diffusion.o \
		diseq_rate.o \
		diseq_timescale.o \
		diseq_calc.o \
		OpacityFITS.o \
		SetupStructure.o \
		SetupOpacities.o \
		Raytrace.o \
		WriteOutput.o \
		AdjustParameters.o \
		ReadData.o \
		Subroutines.o \
		Voigt.o \
		TIPS_2011_v1p0_sub.o \
		CIA.o \
		ComputePart.o \
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
		KuruczData.o \
		MCRad.o \
		polyPartition.o \
		Retrieval.o \
		Genetic.o \
		MakeAI.o \
		dlsei.o \
		Lapack.o \
		writeFITS.o \
		params_multinest.o \
		MultiNestARCiS.o \
		mrqmin.o \
		truncated_normal.o \
		amoeba.o \
		GGchemARCiS.o \
		nasa_polynomial.o \
		GGchem_linpack_q.o \
		GGchem_is_nan.o \
		PostEqualWeights.o \
		TrendCompute.o \
		MCComputeT.o \
		LightCurve.o \
		Run3D.o

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


