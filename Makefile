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
    DEBUGGING = -check all -traceback -check bounds -O0 -g
  endif
endif

# enforce multi core compilation with:
# cl> make multi=true
ifeq ($(multi),true)
	MULTICORE = -openmp -fp-model strict -DUSE_OPENMP
	ifeq ($(debug),true)
  		MULTICORE = -openmp -DUSE_OPENMP
	endif
endif


# Platform specific compilation options
ifeq ($(gfort),true)
  FLAG_ALL      = -O3 -g $(DEBUGGING) $(MULTICORE) -lgfortran
  FLAG_LINUX    = -cpp
  FLAG_MAC      = -m64 -ffixed-line-length-132 -cpp
else
  FLAG_ALL      = -O3 -g -extend-source -zero -prec-div $(DEBUGGING) $(MULTICORE) -assume buffered_io -fp-model precise
  FLAG_LINUX    = -xHOST -fpp
  FLAG_MAC      = -xHOST -opt-prefetch -static-intel -fpp 
endif

LIBS_FITS		= -lcfitsio

ifeq ($(shell uname),Linux)
  FFLAGS   = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS) -diag-disable vec
  LDFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS)
  LIBS     = -L$(HOME)/lib -lm $(LIBS_FITS)
else
  FFLAGS  = $(FLAG_ALL) $(FLAG_MAC) $(FLAG_FITS)
  LDFLAGS = $(FLAG_ALL) $(FLAG_MAC) $(FLAG_FITS)
  LIBS    =  -L/usr/local/lib $(LIBS_FITS)
endif


# files to make
OBJS	= Modules.o \
		InputOutput.o \
		Main.o \
		Init.o \
		SetupStructure.o \
		SetupOpacities.o \
		Raytrace.o \
		WriteOutput.o \
		Version.o \
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
		IronData.o \
		KuruczData.o \
		MCRad.o \
		OpacityFITS.o \
		ComputeT.o \
		polyPartition.o \
		Retrieval.o \
		Genetic.o \
		MakeAI.o \
		dlsei.o \
		Lapack.o \
		easy_chem_extra.o \
		easy_chem.o

# program name and install location
PROGRAM       = SPARC
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
	$(FC) $(LDFLAGS) -c $<

.f90.o:
	$(FC) $(LDFLAGS) -c $<

.F.o:
	$(FC) $(LDFLAGS) -c $<

$(PROGRAM):     version  $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

# recompile everything if Modules.f has changed 
$(OBJS):	Modules.f


