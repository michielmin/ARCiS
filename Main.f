	program main
	use GlobalSetup
	IMPLICIT NONE
	character*500 VersionGIT
	integer i
	real*8 starttime,stoptime

	call SetOutputMode(.true.)
	
	call cpu_time(starttime)

	call GetOutputDir
	open(unit=9,file=trim(outputdir) // "log.dat",RECL=6000)
	call output("Output dir: " // trim(outputdir))

	call output("==================================================================")
	call output("         ARtful modelling code for exoplanet Science - ARCiS")
	call output("==================================================================")
c terms of use
	call output("By using ARCiS you agree to the terms of use.")
	call output("It basically means you offer us co-author rights on any paper")
	call output("that uses results computed with ARCiS.")

	call output("==================================================================")
	call output("Let's get the show on the road!!")
	call output("ARCiS version "//trim(VersionGIT()))
	call output("==================================================================")

	call Init()

	if(retrieval) call ReadObs()

	call cpu_time(stoptime)
	call output("Initialisation time: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	call output("==================================================================")

	if(opacitymode) then
		if(dochemistry) then
			call OnlyChemCompute
		else
			do i=1,nmol
				mixrat_r(:,i)=mixrat(i)
			enddo
		endif
		call SetupOpacities()
		call WriteOpacityFITS()
	else if(domakeai) then
		call MakeAI()
	else if(retrieval) then
		call DoRetrieval()
	else
		call ComputeModel(.true.)
		call WriteStructure()
		call WriteOutput()
	endif

	call cpu_time(stoptime)
	call output("==================================================================")
	call output("Total runtime:       " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	call output("All done!")
	call output("==================================================================")

	end
	

	subroutine ComputeModel(recomputeopacities)
	use GlobalSetup
	IMPLICIT NONE
	logical Tconverged
	integer nTiter
	real*8 starttime,stoptime,starttime_w,stoptime_w,omp_get_wtime
	logical,intent(in),optional :: recomputeopacities
	logical computeopac
	
	if(present(recomputeopacities)) then
		computeopac=recomputeopacities
	else
		computeopac=.true.
	endif

	call cpu_time(starttime)
#if USE_OPENMP
	starttime_w=omp_get_wtime()
#ENDIF
	Tconverged=.false.
	call SetupStructure(computeopac)
	if(domakeai.and..not.modelsucces) return
	if(computeopac) call SetupOpacities()
	if(computeT) then
		nTiter=0
		do while(.not.Tconverged.and.nTiter.lt.maxiter)
			call DoComputeT(Tconverged)
			par_tprofile=.false.
			call SetupStructure(.true.)
			call SetupOpacities()
			nTiter=nTiter+1
		enddo
	endif
	call cpu_time(stoptime)
	call output("Opacity computation: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
#if USE_OPENMP
	stoptime_w=omp_get_wtime()
	call output("Walltime:            " // trim(dbl2string((stoptime_w-starttime_w),'(f10.2)')) // " s")
#ENDIF
	call Raytrace()
	call cpu_time(stoptime)
	call output("Model runtime:       " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")

	return
	end
