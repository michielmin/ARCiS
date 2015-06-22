	program main
	use GlobalSetup
	IMPLICIT NONE
	character*500 VersionGIT
	logical converged
	integer i
	real*8 starttime,starttime0,stoptime

	call cpu_time(starttime)
	starttime0=starttime

	call GetOutputDir
	open(unit=9,file=trim(outputdir) // "log.dat",RECL=6000)
	call output("Output dir: " // trim(outputdir))

	call output("==================================================================")
	call output("         SRON Planetary Atmosphere Retrieval Code - SPARC")
	call output("==================================================================")
c terms of use
	call output("By using SPARC you agree to the terms of use.")
	call output("It basically means you offer us co-author rights on any paper.")
	call output("that uses results computed with SPARC.")

	call output("==================================================================")
	call output("Let's get the show on the road!!")
	call output("SPARC version "//trim(VersionGIT()))
	call output("==================================================================")

	call Init()

	if(retrieval) call ReadObs()

	call cpu_time(stoptime)
	call output("Initialisation time: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	call output("==================================================================")
	starttime=stoptime

	converged=.false.
	do while(.not.converged)
		call SetupStructure()
		call SetupOpacities()
		call cpu_time(stoptime)
		call output("Opacity computation: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
		do i=1,nobs
			call Raytrace(i)
		enddo
		if(retrieval) then
			call AdjustParameters(converged)
		else
			converged=.true.
		endif
		call cpu_time(stoptime)
		call output("Model runtime:       " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
		starttime=stoptime
	enddo

	call WriteOutput()

	call cpu_time(stoptime)
	call output("==================================================================")
	call output("Total runtime:       " // trim(dbl2string((stoptime-starttime0),'(f10.2)')) // " s")
	call output("All done!")
	call output("==================================================================")

	end
	
