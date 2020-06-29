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
	call output("=!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=")
	call output("==================================================================")
	call output("THIS VERSION OF ARCIS IS MEANT FOR REPRODUCTION OF THE RESULTS IN")
	call output("")
	call output("                     MIN ET AL. 2020")
	call output("      The ARCiS framework for Exoplanet Atmospheres")
	call output("           Modelling Philosophy and Retrieval")
	call output("")
	call output("THE AUTHORS STRONGLY ADVICE CONTACTING US FOR USE OF THE CODE IN")
	call output("ANY OTHER WAY. WE ACCEPT NO RESPONSIBILITY FOR RESULTS PUBLISHED")
	call output("WITH THIS CODE WITHOUT CONSULTATION OF THE AUTHORS.")
	call output("")
	call output("                                        Michiel Min: M.Min@sron.nl")
	call output("==================================================================")
	call output("=!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=")
	call output("==================================================================")
	call output("Let's get the show on the road!!")
	call output("ARCiS version "//trim(VersionGIT()))
	call output("==================================================================")

	call Init()

	if(nobs.ne.0) call ReadObs()

	call cpu_time(stoptime)
	call output("Initialisation time: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	call output("==================================================================")

	if(dopostequalweights) then
		call PostEqualWeights()
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
	logical recomputeopacities
	
	call ComputeModel1D(recomputeopacities)
	
	return
	end


	subroutine ComputeModel1D(recomputeopacities)
	use GlobalSetup
	IMPLICIT NONE
	logical Tconverged
	real*8 starttime,stoptime,starttime_w,stoptime_w,omp_get_wtime,f
	logical recomputeopacities
	logical computeopac,temp
	
	computeopac=recomputeopacities

	call cpu_time(starttime)
	Tconverged=.false.
	nTiter=0
	call SetupStructure(computeopac)
	if(domakeai.and..not.modelsucces) return
	if(computeopac) call SetupOpacities()
	call cpu_time(stoptime)
	call output("Opacity computation: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	if(.not.do3D) call Raytrace()
	call cpu_time(stoptime)
	call output("Model runtime:       " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")

	return
	end
