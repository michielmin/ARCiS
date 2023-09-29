	subroutine ComputeModel(recomputeopacities)
	use GlobalSetup
	IMPLICIT NONE
	logical recomputeopacities
	
	if(do3D) then
		call Run3D(recomputeopacities)
	else
		call ComputeModel1D(recomputeopacities)
	endif
	
	return
	end

	
	subroutine ComputeModel1D(recomputeopacities)
	use GlobalSetup
	use TimingModule
	use Constants
	IMPLICIT NONE
	logical Tconverged
	real*8 starttime,stoptime,starttime_w,stoptime_w,omp_get_wtime,f,tot,Planck
	logical recomputeopacities
	logical computeopac,temp
	integer i,ilam
	
	computeopac=recomputeopacities

	call cpu_time(starttime)
	Tconverged=.false.
	nTiter=0
	Tsurface=sqrt(Rstar/(2d0*Dplanet))*Tstar
	if(computeT.and.computeopac) then
		temp=par_tprofile
		par_tprofile=.false.
		f=1d0
		computelam=RTgridpoint
		do nTiter=1,maxiter
			call output("Temperature computation (" // trim(int2string(nTiter,'(i3)')) // " of " 
     &					// trim(int2string(maxiter,'(i3)')) // ")")

			if(fixnight2day.and..not.do3D) call ComputeNight2Day((nTiter.le.1))
			call SetupStructure(.true.)
			call SetupOpacities()
			if(nTiter.eq.1) then
				f=1d0
			else
				f=0.5d0
				if(forceEbalance) f=f+0.5d0*exp(-real(maxiter-nTiter)/5d0)
				if(WaterWorld) f=0.5d0
			endif
			if(f.gt.1d0) f=1d0
			call DoComputeT(Tconverged,f)
			if(Tconverged.and.nTiter.ge.miniter) exit
c			call SetoutputMode(.true.)
c			call output("Chemistry cpu time: " // trim(dbl2string(timechem,'(f10.4)')) // " s")
c			call output("Chemistry walltime: " // trim(dbl2string(dble(itimechem)/dble(rate),'(f10.4)')) // " s")
c			call output("Number of chemistry calls: " // trim(int2string(ctimechem,'(i5)')))
c			call output("Cloud cpu time:    " // trim(dbl2string(timecloud,'(f10.4)')) // " s")
c			call output("Cloud walltime:    " // trim(dbl2string(dble(itimecloud)/dble(rate),'(f10.4)')) // " s")
c			call output("Number of cloud calls:     " // trim(int2string(ctimecloud,'(i5)')))
c			call output("PTstruct cpu time: " // trim(dbl2string(timetemp,'(f10.4)')) // " s")
c			call output("PTstruct walltime: " // trim(dbl2string(dble(itimetemp)/dble(rate),'(f10.4)')) // " s")
c			call output("Number of PTstruct calls:  " // trim(int2string(ctimetemp,'(i5)')))
c			call SetoutputMode(.false.)
		enddo
		if(.not.do3D.or..not.init3D) then
			computelam=.not.RTgridpoint
			if(forceEbalance) computelam=.true.
			call SetupStructure(.true.)
			call SetupOpacities()
			if(forceEbalance) then
				f=1d0
				nTiter=1
				call DoComputeT(Tconverged,f)
				computelam=.not.RTgridpoint
			endif
			par_tprofile=temp
		endif
	else
		call SetupStructure(computeopac)
		if(domakeai.and..not.modelsucces) return
		if(computeopac) call SetupOpacities()
	endif
	call cpu_time(stoptime)
	call output("Opacity computation: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	if(.not.do3D) then
		call Raytrace()
	endif
	if(emisspec.and..not.useobsgrid.and..not.retrieval) then
		Lplanet=0d0
		do i=1,nlam
			if(computelam(i).and.lamemis(i)) Lplanet=Lplanet+dfreq(i)*flux(0,i)
		enddo
		TeffPoutput=1000d0
		do i=1,10
			tot=0d0
			do ilam=1,nlam
				tot=tot+dfreq(ilam)*Planck(TeffPoutput,freq(ilam))
			enddo
			tot=tot*pi*Rplanet**2*1d23/distance**2
			TeffPoutput=TeffPoutput*(Lplanet/tot)**0.25
		enddo
		call output("Teff: " // dbl2string(TeffPoutput,'(f10.2)') // "K" )
	endif

	call cpu_time(stoptime)
	call output("Model runtime:       " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")

	return
	end
