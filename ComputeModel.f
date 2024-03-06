	subroutine ComputeModel(recomputeopacities)
	use GlobalSetup
	IMPLICIT NONE
	logical recomputeopacities
	integer i
	
	TeffPoutput=TeffP
	modelfail=.false.
	if(do3D) then
		call Run3D(recomputeopacities)
	else
		call ComputeModel1D(recomputeopacities)
	endif

	if(doRing) then
	do i=1,nlam
		call RingFlux(Rplanet,Rstar,Tstar,TeffPoutput,tauRing,Rring,dRring,Dplanet,distance,lam(i),FRing(i))
		flux(0,i)=flux(0,i)+FRing(i)
	enddo
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
	if(doRing) computelam=.true.

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
			if(modelfail) return
			call SetupOpacities()
			if(modelfail) return
			if(nTiter.eq.1) then
				f=1d0
			else
				f=0.5d0
				if(forceEbalance) f=f+(1d0-f)*exp(-real(maxiter-nTiter)/5d0)
c				if(WaterWorld) f=f*(1d0-exp(-real(maxiter-nTiter)*3d0/real(maxiter)))
			endif
			if(f.gt.1d0) f=1d0
			call DoComputeT(Tconverged,f)
			if(modelfail) return
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
		if(.not.init3D) then
			computelam=.not.RTgridpoint
			if(forceEbalance) computelam=.true.
			call SetupStructure(.true.)
			if(modelfail) return
			call SetupOpacities()
			if(modelfail) return
			if(forceEbalance) then
				f=1d0
				nTiter=1
				call DoComputeT(Tconverged,f)
				if(modelfail) return
				computelam=.not.RTgridpoint
			endif
			par_tprofile=temp
		else
			return
		endif
	else
		call SetupStructure(computeopac)
		if(modelfail) return
		if(domakeai.and..not.modelsucces) return
		if(computeopac) call SetupOpacities()
		if(modelfail) return
	endif
	call cpu_time(stoptime)
	call output("Opacity computation: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
	if(.not.do3D) then
		call Raytrace()
	endif
	if(.not.do3D.and.(doRing.or.computeT)) then
		Lplanet=0d0
		do ilam=1,nlam
			if(RTgridpoint(ilam)) Lplanet=Lplanet+dfreq(ilam)*phase(1,0,ilam)
		enddo
		TeffPoutput=(Lplanet*distance**2*1e-23/(pi*Rplanet**2*((2d0*(pi*kb)**4)/(15d0*hplanck**3*clight**3))))**0.25d0
		call output("Teff: " // dbl2string(TeffPoutput,'(f10.2)') // "K" )
	else if(.not.do3D.and.emisspec.and..not.useobsgrid.and..not.retrieval) then
		Lplanet=0d0
		do i=1,nlam
			if(computelam(i).and.lamemis(i)) Lplanet=Lplanet+dfreq(i)*flux(0,i)
		enddo
		TeffPoutput=1000d0
		do i=1,10
			tot=0d0
			do ilam=1,nlam
				if(computelam(ilam).and.lamemis(ilam)) tot=tot+dfreq(ilam)*Planck(TeffPoutput,freq(ilam))
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



c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine RingFlux(Rp,Rs,Ts,Tp,tau,Rin,dR,Dp,Ds,lam,flux)
	use Constants
	IMPLICIT NONE
	real*8 Rp,Rs,Ts,Tp,tau,Rin,dR,Dp,Ds,flux
	real*8 T,lam,Planck,nu,A,Ra,Wp,Ws
	integer i,n
	parameter(n=100)
	real*8 R(n)

	nu=1d0/lam
	do i=1,n
		R(i)=exp(log(Rin*Rp)+log((Rin+dR)/Rin)*real(i-1)/real(n-1))
	enddo
	flux=0d0
	Ws=(1d0-sqrt(1d0-(Rs/Dp)**2))
	do i=1,n-1
		A=pi*(R(i+1)**2-R(i)**2)
		Ra=sqrt(R(i+1)*R(i))
		Wp=(1d0-sqrt(1d0-(Rp/Ra)**2))
		T=(Ts**4*Ws+Tp**4*Wp)**0.25
		flux=flux+A*Planck(T,nu)
	enddo		
	flux=flux*(1d0-exp(-tau))*1d23/Ds**2
	
	return
	end
	
