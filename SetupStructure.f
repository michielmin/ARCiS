	subroutine SetupStructure(compute_mixrat)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 g,dp,dz,dlogp,RgasBar,Mtot,Pb(nr+1),tot
	parameter(RgasBar=82.05736*1.01325)
	integer i,imol,nmix,j,niter
	logical ini,compute_mixrat
	character*40 cloudspecies(max(nclouds,1))
	
	do i=1,nclouds
		cloudspecies(i)=Cloud(i)%species
	enddo

	ini = .TRUE.

	niter=1
	if(dochemistry) niter=2
	if(par_tprofile) niter=2

	Pb(1)=P(1)
	do i=2,nr
		Pb(i)=sqrt(P(i-1)*P(i))
	enddo
	Pb(nr+1)=P(nr)

	do j=1,niter

	do i=1,nr
		tot=0d0
		do imol=1,nmol
			if(mixrat_r(i,imol).gt.0d0) tot=tot+mixrat_r(i,imol)
		enddo
		mixrat_r(i,1:nmol)=mixrat_r(i,1:nmol)/tot
	enddo

	Mtot=Mplanet
	g=Ggrav*Mtot/(Rplanet)**2
	call output("log(g) [cgs]: " // dbl2string(log10(g),'(f8.3)'))

	R(1)=Rplanet
	mu=0d0
	do imol=1,nmol
		if(mixrat_r(1,imol).gt.0d0) mu=mu+mixrat_r(1,imol)*Mmol(imol)
	enddo
	call output("Mean molecular weight: " // dbl2string(mu,'(f8.3)'))
	do i=1,nr
		mu=0d0
		do imol=1,nmol
			if(mixrat_r(i,imol).gt.0d0) mu=mu+mixrat_r(i,imol)*Mmol(imol)
		enddo

		g=Ggrav*Mtot/R(i)**2
		
		dp=Pb(i)-Pb(i+1)
		dlogp=log(Pb(i)/Pb(i+1))

		Ndens(i)=P(i)*Avogadro/(RgasBar*T(i))
		Hp(i)=(T(i)*kb)/(g*mp*mu)
		dz=dlogp*Hp(i)
		dens(i)=Ndens(i)*mp*mu
		R(i+1)=R(i)+dz
		Mtot=Mtot+dens(i)*(R(i+1)**3-R(i)**3)*4d0*pi/3d0
	enddo

	if(par_tprofile) call ComputeParamT(T)
	do i=1,nr
		if(T(i).gt.3000d0) T(i)=3000d0
		if(T(i).lt.100d0) T(i)=100d0
	enddo

	if(dochemistry.and.j.ne.niter) then
	if(compute_mixrat) then
		if(ini) call easy_chem_set_molfracs_atoms(COratio,metallicity)
		call output("==================================================================")
		call output("Computing chemistry using easy_chem by Paul Molliere")
		do i=1,nr
			call tellertje(i,nr)
			if(P(i).ge.mixP.or.i.eq.1) then
				call call_easy_chem(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,.false.,cloudspecies,XeqCloud(i,1:nclouds),nclouds)
			else
				mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
			endif
		enddo
		call output("==================================================================")
		mixrat_old_r=mixrat_r
	else
		mixrat_r=mixrat_old_r
	endif
	endif

	enddo

	if(dochemistry) then
	if(compute_mixrat) then
		ini=.true.
		call easy_chem_set_molfracs_atoms(COratio,metallicity)
		call output("==================================================================")
		call output("Computing chemistry using easy_chem by Paul Molliere")
		do i=1,nr
			call tellertje(i,nr)
c			do j=1,nmol
c				call MorleyChemistry(mixrat_r(i,j),T(i),P(i),molname(j),metallicity)
c			enddo
			if(P(i).ge.mixP.or.i.eq.1) then
				call call_easy_chem(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,XeqCloud(i,1:nclouds),nclouds)
			else
				if(cloudcompute) call call_easy_chem(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,XeqCloud(i,1:nclouds),nclouds)
				mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
			endif
		enddo
		call output("==================================================================")
		mixrat_old_r=mixrat_r
	else
		mixrat_r=mixrat_old_r
	endif
	endif

	do i=1,nclouds
		call SetupCloud(i)
	enddo

	return
	end

	subroutine ComputeParamT(x)
	use GlobalSetup
	IMPLICIT NONE
	real*8 x(nr),tau,Tirr,eta,expint
	integer i

	tau=0d0
	Tirr=betaT*sqrt(Rstar/(2d0*Dplanet))*Tstar
	do i=nr,1,-1
		x(i)=(3d0*TeffP**4/4d0)*(2d0/3d0+tau)
		if(tau.lt.100d0) then
			eta=2d0/3d0+(2d0/(3d0*gammaT1))*(1d0+(gammaT1*tau/2d0-1)*exp(-gammaT1*tau))
     &					+(2d0*gammaT1/3d0)*(1d0-tau**2/2d0)*expint(2,gammaT1*tau)
		else
			eta=2d0/3d0+(2d0/(3d0*gammaT1))
		endif
		x(i)=x(i)+(3d0*Tirr**4/4d0)*(1d0-alphaT)*eta
		if(tau.lt.100d0) then
			eta=2d0/3d0+(2d0/(3d0*gammaT2))*(1d0+(gammaT2*tau/2d0-1)*exp(-gammaT2*tau))
     &					+(2d0*gammaT2/3d0)*(1d0-tau**2/2d0)*expint(2,gammaT2*tau)
		else
			eta=2d0/3d0+(2d0/(3d0*gammaT2))
		endif
		x(i)=x(i)+(3d0*Tirr**4/4d0)*alphaT*eta
		x(i)=x(i)**0.25d0
		if(x(i).gt.10000d0) x(i)=10000d0
		tau=tau+kappaT*dens(i)*(R(i+1)-R(i))
	enddo

	return
	end

	subroutine WriteStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 mix(nr,nmol)
	integer i,imol,nmix
	character*1000 form
	character*10 namemix(nmol)
	
	open(unit=50,file=trim(outputdir) // 'densityprofile.dat',RECL=100)
	write(50,'("#",a14,a15,a15,a13,a10,a10)') "radius [cm]","height [cm]","dens [g/cm^3]","N [1/cm^3]","T [K]","P [Ba]"
	do i=1,nr
		write(50,'(es15.7,es15.4,es15.4,es13.4,f10.3,es10.3)') sqrt(R(i)*R(i+1)),sqrt(R(i)*R(i+1))-Rplanet
     &			,dens(i),Ndens(i),T(i),P(i)
	enddo
	close(unit=50)

	open(unit=50,file=trim(outputdir) // 'density.dat',RECL=1000)
	do i=1,nr
		write(50,*) sqrt(R(i)*R(i+1)),sqrt(R(i)*R(i+1))-Rplanet
     &			,dens(i),cloud_dens(i,1:nclouds)
	enddo
	close(unit=50)

	nmix=0
	do imol=1,nmol
		if(includemol(imol)) then
			nmix=nmix+1
			mix(1:nr,nmix)=mixrat_r(1:nr,imol)
			namemix(nmix)=molname(imol)
		endif
	enddo

	open(unit=50,file=trim(outputdir) // 'mixingratios.dat',RECL=1000)
	form='("#",a9,a10,' // trim(int2string(nmix,'(i2)')) // 'a15)'
	write(50,form) "T [K]","P [Ba]",namemix(1:nmix)
	form='(f10.3,es10.3,' // trim(int2string(nmix,'(i2)')) // 'es15.4)'
	do i=1,nr
		write(50,form) T(i),P(i),mix(i,1:nmix)
	enddo
	close(unit=50)

	return
	end
	
	

	subroutine SetupCloud(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 pp,tot,column
	integer ii,i,j,nsubr
	real*8 Xc,Xc1,lambdaC
	
	if(.not.cloudcompute) then
		column=0d0
		nsubr=100
		do i=1,nr
			do j=1,nsubr
				if(i.ne.nr) then
					pp=(log10(P(i))+log10(P(i+1)/P(i))*real(j)/real(nsubr+1))
				else
					pp=log10(P(i))
				endif
				cloud_dens(i,ii)=cloud_dens(i,ii)+exp(-(abs(pp-log10(Cloud(ii)%P))/log10(Cloud(ii)%dP))**Cloud(ii)%s/2d0)
			enddo
			column=column+cloud_dens(i,ii)*(R(i+1)-R(i))
		enddo
		cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*Cloud(ii)%column/column
	else
c use Ackerman & Marley 2001 cloud computation
		column=0d0
		nsubr=1
		Xc1=0d0
		cloud_dens(1,ii)=Xc1*dens(1)
		do i=2,nr
			lambdaC=0.1d0
			do j=1,nsubr
				Xc=(Xc1+(XeqCloud(i,ii)-XeqCloud(i-1,ii))/real(nsubr))/(1d0-Cloud(ii)%frain*((P(i)-P(i-1))/real(nsubr))/(lambdaC*P(i)))
				Xc=max(Xc,0d0)
				Xc1=Xc
			enddo
			cloud_dens(i,ii)=Xc*dens(i)
			write(90,*) P(i),Xc,cloud_dens(i,ii)
		enddo
	endif

	tot=0d0
	do i=1,Cloud(ii)%nsize
		Cloud(ii)%w(i)=(1q4*Cloud(ii)%rv(i))**(1q0+(1q0-3q0*Cloud(ii)%veff)/Cloud(ii)%veff)*
     &					exp(-1q4*Cloud(ii)%rv(i)/(Cloud(ii)%reff*Cloud(ii)%veff))
		Cloud(ii)%w(i)=Cloud(ii)%w(i)*Cloud(ii)%rv(i)**3
		tot=tot+Cloud(ii)%w(i)
	enddo
	Cloud(ii)%w=Cloud(ii)%w/tot

	return
	end	


	real*8 function findpressure(P0)
	use GlobalSetup
	IMPLICIT NONE
	real*8 P0,r1,r2
	integer i

	do i=1,nr-1
		if(P(i).gt.P0.and.P(i+1).lt.P0) then
			r1=log10(sqrt(R(i)*R(i+1)))
			r2=log10(sqrt(R(i+1)*R(i+2)))
			findpressure=10d0**(r1+(r2-r1)*(P(i)-P0)/(P(i)-P(i+1)))
			return
		endif
	enddo
	call output("Error in finding pressure...")
	findpressure=R(nr+1)
	
	return
	end



	