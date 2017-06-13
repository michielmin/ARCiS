	subroutine SetupStructure(compute_mixrat)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 dp,dz,dlogp,RgasBar,Mtot,Pb(nr+1),tot,met_r,dens1bar,minZ,Tc
	real*8 Otot,Ctot,Htot,vescape,vtherm,RHill
	parameter(RgasBar=82.05736*1.01325)
	integer i,imol,nmix,j,niter,k
	logical ini,compute_mixrat
	character*500 cloudspecies(max(nclouds,1))
	
	real*8 starttime,stoptime,chemtime
	chemtime=0d0
	
	do i=1,nclouds
		cloudspecies(i)=Cloud(i)%species
	enddo

	ini = .TRUE.

	if(PTchemAbun) then
		call easy_chem_set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
		call call_easy_chem(Tchem,Pchem,mixrat_r(1,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &					XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),.false.)
		do i=2,nr
			mixrat_r(i,1:nmol)=mixrat_r(1,1:nmol)
		enddo
		mixrat(1:nmol)=mixrat_r(1,1:nmol)		
   	endif

	minZ=-5d0

	niter=1
	if(dochemistry) niter=2
	if(par_tprofile) niter=2

	Pb(1)=P(1)
	do i=2,nr
		Pb(i)=sqrt(P(i-1)*P(i))
	enddo
	Pb(nr+1)=P(nr)

	if(compute_mixrat) nabla_ad=2d0/7d0
	grav=Ggrav*Mplanet/(Rplanet)**2
	if(par_tprofile) call ComputeParamT(T)

	do j=1,niter

	Mtot=Mplanet
	grav=Ggrav*Mtot/(Rplanet)**2

	if(.not.mixratfile.and.(.not.dochemistry.or.j.eq.1)) then
	do i=1,nr
		tot=0d0
		do imol=1,nmol
			if(mixrat_r(i,imol).gt.0d0) tot=tot+mixrat_r(i,imol)
		enddo
		mixrat_r(i,1:nmol)=mixrat_r(i,1:nmol)/tot
	enddo
	endif

	call output("log(g) [cgs]: " // dbl2string(log10(grav(1)),'(f8.3)'))

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

		grav(i)=Ggrav*Mtot/R(i)**2

		if(i.lt.nr) then
			if(P(i).ge.1d0.and.P(i+1).lt.1d0.or.i.eq.1) then
				dens1bar=Avogadro*mp*mu/(RgasBar*T(i))
			endif
		endif
		
		dp=Pb(i)-Pb(i+1)
		dlogp=log(Pb(i)/Pb(i+1))

		Ndens(i)=P(i)*Avogadro/(RgasBar*T(i))
		Hp(i)=(T(i)*kb)/(grav(i)*mp*mu)
		dz=dlogp*Hp(i)
		dens(i)=Ndens(i)*mp*mu
		R(i+1)=R(i)+dz
		Mtot=Mtot+dens(i)*(R(i+1)**3-R(i)**3)*4d0*pi/3d0
	enddo
	Mtot=Mplanet
	do i=1,nr
		RHill=(Dplanet*(Mtot/(3d0*Mstar))**(1d0/3d0))
		if(R(i+1).gt.RHill) then
			print*,'layer',P(i),'is beyond the Hill Sphere'
			print*,'adjusting radius'
			R(i+1)=sqrt(R(i)*RHill)
		endif
		Mtot=Mtot+dens(i)*(R(i+1)**3-R(i)**3)*4d0*pi/3d0
	enddo
	do i=nr,1,-1
		vescape=sqrt(2d0*Ggrav*Mplanet/R(i))
		vtherm=sqrt(3d0*kb*T(i)/(mp*mu))
		Mtot=Mtot-dens(i)*(R(i+1)**3-R(i)**3)*4d0*pi/3d0
		if(vtherm.gt.vescape) then
			Ndens(i)=1d-20
			dens(i)=Ndens(i)*mp*mu
			print*,'layer',P(i),'escapes to space'
c			modelsucces=.false.
c			if(domakeai) return
		else
			exit
		endif
	enddo

	if(par_tprofile) call ComputeParamT(T)
	do i=1,nr
		if(T(i).gt.maxTprofile) T(i)=maxTprofile
		if(T(i).lt.3d0) T(i)=3d0
	enddo

	if(useDRIFT.and.domakeai) then
		modelsucces=.false.
		do i=1,nr
			if(T(i).lt.Tmax.and.Tmax.gt.0d0) modelsucces=.true.
		enddo
		do i=1,nr
			if(T(i).lt.Tmin.and.Tmin.gt.0d0) modelsucces=.false.
		enddo
		if(.not.modelsucces) return
	endif


	if(dochemistry.and.j.ne.niter) then
	if(compute_mixrat) then
		call easy_chem_set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
		if(Tform.gt.0d0) then
			call call_easy_chem(Tform,Pform,mixrat_r(1,1:nmol),molname(1:nmol),nmol,ini,.true.,cloudspecies,
     &					XeqCloud(1,1:nclouds),nclouds,nabla_ad(1),.true.)
    		ini=.true.
    	endif
		call output("==================================================================")
		call output("Computing chemistry using easy_chem by Paul Molliere")
		do i=1,nr
			call tellertje(i,nr)
			if(chemtime.gt.maxchemtime.and.domakeai) then
				modelsucces=.false.
				return
			endif
			met_r=metallicity
			if(sinkZ) then
				met_r=metallicity+log10((dens(i)/dens1bar)**(1d0/alphaZ**2-1d0))
				call easy_chem_set_molfracs_atoms(COratio,met_r,TiScale,enhancecarbon)
			endif
			if(P(i).ge.mixP.or.i.eq.1) then
				if(met_r.gt.minZ) then
					Tc=max(min(T(i),3000d0),100d0)
c					Tc=T(i)
					call cpu_time(starttime)
					call call_easy_chem(Tc,P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,.false.,cloudspecies,
     &						XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),.false.)
					call cpu_time(stoptime)
					chemtime=chemtime+stoptime-starttime
				else if(i.gt.1) then
					mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
					nabla_ad(i)=nabla_ad(i-1)
				else
					mixrat_r(i,1:nmol)=0d0
					mixrat_r(i,45)=0.85453462
					mixrat_r(i,48)=1d0-mixrat_r(i,45)
					nabla_ad(i)=2d0/7d0
				endif
			else
				mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
				nabla_ad(i)=nabla_ad(i-1)
			endif
			do imol=1,nmol
				if(includemol(imol)) then
					if(IsNaN(mixrat_r(i,imol))) then
						if(i.gt.1) then
							mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
						else
							mixrat_r(i,1:nmol)=0d0
							mixrat_r(i,45)=0.85453462
							mixrat_r(i,48)=1d0-mixrat_r(i,45)
						endif
					endif
				endif
			enddo
		enddo
		call output("==================================================================")
	else
		mixrat_r=mixrat_old_r
		XeqCloud=XeqCloud_old
	endif
	endif

	enddo

	if(dochemistry) then
	if(compute_mixrat) then
		ini=.true.
		call easy_chem_set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
		if(Tform.gt.0d0) then
			call call_easy_chem(Tform,Pform,mixrat_r(1,1:nmol),molname(1:nmol),nmol,ini,.true.,cloudspecies,
     &					XeqCloud(1,1:nclouds),nclouds,nabla_ad(1),.true.)
    		ini=.true.
    	endif
		call output("==================================================================")
		call output("Computing chemistry using easy_chem by Paul Molliere")
		do i=1,nr
			call tellertje(i,nr)
			if(chemtime.gt.maxchemtime.and.domakeai) then
				modelsucces=.false.
				return
			endif
			met_r=metallicity
			if(sinkZ) then
				met_r=metallicity+log10(dens(i)/dens1bar)*(1d0/alphaZ**2-1d0)
				call easy_chem_set_molfracs_atoms(COratio,met_r,TiScale,enhancecarbon)
			endif
			if(P(i).ge.mixP.or.i.eq.1) then
				if(met_r.gt.minZ) then
					Tc=max(min(T(i),3000d0),100d0)
					Tc=T(i)
				call cpu_time(starttime)
					call call_easy_chem(Tc,P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &					XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),.false.)
				call cpu_time(stoptime)
				chemtime=chemtime+stoptime-starttime
				else if(i.gt.1) then
					mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
					nabla_ad(i)=nabla_ad(i-1)
				else
					mixrat_r(i,1:nmol)=0d0
					mixrat_r(i,45)=0.85453462
					mixrat_r(i,48)=1d0-mixrat_r(i,45)
					nabla_ad(i)=2d0/7d0
				endif
			else
				Tc=max(min(T(i),3000d0),100d0)
					Tc=T(i)
				call cpu_time(starttime)
				if(cloudcompute) call call_easy_chem(Tc,P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &					XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),.false.)
				call cpu_time(stoptime)
				chemtime=chemtime+stoptime-starttime
				mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
				nabla_ad(i)=nabla_ad(i-1)
			endif
			do imol=1,nmol
				if(includemol(imol)) then
					if(IsNaN(mixrat_r(i,imol))) then
						print*,imol," is a NaN...",met_r
						if(i.gt.1) then
							mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
							nabla_ad(i)=nabla_ad(i-1)
						else
							mixrat_r(i,1:nmol)=0d0
							mixrat_r(i,45)=0.85453462
							mixrat_r(i,48)=1d0-mixrat_r(i,45)
							nabla_ad(i)=2d0/7d0
						endif
					endif
				endif
			enddo
		enddo
		call output("==================================================================")
	else
		mixrat_r=mixrat_old_r
		XeqCloud=XeqCloud_old
	endif
	endif

	do i=1,nclouds
		call SetupCloud(i)
	enddo

	Otot=0d0
	Ctot=0d0
	Htot=0d0
	do i=1,nr
		do imol=1,nmol
			if(includemol(imol)) then
				Otot=Otot+Ndens(i)*mixrat_r(i,imol)*real(Oatoms(imol))
				Ctot=Ctot+Ndens(i)*mixrat_r(i,imol)*real(Catoms(imol))
				Htot=Htot+Ndens(i)*mixrat_r(i,imol)*real(Hatoms(imol))
			endif
		enddo
	enddo
	COret=Ctot/Otot
	call output("C/O: " // dbl2string(COret,'(f8.3)'))
	call output("[O]: " // dbl2string(log10(Otot/Htot)-log10(0.0004509658/0.9207539305),'(f8.3)'))
	call output("[C]: " // dbl2string(log10(Ctot/Htot)-log10(0.0002478241/0.9207539305),'(f8.3)'))

	if(.not.PTchemAbun.and..not.dochemistry) then
		COratio=COret
		metallicity=log10((Ctot+Otot)/Htot)-log10((0.0002478241+0.0004509658)/0.9207539305)
	endif

	open(unit=50,file=trim(outputdir) // 'COprofile.dat',RECL=100)
	write(50,'("#",a14,3a10)') "P [bar]","C/O","[O]","[C]"
	do i=1,nr
		Otot=0d0
		Ctot=0d0
		Htot=0d0
		do imol=1,nmol
			if(includemol(imol)) then
				Otot=Otot+Ndens(i)*mixrat_r(i,imol)*real(Oatoms(imol))
				Ctot=Ctot+Ndens(i)*mixrat_r(i,imol)*real(Catoms(imol))
				Htot=Htot+Ndens(i)*mixrat_r(i,imol)*real(Hatoms(imol))
			endif
		enddo
		write(50,'(es15.4,3f10.3)') P(i),Ctot/Otot,
     &			log10(Otot/Htot)-log10(0.0004509658/0.9207539305),
     &			log10(Ctot/Htot)-log10(0.0002478241/0.9207539305)
	enddo
	close(unit=50)


	call output("Chemistry runtime:  " // trim(dbl2string((chemtime),'(f10.2)')) // " s")
	
	return
	end

	subroutine ComputeParamT(x)
	use GlobalSetup
	IMPLICIT NONE
	real*8 x(nr),tau,Tirr,eta,expint,dlnT,dlnP
	integer i

	tau=0d0
	Tirr=betaT*sqrt(Rstar/(2d0*Dplanet))*Tstar
	do i=nr,1,-1
		tau=kappaT*1d6*P(i)/grav(i)
		if(tau.lt.0d0) tau=0d0
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
		if(IsNaN(x(i))) then
			call output("NaN in temperature structure")
			if(i.gt.1) then
				x(i)=x(i-1)
			else
				x(i)=TeffP
			endif
		endif
		if(i.lt.nr) then
			dlnP=log(P(i+1)/P(i))
			dlnT=log(x(i+1)/x(i))
			if((dlnT/dlnP).gt.(nabla_ad(i)).and.adiabatic_tprofile) then
				dlnT=(nabla_ad(i))*dlnP
				x(i)=x(i+1)/exp(dlnT)
			endif
		endif
	enddo

	return
	end

	subroutine WriteStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 mix(nr,nmol)
	integer i,imol,nmix,ii
	character*1000 form
	character*10 namemix(nmol)
	
	open(unit=50,file=trim(outputdir) // 'densityprofile.dat',RECL=100)
	write(50,'("#",a14,a15,a15,a13,a10,a10)') "radius [cm]","height [cm]","dens [g/cm^3]","N [1/cm^3]","T [K]","P [bar]"
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
	write(50,form) "T [K]","P [bar]",namemix(1:nmix)
	form='(f10.3,es10.3,' // trim(int2string(nmix,'(i2)')) // 'es15.4)'
	do i=1,nr
		write(50,form) T(i),P(i),mix(i,1:nmix)
	enddo
	close(unit=50)

	do ii=1,nclouds
		open(unit=50,file=trim(outputdir) // 'clouddens' // trim(int2string(ii,'(i0.2)')) // '.dat',RECL=1000)
		form='("#",a9,a15,a15,a15)'
		write(50,form) "P [bar]","dens [g/cm^3]","Eq. dens","gas dens"
		form='(es10.3,es15.3,es15.3,es15.3)'
		do i=1,nr
			write(50,form) P(i),cloud_dens(i,ii),dens(i)*XeqCloud(i,ii),dens(i)
		enddo
		close(unit=50)
	enddo

	return
	end
	

	subroutine RunDRIFT(ii)
	use GlobalSetup
	use Constants
	use AtomsModule
	IMPLICIT NONE
	integer i,ii,j
	real*8 tmix,frac(nr,10),temp(nr,3)
	real*8 elabun(nr,7)
	character*500 command
	character*500 filename
	logical ini
	character*500 cloudspecies(max(nclouds,1))
	
	call output("Running DRIFT cloud formation model")

	open(unit=25,file=trim(outputdir) // 'SPARCtoDRIFT.dat',RECL=1000)
	write(25,'("#Elemental abundances")')
	call easy_chem_set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
	do i=1,18
		write(25,'(se18.6,"   ",a5)') molfracs_atoms(i),names_atoms(i)
	enddo

c	modelsucces=.false.
c	do i=1,nr
c		if(T(i).lt.3000d0) modelsucces=.true.
c	enddo
c	if(.not.modelsucces) then
c		close(unit=25)
c		return
c	endif
	modelsucces=.true.

	write(25,'("#Density setup")')
	write(25,'(i5)') nr
	do i=nr,1,-1
		tmix=Cloud(ii)%tmix*P(i)**(-Cloud(ii)%betamix)
		if(T(i).lt.Tmin.and.domakeai) then
			modelsucces=.false.
			close(unit=25)
			return
		endif
   	 	write(25,*) T(i), P(i) , R(i), dens(i), grav(i), 1d0/tmix
	enddo
	close(unit=25)
	command="rm -rf " // trim(outputdir) // "restart.dat"
	call system(command)
	command="cd " // trim(outputdir) // "; gtimeout 900s nohup static_weather12 4 1d-3"
	call system(command)

	inquire(file=trim(outputdir) // "done",exist=modelsucces)
	if(modelsucces) then
		open(unit=25,file=trim(outputdir) // "done")
		read(25,*) modelsucces
	endif
	command="rm -rf " // trim(outputdir) // "nohup.out"
	call system(command)

	if(.not.modelsucces.and.domakeai) return

	filename=trim(outputdir) // "out3_dust.dat"
	call regridN(filename,P*1d6,frac,nr,2,9,10,4,.true.,.true.)
	Cloud(ii)%frac(1:nr,1)=frac(1:nr,1)/3d0
	Cloud(ii)%frac(1:nr,2)=frac(1:nr,1)/3d0
	Cloud(ii)%frac(1:nr,3)=frac(1:nr,1)/3d0

	Cloud(ii)%frac(1:nr,4)=frac(1:nr,2)/3d0
	Cloud(ii)%frac(1:nr,5)=frac(1:nr,2)/3d0
	Cloud(ii)%frac(1:nr,6)=frac(1:nr,2)/3d0

	Cloud(ii)%frac(1:nr,7)=frac(1:nr,3)
	Cloud(ii)%frac(1:nr,8)=frac(1:nr,4)
	Cloud(ii)%frac(1:nr,9)=frac(1:nr,5)
	Cloud(ii)%frac(1:nr,10)=frac(1:nr,6)
	Cloud(ii)%frac(1:nr,11)=frac(1:nr,7)
	Cloud(ii)%frac(1:nr,12)=frac(1:nr,8)

	Cloud(ii)%frac(1:nr,13)=frac(1:nr,9)/3d0
	Cloud(ii)%frac(1:nr,14)=frac(1:nr,9)/3d0
	Cloud(ii)%frac(1:nr,15)=frac(1:nr,9)/3d0

	Cloud(ii)%frac(1:nr,16)=frac(1:nr,10)

	call regridN(filename,P*1d6,cloud_dens(1:nr,ii),nr,2,6,1,4,.true.,.false.)
	call regridN(filename,P*1d6,elabun,nr,2,36,7,4,.true.,.true.)

	do i=1,nclouds
		cloudspecies(i)=Cloud(i)%species
	enddo
	if(mixratfile) then
		call regridN(TPfile,P,mixrat_r,nr,2,3,nmol,2,.true.,.false.)
	else
		call output("Computing chemistry after dust condensation")
		do i=1,nr
			call tellertje(i,nr)
			ini=.true.
			if(cloud_dens(i,ii).lt.1d-25) then
				call easy_chem_set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
			else
				call easy_chem_set_molfracs_atoms_elabun(COratio,metallicity,elabun(i,1:7))
			endif
			call call_easy_chem(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,.false.,cloudspecies,
     &						XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),.false.)
		enddo
	endif

	filename=trim(outputdir) // "out3_dist.dat"
	call regridN(filename,P,temp,nr,1,4,3,4,.true.,.false.)

	cloud_dens(1:nr,ii)=temp(1:nr,1)
	if(.not.allocated(Cloud(ii)%rv)) then
		allocate(Cloud(ii)%rv(nr))
		allocate(Cloud(ii)%sigma(nr))
	endif
	Cloud(ii)%rv(1:nr)=temp(1:nr,2)*1d4
	Cloud(ii)%sigma(1:nr)=temp(1:nr,3)*1d4

	call output("Computing inhomogeneous cloud particles")

	call SetupPartCloud(ii)

	cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*Cloud(ii)%M(1:nr)

	open(unit=25,file=trim(outputdir) // "DRIFTcloud" // trim(int2string(ii,'(i0.4)')),recl=6000)
	do i=1,nr
		write(25,*) P(i),T(i),Cloud(ii)%rv(i),Cloud(ii)%sigma(i),cloud_dens(i,ii),dens(i),
     &			Cloud(ii)%frac(i,1)*3d0,Cloud(ii)%frac(i,4)*3d0,Cloud(ii)%frac(i,7:12),
     &			Cloud(ii)%frac(i,13)*3d0,Cloud(ii)%frac(i,16)
	enddo
	close(unit=25)

	return
	end
	

	subroutine SetupCloud(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 pp,tot,column,tt,z
	integer ii,i,j,nsubr,isize,ilam
	real*8 Xc,Xc1,lambdaC,Ca,Cs,tau,P_SI1,P_SI2
	logical cl
	
	if(useDRIFT) then
		if(cloudcompute) then
			if(.not.allocated(Cloud(ii)%rv)) then
				allocate(Cloud(ii)%rv(nr))
				allocate(Cloud(ii)%sigma(nr))
			endif
			call DiffuseCloud(ii)
			Cloud(ii)%rv=Cloud(ii)%rv*1d4
			Cloud(ii)%sigma=1d-10
			Cloud(ii)%frac(1:nr,1:16)=1d-10
			Cloud(ii)%frac(1:nr,13:15)=1d0/3d0
			call output("Computing inhomogeneous cloud particles")

			call SetupPartCloud(ii)

			open(unit=25,file=trim(outputdir) // "DRIFTcloud" // trim(int2string(ii,'(i0.4)')),recl=6000)
			do i=1,nr
				write(25,*) P(i),T(i),Cloud(ii)%rv(i),Cloud(ii)%sigma(i),cloud_dens(i,ii),dens(i),
     &			Cloud(ii)%frac(i,1)*3d0,Cloud(ii)%frac(i,4)*3d0,Cloud(ii)%frac(i,7:12),
     &			Cloud(ii)%frac(i,13)*3d0,Cloud(ii)%frac(i,16)
			enddo
			close(unit=25)
		else if(Cloud(ii)%file.ne.' ') then
			call regridN(Cloud(ii)%file,P*1d6,cloud_dens(1:nr,ii),nr,2,6,1,3,.false.,.false.)
			cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)
			if(.not.allocated(Cloud(ii)%rv)) then
				allocate(Cloud(ii)%rv(nr))
				allocate(Cloud(ii)%sigma(nr))
			endif
			call regridN(Cloud(ii)%file,P*1d6,Cloud(ii)%rv(1:nr),nr,2,13,1,3,.false.,.true.)
			Cloud(ii)%rv=Cloud(ii)%rv*1d4
			Cloud(ii)%sigma=1d-10
			Cloud(ii)%frac(1:nr,1:16)=1d-10
			Cloud(ii)%frac(1:nr,13:15)=1d0/3d0
			call output("Computing inhomogeneous cloud particles")

			call SetupPartCloud(ii)

			open(unit=25,file=trim(outputdir) // "DRIFTcloud" // trim(int2string(ii,'(i0.4)')),recl=6000)
			do i=1,nr
				write(25,*) P(i),T(i),Cloud(ii)%rv(i),Cloud(ii)%sigma(i),cloud_dens(i,ii),dens(i),
     &			Cloud(ii)%frac(i,1)*3d0,Cloud(ii)%frac(i,4)*3d0,Cloud(ii)%frac(i,7:12),
     &			Cloud(ii)%frac(i,13)*3d0,Cloud(ii)%frac(i,16)
			enddo
			close(unit=25)
		else
			call RunDRIFT(ii)
		endif
		return
	endif

	if(.not.cloudcompute) then
		if(Cloud(ii)%haze) then
			column=0d0
			do i=1,nr
				cloud_dens(i,ii)=Cloud(ii)%mixrat*dens(i)
				column=column+cloud_dens(i,ii)*(R(i+1)-R(i))
			enddo
			Cloud(ii)%column=column
		else
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
		endif
	else
c use Ackerman & Marley 2001 cloud computation
		column=0d0
		nsubr=1
		if(Cloud(ii)%haze) then
			Xc=Cloud(ii)%fcond*XeqCloud(1,ii)
		else
			Xc=0d0
		endif
		Xc1=Xc
		cloud_dens(1,ii)=Xc*dens(1)
		cl=.false.
		do i=2,nr
			lambdaC=max(0.1d0,(log(T(i)/T(i-1))/log(P(i)/P(i-1)))/nabla_ad(i))
			do j=1,nsubr
				if(Cloud(ii)%haze) then
					Xc=Cloud(ii)%fcond*XeqCloud(i,ii)
				else
					Xc=(Xc1+(XeqCloud(i,ii)-XeqCloud(i-1,ii))/real(nsubr))/(1d0-Cloud(ii)%frain*((P(i)-P(i-1))/real(nsubr))/(lambdaC*P(i)))
				endif
				if(Xc.le.0d0.and.cl) then
					cloud_dens(i,ii)=0d0
					goto 2
				endif
				Xc=max(Xc,0d0)
				if(Xc.gt.0d0) cl=.true.
				Xc1=Xc
			enddo
			cloud_dens(i,ii)=Xc*dens(i)
		enddo
2 		continue
	endif

	j=0
1	tot=0d0
	do i=1,Cloud(ii)%nr
		Cloud(ii)%w(i)=(1q4*Cloud(ii)%rv(i))**(1q0+(1q0-3q0*Cloud(ii)%veff)/Cloud(ii)%veff)*
     &					exp(-1q4*Cloud(ii)%rv(i)/(Cloud(ii)%reff*Cloud(ii)%veff))
		Cloud(ii)%w(i)=Cloud(ii)%w(i)*Cloud(ii)%rv(i)**3
		tot=tot+Cloud(ii)%w(i)
	enddo
	if(.not.tot.gt.0d0) then
		if(j.lt.10) then
			Cloud(ii)%veff=Cloud(ii)%veff*2d0
			call output("increasing veff")
			j=j+1
			goto 1
		else
			tot=1d0
		endif
	endif
	Cloud(ii)%w=Cloud(ii)%w/tot
	
	if(Cloud(ii)%tau.gt.0d0) then
		tau=0d0
		do ilam=nlam-1,1,-1
			if(lam(ilam).gt.Cloud(ii)%lam.and.lam(ilam+1).le.Cloud(ii)%lam) exit
		enddo
		if(ilam.lt.1) ilam=1
		do i=nr,2,-1
			Ca=0d0
			Cs=0d0
			do isize=1,Cloud(ii)%nr
				Ca=Ca+
     &		Cloud(ii)%Kabs(isize,ilam)*Cloud(ii)%w(isize)*cloud_dens(i,ii)
				Cs=Cs+
     &		Cloud(ii)%Ksca(isize,ilam)*Cloud(ii)%w(isize)*cloud_dens(i,ii)
			enddo
			tau=tau+(R(i)-R(i-1))*(Ca+Cs)
		enddo
		cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*Cloud(ii)%tau/tau
	endif

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


	subroutine OnlyChemCompute()
	use GlobalSetup
	IMPLICIT NONE
	logical ini
	integer i
	character*500 cloudspecies(max(nclouds,1))
	
	do i=1,nclouds
		cloudspecies(i)=Cloud(i)%species
	enddo
	call easy_chem_set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
	ini=.true.
	do i=1,nr
		call call_easy_chem(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &					XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),.false.)
	enddo

	return
	end
	


	