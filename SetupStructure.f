	subroutine SetupStructure(compute_mixrat)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 dp,dz,dlogp,RgasBar,Mtot,Pb(nr+1),tot,met_r,dens1bar,minZ,Tc,Rscale
	real*8 Otot,Ctot,Htot,vescape,vtherm,RHill,MMW_form,P0,R0
	parameter(RgasBar=82.05736*1.01325)
	integer i,imol,nmix,j,niter,k,i1,i2,di,ii,i3
	logical ini,compute_mixrat
	character*500 cloudspecies(max(nclouds,1))
	
	real*8 starttime,stoptime,chemtime
	chemtime=0d0
		
	do i=1,nclouds
		cloudspecies(i)=Cloud(i)%species
	enddo

	ini = .TRUE.

	if(Tform.gt.0d0) then
		call FormAbun(Tform,f_dry,f_wet,COratio,metallicity0,metallicity)
	endif

	if(PTchemAbun) then
		call set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
		call call_chemistry(Tchem,Pchem,mixrat_r(1,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW_form,didcondens_chem,includemol)
		do i=2,nr
			mixrat_r(i,1:nmol)=mixrat_r(1,1:nmol)
		enddo
		mixrat(1:nmol)=mixrat_r(1,1:nmol)		
   	endif

	minZ=-5d0

	niter=3
	if(dochemistry) niter=3
	if(par_tprofile) niter=3

	Pb(1)=P(1)
	do i=2,nr
		Pb(i)=sqrt(P(i-1)*P(i))
	enddo
	Pb(nr+1)=P(nr)

	if(compute_mixrat) nabla_ad=2d0/7d0
	grav=Ggrav*Mplanet/(Rplanet)**2
	if(par_tprofile) call ComputeParamT(T)

	Rscale=1d0

	grav=Ggrav*Mplanet/Rplanet**2

	do j=1,niter

	if(.not.mixratfile.and.(.not.dochemistry.or.j.eq.1)) then
	do i=1,nr
		tot=0d0
		do imol=1,nmol
			if(mixrat_r(i,imol).gt.0d0) tot=tot+mixrat_r(i,imol)
		enddo
		mixrat_r(i,1:nmol)=mixrat_r(i,1:nmol)/tot
	enddo
	endif

	call output("log(g) [cgs]: " // dbl2string(log10(Ggrav*Mplanet/Rplanet**2),'(f8.3)'))

	Mtot=Mplanet

	if(j.eq.1) then
		MMW=0d0
		do imol=1,nmol
			if(mixrat_r(1,imol).gt.0d0) MMW=MMW+mixrat_r(1,imol)*Mmol(imol)
		enddo
	endif
	call output("Mean molecular weight: " // dbl2string(MMW(1),'(f8.3)'))
	i1=1
	do i=1,nr
		if(Pb(i).ge.Pplanet.and.Pb(i+1).lt.Pplanet) then
			i1=i
		endif
	enddo
	i2=nr
	di=1
	do ii=1,2
	P0=Pplanet
	R0=Rplanet
	do i=i1,i2,di
		if(di.gt.0) then
			i3=i
		else
			i3=i-1
		endif

		if(.not.dochemistry) then
			MMW(i3)=0d0
			do imol=1,nmol
				if(mixrat_r(i3,imol).gt.0d0) MMW(i3)=MMW(i3)+mixrat_r(i3,imol)*Mmol(imol)
			enddo
		endif

		if(Pb(i3).ge.1d0.and.Pb(i3+1).lt.1d0) then
			dens1bar=Avogadro*mp*MMW(i3)/(RgasBar*T(i3))
		endif
		
		dp=P0-Pb(i+di)
		dlogp=log(P0/Pb(i+di))

		Ndens(i3)=P(i3)*Avogadro/(RgasBar*T(i3))
		Hp(i3)=(T(i3)*kb)/(grav(i3)*mp*MMW(i3))

		dz=dlogp*Hp(i3)
		dens(i3)=Ndens(i3)*mp*MMW(i3)

		R(i+di)=R0+dz
		if(R(i+di).lt.R(i)/5d0) R(i+di)=R(i)/5d0
		R0=R(i+di)
		P0=Pb(i+di)
	enddo
	i1=i1+1
	i2=2
	di=-1
	R0=Rplanet
	P0=Pplanet
	enddo
	Mtot=Mplanet
	do i=1,nr
		RHill=(Dplanet*(Mtot/(3d0*Mstar))**(1d0/3d0))
		if(R(i+1).gt.RHill) then
			call output("layer" // dbl2string(P(i),'(es10.3E3)') // "is beyond the Hill Sphere")
c			print*,'adjusting radius'
c			R(i+1)=sqrt(R(i)*RHill)
		endif
		Mtot=Mtot+dens(i)*(R(i+1)**3-R(i)**3)*4d0*pi/3d0
		grav(i)=Ggrav*Mtot/(R(i)*R(i+1))
	enddo
	do i=nr,1,-1
		vescape=sqrt(2d0*Ggrav*Mplanet/R(i))
		vtherm=sqrt(3d0*kb*T(i)/(mp*MMW(i)))
		Mtot=Mtot-dens(i)*(R(i+1)**3-R(i)**3)*4d0*pi/3d0
		if(vtherm.gt.vescape) then
			Ndens(i)=1d-20
			dens(i)=Ndens(i)*mp*MMW(i)
			call output("layer" // dbl2string(P(i),'(es10.3E3)') // "escapes to space")
			modelsucces=.false.
c			if(domakeai.or.retrieval) return
		else
			exit
		endif
	enddo

	if(par_tprofile) call ComputeParamT(T)
	do i=1,nr
		if(T(i).gt.maxTprofile) T(i)=maxTprofile
		if(T(i).lt.3d0) T(i)=3d0
	enddo

c	if(useDRIFT.and.domakeai) then
	if(domakeai) then
		modelsucces=.false.
		do i=1,nr
			if(T(i).lt.Tmax.and.Tmax.gt.0d0) modelsucces=.true.
		enddo
		do i=1,nr
			if(T(i).lt.Tmin.and.Tmin.gt.0d0) modelsucces=.false.
		enddo
		if(.not.modelsucces) return
	endif


	if(dochemistry.and.j.eq.1) then
	if(compute_mixrat) then
		if(Tform.gt.0d0) then
			call FormAbun(Tform,f_dry,f_wet,COratio,metallicity0,metallicity)
		else
			call set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
    	endif
		call output("==================================================================")
c		call output("Computing chemistry using easy_chem by Paul Molliere")
		call output("Computing chemistry using GGchem by Peter Woitke")
		do i=1,nr
			call tellertje(i,nr)
			if(chemtime.gt.maxchemtime.and.domakeai) then
				modelsucces=.false.
				return
			endif
			met_r=metallicity
			if(sinkZ) then
				met_r=metallicity+log10((dens(i)/dens1bar)**(1d0/alphaZ**2-1d0))
				call set_molfracs_atoms(COratio,met_r,TiScale,enhancecarbon)
			endif
			if(P(i).ge.mixP.or.i.eq.1) then
				if(met_r.gt.minZ) then
					Tc=max(min(T(i),20000d0),100d0)
c					Tc=T(i)
					call cpu_time(starttime)
					call call_chemistry(Tc,P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &					XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol)
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

	if(dochemistry.and..false.) then
	if(compute_mixrat) then
		ini=.true.
		if(Tform.gt.0d0) then
			call FormAbun(Tform,f_dry,f_wet,COratio,metallicity0,metallicity)
		else
			call set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
    	endif
		call output("==================================================================")
c		call output("Computing chemistry using easy_chem by Paul Molliere")
		call output("Computing chemistry using GGchem by Peter Woitke")
		do i=1,nr
			call tellertje(i,nr)
			if(chemtime.gt.maxchemtime.and.domakeai) then
				modelsucces=.false.
				return
			endif
			met_r=metallicity
			if(sinkZ) then
				met_r=metallicity+log10(dens(i)/dens1bar)*(1d0/alphaZ**2-1d0)
				call set_molfracs_atoms(COratio,met_r,TiScale,enhancecarbon)
			endif
			if(P(i).ge.mixP.or.i.eq.1) then
				if(met_r.gt.minZ) then
					Tc=max(min(T(i),3000d0),100d0)
					Tc=T(i)
				call cpu_time(starttime)
					call call_chemistry(Tc,P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol)
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
				if(cloudcompute) call call_chemistry(Tc,P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol)
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

	if(compute_mixrat) then
		do i=1,nclouds
			call SetupCloud(i)
		enddo
	endif

c bug needs to be fixed!!!!!
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
c bug needs to be fixed!!!!!
	COret=Ctot/Otot
	if(Tform.gt.0d0) COret=COratio
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


	open(unit=20,file='output_struct.dat',RECL=1000)
	do i=1,nr
		write(20,*) P(i),MMW(i),grav(i),sqrt(R(i)*R(i+1))/Rjup,T(i),Hp(i)
	enddo
	close(unit=20)
	
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
	character*10 side
	
	if(i2d.eq.0) then
		side=" "
	else
		write(side,'("_",i0.2)') i2d
	endif
	
	open(unit=50,file=trim(outputdir) // 'densityprofile' // trim(side) // '.dat',RECL=100)
	write(50,'("#",a14,a15,a15,a13,a10,a11)') "radius [cm]","height [cm]","dens [g/cm^3]","N [1/cm^3]","T [K]","P [bar]"
	do i=1,nr
		write(50,'(es15.7E3,es15.4E3,es15.4E3,es13.4E3,f10.3,es11.3E3)') sqrt(R(i)*R(i+1)),sqrt(R(i)*R(i+1))-Rplanet
     &			,dens(i),Ndens(i),T(i),P(i)
	enddo
	close(unit=50)

	open(unit=50,file=trim(outputdir) // 'density' // trim(side) // '.dat',RECL=1000)
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

	open(unit=50,file=trim(outputdir) // 'mixingratios' // trim(side) // '.dat',RECL=1000)
	form='("#",a9,a13,' // trim(int2string(nmix,'(i2)')) // 'a15)'
	write(50,form) "T [K]","P [bar]",namemix(1:nmix)
	form='(f10.3,es13.3E3,' // trim(int2string(nmix,'(i2)')) // 'es15.4E3)'
	do i=1,nr
		write(50,form) T(i),P(i),mix(i,1:nmix)
	enddo
	close(unit=50)

	do ii=1,nclouds
		open(unit=50,file=trim(outputdir) // 'clouddens' // trim(int2string(ii,'(i0.2)')) // trim(side) // '.dat',RECL=1000)
		form='("#",a12,a15,a15,a15)'
		write(50,form) "P [bar]","dens [g/cm^3]","Eq. dens","gas dens"
		form='(es13.3,es15.3E3,es15.3E3,es15.3E3)'
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
	real*8 elabun(nr,7),MMW_form,Kzz(nr),Nmfp
	character*500 command
	character*500 filename
	logical ini
	character*500 cloudspecies(max(nclouds,1))
	
	call output("Running DRIFT cloud formation model")

	open(unit=25,file=trim(outputdir) // 'SPARCtoDRIFT_CO11.dat',RECL=1000)
	write(25,'("#Elemental abundances")')
  	ini=.true.
	if(Tform.gt.0d0) then
		call FormAbun(Tform,f_dry,f_wet,COratio,metallicity0,metallicity)
	else
		call set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
   	endif
	do i=1,18
		write(25,'(se18.6E3,"   ",a5)') molfracs_atoms(i),names_atoms(i)
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
	if(Cloud(ii)%Kzzfile.ne.' ') then
		call regridN(Cloud(ii)%Kzzfile,P,Kzz,nr,1,2,1,0,.true.,.true.)
		Kzz=Kzz*Cloud(ii)%Kscale
c		do i=nr,1,-1
c			if(didcondens(i)) then
c				do j=i,1,-1
c					if(.not.didcondens(j)) exit
c				enddo
c				if(j.lt.1) j=1
c				Nmfp=log(P(j)/P(i))
c				j=j+real(i-j)*0.7d0
c				tmix=Nmfp**2*Hp(j)**2/Kzz(j)
c				Kzz(i)=Hp(i)**2/tmix
c			endif
c		enddo
		do i=nr,1,-1
c			write(82,*) P(i), Kzz(i)
			tmix=Hp(i)**2/abs(Kzz(i))
			if(T(i).lt.Tmin.and.domakeai) then
				modelsucces=.false.
				close(unit=25)
				return
			endif
   		 	write(25,*) T(i), P(i) , R(i), dens(i), grav(i), 1d0/tmix
		enddo
	else
		do i=nr,1,-1
			tmix=Cloud(ii)%tmix*P(i)**(-Cloud(ii)%betamix)
			tmix=Hp(i)**2/abs(Cloud(ii)%Kzz)
			if(didcondens(i)) then
				do j=i,1,-1
					if(.not.didcondens(j)) exit
				enddo
				if(j.lt.1) j=1
				Nmfp=log(P(j)/P(i))
				j=j+real(i-j)*0.7d0
				tmix=tmix+Nmfp**2*Hp(j)**2/Cloud(ii)%Kzz
			endif
			if(T(i).lt.Tmin.and.domakeai) then
				modelsucces=.false.
				close(unit=25)
				return
			endif
   		 	write(25,*) T(i), P(i) , R(i), dens(i), grav(i), 1d0/tmix
		enddo
	endif
	close(unit=25)
	command="rm -rf " // trim(outputdir) // "restart.dat"
	call system(command)
	command="cd " // trim(outputdir) // "; gtimeout 3600s nohup static_weather13 4 1d-3"
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

	Cloud(ii)%frac(1:nr,17)=0d0
	Cloud(ii)%frac(1:nr,18)=0d0

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
				if(Tform.gt.0d0) then
					call FormAbun(Tform,f_dry,f_wet,COratio,metallicity0,metallicity)
				else
					call set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
		    	endif
			else
				call set_molfracs_atoms_elabun(COratio,metallicity,elabun(i,1:7))
			endif
			call call_chemistry(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,.false.,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol)
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
			call regridN(Cloud(ii)%file,P*1d6,Cloud(ii)%rv(1:nr),nr,2,11,1,3,.false.,.true.)
			Cloud(ii)%rv=Cloud(ii)%rv*1d4
			Cloud(ii)%sigma=1d-10
			Cloud(ii)%frac(1:nr,1:18)=1d-10
c 10% iron
			Cloud(ii)%frac(1:nr,9)=0.1d0
c 90% MgSiO3
			Cloud(ii)%frac(1:nr,13:15)=0.9d0/3d0
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
	call set_molfracs_atoms(COratio,metallicity,TiScale,enhancecarbon)
	ini=.true.
	do i=1,nr
		call call_chemistry(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol)
	enddo

	return
	end
	

	subroutine FormAbun(T,f_dry,f_wet,COratio,metallicity0,metallicity)
	use AtomsModule
	use OutputModule
	IMPLICIT NONE
	character*10 species(20)
	real*8 sil_atoms(N_atoms),tot_atoms(N_atoms),COratio,metallicity
	real*8 Tmax(20),abun(20),max,maxf(20),T,f_dry,f_wet,metallicity0,tot,atoms(20,N_atoms)
	real*8 dry_atoms(N_atoms),wet_atoms(N_atoms)
	integer limit(20,N_atoms),nspecies,i,j,k
	character*100 temp,filename
	logical molecule(20)

	names_atoms(1) = 'H'
	names_atoms(2) = 'He'
	names_atoms(3) = 'C'
	names_atoms(4) = 'N'
	names_atoms(5) = 'O'
	names_atoms(6) = 'Na'
	names_atoms(7) = 'Mg'
	names_atoms(8) = 'Al'
	names_atoms(9) = 'Si'
	names_atoms(10) = 'P'
	names_atoms(11) = 'S'
	names_atoms(12) = 'Cl'
	names_atoms(13) = 'K'
	names_atoms(14) = 'Ca'
	names_atoms(15) = 'Ti'
	names_atoms(16) = 'V'
	names_atoms(17) = 'Fe'
	names_atoms(18) = 'Ni'

	molfracs_atoms = (/ 0.9207539305,
     &  0.0783688694,
     &  0.0002478241, 
     &  6.22506056949881e-05, 
     &  0.0004509658, 
     &  1.60008694353205e-06, 
     &  3.66558742055362e-05, 
     &  2.595e-06, 
     &  2.9795e-05, 
     &  2.36670201997668e-07, 
     &  1.2137900734604e-05, 
     &  2.91167958499589e-07, 
     &  9.86605611925677e-08, 
     &  2.01439011429255e-06, 
     &  8.20622804366359e-08,
     &  7.83688694089992e-09,
     &  2.91167958499589e-05,
     &  1.52807116806281e-06
     &  /)

	atoms=0
	maxf=1d0
	molecule=.false.
	i=1
	species(i)='C'
		atoms(i,3)=1
		Tmax(i)=1000d0
		maxf(i)=0.5d0
	i=i+1
	species(i)='CO'
		atoms(i,3)=1
		atoms(i,5)=1
		Tmax(i)=20d0
		molecule(i)=.true.
	i=i+1
	species(i)='TiO2'
		atoms(i,5)=2
		atoms(i,15)=1
		Tmax(i)=1700d0
	i=i+1
	species(i)='VO'
		atoms(i,5)=1
		atoms(i,16)=1
		Tmax(i)=1700d0
	i=i+1
	species(i)='Silicates'
		atoms(i,9)=1
		atoms(i,6)=min(molfracs_atoms(6),molfracs_atoms(8))/molfracs_atoms(9)
		atoms(i,8)=min(molfracs_atoms(6),molfracs_atoms(8))/molfracs_atoms(9)
		atoms(i,7)=molfracs_atoms(7)/molfracs_atoms(9)
		atoms(i,13)=molfracs_atoms(13)/molfracs_atoms(9)
		atoms(i,14)=molfracs_atoms(14)/molfracs_atoms(9)
		atoms(i,5)=atoms(i,6)+atoms(i,7)+atoms(i,8)+atoms(i,14)+atoms(i,13)+2d0
		Tmax(i)=1700d0
c		write(*,'("Al",f3.1,"Na",f3.1,"Mg",f3.1,"SiO",f3.1)') atoms(i,8),atoms(i,6),atoms(i,7),atoms(i,5)
	i=i+1
	species(i)='SiO2'
		atoms(i,9)=1
		atoms(i,5)=2
		Tmax(i)=1700d0
	i=i+1
	species(i)='FeS'
		atoms(i,11)=1
		atoms(i,17)=1
		Tmax(i)=700d0
	i=i+1
	species(i)='Fe'
		atoms(i,17)=1
		Tmax(i)=1700d0
	i=i+1
	species(i)='SiC'
		atoms(i,9)=1
		atoms(i,3)=1
		Tmax(i)=1700d0
	i=i+1
	species(i)='Al2O3'
		atoms(i,5)=3
		atoms(i,8)=2
		Tmax(i)=1700d0
	i=i+1
	species(i)='CO2'
		atoms(i,3)=1
		atoms(i,5)=2
		Tmax(i)=47d0
	i=i+1
	species(i)='H2O'
		atoms(i,1)=2
		atoms(i,5)=1
		Tmax(i)=135d0
		
	nspecies=i

	gas_atoms=molfracs_atoms
	tot_atoms=0d0
	solid_atoms=0d0
	abun=0d0
	do i=1,nspecies
		max=1d200
		do j=1,n_atoms
			if(atoms(i,j).gt.0) then
				if(maxf(i)*gas_atoms(j)/real(atoms(i,j)).lt.max) max=maxf(i)*gas_atoms(j)/real(atoms(i,j))
			endif
		enddo
		if(T.lt.Tmax(i)) then
			abun(i)=max
			do j=1,n_atoms
				gas_atoms(j)=gas_atoms(j)-abun(i)*real(atoms(i,j))
				solid_atoms(j)=solid_atoms(j)+abun(i)*real(atoms(i,j))
			enddo
		else if(molecule(i)) then
			abun(i)=max
			do j=1,n_atoms
				gas_atoms(j)=gas_atoms(j)-abun(i)*real(atoms(i,j))
				tot_atoms(j)=tot_atoms(j)+abun(i)*real(atoms(i,j))
			enddo
			abun(i)=0d0
		endif
	enddo
	gas_atoms=gas_atoms+tot_atoms

	dry_atoms(1:N_atoms)=solid_atoms(1:N_atoms)-abun(nspecies)*atoms(nspecies,1:N_atoms)
	wet_atoms=solid_atoms-dry_atoms

	tot_atoms=gas_atoms+f_dry*dry_atoms+f_wet*wet_atoms
	tot_atoms(1:2)=tot_atoms(1:2)*10d0**-metallicity0
	tot_atoms=tot_atoms/sum(tot_atoms(1:N_atoms))

	COratio=tot_atoms(3)/tot_atoms(5)
	call output("C/O: " // dbl2string(tot_atoms(3)/tot_atoms(5),'(f6.2)'))
	call output("Si/O:" // dbl2string(tot_atoms(9)/tot_atoms(5),'(f6.2)'))

	metallicity=log10((sum(tot_atoms(3:n_atoms))/sum(tot_atoms(1:n_atoms)))/
     &			(sum(molfracs_atoms(3:n_atoms))/sum(molfracs_atoms(1:n_atoms))))
	call output("[Z]: " // dbl2string(metallicity,'(f6.2)'))

	tot=sum(tot_atoms(1:N_atoms))
	molfracs_atoms=tot_atoms/tot

	open(unit=50,file='atomic.dat')
	do i=1,18
		write(50,'(a5,se18.6)') names_atoms(i),molfracs_atoms(i)
	enddo
	close(unit=50)

	
	return
	end
	

	subroutine set_molfracs_atoms(CO,Z,TiScale,enhancecarbon)
	use AtomsModule
	implicit none
	real*8 CO,Z,tot,TiScale
	integer i
	logical enhancecarbon

	names_atoms(1) = 'H'
	names_atoms(2) = 'He'
	names_atoms(3) = 'C'
	names_atoms(4) = 'N'
	names_atoms(5) = 'O'
	names_atoms(6) = 'Na'
	names_atoms(7) = 'Mg'
	names_atoms(8) = 'Al'
	names_atoms(9) = 'Si'
	names_atoms(10) = 'P'
	names_atoms(11) = 'S'
	names_atoms(12) = 'Cl'
	names_atoms(13) = 'K'
	names_atoms(14) = 'Ca'
	names_atoms(15) = 'Ti'
	names_atoms(16) = 'V'
	names_atoms(17) = 'Fe'
	names_atoms(18) = 'Ni'

	molfracs_atoms = (/ 0.9207539305,
     &	0.0783688694,
     &	0.0002478241,
     &	6.22506056949881e-05,
     &	0.0004509658,
     &	1.60008694353205e-06,
     &	3.66558742055362e-05,
     &	2.595e-06,
     &	2.9795e-05,
     &	2.36670201997668e-07,
     &	1.2137900734604e-05,
     &	2.91167958499589e-07,
     &	9.86605611925677e-08,
     &	2.01439011429255e-06,
     &	8.20622804366359e-08,
     &	7.83688694089992e-09,
     &	2.91167958499589e-05,
     &	1.52807116806281e-06
     &  /)

	molfracs_atoms(1)=molfracs_atoms(1)/(10d0**Z)
	molfracs_atoms(2)=molfracs_atoms(2)/(10d0**Z)
	if(enhancecarbon) then
		molfracs_atoms(3)=molfracs_atoms(5)*CO
	else
		molfracs_atoms(5)=molfracs_atoms(3)/CO
	endif
	
	molfracs_atoms(15)=molfracs_atoms(15)*TiScale

	tot=sum(molfracs_atoms(1:N_atoms))
	molfracs_atoms=molfracs_atoms/tot

	open(unit=50,file='atomic.dat')
	do i=1,18
		write(50,'(a5,se18.6)') names_atoms(i),molfracs_atoms(i)
	enddo
	close(unit=50)


	return
	end



	subroutine set_molfracs_atoms_elabun(CO,Z,elabun)
	use AtomsModule
	implicit none
	real*8 CO,Z,tot,elabun(7)
	integer i

	names_atoms(1) = 'H'
	names_atoms(2) = 'He'
	names_atoms(3) = 'C'
	names_atoms(4) = 'N'
	names_atoms(5) = 'O'
	names_atoms(6) = 'Na'
	names_atoms(7) = 'Mg'
	names_atoms(8) = 'Al'
	names_atoms(9) = 'Si'
	names_atoms(10) = 'P'
	names_atoms(11) = 'S'
	names_atoms(12) = 'Cl'
	names_atoms(13) = 'K'
	names_atoms(14) = 'Ca'
	names_atoms(15) = 'Ti'
	names_atoms(16) = 'V'
	names_atoms(17) = 'Fe'
	names_atoms(18) = 'Ni'

	molfracs_atoms = (/ 0.9207539305,
     &	0.0783688694,
     &	0.0002478241,
     &	6.22506056949881e-05,
     &	0.0004509658,
     &	1.60008694353205e-06,
     &	3.66558742055362e-05,
     &	2.595e-06,
     &	2.9795e-05,
     &	2.36670201997668e-07,
     &	1.2137900734604e-05,
     &	2.91167958499589e-07,
     &	9.86605611925677e-08,
     &	2.01439011429255e-06,
     &	8.20622804366359e-08,
     &	7.83688694089992e-09,
     &	2.91167958499589e-05,
     &	1.52807116806281e-06
     &  /)

	molfracs_atoms(1)=molfracs_atoms(1)/(10d0**Z)
	molfracs_atoms(2)=molfracs_atoms(2)/(10d0**Z)
	molfracs_atoms(5)=molfracs_atoms(3)/CO

	tot=sum(molfracs_atoms(1:N_atoms))
	molfracs_atoms=molfracs_atoms/tot

	molfracs_atoms( 7)=molfracs_atoms(1)*elabun(1)
	molfracs_atoms( 9)=molfracs_atoms(1)*elabun(2)
	molfracs_atoms(15)=molfracs_atoms(1)*elabun(3)
	molfracs_atoms( 5)=molfracs_atoms(1)*elabun(4)
	molfracs_atoms(17)=molfracs_atoms(1)*elabun(5)
	molfracs_atoms( 8)=molfracs_atoms(1)*elabun(6)
	molfracs_atoms( 3)=molfracs_atoms(1)*elabun(7)

	tot=sum(molfracs_atoms(1:N_atoms))
	molfracs_atoms=molfracs_atoms/tot

	return
	end



	subroutine call_chemistry(Tin,Pin,mol_abun,mol_names,nmol,ini,condensates,
     &		cloudspecies,Xcloud,Ncloud,nabla_ad,MMW,didcondens,includemol)
	use AtomsModule
	IMPLICIT NONE
	integer nmol
	real*8 Tin,Pin,mol_abun(nmol),nabla_ad
	character*10 mol_names(nmol)
	logical includemol(nmol),didcondens,ini,condensates
	integer Ncloud
	real*8 Xcloud(max(Ncloud,1)),MMW,Tg
	character*500 cloudspecies(max(Ncloud,1)),namecloud

	Tg=min(max(Tin,70d0),30000d0)

c	goto 1
	call call_GGchem(Tg,Pin,names_atoms,molfracs_atoms,N_atoms,mol_names,mol_abun,nmol,MMW,condensates)

	nabla_ad=2d0/7d0

	return

1	continue
	call call_easy_chem(Tg,Pin,mol_abun,mol_names,nmol,ini,condensates,
     &		cloudspecies,Xcloud,Ncloud,nabla_ad,MMW,didcondens,includemol)


	return
	end

