	subroutine SetupStructure(compute_mixrat)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 dp,dz,dlogp,RgasBar,Mtot,Pb(nr+1),tot,met_r,dens1bar,minZ,Tc,Rscale
	real*8 Otot,Ctot,Htot,vescape,vtherm,RHill,MMW_form,P0,R0,Kzz_r(nr)
	parameter(RgasBar=82.05736*1.01325)
	integer i,imol,nmix,j,niter,k,i1,i2,di,ii,i3,ir
	logical ini,compute_mixrat
	character*500 cloudspecies(max(nclouds,1))
	real*8 ComputeKzz
		
	do i=1,nclouds
		cloudspecies(i)=Cloud(i)%species
	enddo

	ini = .TRUE.

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
	if(par_tprofile.or.(computeT.and.nTiter.le.1)) call ComputeParamT(T)
	if(free_tprofile) call MakePTstruct

	call SetAbun
	nabla_ad=2d0/7d0

	Rscale=1d0

	grav=Ggrav*Mplanet/Rplanet**2

	do j=1,niter

	if(.not.mixratfile.and.(.not.dochemistry.or.j.eq.1)) then
	do i=1,nr
		tot=0d0
		mixrat_r(i,1:nmol)=mixrat(1:nmol)
		do imol=1,nmol
			if(mixrat_r(i,imol).gt.0d0) tot=tot+mixrat_r(i,imol)
		enddo
		if(tot.gt.0d0) mixrat_r(i,1:nmol)=mixrat_r(i,1:nmol)/tot
	enddo
	endif

	call output("log(g) [cgs]: " // dbl2string(log10(Ggrav*Mplanet/Rplanet**2),'(f8.3)'))

	Mtot=Mplanet

	if(j.eq.1) then
		MMW=0d0
		do imol=1,nmol
			if(includemol(imol)) MMW=MMW+mixrat_r(1,imol)*Mmol(imol)
		enddo
	endif
	call output("Mean molecular weight: " // dbl2string(MMW(1),'(f8.3)'))
	i1=1
	if(Pb(1).eq.Pplanet) Pb(1)=Pb(1)*1.001
	do i=2,nr
		if(Pb(i).eq.Pplanet) Pb(i)=Pb(i)**0.99*P(i)**0.01
	enddo
	if(Pb(nr+1).eq.Pplanet) Pb(nr+1)=Pb(nr+1)/1.001
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
				if(includemol(imol)) MMW(i3)=MMW(i3)+mixrat_r(i3,imol)*Mmol(imol)
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
		if(R(i+di).lt.R0/5d0) R(i+di)=R0/5d0
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
c		if(R(i+1).gt.RHill) then
c			call output("layer" // dbl2string(P(i),'(es10.3E3)') // "is beyond the Hill Sphere")
c			R(i+1)=sqrt(R(i)*RHill)
c		endif
		Mtot=Mtot+dens(i)*(R(i+1)**3-R(i)**3)*4d0*pi/3d0
		grav(i)=Ggrav*Mtot/(R(i)*R(i+1))
	enddo
c	if(constant_g) grav=Ggrav*Mplanet/(Rplanet)**2
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

	if(par_tprofile.or.(computeT.and.nTiter.le.1)) call ComputeParamT(T)
	if(free_tprofile) call MakePTstruct
	do i=1,nr
		if(T(i).gt.maxTprofile) T(i)=maxTprofile
		if(T(i).lt.3d0) T(i)=3d0
		if(.not.T(i).gt.3d0) T(i)=0.25d0*sqrt(Rstar/(2d0*Dplanet))*Tstar
	enddo

	if(domakeai) then
		modelsucces=.true.
		if(.not.modelsucces) return
	endif


	if(dochemistry.and.j.eq.1) then
	if(compute_mixrat) then
		call SetAbun
		call output("==================================================================")
c		call output("Computing chemistry using easy_chem by Paul Molliere")
		call output("Computing chemistry using GGchem by Peter Woitke")
		do i=1,nr
			call tellertje(i,nr)
			Tc=max(min(T(i),20000d0),100d0)
			if(P(i).ge.mixP.or.i.eq.1) then
				call call_chemistry(Tc,P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &			XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol)
    		else
    			mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
    			XeqCloud(i,1:nclouds)=XeqCloud(i-1,1:nclouds)
    			nabla_ad(i)=nabla_ad(i-1)
    			MMW(i)=MMW(i-1)
    			didcondens(i)=didcondens(i-1)
    		endif
		enddo
		if(disequilibrium) then
c call disequilibrium code
c input: 	R(1:nr+1) : These are the radial boundaries of the layers (bottom to top)
c			P(1:nr),T(1:nr) : These are the pressure and temperature inside the layers
c			molname(1:nmol) : names of the molecules included
c			Kzz_r(1:nr) : Diffusion coefficient
c input/output:	mixrat_r(1:nr,1:nmol) : number densities inside each layer. Now set to equilibrium abundances.
		   call output("==================================================================")
		   call output("Computing disequilibrium chemistry")
			do i=1,nr
				Kzz_r(i)=ComputeKzz(P(i),T(i),dens(i),complexKzz)
			enddo
		   call diseq_calc(nr,R(1:nr+1),P(1:nr),T(1:nr),nmol,molname(1:nmol),mixrat_r(1:nr, 1:nmol),COratio,Kzz_r(1:nr))
		   
		endif
		if(nfixmol.gt.0) then
			do i=1,nfixmol
				mixrat_r(1:nr,ifixmol(i))=fixmol_abun(i)
			enddo
		endif
		do i=1,nr
			do imol=1,nmol
				if(.not.mixrat_r(i,imol).gt.0d0) mixrat_r(i,imol)=0d0
			enddo
			if(.not.sum(mixrat_r(i,1:nmol)).gt.0d0.and.nmol.ge.48) then
				mixrat_r(i,1:nmol)=0d0
				mixrat_r(i,45)=0.85453462
				mixrat_r(i,48)=1d0-mixrat_r(i,45)
			endif
		enddo
		call output("==================================================================")
	else
		mixrat_r=mixrat_old_r
		XeqCloud=XeqCloud_old
	endif
	endif

	enddo

	if(compute_mixrat) then
		do i=1,nclouds
			call SetupCloud(i)
		enddo
	endif

	if(.not.dochemistry) then
		Otot=0d0
		Ctot=0d0
		Htot=0d0
		metallicity=0d0
		do i=1,nr
			do imol=1,nmol
				if(includemol(imol)) then
					Otot=Otot+Ndens(i)*mixrat_r(i,imol)*real(Oatoms(imol))
					Ctot=Ctot+Ndens(i)*mixrat_r(i,imol)*real(Catoms(imol))
					Htot=Htot+Ndens(i)*mixrat_r(i,imol)*real(Hatoms(imol))
					metallicity=metallicity+Ndens(i)*mixrat_r(i,imol)*real(tot_atoms(imol)-Hatoms(imol))
				endif
			enddo
		enddo
		COret=Ctot/Otot
		if(Tform.gt.0d0) COret=COratio
		call output("C/O: " // dbl2string(COret,'(f8.3)'))
		call output("[O]: " // dbl2string(log10(Otot/Htot)-log10(0.0004509658/0.9207539305),'(f8.3)'))
		call output("[C]: " // dbl2string(log10(Ctot/Htot)-log10(0.0002478241/0.9207539305),'(f8.3)'))

		COratio=COret
		if(includemol(48)) then
			do i=1,nr
				Htot=Htot+Ndens(i)*mixrat_r(i,48)
				metallicity=metallicity-Ndens(i)*mixrat_r(i,48)
			enddo
		endif
		metallicity=log10(metallicity/Htot)+3.0565202503263760
	endif

	if(.not.retrieval) then
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
	endif

	return
	end

	subroutine ComputeParamT(x)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x(nr),tau,Tirr,eta,expint,dlnT,dlnP,beta_used,max
	integer i

	if(.not.gamma_equal) then
		call ComputeParamT2(x)
		Tsurface=x(1)
		return
	endif

	beta_used=betaT
	if(deepredist) beta_used=f_deepredist

	tau=0d0
	Tirr=sqrt(Rstar/(Dplanet))*Tstar
	call output("Irradiation temperature: " // dbl2string(Tirr,'(f8.1)'))
c	if(computeT) then
c		x=(Tirr**4/sqrt(2d0)+TeffP**4)**0.25
c		return
c	endif
	
	do i=nr,1,-1
		if(constant_g) then
			tau=kappaT*1d6*P(i)/(Ggrav*Mplanet/(Rplanet)**2)
		else
			tau=kappaT*1d6*P(i)/grav(i)
		endif
		if(tau.lt.0d0) tau=0d0
		x(i)=(3d0*TeffP**4/4d0)*(2d0/3d0+tau)
c		if(deepRedist) then
c			print*,'warning! deep redistribution does not work properly with parameterised T-structure'
c			x(i)=x(i)+(3d0*Tirr**4*exp(-tau)/4d0/sqrt(2d0))*(2d0/3d0+tau)*(f_deepredist-beta_used)
c		endif
		x(i)=x(i)+(3d0*Tirr**4/4d0)*beta_used*
     &	(2d0/3d0+1d0/(sqrt(3d0)*gammaT1)+(gammaT1/sqrt(3d0)-1d0/(sqrt(3d0)*gammaT1))*exp(-gammaT1*tau*sqrt(3d0)))

		x(i)=x(i)**0.25d0

		if(x(i).gt.10000d0) x(i)=10000d0
		if(IsNaN(x(i))) then
			call output("NaN in temperature structure")
			if(i.lt.nr) then
				x(i)=x(i+1)
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
	Tsurface=x(1)

	return
	end

	subroutine ComputeParamT2(x)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x(nr),tau,Tirr,eta,expint,dlnT,dlnP,beta_used,max
	integer i

	beta_used=betaT
	if(deepredist) beta_used=f_deepredist

	tau=0d0
	Tirr=(beta_used*4d0)**0.25*sqrt(Rstar/(2d0*Dplanet))*Tstar
	call output("Irradiation temperature: " // dbl2string(Tirr,'(f8.1)'))
	if(computeT) then
		x=(Tirr**4+TeffP**4)**0.25
		return
	endif
	
	do i=nr,1,-1
		if(constant_g) then
			tau=kappaT*1d6*P(i)/(Ggrav*Mplanet/(Rplanet)**2)
		else
			tau=kappaT*1d6*P(i)/grav(i)
		endif
		if(tau.lt.0d0) tau=0d0
		x(i)=(3d0*TeffP**4/4d0)*(2d0/3d0+tau)
		if(tau.lt.100d0) then
			eta=2d0/3d0+(2d0/(3d0*gammaT1))*(1d0+(gammaT1*tau/2d0-1d0)*exp(-gammaT1*tau))
     &					+(2d0*gammaT1/3d0)*(1d0-tau**2/2d0)*expint(2,gammaT1*tau)
		else
			eta=2d0/3d0+(2d0/(3d0*gammaT1))
		endif
		x(i)=x(i)+(3d0*Tirr**4/4d0)*(1d0-alphaT)*eta
		if(tau.lt.100d0) then
			eta=2d0/3d0+(2d0/(3d0*gammaT2))*(1d0+(gammaT2*tau/2d0-1d0)*exp(-gammaT2*tau))
     &					+(2d0*gammaT2/3d0)*(1d0-tau**2/2d0)*expint(2,gammaT2*tau)
		else
			eta=2d0/3d0+(2d0/(3d0*gammaT2))
		endif
		x(i)=x(i)+(3d0*Tirr**4/4d0)*alphaT*eta
		x(i)=x(i)**0.25d0
		if(x(i).gt.10000d0) x(i)=10000d0
		if(IsNaN(x(i))) then
			call output("NaN in temperature structure")
			if(i.lt.nr) then
				x(i)=x(i+1)
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
	real*8 mix(nr,nmol),ComputeKzz
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
	write(50,'("#",a14,a15,a15,a13,a10,a11,a10)') "radius [cm]","height [cm]","dens [g/cm^3]","N [1/cm^3]","T [K]",
     & 	"P [bar]","g [cm/s^2]"
	do i=1,nr
		write(50,'(es15.7E3,es15.4E3,es15.4E3,es13.4E3,f10.3,es11.3E3,f10.3)') sqrt(R(i)*R(i+1)),sqrt(R(i)*R(i+1))-Rplanet
     &			,dens(i),Ndens(i),T(i),P(i),grav(i)
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

	open(unit=50,file=trim(outputdir) // 'mixingratios' // trim(side) // '.dat',RECL=6000)
	form='("#",a9,a13,' // trim(int2string(nmix,'(i2)')) // 'a15)'
	write(50,form) "T [K]","P [bar]",namemix(1:nmix)
	form='(f10.3,es13.3E3,' // trim(int2string(nmix,'(i2)')) // 'es15.4E3)'
	do i=1,nr
		write(50,form) T(i),P(i),mix(i,1:nmix)
	enddo
	close(unit=50)

	open(unit=50,file=trim(outputdir) // 'Kzz' // trim(side) // '.dat',RECL=6000)
	form='("#",a12,a13)'
	write(50,form) "Kzz [cm^2/s]","P [bar]"
	form='(es13.3E3,es13.3E3)'
	do i=1,nr
		write(50,form) ComputeKzz(P(i),T(i),dens(i),complexKzz),P(i)
	enddo
	close(unit=50)

	do ii=1,nclouds
		open(unit=50,file=trim(outputdir) // 'clouddens' // trim(int2string(ii,'(i0.2)')) // trim(side) // '.dat',RECL=6000)
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
	

	subroutine SetupCloud(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 pp,tot,column,tt,z,densdust,CloudMass
	integer ii,i,j,nsubr,isize,ilam
	real*8 Xc,Xc1,lambdaC,Ca,Cs,tau,P_SI1,P_SI2,veff,frac(nr,10)
	logical cl

	if(.not.allocated(Cloud(ii)%rv)) then
		allocate(Cloud(ii)%rv(nr))
		allocate(Cloud(ii)%sigma(nr))
	endif
	if(Cloud(ii)%simplecloud) then
		Cloud(ii)%ptype='SIMPLE'
		call SetupPartCloud(ii)
		Cloud(ii)%w=1d0
		return
	else if(cloudcompute) then
		call DiffuseCloud(ii)
		Cloud(ii)%rv=Cloud(ii)%rv*1d4
		Cloud(ii)%sigma=1d-10
		call output("Computing inhomogeneous cloud particles")
		call SetupPartCloud(ii)
		if(.not.retrieval) then
			open(unit=25,file=trim(outputdir) // "DRIFTcloud" // trim(int2string(ii,'(i0.4)')),recl=6000)
			do i=1,nr
				write(25,*) P(i),T(i),Cloud(ii)%rv(i),Cloud(ii)%sigma(i),cloud_dens(i,ii),dens(i),
     &			Cloud(ii)%frac(i,1)*3d0,Cloud(ii)%frac(i,4)*3d0,Cloud(ii)%frac(i,7:12),
     &			Cloud(ii)%frac(i,13)*3d0,Cloud(ii)%frac(i,16)
			enddo
			close(unit=25)
		endif
	else if(Cloud(ii)%file.ne.' ') then
		if(useDRIFT) then
		
		call regridN(Cloud(ii)%file,P*1d6,frac,nr,2,9,10,4,.true.,.true.)
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

		call regridN(Cloud(ii)%file,P*1d6,cloud_dens(1:nr,ii),nr,2,43,1,4,.true.,.false.)
		cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*dens(1:nr)

		call regridN(Cloud(ii)%file,P*1d6,Cloud(ii)%rv(1:nr),nr,2,7,1,4,.true.,.false.)
		
		else

		call regridN(Cloud(ii)%file,P,cloud_dens(1:nr,ii),nr,2,6,1,1,.false.,.false.)
		call regridN(Cloud(ii)%file,P,Cloud(ii)%rv(1:nr),nr,2,5,1,1,.false.,.true.)


		Cloud(ii)%frac(1:nr,1:19)=1d-10
		select case(Cloud(ii)%hazetype)
			case("SOOT","soot","THOLIN","tholin")
				Cloud(ii)%frac(1:nr,19)=1d0
				densdust=1.00
				Cloud(ii)%haze=.true.
			case("CHRIS")
c 10% iron
				Cloud(ii)%frac(1:nr,9)=0.1d0
c 90% MgSiO3
				Cloud(ii)%frac(1:nr,13:15)=0.9d0/3d0
				densdust=3.7
				Cloud(ii)%hazetype='THOLIN'
			case default
				call output("Cloud species unknown for file readin")
				stop
		end select
		
		cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*(4d0*pi*densdust*Cloud(ii)%rv(1:nr)**3)/3d0
		endif
		
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
	else if(.true.) then
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
			cloud_dens(1:nr,ii)=0d0
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
					Xc=(Xc1+(XeqCloud(i,ii)-XeqCloud(i-1,ii))/real(nsubr))/(1d0-Cloud(ii)%frain*((P(i)-P(i-1))
     &						/real(nsubr))/(lambdaC*P(i)))
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
	veff=Cloud(ii)%veff
1	tot=0d0
	do i=1,Cloud(ii)%nr
		Cloud(ii)%w(i)=(1d4*Cloud(ii)%rv(i))**(1d0+(1d0-3d0*veff)/veff)*
     &					exp(-1d0*Cloud(ii)%rv(i)/(Cloud(ii)%reff*veff))
		Cloud(ii)%w(i)=Cloud(ii)%w(i)*Cloud(ii)%rv(i)**3
		tot=tot+Cloud(ii)%w(i)
	enddo
	if(.not.tot.gt.0d0) then
		if(j.lt.10) then
			veff=veff*2d0
			call output("increasing veff")
			j=j+1
			goto 1
		else
			tot=1d0
			Cloud(ii)%w(1:Cloud(ii)%nr)=1d0
		endif
	endif
	Cloud(ii)%w=Cloud(ii)%w/tot
	
	if(Cloud(ii)%tau.gt.0d0.and..not.cloudcompute) then
		tau=0d0
		do ilam=nlam-1,1,-1
			if((lam(ilam)*1d4).lt.Cloud(ii)%lam.and.(lam(ilam+1)*1d4).ge.Cloud(ii)%lam) exit
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

	CloudMass=0d0
	do i=1,nr
		CloudMass=CloudMass+cloud_dens(i,ii)*4d0*pi*abs(R(i+1)**3-R(i)**3)/3d0
	enddo
	call output("Cloud mass: " // dbl2string(CloudMass/Mearth,'(es13.4E3)') // " Mearth")

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
	call set_molfracs_atoms(COratio,SiOratio,NOratio,SOratio,metallicity)
	ini=.true.
	do i=1,nr
		call call_chemistry(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol)
	enddo

	return
	end
	
	
	subroutine SetAbun()
	use GlobalSetup
	use Constants
	use AtomsModule
	IMPLICIT NONE
	real*8 Zout,abun_dust(N_atoms)
	integer i

	if(secondary_atmosphere) then
		if(.not.computeT) Tsurface=T(1)
		if(Toutgas.eq.Tsurface.and.Poutgas.eq.P(1)) then
			molfracs_atoms(1:N_atoms)=molfracs_atoms_outgas(1:N_atoms)
			return
		endif
	endif

	if(element_abun_file.ne.' ') then
		call read_molfracs_atoms(element_abun_file,COratio,metallicity)
	else
		call set_molfracs_atoms(COratio,SiOratio,NOratio,SOratio,metallicity)
	endif

	if(secondary_atmosphere) then
		Toutgas=Tsurface
		Poutgas=P(1)
		call SplitGasDust(Toutgas,Poutgas,molfracs_atoms,molfracs_atoms_outgas,abun_dust)
	endif
	
	return
	end


	subroutine SplitGasDust(Temp,Pres,abun_total,abun_gas,abun_dust)
	use GlobalSetup
	use Constants
	use AtomsModule
	IMPLICIT NONE
	real*8 mutemp,mol_abun(nmol),abun_total(N_atoms),abun_gas(N_atoms),abun_dust(N_atoms)
	real*8 Pres,Temp
	integer methGGchem
	
	methGGchem=0

	call call_GGchem(Temp,Pres,names_atoms,abun_total,N_atoms,molname(1:nmol),
     &			mol_abun,nmol,mutemp,.true.,abun_gas,methGGchem)
	abun_dust=abun_total-abun_gas	

	return
	end
	

	subroutine set_molfracs_atoms(CO,SiO,NO,SO,Z)
	use GlobalSetup
	use Constants
	use AtomsModule
	implicit none
	real*8 CO,Z,tot,SiO,Z0,scale,NO,CO0,SO
	integer i
	character*500 command,homedir
	character*10 name
	real*8 abun,Mcore,Rstart

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

	molfracs_atoms = (/ 0.9207539305, 	!H
     &	0.0783688694,					!He
     &	0.0002478241,					!C
     &	6.22506056949881e-05,			!N
     &	0.0004509658,					!O
     &	1.60008694353205e-06,			!Na
     &	3.66558742055362e-05,			!Mg
     &	2.595e-06,						!Al
     &	2.9795e-05,						!Si
     &	2.36670201997668e-07,			!P
     &	1.2137900734604e-05,			!S
     &	2.91167958499589e-07,			!Cl
     &	9.86605611925677e-08,			!K
     &	2.01439011429255e-06,			!Ca
     &	8.20622804366359e-08,			!Ti
     &	7.83688694089992e-09,			!V
     &	2.91167958499589e-05,			!Fe
     &	1.52807116806281e-06			!Ni
     &  /)


	if(planetform) then

	if(.false.) then	!run the python version

	Z0=sum(molfracs_atoms(3:N_atoms))/sum(molfracs_atoms(1:2))
	
	if((planetform_fdust+planetform_fplan).ge.1d0) then	
		scale=planetform_fdust+planetform_fplan
		planetform_fdust=planetform_fdust/scale
		planetform_fplan=planetform_fplan/scale
	endif
	
	call getenv('HOME',homedir) 
	write(command,'("python3 ",a,"/ARCiS/Data/planetform/elements.py ",a,"/atomic.dat ",6es15.4)') 
     &	trim(homedir),trim(outputdir),Dplanet/AU,Mplanet/Mearth,planetform_Rstart,planetform_Mstart,planetform_fdust,planetform_fplan
	if(trim(command).ne.trim(formationcommand)) then
		call system(command)
		formationcommand=command
	endif

	molfracs_atoms=1d-200
	open(unit=43,file=trim(outputdir) // 'atomic.dat')
1	read(43,*,err=1,end=2) name,abun
	do i=1,N_atoms
		if(trim(name).eq.names_atoms(i)) then
			molfracs_atoms(i)=abun
			exit
		endif
	enddo
	goto 1
2	close(unit=43)

	else	! run the fortran version

c Setup names and weights of the elements
	Z0=sum(molfracs_atoms(3:N_atoms))/sum(molfracs_atoms(1:2))

	CO=molfracs_atoms(3)/molfracs_atoms(5)
	SiO=molfracs_atoms(9)/molfracs_atoms(5)
	NO=molfracs_atoms(4)/molfracs_atoms(5)
	SO=molfracs_atoms(11)/molfracs_atoms(5)
	Z=log10(sum(molfracs_atoms(3:N_atoms))/sum(molfracs_atoms(1:2))/Z0)

	call output("Solar abundances:")
	call output("C/O: " // dbl2string(CO,'(f7.3)'))
	call output("Si/O:" // dbl2string(SiO,'(f7.3)'))
	call output("N/O: " // dbl2string(NO,'(f7.3)'))
	call output("S/O: " // dbl2string(SO,'(f7.3)'))
	call output("[Z]: " // dbl2string(Z,'(f7.3)'))

	Mcore=planetform_Mstart*Mearth
	Rstart=planetform_Rstart*AU

	call InitFormation(Mstar,Tstar,Rstar,planetform_SolidC,planetform_Macc)
	call Formation(Mplanet,Mcore,Rstart,Dplanet,planetform_fdust,planetform_fplan,simAb_converge,MSimAb)
	if(.not.simAb_converge) then
		call output("WARNING: SimAb did not converge!")
		call output("WARNING: final mass of the planet:" // dbl2string(MSimAb/Mjup,'(se10.3)') // " Mjup")
	endif

	endif

	molfracs_atoms=molfracs_atoms/sum(molfracs_atoms(1:N_atoms))
	do i=1,N_atoms
		if(molfracs_atoms(i).lt.1d-50) molfracs_atoms(i)=1d-50
	enddo

	CO=molfracs_atoms(3)/molfracs_atoms(5)
	SiO=molfracs_atoms(9)/molfracs_atoms(5)
	NO=molfracs_atoms(4)/molfracs_atoms(5)
	SO=molfracs_atoms(11)/molfracs_atoms(5)
	Z=log10(sum(molfracs_atoms(3:N_atoms))/sum(molfracs_atoms(1:2))/Z0)

	call output("SimAB abundances:")
	call output("C/O: " // dbl2string(CO,'(f7.3)'))
	call output("Si/O:" // dbl2string(SiO,'(f7.3)'))
	call output("N/O: " // dbl2string(NO,'(f7.3)'))
	call output("S/O: " // dbl2string(SO,'(f7.3)'))
	call output("[Z]: " // dbl2string(Z,'(f7.3)'))

	if(.not.retrieval) then
	open(unit=50,file=trim(outputdir) // 'atomic.dat')
	write(50,'("SimAB abundances:")')
	write(50,'(a6,f7.3)') "C/O:",CO
	write(50,'(a6,f7.3)') "Si/O:",SiO
	write(50,'(a6,f7.3)') "N/O:",NO
	write(50,'(a6,f7.3)') "S/O:",SO
	write(50,'(a6,f7.3)') "[Z]:",Z
	write(50,'("=============")')
	do i=1,18
		write(50,'(a5,se18.6)') names_atoms(i),molfracs_atoms(i)
	enddo
	close(unit=50)
	endif

	return

	else

	Z0=sum(molfracs_atoms(3:N_atoms))/sum(molfracs_atoms(1:2))

c	adjust C/O ratio
	CO0=molfracs_atoms(3)/molfracs_atoms(5)
	if((.not.inverseCOratio.and.CO.lt.CO0).or.(inverseCOratio.and.CO.gt.CO0)) then
		molfracs_atoms(3)=molfracs_atoms(5)*CO
	else
		molfracs_atoms(5)=molfracs_atoms(3)/CO
	endif

c	adjust N/O ratio
	molfracs_atoms(4)=molfracs_atoms(5)*NO

c	adjust Si/O ratio
	scale=molfracs_atoms(5)*SiO/molfracs_atoms(9)
	molfracs_atoms(6:18)=molfracs_atoms(6:18)*scale

c	adjust S/O ratio
	molfracs_atoms(11)=molfracs_atoms(5)*SO

c	adjust metallicity
	scale=sum(molfracs_atoms(3:N_atoms))/(sum(molfracs_atoms(1:2))*Z0*10d0**Z)
	molfracs_atoms(1:2)=molfracs_atoms(1:2)*scale

	tot=sum(molfracs_atoms(1:N_atoms))
	molfracs_atoms=molfracs_atoms/tot

c	open(unit=50,file='atomic.dat')
c	do i=1,18
c		write(50,'(a5,se18.6)') names_atoms(i),molfracs_atoms(i)
c	enddo
c	close(unit=50)

	endif

	return
	end

	
	subroutine read_molfracs_atoms(filename,COratio,metallicity)
	use AtomsModule
	use OutputModule
	IMPLICIT NONE
	character*500 filename
	real*8 molfracs_atoms_solar(N_atoms),abun,COratio,metallicity
	character*40 name
	logical readin(N_atoms)
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

	molfracs_atoms_solar = (/ 0.9207539305,
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

	readin=.false.
	molfracs_atoms=molfracs_atoms_solar

	open(unit=43,file=filename)
1	read(43,*,err=1,end=2) name,abun
	do i=1,N_atoms
		if(trim(name).eq.names_atoms(i)) then
			molfracs_atoms(i)=abun
			readin(i)=.true.
		endif
	enddo
	goto 1
2	close(unit=43)

	molfracs_atoms=molfracs_atoms/sum(molfracs_atoms(1:N_atoms))

	molfracs_atoms_solar=molfracs_atoms_solar/sum(molfracs_atoms_solar(1:N_atoms))

	molfracs_atoms=molfracs_atoms/sum(molfracs_atoms(1:N_atoms))

	COratio=molfracs_atoms(3)/molfracs_atoms(5)
	call output("C/O: " // dbl2string(molfracs_atoms(3)/molfracs_atoms(5),'(f6.2)'))
	call output("Si/O:" // dbl2string(molfracs_atoms(9)/molfracs_atoms(5),'(f6.2)'))

	metallicity=log10((sum(molfracs_atoms(3:n_atoms))/sum(molfracs_atoms(1:n_atoms)))/
     &			(sum(molfracs_atoms_solar(3:n_atoms))/sum(molfracs_atoms_solar(1:n_atoms))))
	call output("[Z]: " // dbl2string(metallicity,'(f6.2)'))

	open(unit=50,file='atomic.dat')
	do i=1,18
		write(50,'(a5,se18.6)') names_atoms(i),molfracs_atoms(i)
	enddo
	close(unit=50)

	return
	end


	subroutine set_molfracs_atoms_old(CO,Z,TiScale,enhancecarbon)
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
	molfracs_atoms(16)=molfracs_atoms(16)*TiScale

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
	use TimingModule
	IMPLICIT NONE
	integer nmol
	real*8 Tin,Pin,mol_abun(nmol),nabla_ad
	character*10 mol_names(nmol)
	logical includemol(nmol),didcondens,ini,condensates
	integer Ncloud,i,imol,methGGchem
	real*8 Xcloud(max(Ncloud,1)),MMW,Tg
	character*500 cloudspecies(max(Ncloud,1)),namecloud
	real*8 P1,P2,abun_temp(nmol),M,gas_atoms(N_atoms)
	real*8 time
	integer itime

	if(ini) then
		methGGchem=0
	else
		methGGchem=2
	endif
	ini=.false.

	call cpu_time(time)
	timechem=timechem-time
	call system_clock(itime)
	itimechem=itimechem-itime
	ctimechem=ctimechem+1

	Tg=min(max(Tin,100d0),30000d0)

	mol_abun=0d0
	Xcloud=0d0
	call call_GGchem(Tg,Pin,names_atoms,molfracs_atoms,N_atoms,mol_names,mol_abun,nmol,MMW,condensates,gas_atoms,methGGchem)

c	call readBaud(mol_abun,nmol,Pin,MMW)

c	mol_abun=1d-4
c	MMW=2.2

	do i=1,nmol
		if(.not.mol_abun(i).gt.0d0) mol_abun(i)=0d0
	enddo

	nabla_ad=2d0/7d0

	call cpu_time(time)
	timechem=timechem+time
	call system_clock(itime)
	itimechem=itimechem+itime

	return

	call call_easy_chem(Tin,Pin,mol_abun,mol_names,nmol,ini,condensates,
     &		cloudspecies,Xcloud,Ncloud,nabla_ad,MMW,didcondens,includemol)
	
	return
	end

	subroutine readBaud(mf,nm,Pin,mm)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nm,i
	real*8 mf(nm),Pin,P0,mm
	
	open(unit=20,file='dbf26/a_p_G4___.tsv')
1	read(20,*,end=2) P0,mf(1),mf(6),mf(11),mf(5),mf(2),mf(28),mf(56),mf(57)
	if(P0.ge.Pin) then
		mm=0d0
		do i=1,nm
			if(includemol(i)) mm=mm+mf(i)*Mmol(i)
		enddo
		mm=mm/sum(mf(1:nm))
		close(unit=20)
		return
	endif
	goto 1
2	continue
	close(unit=20)
	return
	end
	


	subroutine ComputeBeta(theta,t,beta)
	IMPLICIT NONE
	integer n
	real*8 pi
	parameter(n=1000)
	parameter(pi=3.1415926536)
	real*8 theta,t,beta,b(n),b0(n),th(n),diff,tot,tot0
	integer i,j

	if(t.lt.0d0) then
		beta=1d0
		return
	endif

	tot=0d0
	tot0=0d0
	b=0d0
	do i=1,n
		th(i)=2d0*pi*real(i-1)/real(n)
		b0(i)=cos(th(i))
		if(b0(i).lt.0d0) b0(i)=0d0
		tot0=tot0+b0(i)
	enddo
	do i=1,n
		do j=1,n
			diff=(th(i)-th(j))
			if(diff.lt.0d0) diff=diff+2d0*pi
			diff=diff*t/(2d0*pi)
			b(i)=b(i)+b0(j)*exp(-diff)
		enddo
		tot=tot+b(i)
	enddo
	b=b*tot0/tot

	beta=0d0
	do i=1,n
		diff=(theta*pi/180d0-th(i))
		if(diff.lt.0d0) diff=diff+2d0*pi
		diff=diff*t/(2d0*pi)
		beta=beta+b0(i)*exp(-diff)
	enddo
	beta=beta*tot0/tot

	return
	end
	

	subroutine MakePTstruct()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 Tc,Ppoint(nVpoints+nIRpoints),Tpoint(nVpoints+nIRpoints),PrefTpoint,Gplan
	real*8 TV(nr),TIR(nr)
	integer i

	TV=0d0
	TIR=0d0
	
	Gplan=(Ggrav*Mplanet/(Rplanet)**2)
	PrefTpoint=Gplan/(kappaT*1d6)
	
	if(nVpoints.gt.0) then
		do i=1,nVpoints
			Ppoint(i)=betaT*tau_Vpoint(i)*Gplan/(kappaT*gammaT1*1d6)
			Tpoint(i)=dT_Vpoint(i)
		enddo
		Tc=(Tstar**4*(Rstar/Dplanet)**2*(betaT*gammaT2))**0.25
		call MakePTstruct_dT(P,TV,nr,Ppoint,Tpoint,nVpoints,Tc,PrefTpoint)
	endif
	if(nIRpoints.gt.0) then
		do i=1,nIRpoints
			Ppoint(i)=tau_IRpoint(i)*Gplan/(kappaT*1d6)
			Tpoint(i)=dT_IRpoint(i)
		enddo
		Tc=TeffP
		call MakePTstruct_dT(P,TIR,nr,Ppoint,Tpoint,nIRpoints,Tc,PrefTpoint)
	endif
	T=(TV**4+TIR**4)**0.25
	
	return
	end
	

	subroutine MakePTstruct_dT(P,T,np,Pp_in,dTp_in,nT_in,T0,P0)
	IMPLICIT NONE
	integer np,i,nT,nT_in,j,i0
	real*8 P(np),T(np),Pp_in(nT_in),dTp_in(nT_in),dTp(nT_in),P0,Pp(nT_in)
	real*8 logPp(nT_in),logP(np),logT(np),logT0,T0,logTp(nT_in),dT0,Tp(nT_in),logP0
	real*8 a,b,c

	logPp=log(Pp_in+1d-30)
	dTp=dTp_in

	nT=nT_in
	
	call sortw(logPp,dTp,nT)
1	continue
	do i=1,nT-1
		if(logPp(i).eq.logPp(i+1)) then
			logPp(i)=logPp(i)-1d-4
			goto 1
		endif
	enddo

	logPp=-logPp
	call sortw(logPp,dTp,nT)
	logPp=-logPp
	Pp=exp(logPp)

	logP0=log(P0)
	logT0=log(T0)
	logP=log(P)

	if(logP0.le.logPp(nT)) then
		a=dTp(nT)/Pp(nT)
		b=logT0-a*P0
		logTp(nT)=a*Pp(nT)+b
		do i=nT-1,1,-1
			a=(dTp(i+1)-dTp(i))/(2d0*(logPp(i+1)-logPp(i)))
			b=dTp(i+1)-2d0*a*logPp(i+1)
			c=logTp(i+1)-b*logPp(i+1)-a*logPp(i+1)**2
			logTp(i)=a*logPp(i)**2+b*logPp(i)+c
		enddo
	else if(logP0.ge.logPp(1)) then
		a=dTp(1)
		b=logT0-a*logP0
		logTp(1)=a*logPp(1)+b
		do i=2,nT
			a=(dTp(i-1)-dTp(i))/(2d0*(logPp(i-1)-logPp(i)))
			b=dTp(i-1)-2d0*a*logPp(i-1)
			c=logTp(i-1)-b*logPp(i-1)-a*logPp(i-1)**2
			logTp(i)=a*logPp(i)**2+b*logPp(i)+c
		enddo
	else
		do i0=1,nT-1
			if(logP0.le.logPp(i0).and.logP0.gt.logPp(i0+1)) exit
		enddo
		a=(dTp(i0)-dTp(i0+1))/(2d0*(logPp(i0)-logPp(i0+1)))
		b=dTp(i0)-2d0*a*logPp(i0)
		c=logT0-b*logP0-a*logP0**2
		logTp(i0)=a*logPp(i0)**2+b*logPp(i0)+c
		do i=i0-1,1,-1
			a=(dTp(i+1)-dTp(i))/(2d0*(logPp(i+1)-logPp(i)))
			b=dTp(i+1)-2d0*a*logPp(i+1)
			c=logTp(i+1)-b*logPp(i+1)-a*logPp(i+1)**2
			logTp(i)=a*logPp(i)**2+b*logPp(i)+c
		enddo
		do i=i0+1,nT
			a=(dTp(i-1)-dTp(i))/(2d0*(logPp(i-1)-logPp(i)))
			b=dTp(i-1)-2d0*a*logPp(i-1)
			c=logTp(i-1)-b*logPp(i-1)-a*logPp(i-1)**2
			logTp(i)=a*logPp(i)**2+b*logPp(i)+c
		enddo
	endif

	do i=np,1,-1
		if(logP(i).le.logPp(nT)) then
			a=dTp(nT)/Pp(nT)
			b=logTp(nT)-a*Pp(nT)
			logT(i)=a*P(i)+b
		else if(logP(i).ge.logPp(1)) then
			a=dTp(1)
			b=logTp(1)-a*logPp(1)
			logT(i)=a*logP(i)+b
		else
			do i0=1,nT-1
				if(logP(i).le.logPp(i0).and.logP(i).gt.logPp(i0+1)) exit
			enddo
			a=(dTp(i0)-dTp(i0+1))/(2d0*(logPp(i0)-logPp(i0+1)))
			b=dTp(i0)-2d0*a*logPp(i0)
			c=logTp(i0)-b*logPp(i0)-a*logPp(i0)**2
			logT(i)=a*logP(i)**2+b*logP(i)+c
		endif
		T(i)=exp(logT(i))
	enddo
	
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	real*8 function erfinv(x)
	IMPLICIT NONE
	real*8 x,f,eps,y,pi,ck(0:100)
	integer k,m
	parameter(pi=3.141592653589793238462643383279502884)
	parameter(eps=1d-8)
	y=sqrt(pi)*x/2d0

	erfinv=0d0
	k=0
	ck(0)=1d0
1	continue
	f=ck(k)*(y**(2d0*real(k)+1d0))/(2d0*real(k)+1d0)
	erfinv=erfinv+f
	if(abs(f).gt.eps.and.k.lt.100) then
		k=k+1
		ck(k)=0d0
		do m=0,k-1
			ck(k)=ck(k)+ck(m)*ck(k-1-m)/real((m+1)*(2*m+1))
		enddo
		goto 1
	endif

	return
	end



	real*8 function ComputeKzz(x,Tg,rhog,addmicro)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,Tg,rhog,lmfp,vth,sigmamol,Kmax,Kmin,Kp
	logical addmicro
	sigmamol=8d-15

	if(Kzz_deep.gt.0d0.and.Kzz_1bar.gt.0d0) then
		if(Kzz_contrast.gt.1d0) then
			Kmax=Kzz_deep*Kzz_contrast
			Kmin=Kzz_deep
			Kp=abs(Kzz_P)
		else
			Kmax=Kzz_deep
			Kmin=Kzz_deep*Kzz_contrast
			Kp=-abs(Kzz_P)
		endif
c		ComputeKzz=1d0/(1d0/Kmax+1d0/(Kmin+Kzz_1bar/x**Kp))
		ComputeKzz=min(max(Kzz_1bar/x**Kp,Kmin),Kmax)
	else
		ComputeKzz=Kzz
	endif
	if(addmicro) then
		lmfp=2.3d0*mp/(sqrt(2d0)*rhog*sigmamol)
		vth=sqrt(8d0*kb*Tg/(pi*2.3d0*mp))
		ComputeKzz=ComputeKzz+lmfp*vth/3d0
	endif
	
	return
	end
	