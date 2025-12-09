	subroutine SetupStructure(compute_mixrat)
	use GlobalSetup
	use Constants
	use AtomsModule
	use CloudModule
	IMPLICIT NONE
	real*8 dp,dz,dlogp,RgasBar,Mtot,Pb(nr+1),tot,tot2,met_r,dens1bar,minZ,Tc,Rscale
	real*8 Otot,Ctot,Htot,vescape,vtherm,RHill,MMW_form,P0,R0,Kzz_g_in(nr),Kzz_b_in(nr)
	parameter(RgasBar=82.05736d0*1.01325d0)
	integer i,imol,nmix,j,niter,k,i1,i2,di,ii,i3,ir,icc
	logical ini,compute_mixrat
	character*500 cloudspecies(max(nclouds,1))
	real*8 ComputeKzz,molfracs_atoms0(N_atoms)
	real*8,allocatable :: Otot_r(:),Ctot_r(:),Ntot_r(:)
	logical dochemR(nr)
	
	Kzz_g_in=Kzz_g
	Kzz_b_in=Kzz_b
	do i=1,nclouds
		cloudspecies(i)=Cloud(i)%species
	enddo
	
	if(ComputeTeff) then
c	Use Thorngren & Fortney (2018)
		Tc=sqrt(Rstar/(Dplanet))*Tstar
		TeffP=0.39*Tc*exp(-(log10(4d-9*sigma*Tc**4)-0.14)**2/1.095)
c	Add a minimum of 85K for compatibility with cold/old gas Giants
		TeffP=(TeffP**4+85d0**4)**0.25
		if(.not.retrieval) call output("Internal temperature: " // dbl2string(TeffP,'(f6.1)') // "K")
	endif
	
	if(.not.do3D) betaF=betaT

	ini = .TRUE.

	minZ=-5d0
	mixrat_optEC_r=0d0

	niter=3
	if(dochemistry) niter=3
	if(par_tprofile) niter=3

	Pb(1)=P(1)
	do i=2,nr
		Pb(i)=sqrt(P(i-1)*P(i))
	enddo
	Pb(nr+1)=P(nr)

	nabla_ad=(1d0-1d0/exp_ad)	!2d0/7d0
	grav=Ggrav*Mplanet/(Rplanet)**2
	if(((par_tprofile.and..not.computeT).or.(par_tprofile.and.computeT.and.nTiter.le.1)).and.i_alb.le.1) 
     &			call ComputeParamT(T)
	if(free_tprofile.and.(.not.computeT.or.nTiter.le.1).and.i_alb.le.1) call MakePTstruct
	if(WaterWorld.and.(.not.do3D.or.setsurfpressure)) call doWaterWorld()

	call SetAbun

	Rscale=1d0

	grav=Ggrav*Mplanet/Rplanet**2

	do j=1,niter

	if(.not.mixratfile.and.(.not.dochemistry.or.j.eq.1)) then
		do i=1,nr
			tot=0d0
			mixrat_r(i,1:nmol)=mixrat(1:nmol)
			do imol=1,nmol
				if(P(i).lt.Pswitch_mol(imol)) mixrat_r(i,imol)=abun_switch_mol(imol)
			enddo
			if(dobackgroundgas) then
				do imol=1,nmol
					if(mixrat_r(i,imol).gt.0d0.and..not.backgroundgas(imol)) tot=tot+mixrat_r(i,imol)
				enddo
				if(tot.gt.1d0) then
					modelfail=.true.
					return
				else
					tot2=0d0
					do imol=1,nmol
						if(mixrat_r(i,imol).gt.0d0.and.backgroundgas(imol)) tot2=tot2+mixrat_r(i,imol)
					enddo
					do imol=1,nmol
						if(mixrat_r(i,imol).gt.0d0.and.backgroundgas(imol)) mixrat_r(i,imol)=mixrat_r(i,imol)*(1d0-tot)/tot2
					enddo
				endif
			else
				do imol=1,nmol
					if(mixrat_r(i,imol).gt.0d0) tot=tot+mixrat_r(i,imol)
				enddo
				if(tot.gt.0d0) mixrat_r(i,1:nmol)=mixrat_r(i,1:nmol)/tot
			endif
		enddo
		mixrat_optEC_r=0d0
		call doPhotoChemMol()
		do imol=1,nmol
			if(isotope(imol).gt.0) then
				do i=1,nr
					mixrat_r(i,isotope(imol))=mixrat_r(i,imol)/f_isotope(imol)
					mixrat_r(i,imol)=mixrat_r(i,imol)*(1d0-1d0/f_isotope(imol))
				enddo
			endif
		enddo
	endif

	call output("log(g) [cgs]: " // dbl2string(log10(Ggrav*Mplanet/Rplanet**2),'(f8.3)'))

	Mtot=Mplanet

	if(j.eq.1) then
		if(fixMMW) then
			MMW=MMW0
		else
			MMW=0d0
			do imol=1,nmol
				if(includemol(imol)) MMW=MMW+mixrat_r(1,imol)*Mmol(imol)
			enddo
		endif
	endif
	call output("Mean molecular weight: " // dbl2string(MMW(1),'(f8.3)'))
	i1=1
	if(Pb(1).eq.Pplanet) Pb(1)=Pb(1)*1.001
	do i=2,nr
		if(Pb(i).eq.Pplanet) Pb(i)=Pb(i)**0.99*P(i)**0.01
	enddo
	if(Pb(nr+1).eq.Pplanet) Pb(nr+1)=Pb(nr+1)/1.001
	if(.not.outflow) then
	P0=min(Pplanet,Pb(1))
	do i=1,nr
		if(Pb(i).ge.P0.and.Pb(i+1).lt.P0) then
			i1=i
		endif
	enddo
	i2=nr
	di=1
	do ii=1,2
	P0=min(Pplanet,Pb(1))
	R0=Rplanet
	do i=i1,i2,di
		if(di.gt.0) then
			i3=i
		else
			i3=i-1
		endif

		if(.not.dochemistry.and..not.fixMMW) then
			MMW(i3)=0d0
			do imol=1,nmol
				if(includemol(imol)) MMW(i3)=MMW(i3)+mixrat_r(i3,imol)*Mmol(imol)
			enddo
			if(mixratfile) then
				MMW(i3)=MMW(i3)+(1d0-sum(mixrat_r(i3,1:nmol)))
			endif
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
	P0=min(Pplanet,Pb(1))
	enddo
	else
	R(1)=Rplanet
	do i=1,nr
		if(.not.dochemistry.and..not.fixMMW) then
			MMW(i)=0d0
			do imol=1,nmol
				if(includemol(imol)) MMW(i)=MMW(i)+mixrat_r(i,imol)*Mmol(imol)
			enddo
			if(mixratfile) then
				MMW(i)=MMW(i)+(1d0-sum(mixrat_r(i,1:nmol)))
			endif
		endif
				
		Ndens(i)=P(i)*Avogadro/(RgasBar*T(i))
		Hp(i3)=(T(i)*kb)/(grav(i)*mp*MMW(i))

		dz=dlogp*Hp(i)
		dens(i)=Ndens(i)*mp*MMW(i)

		if(i.lt.nr) then
			if(Pb(i).ge.1d0.and.Pb(i+1).lt.1d0) then
				dens1bar=Avogadro*mp*MMW(i)/(RgasBar*T(i))
			endif
		endif
		if(i.gt.1.and.i.le.nr) then
			R(i)=sqrt(R(i-1)**2*dens(i-1)/dens(i))
		endif
	enddo
	R(nr+1)=sqrt(R(nr)**2*dens(nr-1)/dens(nr))
	endif
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
	if(constant_g) grav=Ggrav*Mplanet/(Rplanet)**2

	if(((par_tprofile.and..not.computeT).or.(par_tprofile.and.computeT.and.nTiter.le.1)).and.i_alb.le.1) 
     &			call ComputeParamT(T)
	if(free_tprofile.and.(.not.computeT.or.nTiter.le.1).and.i_alb.le.1) call MakePTstruct
	do i=1,nr
		if(T(i).gt.maxTprofile) T(i)=maxTprofile
		if(T(i).lt.3d0) T(i)=3d0
		if(.not.T(i).gt.3d0) T(i)=0.25d0*sqrt(Rstar/(2d0*Dplanet))*Tstar
	enddo

	if(domakeai) then
		modelsucces=.true.
		if(.not.modelsucces) return
	endif

	do i=1,nr
		if(nTiter.le.1) then
			Kzz_b(i)=ComputeKzz(P(i),T(i),dens(i),Hp(i),.false.)+Kzz_convect(i)
			if(complexKzz) then
				Kzz_g(i)=ComputeKzz(P(i),T(i),dens(i),Hp(i),complexKzz)+Kzz_convect(i)
			else
				Kzz_g(i)=Kzz_b(i)
			endif
		else
			Kzz_b(i)=exp((log(Kzz_b_in(i))*real(nTiter-1)+log(ComputeKzz(P(i),T(i),dens(i),Hp(i),.false.)+Kzz_convect(i)))/real(nTiter))
			if(complexKzz) then
				Kzz_g(i)=exp((log(Kzz_g_in(i))*real(nTiter-1)+log(ComputeKzz(P(i),T(i),dens(i),Hp(i),.true.)+Kzz_convect(i)))/real(nTiter))
			else
				Kzz_g(i)=Kzz_b(i)
			endif
		endif		
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

	if(dochemistry.and.j.eq.1) then
	if(compute_mixrat) then
		call SetAbun
		call output("==================================================================")
c		call output("Computing chemistry using easy_chem by Paul Molliere")
		call output("Computing chemistry using GGchem by Peter Woitke")

		dochemR=.false.
		dochemR(1)=.true.
		dochemR(nr)=.true.
		do i=1,nr,nrstepchem
			dochemR(i)=.true.
		enddo

		do i=1,nr
			call tellertje(i,nr)
			if(nPhotoReacts.gt.0) then
				molfracs_atoms0=molfracs_atoms
				call doPhotoChemAtom(i)
			endif
			if(dochemR(i)) then
			Tc=max(min(T(i),20000d0),100d0)
			do ii=1,nclouds
				if(nTiter.gt.Cloud(ii)%fixcloud.and.Cloud(ii)%fixcloud.gt.0) molfracs_atoms(1:N_atoms)=molfracs_atoms_cloud(1:N_atoms,i)
			enddo
			if(P(i).ge.mixP.or.i.eq.1) then
				call call_chemistry(Tc,P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &			XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol,dosimplerainout,useEOS,x_el(i))
    		else
    			mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
    			XeqCloud(i,1:nclouds)=XeqCloud(i-1,1:nclouds)
    			nabla_ad(i)=nabla_ad(i-1)
    			MMW(i)=MMW(i-1)
    			didcondens(i)=didcondens(i-1)
    		endif
			endif
			if(nPhotoReacts.gt.0) then
				molfracs_atoms=molfracs_atoms0
			endif
		enddo
		if(nrstepchem.ne.1) then
			do i=1,nmol
				if(includemol(i).or.diseqmol(i)) then
					call fillblanks(P,mixrat_r(1:nr,i),nr,dochemR,.true.)
				endif
			enddo
			call fillblanks(P,MMW,nr,dochemR,.true.)
			call fillblanks(P,nabla_ad,nr,dochemR,.true.)
		endif

		if(fixMMW) MMW=MMW0
		if(disequilibrium) then
			allocate(Otot_r(nr),Ctot_r(nr),Ntot_r(nr))
			call AddDiseqAtoms(Otot_r,Ctot_r,Ntot_r)
c call disequilibrium code
c input: 	R(1:nr+1) : These are the radial boundaries of the layers (bottom to top)
c			P(1:nr),T(1:nr) : These are the pressure and temperature inside the layers
c			molname(1:nmol) : names of the molecules included
c			Kzz_g(1:nr) : Diffusion coefficient
c input/output:	mixrat_r(1:nr,1:nmol) : number densities inside each layer. Now set to equilibrium abundances.
		   call output("==================================================================")
		   call output("Computing disequilibrium chemistry")
		   call diseq_calc(nr,R(1:nr+1),P(1:nr),T(1:nr),nmol,molname(1:nmol),mixrat_r(1:nr, 1:nmol),COratio,Kzz_g(1:nr))		   
			call CorrectDiseqAtoms(Otot_r,Ctot_r,Ntot_r)
			deallocate(Otot_r,Ctot_r,Ntot_r)
		endif
		if(nfixmol.gt.0) then
			do i=1,nfixmol
				mixrat_r(1:nr,ifixmol(i))=fixmol_abun(i)*exp(-P(1:nr)/fixmol_P(i))
			enddo
		endif
		do i=1,nr
			do imol=1,nmol
				if(P(i).lt.Pswitch_mol(imol)) mixrat_r(i,imol)=abun_switch_mol(imol)
			enddo
		enddo
		call doPhotoChemMol()
		do imol=1,nmol
			if(isotope(imol).gt.0) then
				do i=1,nr
					mixrat_r(i,isotope(imol))=mixrat_r(i,imol)/f_isotope(imol)
					mixrat_r(i,imol)=mixrat_r(i,imol)*(1d0-1d0/f_isotope(imol))
				enddo
			endif
		enddo
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
			if(nTiter.le.Cloud(i)%fixcloud.or.Cloud(i)%fixcloud.le.0) call SetupCloud(i)
		enddo
		if(cloudoverlap) then
			docloud=.false.
			do icc=2,ncc
				docloud(icc,1:nclouds)=docloud(icc-1,1:nclouds)
				i=0
10				i=i+1
				docloud(icc,i)=.not.docloud(icc,i)
				if(.not.docloud(icc,i)) goto 10
			enddo
			do icc=1,ncc
				cloudfrac(icc)=1d0
				do i=1,nclouds
					if(docloud(icc,i)) then
						cloudfrac(icc)=cloudfrac(icc)*Cloud(i)%coverage
					else
						cloudfrac(icc)=cloudfrac(icc)*(1d0-Cloud(i)%coverage)
					endif
				enddo
			enddo
		else
			docloud=.false.
			do icc=2,ncc
				docloud(icc,icc-1)=.true.
				cloudfrac(icc)=Cloud(icc-1)%coverage
			enddo
			cloudfrac(1)=1d0-sum(cloudfrac(2:ncc))
			if(cloudfrac(1).lt.0d0) then
				print*,"#####################################"
				print*," CLOUD FRACTIONS ADD UP TO > 1"
				print*," Normalizing to fully clouded planet"
				cloudfrac(2:ncc)=cloudfrac(2:ncc)/sum(cloudfrac(2:ncc))
				cloudfrac(1)=0d0
				print*," Cloud coverage used:"
				do i=1,nclouds
					print*,"  cloud layer ",trim(int2string(i,'(i3)')),": ",trim(dbl2string(cloudfrac(i+1),'(f4.2)'))
				enddo
				print*,"#####################################"
			endif
		endif

		if(ncc.eq.1) then
			docloud=.true.
			cloudfrac=1d0
		endif
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
		COret=0d0
		if(Otot.gt.0d0) COret=Ctot/Otot
		if(Tform.gt.0d0) COret=COratio
		if(Otot.gt.0d0) call output("C/O: " // dbl2string(COret,'(f8.3)'))
		if(Otot.gt.0d0.and.Htot.gt.0d0) call output("[O]: " // 
     &				dbl2string(log10(Otot/Htot)-log10(0.0004509658/0.9207539305),'(f8.3)'))
		if(Ctot.gt.0d0.and.Htot.gt.0d0) call output("[C]: " // 
     &				dbl2string(log10(Ctot/Htot)-log10(0.0002478241/0.9207539305),'(f8.3)'))

		COratio=COret
		if(nmol.ge.48) then
		if(includemol(48)) then
			do i=1,nr
				Htot=Htot+Ndens(i)*mixrat_r(i,48)
				metallicity=metallicity-Ndens(i)*mixrat_r(i,48)
			enddo
		endif
		endif
		if(Htot.gt.0d0) metallicity=log10(metallicity/Htot)+3.0565202503263760
	endif

	if(writefiles) then
		open(unit=50,file=trim(outputdir) // 'COprofile.dat',FORM="FORMATTED",ACCESS="STREAM")
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
			if(Otot.gt.0d0.and.Ctot.gt.0d0.and.Htot.gt.0d0) write(50,'(es15.4,3f10.3)') P(i),Ctot/Otot,
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

	beta_used=abs(betaT)
	if(deepredist) beta_used=betaF

	tau=0d0
	Tirr=sqrt(Rstar/(Dplanet))*Tstar
	call output("Irradiation temperature: " // dbl2string(Tirr,'(f8.1)'))
	if(grey_isoT) then
		x=(45d0**4+Tirr**4*beta_used)**0.25
		return
	endif
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
	if(deepredist) beta_used=betaF

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
	
	open(unit=50,file=trim(outputdir) // 'densityprofile' // trim(side) // '.dat',FORM="FORMATTED",ACCESS="STREAM")
	write(50,'("#",a14,a15,a15,a13,a10,a11,a10,a10)') "radius [cm]","height [cm]","dens [g/cm^3]","N [1/cm^3]","T [K]",
     & 	"P [bar]","g [cm/s^2]","MMW"
	do i=1,nr
		write(50,'(es15.7E3,es15.4E3,es15.4E3,es13.4E3,f10.3,es11.3E3,f10.3,f10.3)') sqrt(R(i)*R(i+1)),sqrt(R(i)*R(i+1))-Rplanet
     &			,dens(i),Ndens(i),T(i),P(i),grav(i),MMW(i)
	enddo
	close(unit=50)

	open(unit=50,file=trim(outputdir) // 'density' // trim(side) // '.dat',FORM="FORMATTED",ACCESS="STREAM")
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

	open(unit=50,file=trim(outputdir) // 'mixingratios' // trim(side) // '.dat',FORM="FORMATTED",ACCESS="STREAM")
	form='("#",a9,a13,' // trim(int2string(nmix,'(i2)')) // 'a15)'
	write(50,trim(form)) "T [K]","P [bar]",namemix(1:nmix)
	form='(f10.3,es13.3E3,' // trim(int2string(nmix,'(i2)')) // 'es15.4E3)'
	do i=1,nr
		write(50,trim(form)) T(i),P(i),mix(i,1:nmix)
	enddo
	close(unit=50)

	open(unit=50,file=trim(outputdir) // 'Kzz' // trim(side) // '.dat',FORM="FORMATTED",ACCESS="STREAM")
	form='("#",a12,a13)'
	write(50,trim(form)) "Kzz [cm^2/s]","P [bar]"
	form='(es13.3E3,es13.3E3)'
	do i=1,nr
		write(50,trim(form)) Kzz_g(i),P(i)
	enddo
	close(unit=50)

	do ii=1,nclouds
		open(unit=50,file=trim(outputdir) // 'clouddens' // trim(int2string(ii,'(i0.2)')) // trim(side) // '.dat',
     &					FORM="FORMATTED",ACCESS="STREAM")
		form='("#",a12,a15,a15,a15)'
		write(50,trim(form)) "P [bar]","dens [g/cm^3]","Eq. dens","gas dens"
		form='(es13.3,es15.3E3,es15.3E3,es15.3E3)'
		do i=1,nr
			write(50,trim(form)) P(i),cloud_dens(i,ii),dens(i)*XeqCloud(i,ii),dens(i)
		enddo
		close(unit=50)
	enddo

	if(mixrat_optEC0+sum(mixrat_optEC_r(1:nr)).gt.0d0) then
		open(unit=50,file=trim(outputdir) // 'optEC' // trim(side) // '.dat',FORM="FORMATTED",ACCESS="STREAM")
		form='("#",a12,a13)'
		write(50,trim(form)) "VMR C-atoms","P [bar]"
		form='(es13.3E3,es13.3E3)'
		do i=1,nr
			write(50,trim(form)) mixrat_optEC0+mixrat_optEC_r(i),P(i)
		enddo
		close(unit=50)
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
	call set_molfracs_atoms(COratio,SiOratio,NOratio,SOratio,metallicity)
	ini=.true.
	do i=1,nr
		call call_chemistry(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,condensates,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol,dosimplerainout,useEOS,x_el(i))
	enddo

	return
	end
	
	
	subroutine SetAbun()
	use GlobalSetup
	use Constants
	use AtomsModule
	IMPLICIT NONE
	real*8 Zout,abun_dust(N_atoms),tot
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
		tot=sum(molfracs_atoms_outgas(1:N_atoms))
		molfracs_atoms=molfracs_atoms_outgas/tot
	endif

	do i=1,N_atoms
		if(molfracs_atoms(i).lt.1d-40) molfracs_atoms(i)=1d-40
	enddo
	
	return
	end


	subroutine SplitGasDust(Temp,Pres,abun_total,abun_gas,abun_dust)
	use GlobalSetup
	use Constants
	use AtomsModule
	IMPLICIT NONE
	real*8 mutemp,mol_abun(nmol),abun_total(N_atoms),abun_gas(N_atoms),abun_dust(N_atoms)
	real*8 Pres,Temp,xtemp
	integer methGGchem
	
	methGGchem=0

	call call_GGchem(Temp,Pres,names_atoms,abun_total,N_atoms,molname(1:nmol),
     &			mol_abun,nmol,mutemp,.true.,abun_gas,methGGchem,xtemp)
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
	real*8 abun,Mcore,Rstart,Rend

	call SetupAtoms

	if(planetform) then

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
	Rstart=(planetform_Rend+planetform_Dmigrate)*AU
	Rend=planetform_Rend*AU
	if(Rend.le.0d0) Rend=Dplanet

	call InitFormation(Mstar,Tstar,Rstar,planetform_SolidC,planetform_Macc)
	call Formation(Mplanet,Mcore,Rstart,Rend,planetform_fdust,planetform_fplan,simAb_converge,MSimAb)
	if(.not.simAb_converge) then
		call output("WARNING: SimAb did not converge!")
		call output("WARNING: final mass of the planet:" // dbl2string(MSimAb/Mjup,'(se10.3)') // " Mjup")
	endif

	molfracs_atoms(1)=molfracs_atoms(1)*Hydrogenloss
	molfracs_atoms(2)=molfracs_atoms(2)*Hydrogenloss
	
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

	if(writefiles) then
	open(unit=50,file=trim(outputdir) // 'atomic.dat',FORM="FORMATTED",ACCESS="STREAM")
	write(50,'("SimAB abundances:")')
	write(50,'(a6,f7.3)') "C/O:",CO
	write(50,'(a6,f7.3)') "Si/O:",SiO
	write(50,'(a6,f7.3)') "N/O:",NO
	write(50,'(a6,f7.3)') "S/O:",SO
	write(50,'(a6,f7.3)') "[Z]:",Z
	write(50,'("=============")')
	do i=1,N_atoms
		write(50,'(a5,se18.6)') names_atoms(i),molfracs_atoms(i)
	enddo
	close(unit=50)
	endif

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
	molfracs_atoms(6:N_atoms)=molfracs_atoms(6:N_atoms)*scale

c	adjust S/O ratio
	molfracs_atoms(11)=molfracs_atoms(5)*SO

c	adjust metallicity
	scale=sum(molfracs_atoms(3:N_atoms))/(sum(molfracs_atoms(1:2))*Z0*10d0**Z)
	molfracs_atoms(1:2)=molfracs_atoms(1:2)*scale

	tot=sum(molfracs_atoms(1:N_atoms))
	molfracs_atoms=molfracs_atoms/tot

c	open(unit=50,file='atomic.dat',FORM="FORMATTED",ACCESS="STREAM")
c	do i=1,N_atoms
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

	call SetupAtoms

	readin=.false.
	molfracs_atoms_solar=molfracs_atoms
	molfracs_atoms=1d-50

	open(unit=43,file=filename,FORM="FORMATTED",ACCESS="STREAM")
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

	COratio=molfracs_atoms(3)/molfracs_atoms(5)
	call output("C/O: " // dbl2string(molfracs_atoms(3)/molfracs_atoms(5),'(f6.2)'))
	call output("Si/O:" // dbl2string(molfracs_atoms(9)/molfracs_atoms(5),'(f6.2)'))

	metallicity=log10((sum(molfracs_atoms(3:n_atoms))/sum(molfracs_atoms(1:n_atoms)))/
     &			(sum(molfracs_atoms_solar(3:n_atoms))/sum(molfracs_atoms_solar(1:n_atoms))))
	call output("[Z]: " // dbl2string(metallicity,'(f6.2)'))

c	open(unit=50,file='atomic.dat',FORM="FORMATTED",ACCESS="STREAM")
c	do i=1,N_atoms
c		write(50,'(a5,se18.6)') names_atoms(i),molfracs_atoms(i)
c	enddo
c	close(unit=50)

	return
	end


	subroutine call_chemistry(Tin,Pin,mol_abun,mol_names,nmol,ini,condensates,
     &		cloudspecies,Xcloud,Ncloud,nabla_ad,MMW,didcondens,includemol,rainout,useEOS,x_el)
	use AtomsModule
	use TimingModule
	IMPLICIT NONE
	integer nmol
	real*8 Tin,Pin,mol_abun(nmol),nabla_ad,x_el
	character*10 mol_names(nmol)
	logical includemol(nmol),didcondens,ini,condensates,rainout,useEOS
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

	if(useEOS) then
		call GetNablaEOS(Pin,Tg,molfracs_atoms(1)/(molfracs_atoms(1)+molfracs_atoms(2)),nabla_ad)
	endif

	mol_abun=0d0
	Xcloud=0d0
	call call_GGchem(Tg,Pin,names_atoms,molfracs_atoms,N_atoms,mol_names,mol_abun,nmol,MMW,condensates,gas_atoms,methGGchem,x_el)
c	call call_easy_chem(Tg,Pin,mol_abun,mol_names,nmol,ini,condensates,
c     &		cloudspecies,Xcloud,Ncloud,nabla_ad,MMW,didcondens,includemol)

	if(condensates.and.rainout) then
		molfracs_atoms=gas_atoms/sum(gas_atoms(1:N_atoms))
	endif

c	call readBaud(mol_abun,nmol,Pin,MMW)

c	mol_abun=1d-4
c	MMW=2.2

	do i=1,nmol
		if(.not.mol_abun(i).gt.0d0) mol_abun(i)=0d0
	enddo

	call cpu_time(time)
	timechem=timechem+time
	call system_clock(itime)
	itimechem=itimechem+itime

	return

	call call_easy_chem(Tg,Pin,mol_abun,mol_names,nmol,ini,condensates,
     &		cloudspecies,Xcloud,Ncloud,nabla_ad,MMW,didcondens,includemol)
	
	return
	end

	subroutine readBaud(mf,nm,Pin,mm)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nm,i
	real*8 mf(nm),Pin,P0,mm
	
	open(unit=20,file='dbf26/a_p_G4___.tsv',FORM="FORMATTED",ACCESS="STREAM")
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
	real*8 Ppoint0(nTpoints+nVpoints+nIRpoints),Tpoint0(nTpoints+nVpoints+nIRpoints),PrefTp0
	real*8 Tc,Gplan,TV(nr),TIR(nr),Tfree(nr)
	integer i,j,nj

	TV=0d0
	TIR=0d0
	Tfree=0d0
	
	Gplan=(Ggrav*Mplanet/(Rplanet)**2)
	PrefTp0=Gplan/(kappaT*1d6)
	
	if(nVpoints.gt.0) then
		do i=1,nVpoints
			Ppoint0(i)=betaT*tau_Vpoint(i)*Gplan/(kappaT*gammaT1*1d6)
			Tpoint0(i)=dT_Vpoint(i)
		enddo
		Tc=(Tstar**4*(Rstar/Dplanet)**2*(betaT*gammaT2))**0.25
		call MakePTstruct_dT(P,TV,nr,Ppoint0,Tpoint0,nVpoints,Tc,PrefTp0)
c		call MakePTstruct_T(P,TV,nr,Ppoint0,Tpoint0,nVpoints)
	endif
	if(nIRpoints.gt.0) then
		do i=1,nIRpoints
			Ppoint0(i)=tau_IRpoint(i)*Gplan/(kappaT*1d6)
			Tpoint0(i)=dT_IRpoint(i)
		enddo
		Tc=TeffP
		call MakePTstruct_dT(P,TIR,nr,Ppoint0,Tpoint0,nIRpoints,Tc,PrefTp0)
c		call MakePTstruct_T(P,TIR,nr,Ppoint0,Tpoint0,nIRpoints)
	endif
	if(nTpoints.gt.0) then
		if(freePT_fitT) then
			call MakePTstruct_T(P,Tfree,nr,Ppoint,dTpoint,nTpoints)
		else
			Tc=(0.5d0*TeffP**4+0.5d0*Tstar**4*(Rstar/Dplanet)**2*(betaT*gammaT1))**0.25
			call MakePTstruct_dT(P,Tfree,nr,Ppoint,dTpoint,nTpoints,Tc,PrefTpoint)
		endif
	endif

	T=(TV**4+TIR**4+Tfree**4+Tmin**4)**0.25

	if(taurexprofile) then
		call regridarray(log10(Ppoint),log10(dTpoint),nTpoints,log10(P),T,nr)
		T=10d0**T
		TV=T
		do i=1,nr
			nj=0
			Tc=0d0
			do j=1,nr
				if(abs(log10(P(i)/P(j))).le.(taurexsmooth/2d0)) then
					Tc=Tc+TV(j)
					nj=nj+1
				endif
			enddo
			if(nj.gt.0) T(i)=Tc/real(nj)
		enddo
	endif
	
	return
	end
	

	subroutine MakePTstruct_T(P,T,np,Pp_in,Tp_in,nT_in)
	IMPLICIT NONE
	integer np,i,nT,nT_in,j,i0
	real*8 P(np),T(np),Pp_in(nT_in),Tp_in(nT_in),dTp(nT_in)
	real*8 logPp(nT_in),logP(np),logT(np),logTp(nT_in)
	logical SKIP
	integer INCFD,IERR
	
	logPp=-log(Pp_in+1d-30)
	logTp=log(Tp_in)
	call sortw(logPp,logTp,nT_in)
	logP=-log(P)

	SKIP=.false.
	INCFD=1
	call DPCHIM(nT_in,logPp,logTp,dTp,INCFD)
	call DPCHFE(nT_in,logPp,logTp,dTp,INCFD,SKIP,np,logP,logT,IERR)

c	do i=1,np
c		if(logP(i).lt.logPp(1)) then
c			logT(i)=logTp(1)+dTp(1)*(logP(i)-logPp(1))
c		endif
c		if(logP(i).gt.logPp(nT_in)) then
c			logT(i)=logTp(nT_in)+dTp(nT_in)*(logP(i)-logPp(nT_in))
c		endif
c	enddo

	T(1:np)=exp(logT(1:np))

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


	subroutine doPhotoChemMol()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 nmax,nreac
	integer i,ir,imol
	
	do ir=1,nr
		if(kappaUV0.gt.0d0) then
			tauUV(ir)=exp(-kappaUV0*1d6*P(ir)/grav(ir))
		else if(tauUV(ir).lt.0d0) then
			tauUV(ir)=exp(-1d6*P(ir)/grav(ir))
		endif
	enddo

	do i=1,nPhotoReacts
		if(.not.PhotoReacts(i)%atomic) then
		do ir=1,nr
			nmax=1d200
			do imol=1,nmol
				if(includemol(imol).and.PhotoReacts(i)%react(imol).gt.0d0) then
					nreac=mixrat_r(ir,imol)/PhotoReacts(i)%react(imol)
					if(nreac.lt.nmax) nmax=nreac
				endif
			enddo
			nreac=nmax*min(1d0,PhotoReacts(i)%f_eff*tauUV(ir)**(PhotoReacts(i)%scaleKappa))
			do imol=1,nmol
				if(includemol(imol)) then
					if(PhotoReacts(i)%react(imol).gt.0d0) then
						mixrat_r(ir,imol)=mixrat_r(ir,imol)-nreac*PhotoReacts(i)%react(imol)
					endif
					if(PhotoReacts(i)%product(imol).gt.0d0) then
						mixrat_r(ir,imol)=mixrat_r(ir,imol)+nreac*PhotoReacts(i)%product(imol)
					endif
				endif
			enddo
			mixrat_optEC_r(ir)=nreac*PhotoReacts(i)%haze
		enddo
		else
		do ir=1,nr
			do imol=1,nmol
				if(includemol(imol)) then
					mixrat_r(ir,imol)=mixrat_r(ir,imol)+PhotoReacts(i)%abun(ir,imol)
				endif
			enddo
		enddo
		endif
	enddo
	do ir=1,nr
		do imol=1,nmol
			mixrat_r(ir,imol)=mixrat_r(ir,imol)*exp(-(PphotMol(imol)/P(ir))**2)
		enddo
	enddo

	return
	end

	subroutine doPhotoChemAtom(ir)
	use GlobalSetup
	use AtomsModule
	use Constants
	IMPLICIT NONE
	real*8 nmax,nreac
	integer i,ir,imol
	
	if(kappaUV0.gt.0d0) then
		tauUV(ir)=exp(-kappaUV0*1d6*P(ir)/grav(ir))
	else if(tauUV(ir).lt.0d0) then
		tauUV(ir)=exp(-1d6*P(ir)/grav(ir))
	endif
		
	mixrat_optEC_r=0d0
	do i=1,nPhotoReacts
		PhotoReacts(i)%abun(ir,1:nmol)=0d0
		if(PhotoReacts(i)%atomic) then
			nmax=1d200
			do imol=1,N_atoms
				if(PhotoReacts(i)%react(imol).gt.0d0) then
					nreac=molfracs_atoms(imol)/PhotoReacts(i)%react(imol)
					if(nreac.lt.nmax) nmax=nreac
				endif
			enddo
			nreac=nmax*min(1d0,PhotoReacts(i)%f_eff*tauUV(ir)**(PhotoReacts(i)%scaleKappa))
			do imol=1,N_atoms
				if(PhotoReacts(i)%react(imol).gt.0d0) then
					molfracs_atoms(imol)=molfracs_atoms(imol)-nreac*PhotoReacts(i)%react(imol)
				endif
			enddo
			do imol=1,N_atoms
				if(includemol(imol)) then
					if(PhotoReacts(i)%product(imol).gt.0d0) then
						PhotoReacts(i)%abun(ir,imol)=PhotoReacts(i)%abun(ir,imol)+nreac*PhotoReacts(i)%product(imol)
					endif
				endif
			enddo
			mixrat_optEC_r(ir)=nreac*PhotoReacts(i)%haze
		endif
	enddo

	return
	end


	real*8 function ComputeKzz(x,Tg,rhog,H,addmicro)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 x,Tg,rhog,lmfp,vth,sigmamol,Kmax,Kmin,Kp,H,Te
	integer i
	logical addmicro
	sigmamol=8d-15

	if(SCKzz) then
c use parametrization from Moses et al. (2022)
		Te=(TeffP**4+(Rstar/(Dplanet))**2*Tstar**4)**0.25
		ComputeKzz=(Kzz/sqrt(x))*(H/620d5)*(Te/1450d0)**4
c put some limits on the Kzz
		ComputeKzz=ComputeKzz+Kzz_offset
		ComputeKzz=1d0/(1d0/ComputeKzz+1d0/Kzz_max)
	else if(Kzz_deep.gt.0d0.and.Kzz_1bar.gt.0d0) then
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

	subroutine doWaterWorld()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 c0,c1,c2,c3,c4,Pm,PH2Omax,mutot,mixrat_tot,Pold(nr+1),Told(nr),fact
	integer i
	logical liquid

	c0=2.98605d+01
	c1=-3.15220d+03
	c2=-7.30370d+00
	c3=2.42470d-09
	c4=1.80900d-06

	PH2Omax=fH2O*Ggrav*Mplanet**2/(4d0*pi*Rplanet**4*1d6)

	if(Tsurface.gt.647.096) then
		call output("Warning! Waterworld ocean reaches critical temperature!")
		print*,"Warning! Waterworld ocean reaches critical temperature!"
		Pmax=PH2Omax
		liquid=.true.
	else
		call PvapH2O(Tsurface,Pmax,liquid)
		Pmax=f_water*Pmax
		if(Pmax.gt.PH2Omax) Pmax=PH2Omax
	endif
	if(Pmax.le.1d-20) Pmax=1d-20

	mutot=0d0
	mixrat_tot=0d0
	do i=1,nmol
		if(includemol(i)) then
			mutot=mutot+mixrat(i)*Mmol(i)
			mixrat_tot=mixrat_tot+mixrat(i)
		endif
	enddo
	mutot=mutot/mixrat_tot
	
	if(setsurfpressure) then
		fact=Pmax/(mixrat(1))
c		fact=fact**2
		mixrat(1)=mixrat(1)*sqrt(fact)
		Pmax=0d0
		do i=1,nmol
			Pmax=Pmax+mixrat(i)
		enddo
	else
		mixrat(1:nmol)=mixrat(1:nmol)/mixrat_tot
		fact=Pmax/(P(1)*mixrat(1))
c		fact=fact**2
		Pmax=P(1)*sqrt(fact)
	endif
	Pm=Pmin
	if(Pm.gt.Pmax/100d0) Pm=Pmax/100d0
	Pplanet=Pmax

	call output("Surface pressure: " // dbl2string(Pmax,'(es8.2)'))
	if(.not.liquid) then
		call output("Ice world!")
		print*,"Ice world!"
	endif

	Pold=P
	Told=T
	do i=1,nr
		P(i)=10d0**(log10(Pmax)+(log10(Pm/Pmax)*real(i-1)/real(nr-1)))
	enddo
	return
	call regridarray(-log(Pold(1:nr)),Told,nr,-log(P(1:nr)),T,nr)
	do i=1,nr
		if(P(i).gt.Pold(1)) then
			T(i)=Told(1)/exp(log(Pold(1)/P(i))*(log(Told(1)/Told(2))/log(Pold(1)/Pold(2))))
		endif
	enddo

	return
	end
	