	subroutine SetupOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE

	call ReadOpacities()

	return
	end


	subroutine ReadOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer imol
	real*8 nu1,nu2,tanscale,ll,tot,tot2
	real*8 x1,x2,rr,gasdev,random,dnu,Saver,starttime,stoptime,cwg(ng),w1
	real*8,allocatable :: nu_line(:),dnu_line(:),mixrat_tmp(:)
	real*8,allocatable :: opac_tot(:,:),cont_tot(:),kaver(:),kappa_mol(:,:,:)
	logical,allocatable,save :: fulladd(:)
	real*8,allocatable,save :: k_line(:),ktemp(:),kappa(:),w_line(:),kappa_tot(:)
!$OMP THREADPRIVATE(fulladd,k_line,ktemp,kappa,w_line,kappa_tot)
	integer n_nu_line,iT
	integer i,j,ir,k,nl,ig,ig_c,imol0
	integer,allocatable :: inu1(:),inu2(:)
	character*500 filename

	allocate(cont_tot(nlam))
	allocate(kaver(nlam))
	allocate(opac_tot(nlam,ng))
	allocate(kappa_mol(ng,nlam,nmol))
	allocate(mixrat_tmp(nmol))

	j=0
	do i=1,nmol
		if(opacitymol(i)) j=j+1
	enddo
	n_nu_line=ng*ng
c	n_nu_line=ng*min(j,4)
	if(.not.emisspec.and..not.computeT) n_nu_line=ng
	
	allocate(nu_line(n_nu_line))

	cwg(1)=wgg(1)
	do ig=2,ng
		cwg(ig)=wgg(ig)+cwg(ig-1)
	enddo
	cwg(1:ng)=cwg(1:ng)/cwg(ng)

	if(.not.allocated(ig_comp)) then
		ng_comp=1d6
		allocate(ig_comp(ng_comp))
		do i=1,ng_comp
			rr=random(idum)
			call hunt(cwg,ng,rr,ig)
			ig_comp(i)=ig+1
		enddo
	endif
	
	opac_tot=0d0

	nu1=0d0
	nu2=1d0
	do i=1,n_nu_line
		nu_line(i)=1d0-real(i-1)/real(n_nu_line-1)
	enddo

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(n_nu_line,ng,nmol)
		allocate(k_line(n_nu_line))
		allocate(ktemp(ng))
		allocate(kappa(ng))
		allocate(w_line(n_nu_line))
		allocate(fulladd(nmol))
		allocate(kappa_tot(0:nmol))
!$OMP END PARALLEL

	do ir=nr,1,-1
		call tellertje(nr-ir+1,nr)
		cont_tot=0d0

c===============
c UV cross sections of CO2 from Venot et al.
		call CO2_UV_cross(lam,cont_tot,nlam,min(T(ir),800d0))
		cont_tot=cont_tot*mixrat_r(ir,2)
		do i=1,nlam
			if(lam(i).gt.0.3d-4) cont_tot(i)=0d0
		enddo
c===============

		if(P(ir).gt.psimplecloud) then
			cont_tot(1:nlam)=1d0/Ndens(ir)
		endif
		mixrat_tmp(1:nmol)=mixrat_r(ir,1:nmol)
		do i=1,ncia
			if(T(ir).lt.CIA(i)%T(1)) then
				iT=1
			else if(T(ir).gt.CIA(i)%T(CIA(i)%nT)) then
				iT=CIA(i)%nT-1
			else
				do iT=1,CIA(i)%nT-1
					if(T(ir).ge.CIA(i)%T(iT).and.T(ir).le.CIA(i)%T(iT+1)) exit
				enddo
			endif
			cont_tot(1:nlam)=cont_tot(1:nlam)+CIA(i)%Cabs(iT,1:nlam)*Ndens(ir)*mixrat_tmp(CIA(i)%imol1)*mixrat_tmp(CIA(i)%imol2)
		enddo
		kappa_mol=0d0
		do imol=1,nmol
			if(includemol(imol)) then
				call ReadOpacityFITS(kappa_mol,imol,ir)
			endif
		enddo
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,j,imol,ig,ig_c,tot,tot2,imol0,w1)
!$OMP& SHARED(nlam,n_nu_line,nmol,mixrat_tmp,ng,ir,kappa_mol,cont_tot,Cabs,Csca,opac_tot,Ndens,R,computelam,
!$OMP&        ig_comp,retrieval,domakeai,gg,wgg,ng_comp,opacitymol,emisspec,computeT,lamemis,useobsgrid,
!$OMP&        RTgridpoint)
!$OMP DO SCHEDULE(DYNAMIC,1)
		do i=1,nlam
			if(computelam(i).and.(emisspec.or.computeT).and.(.not.useobsgrid.or.lamemis(i).or.RTgridpoint(i))) then
			kappa_tot(0:nmol)=0d0
			kappa_tot(0)=cont_tot(i)
			do imol=1,nmol
				if(opacitymol(imol)) then
					do ig=1,ng
						kappa_tot(imol)=kappa_tot(imol)+wgg(ig)*kappa_mol(ig,i,imol)*mixrat_tmp(imol)
					enddo
					kappa_tot(0)=kappa_tot(0)+kappa_tot(imol)
				endif
			enddo
			do imol=1,nmol
				if(opacitymol(imol)) then
					fulladd(imol)=(kappa_tot(imol).gt.0.01*kappa_tot(0))
				endif
			enddo
			tot=0d0
			kappa(1:ng)=0d0
			do imol=1,nmol
				if(opacitymol(imol).and.fulladd(imol)) then
					kappa(1:ng)=kappa_mol(1:ng,i,imol)*mixrat_tmp(imol)
					tot=tot+sum(kappa(1:ng)*wgg(1:ng))
					exit
				endif
			enddo
			imol0=imol
			if(imol0.lt.nmol) then
			do imol=1,nmol
				if(opacitymol(imol).and.fulladd(imol).and.imol.ne.imol0) then
					ktemp(1:ng)=kappa_mol(1:ng,i,imol)*mixrat_tmp(imol)
					ig_c=0
					tot=tot+sum(ktemp(1:ng)*wgg(1:ng))
					do ig=1,ng
						do j=1,ng
							ig_c=ig_c+1
							k_line(ig_c)=kappa(ig)+ktemp(j)
							w_line(ig_c)=wgg(ig)*wgg(j)
						enddo
					enddo
					n_nu_line=ig_c
					do ig=1,n_nu_line
						if(.not.k_line(ig).ge.1d-80) k_line(ig)=1d-80
					enddo
					call sortw(k_line,w_line,n_nu_line)
					do ig=2,n_nu_line
						w_line(ig)=w_line(ig)+w_line(ig-1)
					enddo
					w_line(1:n_nu_line)=w_line(1:n_nu_line)/w_line(n_nu_line)
					do ig=1,ng
						call hunt(w_line,n_nu_line,gg(ig),ig_c)
						if(ig_c.le.1) then
							kappa(ig)=k_line(1)
						else
							w1=(gg(ig)-w_line(ig_c+1))/(w_line(ig_c)-w_line(ig_c+1))
							kappa(ig)=k_line(ig_c)*w1+k_line(ig_c+1)*(1d0-w1)
						endif
					enddo
				endif
			enddo
			endif
			do imol=1,nmol
				if(opacitymol(imol).and..not.fulladd(imol)) then
					ktemp(1:ng)=kappa_mol(1:ng,i,imol)
					tot=tot+sum(ktemp(1:ng)*wgg(1:ng))*mixrat_tmp(imol)
					kappa(1:ng)=kappa(1:ng)+ktemp(1:ng)*mixrat_tmp(imol)
				endif
			enddo
			else
			kappa(1:ng)=0d0
			tot=0d0
			do imol=1,nmol
				if(opacitymol(imol)) then
					ktemp(1:ng)=kappa_mol(1:ng,i,imol)
					tot=tot+sum(ktemp(1:ng)*wgg(1:ng))*mixrat_tmp(imol)
					kappa(1:ng)=kappa(1:ng)+ktemp(1:ng)*mixrat_tmp(imol)
				endif
			enddo
			endif

			tot2=sum(kappa(1:ng)*wgg(1:ng))+cont_tot(i)
			tot=tot+cont_tot(i)
			kappa=kappa+cont_tot(i)
			if(tot2.gt.0d0) then
				kappa=kappa*tot/tot2
			endif
			Cabs(ir,i,1:ng)=kappa(1:ng)
			call RayleighScattering(Csca(ir,i),ir,i)
			do ig=1,ng
				if(Cabs(ir,i,ig).lt.1d-6*Csca(ir,i)) Cabs(ir,i,ig)=1d-6*Csca(ir,i)
			enddo
			opac_tot(i,1:ng)=opac_tot(i,1:ng)+(Cabs(ir,i,1:ng)+Csca(ir,i))*Ndens(ir)*(R(ir+1)-R(ir))
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		Cext_cont(ir,1:nlam)=(cont_tot(1:nlam)+Csca(ir,1:nlam))*Ndens(ir)
		do imol=1,nmol
			if(includemol(imol)) then
				do i=1,nlam
					do j=1,ng
						Cabs_mol(j,i,imol,ir)=kappa_mol(j,i,imol)*Ndens(ir)*mixrat_r(ir,imol)
					enddo
				enddo
			endif
		enddo
		if(outputopacity) then
			call WriteOpacity(ir,"ktab",freq,Cabs(ir,1:nlam,1:ng),nlam,ng)
			do i=1,nlam
				kaver(i)=0d0
				do j=1,ng
					kaver(i)=kaver(i)+wgg(j)*Cabs(ir,i,j)
				enddo
			enddo
			call WriteOpacity(ir,"aver",freq,kaver(1:nlam),nlam,1)
			call WriteOpacity(ir,"scat",freq,Csca(ir,1:nlam)*Ndens(ir)/dens(ir),nlam,1)
		endif
	enddo

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
		deallocate(k_line)
		deallocate(w_line)
		deallocate(ktemp)
		deallocate(kappa,fulladd,kappa_tot)
!$OMP END PARALLEL

	if(.not.retrieval) then
		open(unit=30,file=trim(outputdir) // "opticaldepth.dat",RECL=6000)
		write(30,'("#",a13,a19)') "lambda [mu]","total average tau"
		do i=1,nlam
			write(30,'(f12.6,e19.7)') lam(i)/micron,sum(opac_tot(i,1:ng)*wgg(1:ng))
		enddo
		close(unit=30)
	endif
	
	deallocate(cont_tot)
	deallocate(kaver)
	deallocate(opac_tot)
	deallocate(kappa_mol)
	deallocate(nu_line)
	deallocate(mixrat_tmp)

	return
	end



	subroutine RayleighScattering(Cs,ir,i)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 ll,Cs
	integer ir,i,j

	ll=1d0/lam(i)**2
c Rayleigh cross sections from Dalgarno & Williams (1962)
c For other than H2 from Sneep & Ubachs (2005)
	Cs=0d0
	do j=1,nmol
		if(mixrat_r(ir,j).gt.0d0) then
			select case(j)
				case(2) !CO2
c Sneep & Ubachs (2005)
					Cs=Cs+12.4*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
				case(4) !N2O
c Sneep & Ubachs (2005)
					Cs=Cs+15.9*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
				case(5) !CO
c Sneep & Ubachs (2005)
					Cs=Cs+6.19*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
				case(6) !CH4
c Sneep & Ubachs (2005)
					Cs=Cs+12.47*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
				case(22) !N2
c Sneep & Ubachs (2005)
					Cs=Cs+5.1*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
				case(30) !SF6
c Sneep & Ubachs (2005)
					Cs=Cs+32.3*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
				case(45) ! H2
c Dalgarno & Williams (1962)
					Cs=Cs+mixrat_r(ir,j)*(8.14d-45*ll**2+1.28d-54*ll**3+1.61d-64*ll**4)
				case(1) !H2O
					Cs=Cs+4.43*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
			end select
		endif
	enddo
c============================
c old parameterization
	Cs=Cs+mixratHaze*(8.14d-45*ll**2+1.28d-54*ll**3+1.61d-64*ll**4+kappaHaze*2d0*mp)*exp(-(abs(log10(P(ir)/PHaze))
     &					/log10(dPHaze))**2/2d0)
c============================

c============================
c Pinhas 2019
c	Cs=Cs+mixratHaze*5.31e-27*(sqrt(ll)*0.35e-4)**4
c============================

	return
	end


	subroutine ReadOpacityFITS(kappa_mol,imol,ir)
	use GlobalSetup
	use OpacityFITSdata
	implicit none
	character*80 comment,errmessage
	character*30 errtext
	integer ig,ilam,iT,iP,imol,i,j,ir,ngF,i1,i2
	real*8 kappa_mol(ng,nlam,nmol),wP1,wP2,wT1,wT2,x1,x2,tot,tot2,random,w1,ww
	real*8,allocatable :: temp(:),wtemp(:)
	real*8,dimension(:,:),pointer :: tab1,tab2,tab3,tab4
	type(databaseKtable),pointer :: Ktab

	Ktab => Ktable(imol)
	if(.not.Ktab%available.or..not.includemol(imol)) then
		kappa_mol(1:ng,1:nlam,imol)=0d0
		return
	endif

	if(P(ir).lt.Ktab%P1) then
		iP=1
		wP1=1d0
		wP2=0d0
	else if(P(ir).gt.Ktab%P2) then
		iP=Ktab%nP-1
		wP1=0d0
		wP2=1d0
	else
		call hunt(Ktab%P,Ktab%nP,P(ir),iP)
c		do iP=1,Ktab%nP
c			if(P(ir).ge.Ktab%P(iP).and.P(ir).lt.Ktab%P(iP+1)) exit
c		enddo
		wP1=1d0-log10(P(ir)/Ktab%P(iP))/log10(Ktab%P(iP+1)/Ktab%P(iP))
		wP2=1d0-wP1
	endif

	if(T(ir).lt.Ktab%T1) then
		iT=1
		wT1=1d0
		wT2=0d0
	else if(T(ir).gt.Ktab%T2) then
		iT=Ktab%nT-1
		wT1=0d0
		wT2=1d0
	else
		call hunt(Ktab%T,Ktab%nT,T(ir),iT)
c		do iT=1,Ktab%nT
c			if(T(ir).ge.Ktab%T(iT).and.T(ir).lt.Ktab%T(iT+1)) exit
c		enddo
		wT1=1d0-log10(T(ir)/Ktab%T(iT))/log10(Ktab%T(iT+1)/Ktab%T(iT))
		wT2=1d0-wT1
	endif
 
 	tab1 => Ktab%ktable(1:ng,1:nlam,iT,iP)
 	tab2 => Ktab%ktable(1:ng,1:nlam,iT+1,iP)
 	tab3 => Ktab%ktable(1:ng,1:nlam,iT,iP+1)
 	tab4 => Ktab%ktable(1:ng,1:nlam,iT+1,iP+1)
 
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ilam,ig)
!$OMP& SHARED(nlam,ng,kappa_mol,imol,tab1,tab2,tab3,tab4,wT1,wT2,wP1,wP2,iT,iP,computelam)
!$OMP DO
	do ilam=1,nlam
		if(computelam(ilam)) then
			kappa_mol(1:ng,ilam,imol)=tab1(1:ng,ilam)*wT1*wP1+tab2(1:ng,ilam)*wT2*wP1+tab3(1:ng,ilam)*wT1*wP2+tab4(1:ng,ilam)*wT2*wP2
		endif
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	return
	end


	
	subroutine CO2_UV_cross(lam,C,nlam,T)
c Venot et al. 2018
	IMPLICIT NONE
	integer nlam,i
	real*8 lam(nlam),C(nlam),T
	real*8 s1,s2,s3,A1,A2,A3,nu,nu1,nu2,nu3
	
	s1=877.36+10947.81*exp(-1382.63/T)
	s2=T*(2.78+49.52*exp(-0.00654*T))
	s3=T*(8.17+46.12*exp(-0.00813*T))
	
	A1=50d-19
	A2=(3.58+9.18*exp(-580.92/T))*1d-19
	A3=(4.09+0.0022*T)*1d-19

	nu1=88574.0
	nu2=76000.0
	nu3=68000.0
	
	do i=1,nlam
		nu=1d0/lam(i)
		C(i)=A1*exp(-(nu-nu1)**2/(2d0*s1**2))
     &		+A2*exp(-(nu-nu2)**2/(2d0*s2**2))
     &		+A3*exp(-(nu-nu3)**2/(2d0*s3**2))
	enddo
	
	return
	end
		
