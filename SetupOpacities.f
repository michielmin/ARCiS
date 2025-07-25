	subroutine SetupOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE

	call ComputeSurface()
	call ReadOpacities()

	return
	end


	subroutine ReadOpacities()
	use GlobalSetup
	use Constants
	use OpacityFITSdata
	IMPLICIT NONE
	integer imol
	real*8 nu1,nu2,tanscale,ll,tot,tot2
	real*8 x1,x2,rr,gasdev,random,dnu,Saver,starttime,stoptime,cwg(ng),w1
	real*8,allocatable,save :: nu_line(:),mixrat_tmp(:),kappa_omp(:,:)
	real*8,allocatable,save :: opac_tot(:,:),cont_tot(:),kaver(:),kappa_mol(:,:,:)
	integer,allocatable,save :: ifull(:),ifast(:)
	real*8,allocatable,save :: k_line(:),ktemp(:),kappa(:),w_line(:),kappa_tot(:),work1(:),work2(:),work3(:)
!$OMP THREADPRIVATE(ifull,ifast,k_line,ktemp,w_line,kappa_tot,work1,work2,work3,kappa)
	integer n_nu_line,iT,nfull,nfast,ivel
	integer i,j,ir,k,nl,ig,ig_c,imol0,jg,nmap,imap(nmol)
	character*500 filename
	logical,save :: first_entry=.true.
	real*8,allocatable :: Cabs_optEC(:),Csca_optEC(:)
	logical do_optEC
	
	do_optEC=.false.
	if(mixrat_optEC0.gt.0d0) then
		do_optEC=.true.
	else
		do ir=1,nr
			if(mixrat_optEC_r(ir).gt.0d0) then
				do_optEC=.true.
				exit
			endif
		enddo
	endif
	if(do_optEC) then
		allocate(Cabs_optEC(nlam),Csca_optEC(nlam))
		Cabs_optEC=0d0
		Csca_optEC=0d0
		mixrat_optEC=1d0
		call Compute_optEC(Cabs_optEC,Csca_optEC,computelam)
	endif

	j=0
	do i=1,nmol
		if(opacitymol(i)) j=j+1
	enddo
	n_nu_line=ng*ng
c	n_nu_line=ng*min(j,4)
	if(.not.emisspec.and..not.computeT.and..not.doRing) n_nu_line=ng
	
	if(first_entry) then
		allocate(cont_tot(nlam))
		allocate(kaver(nlam))
		allocate(opac_tot(nlam,ng))
		allocate(kappa_mol(ng,nlam,nmol))
		allocate(mixrat_tmp(nmol))
		allocate(nu_line(n_nu_line))
!$OMP PARALLEL IF(useomp)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(n_nu_line,ng,nmol)
		allocate(kappa(ng))
		allocate(k_line(n_nu_line))
		allocate(ktemp(ng))
		allocate(w_line(n_nu_line))
		allocate(ifull(nmol))
		allocate(ifast(nmol))
		allocate(kappa_tot(0:nmol))
		allocate(work1(n_nu_line))
		allocate(work2(n_nu_line+1))
		allocate(work3(n_nu_line))
!$OMP END PARALLEL
		first_entry=.false.
	endif

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

	do ir=nr,1,-1
		call tellertje(nr-ir+1,nr)
		do i=1,nlam
			cont_tot(i)=0d0
		enddo

c===============
c UV cross sections of CO2 from Venot et al.
		if(includemol(2)) then
			call CO2_UV_cross(lam,cont_tot,nlam,min(T(ir),800d0))
			cont_tot=cont_tot*mixrat_r(ir,2)
			do i=1,nlam
				if(lam(i).gt.0.3d-4) cont_tot(i)=0d0
			enddo
		endif
c===============

		mixrat_tmp(1:nmol)=mixrat_r(ir,1:nmol)
		do i=1,ncia
			if(includemol(CIA(i)%imol1).and.includemol(CIA(i)%imol2)) then
				if(T(ir).lt.CIA(i)%T(1)) then
					iT=1
				else if(T(ir).gt.CIA(i)%T(CIA(i)%nT)) then
					iT=CIA(i)%nT
				else
					do iT=1,CIA(i)%nT-1
						if(T(ir).ge.CIA(i)%T(iT).and.T(ir).le.CIA(i)%T(iT+1)) exit
					enddo
				endif
				if(T(ir).lt.CIA(i)%T(1)) then
					cont_tot(1:nlam)=cont_tot(1:nlam)+CIA(i)%Cabs(1,1:nlam)*Ndens(ir)*
     &						mixrat_tmp(CIA(i)%imol1)*mixrat_tmp(CIA(i)%imol2)
				else if(iT.lt.CIA(i)%nT) then
					w1=(CIA(i)%T(iT+1)-T(ir))/(CIA(i)%T(iT+1)-CIA(i)%T(iT))
					cont_tot(1:nlam)=cont_tot(1:nlam)+(w1*CIA(i)%Cabs(iT,1:nlam)+(1d0-w1)*CIA(i)%Cabs(iT+1,1:nlam))*
     &									Ndens(ir)*mixrat_tmp(CIA(i)%imol1)*mixrat_tmp(CIA(i)%imol2)
				else
					cont_tot(1:nlam)=cont_tot(1:nlam)+CIA(i)%Cabs(CIA(i)%nT,1:nlam)*Ndens(ir)*
     &						mixrat_tmp(CIA(i)%imol1)*mixrat_tmp(CIA(i)%imol2)
				endif
			endif
		enddo
		Csca(ir,1:nlam)=0d0
		
		nmap=0
		do imol=1,nmol
			if(opacitymol(imol)) then
				nmap=nmap+1
				imap(nmap)=imol
			endif
		enddo

		if(do_rayleigh) then
			do i=1,nlam
				if(computelam(i)) call RayleighScattering(Csca(ir,i),ir,i)
			enddo
		endif

		if(mixrat_PAH.gt.0d0) call ComputePAH(cont_tot,Csca(ir,1:nlam),computelam)
		mixrat_optEC=mixrat_optEC0+mixrat_optEC_r(ir)
		if(do_optEC) then
			cont_tot=cont_tot+Cabs_optEC*mixrat_optEC
			Csca(ir,1:nlam)=Csca(ir,1:nlam)+Csca_optEC(1:nlam)*mixrat_optEC
		endif

		do ivel=-nvel,nvel
		do i=1,nlam
			kappa_mol(1:ng,i,1:nmol)=0d0
		enddo

		do imol=1,nmol
			if(includemol(imol)) then
				call ReadOpacityFITS(kappa_mol,imol,ir,ivel)
			endif
		enddo

!$OMP PARALLEL IF(useomp)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,j,imol,ig,jg,ig_c,imol0,w1,nfull,nfast)
!$OMP& SHARED(nlam,n_nu_line,nmol,mixrat_tmp,ng,ir,kappa_mol,cont_tot,Cabs,Csca,opac_tot,Ndens,R,computelam,
!$OMP&        ig_comp,retrieval,domakeai,gg,wgg,ng_comp,opacitymol,emisspec,computeT,doRing,lamemis,useobsgrid,
!$OMP&        RTgridpoint,includemol,do_rayleigh,mixrat_PAH,mixrat_optEC0,mixrat_optEC,mixrat_optEC_r,nmap,imap,
!$OMP&        do_optEC,Cabs_optEC,Csca_optEC,ivel)
!$OMP DO SCHEDULE(DYNAMIC)
		do i=1,nlam
			if(computelam(i).and.(emisspec.or.computeT.or.doRing).and.(.not.useobsgrid.or.lamemis(i).or.RTgridpoint(i))) then
			kappa_tot(0)=cont_tot(i)
			do imol=1,nmap
				kappa_tot(imol)=dot_product(wgg(1:ng),kappa_mol(1:ng,i,imap(imol)))
				kappa_tot(imol)=kappa_tot(imol)*mixrat_tmp(imap(imol))
				kappa_tot(0)=kappa_tot(0)+kappa_tot(imol)
			enddo
			nfull=0
			nfast=0d0
			do imol=1,nmap
				if(kappa_tot(imol).ge.0.01*kappa_tot(0)) then
					nfull=nfull+1
					ifull(nfull)=imap(imol)
					kappa_tot(nfull)=-kappa_tot(imol)
				else
					nfast=nfast+1
					ifast(nfast)=imap(imol)
				endif
			enddo
			call dpquicksort_indx(kappa_tot(1:nfull),ifull(1:nfull),nfull)
c			call sortidx_2(kappa_tot(1:nfull),ifull(1:nfull),nfull)
			kappa(1:ng)=0d0
			if(nfull.gt.0) then
				imol=ifull(1)
				kappa(1:ng)=kappa_mol(1:ng,i,imol)*mixrat_tmp(imol)
				do j=2,nfull
					imol=ifull(j)
					ktemp(1:ng)=kappa_mol(1:ng,i,imol)*mixrat_tmp(imol)
					ig_c=0
					do ig=1,ng
						do jg=1,ng
							ig_c=ig_c+1
							k_line(ig_c)=kappa(ig)+ktemp(jg)
							w_line(ig_c)=wgg(ig)*wgg(jg)
						enddo
					enddo
					n_nu_line=ig_c
					do ig=1,n_nu_line
						if(.not.k_line(ig).ge.1d-80) k_line(ig)=1d-80
					enddo
					call regridKtable(k_line,w_line,n_nu_line,gg,kappa(1:ng),wgg,ng,work1,work2,work3)
				enddo
			endif
			do j=1,nfast
				imol=ifast(j)
				kappa(1:ng)=kappa(1:ng)+kappa_mol(1:ng,i,imol)*mixrat_tmp(imol)
			enddo
			else
			kappa(1:ng)=0d0
			do imol=1,nmol
				if(opacitymol(imol)) then
					kappa(1:ng)=kappa(1:ng)+kappa_mol(1:ng,i,imol)*mixrat_tmp(imol)
				endif
			enddo
			endif

			kappa(1:ng)=kappa(1:ng)+cont_tot(i)
			Cabs(ir,i,1:ng,ivel)=kappa(1:ng)
			do ig=1,ng
				if(Cabs(ir,i,ig,ivel).lt.1d-6*Csca(ir,i)) Cabs(ir,i,ig,ivel)=1d-6*Csca(ir,i)
			enddo
			if(.not.retrieval.and.ivel.eq.0) opac_tot(i,1:ng)=opac_tot(i,1:ng)+(Cabs(ir,i,1:ng,ivel)+Csca(ir,i))*Ndens(ir)*(R(ir+1)-R(ir))
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		if(ivel.eq.0) then
		do imol=1,nmol
			if(includemol(imol)) then
				do i=1,nlam
					do j=1,ng
						Cabs_mol(j,i,imol,ir)=kappa_mol(j,i,imol)*Ndens(ir)*mixrat_r(ir,imol)
					enddo
				enddo
			endif
		enddo
		endif
		enddo
		Cext_cont(ir,1:nlam)=(cont_tot(1:nlam)+Csca(ir,1:nlam))*Ndens(ir)
		if(outputopacity) then
			call WriteOpacity(ir,"ktab",freq,Cabs(ir,1:nlam,1:ng,0),nlam,ng)
			do i=1,nlam
				kaver(i)=0d0
				do j=1,ng
					kaver(i)=kaver(i)+wgg(j)*Cabs(ir,i,j,0)
				enddo
			enddo
			call WriteOpacity(ir,"aver",freq,kaver(1:nlam),nlam,1)
			call WriteOpacity(ir,"scat",freq,Csca(ir,1:nlam)*Ndens(ir)/dens(ir),nlam,1)
		endif
	enddo

	if(writefiles) then
		open(unit=30,file=trim(outputdir) // "opticaldepth.dat",FORM="FORMATTED",ACCESS="STREAM")
		write(30,'("#",a13,a19)') "lambda [mu]","total average tau"
		do i=1,nlam
			if(computelam(i)) write(30,'(f12.6,e19.7)') lam(i)/micron,sum(opac_tot(i,1:ng)*wgg(1:ng))
		enddo
		close(unit=30)
	endif
	
	if(do_optEC) deallocate(Cabs_optEC,Csca_optEC)


	return
	end



	subroutine RayleighScattering(Cs,ir,i)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 ll,Cs,ll2,RayleighHatom
	integer ir,i,j

	ll=1d0/lam(i)**2
	ll2=ll*ll
c Rayleigh cross sections from Dalgarno & Williams (1962)
c For other than H2 from Sneep & Ubachs (2005)
	Cs=0d0
	do j=1,nmol
		if(mixrat_r(ir,j).gt.0d0) then
			select case(j)
				case(2) !CO2
c Sneep & Ubachs (2005)
					Cs=Cs+28.499d-45*mixrat_r(ir,j)*ll2**(4.1343/4.0)
				case(4) !N2O
c Sneep & Ubachs (2005)
					Cs=Cs+15.90e-27*mixrat_r(ir,j)*(ll2/18788.4**4)
				case(5) !CO
c Sneep & Ubachs (2005)
					Cs=Cs+6.19e-27*mixrat_r(ir,j)*(ll2/18788.4**4)
				case(6) !CH4
c Sneep & Ubachs (2005)
					Cs=Cs+12.47e-27*mixrat_r(ir,j)*(ll2/18788.4**4)
				case(22) !N2
c Sneep & Ubachs (2005)
					Cs=Cs+16.31d-45*mixrat_r(ir,j)*ll2**(4.0974/4.0)
				case(30) !SF6
c Sneep & Ubachs (2005)
					Cs=Cs+32.3e-27*mixrat_r(ir,j)*(ll2/18788.4**4)
				case(45) ! H2
c Dalgarno & Williams (1962)
					Cs=Cs+mixrat_r(ir,j)*(8.14d-45*ll2+1.28d-54*ll2*ll+1.61d-64*ll2*ll2)
				case(1) !H2O
					Cs=Cs+3.555040361134916e-44*mixrat_r(ir,j)*ll2
				case(48) ! He
c https://www.climate-policy-watcher.org/surface-temperature/scattering-by-molecules-rayleigh-scattering.html
					Cs=Cs+0.0641*mixrat_r(ir,j)*(8.14d-45*ll2+1.28d-54*ll2*ll+1.61d-64*ll2*ll2)
				case(64) ! H
c from Helios source code
					Cs=Cs+mixrat_r(ir,j)*RayleighHatom(lam(i))
				case(11) ! NH3
c https://www.climate-policy-watcher.org/surface-temperature/scattering-by-molecules-rayleigh-scattering.html
					Cs=Cs+7.3427*mixrat_r(ir,j)*(8.14d-45*ll2+1.28d-54*ll2*ll+1.61d-64*ll2*ll2)
				case(13) ! OH
c Tarafdar and Vardya 1973
					Cs=Cs+(12.41**2)*mixrat_r(ir,j)*13056.839884384113d-50*ll2
				case(54) ! SiO
c Tarafdar and Vardya 1973
					Cs=Cs+(26.41**2)*mixrat_r(ir,j)*13056.839884384113d-50*ll2
				case(23) ! HCN
c Tarafdar and Vardya 1973
					Cs=Cs+(25.54**2)*mixrat_r(ir,j)*13056.839884384113d-50*ll2
				case(31) ! H2S
c Tarafdar and Vardya 1973
					Cs=Cs+(32.14**2)*mixrat_r(ir,j)*13056.839884384113d-50*ll2
				case(38) ! C2H4
c Tarafdar and Vardya 1973
					Cs=Cs+(42.61**2)*mixrat_r(ir,j)*13056.839884384113d-50*ll2
c				case(9) ! SO2 (but taken the value from CO2 as a poor-man's solution)
c Sneep & Ubachs (2005)
c					Cs=Cs+28.499d-47*mixrat_r(ir,j)*ll2**(4.1343/4.0)
			end select
		endif
	enddo

	return
	end
	
	real*8 function RayleighHatom(lam)
	IMPLICIT NONE
	real*8 lam,Cp(10),sigma_T,lam_l,sum_term
	parameter(Cp = (/1.26563d0, 3.73828125d0, 8.813930935d0, 19.15379502d0, 
     &		39.92303232d0, 81.10881152d0, 161.9089166d0, 319.0231631d0, 622.2679809d0, 1203.891509d0/))	
	parameter(sigma_T = 0.665d-24)
	parameter(lam_l = 91.2d-7)
	integer i

	sum_term = 0d0
	do i=1,10
		sum_term=sum_term+cp(i)*(lam_l/lam)**(2d0*real(i))
	enddo
	RayleighHatom = sigma_T * (lam_l / lam)**4 * sum_term

	return
	end

	subroutine ReadOpacityFITS(kappa_mol,imol,ir,ivel)
	use GlobalSetup
	use OpacityFITSdata
	implicit none
	character*80 comment,errmessage
	character*30 errtext
	integer ig,ilam,iT,iP,imol,i,j,ir,ngF,i1,i2,ivel
	real*8 kappa_mol(ng,nlam,nmol),wP1,wP2,wT1,wT2,x1,x2,tot,tot2,random,w1,ww
	type(databaseKtable),pointer :: Ktab

	Ktab => Ktable(imol)
	if(.not.Ktab%available.or..not.includemol(imol)) then
		kappa_mol(1:ng,1:nlam,imol)=0d0
		return
	endif

	if(P(ir).le.Ktab%P1) then
		iP=1
		wP1=1d0
		wP2=0d0
	else if(P(ir).ge.Ktab%P2) then
		iP=Ktab%nP-1
		wP1=0d0
		wP2=1d0
	else
		call hunt(Ktab%P,Ktab%nP,P(ir),iP)
c		do iP=1,Ktab%nP
c			if(P(ir).ge.Ktab%P(iP).and.P(ir).lt.Ktab%P(iP+1)) exit
c		enddo
		if(iP.lt.1) then
			iP=1
			wP1=1d0
			wP2=0d0
		else if(iP.ge.Ktab%nP) then
			iP=Ktab%nP-1
			wP1=0d0
			wP2=1d0
		else
			wP1=1d0-log10(P(ir)/Ktab%P(iP))/log10(Ktab%P(iP+1)/Ktab%P(iP))
			wP2=1d0-wP1
		endif
	endif

	if(T(ir).le.Ktab%T1) then
		iT=1
		wT1=1d0
		wT2=0d0
	else if(T(ir).ge.Ktab%T2) then
		iT=Ktab%nT-1
		wT1=0d0
		wT2=1d0
	else
		call hunt(Ktab%T,Ktab%nT,T(ir),iT)
c		do iT=1,Ktab%nT
c			if(T(ir).ge.Ktab%T(iT).and.T(ir).lt.Ktab%T(iT+1)) exit
c		enddo
		if(iT.lt.1) then
			iT=1
			wT1=1d0
			wT2=0d0
		else if(iT.ge.Ktab%nT) then
			iT=Ktab%nT-1
			wT1=0d0
			wT2=1d0
		else
			wT1=1d0-log10(T(ir)/Ktab%T(iT))/log10(Ktab%T(iT+1)/Ktab%T(iT))
			wT2=1d0-wT1
		endif
	endif
 
	do ilam=1,nlam
		if(computelam(ilam)) then
			kappa_mol(1:ng,ilam,imol)=Ktab%ktable(1:ng,ilam,iT,iP,ivel)*wT1*wP1+
     &			Ktab%ktable(1:ng,ilam,iT+1,iP,ivel)*wT2*wP1+
     &			Ktab%ktable(1:ng,ilam,iT,iP+1,ivel)*wT1*wP2+
     &			Ktab%ktable(1:ng,ilam,iT+1,iP+1,ivel)*wT2*wP2
		endif
	enddo

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
		
		
	subroutine regridKtable(k0,w0,n0,g1,k1,w1,n1, g0,b0,gg1)
	IMPLICIT NONE
	integer,intent(in) :: n0,n1
	real*8,intent(in) :: g1(n1),w1(n1)
	real*8,intent(out) :: k1(n1),g0(n0),b0(n0+1),gg1(n1)
	real*8,intent(inout) :: k0(n0),w0(n0)
	integer ig,j,j1,j2
	real*8 tot0,tot1,ww,bg0,bg1
	
	tot0=dot_product(k0,w0)
	tot0=tot0/sum(w0(1:n0))
	call dpquicksort_w(k0,w0,n0)
c	call sortw_2(k0,w0,n0)
	if(n1.eq.1) then
		k1(1)=tot0
	else
		g0=w0
		do ig=2,n0
			g0(ig)=g0(ig)+g0(ig-1)
		enddo
		g0(1:n0)=g0(1:n0)/g0(n0)
		w0=w0/sum(w0(1:n0))
		gg1=w1
		do ig=2,n1
			gg1(ig)=gg1(ig)+gg1(ig-1)
		enddo
		gg1(1:n1)=gg1(1:n1)/gg1(n1)
		b0(1)=0d0
		do ig=1,n0
			b0(ig+1)=g0(ig)
		enddo
		bg0=0d0
		do ig=1,n1
			bg1=gg1(ig)
			k1(ig)=0d0
			do j1=1,n0
				if(b0(j1+1).ge.bg0) exit
			enddo
c			call hunt(b0,n0,bg0,j1)
			do j2=j1,n0
				if(b0(j2+1).ge.bg1) exit
			enddo
c			call hunt(b0(j1:n0),n0-j1+1,bg1,j2)
c			j2=j2+j1-1
			if(j1.gt.n0) j1=n0
			if(j2.gt.n0) j2=n0
			if(j1.ge.j2) then
				k1(ig)=k0(j1)
			else if((j2-j1).eq.1) then
				k1(ig)=(k0(j1)*(b0(j1+1)-bg0)+k0(j2)*(bg1-b0(j2)))/(bg1-bg0)
			else
			k1(ig)=(k0(j1)*(b0(j1+1)-bg0)+k0(j2)*(bg1-b0(j2)))
				do j=j1+1,j2-1
					k1(ig)=k1(ig)+k0(j)*(b0(j+1)-b0(j))
				enddo
				k1(ig)=k1(ig)/(bg1-bg0)
			endif
			bg0=bg1
		enddo
	endif
	
	return
	end
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
      SUBROUTINE sortw_2(arr,brr,n)
      INTEGER n,M,NSTACK
      REAL*8 arr(n),brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,b,temp
	if(n.le.1) return
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=l-1
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
          temp=brr(l)
          brr(l)=brr(l+1)
          brr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
        b=brr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        brr(l+1)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
      SUBROUTINE sortidx_2(arr,brr,n)
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      integer brr(n),b,itemp
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,temp
	if(n.le.1) return
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=l-1
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        itemp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=itemp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          itemp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=itemp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          itemp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=itemp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
          itemp=brr(l)
          brr(l)=brr(l+1)
          brr(l+1)=itemp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
        b=brr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        itemp=brr(i)
        brr(i)=brr(j)
        brr(j)=itemp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        brr(l+1)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END



  ! dual pivot quicksort
	recursive subroutine dpquicksort_w(array,warray,last)
	IMPLICIT NONE
	integer, intent(in) :: last
	real*8, intent(inout) :: array(last),warray(last)
	real*8 :: temp,p1,p2,wtemp,wp1,wp2
	integer :: i,j,l,k,g
	
	if (last.lt.40) then ! use insertion sort on small arrays
		do i=2,last
			temp=array(i)
			wtemp=warray(i)
			do j=i-1,1,-1
				if (array(j).le.temp) exit
				array(j+1)=array(j)
				warray(j+1)=warray(j)
			enddo
			array(j+1)=temp
			warray(j+1)=wtemp
		enddo
		return
	endif
	p1=array(last/3)
	p2=array(2*last/3)
	wp1=warray(last/3)
	wp2=warray(2*last/3)
	if (p2.lt.p1) then
		temp=p1
		p1=p2
		p2=temp
		wtemp=wp1
		wp1=wp2
		wp2=wtemp
	endif
	array(last/3)=array(1)
	array(1)=p1
	array(2*last/3)=array(last)
	array(last)=p2
	
	warray(last/3)=warray(1)
	warray(1)=wp1
	warray(2*last/3)=warray(last)
	warray(last)=wp2
	
	g=last
	l=2
	do while (array(l).lt.p1)
		l=l+1
	enddo
	k=l
	
	do while(k.lt.g)
		temp=array(k)
		wtemp=warray(k)
		if (temp.lt.p1) then
			array(k)=array(l)
			array(l)=temp
			warray(k)=warray(l)
			warray(l)=wtemp
			l=l+1
		else if (temp.gt.p2) then
			do while(array(g-1).gt.p2)
				g=g-1
			enddo
			if (k.ge.g) exit
			g=g-1
			if (array(g).lt.p1) then
				array(k)=array(l)
				array(l)=array(g)
				array(g)=temp
				warray(k)=warray(l)
				warray(l)=warray(g)
				warray(g)=wtemp
				l=l+1
			else
				array(k)=array(g)
				array(g)=temp
				warray(k)=warray(g)
				warray(g)=wtemp
			endif
		endif
		k=k+1
	enddo
	if (l.gt.2) then
		array(1)=array(l-1)
		array(l-1)=p1
		warray(1)=warray(l-1)
		warray(l-1)=wp1
		call dpquicksort_w(array(1:l-2),warray(1:l-2),l-2)
	endif
	call dpquicksort_w(array(l:g-1),warray(l:g-1),g-l)
	if (g.lt.last) then
		array(last)=array(g)
		array(g)=p2
		warray(last)=warray(g)
		warray(g)=wp2
		call dpquicksort_w(array(g+1:last),warray(g+1:last),last-g)
	endif
	
	return
	end



  ! dual pivot quicksort
	recursive subroutine dpquicksort_indx(array,warray,last)
	IMPLICIT NONE
	integer, intent(in) :: last
	real*8, intent(inout) :: array(last)
	integer, intent(inout) :: warray(last)
	real*8 :: temp,p1,p2
	integer :: i,j,l,k,g,wtemp,wp1,wp2
	
	if (last.lt.40) then ! use insertion sort on small arrays
		do i=2,last
			temp=array(i)
			wtemp=warray(i)
			do j=i-1,1,-1
				if (array(j).le.temp) exit
				array(j+1)=array(j)
				warray(j+1)=warray(j)
			enddo
			array(j+1)=temp
			warray(j+1)=wtemp
		enddo
		return
	endif
	p1=array(last/3)
	p2=array(2*last/3)
	wp1=warray(last/3)
	wp2=warray(2*last/3)
	if (p2.lt.p1) then
		temp=p1
		p1=p2
		p2=temp
		wtemp=wp1
		wp1=wp2
		wp2=wtemp
	endif
	array(last/3)=array(1)
	array(1)=p1
	array(2*last/3)=array(last)
	array(last)=p2
	
	warray(last/3)=warray(1)
	warray(1)=wp1
	warray(2*last/3)=warray(last)
	warray(last)=wp2
	
	g=last
	l=2
	do while (array(l).lt.p1)
		l=l+1
	enddo
	k=l
	
	do while(k.lt.g)
		temp=array(k)
		wtemp=warray(k)
		if (temp.lt.p1) then
			array(k)=array(l)
			array(l)=temp
			warray(k)=warray(l)
			warray(l)=wtemp
			l=l+1
		else if (temp.gt.p2) then
			do while(array(g-1).gt.p2)
				g=g-1
			enddo
			if (k.ge.g) exit
			g=g-1
			if (array(g).lt.p1) then
				array(k)=array(l)
				array(l)=array(g)
				array(g)=temp
				warray(k)=warray(l)
				warray(l)=warray(g)
				warray(g)=wtemp
				l=l+1
			else
				array(k)=array(g)
				array(g)=temp
				warray(k)=warray(g)
				warray(g)=wtemp
			endif
		endif
		k=k+1
	enddo
	if (l.gt.2) then
		array(1)=array(l-1)
		array(l-1)=p1
		warray(1)=warray(l-1)
		warray(l-1)=wp1
		call dpquicksort_indx(array(1:l-2),warray(1:l-2),l-2)
	endif
	call dpquicksort_indx(array(l:g-1),warray(l:g-1),g-l)
	if (g.lt.last) then
		array(last)=array(g)
		array(g)=p2
		warray(last)=warray(g)
		warray(g)=wp2
		call dpquicksort_indx(array(g+1:last),warray(g+1:last),last-g)
	endif
	
	return
	end

