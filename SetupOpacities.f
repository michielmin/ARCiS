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
	real*8,allocatable :: nu_line(:),dnu_line(:),mixrat_tmp(:)
	real*8,allocatable :: opac_tot(:,:),cont_tot(:),kaver(:),kappa_mol(:,:,:)
	integer,allocatable,save :: ifull(:),ifast(:)
	real*8,allocatable,save :: k_line(:),ktemp(:),kappa(:),w_line(:),kappa_tot(:),work1(:),work2(:),work3(:)
!$OMP THREADPRIVATE(ifull,ifast,k_line,ktemp,kappa,w_line,kappa_tot,work1,work2,work3)
	integer n_nu_line,iT,nfull,nfast
	integer i,j,ir,k,nl,ig,ig_c,imol0,jg
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
	if(.not.allocated(k_line)) allocate(k_line(n_nu_line))
	if(.not.allocated(ktemp)) allocate(ktemp(ng))
	if(.not.allocated(kappa)) allocate(kappa(ng))
	if(.not.allocated(w_line)) allocate(w_line(n_nu_line))
	if(.not.allocated(ifull)) allocate(ifull(nmol))
	if(.not.allocated(ifast)) allocate(ifast(nmol))
	if(.not.allocated(kappa_tot)) allocate(kappa_tot(0:nmol))
	if(.not.allocated(work1)) allocate(work1(n_nu_line))
	if(.not.allocated(work2)) allocate(work2(n_nu_line))
	if(.not.allocated(work3)) allocate(work3(n_nu_line))
!$OMP END PARALLEL

	do ir=nr,1,-1
		call tellertje(nr-ir+1,nr)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(cont_tot,nlam)
!$OMP DO SCHEDULE(STATIC)
		do i=1,nlam
			cont_tot(i)=0d0
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

c===============
c UV cross sections of CO2 from Venot et al.
		call CO2_UV_cross(lam,cont_tot,nlam,min(T(ir),800d0))
		cont_tot=cont_tot*mixrat_r(ir,2)
		do i=1,nlam
			if(lam(i).gt.0.3d-4) cont_tot(i)=0d0
		enddo
c===============

		mixrat_tmp(1:nmol)=mixrat_r(ir,1:nmol)
		do i=1,ncia
			if(includemol(CIA(i)%imol1).and.includemol(CIA(i)%imol2)) then
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
			endif
		enddo
		Csca(ir,1:nlam)=0d0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,j,imol,ig,jg,ig_c,imol0,w1,nfull,nfast)
!$OMP& SHARED(nlam,n_nu_line,nmol,mixrat_tmp,ng,ir,kappa_mol,cont_tot,Cabs,Csca,opac_tot,Ndens,R,computelam,
!$OMP&        ig_comp,retrieval,domakeai,gg,wgg,ng_comp,opacitymol,emisspec,computeT,lamemis,useobsgrid,
!$OMP&        RTgridpoint,includemol,do_rayleigh,mixrat_PAH,mixrat_optEC)

		if(do_rayleigh) then
!$OMP DO SCHEDULE(STATIC)
			do i=1,nlam
				if(computelam(i)) call RayleighScattering(Csca(ir,i),ir,i)
			enddo
!$OMP END DO
!$OMP FLUSH
		endif

		if(mixrat_PAH.gt.0d0) call ComputePAH(cont_tot,Csca(ir,1:nlam),computelam)
		if(mixrat_optEC.gt.0d0) call Compute_optEC(cont_tot,Csca(ir,1:nlam),computelam)

!$OMP DO SCHEDULE(STATIC)
		do i=1,nlam
			kappa_mol(1:ng,i,1:nmol)=0d0
		enddo
!$OMP END DO
!$OMP FLUSH

!$OMP DO SCHEDULE(DYNAMIC)
		do imol=1,nmol
			if(includemol(imol)) then
				call ReadOpacityFITS(kappa_mol,imol,ir)
			endif
		enddo
!$OMP END DO
!$OMP FLUSH

!$OMP DO SCHEDULE(DYNAMIC)
		do i=1,nlam
			if(computelam(i).and.(emisspec.or.computeT).and.(.not.useobsgrid.or.lamemis(i).or.RTgridpoint(i))) then
			kappa_tot(0)=cont_tot(i)
			do imol=1,nmol
				if(opacitymol(imol)) then
					kappa_tot(imol)=0d0
					do ig=1,ng
						kappa_tot(imol)=kappa_tot(imol)+wgg(ig)*kappa_mol(ig,i,imol)
					enddo
					kappa_tot(imol)=kappa_tot(imol)*mixrat_tmp(imol)
					kappa_tot(0)=kappa_tot(0)+kappa_tot(imol)
				endif
			enddo
			nfull=0
			nfast=0d0
			do imol=1,nmol
				if(opacitymol(imol)) then
					if(kappa_tot(imol).ge.0.01*kappa_tot(0)) then
						nfull=nfull+1
						ifull(nfull)=imol
						kappa_tot(nfull)=-kappa_tot(imol)
					else
						nfast=nfast+1
						ifast(nfast)=imol
					endif
				endif
			enddo
			call dpquicksort_indx(kappa_tot(1:nfull),ifull(1:nfull),nfull)
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
					call regridKtable(k_line,w_line,n_nu_line,gg,kappa,wgg,ng,work1,work2,work3)
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

			kappa=kappa+cont_tot(i)
			Cabs(ir,i,1:ng)=kappa(1:ng)
			do ig=1,ng
				if(Cabs(ir,i,ig).lt.1d-6*Csca(ir,i)) Cabs(ir,i,ig)=1d-6*Csca(ir,i)
			enddo
			if(.not.retrieval) opac_tot(i,1:ng)=opac_tot(i,1:ng)+(Cabs(ir,i,1:ng)+Csca(ir,i))*Ndens(ir)*(R(ir+1)-R(ir))
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
	real*8 ll,Cs,ll2
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
c					Cs=Cs+12.4*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
					Cs=Cs+9.950903042454394e-44*mixrat_r(ir,j)*ll2
				case(4) !N2O
c Sneep & Ubachs (2005)
c					Cs=Cs+15.9*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
					Cs=Cs+1.2759625675405231e-43*mixrat_r(ir,j)*ll2
				case(5) !CO
c Sneep & Ubachs (2005)
c					Cs=Cs+6.19*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
					Cs=Cs+4.967426599418766e-44*mixrat_r(ir,j)*ll2
				case(6) !CH4
c Sneep & Ubachs (2005)
c					Cs=Cs+12.47*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
					Cs=Cs+1.000707749511341e-43*mixrat_r(ir,j)*ll2
				case(22) !N2
c Sneep & Ubachs (2005)
c					Cs=Cs+5.1*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
					Cs=Cs+4.0927101222997906e-44*mixrat_r(ir,j)*ll2
				case(30) !SF6
c Sneep & Ubachs (2005)
c					Cs=Cs+32.3*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
					Cs=Cs+2.5920497441232007e-43*mixrat_r(ir,j)*ll2
				case(45) ! H2
c Dalgarno & Williams (1962)
					Cs=Cs+mixrat_r(ir,j)*(8.14d-45*ll2+1.28d-54*ll2*ll+1.61d-64*ll2*ll2)
				case(1) !H2O
c					Cs=Cs+4.43*mixrat_r(ir,j)*1e-27*(sqrt(ll)/18788.4)**4
					Cs=Cs+3.555040361134916e-44*mixrat_r(ir,j)*ll2
				case(48) ! He
c https://www.climate-policy-watcher.org/surface-temperature/scattering-by-molecules-rayleigh-scattering.html
					Cs=Cs+0.0641*mixrat_r(ir,j)*(8.14d-45*ll2+1.28d-54*ll2*ll+1.61d-64*ll2*ll2)
				case(11) ! NH3
c https://www.climate-policy-watcher.org/surface-temperature/scattering-by-molecules-rayleigh-scattering.html
					Cs=Cs+7.3427*mixrat_r(ir,j)*(8.14d-45*ll2+1.28d-54*ll2*ll+1.61d-64*ll2*ll2)
			end select
		endif
	enddo

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
 
	do ilam=1,nlam
		if(computelam(ilam)) then
			kappa_mol(1:ng,ilam,imol)=tab1(1:ng,ilam)*wT1*wP1+tab2(1:ng,ilam)*wT2*wP1+tab3(1:ng,ilam)*wT1*wP2+tab4(1:ng,ilam)*wT2*wP2
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
	
	tot0=0d0
	do ig=1,n0
		tot0=tot0+k0(ig)*w0(ig)
	enddo
	tot0=tot0/sum(w0(1:n0))
c	call sortw_2(k0,w0,n0)
	call dpquicksort_w(k0,w0,n0)
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

