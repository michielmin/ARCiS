	subroutine SetupOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer imol
	real*8 kappa(ng),nu1,nu2,opac_tot(nlam,ng),tanscale
	real*8 x1,x2,rr,gasdev,random,dnu,cia_tot(nlam)
	real*8,allocatable :: k_line(:),nu_line(:),dnu_line(:)
	integer n_nu_line,iT
	integer i,j,ir,inu1(nlam),inu2(nlam)
	character*500 filename

	n_voigt=1d6
	allocate(a_therm(n_voigt))
	allocate(a_press(n_voigt))
	tanscale=atan(cutoff_lor)
	do i=1,n_voigt
		x1=gasdev(idum)/sqrt(2d0)
		x2=tan((random(idum)-0.5d0)*2d0*tanscale)
		rr=random(idum)
		if(rr.lt.0.25) then
			x1=-x1
		else if(rr.lt.0.5) then
			x2=-x2
		else if(rr.lt.0.75) then
			x1=-x1
			x2=-x2
		endif
		a_therm(i)=x1
		a_press(i)=x2
	enddo

	opac_tot=0d0

	do ir=1,nr
		call output("Opacities for layer: " // 
     &		trim(int2string(ir,'(i4)')) // " of " // trim(int2string(nr,'(i4)')))
		call output("T = " // trim(dbl2string(T(ir),'(f8.2)')) // " K")
		call output("P = " // trim(dbl2string(P(ir),'(es8.2)')) // " Ba")
		call LineStrengthWidth(ir,dnu,freq(nlam),freq(1))
		dnu=dnu/3d0
		n_nu_line=abs(freq(1)/freq(nlam))/dnu
		nu1=freq(1)
		n_nu_line=1
		do while(nu1.gt.freq(nlam))
			nu1=nu1-nu1*dnu
			n_nu_line=n_nu_line+1
		enddo
		allocate(k_line(n_nu_line))
		allocate(nu_line(n_nu_line))
		allocate(dnu_line(n_nu_line))
		do i=1,n_nu_line
			nu_line(i)=exp(log(freq(1))+log(freq(nlam)/freq(1))*real(i-1)/real(n_nu_line-1))
			dnu_line(i)=nu_line(i)*dnu
		enddo
		cia_tot=0d0
		do i=1,ncia
			if(T(ir).lt.CIA(i)%T(1)) then
				iT=1
			else if(T(iR).gt.CIA(i)%T(CIA(i)%nT)) then
				iT=CIA(i)%nT
			else
				do iT=1,CIA(i)%nT-1
					if(T(ir).ge.CIA(i)%T(iT).and.T(ir).le.CIA(i)%T(iT+1)) exit
				enddo
			endif
			cia_tot(1:nlam)=cia_tot(1:nlam)+CIA(i)%Cabs(iT,1:nlam)*Ndens(ir)*cia_mixrat(CIA(i)%imol1)*cia_mixrat(CIA(i)%imol2)
		enddo
		call ComputeKline(ir,nu_line,k_line,n_nu_line,dnu_line)
c		call ComputeKline_exact(ir,nu_line,k_line,n_nu_line,dnu_line)
		if(outputopacity) call WriteOpacity(ir,"line",nu_line,k_line,n_nu_line,1)
		call tellertje(1,nlam-1)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,nu1,nu2,kappa)
!$OMP& SHARED(nlam,freq,ir,nu_line,k_line,n_nu_line,cia_tot,opac,opac_tot,Ndens,R,ng)
!$OMP DO SCHEDULE(STATIC, 1)
		do i=1,nlam-1
			call tellertje(i+1,nlam+1)
			nu1=freq(i+1)
			nu2=freq(i)
			call ComputeKtable(ir,nu1,nu2,nu_line,k_line,n_nu_line,kappa,cia_tot(i))
			opac(ir,i,1:ng)=kappa(1:ng)
			opac_tot(i,1:ng)=opac_tot(i,1:ng)+opac(ir,i,1:ng)*Ndens(ir)*(R(ir+1)-R(ir))
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		if(outputopacity) call WriteOpacity(ir,"aver",freq,opac(ir,1:nlam,1:ng),nlam,ng)
		call tellertje(nlam,nlam)
	
	open(unit=30,file=trim(outputdir) // "opticaldepth.dat",RECL=6000)
	write(30,'("#",a13,a19)') "lambda [mu]","total average tau"
	do i=1,nlam-1
		write(30,'(f12.6,e19.7)') sqrt(lam(i)*lam(i+1))/micron,sum(opac_tot(i,1:ng))/real(ng)
	enddo
	close(unit=30)

		deallocate(k_line)
		deallocate(nu_line)
		deallocate(dnu_line)
	enddo

	return
	end
	
	subroutine LineStrengthWidth(ir,minw,nu1,nu2)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 w,x1,x2,x3,x4,minw,nu1,nu2
	integer imol,iT,i,ir,iiso

	minw=0.1d0
	call hunt(TZ,nTZ,T(ir),iT)

	do i=1,nlines
		imol=Lines(i)%imol
		iiso=Lines(i)%iiso
c thermal broadening
		w=(sqrt(2d0*kb*T(ir)/(Mmol(imol)*mp)))
		Lines(i)%a_therm=w*Lines(i)%freq/clight
c pressure broadening
		Lines(i)%a_press=Lines(i)%gamma_air*P(ir)*(1d0-mixrat(imol))/atm
		Lines(i)%a_press=Lines(i)%a_press+Lines(i)%gamma_self*P(ir)*mixrat(imol)/atm
		Lines(i)%a_press=Lines(i)%a_press*(296d0/T(ir))**Lines(i)%n
c line strength
		x1=exp(-hplanck*clight*Lines(i)%Elow/(kb*T(ir)))
		x2=exp(-hplanck*clight*Lines(i)%freq/(kb*T(ir)))
		x3=exp(-hplanck*clight*Lines(i)%Elow/(kb*296d0))
		x4=exp(-hplanck*clight*Lines(i)%freq/(kb*296d0))

		Lines(i)%S=Lines(i)%S0*(x1*(1d0-x2))/(x3*ZZ(imol,iiso,iT)*(1d0-x4))

		if((Lines(i)%freq).gt.nu1.and.(Lines(i)%freq).lt.nu2) then
			if(sqrt((Lines(i)%a_press*4d0)**2+Lines(i)%a_therm**2)/Lines(i)%freq.lt.minw) then
				minw=sqrt((Lines(i)%a_press*4d0)**2+Lines(i)%a_therm**2)/Lines(i)%freq
			endif
		endif
	enddo
	
	return
	end
	
	
	subroutine ComputeKtable(ir,nu1,nu2,nu,kline,nnu,kappa,Ccont)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nnu,nkdis
	parameter(nkdis=100)
	real*8 nu1,nu2,kappa(ng),g(ng),w,dnu,gamma,fact,eps
	real*8 nu(nnu),kline(nnu),Ccont,kap
	real*8 kdis(nkdis),dis(nkdis),gb(ng+1)
	real*8 Eu,El,A,x,kmax,kmin,V,scale,x1,x2,gasdev,random,rr,gu,gl
	integer inu,iT,imol,i,ju,jl,j,NV,nl,k,iiso,ir,NV0,iter,maxiter
	integer i_therm,i_press,inu1,inu2,ig(ng+1)
	logical converged
	real*8 f,a_t,a_p
	
	do i=1,ng
		g(i)=(real(i)-0.5)/real(ng)
	enddo
	do i=1,ng+1
		gb(i)=real(i-1)/real(ng)
	enddo

	call hunt(nu,nnu,nu1,inu2)
	call hunt(nu,nnu,nu2,inu1)
	inu1=inu1+1
	if(inu1.gt.nnu) inu1=nnu
	if(inu2.gt.nnu) inu2=nnu

	kmin=1d200
	kmax=Ccont
	do inu=inu1+1,inu2-1
		kap=kline(inu)+Ccont
		if(kap.gt.kmax) kmax=kap
		if(kap.lt.kmin) kmin=kap
	enddo
	if(kmax.eq.0d0) then
		kappa=0d0
		return
	endif
	if(kmin.lt.kmax/1d10) kmin=kmax/1d10
	kmin=log10(kmin)
	kmax=log10(kmax)

	do i=1,nkdis
		kdis(i)=10d0**(kmin+real(i-1)*(kmax-kmin)/real(nkdis-1))
	enddo
	dis=0d0
	do inu=inu1,inu2
		do i=nkdis,1,-1
			if((kline(inu)+Ccont).le.kdis(i)) then
				dis(i)=dis(i)+1d0
			else
				exit
			endif
		enddo
	enddo

	dis(1:nkdis-1)=dis(1:nkdis-1)/dis(nkdis)
	dis(nkdis)=1d0

	do i=1,ng+1
		call hunt(dis,nkdis,gb(i),ig(i))
		if(ig(i).lt.1) ig(i)=1
		if(ig(i).gt.nkdis) ig(i)=nkdis
	enddo
	do i=1,ng
		kappa(i)=sqrt(kdis(ig(i))*kdis(ig(i+1)))
	enddo

	return
	end



	subroutine ComputeKline(ir,nu,kline,nnu,dnu)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nnu
	real*8 w,gamma,fact
	real*8 nu(nnu),dnu(nnu)
	real*8,target :: kline(nnu)
	real*8 Eu,El,A,x,kmax,kmin,V,scale,x1,x2,gasdev,random,rr,gu,gl
	integer iT,imol,i,ju,jl,j,nkdis,NV,nl,k,iiso,ir,NV0,iter,maxiter
	integer i_therm,i_press,il,idnu,inu1,inu2,inu
	real*8 f,a_t,a_p
	
	fact=50d0
	kline=0d0

	nl=0
	do i=1,nlines
		gamma=fact*sqrt(Lines(i)%a_therm**2+Lines(i)%a_press**2)
		if((Lines(i)%freq+gamma).gt.nu(nnu).and.(Lines(i)%freq-gamma).lt.nu(1)) nl=nl+1
	enddo

	scale=real(nnu-1)/log(nu(nnu)/nu(1))

	NV=real(nnu)*100d0/real(nl+1)+25d0

	il=0
	call hunt(TZ,nTZ,T(ir),iT)

	call tellertje(1,nl)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,imol,iiso,gamma,A,a_t,a_p,f,x1,x2,rr,x,inu,i_press,i_therm,idnu,inu1,inu2)
!$OMP& SHARED(fact,Lines,mixrat,scale,NV,kline,nnu,nu,nlines,a_therm,a_press,il,nl,n_voigt,P,ir,cutoff_abs)
!$OMP DO SCHEDULE(STATIC, 1)
	do i=1,nlines
		gamma=sqrt(Lines(i)%a_therm**2+Lines(i)%a_press**2)
		if((Lines(i)%freq+fact*gamma).gt.nu(nnu).and.(Lines(i)%freq-fact*gamma).lt.nu(1)) then
			il=il+1
			call tellertje(il+1,nl+2)
			imol=Lines(i)%imol
			iiso=Lines(i)%iiso
			A=Lines(i)%S*mixrat(imol)
			a_t=Lines(i)%a_therm
			a_p=Lines(i)%a_press
			f=Lines(i)%freq
c	Random sampling of the Voigt profile
			A=A/real(NV)
			i_therm=random(idum)*real(n_voigt)
			i_press=random(idum)*real(n_voigt)
			do j=1,NV
				i_therm=i_therm+1
				i_press=i_press+1
				if(i_therm.gt.n_voigt) i_therm=1
				if(i_press.gt.n_voigt) i_press=1
				x1=a_therm(i_therm)
				x1=x1*a_t
				x2=a_press(i_press)
				x2=x2*a_p
				x=x1+x2
				if(abs(x).lt.cutoff_abs) then
					idnu=abs(x/gamma)
					x=x+f
					x=log(x/nu(1))*scale
					inu=int(x)+1
					inu1=inu-idnu
					inu2=inu+idnu
					if(inu1.le.nnu.and.inu2.ge.1) then
						if(inu1.lt.1) inu1=1
						if(inu2.gt.nnu) inu2=nnu
						kline(inu1:inu2)=kline(inu1:inu2)+A/real(inu2-inu1+1)
					endif
				endif
			enddo
		endif
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(nl,nl)

	kline=kline/dnu

	return
	end






	subroutine ComputeKline_exact(ir,nu,kline,nnu,dnu)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nnu
	real*8 w,gamma,fact
	real*8 nu(nnu),dnu(nnu)
	real*8,target :: kline(nnu)
	real*8 Eu,El,A,x,kmax,kmin,V,scale,x1,x2,gasdev,random,rr,gu,gl
	integer iT,imol,i,ju,jl,j,nkdis,NV,nl,k,iiso,ir,NV0,iter,maxiter
	integer i_therm,i_press,il,inu
	real*8 f,a_t,a_p,y
	real*8,allocatable :: psi(:),x_psi(:)
	integer n_psi
	
	n_psi=1d6
	allocate(psi(n_psi))
	allocate(x_psi(n_psi))
	do i=1,n_psi
		x_psi(i)=real(i-1)/100d0
		psi(i)=1d0/(x_psi(i)**2+1d0)/sqrt(3.1415926536)
	enddo
	
	fact=50d0
	kline=0d0

	nl=0
	do i=1,nlines
		gamma=fact*sqrt(Lines(i)%a_therm**2+Lines(i)%a_press**2)
		if((Lines(i)%freq+gamma).gt.nu(nnu).and.(Lines(i)%freq-gamma).lt.nu(1)) nl=nl+1
	enddo

	scale=real(nnu-1)/log(nu(nnu)/nu(1))

	NV=real(nnu)*100d0/real(nl+1)+10d0

	il=0
	call hunt(TZ,nTZ,T(ir),iT)

	call tellertje(1,nl)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,imol,iiso,gamma,A,a_t,a_p,f,x1,x2,rr,x,inu,i_press,i_therm,y,j)
!$OMP& SHARED(fact,Lines,mixrat,scale,NV,kline,nnu,nu,nlines,a_therm,a_press,il,nl,n_voigt,P,ir,n_psi,psi)
!$OMP DO SCHEDULE(STATIC, 1)
	do i=1,nlines
		gamma=fact*sqrt(Lines(i)%a_therm**2+Lines(i)%a_press**2)
		if((Lines(i)%freq+gamma).gt.nu(nnu).and.(Lines(i)%freq-gamma).lt.nu(1)) then
			il=il+1
			call tellertje(il+1,nl+2)
			imol=Lines(i)%imol
			iiso=Lines(i)%iiso
			A=Lines(i)%S*mixrat(imol)
			a_t=Lines(i)%a_therm
			a_p=Lines(i)%a_press
			f=Lines(i)%freq
			y=a_p/a_t
			do inu=1,nnu
				x=abs(nu(inu)-f)
				x=x/a_t
				if(x.gt.8d0) then
					j=x*100d0+1
					if(j.lt.n_psi.and.j.gt.0) then
						kline(inu)=kline(inu)+psi(j)*A/y
					endif
				endif
			enddo
			A=A/real(NV)
			i_therm=random(idum)*real(n_voigt)
			i_press=random(idum)*real(n_voigt)
			do j=1,NV
				i_therm=i_therm+1
				i_press=i_press+1
				if(i_therm.gt.n_voigt) i_therm=1
				if(i_press.gt.n_voigt) i_press=1
				x1=a_therm(i_therm)
				x1=x1*a_t
				x2=a_press(i_press)
				x2=x2*a_p
				if(abs(x1+x2).lt.8d0*a_t) then
					x=x1+x2+f
					x=log(x/nu(1))*scale
					inu=int(x)+1
					if(inu.ge.1.and.inu.le.nnu) kline(inu)=kline(inu)+A
				endif
			enddo
		endif
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(nl,nl)

	kline=kline/dnu

	return
	end


