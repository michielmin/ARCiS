	subroutine SetupOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer imol
	real*8 kappa(ng),nu1,nu2,opac_tot(nlam,ng)
	real*8 x1,x2,rr,gasdev,random
	real*8,allocatable :: k_line(:),nu_line(:)
	integer n_nu_line
	integer i,j,ir

	n_voigt=1d6
	allocate(a_therm(n_voigt))
	allocate(a_press(n_voigt))
	do i=1,n_voigt
		x1=gasdev(idum)/sqrt(2d0)
		x2=tan((random(idum)-0.5d0)*pi)
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
		call LineStrengthWidth(ir)
		do i=1,nlam-1
			call tellertje(i,nlam-1)
			nu1=freq(i+1)
			nu2=freq(i)
			call ComputeKtable(ir,nu1,nu2,kappa,epsCk)
			opac(ir,i,1:ng)=kappa(1:ng)
			opac_tot(i,1:ng)=opac_tot(i,1:ng)+opac(ir,i,1:ng)*Ndens(ir)*(R(ir+1)-R(ir))
		enddo
	enddo
	
	open(unit=30,file=trim(outputdir) // "opticaldepth.dat",RECL=6000)
	do i=1,nlam-1
		write(30,*) sqrt(lam(i)*lam(i+1))/micron,sum(opac_tot(i,1:ng))/real(ng),opac_tot(i,1:ng)
	enddo
	close(unit=30)

	return
	end
	
	subroutine LineStrengthWidth(ir)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 w,x1,x2,x3,x4
	integer imol,iT,i,ir,iiso

	call hunt(TZ,nTZ,T(ir),iT)

	do i=1,nlines
		imol=Lines(i)%imol
		iiso=Lines(i)%iiso
c thermal broadening
		w=(sqrt(2d0*kb*T(ir)/(Mmol(imol)*mp)))
		Lines(i)%a_therm=w*Lines(i)%freq/clight
c pressure broadening
		Lines(i)%a_press=Lines(i)%gamma_air*P(ir)*(1d0-mixrat(imol))
		Lines(i)%a_press=Lines(i)%a_press+Lines(i)%gamma_air*P(ir)*mixrat(imol)
		Lines(i)%a_press=Lines(i)%a_press*(296d0/T(ir))**Lines(i)%n
c line strength
		x1=exp(-hplanck*clight*Lines(i)%Elow/(kb*T(ir)))
		x2=exp(-hplanck*clight*Lines(i)%freq/(kb*T(ir)))
		x3=exp(-hplanck*clight*Lines(i)%Elow/(kb*296d0))
		x4=exp(-hplanck*clight*Lines(i)%freq/(kb*296d0))
		Lines(i)%S=Lines(i)%S0*(x1*(1d0-x2))/(x3*ZZ(imol,iiso,iT)*(1d0-x4))
	enddo
	
	return
	end
	
	
	subroutine ComputeKtable(ir,nu1,nu2,kappa,eps)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 nu1,nu2,kappa(ng),g(ng),w,dnu,gamma,fact,kappa0(ng),eps
	real*8,allocatable :: nu(:),kline(:),kdis(:),dis(:),kline0(:)
	real*8 Eu,El,A,x,kmax,kmin,V,scale,x1,x2,gasdev,random,rr,gu,gl
	integer nnu,inu,iT,imol,i,ju,jl,j,nkdis,NV,nl,k,iiso,ir,NV0,iter,maxiter
	integer i_therm,i_press
	logical converged
	real*8 f,a_t,a_p
	
	do i=1,ng
		g(i)=(real(i)-0.5)/real(ng)
	enddo

	maxiter=5

	fact=50d0

	dnu=0.01d0
	nl=0
	do i=1,nlines
		gamma=fact*sqrt(Lines(i)%a_therm**2+Lines(i)%a_press**2)
		if((Lines(i)%freq+gamma).gt.nu1.and.(Lines(i)%freq-gamma).lt.nu2) nl=nl+1
	enddo

	nnu=abs(nu2-nu1)/dnu

	allocate(nu(nnu))
	allocate(kline(nnu))
	allocate(kline0(nnu))
	
	do inu=1,nnu
		nu(inu)=nu1+real(inu-1)*(nu2-nu1)/real(nnu-1)
	enddo
	scale=real(nnu)/(nu2-nu1)

	NV=nnu*100/(nl+1)+50

	kline=0d0
	call hunt(TZ,nTZ,T(ir),iT)

	nkdis=1000
	allocate(dis(nkdis))
	allocate(kdis(nkdis))

	converged=.false.
	iter=0
	do while(.not.converged)

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,imol,iiso,gamma,A,a_t,a_p,f,x1,x2,rr,x,inu)
!$OMP& SHARED(fact,Lines,mixrat,scale,NV,kline,nnu,nu1,nu2,nlines)
!$OMP DO SCHEDULE(STATIC,1)
	do i=1,nlines
		gamma=fact*sqrt(Lines(i)%a_therm**2+Lines(i)%a_press**2)
		if((Lines(i)%freq+gamma).gt.nu1.and.(Lines(i)%freq-gamma).lt.nu2) then
			imol=Lines(i)%imol
			iiso=Lines(i)%iiso
			A=Lines(i)%S*mixrat(imol)*scale
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
c				x1=gasdev(idum)/sqrt(2d0)
				x1=a_therm(i_therm)
				x1=x1*a_t
c				x2=tan((random(idum)-0.5d0)*pi)
				x2=a_press(i_press)
				x2=x2*a_p
				x=x1+x2
c				rr=random(idum)
c				if(rr.lt.0.25) then
c					x=x1+x2
c				else if(rr.lt.0.5) then
c					x=-x1+x2
c				else if(rr.lt.0.75) then
c					x=x1-x2
c				else
c					x=-x1-x2
c				endif
				x=(x+f-nu1)*scale
				inu=int(x)
				if(inu.ge.1.and.inu.le.nnu) then
					kline(inu)=kline(inu)+A
				endif
			enddo
		endif
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	if(iter.gt.0) then
		kline=(real(NV)*kline+real(NV0)*kline0)/real(NV+NV0)
	endif

	kmin=1d200
	kmax=0d0
	do inu=1,nnu	
		if(kline(inu).gt.kmax) kmax=kline(inu)
		if(kline(inu).lt.kmin) kmin=kline(inu)
	enddo
	if(kmin.lt.kmax/1d30) kmin=kmax/1d30
	kmin=log10(kmin)
	kmax=log10(kmax)

	do i=1,nkdis
		kdis(i)=10d0**(kmin+real(i-1)*(kmax-kmin)/real(nkdis-1))
	enddo
	dis=0d0
	do inu=1,nnu
		do i=nkdis,1,-1
			if(kline(inu).le.kdis(i)) then
				dis(i)=dis(i)+1d0
			else
				exit
			endif
		enddo
	enddo

	dis(1:nkdis-1)=dis(1:nkdis-1)/dis(nkdis)
	dis(nkdis)=1d0

	do i=1,ng
		if(g(i).lt.dis(1)) then
			kappa(i)=kdis(1)
		else
			call hunt(dis,nkdis,g(i),j)
			if(j.lt.nkdis) then
				if(abs(dis(j+1)-g(i)).lt.abs(dis(j)-g(i))) j=j+1
			endif
			kappa(i)=kdis(j)
		endif
	enddo
	kappa(ng)=10d0**kmax

	if(iter.gt.maxiter) then
		converged=.true.
	else if(iter.eq.0) then
		kappa0=kappa
		kline0=kline
		NV0=NV
		converged=.false.
		iter=iter+1
	else if(maxval(abs(kappa0-kappa)/(kappa0+kappa)).gt.eps) then
		kappa0=kappa
		kline0=kline
		if(iter.gt.0) then
			NV=NV0*2
		endif
		NV0=NV
		converged=.false.
		iter=iter+1
	else
		converged=.true.
	endif

	enddo

	deallocate(dis)
	deallocate(kdis)

	deallocate(nu)
	deallocate(kline)

	return
	end


