	subroutine SetupOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer imol
	real*8 kappa(100),g(100),nu1,nu2,Temp,dens0
	integer ng,i,j
	
	call ReadHITRAN()

	Temp=200d0
	dens0=1d-10
	call LineStrengthWidth(dens0,Temp)

	ng=100
	do i=1,nlam-1
		call tellertje(i,nlam-1)
		nu1=freq(i+1)
		nu2=freq(i)
		call ComputeKtable(dens0,Temp,nu1,nu2,kappa,g,ng)
		do j=1,ng
			write(90,*) sqrt(lam(i)*lam(i+1))/micron,sum(kappa)/real(ng),kappa(j)
		enddo
	enddo
	
	return
	end
	
	subroutine LineStrengthWidth(dens0,Temp)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 dens0,Temp,w
	integer imol,iT,i
	
	do i=1,nlines
		call tellertje(i,nlines)
		imol=Lines(i)%imol
		w=sqrt(2d0*kb*Temp/(Mmol(imol)*mp))
		Lines(i)%a_therm=w*Lines(i)%freq/clight
		Lines(i)%a_press=Lines(i)%a_therm*500d0
	enddo
	
	return
	end
	
	
	subroutine ComputeKtable(dens0,Temp,nu1,nu2,kappa,g,ng)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ng
	real*8 dens0,Temp,nu1,nu2,kappa(ng),g(ng),w,dnu,gamma
	real*8,allocatable :: nu(:),kline(:),kdis(:),dis(:)
	real*8 Eu,El,A,x,kmax,kmin,V,scale,x1,x2,gasdev,random,rr,gu,gl
	integer nnu,inu,iT,imol,i,ju,jl,j,nkdis,NV,nl,k,iiso
	
	do i=1,ng
		g(i)=(real(i)-0.5)/real(ng)
	enddo

	nnu=10d0*(nu2/clight-nu1/clight)
	allocate(nu(nnu))
	allocate(kline(nnu))
	do inu=1,nnu
		nu(inu)=nu1+real(inu-1)*(nu2-nu1)/real(nnu-1)
	enddo
	scale=real(nnu)/(nu2-nu1)

	nl=0
	do i=1,nlines
		gamma=50d0*sqrt(Lines(i)%a_therm**2+Lines(i)%a_press**2)
		if((Lines(i)%freq+gamma).gt.nu1.and.(Lines(i)%freq-gamma).lt.nu2) nl=nl+1
	enddo
	NV=nnu*100/(nl+1)+100

	kline=0d0
	call hunt(TZ,nTZ,Temp,iT)

c This part needs some serious attention!!!

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,gamma,ju,jl,Eu,El,A,x1,x2,x,inu,j,V)
!$OMP& SHARED(Mol,Temp,iT,NV,imol,nu1,nu2,scale,nnu,kline,w)
!$OMP DO SCHEDULE(STATIC,j)
	do i=1,nlines
		gamma=50d0*sqrt(Lines(i)%a_therm**2+Lines(i)%a_press**2)
		if((Lines(i)%freq+gamma).gt.nu1.and.(Lines(i)%freq-gamma).lt.nu2) then
			imol=Lines(i)%imol
			iiso=Lines(i)%iiso
			gu=Lines(i)%gu
			gl=Lines(i)%gl
			Eu=Lines(i)%Eup
			El=Eu-Lines(i)%freq*hplanck
			A=gu*Lines(i)%Aul*(exp(-El/Temp)-exp(-Eu/Temp))
			A=A/(Lines(i)%freq**3*ZZ(imol,iiso,iT))/(Lines(i)%a_therm*sqrt(pi))
			A=A*1d50*mixrat(imol)

c	Random sampling of the Voigt profile
			A=A/real(NV)
			do j=1,NV
				x1=gasdev(idum)/sqrt(2d0)
				x1=x1*Lines(i)%a_therm
				x2=tan((random(idum)-0.5d0)*pi)
				x2=x2*Lines(i)%a_press
				rr=random(idum)
				if(rr.lt.0.25) then
					x=x1+x2
				else if(rr.lt.0.5) then
					x=-x1+x2
				else if(rr.lt.0.75) then
					x=x1-x2
				else
					x=-x1-x2
				endif
				x=(x+Lines(i)%freq-nu1)*scale
				inu=int(x)
				if(inu.ge.1.and.inu.le.nnu) then
					kline(inu)=kline(inu)+A
				endif
			enddo
c	exact computation of the Voigt profile
c				do inu=1,nnu
c					dnu=Lines(i)%freq-nu(inu)
c					x=dnu/gamma
cc					if(x.lt.500d0*100d0/gamma) then
c						call voigt(500d0,x,V)
c						kline(inu)=kline(inu)+A*V
cc					endif
c				enddo

		endif
	enddo
c Above is crap, just to figure out if stuff works
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	

	kmin=1d200
	kmax=0d0
	do inu=1,nnu	
		if(kline(inu).gt.kmax) kmax=kline(inu)
		if(kline(inu).lt.kmin) kmin=kline(inu)
	enddo
	if(kmin.lt.kmax/1d30) kmin=kmax/1d30
	kmin=log10(kmin)
	kmax=log10(kmax)

	nkdis=1000
	allocate(dis(nkdis))
	allocate(kdis(nkdis))
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
			kappa(i)=kdis(j)
		endif
	enddo

	deallocate(nu)
	deallocate(kline)

	return
	end


