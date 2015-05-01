	subroutine SetupOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer imol
	real*8 kappa(100),g(100),nu1,nu2,Temp,dens0
	integer ng
	
	do imol=1,nmol
		select case(Mol(imol)%filetype)
			case("LAMBDA")
				call ReadLambdaFiles(imol)
			case default
				call output("Unknown filetype")
		end select
	enddo

	Temp=800d0
	dens0=1d-10
	nu1=800d0*clight
	nu2=1200d0*clight
	ng=100
	call ComputeKtable(dens0,Temp,nu1,nu2,kappa,g,ng)

	
	return
	end
	
	
	subroutine ComputeKtable(dens0,Temp,nu1,nu2,kappa,g,ng)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ng
	real*8 dens0,Temp,nu1,nu2,kappa(ng),g(ng),w,dnu,gamma
	real*8,allocatable :: nu(:),kline(:),kdis(:),dis(:)
	real*8 Eu,El,A,x,kmax,kmin,V,scale,x1,x2,gasdev,random
	integer nnu,inu,iT,imol,i,ju,jl,j,nkdis,NV,nl,k
	
	nnu=10d0*(nu2/clight-nu1/clight)
	allocate(nu(nnu))
	allocate(kline(nnu))
	do inu=1,nnu
		nu(inu)=nu1+real(inu-1)*(nu2-nu1)/real(nnu-1)
	enddo
	scale=real(nnu)/(nu2-nu1)

	nl=0
	do imol=1,nmol
		do i=1,Mol(imol)%nlines
			if(Mol(imol)%L(i)%freq.gt.nu1.and.Mol(imol)%L(i)%freq.lt.nu2) nl=nl+1
		enddo
	enddo
	NV=nnu*100/nl+100
	print*,nl,NV


	kline=0d0
	do k=1,1000
	do imol=1,nmol
		call hunt(Mol(imol)%T,Mol(imol)%nT,Temp,iT)

c This part needs some serious attention!!!
		w=sqrt(2d0*kb*Temp/(Mol(imol)%M*mp))

!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,gamma,ju,jl,Eu,El,A,x1,x2,x,inu,j,V)
!$OMP& SHARED(Mol,Temp,iT,NV,imol,nu1,nu2,scale,nnu,kline,w)
!$OMP DO SCHEDULE(STATIC,j)
		do i=1,Mol(imol)%nlines
			gamma=w*Mol(imol)%L(i)%freq/clight
			if(Mol(imol)%L(i)%freq.gt.nu1.and.Mol(imol)%L(i)%freq.lt.nu2) then
				ju=Mol(imol)%L(i)%jup
				jl=Mol(imol)%L(i)%jlow
				Eu=Mol(imol)%E(ju)
				El=Mol(imol)%E(jl)
				A=Mol(imol)%g(ju)*Mol(imol)%L(i)%Aul*(exp(-El/Temp)-exp(-Eu/Temp))
				A=A/(Mol(imol)%L(i)%freq**3*Mol(imol)%Z(iT))/(gamma*sqrt(pi))
				A=A*1d50

c	Random sampling of the Voigt profile
				A=A/real(NV)
				do j=1,NV
					x1=gasdev(idum)/sqrt(2d0)
					x2=tan((random(idum)-0.5d0)*pi)
					x=((x1+500*x2)*gamma+Mol(imol)%L(i)%freq-nu1)*scale
					inu=int(x)
					if(inu.ge.1.and.inu.le.nnu) then
						kline(inu)=kline(inu)+A
					endif
				enddo
c	exact computation of the Voigt profile
c				do inu=1,nnu
c					dnu=Mol(imol)%L(i)%freq-nu(inu)
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

	enddo
	enddo
	
	do i=1,nnu
		write(90,*) nu(i),kline(i)
	enddo


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
		g(i)=(real(i)-0.5)/real(ng)
		if(g(i).lt.dis(1)) then
			kappa(i)=kdis(1)
		else
			call hunt(dis,nkdis,g(i),j)
			kappa(i)=kdis(j)
		endif
	enddo

	do i=1,ng
		write(91,*) g(i),kappa(i)
	enddo

	return
	end


