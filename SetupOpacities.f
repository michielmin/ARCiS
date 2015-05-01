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
	nu1=900d0*clight
	nu2=1100d0*clight
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
	real*8 Eu,El,A,x,kmax,kmin,V
	integer nnu,inu,iT,imol,i,ju,jl,j,nkdis
	
	nnu=1000d0*(nu2/clight-nu1/clight)
	allocate(nu(nnu))
	allocate(kline(nnu))
	do inu=1,nnu
		nu(inu)=nu1+real(inu-1)*(nu2-nu1)/real(nnu-1)
	enddo

	kline=0d0
	do imol=1,nmol
		call hunt(Mol(imol)%T,Mol(imol)%nT,Temp,iT)

c This part needs some serious attention!!!
		w=sqrt(2d0*kb*Temp/(Mol(imol)%M*mp))
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
				do inu=1,nnu
					dnu=Mol(imol)%L(i)%freq-nu(inu)
					x=dnu**2/gamma**2
					call voigt(5000d0,x,V)
					kline(inu)=kline(inu)+A*V
				enddo
			endif
		enddo
c Above is crap, just to figure out if stuff works

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
	do i=1,nkdis
		write(90,*) kdis(i),dis(i)
	enddo

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


