	subroutine DiffuseCloud(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8,allocatable :: x(:),vsed(:),z(:),cond(:),xtot(:)
	real*8,allocatable :: Nseeds(:)
	real*8,allocatable :: A(:,:),Ainv(:,:),ALU(:,:),y(:,:)
	real*8 K,dz,f,z12,z13,z12_2,z13_2,g,mumol,Pv,Atherm,Btherm,frac,rhodust,rr
	integer ii,info,i,j,iter,NN,NRHS,niter
	real*8 cs,err,maxerr,eps
	integer,allocatable :: IWORK(:)
	character*500 tmp,input

	eps=1d-3
	
	niter=10000

	mumol=86.3d0
	rhodust=2.8d0
	
	Atherm=6.89d4
	Btherm=37.8
	frac=1d-2

	NN=3*nr
	allocate(z(nr))
	allocate(Nseeds(nr))
	allocate(x(NN))
	allocate(xtot(NN))
	allocate(vsed(nr))
	allocate(cond(NN))
	allocate(A(NN,NN))
	allocate(Ainv(NN,NN))
	allocate(ALU(NN,NN))
	allocate(IWORK(10*NN*NN))

	z=R
	K=Cloud(ii)%Kzz

	Cloud(ii)%rv=0.01d0*1d-4	! cm

	Nseeds=Ndens*8.20622804366359e-08/(1000d0**3)

	f=0.5
	x=0d0
	xtot=0d0

	do iter=1,niter
	call tellertje(iter,niter)

	vsed=0d0
	do i=1,nr
		cs=sqrt(kb*T(i)/(2.3*mp))
		vsed(i)=Cloud(ii)%rv(i)*rhodust*Ggrav*Mplanet/(dens(i)*cs*Rplanet**2)
	enddo

	A=0d0
	do i=2,nr-1
		if(i.eq.1) then
			dz=(z(i+1)-z(i))
		else if(i.eq.nr) then
			dz=(z(i)-z(i-1))
		else
			dz=z(i+1)-z(i-1)
		endif

c	quadratic derivative
		z12=(z(i-1)-z(i))
		z13=(z(i-1)-z(i+1))
		z12_2=(z(i-1)**2-z(i)**2)
		z13_2=(z(i-1)**2-z(i+1)**2)
		g=((z13)-((z12)*(z13_2))/(z12_2))
		A(i,i-1)=A(i,i-1)-dens(i-1)*vsed(i-1)*((2d0*z(i)/z12_2)*(1d0-(1d0-z13_2/z12_2)*z12/g) + (1d0-z13_2/z12_2)/g)
		A(i,i)=A(i,i)-dens(i)*vsed(i)*((2d0*z(i)/z12_2)*(-1d0-z13_2*z12/(z12_2*g)) + z13_2/(z12_2*g))
		A(i,i+1)=A(i,i+1)-dens(i+1)*vsed(i+1)*((2d0*z(i)/z12_2)*z12/g - 1d0/g)

		A(i,i-1)=A(i,i-1)-dens(i)*K*2d0*(1d0-((1d0-z13_2/z12_2)/g)*z12)/z12_2
		A(i,i)=A(i,i)-dens(i)*K*2d0*(-1d0-((z13_2/z12_2)/g)*z12)/z12_2
		A(i,i+1)=A(i,i+1)-dens(i)*K*(2d0/g)*z12/z12_2

		A(i,i-1)=A(i,i-1)-((dens(i+1)-dens(i-1))/dz)*K*((2d0*z(i)/z12_2)*(1d0-(1d0-z13_2/z12_2)*z12/g) + (1d0-z13_2/z12_2)/g)
		A(i,i)=A(i,i)-((dens(i+1)-dens(i-1))/dz)*K*((2d0*z(i)/z12_2)*(-1d0-z13_2*z12/(z12_2*g)) + z13_2/(z12_2*g))
		A(i,i+1)=A(i,i+1)-((dens(i+1)-dens(i-1))/dz)*K*((2d0*z(i)/z12_2)*z12/g - 1d0/g)

		j=i+nr

		A(j,j-1)=A(j,j-1)-dens(i)*K*2d0*(1d0-((1d0-z13_2/z12_2)/g)*z12)/z12_2
		A(j,j)=A(j,j)-dens(i)*K*2d0*(-1d0-((z13_2/z12_2)/g)*z12)/z12_2
		A(j,j+1)=A(j,j+1)-dens(i)*K*(2d0/g)*z12/z12_2

		A(j,j-1)=A(j,j-1)-((dens(i+1)-dens(i-1))/dz)*K*((2d0*z(i)/z12_2)*(1d0-(1d0-z13_2/z12_2)*z12/g) + (1d0-z13_2/z12_2)/g)
		A(j,j)=A(j,j)-((dens(i+1)-dens(i-1))/dz)*K*((2d0*z(i)/z12_2)*(-1d0-z13_2*z12/(z12_2*g)) + z13_2/(z12_2*g))
		A(j,j+1)=A(j,j+1)-((dens(i+1)-dens(i-1))/dz)*K*((2d0*z(i)/z12_2)*z12/g - 1d0/g)

		j=i+nr*2

		A(j,j-1)=A(j,j-1)-dens(i-1)*vsed(i-1)*((2d0*z(i)/z12_2)*(1d0-(1d0-z13_2/z12_2)*z12/g) + (1d0-z13_2/z12_2)/g)
		A(j,j)=A(j,j)-dens(i)*vsed(i)*((2d0*z(i)/z12_2)*(-1d0-z13_2*z12/(z12_2*g)) + z13_2/(z12_2*g))
		A(j,j+1)=A(j,j+1)-dens(i+1)*vsed(i+1)*((2d0*z(i)/z12_2)*z12/g - 1d0/g)

		A(j,j-1)=A(j,j-1)-dens(i)*K*2d0*(1d0-((1d0-z13_2/z12_2)/g)*z12)/z12_2
		A(j,j)=A(j,j)-dens(i)*K*2d0*(-1d0-((z13_2/z12_2)/g)*z12)/z12_2
		A(j,j+1)=A(j,j+1)-dens(i)*K*(2d0/g)*z12/z12_2

		A(j,j-1)=A(j,j-1)-((dens(i+1)-dens(i-1))/dz)*K*((2d0*z(i)/z12_2)*(1d0-(1d0-z13_2/z12_2)*z12/g) + (1d0-z13_2/z12_2)/g)
		A(j,j)=A(j,j)-((dens(i+1)-dens(i-1))/dz)*K*((2d0*z(i)/z12_2)*(-1d0-z13_2*z12/(z12_2*g)) + z13_2/(z12_2*g))
		A(j,j+1)=A(j,j+1)-((dens(i+1)-dens(i-1))/dz)*K*((2d0*z(i)/z12_2)*z12/g - 1d0/g)


c	linear derivative
c		A(i,i-1)=A(i,i-1)+dens(i-1)*vsed(i-1)/dz
c		A(i,i+1)=A(i,i+1)-dens(i+1)*vsed(i+1)/dz

c		A(i,i-1)=A(i,i-1)-0.5d0*(dens(i-1)+dens(i))*K/(z(i)-z(i-1))/dz
c		A(i,i)=A(i,i)+0.5d0*(dens(i-1)+dens(i))*K/(z(i)-z(i-1))/dz
c		A(i,i+1)=A(i,i+1)-0.5d0*(dens(i+1)+dens(i))*K/(z(i+1)-z(i))/dz
c		A(i,i)=A(i,i)+0.5d0*(dens(i+1)+dens(i))*K/(z(i+1)-z(i))/dz

c		j=i+nr

c		A(j,j-1)=A(j,j-1)-0.5d0*(dens(i-1)+dens(i))*K/(z(i)-z(i-1))/dz
c		A(j,j)=A(j,j)+0.5d0*(dens(i-1)+dens(i))*K/(z(i)-z(i-1))/dz
c		A(j,j+1)=A(j,j+1)-0.5d0*(dens(i+1)+dens(i))*K/(z(i+1)-z(i))/dz
c		A(j,j)=A(j,j)+0.5d0*(dens(i+1)+dens(i))*K/(z(i+1)-z(i))/dz

c		j=i+nr*2

c		A(j,j-1)=A(j,j-1)+dens(i-1)*vsed(i-1)/dz
c		A(j,j+1)=A(j,j+1)-dens(i+1)*vsed(i+1)/dz
c		A(j,j-1)=A(j,j-1)-0.5d0*(dens(i-1)+dens(i))*K/(z(i)-z(i-1))/dz
c		A(j,j)=A(j,j)+0.5d0*(dens(i-1)+dens(i))*K/(z(i)-z(i-1))/dz
c		A(j,j+1)=A(j,j+1)-0.5d0*(dens(i+1)+dens(i))*K/(z(i+1)-z(i))/dz
c		A(j,j)=A(j,j)+0.5d0*(dens(i+1)+dens(i))*K/(z(i+1)-z(i))/dz
		

c sources and sinks
		j=nr+i

		Pv=exp(-Atherm/T(i)+Btherm)*1d6
		A(i,i)=A(i,i)+dens(i)*(3d0*Pv*mp*mumol/(rhodust*Cloud(ii)%rv(i)*sqrt(2d0*pi*mp*mumol*kb*T(i))))
		A(j,i)=A(j,i)-dens(i)*(3d0*Pv*mp*mumol/(rhodust*Cloud(ii)%rv(i)*sqrt(2d0*pi*mp*mumol*kb*T(i))))

		A(i,j)=A(i,j)-dens(i)*(4d0*pi*Cloud(ii)%rv(i)**2*Nseeds(i)*Rgas*T(i)/(sqrt(2d0*pi*mp*mumol*kb*T(i))*Avogadro))
		A(j,j)=A(j,j)+dens(i)*(4d0*pi*Cloud(ii)%rv(i)**2*Nseeds(i)*Rgas*T(i)/(sqrt(2d0*pi*mp*mumol*kb*T(i))*Avogadro))

	enddo
	
c boundary conditions

	cond=0d0
	cond(nr+1)=frac
	cond(nr*2)=0d0!frac
	cond(nr*2+1)=Ndens(1)*8.20622804366359e-08/(1000d0**3)/dens(1)
	cond(nr*3)=0d0!Ndens(nr)*8.20622804366359e-08/(1000d0**3)/dens(nr)

c	cond(nr*2+1:nr*3-1)=Ndens(1)*8.20622804366359e-08/(1000d0**3)/dens(1)/1e11


	A(1,1)=1d0
	i=nr
	A(i,i-1)=dens(i)*vsed(i)-K*dens(i)/(z(i)-z(i-1))
	A(i,i)=K*dens(i)/(z(i)-z(i-1))

	A(nr+1,nr+1)=1d0
	j=nr*2
	i=nr
	A(j,j-1)=-K*dens(i)/(z(i)-z(i-1))
	A(j,j)=K*dens(i)/(z(i)-z(i-1))

	A(nr*2+1,nr*2+1)=1d0
	j=nr*3
	i=nr
	A(j,j-1)=dens(i)*vsed(i)-K*dens(i)/(z(i)-z(i-1))
	A(j,j)=K*dens(i)/(z(i)-z(i-1))

	NRHS=1
	call DGESV( NN, NRHS, A, NN, IWORK, cond, NN, info )
	x=cond
	
	do i=1,NN
		if(.not.x(i).gt.0d0) x(i)=0d0
		xtot(i)=xtot(i)*(1d0-f)+x(i)*f
		if(iter.eq.1) xtot(i)=x(i)
	enddo
	do i=1,nr
		Nseeds(i)=dens(i)*xtot(nr*2+i)
	enddo
	do i=1,nr
		rr=(3d0*xtot(i)*dens(i)/(4d0*pi*Nseeds(i)*rhodust))**(1d0/3d0)
		if(.not.rr.gt.1d-6) rr=1d-6
		Cloud(ii)%rv(i)=10d0**(log10(Cloud(ii)%rv(i))*(1d0-f)+log10(rr)*f)
		if(iter.eq.1) Cloud(ii)%rv(i)=rr
	enddo
	maxerr=0d0
	do i=1,NN
		err=log10(xtot(i)/x(i))/log10(xtot(i)*x(i))
		if(err.gt.maxerr) maxerr=err
	enddo
	if(iter.gt.10.and.maxerr.lt.eps) exit
c	xtot=x

	enddo
	print*,iter

	cloud_dens(1:nr,ii)=xtot(1:nr)*dens(1:nr)

	deallocate(z)
	deallocate(Nseeds)
	deallocate(x)
	deallocate(xtot)
	deallocate(vsed)
	deallocate(cond)
	deallocate(A)
	deallocate(Ainv)
	deallocate(ALU)
	deallocate(IWORK)

	return
	end
	
	
