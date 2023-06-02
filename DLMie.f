	module AImod
	IMPLICIT NONE
	integer nlay,nn
	real*8,allocatable :: w(:,:,:),b(:,:),a(:,:)
!$OMP THREADPRIVATE(a)
	real*8 DLMie_e1max,DLMie_e2max
	real*8 DLMie_xmin,DLMie_xmax
	end module AImod
	
	
	subroutine InitDLMie()
	use AImod
	IMPLICIT NONE
	integer i,j,n1,n2,k,l
	character*1000 file,datadir
	logical check

	call getenv('HOME',datadir)
	datadir = trim(datadir) // "/ARCiS/Data/DLMie/"

	DLMie_e1max=10d0
	DLMie_e2max=10d0
	DLMie_xmin=0.0005
	DLMie_xmax=5000d0
	
	nlay=6
	nn=130
	allocate(w(nn,nn,nlay+1),b(nn,nlay+1))
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(nn,nlay)
	allocate(a(nn,nlay))
!$OMP END PARALLEL
	do i=1,nlay+1
		if(i.eq.1) then
			n1=4
			n2=nn
		else if(i.eq.nlay+1) then
			n1=nn
			n2=2
		else
			n1=nn
			n2=nn
		endif
		write(file,'(a,"layer",i1,"_weightmatrix.dat")') trim(datadir),i+1
		inquire(file=file,exist=check)
		if(.not.check) then
			call output("Datafile for DLMie not found!: " // trim(file))
			stop
		endif
		open(unit=20,file=file,FORM="FORMATTED",ACCESS="STREAM")
		do j=1,n1
			read(20,*) w(1:n2,j,i)
		enddo
		close(unit=20)
		write(file,'(a,"layer",i1,"_biasvector.dat")') trim(datadir),i+1
		inquire(file=file,exist=check)
		if(.not.check) then
			call output("Datafile for DLMie not found!: " // trim(file))
			stop
		endif
		open(unit=20,file=file,FORM="FORMATTED",ACCESS="STREAM")
		do j=1,n2
			read(20,*) b(j,i)
		enddo
		close(unit=20)
	enddo

	return
	end
	
	subroutine DLMie(rad,lam,e1,e2,fmax,Cabs,Csca)
	use AImod
	IMPLICIT NONE
	real*8 rad,lam,e1,e2,fmax,Cabs,Csca,x,pi,scaleA,scaleS
	integer i,j
	parameter(pi=3.1415926536)
	x=2d0*pi*rad/lam
	
	scaleA=1d0
	scaleS=1d0
	if(x.lt.DLMie_xmin) then
		scaleA=x/DLMie_xmin
		scaleS=scaleA**4
		x=DLMie_xmin
	else if(x.gt.DLMie_xmax) then
		x=DLMie_xmax
	endif
	
	a(1:nn,1)=b(1:nn,1)
	a(1:nn,1)=a(1:nn,1)+w(1:nn,1,1)*x
	a(1:nn,1)=a(1:nn,1)+w(1:nn,2,1)*e1
	a(1:nn,1)=a(1:nn,1)+w(1:nn,3,1)*e2
	a(1:nn,1)=a(1:nn,1)+w(1:nn,4,1)*fmax
	do j=1,nn
		a(j,1)=max(0d0,a(j,1))
	enddo
	do i=2,nlay
		a(1:nn,i)=b(1:nn,i)
		do j=1,nn
			a(1:nn,i)=a(1:nn,i)+w(1:nn,j,i)*a(j,i-1)
		enddo
		do j=1,nn
			a(j,i)=max(0d0,a(j,i))
		enddo
	enddo
	Cabs=b(1,nlay+1)
	Csca=b(2,nlay+1)
	do j=1,nn
		Cabs=Cabs+w(1,j,nlay+1)*a(j,nlay)
		Csca=Csca+w(2,j,nlay+1)*a(j,nlay)
	enddo
	Cabs=10.**Cabs
	Csca=10.**Csca
	Cabs=scaleA*Cabs*pi*rad**2
	Csca=scaleS*Csca*pi*rad**2

	return
	end
	