	subroutine ReadLambdaFiles(imol)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,i_low,i_up,imol
	real*8 v1,v2
	
	call output("Reading Lambda file: " // trim(Mol(imol)%filename))
	open(unit=80,file=Mol(imol)%filename,RECL=6000)
	read(80,*)
	read(80,*) Mol(imol)%name
	read(80,*)
	read(80,*) Mol(imol)%M
	read(80,*)
	read(80,*) Mol(imol)%nlevels
	read(80,*)
	allocate(Mol(imol)%E(Mol(imol)%nlevels))
	allocate(Mol(imol)%g(Mol(imol)%nlevels))
	do i=1,Mol(imol)%nlevels
		read(80,*) j,Mol(imol)%E(i),Mol(imol)%g(i)
		Mol(imol)%E(i)=Mol(imol)%E(i)*11600d0/8065.54d0
	enddo
	read(80,*)
	read(80,*) Mol(imol)%nlines
	read(80,*)

	call output("  number of levels: " // int2string(Mol(imol)%nlevels,'(i14)'))
	call output("  number of lines:  " // int2string(Mol(imol)%nlines,'(i14)'))
	allocate(Mol(imol)%L(Mol(imol)%nlines))
	do i=1,Mol(imol)%nlines
		read(80,*) j,Mol(imol)%L(i)%jup,Mol(imol)%L(i)%jlow,Mol(imol)%L(i)%Aul,
     &				Mol(imol)%L(i)%freq,!Mol(imol)%E(Mol(imol)%L(i)%jup)
     &				Mol(imol)%L(i)%Eup

		Mol(imol)%L(i)%freq=Mol(imol)%L(i)%freq*1d9
		Mol(imol)%L(i)%lam=clight*1d4/(Mol(imol)%L(i)%freq)
	enddo
	close(unit=80)

c compute the partition function
	Mol(imol)%nT=500
	allocate(Mol(imol)%Z(Mol(imol)%nT))
	allocate(Mol(imol)%T(Mol(imol)%nT))
	do i=1,Mol(imol)%nT
		Mol(imol)%T(i)=10d0**(4d0*real(i-1)/real(Mol(imol)%nT-1))
		Mol(imol)%Z(i)=0d0
		do j=1,Mol(imol)%nlevels
			Mol(imol)%Z(i)=Mol(imol)%Z(i)+Mol(imol)%g(j)*exp(-Mol(imol)%E(j)/Mol(imol)%T(i))
		enddo
	enddo

	return
	end
	
