	module EOSdata
	IMPLICIT NONE
	real*8,allocatable :: P_EOS(:),T_EOS(:)
	real*8,allocatable :: S_EOS(:,:,:),nab_EOS(:,:,:)
	real*8,allocatable :: dSdP_EOS(:,:,:),dSdT_EOS(:,:,:)
	integer nP_EOS,nT_EOS
	end module EOSdata

	subroutine InitEOS(dir)
	use EOSdata
	IMPLICIT NONE
	character*500 dir,fileH,fileHe,line
	integer i,j
	real*8 T,P,dum,S,dSdP,dSdT
	
	fileH=trim(dir) // "TABLE_H_CMS_vLS.dat"
	fileHe=trim(dir) // "TABLE_He_CMS_vLS.dat"
	
	open(unit=20,file=fileH,RECL=6000)
	open(unit=21,file=fileHe,RECL=6000)
	i=1
	read(20,*)
	read(20,*)
1	j=1
2	read(20,*,err=3,end=4) T,P,dum,dum,S,dum,dum,dSdP,dSdT
	j=j+1
	goto 2
3	i=i+1
	goto 1
4	nP_EOS=j-1
	nT_EOS=i
	print*,nT_EOS,nP_EOS

	rewind(20)

	allocate(T_EOS(nT_EOS),P_EOS(nP_EOS))
	allocate(S_EOS(nT_EOS,nP_EOS,2))
	allocate(dSdP_EOS(nT_EOS,nP_EOS,2))
	allocate(dSdT_EOS(nT_EOS,nP_EOS,2))
	allocate(nab_EOS(nT_EOS,nP_EOS,2))
	
	read(20,*)
	read(21,*)
	do i=1,nT_EOS
		read(20,*)
		read(21,*) 
		do j=1,nP_EOS
			read(20,*) T_EOS(i),P_EOS(j),
     &			dum,dum,S_EOS(i,j,1),dum,dum,
     &			dSdT_EOS(i,j,1),dSdP_EOS(i,j,1),nab_EOS(i,j,1)
			read(21,*) T_EOS(i),P_EOS(j),
     &			dum,dum,S_EOS(i,j,2),dum,dum,
     &			dSdT_EOS(i,j,2),dSdP_EOS(i,j,2),nab_EOS(i,j,2)
		enddo
	enddo
	close(unit=20)
	close(unit=21)

	return
	end
	
	
	subroutine GetNablaEOS(P,T,H_to_He,nabla)
	use EOSdata
	IMPLICIT NONE
	real*8 P,T,H_to_He,nabla,f(2),nab(2)
	real*8 logP,logT,x,y,wP1,wP2,wT1,wT2
	real*8 logS,dSdP,dSdT
	integer i,j,iP,iT

	logP=log10(P/1d4)
	logT=log10(T)
	
	if(logP.le.P_EOS(1)) then
		iP=1
		wP1=1d0
		wP2=0d0
	else if(logP.ge.P_EOS(nP_EOS)) then
		iP=nP_EOS-1
		wP1=0d0
		wP2=1d0
	else
		call hunt(P_EOS,nP_EOS,logP,iP)
		if(iP.lt.1) then
			iP=1
			wP1=1d0
			wP2=0d0
		else if(iP.ge.nP_EOS) then
			iP=nP_EOS-1
			wP1=0d0
			wP2=1d0
		else
			wP1=1d0-(logP-P_EOS(iP))/(P_EOS(iP+1)-P_EOS(iP))
			wP2=1d0-wP1
		endif
	endif

	if(logT.le.T_EOS(1)) then
		iT=1
		wT1=1d0
		wT2=0d0
	else if(logT.ge.T_EOS(nT_EOS)) then
		iT=nT_EOS-1
		wT1=0d0
		wT2=1d0
	else
		call hunt(T_EOS,nT_EOS,logT,iT)
		if(iT.lt.1) then
			iT=1
			wT1=1d0
			wT2=0d0
		else if(iT.ge.nT_EOS) then
			iT=nT_EOS-1
			wT1=0d0
			wT2=1d0
		else
			wT1=1d0-(logT-T_EOS(iT))/(T_EOS(iT+1)-T_EOS(iT))
			wT2=1d0-wT1
		endif
	endif
 
 	x=0d0
 	y=0d0
 	f(1)=H_to_He
 	f(2)=1d0-H_to_He
 	do i=1,2
		logS=S_EOS(iT,iP,i)*wT1*wP1
		logS=logS+S_EOS(iT+1,iP,i)*wT2*wP1
		logS=logS+S_EOS(iT,iP+1,i)*wT1*wP2
		logS=logS+S_EOS(iT+1,iP+1,i)*wT2*wP2
		dSdP=dSdP_EOS(iT,iP,i)*wT1*wP1
		dSdP=dSdP+dSdP_EOS(iT+1,iP,i)*wT2*wP1
		dSdP=dSdP+dSdP_EOS(iT,iP+1,i)*wT1*wP2
		dSdP=dSdP+dSdP_EOS(iT+1,iP+1,i)*wT2*wP2
		dSdT=dSdT_EOS(iT,iP,i)*wT1*wP1
		dSdT=dSdT+dSdT_EOS(iT+1,iP,i)*wT2*wP1
		dSdT=dSdT+dSdT_EOS(iT,iP+1,i)*wT1*wP2
		dSdT=dSdT+dSdT_EOS(iT+1,iP+1,i)*wT2*wP2
		x=x+f(i)*10.0**logS*dSdP
		y=y+f(i)*10.0**logS*dSdT
		nab(i)=nab_EOS(iT,iP,i)*wT1*wP1+
     &		nab_EOS(iT+1,iP,i)*wT2*wP1+
     &		nab_EOS(iT,iP+1,i)*wT1*wP2+
     &		nab_EOS(iT+1,iP+1,i)*wT2*wP2
	enddo
	nabla=-x/y
	if(.not.nabla.lt.maxval(nab)) nabla=maxval(nab)
	if(.not.nabla.gt.minval(nab)) nabla=minval(nab)

	
	return
	end
	
