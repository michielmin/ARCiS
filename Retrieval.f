	subroutine ReadObs()
	use GlobalSetup
	IMPLICIT NONE
	integer i,j,ilam,nj
	real*8 x,y,dy
	
	do i=1,nobs
		open(unit=20,file=ObsSpec(i)%file,RECL=1000)
		j=1
		ilam=1
1		read(20,*,end=2,err=1) x,y,dy
		if(x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)) ilam=ilam+1
		j=j+1
		goto 1
2		ObsSpec(i)%nlam=ilam-1
		nj=j-1
		rewind(20)
		allocate(ObsSpec(i)%lam(ObsSpec(i)%nlam))
		allocate(ObsSpec(i)%y(ObsSpec(i)%nlam))
		allocate(ObsSpec(i)%dy(ObsSpec(i)%nlam))
		ilam=1
		do j=1,nj
3			read(20,*,err=3) x,y,dy
			if(x.gt.(lam(1)*1d4).and.x.lt.(lam(nlam)*1d4)) then
				ObsSpec(i)%lam(ilam)=x*1d-4
				ObsSpec(i)%y(ilam)=y
				ObsSpec(i)%dy(ilam)=dy
				ilam=ilam+1
			endif
		enddo
	enddo
	
	return
	end



	subroutine DoRetrieval()
	use GlobalSetup
	IMPLICIT NONE
	external ComputeChi2
	real*8 var0(n_ret),dvar0(2,n_ret),ComputeChi2,dvar(n_ret),x(n_ret)
	real*8,allocatable :: y(:),y1(:),y2(:),dy(:,:),yobs(:),dyobs(:)
	real*8 chi2obs(nobs),var(n_ret),chi2,chi2max,gasdev,maxd,error(n_ret)
	real*8,allocatable :: W(:,:),WS(:)
	real*8 chi2_0,chi2_1,chi2_2,x1,x2,minT(nr),maxT(nr),ran1,tot,W0(n_ret,n_ret)
	integer imodel,ny,i,j,iter,iy,k
	integer MDW,ME,MA,MG,MODE
	integer,allocatable :: IP(:)
	real*8 PRGOPT(10),RNORME,RNORML
	
	var0=0.5d0
	dvar0=10d0
	open(unit=31,file=trim(outputdir) // "Wolk.dat",RECL=6000)

	if(ngen.gt.0) then
c first genetic algoritm to make the first estimate
		call Genetic(ComputeChi2,var0,dvar0,n_ret,nobs,npop,ngen,idum,gene_cross)
	endif

	imodel=npop*ngen
	ny=0
	do i=1,nobs
		ny=ny+ObsSpec(i)%nlam
	enddo
	allocate(y(ny))
	allocate(y1(ny))
	allocate(y2(ny))
	allocate(dy(n_ret,ny))
	allocate(yobs(ny))
	allocate(dyobs(ny))
	MDW=ny+2*n_ret
	ME=0
	MA=ny
	MG=n_ret*2
	allocate(W(MDW,n_ret+1))
	allocate(IP(10*(MG+2*n_ret+2)))
	allocate(WS((2*(ME+n_ret)+max(MA+MG,n_ret)+(MG+2)*(n_ret+7))*10))
	iy=1
	do i=1,nobs
		yobs(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%y(1:ObsSpec(i)%nlam)
		dyobs(iy:iy+ObsSpec(i)%nlam-1)=ObsSpec(i)%dy(1:ObsSpec(i)%nlam)
		iy=iy+ObsSpec(i)%nlam
	enddo
	dvar=0.1d0
	var=var0
	chi2_0=1d200
	chi2_1=1d200
	chi2_2=1d200
	chi2max=0d0
	do iter=1,100
		call output("Iteration number" // trim(int2string(iter,'(i4)')))
		j=0
1		continue
		j=j+1
		imodel=imodel+1
		chi2=ComputeChi2(imodel,n_ret,var,nobs,chi2obs,error)
		if(chi2.ge.chi2max) then
			var0=var
			chi2max=chi2
			call output("Updating best fit")
			call WriteStructure()
			call WriteOutput()
			call WriteRetrieval(imodel,1d0/chi2,var(1:n_ret))
			call WritePTlimits(var0,W0)
		else if(j.gt.3.or.iter.eq.1) then
			var0=var
		else
			dvar=dvar/2d0
			var=var0+dvar
			goto 1
		endif
		chi2_2=1d0/chi2
c		if(abs((chi2_0-chi2_2)/(chi2_0+chi2_2)).lt.1d-4) exit
		chi2_0=chi2_1
		chi2_1=chi2_2
		iy=1
		do i=1,nobs
	 		call RemapObs(i,y(iy:iy+ObsSpec(i)%nlam-1))
			iy=iy+ObsSpec(i)%nlam
		enddo
     	dvar=dvar/2d0
		do i=1,n_ret
			if(dvar(i).lt.1d-2) dvar(i)=1d-2
			var=var0
			var(i)=var(i)+dvar(i)
			if(var(i).gt.1d0) var(i)=1d0
			if(var(i).lt.0d0) var(i)=0d0
			x1=var(i)
			imodel=imodel+1
			chi2=ComputeChi2(imodel,n_ret,var,nobs,chi2obs,error)
			if(chi2.gt.chi2max) then
				chi2max=chi2
				call output("Updating best fit")
				call WriteStructure()
				call WriteOutput()
				call WriteRetrieval(imodel,1d0/chi2,var(1:n_ret))
				call WritePTlimits(var,W0)
			endif
			iy=1
			do j=1,nobs
		 		call RemapObs(j,y1(iy:iy+ObsSpec(j)%nlam-1))
				iy=iy+ObsSpec(j)%nlam
			enddo
			var=var0
			var(i)=var(i)-dvar(i)
			if(var(i).gt.1d0) var(i)=1d0
			if(var(i).lt.0d0) var(i)=0d0
			x2=var(i)
			imodel=imodel+1
			chi2=ComputeChi2(imodel,n_ret,var,nobs,chi2obs,error)
			if(chi2.gt.chi2max) then
				chi2max=chi2
				call output("Updating best fit")
				call WriteStructure()
				call WriteOutput()
				call WriteRetrieval(imodel,1d0/chi2,var(1:n_ret))
				call WritePTlimits(var,W0)
			endif
			iy=1
			do j=1,nobs
		 		call RemapObs(j,y2(iy:iy+ObsSpec(j)%nlam-1))
				iy=iy+ObsSpec(j)%nlam
			enddo
			dy(i,1:ny)=(y1(1:ny)-y2(1:ny))/abs(x1-x2)
		enddo
		iy=0
		do i=1,n_ret
			W(1:ny,i)=dy(i,1:ny)/dyobs(1:ny)
			iy=iy+1
			W(ny+iy,1:n_ret)=0d0
			W(ny+iy,i)=-1d0
			W(ny+iy,n_ret+1)=var0(i)-1d0
			iy=iy+1
			W(ny+iy,1:n_ret)=0d0
			W(ny+iy,i)=1d0
			W(ny+iy,n_ret+1)=-var0(i)
		enddo
		W(1:ny,n_ret+1)=(yobs(1:ny)-y(1:ny))/dyobs(1:ny)		

		IP(1)=(2*(ME+n_ret)+max(MA+MG,n_ret)+(MG+2)*(n_ret+7))*10
		IP(2)=(MG+2*n_ret+2)*10
		PRGOPT(1)=2
		PRGOPT(2)=1
		PRGOPT(3)=1
		PRGOPT(4)=1
		MODE=0
		call DLSEI (W, MDW, ME, MA, MG, n_ret, PRGOPT, dvar, RNORME,
     +   RNORML, MODE, WS, IP)
		var=var0+dvar

		W0(1:n_ret,1:n_ret)=W(1:n_ret,1:n_ret)
		do i=1,n_ret
			if(var(i).lt.0d0) var(i)=0d0
			if(var(i).gt.1d0) var(i)=1d0
			error(i)=sqrt(W(i,i))
			do j=1,n_ret
				if(error(i).lt.sqrt(abs(W(i,j)))) then
					error(i)=sqrt(abs(W(i,j)))
				endif
			enddo
		enddo
	enddo

	close(unit=31)
	
	return
	end



	subroutine WritePTlimits(var0,W)
	use GlobalSetup
	IMPLICIT NONE
	real*8 minT(nr),maxT(nr),var(n_ret),error(n_ret),tot,W(n_ret,n_ret)
	real*8 var0(n_ret),ran1
	integer i,j,k

	maxT=0d0
	minT=1d200
	call SetOutputMode(.false.)
	do k=1,1000
2		continue
		tot=0d0
		do j=1,n_ret
			error(j)=2d0*(ran1(idum)-0.5d0)
			tot=tot+error(j)**2
			if(tot.gt.1d0) goto 2
		enddo
		do i=1,n_ret
			var(i)=var0(i)
			do j=1,n_ret
				var(i)=var(i)+sqrt(abs(W(i,j)))*error(j)
			enddo
			if(var(i).lt.0d0) var(i)=0d0
			if(var(i).gt.1d0) var(i)=1d0
		enddo
		error=0d0
		call MapRetrieval(var,error)
		call InitDens()
		call SetupStructure()
		do i=1,nr
			if(T(i).gt.maxT(i)) maxT(i)=T(i)
			if(T(i).lt.minT(i)) minT(i)=T(i)
		enddo				
	enddo
	call SetOutputMode(.true.)
	open(unit=45,file=trim(outputdir) // "limits.dat")
	do i=1,nr
		write(45,*) P(i),minT(i),maxT(i)
	enddo
	close(unit=45)

	return
	end



	subroutine RemapObs(i,spec)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j
	real*8 spec(*)
	real*8 lamobs(nlam-1)

	do j=1,nlam-1
		lamobs(j)=sqrt(lam(j)*lam(j+1))
	enddo
	select case(ObsSpec(i)%type)
		case("trans","transmission")
			call regridarray(lamobs,obsA(0,1:nlam-1)/(pi*Rstar**2),nlam-1,
     &					ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
		case("emisr","emisR")
			call regridarray(lamobs,flux(0,1:nlam-1)/(Fstar(1:nlam-1)*1d23/distance**2),
     &					nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
		case("emisa","emis","emission")
			call regridarray(lamobs,flux(0,1:nlam-1),
     &					nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
	end select


	return
	end
		
	
	real*8 function ComputeChi2(imodel,nvars,var,nobs0,chi2obs,error0)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nvars,nobs0,i,j,imodel
	real*8 var(nvars),chi2obs(nobs0),error(n_ret)
	real*8,allocatable :: spec(:)
	real*8,intent(in),optional :: error0(n_ret)
	error=0d0
	if(present(error0)) error=error0
	
	call MapRetrieval(var,error)
	
	call InitDens()
	call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam)
	Fstar=Fstar*pi*Rstar**2*pi/3.336e11
	call SetOutputMode(.false.)
	call ComputeModel()
	call SetOutputMode(.true.)
	
	ComputeChi2=0d0
	do i=1,nobs
		allocate(spec(ObsSpec(i)%nlam))
		call RemapObs(i,spec)
		chi2obs(i)=0d0
		do j=1,ObsSpec(i)%nlam
			chi2obs(i)=chi2obs(i)+((spec(j)-ObsSpec(i)%y(j))/ObsSpec(i)%dy(j))**2
		enddo
		chi2obs(i)=chi2obs(i)/real(ObsSpec(i)%nlam)
		ComputeChi2=ComputeChi2+chi2obs(i)
		open(unit=40,file='file' // trim(int2string(i,'(i1)')))
		do j=1,ObsSpec(i)%nlam
			write(40,*) ObsSpec(i)%lam(j),spec(j),ObsSpec(i)%y(j)
		enddo
		close(unit=40)
		deallocate(spec)
	enddo

	write(31,*) imodel,ComputeChi2,(trim(dbl2string(RetPar(i)%value,'(es18.7)')),i=1,n_ret)
	flush(31)

	chi2obs=1d0/chi2obs
	ComputeChi2=1d0/ComputeChi2
	
	return
	end
	
	
	subroutine MapRetrieval(var,dvar)
	use GlobalSetup
	use ReadKeywords
	use Constants
	IMPLICIT NONE
	type(SettingKey) key
	character*1000 readline
	integer i
	real*8 var(n_ret),dvar(n_ret),x

	do i=1,n_ret
		if(RetPar(i)%logscale) then
c	log
			x=var(i)
			RetPar(i)%value=10d0**(log10(RetPar(i)%xmin)+log10(RetPar(i)%xmax/RetPar(i)%xmin)*x)
			x=var(i)+dvar(i)
			RetPar(i)%error=10d0**(log10(RetPar(i)%xmin)+log10(RetPar(i)%xmax/RetPar(i)%xmin)*x)
			RetPar(i)%error=RetPar(i)%error/RetPar(i)%value
		else if(RetPar(i)%squarescale) then
c	square
			x=var(i)
			RetPar(i)%value=sqrt(RetPar(i)%xmin**2+(RetPar(i)%xmax**2-RetPar(i)%xmin**2)*x)
			x=var(i)+dvar(i)
			RetPar(i)%error=sqrt(RetPar(i)%xmin**2+(RetPar(i)%xmax**2-RetPar(i)%xmin**2)*x)
			RetPar(i)%error=RetPar(i)%error-RetPar(i)%value
		else
c	linear
			x=var(i)
			RetPar(i)%value=RetPar(i)%xmin+(RetPar(i)%xmax-RetPar(i)%xmin)*x
			x=var(i)+dvar(i)
			RetPar(i)%error=RetPar(i)%xmin+(RetPar(i)%xmax-RetPar(i)%xmin)*x
			RetPar(i)%error=RetPar(i)%error-RetPar(i)%value
		endif
	enddo

	Rplanet=Rplanet/Rjup
	Mplanet=Mplanet/Mjup
	Rstar=Rstar/Rsun
	Dplanet=Dplanet/AU
	lam1=lam1/micron
	lam2=lam2/micron
	distance=distance/parsec
	do i=1,n_ret
		readline=trim(RetPar(i)%keyword) // "=" // trim(dbl2string(RetPar(i)%value,'(es14.7)'))
		call get_key_value(readline,key%key,key%key1,key%key2,key%value,key%nr1,key%nr2)
		call ReadAndSetKey(key)
	enddo
	call ConvertUnits()

	return
	end
	
	
	subroutine regridarray(x0,y0,n0,x1,y1,n1)
	IMPLICIT NONE
	integer i1,i0,n0,n1
	real*8 x0(n0),y0(n0),x1(n1),y1(n1)

	do i1=1,n1
		if(x1(i1).lt.x0(1)) then
			y1(i1)=y0(1)
		else if(x1(i1).gt.x0(n0)) then
			y1(i1)=y0(n0)
		else
			do i0=1,n0-1
				if(x1(i1).ge.x0(i0).and.x1(i1).le.x0(i0+1)) then
					y1(i1)=y0(i0)+(y0(i0+1)-y0(i0))*(x0(i0+1)-x1(i1))/(x0(i0+1)-x0(i0))
				endif
			enddo
		endif
	enddo

	return
	end
	
	subroutine WriteRetrieval(imodel,chi2,var)
	use GlobalSetup
	IMPLICIT NONE
	real*8 var(n_ret),chi2
	integer i,imodel

	open(unit=20,file=trim(outputdir) // "retrieval",RECL=1000)
	write(20,'("Model ",i)') imodel
	write(20,'("chi2=",f14.6)') chi2
	do i=1,n_ret
		if(RetPar(i)%logscale) then
c	log
			write(20,'(a15," = ",es14.7," x/: ",es11.4)') trim(RetPar(i)%keyword),RetPar(i)%value,RetPar(i)%error
 		else
c	linear/squared
			write(20,'(a15," = ",es14.7," +/- ",es11.4)') trim(RetPar(i)%keyword),RetPar(i)%value,RetPar(i)%error
		endif
	enddo
	close(unit=20)
		
		
	return
	end

