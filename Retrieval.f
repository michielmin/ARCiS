	subroutine ReadObs()
	use GlobalSetup
	IMPLICIT NONE
	integer i,j
	real*8 x,y,dy
	
	do i=1,nobs
		open(unit=20,file=ObsSpec(i)%file,RECL=1000)
		j=1
1		read(20,*,end=2,err=1) x,y,dy
		j=j+1
		goto 1
2		ObsSpec(i)%nlam=j-1
		rewind(20)
		allocate(ObsSpec(i)%lam(ObsSpec(i)%nlam))
		allocate(ObsSpec(i)%y(ObsSpec(i)%nlam))
		allocate(ObsSpec(i)%dy(ObsSpec(i)%nlam))
		do j=1,ObsSpec(i)%nlam
3			read(20,*,err=3) ObsSpec(i)%lam(j),ObsSpec(i)%y(j),ObsSpec(i)%dy(j)
			ObsSpec(i)%lam(j)=ObsSpec(i)%lam(j)*1d-4
		enddo
	enddo
	
	return
	end



	subroutine GeneticRetrieval()
	use GlobalSetup
	IMPLICIT NONE
	external ComputeChi2
	real*8 var0(n_ret),dvar0(2,n_ret),ComputeChi2
	
	var0=0.5d0
	dvar0=10d0
	open(unit=31,file=trim(outputdir) // "Wolk.dat",RECL=6000)

	call Genetic(ComputeChi2,var0,dvar0,n_ret,nobs,npop,ngen,idum,gene_cross)	

	close(unit=31)
	
	return
	end
	
	
	real*8 function ComputeChi2(imodel,nvars,var,nobs0,chi2obs)
	use GlobalSetup
	use ReadKeywords
	use Constants
	IMPLICIT NONE
	type(SettingKey) key
	character*1000 readline
	integer nvars,nobs0,i,j,imodel
	real*8 var(nvars),chi2obs(nobs0)
	real*8 lamobs(nlam-1)
	real*8,allocatable :: spec(:)
	
	Rplanet=Rplanet/Rjup
	Mplanet=Mplanet/Mjup
	Rstar=Rstar/Rsun
	Dplanet=Dplanet/AU
	lam1=lam1/micron
	lam2=lam2/micron
	distance=distance/parsec
	do i=1,n_ret
		call MapRetrieval(var)
		readline=trim(RetPar(i)%keyword) // "=" // trim(dbl2string(RetPar(i)%value,'(es14.7)'))
		call get_key_value(readline,key%key,key%key1,key%key2,key%value,key%nr1,key%nr2)
		call ReadAndSetKey(key)
	enddo
	call ConvertUnits()
	
	call InitDens()
	call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam)
	Fstar=Fstar*pi*Rstar**2*pi/3.336e11
	call ComputeModel()
	
	do i=1,nlam-1
		lamobs(i)=sqrt(lam(i)*lam(i+1))
	enddo
	ComputeChi2=0d0
	do i=1,nobs
		allocate(spec(ObsSpec(i)%nlam))
		select case(ObsSpec(i)%type)
			case("trans","transmission")
				call regridarray(lamobs,obsA(0,1:nlam-1)/(pi*Rstar**2),nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
			case("emis","emission")
				call regridarray(lamobs,flux(0,1:nlam-1)/(Fstar(1:nlam-1)*1d23/distance**2),
     &					nlam-1,ObsSpec(i)%lam,spec,ObsSpec(i)%nlam)
		end select
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
	
	
	subroutine MapRetrieval(var)
	use GlobalSetup
	IMPLICIT NONE
	integer i
	real*8 var(n_ret)

	do i=1,n_ret
		if(RetPar(i)%logscale) then
c	log
			RetPar(i)%value=10d0**(log10(RetPar(i)%xmin)+log10(RetPar(i)%xmax/RetPar(i)%xmin)*var(i))
		else
c	linear
			RetPar(i)%value=RetPar(i)%xmin+(RetPar(i)%xmax-RetPar(i)%xmin)*var(i)
		endif
	enddo

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
		write(20,*) trim(RetPar(i)%keyword) // "=" // trim(dbl2string(RetPar(i)%value,'(es14.7)'))
	enddo
	close(unit=20)
		
		
	return
	end

