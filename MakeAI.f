	subroutine MakeAI()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,k
	real*8 var(n_ret),dvar(2,n_ret),chi2,random
	character*500 outputdir0,command
	logical exist,saneplanet
	
	write(outputdir0,'(a)') trim(outputdir)
	do i=1,nai
		call output("Model number: " // int2string(i,'(i0.6)'))
1		chi2=0d0
		modelsucces=.true.
		do j=1,n_ret
			var(j)=random(idum)
		enddo

		write(outputdir,'(a,"model",i0.6,"/")') trim(outputdir0),i
		write(command,'("mkdir -p ",a)') trim(outputdir)
		call system(command)
		inquire(file=trim(outputdir) // "emis",exist=exist)
		if(.not.exist) then
			call MapRetrieval(var,dvar)
			call WriteRetrieval(i,chi2,var,dvar)
			call InitDens()
			call CheckPlanet(saneplanet)
			if(.not.saneplanet) then
				call output("This is an insane planet...")
				call output("Radius: " // dbl2string(Rplanet/Rjup,'(es10.4)'))
				call output("Mass:   " // dbl2string(Mplanet/Mjup,'(es10.4)'))
				goto 1
			endif
			call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam,starfile)
			Fstar=Fstar*pi*Rstar**2
c			call SetOutputMode(.false.)
			call ComputeModel(.true.)
			if(PTchemAbun) then
				do k=1,n_ret
					do j=1,nmol
						if(RetPar(k)%keyword.eq.molname(j)) then
							RetPar(k)%value=mixrat(j)
						endif
					enddo
				enddo
				call MapRetrievalInverse(var)
				do k=1,n_ret
					if(var(k).lt.0d0) var(k)=0d0
					if(var(k).gt.1d0) var(k)=1d0
				enddo
			else if(.not.dochemistry) then
				do k=1,n_ret
					if(RetPar(k)%keyword.eq.'COratio'.or.RetPar(k)%keyword.eq.'coratio') then
						RetPar(k)%value=COratio
					endif
					if(RetPar(k)%keyword.eq.'metallicity') then
						RetPar(k)%value=metallicity
					endif
				enddo
				call MapRetrievalInverse(var)
			endif
			if(modelsucces) then
				call SetOutputMode(.true.)
				call WriteRetrieval(i,chi2,var,dvar)
				call WriteStructure()
				call WriteOutput()
			else
				call SetOutputMode(.true.)
				call output("something is wrong...")
				call output("try different set of parameters")
				goto 1
			endif
		else
			chi2=random(idum)
		endif
	enddo		
	
	return
	end


	subroutine CheckPlanet(saneplanet)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	logical saneplanet
	real*8 RHill
	integer i

	saneplanet=.true.
c	call SetupStructure(.true.)
c	do i=1,nr
c		if(T(i).le.Tmin) saneplanet=.false.
c	enddo
	
	return
	end
	
	