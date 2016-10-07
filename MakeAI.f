	subroutine MakeAI()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j
	real*8 var(n_ret),dvar(2,n_ret),chi2,random
	character*500 outputdir0,command
	logical exist

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
			call InitDens()
			call ReadKurucz(Tstar,logg,1d4*lam,Fstar,nlam,starfile)
			Fstar=Fstar*pi*Rstar**2
			call SetOutputMode(.false.)
			call ComputeModel(.true.)
			if(modelsucces) then
				call SetOutputMode(.true.)
				call WriteRetrieval(i,chi2,var,dvar)
				call WriteStructure()
				call WriteOutput()
			else
				call output("chemistry is taking too long...")
				call output("try different set of parameters")
				goto 1
			endif
		else
			chi2=random(idum)
		endif
	enddo		
	
	return
	end
	