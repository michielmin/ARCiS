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
			if(mapCOratio.and..not.dochemistry) call DoMapCOratio()
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
	RHill=(Dplanet*(Mplanet/(3d0*Mstar))**(1d0/3d0))
	if(Rplanet.gt.RHill) saneplanet=.false.
	
c	call SetupStructure(.true.)
c	do i=1,nr
c		if(T(i).le.Tmin) saneplanet=.false.
c	enddo
	
	return
	end



	subroutine DoMapCOratio()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 random,COrat,tot,met,scale
	real*8 f(nmol),Ctot,Otot,Htot,f0(nmol),eps,prev,eps0
	integer i,n,imol
	
	COrat=-1d0
	f=0d0
	eps=10d0
	prev=10d0
	f0=f
	eps0=1d-4
	do while(eps.gt.eps0)
		Otot=0d0
		Ctot=0d0
		Htot=0d0
		tot=0d0
		do imol=1,nmol
			if(molname(imol).ne.'H2'.and.molname(imol).ne.'He') then
				f(imol)=f(imol)+10d0**(10d0*random(idum)-10d0)!random(idum)
			else if(molname(imol).eq.'H2') then
				f(imol)=0.85d0
			else if(molname(imol).eq.'He') then
				f(imol)=0.15d0
			endif
			if(includemol(imol)) then
				tot=tot+f(imol)
			endif
		enddo
		f=f/tot
		do imol=1,nmol
			if(includemol(imol)) then
				Otot=Otot+f(imol)*real(Oatoms(imol))
				Ctot=Ctot+f(imol)*real(Catoms(imol))
				Htot=Htot+f(imol)*real(Hatoms(imol))
			endif
		enddo
		COrat=Ctot/Otot
		eps=abs((COrat-COratio)/(COrat+COratio))
		if(eps.lt.prev) then
			f0=f
			prev=eps
		else
			f=f0
		endif
	enddo
	Otot=Otot/tot
	Ctot=Ctot/tot
	Htot=Htot/tot

	do imol=1,nmol
		if(molname(imol).ne.'H2'.and.molname(imol).ne.'He'
     &				.and.Oatoms(imol).eq.0.and.Catoms(imol).eq.0) then
			f(imol)=10d0**(10d0*random(idum)-10d0)
		else if(molname(imol).eq.'H2') then
			f(imol)=0.85d0
		else if(molname(imol).eq.'He') then
			f(imol)=0.15d0
		endif
		if(includemol(imol)) then
			Otot=Otot+f(imol)*real(Oatoms(imol))
			Ctot=Ctot+f(imol)*real(Catoms(imol))
			Htot=Htot+f(imol)*real(Hatoms(imol))
		endif
	enddo

	scale=10d0**(metallicity)*(Htot*0.0002478241/0.9207539305)/Ctot
	tot=0d0
	do imol=1,nmol
		if(includemol(imol)) then
			if(molname(imol).ne.'H2'.and.molname(imol).ne.'He') then
				f(imol)=f(imol)*scale
			endif
			tot=tot+f(imol)
		else
			f(imol)=-1d-8
		endif
	enddo
	f=f/tot

	do i=1,nr
		mixrat_r(i,1:nmol)=f(1:nmol)
	enddo
	mixrat(1:nmol)=f(1:nmol)
	
	return
	end
		
