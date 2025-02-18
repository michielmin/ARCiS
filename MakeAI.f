	subroutine MakeAI()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,k,iobs
	real*8 var(n_ret),dvar(2,n_ret),random
	character*500 outputdir0,command
	logical exist,saneplanet,usefile
	real*8 spectemp(nlam),tot,lnew
	real*8,allocatable :: spec(:)
	
	write(outputdir0,'(a)') trim(outputdir)

	if(pargridfile.ne." ") then
		print*,pargridfile
		open(unit=53,file=pargridfile,FORM="FORMATTED",ACCESS="STREAM")
		usefile=.true.
	else
		usefile=.false.
	endif

	i=0
	do while(i.lt.nai.or.usefile)
		i=i+1
		if(usefile) read(53,*,end=2) var(1:n_ret)
		call output("Model number: " // int2string(i,'(i0.6)'))
		write(outputdir,'(a,"model",i0.6,"/")') trim(outputdir0),i
		write(command,'("mkdir -p ",a)') trim(outputdir)
		call system(command)
		inquire(file=trim(outputdir) // "parameters",exist=exist)

1		modelsucces=.true.

		if(.not.usefile) then
			do j=1,n_ret
				var(j)=random(idum)
			enddo
		endif

		if(.not.exist) then
			if(usefile) then
				call MapRetrievalMN(var,dvar)
			else
				call MapRetrieval(var,dvar)
			endif
			if(mapCOratio.and..not.dochemistry) call DoMapCOratio()
			lnew=0d0
			call WriteAI(i,var,lnew)
			call InitDens()
			call CheckPlanet(saneplanet)
			if(.not.saneplanet) then
				call output("This is an insane planet...")
				call output("Radius: " // dbl2string(Rplanet/Rjup,'(es10.4)'))
				call output("Mass:   " // dbl2string(Mplanet/Mjup,'(es10.4)'))
				if(.not.usefile) goto 1
			endif
			if(standardstarname.ne.' ') call StandardStar(standardstarname,Tstar,Rstar,Mstar)
			call StarSpecSetup(Tstar,logg,1d4*lam,1d4*blam,Fstar,nlam,starfile,blackbodystar)
			Fstar=Fstar*pi*Rstar**2
			call SetOutputMode(.false.)
			call ComputeModel(.true.)
			if(.not.dochemistry) then
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
			if(nobs.gt.0) then
				do iobs=1,nobs
					allocate(spec(ObsSpec(iobs)%ndata))
					call RemapObs(iobs,spec,spectemp)
					deallocate(spec)
				enddo
				tot=0d0
				lnew=0d0
				do iobs=1,nobs
					do j=1,ObsSpec(iobs)%ndata
						tot=tot-log(sqrt(2d0*pi)*ObsSpec(iobs)%dy(j))
						lnew=lnew+((ObsSpec(iobs)%model(j)-ObsSpec(iobs)%y(j))/ObsSpec(iobs)%dy(j))**2
					enddo
				enddo
				if(massprior) then
					lnew=lnew+((Mplanet/Mjup-Mp_prior)/dMp_prior)**2
					tot=tot-log(sqrt(2d0*pi)*dMp_prior)
				endif
				lnew=-lnew/2d0+tot
			endif
			if(modelsucces) then
				call SetOutputMode(.true.)
				call WriteAI(i,var,lnew)
				call WriteStructure()
				call WriteOutput()
			else
				call SetOutputMode(.true.)
				call output("something is wrong...")
				call output("try different set of parameters")
				if(.not.usefile) goto 1
			endif
			do iobs=1,nobs
				select case(ObsSpec(iobs)%type)
					case("trans","transmission","emisr","emisR","emisa","emis","emission","transC","phase",
     &							 "phaser","phaseR","transM","transE")
						open(unit=20,file=trim(outputdir) // "obs" // trim(int2string(iobs,'(i0.3)')),FORM="FORMATTED",ACCESS="STREAM")
						do j=1,ObsSpec(iobs)%ndata
							write(20,*) ObsSpec(iobs)%lam(j)*1d4,ObsSpec(iobs)%model(j)/ObsSpec(iobs)%scale,
     &									ObsSpec(iobs)%scale*ObsSpec(iobs)%y(j),
     &									ObsSpec(iobs)%dy(j)
						enddo
						close(unit=20)
				end select
			enddo
		endif
	enddo
2	continue
	if(usefile) close(unit=53)
	
	return
	end


	subroutine CheckPlanet(saneplanet)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	logical saneplanet,MassRadius
	real*8 RHill,SH,Tirr
	integer i

	saneplanet=MassRadius(Mplanet,Rplanet)
	RHill=(Dplanet*(Mplanet/(3d0*Mstar))**(1d0/3d0))
	if(Rplanet.gt.RHill) saneplanet=.false.
	
	SH=(sqrt(Rstar/(2d0*Dplanet))*Tstar*kb)/((Ggrav*Mplanet/(Rplanet**2))*mp*2.3d0)
	if(SH.gt.Rplanet/5d0) saneplanet=.false.

	return
	end

	logical function MassRadius(Mg,Rg)
	use Constants
	IMPLICIT NONE
	real*8 Mg,Rg,Me,Re,Mb,Rb,s
	
	Rb=12.1d0
	Mb=124d0
	if(Mg/Mearth.lt.Mb) then
		s=0.55
	else
		s=0.1
	endif
	Me=Mg/Mearth
	Re=Rb*(Me/Mb)**s
	Re=Re*Rearth
	MassRadius=.false.
	if(abs(log10(Re/Rg)).lt.0.5d0) MassRadius=.true.
	
	return
	end
	


	subroutine DoMapCOratio()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 random,COrat,tot,met,scale
	real*8 f(nmol),Ctot,Otot,Htot,f0(nmol),eps,prev,eps0,minf
	integer i,n,imol
	
	minf=20d0
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
c				f(imol)=f(imol)+10d0**(10d0*random(idum)-10d0)!random(idum)
				f(imol)=(f(imol)-random(idum)*minf)/2d0
			else if(molname(imol).eq.'H2') then
				f(imol)=log10(0.85d0)
			else if(molname(imol).eq.'He') then
				f(imol)=log10(0.15d0)
			endif
			if(includemol(imol)) then
				tot=tot+10d0**f(imol)
			endif
		enddo
c		f=f/tot
		do imol=1,nmol
			if(includemol(imol)) then
				Otot=Otot+10d0**f(imol)*real(Oatoms(imol))/tot
				Ctot=Ctot+10d0**f(imol)*real(Catoms(imol))/tot
				Htot=Htot+10d0**f(imol)*real(Hatoms(imol))/tot
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
	Otot=0d0
	Ctot=0d0
	Htot=0d0

	f=10d0**f

	do imol=1,nmol
		if(molname(imol).ne.'H2'.and.molname(imol).ne.'He'
     &				.and.Oatoms(imol).eq.0.and.Catoms(imol).eq.0) then
			f(imol)=10d0**(minf*random(idum)-minf)
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
			f(imol)=-1d-20
		endif
	enddo
	f=f/tot

	do i=1,nr
		mixrat_r(i,1:nmol)=f(1:nmol)
	enddo
	mixrat(1:nmol)=f(1:nmol)
	
	return
	end



	subroutine WriteAI(imodel,var,like)
	use GlobalSetup
	IMPLICIT NONE
	real*8 var(n_ret),error(2,n_ret),like
	integer i,imodel,j
	character*500 command

c	call MapRetrieval(var,error)

98	open(unit=20,file=trim(outputdir) // "parameters",FORM="FORMATTED",ERR=99,ACCESS="STREAM")
	goto 100
99	write(command,'("mkdir -p ",a)') trim(outputdir)
	call system(command)
	call sleep(10)
	goto 98
100	continue
	write(20,'("Model ",i7)') imodel
	do i=1,n_ret
		if(RetPar(i)%logscale) then
c	log
			write(20,'(a15," = ",es14.7)') trim(RetPar(i)%keyword),RetPar(i)%value
 		else
c	linear/squared
			write(20,'(a15," = ",es14.7)') trim(RetPar(i)%keyword),RetPar(i)%value
		endif
	enddo
	if(.not.dochemistry.or.planetform) then
		write(20,'(a15," = ",es14.7)') 'COratio',COratio
		write(20,'(a15," = ",es14.7)') 'metallicity',metallicity
	endif
	if(mapCOratio) then
		do i=1,nmol
			if(mixrat_r(1,i).gt.0d0) write(20,'(a15," = ",es14.7)') trim(molname(i)),mixrat_r(1,i)
		enddo
	endif
	if(nobs.gt.0.and.like.ne.0d0) write(20,'(a15," = ",es14.7)') 'likelyhood',like

	close(unit=20)

	return
	end

		
