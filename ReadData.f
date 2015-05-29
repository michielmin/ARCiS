	subroutine ReadData()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	character*12 imol,iiso,nu,S,A,gamma_air,gamma_self,E,n,delta,gu,gl
	character*100 dummy,homedir,file
	logical exist,doneHITEMP(48)
	integer i,j,k,maxiiso,it,ifile,nfile,nHITEMP,nHITRAN
	type(Line),pointer :: L
	real*8 scale,f_numin,f_numax
c H2O, CO2, CO, NO, OH
	parameter(nfile=34+20+1+1+1)
	character*30 files(nfile)
	parameter(files = (/'01_0-50_HITEMP2010.par','01_1000-1150_HITEMP2010.par','01_11000-30000_HITEMP2010.par',
     &	'01_1150-1300_HITEMP2010.par','01_1300-1500_HITEMP2010.par','01_150-250_HITEMP2010.par',
     &	'01_1500-1750_HITEMP2010.par','01_1750-2000_HITEMP2010.par','01_2000-2250_HITEMP2010.par',
     &	'01_2250-2500_HITEMP2010.par','01_250-350_HITEMP2010.par','01_2500-2750_HITEMP2010.par',
     &	'01_2750-3000_HITEMP2010.par','01_3000-3250_HITEMP2010.par','01_3250-3500_HITEMP2010.par',
     &	'01_350-500_HITEMP2010.par','01_3500-4150_HITEMP2010.par','01_4150-4500_HITEMP2010.par',
     &	'01_4500-5000_HITEMP2010.par','01_50-150_HITEMP2010.par','01_500-600_HITEMP2010.par',
     &	'01_5000-5500_HITEMP2010.par','01_5500-6000_HITEMP2010.par','01_600-700_HITEMP2010.par',
     &	'01_6000-6500_HITEMP2010.par','01_6500-7000_HITEMP2010.par','01_700-800_HITEMP2010.par',
     &	'01_7000-7500_HITEMP2010.par','01_7500-8000_HITEMP2010.par','01_800-900_HITEMP2010.par',
     &	'01_8000-8500_HITEMP2010.par','01_8500-9000_HITEMP2010.par','01_900-1000_HITEMP2010.par',
     &	'01_9000-11000_HITEMP2010.par','02_0-500_HITEMP2010.par','02_1000-1500_HITEMP2010.par',
     &	'02_1500-2000_HITEMP2010.par','02_2000-2125_HITEMP2010.par','02_2125-2250_HITEMP2010.par',
     &	'02_2250-2500_HITEMP2010.par','02_2500-3000_HITEMP2010.par','02_3000-3250_HITEMP2010.par',
     &	'02_3250-3500_HITEMP2010.par','02_3500-3750_HITEMP2010.par','02_3750-4000_HITEMP2010.par',
     &	'02_4000-4500_HITEMP2010.par','02_4500-5000_HITEMP2010.par','02_500-625_HITEMP2010.par',
     &	'02_5000-5500_HITEMP2010.par','02_5500-6000_HITEMP2010.par','02_6000-6500_HITEMP2010.par',
     &	'02_625-750_HITEMP2010.par','02_6500-12785_HITEMP2010.par','02_750-1000_HITEMP2010.par',
     &	'05_HITEMP2010new.par','08_HITEMP2010.par','13_HITEMP2010.par' /))
	integer ind1,ind2,ind3

	doneHITEMP=.false.
	nlines=0
	if(HITEMP) then
		call output("Counting HITEMP database")
		do ifile=1,nfile
			file=trim(HITEMPdir) // trim(files(ifile))
			read(files(ifile)(1:2),*) j
			doneHITEMP(j)=.true.
			if(j.le.nmol) then
			if(includemol(j)) then
				ind1=index(file,'_')
				ind2=index(file,'-')
				ind3=index(file(ind1+1:len_trim(file)),'_')+ind1
				if(ind2.ne.0) then
					read(file(ind1+1:ind2-1),*) f_numin
					read(file(ind2+1:ind3-1),*) f_numax
				else
					f_numin=0d0
					f_numax=1d200
				endif
				if(freq(1).ge.f_numin.and.freq(nlam).le.f_numax) then
					call output("file: " // trim(file))
					open(unit=30,file=file,RECL=500)

1					read(30,'(a2)',end=2) imol
					if(j.le.nmol) then
						if(includemol(j)) nlines=nlines+1
					endif
					goto 1
2					close(unit=30)
					call output("number of lines so far: " // trim(dbl2string(dble(nlines),'(es7.1)')))
				endif
			endif
			endif
		enddo
	endif
	nHITEMP=nlines

	call output("Counting HITRAN database")
	call output("file: " // trim(HITRANdir) // "HITRAN2012.par")

	open(unit=30,file=trim(HITRANdir) // "HITRAN2012.par",RECL=500)

3	read(30,'(a2)',end=4) imol
	read(imol,*) j
	if(j.le.nmol) then
		if(includemol(j).and..not.doneHITEMP(j)) nlines=nlines+1
	endif
	goto 3
4	close(unit=30)

	nHITRAN=nlines-nHITEMP

	allocate(Lines(nlines+1))

	call output("number of lines: " // trim(dbl2string(dble(nlines),'(es7.1)')))

c done counting, now read it in!

	call output("Reading database")

	i=1
	maxiiso=0
	L => Lines(1)
	if(HITEMP) then
		do ifile=1,nfile
			file=trim(HITEMPdir) // trim(files(ifile))
			read(files(ifile)(1:2),*) j
			doneHITEMP(j)=.true.
			if(j.le.nmol) then
			if(includemol(j)) then
				ind1=index(file,'_')
				ind2=index(file,'-')
				ind3=index(file(ind1+1:len_trim(file)),'_')+ind1
				if(ind2.ne.0) then
					read(file(ind1+1:ind2-1),*) f_numin
					read(file(ind2+1:ind3-1),*) f_numax
				else
					f_numin=0d0
					f_numax=1d200
				endif
				if(freq(1).ge.f_numin.and.freq(nlam).le.f_numax) then
					open(unit=30,file=file,RECL=500)

5					read(30,'(i2,i1,f12.0,f10.0,f10.0,f5.0,f5.0,f10.0,f4.0,a87,f7.0,f7.0)',end=6) 
     &					L%imol,L%iiso,L%freq,L%S0,L%Aul,L%gamma_air,L%gamma_self,L%Elow,L%n,dummy,L%gu,L%gl
					j=L%imol
					if(j.le.nmol) then
						if(includemol(j)) then
							call tellertje(i,nlines)
							if(L%imol.gt.nmol) nmol=L%imol
							if(L%iiso.gt.maxiiso) maxiiso=L%iiso
							i=i+1
							L => Lines(i)
						endif
					endif
					goto 5
6					close(unit=30)
				endif
			endif
			endif
		enddo
	endif

	if(nHITRAN.gt.0) then
		open(unit=30,file=trim(HITRANdir) // "HITRAN2012.par",RECL=500)
7		read(30,'(i2,i1,f12.0,f10.0,f10.0,f5.0,f5.0,f10.0,f4.0,a87,f7.0,f7.0)',end=8) 
     &			L%imol,L%iiso,L%freq,L%S0,L%Aul,L%gamma_air,L%gamma_self,L%Elow,L%n,dummy,L%gu,L%gl
		j=L%imol
		if(j.le.nmol) then
			if(includemol(j).and..not.doneHITEMP(j)) then
				call tellertje(i,nlines)
				if(L%imol.gt.nmol) nmol=L%imol
				if(L%iiso.gt.maxiiso) maxiiso=L%iiso
				i=i+1
				L => Lines(i)
			endif
		endif
		goto 7
8		close(unit=30)
	endif

	allocate(niso(nmol))
	niso=0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i)
!$OMP& SHARED(nlines,Lines,niso)
!$OMP DO
	do i=1,nlines
		if(Lines(i)%iiso.gt.niso(Lines(i)%imol)) niso(Lines(i)%imol)=Lines(i)%iiso
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	call output("Partition functions")

	nTZ=5000
	allocate(ZZ(nmol,maxiiso,nTZ))
	allocate(TZ(nTZ))
	ZZ=0d0
	do it=1,nTZ
		TZ(it)=exp(log(71d0)+log(2900d0/71d0)*real(it-1)/real(nTZ-1))
	enddo
	do j=1,nmol
		do k=1,niso(j)
			call TIPS_2011(j,k,296d0,scale)
			do it=1,nTZ
				call TIPS_2011(j,k,TZ(it),ZZ(j,k,it))
				ZZ(j,k,it)=ZZ(j,k,it)/scale
			enddo
		enddo
	enddo

	if(ncia.gt.0) call output("Reading CIA opacities")
	do i=1,ncia
		call InitCIA(i)
		call output("CIA: " // trim(molname(CIA(i)%imol1)) // "-" // trim(molname(CIA(i)%imol2)))
	enddo

	cia_mixrat=-1d0
	do i=1,nmol
		cia_mixrat(i)=mixrat(i)
	enddo
c add Helium (arbitrary value for now...)
	cia_mixrat(48)=0.1
c set default for H2 to 1.0
	if(cia_mixrat(45).lt.0d0) cia_mixrat(45)=1d0
	
	return
	end
	