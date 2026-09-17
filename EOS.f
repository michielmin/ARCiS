c===================================================================================
c Test EOS implementation with separate H, He, and metal (AQUA) tables.
c
c The H and He tables use a common grid, while TABLE_AQUA_vLS.dat has its
c own T-P grid.  Each pure-component entropy and its derivatives are
c therefore interpolated on its native grid at the requested state before
c the mixture is formed with
c
c     fH  = X
c     fHe = Y
c     fZ  = Z .
c
c The X, Y, and Z values are the same hardcoded metallicity-dependent mass
c fractions as in EOS.f.  H_to_He is retained in GetNablaEOS so this file
c can be substituted for EOS.f without changing the existing call sites.
c===================================================================================
	module EOSdata_Z_test
	IMPLICIT NONE
	real*8,allocatable :: P_HHe_EOS(:),T_HHe_EOS(:)
	real*8,allocatable :: S_HHe_EOS(:,:,:)
	real*8,allocatable :: nab_HHe_EOS(:,:,:)
	real*8,allocatable :: dSdP_HHe_EOS(:,:,:)
	real*8,allocatable :: dSdT_HHe_EOS(:,:,:)
	real*8,allocatable :: P_Z_EOS(:),T_Z_EOS(:)
	real*8,allocatable :: S_Z_EOS(:,:),nab_Z_EOS(:,:)
	real*8,allocatable :: dSdP_Z_EOS(:,:),dSdT_Z_EOS(:,:)
	integer nP_HHe_EOS,nT_HHe_EOS,nP_Z_EOS,nT_Z_EOS
	end module EOSdata_Z_test


c===================================================================================
c Count the rectangular grid in an EOS table.  The first line is the column
c header and every temperature block starts with another line beginning '#'.
c===================================================================================
	subroutine CountEOSTable(file,nT,nP)
	IMPLICIT NONE
	character*(*) file
	character*6000 line
	integer nT,nP,nrow,ios

	nT=0
	nP=0
	nrow=0
	open(unit=29,file=file,status='old',action='read',iostat=ios,
     &	RECL=6000)
	if(ios.ne.0) then
		write(*,'(A,A)') 'EOS ERROR: cannot open ',trim(file)
		stop 1
	endif

c Skip the global column header.
	read(29,'(A)',iostat=ios) line
	if(ios.ne.0) then
		write(*,'(A,A)') 'EOS ERROR: empty EOS table ',trim(file)
		stop 1
	endif

10	continue
	read(29,'(A)',iostat=ios) line
	if(ios.lt.0) goto 20
	if(ios.gt.0) then
		write(*,'(A,A)') 'EOS ERROR: failed reading ',trim(file)
		stop 1
	endif
	if(len_trim(line).eq.0) goto 10
	if(line(1:1).eq.'#') then
		nT=nT+1
	else
		nrow=nrow+1
	endif
	goto 10

20	continue
	close(unit=29)
	if(nT.lt.2.or.nrow.lt.2.or.mod(nrow,nT).ne.0) then
		write(*,'(A,A)') 'EOS ERROR: non-rectangular table ',
     &		trim(file)
		stop 1
	endif
	nP=nrow/nT
	if(nP.lt.2) then
		write(*,'(A,A)') 'EOS ERROR: EOS pressure grid too small ',
     &		trim(file)
		stop 1
	endif

	return
	end


	subroutine InitEOS(dir)
	use EOSdata_Z_test
	IMPLICIT NONE
	character*500 dir,fileH,fileHe,fileZ
	character*6000 line
	integer i,j,ios,nT_check,nP_check
	real*8 Tval,Pval,Tcheck,Pcheck,dum,tol
	parameter(tol=1d-8)

	fileH=trim(dir) // 'TABLE_H_CMS_vLS.dat'
	fileHe=trim(dir) // 'TABLE_He_CMS_vLS.dat'
	fileZ=trim(dir) // 'TABLE_AQUA_vLS.dat'

	call CountEOSTable(fileH,nT_HHe_EOS,nP_HHe_EOS)
	call CountEOSTable(fileHe,nT_check,nP_check)
	if(nT_check.ne.nT_HHe_EOS.or.nP_check.ne.nP_HHe_EOS) then
		write(*,'(A)')
     &		'EOS ERROR: H and He EOS grids have different sizes.'
		stop 1
	endif
	call CountEOSTable(fileZ,nT_Z_EOS,nP_Z_EOS)

c Permit InitEOS to be called again in a test process.
	if(allocated(P_HHe_EOS)) deallocate(P_HHe_EOS)
	if(allocated(T_HHe_EOS)) deallocate(T_HHe_EOS)
	if(allocated(S_HHe_EOS)) deallocate(S_HHe_EOS)
	if(allocated(nab_HHe_EOS)) deallocate(nab_HHe_EOS)
	if(allocated(dSdP_HHe_EOS)) deallocate(dSdP_HHe_EOS)
	if(allocated(dSdT_HHe_EOS)) deallocate(dSdT_HHe_EOS)
	if(allocated(P_Z_EOS)) deallocate(P_Z_EOS)
	if(allocated(T_Z_EOS)) deallocate(T_Z_EOS)
	if(allocated(S_Z_EOS)) deallocate(S_Z_EOS)
	if(allocated(nab_Z_EOS)) deallocate(nab_Z_EOS)
	if(allocated(dSdP_Z_EOS)) deallocate(dSdP_Z_EOS)
	if(allocated(dSdT_Z_EOS)) deallocate(dSdT_Z_EOS)

	allocate(P_HHe_EOS(nP_HHe_EOS),T_HHe_EOS(nT_HHe_EOS))
	allocate(S_HHe_EOS(nT_HHe_EOS,nP_HHe_EOS,2))
	allocate(nab_HHe_EOS(nT_HHe_EOS,nP_HHe_EOS,2))
	allocate(dSdP_HHe_EOS(nT_HHe_EOS,nP_HHe_EOS,2))
	allocate(dSdT_HHe_EOS(nT_HHe_EOS,nP_HHe_EOS,2))
	allocate(P_Z_EOS(nP_Z_EOS),T_Z_EOS(nT_Z_EOS))
	allocate(S_Z_EOS(nT_Z_EOS,nP_Z_EOS))
	allocate(nab_Z_EOS(nT_Z_EOS,nP_Z_EOS))
	allocate(dSdP_Z_EOS(nT_Z_EOS,nP_Z_EOS))
	allocate(dSdT_Z_EOS(nT_Z_EOS,nP_Z_EOS))

	open(unit=20,file=fileH,status='old',action='read',iostat=ios,
     &	RECL=6000)
	if(ios.ne.0) goto 900
	open(unit=21,file=fileHe,status='old',action='read',iostat=ios,
     &	RECL=6000)
	if(ios.ne.0) goto 900

c Skip the global column headers.
	read(20,'(A)',iostat=ios) line
	if(ios.ne.0) goto 900
	read(21,'(A)',iostat=ios) line
	if(ios.ne.0) goto 900

	do i=1,nT_HHe_EOS
c Skip the temperature-block headers.
		read(20,'(A)',iostat=ios) line
		if(ios.ne.0) goto 900
		read(21,'(A)',iostat=ios) line
		if(ios.ne.0) goto 900
		do j=1,nP_HHe_EOS
			read(20,*,iostat=ios) Tval,Pval,dum,dum,
     &			S_HHe_EOS(i,j,1),dum,dum,
     &			dSdT_HHe_EOS(i,j,1),dSdP_HHe_EOS(i,j,1),
     &			nab_HHe_EOS(i,j,1)
			if(ios.ne.0) goto 900
			read(21,*,iostat=ios) Tcheck,Pcheck,dum,dum,
     &			S_HHe_EOS(i,j,2),dum,dum,
     &			dSdT_HHe_EOS(i,j,2),dSdP_HHe_EOS(i,j,2),
     &			nab_HHe_EOS(i,j,2)
			if(ios.ne.0) goto 900

			if(abs(Tval-Tcheck).gt.tol.or.
     &			abs(Pval-Pcheck).gt.tol) then
				write(*,'(A,2I8)')
     &				'EOS ERROR: H/He grid mismatch at ',i,j
				stop 1
			endif
			if(j.eq.1) then
				T_HHe_EOS(i)=Tval
			else if(abs(Tval-T_HHe_EOS(i)).gt.tol) then
				write(*,'(A)')
     &				'EOS ERROR: inconsistent H/He temperature block.'
				stop 1
			endif
			if(i.eq.1) then
				P_HHe_EOS(j)=Pval
			else if(abs(Pval-P_HHe_EOS(j)).gt.tol) then
				write(*,'(A)')
     &				'EOS ERROR: inconsistent H/He pressure grid.'
				stop 1
			endif
		enddo
	enddo
	close(unit=20)
	close(unit=21)

	open(unit=22,file=fileZ,status='old',action='read',iostat=ios,
     &	RECL=6000)
	if(ios.ne.0) goto 910
	read(22,'(A)',iostat=ios) line
	if(ios.ne.0) goto 910
	do i=1,nT_Z_EOS
		read(22,'(A)',iostat=ios) line
		if(ios.ne.0) goto 910
		do j=1,nP_Z_EOS
			read(22,*,iostat=ios) Tval,Pval,dum,dum,
     &			S_Z_EOS(i,j),dum,dum,dSdT_Z_EOS(i,j),
     &			dSdP_Z_EOS(i,j),nab_Z_EOS(i,j)
			if(ios.ne.0) goto 910
			if(j.eq.1) then
				T_Z_EOS(i)=Tval
			else if(abs(Tval-T_Z_EOS(i)).gt.tol) then
				write(*,'(A)')
     &				'EOS ERROR: inconsistent metal temperature block.'
				stop 1
			endif
			if(i.eq.1) then
				P_Z_EOS(j)=Pval
			else if(abs(Pval-P_Z_EOS(j)).gt.tol) then
				write(*,'(A)')
     &				'EOS ERROR: inconsistent metal pressure grid.'
				stop 1
			endif
		enddo
	enddo
	close(unit=22)

	write(*,'(A,2I7)') 'H/He EOS grid (nT,nP): ',
     &	nT_HHe_EOS,nP_HHe_EOS
	write(*,'(A,2I7)') 'Metal EOS grid (nT,nP): ',
     &	nT_Z_EOS,nP_Z_EOS

	return

900	continue
	write(*,'(A)') 'EOS ERROR: failed while reading H/He EOS tables.'
	stop 1
910	continue
	write(*,'(A)') 'EOS ERROR: failed while reading metal EOS table.'
	stop 1
	end


c===================================================================================
c Return the hardcoded bulk mass fractions for a supported metallicity.
c The nearest table entry is selected, but it must agree within mh_tol.
c===================================================================================
	subroutine GetHardcodedEOSComposition(mh,Xmass,Ymass,Zmass)
	IMPLICIT NONE
	integer,parameter :: ncomp=12
	integer iz
	real*8 mh,Xmass,Ymass,Zmass,mh_tol
	real*8 mh_grid(ncomp),X_grid(ncomp),Y_grid(ncomp),Z_grid(ncomp)
	parameter(mh_tol=1d-6)
	data mh_grid /
     &	-0.5d0,-0.3d0,0.0d0,0.3d0,0.5d0,0.7d0,
     &	 1.0d0, 1.5d0,1.7d0,1.75d0,2.0d0,2.21d0 /
	data Z_grid /
     &	0.0044d0,0.0070d0,0.0139d0,0.0273d0,
     &	0.0426d0,0.0659d0,0.1234d0,0.3081d0,
     &	0.4138d0,0.4419d0,0.5847d0,0.6955d0 /
	data Y_grid /
     &	0.2446d0,0.2440d0,0.2423d0,0.2390d0,
     &	0.2352d0,0.2295d0,0.2154d0,0.1700d0,
     &	0.1440d0,0.1371d0,0.1020d0,0.0748d0 /
	data X_grid /
     &	0.7509d0,0.7490d0,0.7438d0,0.7337d0,
     &	0.7221d0,0.7046d0,0.6612d0,0.5219d0,
     &	0.4422d0,0.4210d0,0.3132d0,0.2297d0 /

	iz=minloc(abs(mh_grid-mh),dim=1)
	if(abs(mh_grid(iz)-mh).gt.mh_tol) then
		write(*,'(A,ES14.6)')
     &		'EOS ERROR: unsupported [M/H] = ',mh
		write(*,'(A)')
     &		'Add this metallicity to the EOS composition table.'
		stop 1
	endif

	Xmass=X_grid(iz)
	Ymass=Y_grid(iz)
	Zmass=Z_grid(iz)

	return
	end


c===================================================================================
c Find a grid cell and linear interpolation weights.  Values outside a
c component table are clamped to that component's nearest grid boundary,
c matching the behavior of EOS.f.
c===================================================================================
	subroutine GetEOSBracket(grid,n,value,idx,w1,w2)
	IMPLICIT NONE
	integer n,idx
	real*8 grid(n),value,w1,w2

	idx=0
	if(value.le.grid(1)) then
		idx=1
		w1=1d0
		w2=0d0
	else if(value.ge.grid(n)) then
		idx=n-1
		w1=0d0
		w2=1d0
	else
		call hunt(grid,n,value,idx)
		if(idx.lt.1) then
			idx=1
			w1=1d0
			w2=0d0
		else if(idx.ge.n) then
			idx=n-1
			w1=0d0
			w2=1d0
		else
			w1=1d0-(value-grid(idx))/(grid(idx+1)-grid(idx))
			w2=1d0-w1
		endif
	endif

	return
	end


c===================================================================================
c Bilinearly interpolate one pure-component EOS table.
c===================================================================================
	subroutine InterpolateEOSComponent(logP,logT,Pgrid,nP,Tgrid,nT,
     &	Stab,dPdtab,dTdtab,nabtab,logS,dSdP,dSdT,nab)
	IMPLICIT NONE
	integer nP,nT,iP,iT
	real*8 logP,logT,Pgrid(nP),Tgrid(nT)
	real*8 Stab(nT,nP),dPdtab(nT,nP),dTdtab(nT,nP)
	real*8 nabtab(nT,nP),logS,dSdP,dSdT,nab
	real*8 wP1,wP2,wT1,wT2

	call GetEOSBracket(Pgrid,nP,logP,iP,wP1,wP2)
	call GetEOSBracket(Tgrid,nT,logT,iT,wT1,wT2)

	logS=Stab(iT,iP)*wT1*wP1+
     &	Stab(iT+1,iP)*wT2*wP1+
     &	Stab(iT,iP+1)*wT1*wP2+
     &	Stab(iT+1,iP+1)*wT2*wP2
	dSdP=dPdtab(iT,iP)*wT1*wP1+
     &	dPdtab(iT+1,iP)*wT2*wP1+
     &	dPdtab(iT,iP+1)*wT1*wP2+
     &	dPdtab(iT+1,iP+1)*wT2*wP2
	dSdT=dTdtab(iT,iP)*wT1*wP1+
     &	dTdtab(iT+1,iP)*wT2*wP1+
     &	dTdtab(iT,iP+1)*wT1*wP2+
     &	dTdtab(iT+1,iP+1)*wT2*wP2
	nab=nabtab(iT,iP)*wT1*wP1+
     &	nabtab(iT+1,iP)*wT2*wP1+
     &	nabtab(iT,iP+1)*wT1*wP2+
     &	nabtab(iT+1,iP+1)*wT2*wP2

	return
	end


	subroutine GetNablaEOS(P,T,ir,nabla)
	use EOSdata_Z_test
	use GlobalSetup,only : nmol, mixrat_r,Mmol,Hatoms,dochemistry
	use AtomsModule
	IMPLICIT NONE
	real*8 P,T,H_to_He,nabla
	real*8 logP,logT,x,y,Sref,fsum
	real*8 Xmass,Ymass,Zmass,f(3),logS(3)
	real*8 dSdP(3),dSdT(3),nab(3)
	integer i,ir

c Outside the lower pressure/temperature bounds of the H/He tables, use
c the adiabatic gradient of an ideal diatomic gas, as in EOS.f.
	if(P.lt.1d-5.or.T.lt.100d0) then
		nabla=2d0/7d0
		return
	endif

c ARCiS supplies pressure in bar.  The EOS pressure coordinate is log10(GPa).
	logP=log10(P/1d4)
	logT=log10(T)

c	call GetHardcodedEOSComposition(metallicity,Xmass,Ymass,Zmass)
	if(dochemistry) then
		Xmass=molfracs_atoms(1)*mass_atoms(1)
		Ymass=molfracs_atoms(2)*mass_atoms(2)
		Zmass=0d0
		do i=3,N_atoms
			Zmass=Zmass+molfracs_atoms(i)*mass_atoms(i)
		enddo
	else
		Xmass=0d0
		Ymass=mixrat_r(ir,48)*Mmol(48)
		Zmass=0d0
		do i=1,nmol
			Xmass=Xmass+mixrat_r(ir,i)*Hatoms(i)*1.0079
			Zmass=Zmass+mixrat_r(ir,i)*(Mmol(i)-Hatoms(i)*1.0079)
		enddo
	endif
	fsum=Xmass+Ymass+Zmass
	if(Xmass.lt.0d0.or.Ymass.lt.0d0.or.Zmass.lt.0d0.or.
     &	fsum.le.0d0) then
		write(*,'(A,3ES14.6)')
     &		'EOS ERROR: invalid X, Y, Z fractions: ',
     &		Xmass,Ymass,Zmass
		stop 1
	endif
	f(1)=Xmass/fsum
	f(2)=Ymass/fsum
	f(3)=Zmass/fsum

c H and He are interpolated on their shared native grid.
	call InterpolateEOSComponent(logP,logT,P_HHe_EOS,
     &	nP_HHe_EOS,T_HHe_EOS,nT_HHe_EOS,S_HHe_EOS(:,:,1),
     &	dSdP_HHe_EOS(:,:,1),dSdT_HHe_EOS(:,:,1),
     &	nab_HHe_EOS(:,:,1),logS(1),dSdP(1),dSdT(1),nab(1))
	call InterpolateEOSComponent(logP,logT,P_HHe_EOS,
     &	nP_HHe_EOS,T_HHe_EOS,nT_HHe_EOS,S_HHe_EOS(:,:,2),
     &	dSdP_HHe_EOS(:,:,2),dSdT_HHe_EOS(:,:,2),
     &	nab_HHe_EOS(:,:,2),logS(2),dSdP(2),dSdT(2),nab(2))

c Metals are independently interpolated on the AQUA table's native grid.
	call InterpolateEOSComponent(logP,logT,P_Z_EOS,nP_Z_EOS,
     &	T_Z_EOS,nT_Z_EOS,S_Z_EOS,dSdP_Z_EOS,dSdT_Z_EOS,
     &	nab_Z_EOS,logS(3),dSdP(3),dSdT(3),nab(3))

c The additive entropy-of-mixing term is independent of P and T and hence
c does not enter these derivatives.  Shift the logarithmic entropies before
c exponentiation to avoid unnecessary overflow or underflow.
	Sref=maxval(logS)
	x=0d0
	y=0d0
	do i=1,3
		x=x+f(i)*10d0**(logS(i)-Sref)*dSdP(i)
		y=y+f(i)*10d0**(logS(i)-Sref)*dSdT(i)
	enddo

	if(abs(y).le.tiny(1d0)) then
		write(*,'(A,2ES14.6)')
     &		'EOS ERROR: zero entropy T derivative at P,T = ',P,T
		stop 1
	endif
	nabla=-x/y

c Preserve the pure-component bounds used by EOS.f, now including Z.
	if(.not.nabla.lt.maxval(nab)) nabla=maxval(nab)
	if(.not.nabla.gt.minval(nab)) nabla=minval(nab)

	!write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)')
    ! &		'New EOS mixed table: P [bar] = ',P,
    ! &		', T [K] = ',T,
    ! &		', nabla_ad = ',nabla
	!call flush(6)

	return
	end
