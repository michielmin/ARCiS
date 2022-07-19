c===================================================================================
c First the modules that are also available in ARCiS
c===================================================================================

c===================================================================================
c Setup the names of the atoms and set the molfracs_atoms to Solar abundances
c===================================================================================
	subroutine SetupAtoms
	use AtomsModule
	IMPLICIT NONE

	names_atoms(1) = 'H'
	names_atoms(2) = 'He'
	names_atoms(3) = 'C'
	names_atoms(4) = 'N'
	names_atoms(5) = 'O'
	names_atoms(6) = 'Na'
	names_atoms(7) = 'Mg'
	names_atoms(8) = 'Al'
	names_atoms(9) = 'Si'
	names_atoms(10) = 'P'
	names_atoms(11) = 'S'
	names_atoms(12) = 'Cl'
	names_atoms(13) = 'K'
	names_atoms(14) = 'Ca'
	names_atoms(15) = 'Ti'
	names_atoms(16) = 'V'
	names_atoms(17) = 'Fe'
	names_atoms(18) = 'Ni'

	mass_atoms = (/ 1.00784,
     &	4.002602,
     &	12.0107,
     &	14.0067,
     &	15.999,
     &	22.989769,
     &	24.305,
     &	26.981539,
     &	28.0855,
     &	30.973762,
     &	32.065,
     &	35.453,
     &	39.0983,
     &	40.078,
     &	47.867,
     &	50.9415,
     &	55.845,
     &	58.6934
     &  /)

	molfracs_atoms = (/ 0.9207539305,
     &	0.0783688694,
     &	0.0002478241,
     &	6.22506056949881e-05,
     &	0.0004509658,
     &	1.60008694353205e-06,
     &	3.66558742055362e-05,
     &	2.595e-06,
     &	2.9795e-05,
     &	2.36670201997668e-07,
     &	1.2137900734604e-05,
     &	2.91167958499589e-07,
     &	9.86605611925677e-08,
     &	2.01439011429255e-06,
     &	8.20622804366359e-08,
     &	7.83688694089992e-09,
     &	2.91167958499589e-05,
     &	1.52807116806281e-06
     &  /)


	return
	end
c===================================================================================
c===================================================================================

c===================================================================================
c Subroutine to setup the properties of the pp disk in which the planet forms
c===================================================================================
	subroutine SetupPPdisk(frac_SolidC,Macc_in)
	use FormationModule
	use Constants
	use AtomsModule
	IMPLICIT NONE
	real*8 el_abun(100),gasmass,dustmass,planetmass,const
	real*8 const2,r1,r2,frac_SolidC,Macc_in
	integer i,j,k
	mu=2.5
	d2g_T=0.01
	kappa_r=508d0	! cm^2/g
	alpha_disk=0.1d0
	M_acc=Macc_in*Msun/year

c Compute the abudances in gas phase, solid phase and planetessimals
c	...
	Nr_disk=15

c dtg precalculated, Nr_disk = regions
	if(.not.allocated(R_disk)) then
		allocate(R_disk(Nr_disk),T_disk(Nr_disk),d2g_disk(Nr_disk),p2g_disk(Nr_disk))
		allocate(abun_gas(N_atoms,Nr_disk),abun_dust(N_atoms,Nr_disk),abun_planet(N_atoms,Nr_disk))
	endif
	abun_dust=0d0
	abun_planet=0d0
	abun_gas=0d0

	T_disk=(/3,20,47,120,131,500,650,800,850,1150,1200,1250,1400,1550,1600/)
c   Temperatures adjusted to be just below the temperature boundaries to ensure the composition is for the region beyond that.
	T_disk=T_disk - 1d-4

	do i=1,Nr_disk
		call SplitAtoms(T_disk(i),frac_SolidC)
		abun_gas(1:N_atoms,i)=molfracs_atoms(1:N_atoms)
	enddo
	do i=1,Nr_disk-1
		abun_dust(1:N_atoms,i)=abun_gas(1:N_atoms,Nr_disk)-abun_gas(1:N_atoms,i)
		abun_planet(1:N_atoms,i)=abun_dust(1:N_atoms,i)
	enddo
	do i=1,Nr_disk
		gasmass=sum(abun_gas(1:N_atoms,i)*mass_atoms(1:N_atoms))
		dustmass=sum(abun_dust(1:N_atoms,i)*mass_atoms(1:N_atoms))
		planetmass=sum(abun_planet(1:N_atoms,i)*mass_atoms(1:N_atoms))
		d2g_disk(i)=dustmass/gasmass
		p2g_disk(i)=planetmass/gasmass
		
		abun_gas(1:N_atoms,i)=abun_gas(1:N_atoms,i)*mass_atoms(1:N_atoms)/gasmass
		abun_dust(1:N_atoms,i)=abun_dust(1:N_atoms,i)*mass_atoms(1:N_atoms)/dustmass
		abun_planet(1:N_atoms,i)=abun_planet(1:N_atoms,i)*mass_atoms(1:N_atoms)/planetmass

		do j=1,N_atoms
			if(.not.abun_dust(j,i).ge.0d0) abun_dust(j,i)=abun_gas(j,i)
			if(.not.abun_planet(j,i).ge.0d0) abun_planet(j,i)=abun_gas(j,i)
		enddo

		const = (3.*mu*mp*kappa_r)/(128.*pi**2.*alpha_disk*kb*sigma)
		const = const*d2g_T
      const2 = 150.*(Lstar/Lsun)**(2./7)*(Mstar/Msun)**(-1./7)*(1/AU)**(-3./7)
      r1 = ((const*M_acc**2/T_disk(i)**5)**(2./3)*Ggrav*Mstar)**(1./3)
      r2 = (const2/T_disk(i))**(7./3)
		
		
		
		R_disk(i) = max(r1,r2)
	enddo
		
	return
	end
c===================================================================================
c===================================================================================



c===================================================================================
c Main subroutine computing the formation of the atmosphere during accretion
c===================================================================================
	subroutine Formation(Mplanet,Mcore,Rstart,Rend,f_dust,f_planet,flag_converge)!,Mtot,R,alpha_d)
	use FormationModule
	use Constants
	use AtomsModule
	IMPLICIT NONE
	real*8 Mplanet,Rstart,Rend,f_dust,f_planet,Mcore,R,Rhill,Rmax,Rmin
	real*8 surfacedens 
	real*8 Mtot,eps,mass_gas(Nr_disk),mass_planet(Nr_disk),dMdt,dM,dR,dRdt,dt
	real*8 alpha_d_max,alpha_d_min,alpha_d	
	integer i,j,Nsteps,izone,counter
	logical flag_converge
	parameter(eps=1d-4)

c This is where the magic happens!

	Mtot=0d0
	alpha_d_max = 1d0
	alpha_d_min = 1d-20
	alpha_d = (alpha_d_max*alpha_d_min)**(1./2)
	
	j=0
	counter = 100
	Nsteps=1000
	do while(j.lt.counter)
		Mtot=Mcore
		dR=(Rstart-Rend)/real(Nsteps)
		R=Rstart
		izone=1
		do while(R.lt.R_disk(izone).and.izone.lt.Nr_disk) 
			izone=izone+1
		enddo
		call evolve(R,izone,alpha_d,Mtot,surfacedens,drdt,dMdt,Rhill)
		Rmax=Rstart-Rhill
		do i=1,Nsteps
c			accrete stuff, compute:
c				- Mtot
		
			call evolve(R,izone,alpha_d,Mtot,surfacedens,drdt,dMdt,Rhill)
			dt=dr/drdt
			mass_gas(izone)=mass_gas(izone)+dMdt*dt
			Mtot=Mtot+dMdt*dt
			Rmin=(R)-Rhill
			dM=surfacedens*p2g_disk(izone)*f_planet*max(Rmax**2-Rmin**2,0d0)*pi
			mass_planet(izone)=mass_planet(izone)+dM
			Mtot=Mtot+dM
			Rmax=Rmin
			R=R-dr
			izone=1
			do while(R.lt.R_disk(izone).and.izone.lt.Nr_disk) 
				izone=izone+1
			enddo
		enddo
		if (abs(Mtot-Mplanet)/(Mtot+Mplanet).lt.eps) then
c				alpha_d has a good value
				flag_converge = .true.
				exit
		else if ((j+1).lt.counter) then
			if(Mtot.gt.Mplanet) then
c				Adjust accretion to go faster
				alpha_d_min = alpha_d
				alpha_d = (alpha_d_max*alpha_d_min)**(1./2)
				j=j+1
			else 
c				Adjust accretion to go slower
				alpha_d_max = alpha_d
				alpha_d = (alpha_d_max*alpha_d_min)**(1./2)
				j=j+1
			endif
		else
c			Did not converg!
			flag_converge=.false.
			exit
		endif
	enddo
c	print*,flag_converge,Mtot/Mjup,Mplanet/Mjup,Rstart/AU
	Nsteps=1000
	if (flag_converge) then
		Mtot = Mcore
		mass_gas=0d0
		mass_planet=0d0
		dR=(Rstart-Rend)/real(Nsteps)
		R=Rstart
		izone=1
		do while(R.lt.R_disk(izone).and.izone.lt.Nr_disk) 
			izone=izone+1
		enddo
		call evolve(R,izone,alpha_d,Mtot,surfacedens,drdt,dMdt,Rhill)
		Rmax=Rstart-Rhill
		do i=1,Nsteps
			call evolve(R,izone,alpha_d,Mtot,surfacedens,drdt,dMdt,Rhill)

			dt=dr/drdt
			mass_gas(izone)=mass_gas(izone)+dMdt*dt
			Mtot=Mtot+dMdt*dt

			Rmin=R-Rhill
			dM=surfacedens*p2g_disk(izone)*f_planet*max(Rmax**2-Rmin**2,0d0)*pi
			mass_planet(izone)=mass_planet(izone)+dM
			Mtot=Mtot+dM
			Rmax=Rmin
			R=R-dr
			izone=1
			do while(R.lt.R_disk(izone).and.izone.lt.Nr_disk) 
				izone=izone+1
			enddo
		enddo
		molfracs_atoms=0d0
		do izone=1,Nr_disk
			do i=1,N_atoms
				molfracs_atoms(i)=molfracs_atoms(i)+mass_gas(izone)*(1d0-d2g_disk(izone)*f_dust)*abun_gas(i,izone)
				molfracs_atoms(i)=molfracs_atoms(i)+mass_gas(izone)*d2g_disk(izone)*f_dust*abun_dust(i,izone)
				molfracs_atoms(i)=molfracs_atoms(i)+mass_planet(izone)*abun_planet(i,izone)
			enddo
		enddo
		molfracs_atoms(1:N_atoms)=molfracs_atoms(1:N_atoms)/mass_atoms(1:N_atoms)
		molfracs_atoms(1:N_atoms)=1d12*molfracs_atoms(1:N_atoms)/molfracs_atoms(1)
	endif
	return
	end
c===================================================================================
c===================================================================================

	
	subroutine evolve(R,izone,alpha_d,Mtot,surfacedens,drdt,dMdt,Rhill)
	use FormationModule
	use Constants
	IMPLICIT NONE
	real*8 R,alpha_d,Mtot,surfacedens,drdt,dMdt,Rhill
	real*8 T_local,Cs,scalehight,ang_v,T_irr,T_vis
	real*8 const
	integer izone
	const=(3.*mu*mp*kappa_r*d2g_T)/(128.*pi**(2.)*alpha_disk*kb*sigma)
	T_vis = (const*M_acc**2*(Ggrav*Mstar/R**3)**(3./2))**(1./5)
	T_irr = 150.*(Lstar/Lsun)**(2./7)*(Mstar/Msun)**(-1./7)*(r/AU)**(-3./7)
	T_local =max(T_irr,T_vis)
	Cs =(kb*T_local/(mu*mp))**(1./2)
	ang_v = (Ggrav*Mstar/R**3.)**(1./2)
	scalehight=Cs/ang_v
	const=mu*mp/(3*pi*alpha_disk*kb)
	surfacedens=(const/T_local)*M_acc*(Ggrav*Mstar/R**3)**(1./2)

	Rhill = (Mtot/(3.*Mstar))**(1./3)*R
	drdt = (3.*alpha_d*Cs**2.)/(2.*R*ang_v)
	dMdt = scalehight**2*ang_v*surfacedens*min(3.*pi*alpha_disk,(Rhill/scalehight)**(9./2))

	return
	end

	
c===================================================================================
c===================================================================================

	
	subroutine SplitAtoms(T,frac_solidC)
	use AtomsModule
	IMPLICIT NONE
	real*8 T,limit(N_atoms),lockCO2,lockCO,mol(N_atoms)
	real*8 solidC,frac_solidC,removeH2O
	integer H,He,C,N,O,Na,Mg,Al,Si,P,S,Cl,K,Ca,Ti,V,Fe,Ni
	parameter(H=1)
	parameter(He=2)
	parameter(C=3)
	parameter(N=4)
	parameter(O=5)
	parameter(Na=6)
	parameter(Mg=7)
	parameter(Al=8)
	parameter(Si=9)
	parameter(P=10)
	parameter(S=11)
	parameter(Cl=12)
	parameter(K=13)
	parameter(Ca=14)
	parameter(Ti=15)
	parameter(V=16)
	parameter(Fe=17)
	parameter(Ni=18)
	
c	frac_solidC=0d0
	
c First set Solar abundances
	call SetupAtoms
	molfracs_atoms=1e12*molfracs_atoms/molfracs_atoms(1)

	solidC=frac_solidC*molfracs_atoms(C)
		
c Follow recipe according to Khorshid+ 2022
	if(T.le.850d0) then
		mol=0
		limit=1d0
		mol(Na)=1
		mol(Al)=1
		mol(Si)=3
		mol(O)=8
		call remove(mol,limit)
	endif
	if(T.le.850d0) then
		mol=0
		limit=1d0
		mol(K)=1
		mol(Al)=1
		mol(Si)=3
		mol(O)=8
		call remove(mol,limit)
	endif
	if(T.le.1250d0) then
		mol=0
		limit=1d0
		mol(Mg)=1
		mol(Al)=2
		mol(O)=4
		limit(Al)=0.5
		call remove(mol,limit)
	endif
	if(T.le.1550d0) then
		mol=0
		limit=1d0
		mol(Al)=2
		mol(O)=3
		call remove(mol,limit)
	endif
c Additionally remove TiO2 and VO
	if(T.le.1550d0) then
		mol=0
		limit=1d0
		mol(Ti)=1
		mol(O)=2
		call remove(mol,limit)
	endif
	if(T.le.1550d0) then
		mol=0
		limit=1d0
		mol(V)=1
		mol(O)=1
		call remove(mol,limit)
	endif
	if(T.le.1150d0) then
		mol=0
		limit=1d0
		mol(Mg)=1
		mol(Si)=1
		mol(O)=3
		limit(Si)=0.83
		call remove(mol,limit)
	endif
	if(T.le.1200d0) then
		mol=0
		limit=1d0
		mol(Mg)=2
		mol(Si)=1
		mol(O)=4
		call remove(mol,limit)
	endif
	if(T.le.1400d0) then
		mol=0
		limit=1d0
		mol(Si)=1
		mol(O)=2
		call remove(mol,limit)
	endif
	if(T.le.650d0) then
		mol=0
		limit=1d0
		mol(Fe)=1
		mol(S)=1
		call remove(mol,limit)
	endif
	if(T.le.1200d0) then
		mol=0
		limit=1d0
		mol(Fe)=1
		call remove(mol,limit)
	endif
	
c Soot line: remove solidC fraction of the carbon atoms into soot when T<800K
	if(T.le.800d0) then
c Meteoritic abundances
		if(T.le.500d0) then
			solidC=max(solidC,10.0**7.39)
			molfracs_atoms(6:N_atoms)=0d0
			molfracs_atoms(He)=max(0d0,molfracs_atoms(He)-10.0**1.29)
			molfracs_atoms(N)=max(0d0,molfracs_atoms(N)-10.0**6.26)
			removeH2O=10.0**8.4-(10.0**8.69-molfracs_atoms(O))
			molfracs_atoms(H)=max(0d0,molfracs_atoms(H)-removeH2O*2d0)
		endif
		molfracs_atoms(C)=max(0d0,molfracs_atoms(C)-solidC)
	endif

c Adjust lockCO2 and lockCO according to limiting values to make sure abundances do not drop below 0
	lockCO2=min(0.15*molfracs_atoms(C),2d0*molfracs_atoms(O))
	molfracs_atoms(C)=molfracs_atoms(C)-lockCO2
	molfracs_atoms(O)=molfracs_atoms(O)-2d0*lockCO2

	lockCO=min(molfracs_atoms(C),molfracs_atoms(O))
	molfracs_atoms(C)=molfracs_atoms(C)-lockCO
	molfracs_atoms(O)=molfracs_atoms(O)-lockCO

	if(T.le.120d0) then
		mol=0
		limit=1d0
		mol(H)=2
		mol(O)=1
		call remove(mol,limit)
	endif

c Possibility to add NH3 snowline (this might need to be adjusted for N2 iceline if needed)
	if(T.le.131d0) then
		mol=0
		limit=1d0
		mol(H)=3
		mol(N)=1
		call remove(mol,limit)
	endif

c Add CO and/or CO2 back to the gas phase if temperature is too high
	if(T.gt.47d0) then
		molfracs_atoms(C)=molfracs_atoms(C)+lockCO2
		molfracs_atoms(O)=molfracs_atoms(O)+2d0*lockCO2
	endif
	if(T.gt.20d0) then
		molfracs_atoms(C)=molfracs_atoms(C)+lockCO
		molfracs_atoms(O)=molfracs_atoms(O)+lockCO
	endif
	
	
	return
	end
	
	
	subroutine remove(mol,limit)
	use AtomsModule
	IMPLICIT NONE
	real*8 f,mol(N_atoms),limit(N_atoms)
	
	f=minval(molfracs_atoms*limit/mol)
	molfracs_atoms=molfracs_atoms-mol*f
	
	return
	end
	


