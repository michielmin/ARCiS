	subroutine TrendCompute(Mp,Rp,Tp,Tstar,Zstar,idum,COratio,NOratio,SiOratio,Z,TiScale)
c	Mp,Rp,Tp,Zstar,idum				are input
c	COratio,NOratio,SiOratio,Z		are output
c	molfracs_atoms 					is set
	use AtomsModule
	IMPLICIT NONE
	real*8 COratio,Z,TiScale,fdry,fwet,Z0,COratio0
	logical enhancecarbon
	real*8 NOratio,SiOratio
	integer i,idum,N
	real*8 Mp,Rp,Zstar,gasdev,Tp,random,Nscale,Tstar

	idum=sign(idum,-1)

	TiScale=1d0
	enhancecarbon=.false.
	COratio0=0.0002478241d0/0.0004509658d0

	Z=0d0
	fdry=1d0
	fwet=1d0

c ========================================================================================
c This is for now: no trend inserted, just solar composition
c ========================================================================================
		TiScale=1d0
		if(Tp.lt.1800d0) TiScale=1d-8

		call set_molfracs_atoms_old(COratio0,Z,TiScale,enhancecarbon)

		COratio= molfracs_atoms(3)/molfracs_atoms(5)
		NOratio= molfracs_atoms(4)/molfracs_atoms(5)
		SiOratio=molfracs_atoms(9)/molfracs_atoms(5)

		return
c ========================================================================================
c ========================================================================================

	COratio=min(COratio0*max(0.5d0,Mp**0.25),1d0)
	Z=Zstar-1d0+min(3d0,(0.75*Rp/Mp)*(1d0/(1d0+exp(-(Tstar-1000d0)/1000d0))))
	fdry=min(1d0,1d0/Mp**0.25)*min(Rp/Mp,1d0)**0.25
	fwet=min(1d0,1d0/Mp**0.25)
	Nscale=(2d0/(1d0+exp((Tp+gasdev(idum)*50d0-750d0)/500d0)))

	Z=Z+gasdev(idum)*0.05d0
	fdry=max(min(fdry*(1d0+gasdev(idum)*0.05),1d0),0d0)
	fwet=max(min(fwet*(1d0+gasdev(idum)*0.05),1d0),0d0)
	COratio=max(0.1d0,COratio+gasdev(idum)*0.05)

	call set_molfracs_atoms_old(COratio0,0d0,1d0,.false.)
	Z0=sum(molfracs_atoms(3:N_atoms))/sum(molfracs_atoms(1:2))
	
	call set_molfracs_atoms_old(COratio,Z,TiScale,enhancecarbon)
	call removecloud(fdry,fwet)
	molfracs_atoms(4)=molfracs_atoms(4)*Nscale
	molfracs_atoms=molfracs_atoms/sum(molfracs_atoms)
	
	COratio= molfracs_atoms(3)/molfracs_atoms(5)
	NOratio= molfracs_atoms(4)/molfracs_atoms(5)
	SiOratio=molfracs_atoms(9)/molfracs_atoms(5)
	Z=-log10(Z0/(sum(molfracs_atoms(3:N_atoms))/sum(molfracs_atoms(1:2))))

	return
	end


	subroutine removecloud(f_dry,f_wet)
	use AtomsModule
	IMPLICIT NONE
	real*8 COabun,f_dry,f_wet,f
	integer i,iCS,nCS,k
	logical,allocatable :: ice(:)
	character*100,allocatable :: CSname(:)
	real*8,allocatable :: atoms_cloud(:,:),xv_bot(:)

	COabun=min(molfracs_atoms(3),molfracs_atoms(5))
	molfracs_atoms(3)=molfracs_atoms(3)-COabun
	molfracs_atoms(5)=molfracs_atoms(5)-COabun

	nCS=10
	allocate(ice(nCS),CSname(nCS),atoms_cloud(nCS,N_atoms),xv_bot(nCS))

	atoms_cloud=0
	i=0
c TiO2
	i=i+1
	CSname(i)='TiO2'
	atoms_cloud(i,15)=1
	atoms_cloud(i,5)=2
	ice(i)=.false.
c VO
	i=i+1
	CSname(i)='VO'
	atoms_cloud(i,16)=1
	atoms_cloud(i,5)=1
	ice(i)=.false.
c Al2O3
	i=i+1
	CSname(i)='Al2O3'
	atoms_cloud(i,8)=2
	atoms_cloud(i,5)=3
	ice(i)=.false.
c SiO2
	i=i+1
	CSname(i)='SiO2'
	atoms_cloud(i,9)=1
	atoms_cloud(i,5)=2
	ice(i)=.false.
c Silicates
	i=i+1
	atoms_cloud(i,9)=1
c	atoms_cloud(i,6)=min(molfracs_atoms(6),molfracs_atoms(8))/molfracs_atoms(9)
c	atoms_cloud(i,8)=min(molfracs_atoms(6),molfracs_atoms(8))/molfracs_atoms(9)
	atoms_cloud(i,7)=molfracs_atoms(7)/molfracs_atoms(9)
c	atoms_cloud(i,13)=molfracs_atoms(13)/molfracs_atoms(9)
c	atoms_cloud(i,14)=molfracs_atoms(14)/molfracs_atoms(9)
	atoms_cloud(i,5)=atoms_cloud(i,6)+atoms_cloud(i,7)+atoms_cloud(i,8)+atoms_cloud(i,14)+atoms_cloud(i,13)+2d0
c	write(CSname(i),'("Al",f3.1,"Na",f3.1,"Mg",f3.1,"SiO",f3.1)') atoms_cloud(i,8),atoms_cloud(i,6),atoms_cloud(i,7),atoms_cloud(i,5)
	write(CSname(i),'("Mg",f3.1,"SiO",f3.1)') atoms_cloud(i,7),atoms_cloud(i,5)
	atoms_cloud(i,9)=0
	atoms_cloud(i,5)=atoms_cloud(i,5)-2
	ice(i)=.false.
c H2O
	i=i+1
	CSname(i)='H2O'
	atoms_cloud(i,1)=2
	atoms_cloud(i,5)=1
	ice(i)=.true.
c Fe
	i=i+1
	CSname(i)='Fe'
	atoms_cloud(i,17)=1
	ice(i)=.false.
c FeS
	i=i+1
	CSname(i)='FeS'
c	atoms_cloud(i,17)=1
	atoms_cloud(i,11)=1
	ice(i)=.false.
c C
	i=i+1
	CSname(i)='C'
	atoms_cloud(i,3)=1
	ice(i)=.false.
c SiC
	i=i+1
	CSname(i)='SiC'
	atoms_cloud(i,9)=1
c	atoms_cloud(i,3)=1
	ice(i)=.false.

	nCS=i

	xv_bot=1d200
	do iCS=1,nCS
		do k=1,N_atoms
			if(atoms_cloud(iCS,k).gt.0) then
				f=molfracs_atoms(k)/atoms_cloud(iCS,k)
				if(f.lt.xv_bot(iCS)) then
					xv_bot(iCS)=f
				endif
			endif
		enddo
		do k=1,N_atoms
			molfracs_atoms(k)=molfracs_atoms(k)-xv_bot(iCS)*atoms_cloud(iCS,k)
			if(molfracs_atoms(k).lt.0d0) then
				molfracs_atoms(k)=0d0
			endif
		enddo
	enddo

	do iCS=1,nCS
		if(ice(iCS)) then
			xv_bot(iCS)=xv_bot(iCS)*f_wet
		else
			xv_bot(iCS)=xv_bot(iCS)*f_dry
		endif
	enddo

	molfracs_atoms(3)=molfracs_atoms(3)+COabun
	molfracs_atoms(5)=molfracs_atoms(5)+COabun
	do iCS=1,nCS
		do k=1,N_atoms
			molfracs_atoms(k)=molfracs_atoms(k)+xv_bot(iCS)*atoms_cloud(iCS,k)
		enddo
	enddo

	deallocate(ice,CSname,atoms_cloud,xv_bot)

	return
	end




