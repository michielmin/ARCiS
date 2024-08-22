	subroutine CondensationCloud(ii)
	use GlobalSetup
	use Constants
	use AtomsModule
	use CloudModule
	use TimingModule
	IMPLICIT NONE
	real*8,allocatable :: x(:),vsed(:),dx(:),vth(:)
	real*8,allocatable :: Sn(:),mpart(:)
	real*8,allocatable :: y(:,:),CloudHp(:),CloudMMW(:)
	real*8,allocatable :: at_ab(:,:)
	real*8,allocatable,save :: Sc(:,:),vthv(:),IWORK(:),AB(:,:)
	real*8,allocatable :: tcinv(:,:),rho_av(:),Kd(:),Kg(:),Km(:),tcrystinv(:)
	real*8 dz,g,rr,mutot,npart,tot,lambda,tot1,tot2,tot3,nucryst,Tcryst
	integer info,i,j,iter,NN,NRHS,niter,ii,k,ihaze,kl,ku,jCS
	real*8 cs,eps,frac_nuc,m_nuc,tcoaginv,Dp,vmol,f,mm,err,maxerr,dztot
	real*8 Pv,molfracs_atoms0(N_atoms),NKn,Kzz_r(nr),vBM,scale
	integer,allocatable :: ixv(:,:),ixc(:,:),iVL(:,:),ixn(:),icryst(:),iamorph(:),inuc(:),ifit(:)
	real*8 sigmastar,Sigmadot,Pstar,sigmamol,COabun,lmfp,fstick,kappa_cloud,fmin,rho_nuc,Gibbs
	logical ini,Tconverged,haze_nuc
	character*500 cloudspecies(max(nclouds,1)),form
	real*8,allocatable :: xn_iter(:,:),xc_iter(:,:,:),xv_iter(:,:,:)
	logical,allocatable :: docondense(:)
	integer iCS,ir,nrdo,iconv,nconv,iVS,nVS,NStot,ik,nfscale
	real*8 logP(nr),logx(nr),dlogx(nr),St,fsed,fscale,pscale,Jn_temp
	real*8,allocatable :: logCloudP(:),CloudtauUV(:),CloudkappaUV(:),CloudG(:)
	character*10,allocatable :: v_names(:),v_names_out(:)
	logical,allocatable :: v_include(:),c_rainout(:),do_nuc(:)
	integer INCFD,IERR
	logical SKIP,liq
	real*8 time,kp,Otot(nr),Ctot(nr),Ntot(nr),compGibbs,ffrag
	integer itime
	real*8,allocatable :: v_atoms(:,:),muC(:),muV(:),v_cloud(:,:),Sat(:,:),Sat0(:,:),fSat(:,:),v_H2(:)
	real*8,allocatable :: xv_out(:),Jn_xv(:,:),A_J(:),B_j(:),sigma_nuc(:),r0_nuc(:),Nf_nuc(:),Nc_nuc(:,:)
	real*8,allocatable :: bv(:,:),bc(:,:),bH2(:),x0(:)

	logical dochemR(nr)

	call cpu_time(time)
	timecloud=timecloud-time
	call system_clock(itime)
	itimecloud=itimecloud-itime
	ctimecloud=ctimecloud+1

	nVS=16
	allocate(v_names(nVS),v_atoms(nVS,N_atoms),v_include(nVS))
	allocate(bv(nVS,0:4),bH2(0:4))
	bv=0d0

	bH2(0)=5.19096E+04
	bH2(1)=-1.80117E+00
	bH2(2)=8.72246E-02
	bH2(3)=2.56139E-04
	bH2(4)=-5.35403E-09
	
	v_atoms=0d0

	v_names(1)="SiO"
	v_atoms(1,9)=1
	v_atoms(1,5)=1
	bv(1,0)=9.58807E+04
	bv(1,1)=-1.26716E+00
	bv(1,2)=-5.98429E+00
	bv(1,3)=2.61720E-04
	bv(1,4)=-1.67404E-08

	v_names(2)="TiO"
	v_atoms(2,5)=1
	v_atoms(2,15)=1
	bv(2,0)=8.02469E+04
	bv(2,1)=-8.46602E-01
	bv(2,2)=-7.67141E+00
	bv(2,3)=1.51196E-04
	bv(2,4)=-2.20473E-08

	v_names(3)="Mg"
	v_atoms(3,7)=1

	v_names(4)="H2O"
	v_atoms(4,1)=2
	v_atoms(4,5)=1
	bv(4,0)=1.10336E+05
	bv(4,1)=-4.17836E+00
	bv(4,2)=3.17447E+00
	bv(4,3)=9.40647E-04
	bv(4,4)=-4.04825E-08

	v_names(5)="H2S"
	v_atoms(5,1)=2
	v_atoms(5,11)=1
	bv(5,0)=8.71685E+04
	bv(5,1)=-4.03150E+00
	bv(5,2)=3.16992E+00
	bv(5,3)=1.08827E-03
	bv(5,4)=-5.46714E-08

	v_names(6)="Fe"
	v_atoms(6,17)=1

	v_names(7)="Al"
	v_atoms(7,8)=1

	v_names(8)="Na"
	v_atoms(8,6)=1

	v_names(9)="K"
	v_atoms(9,13)=1

	v_names(10)="HCl"
	v_atoms(10,1)=1
	v_atoms(10,12)=1
	bv(10,0)=5.13028E+04
	bv(10,1)=-2.13173E+00
	bv(10,2)=2.87709E+00
	bv(10,3)=4.45605E-04
	bv(10,4)=-1.92144E-08

	v_names(11)="NH3"
	v_atoms(11,1)=3
	v_atoms(11,4)=1
	bv(11,0)=1.39343E+05
	bv(11,1)=-6.39532E+00
	bv(11,2)=4.95981E+00
	bv(11,3)=1.81530E-03
	bv(11,4)=-8.85794E-08

	v_names(12)="Zn"
	v_atoms(12,30)=1

	v_names(13)="Mn"
	v_atoms(13,27)=1

	v_names(14)="Cr"
	v_atoms(14,26)=1

	v_names(15)="CH4"
	v_atoms(15,1)=4
	v_atoms(15,3)=1
	bv(15,0)=1.97846E+05
	bv(15,1)=-8.83168E+00
	bv(15,2)=5.27931E+00
	bv(15,3)=2.75677E-03
	bv(15,4)=-1.39667E-07

	v_names(16)="SO2"
	v_atoms(16,11)=1
	v_atoms(16,5)=2
	bv(16,0)=1.27905E+05
	bv(16,1)=-3.68505E+00
	bv(16,2)=-5.04460E+00
	bv(16,3)=1.38250E-03
	bv(16,4)=-8.66752E-08

	v_include=.false.

	nnr=(nr-1)*nr_cloud+1
	nCS=Cloud(ii)%nmat
	if(Cloud(ii)%haze) nCS=nCS-1
	allocate(Kd(nnr),Kg(nnr),Km(nnr))
	allocate(logCloudP(nnr))
	allocate(CloudMMW(nnr),CloudHp(nnr),Sat(nnr,nCS),Sat0(nnr,nCS),fSat(nnr,nCS),CloudG(nnr))
	
	niter=1000
	nconv=20
	if(computeT) then
		if(nTiter.eq.1) then
			niter=50
			nconv=5
		else if(nTiter.le.3) then
			niter=100
			nconv=10
		endif
	endif
	if(Cloud(ii)%computeJn) niter=niter*4
	
	allocate(docondense(nCS))

	if(.not.allocated(CloudP)) then
		allocate(CloudP(nnr))
		allocate(CloudT(nnr))
		allocate(CloudR(nnr))
		allocate(Clouddens(nnr))
		allocate(xv(nVS,nnr))
		allocate(xc(nCS,nnr))
		allocate(xn(nnr))
		allocate(xnv(nnr))
		allocate(rpart(nnr))
	endif
	if(.not.allocated(ATP)) then
		allocate(ATP(nCS),BTP(nCS))
		allocate(rhodust(nCS))
		allocate(atoms_cloud(nCS,N_atoms))
		allocate(xv_bot(nVS))
		allocate(mu(nCS))
		allocate(CSname(nCS),maxT(nCS))
	endif
	allocate(muC(nCS))
	allocate(muV(nVS))
	allocate(v_cloud(nCS,nVS),iVL(nnr,nCS),v_H2(nCS))
	allocate(icryst(nCS),iamorph(nCS),tcrystinv(nnr))
	allocate(A_J(nCS),B_J(nCS),do_nuc(nCS),Jn_xv(nnr,nCS),sigma_nuc(nCS),r0_nuc(nCS),Nf_nuc(nCS),Nc_nuc(nnr,nCS))
	allocate(inuc(nCS))
	allocate(bc(nCS,0:4),ifit(nCS))
	bc=0d0
	ifit=0

	call SetAbun

	COabun=min(molfracs_atoms(3),molfracs_atoms(5))
	molfracs_atoms(3)=molfracs_atoms(3)-COabun
	molfracs_atoms(5)=molfracs_atoms(5)-COabun

	atoms_cloud=0
	v_cloud=0d0
	v_H2=0d0
	icryst=0
	iamorph=0
	do_nuc=.false.
	maxT=10000d0
	do i=1,nCS
		select case(Cloud(ii)%condensate(i))
			case('SiO2','QUARTZ')
				CSname(i)='SiO2'
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=2
				v_cloud(i,1)=1
				v_cloud(i,4)=1
				v_H2(i)=-1
				v_include(1)=.true.
				v_include(4)=.true.
				rhodust(i)=2.65
				bc(i,0)=2.222506e+05
				bc(i,1)=-5.478967e+00
				bc(i,2)=-1.969170e+01
				bc(i,3)=4.907829e-03
				bc(i,4)=-5.612476e-07
				ifit(i)=0
			case('MgSiO3','ENSTATITE')
				CSname(i)=trim(Cloud(ii)%condensate(i))
				atoms_cloud(i,7)=1
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=3
				v_cloud(i,1)=1
				v_cloud(i,3)=1
				v_cloud(i,4)=2
				v_H2(i)=-2
				v_include(1)=.true.
				v_include(3)=.true.
				v_include(4)=.true.
				rhodust(i)=3.19
				bc(i,0)=3.458775e+05
				bc(i,1)=-7.257035e+00
				bc(i,2)=-4.345839e+01
				bc(i,3)=7.516612e-03
				bc(i,4)=-8.200336e-07
				ifit(i)=0
			case('Mg2SiO4','FORSTERITE')
				CSname(i)=trim(Cloud(ii)%condensate(i))
				atoms_cloud(i,7)=2
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=4
				v_cloud(i,1)=1
				v_cloud(i,3)=2
				v_cloud(i,4)=3
				v_H2(i)=-3
				v_include(1)=.true.
				v_include(3)=.true.
				v_include(4)=.true.
				rhodust(i)=3.21
				bc(i,0)=4.685652e+05
				bc(i,1)=-9.176566e+00
				bc(i,2)=-6.544643e+01
				bc(i,3)=9.627992e-03
				bc(i,4)=-1.038920e-06
				ifit(i)=0
			case('FeSiO3','FERROSILITE')
				CSname(i)=trim(Cloud(ii)%condensate(i))
				atoms_cloud(i,17)=1
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=3
				v_cloud(i,1)=1
				v_cloud(i,6)=1
				v_cloud(i,4)=2
				v_H2(i)=-2
				v_include(1)=.true.
				v_include(6)=.true.
				v_include(4)=.true.
				rhodust(i)=3.5
				bc(i,0)=3.361566e+05
				bc(i,1)=-6.749696e+00
				bc(i,2)=-4.728394e+01
				bc(i,3)=7.277447e-03
				bc(i,4)=-7.573112e-07
				ifit(i)=0
			case('Fe2SiO4','FAYALITE')
				CSname(i)=trim(Cloud(ii)%condensate(i))
				atoms_cloud(i,17)=2
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=4
				v_cloud(i,1)=1
				v_cloud(i,6)=2
				v_cloud(i,4)=3
				v_H2(i)=-3
				v_include(1)=.true.
				v_include(6)=.true.
				v_include(4)=.true.
				rhodust(i)=4.39
				bc(i,0)=4.495737e+05
				bc(i,1)=-9.271411e+00
				bc(i,2)=-6.615000e+01
				bc(i,3)=1.041250e-02
				bc(i,4)=-1.092031e-06
				ifit(i)=0
			case('NaAlSi3O8')
				CSname(i)='NaAlSi3O8'
				atoms_cloud(i,6)=1
				atoms_cloud(i,8)=1
				atoms_cloud(i,9)=3
				atoms_cloud(i,5)=8
				v_cloud(i,1)=3
				v_cloud(i,4)=5
				v_cloud(i,7)=1
				v_cloud(i,8)=1
				v_H2(i)=-5
				v_include(1)=.true.
				v_include(4)=.true.
				v_include(7)=.true.
				v_include(8)=.true.
				rhodust(i)=2.36
				bc(i,0)=9.232283e+05
				bc(i,1)=-2.012507e+01
				bc(i,2)=-1.039484e+02
				bc(i,3)=1.953868e-02
				bc(i,4)=-2.217366e-06
				ifit(i)=0
			case('MgO')
				CSname(i)='MgO'
				atoms_cloud(i,7)=1
				atoms_cloud(i,5)=1
				v_cloud(i,3)=1
				v_cloud(i,4)=1
				v_H2(i)=-1
				v_include(3)=.true.
				v_include(4)=.true.
				rhodust(i)=3.58
				bc(i,0)=1.195366e+05
				bc(i,1)=-2.104910e+00
				bc(i,2)=-2.122681e+01
				bc(i,3)=2.564436e-03
				bc(i,4)=-2.957752e-07
				ifit(i)=0
			case('H2O','WATER')
				CSname(i)='H2O'
				atoms_cloud(i,1)=2
				atoms_cloud(i,5)=1
				v_cloud(i,4)=1
				v_include(4)=.true.
				rhodust(i)=0.93
				maxT(i)=747d0
				do_nuc(i)=Cloud(ii)%ComputeJn
				inuc(i)=4
				sigma_nuc(i)=109.0
				r0_nuc(i)=1.973e-8
				Nf_nuc(i)=1d0
				ifit(i)=-1
			case('Fe','IRON')
				CSname(i)='Fe'
				atoms_cloud(i,17)=1
				v_cloud(i,6)=1
				v_include(6)=.true.
				rhodust(i)=7.87
				do_nuc(i)=.false.!Cloud(ii)%ComputeJn
				inuc(i)=6
				sigma_nuc(i)=1870	! from Brooks et al. 2001
				r0_nuc(i)=3.7e-9
				Nf_nuc(i)=1d0
				bc(i,0)=4.989867e+04
				bc(i,1)=-6.889399e-01
				bc(i,2)=-1.461748e+01
				bc(i,3)=1.161476e-03
				bc(i,4)=-5.259868e-08
				ifit(i)=0
			case('FeS','TROILITE')
				CSname(i)='FeS'
				atoms_cloud(i,17)=1
				atoms_cloud(i,11)=1
				v_cloud(i,5)=1
				v_cloud(i,6)=1
				v_H2(i)=-1
				v_include(5)=.true.
				v_include(6)=.true.
				rhodust(i)=4.83
				bc(i,0)=9.450871e+04
				bc(i,1)=-2.495865e+00
				bc(i,2)=-1.829773e+01
				bc(i,3)=3.174381e-03
				bc(i,4)=-3.069607e-07
				ifit(i)=0
			case('FeO')
				CSname(i)='FeO'
				atoms_cloud(i,17)=1
				atoms_cloud(i,5)=1
				v_cloud(i,6)=1
				v_cloud(i,4)=1
				v_H2(i)=-1
				v_include(6)=.true.
				v_include(4)=.true.
				rhodust(i)=5.9
				bc(i,0)=1.121860e+05
				bc(i,1)=-2.152325e+00
				bc(i,2)=-2.083896e+01
				bc(i,3)=2.956682e-03
				bc(i,4)=-3.223296e-07
				ifit(i)=0
			case('Fe2O3')
				CSname(i)='Fe2O3'
				atoms_cloud(i,17)=2
				atoms_cloud(i,5)=3
				v_cloud(i,6)=2
				v_cloud(i,4)=3
				v_H2(i)=-3
				v_include(6)=.true.
				v_include(4)=.true.
				rhodust(i)=5.24
				bc(i,0)=2.876529e+05
				bc(i,1)=-6.810027e+00
				bc(i,2)=-4.976293e+01
				bc(i,3)=8.724692e-03
				bc(i,4)=-1.038569e-06
				ifit(i)=0
			case('Al2O3','CORRUNDUM')
				CSname(i)='Al2O3'
				atoms_cloud(i,5)=3
				atoms_cloud(i,8)=2
				v_cloud(i,4)=3
				v_cloud(i,7)=2
				v_H2(i)=-3
				v_include(4)=.true.
				v_include(7)=.true.
				rhodust(i)=3.97
				bc(i,0)=3.686315e+05
				bc(i,1)=-8.641215e+00
				bc(i,2)=-3.798714e+01
				bc(i,3)=8.874690e-03
				bc(i,4)=-1.024209e-06
				ifit(i)=0
			case("NaCl")
				CSname(i)='NaCl'
				atoms_cloud(i,6)=1
				atoms_cloud(i,12)=1
				v_cloud(i,8)=1
				v_cloud(i,10)=1
				v_H2(i)=-0.5
				v_include(8)=.true.
				v_include(10)=.true.
				rhodust(i)=2.17
				do_nuc(i)=Cloud(ii)%ComputeJn
				inuc(i)=10
				sigma_nuc(i)=113.3
				r0_nuc(i)=2.205e-8
				Nf_nuc(i)=1d0
				bc(i,0)=7.717022e+04
				bc(i,1)=3.607723e-01
				bc(i,2)=-3.269074e+01
				bc(i,3)=9.594929e-04
				bc(i,4)=1.834797e-08
				ifit(i)=0
			case("KCl")
				CSname(i)='KCl'
				atoms_cloud(i,13)=1
				atoms_cloud(i,12)=1
				v_cloud(i,9)=1
				v_cloud(i,10)=1
				v_H2(i)=-0.5
				v_include(9)=.true.
				v_include(10)=.true.
				rhodust(i)=1.99
				do_nuc(i)=Cloud(ii)%ComputeJn
				inuc(i)=10
				sigma_nuc(i)=100.3
				r0_nuc(i)=2.462e-8
				Nf_nuc(i)=1d0
				bc(i,0)=7.801979e+04
				bc(i,1)=4.147547e-01
				bc(i,2)=-3.263640e+01
				bc(i,3)=1.037051e-03
				bc(i,4)=6.931110e-09
				ifit(i)=0
			case("Na2S")
				CSname(i)='Na2S'
				atoms_cloud(i,6)=1
				atoms_cloud(i,11)=2
				v_cloud(i,8)=2
				v_cloud(i,5)=1
				v_H2(i)=-1
				v_include(8)=.true.
				v_include(5)=.true.
				rhodust(i)=1.86
				bc(i,0)=1.99053E+06
				bc(i,1)=-8.70027E+05
				bc(i,2)=4.06721E+02
				bc(i,3)=-2.78298E-02
				bc(i,4)=0.00000E+00
				ifit(i)=2
			case("NH3","AMONIA")
				CSname(i)='NH3'
				atoms_cloud(i,1)=3
				atoms_cloud(i,4)=1
				v_cloud(i,11)=1
				v_include(11)=.true.
				rhodust(i)=0.87
				maxT(i)=220d0
				do_nuc(i)=Cloud(ii)%ComputeJn
				inuc(i)=11
				sigma_nuc(i)=23.4
				r0_nuc(i)=1.980e-8
				Nf_nuc(i)=1d0
				ifit(i)=-1
			case('TiO2')
				CSname(i)='TiO2'
				atoms_cloud(i,5)=1
				atoms_cloud(i,15)=2
				v_cloud(i,2)=1
				v_cloud(i,4)=1
				v_H2(i)=-1
				v_include(2)=.true.
				v_include(4)=.true.
				rhodust(i)=4.23
				do_nuc(i)=Cloud(ii)%ComputeJn
				A_J(i)=1.112e12
				B_J(i)=-24.0
				inuc(i)=2
				sigma_nuc(i)=480.6
				r0_nuc(i)=1.956e-8
				Nf_nuc(i)=0d0
				bc(i,0)=2.296961e+05
				bc(i,1)=-3.392573e+00
				bc(i,2)=-3.338416e+01
				bc(i,3)=3.426975e-03
				bc(i,4)=-3.572185e-07
				ifit(i)=0
			case('H2SO4')
				CSname(i)='H2SO4'
				atoms_cloud(i,1)=2
				atoms_cloud(i,5)=4
				atoms_cloud(i,11)=1

				v_cloud(i,4)=2
				v_cloud(i,16)=1
				v_H2(i)=-1
				v_include(4)=.true.
				v_include(16)=.true.

c				v_cloud(i,4)=4
c				v_cloud(i,5)=1
c				v_H2(i)=-4
c				v_include(4)=.true.
c				v_include(5)=.true.

				rhodust(i)=1.84d0
				bc(i,0)=9.70368E+05
				bc(i,1)=-2.53825E+06
				bc(i,2)=9.35422E+02
				bc(i,3)=-4.96224E-02
				bc(i,4)=0.00000E+00
				ifit(i)=2
			case('MnS')
				CSname(i)='MnS'
				atoms_cloud(i,27)=1
				atoms_cloud(i,11)=1
				v_cloud(i,13)=1
				v_cloud(i,5)=1
				v_H2(i)=-1
				v_include(13)=.true.
				v_include(5)=.true.
				rhodust(i)=4.08
				bc(i,0)=1.12482E+05
				bc(i,1)=-1.81938E+05
				bc(i,2)=5.87107E+01
				bc(i,3)=8.89360E-05
				bc(i,4)=-4.20876E-09
				ifit(i)=1
			case('Cr')
				CSname(i)='Cr'
				atoms_cloud(i,26)=1
				v_cloud(i,14)=1
				v_include(14)=.true.
				rhodust(i)=7.19d0
				do_nuc(i)=Cloud(ii)%ComputeJn
				inuc(i)=14
				sigma_nuc(i)=3330.0
				r0_nuc(i)=1.421
				Nf_nuc(i)=1d0
				ifit(i)=-1
			case('SiO')
				CSname(i)='SiO'
				atoms_cloud(i,9)=1
				atoms_cloud(i,5)=1
				v_cloud(i,1)=1
				v_include(1)=.true.
				rhodust(i)=2.18
				maxT(i)=5000d0
				do_nuc(i)=Cloud(ii)%computeJn
				A_J(i)=4.4e12
				B_J(i)=1.33
				inuc(i)=1
				sigma_nuc(i)=849.4
				r0_nuc(i)=2.0e-8
				Nf_nuc(i)=1d0
				ifit(i)=-1
			case('ZnS')
				CSname(i)='ZnS'
				atoms_cloud(i,30)=1
				atoms_cloud(i,11)=1
				v_cloud(i,12)=1
				v_cloud(i,5)=1
				v_H2(i)=-1
				v_include(12)=.true.
				v_include(5)=.true.
				rhodust(i)=4.09
				bc(i,0)=7.30636867e+04
				bc(i,1)=-2.72297464e+00
				bc(i,2)=-1.60451251e+01
				bc(i,3)=3.44107646e-03
				bc(i,4)=-2.56383825e-07
				ifit(i)=0
			case('Zn')
				CSname(i)='Zn'
				atoms_cloud(i,30)=1
				v_cloud(i,12)=1
				v_include(12)=.true.
				rhodust(i)=7.14d0
				bc(i,0)=1.56414E+04     
				bc(i,1)=-1.32671E+00
				bc(i,2)=-7.87964E+00
				bc(i,3)=4.59419E-03
				bc(i,4)=-1.46651E-06
				ifit(i)=0
			case default
				call output("Unknown condensate: " // trim(Cloud(ii)%condensate(i)))
				print*,'unknown condensate: ', trim(Cloud(ii)%condensate(i))
				stop
		end select
	enddo
	
	do iCS=1,nCS
		do jCS=1,nCS
			if(Cloud(ii)%condensate(iCS).eq.'FORSTERITE'.and.Cloud(ii)%condensate(jCS).eq.'Mg2SiO4') then
				iamorph(iCS)=jCS
				icryst(jCS)=iCS
			endif
			if(Cloud(ii)%condensate(iCS).eq.'ENSTATITE'.and.Cloud(ii)%condensate(jCS).eq.'MgSiO3') then
				iamorph(iCS)=jCS
				icryst(jCS)=iCS
			endif
			if(Cloud(ii)%condensate(iCS).eq.'FAYALITE'.and.Cloud(ii)%condensate(jCS).eq.'Fe2SiO4') then
				iamorph(iCS)=jCS
				icryst(jCS)=iCS
			endif
			if(Cloud(ii)%condensate(iCS).eq.'FERROSILITE'.and.Cloud(ii)%condensate(jCS).eq.'FeSiO3') then
				iamorph(iCS)=jCS
				icryst(jCS)=iCS
			endif
		enddo
	enddo

	docondense=.true.

	haze_nuc=.false.
	select case(Cloud(ii)%hazetype)
		case("SOOT","soot","Soot")
			rho_nuc=1.00
		case("THOLIN","tholin","Tholin")
			rho_nuc=1.00
		case("optEC")
			rho_nuc=1.50
			ihaze=15
			v_include(15)=(.not.dochemistry)
			haze_nuc=(.not.dochemistry)
		case("SiC")
			rho_nuc=3.22
		case("CARBON","Carbon","carbon")
			rho_nuc=1.80
		case("CORRUNDUM","Corrundum","corrundum","Al2O3")
			rho_nuc=3.97
		case("IRON","Iron","iron","Fe")
			rho_nuc=7.87
		case("SiO")
			rho_nuc=2.18
		case("TiO2")
			rho_nuc=4.23
		case("Enstatite","enstatite","ENSTATITE")
			rho_nuc=3.20
		case default
			call output("hazetype unknown")
			stop
	end select

	do i=1,nCS
		muC(i)=sum(mass_atoms(1:N_atoms)*atoms_cloud(i,1:N_atoms))
	enddo
	do i=1,nVS
		muV(i)=sum(mass_atoms(1:N_atoms)*v_atoms(i,1:N_atoms))
	enddo

	if(dochemistry) then
		mutot=COabun*(mass_atoms(3)+mass_atoms(5))
		xv_bot=1d200
		do i=1,N_atoms
			mutot=mutot+mass_atoms(i)*molfracs_atoms(i)
		enddo
	else
		mutot=0d0
		do i=1,nmol
			if(includemol(i)) then
				mutot=mutot+Mmol(i)*mixrat_r(1,i)
			endif
		enddo
	endif
	do iVS=1,nVS
		if(v_include(iVS)) then
			if(dochemistry) then
				do k=1,N_atoms
					if(v_atoms(iVS,k).gt.0d0) then
						f=molfracs_atoms(k)/v_atoms(iVS,k)
						if(f.lt.xv_bot(iVS)) then
							xv_bot(iVS)=f
						endif
					endif
				enddo
				do k=1,N_atoms
					molfracs_atoms(k)=molfracs_atoms(k)-xv_bot(iVS)*v_atoms(iVS,k)
					if(molfracs_atoms(k).lt.0d0) molfracs_atoms(k)=0d0
				enddo
			else
				do k=1,nmol
					if(molname(k).eq.v_names(iVS)) then
						xv_bot(iVS)=mixrat_r(1,k)
					endif
				enddo
			endif
		endif
	enddo
	
	molfracs_atoms0=molfracs_atoms
	xv_bot=xv_bot*muV/mutot
	do iVS=1,nVS
		if(.not.v_include(iVS)) xv_bot(iVS)=0d0
	enddo

c	print*,xv_bot(1:7)
c	xv_bot(1:7) = 10.0**metallicity*(/ 6.1e-4, 2.4e-6, 4.1e-4, 1.4e-3, 1.9e-4, 7.6e-4, 3.2e-5 /)
c	print*,xv_bot(1:7)

	do iCS=1,nCS
		if(do_nuc(iCS).and.CSname(iCS).ne.'SiO'.and.CSname(iCS).eq.'TiO2') then
			A_J(iCS)=(16.*pi*sigma_nuc(iCS)**3*(muC(iCS)*mp/rhodust(iCS))**2/(3.*kb**3))
			B_J(iCS)=log(sqrt(2.*sigma_nuc(iCS)/(pi*muC(iCS)*mp))*(muC(iCS)*mp/rhodust(iCS)))
		endif
	enddo

	Cloud(ii)%frac=0d0

	sigmastar=0.1
	Pstar=60d-6

	Pstar=Cloud(ii)%P
	sigmastar=log(Cloud(ii)%dP)

	fstick=1d0
	
	sigmamol=8d-15

	Sigmadot=Cloud(ii)%Sigmadot
	

	m_nuc=4d0*pi*Cloud(ii)%rnuc**3*rho_nuc/3d0

	allocate(mpart(nnr))
	allocate(rho_av(nnr))
	allocate(y(nnr,5))
	allocate(Sn(nnr))
	allocate(vth(nnr))
	allocate(tcinv(niter,nnr))
	allocate(xn_iter(niter,nnr))
	allocate(xc_iter(niter,nCS,nnr),xv_iter(niter,nVS,nnr))
	allocate(vsed(nnr))

	allocate(ixv(nVS,nnr))
	allocate(ixc(nCS,nnr))
	allocate(ixn(nnr))

	j=0
	do i=1,nnr
		if(.not.Cloud(ii)%usefsed) then
			j=j+1
			ixn(i)=j
		endif
		do iCS=1,nCS
			j=j+1
			ixc(iCS,i)=j
		enddo
		do iVS=1,nVS
			if(v_include(iVS)) then
				j=j+1
				ixv(iVS,i)=j
			endif
		enddo
		if(i.eq.1) NStot=j
	enddo
	NN=j

	rho_av=sum(rhodust)/real(nCS)

	do i=1,nnr
		CloudP(i)=10d0**(log10(P(1))+log10(P(nr)/P(1))*real(i-1)/real(nnr-1))
	enddo

	logP(1:nr)=-log(P(1:nr))
	logCloudP(1:nnr)=-log(CloudP(1:nnr))

	SKIP=.false.
	INCFD=1
	logx(1:nr)=log(T(1:nr))
	call DPCHIM(nr,logP,logx,dlogx,INCFD)
	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudT,IERR)
	CloudT(1:nnr)=exp(CloudT(1:nnr))

	SKIP=.false.
	INCFD=1
	logx(1:nr)=log(dens(1:nr))
	call DPCHIM(nr,logP,logx,dlogx,INCFD)
	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,Clouddens,IERR)
	Clouddens(1:nnr)=exp(Clouddens(1:nnr))

	SKIP=.false.
	INCFD=1
	logx(1:nr)=log(R(1:nr))
	call DPCHIM(nr,logP,logx,dlogx,INCFD)
	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudR,IERR)
	CloudR(1:nnr)=exp(CloudR(1:nnr))

	SKIP=.false.
	INCFD=1
	logx(1:nr)=MMW(1:nr)
	call DPCHIM(nr,logP,logx,dlogx,INCFD)
	call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudMMW,IERR)

	if(complexKzz.or.Cloud(ii)%usefsed) then
		SKIP=.false.
		INCFD=1
		logx(1:nr)=log(Hp(1:nr))
		call DPCHIM(nr,logP,logx,dlogx,INCFD)
		call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudHp,IERR)
		CloudHp(1:nnr)=exp(CloudHp(1:nnr))
	endif
	if(Cloud(ii)%globalKzz.or.Cloud(ii)%Kzz.le.0d0) then
		SKIP=.false.
		INCFD=1
		logx(1:nr)=log10(Kzz_b(1:nr))
		call DPCHIM(nr,logP,logx,dlogx,INCFD)
		call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,Km,IERR)
		Km=10d0**Km
		do i=1,nnr
			Kd(i)=Km(i)
			if(complexKzz) then
				SKIP=.false.
				INCFD=1
				logx(1:nr)=log10(Kzz_g(1:nr))
				call DPCHIM(nr,logP,logx,dlogx,INCFD)
				call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,Kg,IERR)
				Kg=10d0**Kg
			else
				Kg(i)=Kd(i)
			endif
		enddo
	else
		Km=Cloud(ii)%Kzz
		Kd=Km
		Kg=Kd
	endif

	allocate(x(NN),x0(NN))

	rpart=Cloud(ii)%rnuc
	xn=0d0
	xv=0d0
	xc=0d0
	do j=1,nVS
		do i=1,nnr
			xv(j,i)=xv_bot(j)
		enddo
	enddo
	tcinv=0d0
	Jn_xv=0d0
	Nc_nuc=0d0
	
	docondense=.true.
	if(Cloud(ii)%hazetype.eq.'optEC') then
		allocate(CloudtauUV(nnr),CloudkappaUV(nnr))
		do ir=1,nr
			if(kappaUV0.gt.0d0) then
				tauUV(ir)=exp(-kappaUV0*1d6*P(ir)/grav(ir))
			else if(tauUV(ir).lt.0d0) then
				tauUV(ir)=exp(-1d6*P(ir)/grav(ir))
			endif
		enddo
		SKIP=.false.
		INCFD=1
		logx(1:nr)=tauUV(1:nr)
		call DPCHIM(nr,logP,logx,dlogx,INCFD)
		call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudtauUV,IERR)
		logx(1:nr)=kappaUV(1:nr)
		call DPCHIM(nr,logP,logx,dlogx,INCFD)
		call DPCHFE(nr,logP,logx,dlogx,INCFD,SKIP,nnr,logCloudP,CloudkappaUV,IERR)
		tot=0d0
	endif

	if(.not.retrieval) then
		open(unit=20,file='Tevap.dat',RECL=6000)
		form='("#",a18,' // trim(int2string(nCS,'(i4)')) // 'a23)'
		write(20,form) "T[K]",(trim(CSname(i)),i=1,nCS)
		form='(es19.7E3,' // trim(int2string(nCS,'(i4)')) // 'es23.7E3)'
		do i=10,5000,10
			tot1=real(i)
			do iCS=1,nCS
				select case(CSname(iCS))
					case("H2O")
						call PvapH2O(tot1,Sat(1,iCS),liq)
						Sat(1,iCS)=1d0/Sat(1,iCS)
					case("NH3")
c	Gibbs energy as derived from Eq from GGChem paper does not work at high pressures. Use adjusted Pvap.
						call PvapNH3(tot1,Sat(1,iCS),liq)
						Sat(1,iCS)=1d0/Sat(1,iCS)
					case("SiO")
c	Gibbs energy as derived from Eq from GGChem paper does not work at high pressures. Use adjusted Pvap.
						call PvapSiO(tot1,Sat(1,iCS),liq)
						Sat(1,iCS)=1d0/Sat(1,iCS)
					case("Cr")
c	Gibbs energy as derived from Eq from GGChem paper does not work at high pressures. Use adjusted Pvap.
						call PvapCr(tot1,Sat(1,iCS),liq)
						Sat(1,iCS)=1d0/Sat(1,iCS)
					case default
						Gibbs=compGibbs(tot1,bc(iCS,0:4),ifit(iCS))
						Sat(1,iCS)=Gibbs-v_H2(iCS)*compGibbs(tot1,bH2(0:4),0)
						do iVS=1,nVS
							if(v_cloud(iCS,iVS).gt.0d0.and.v_include(iVS)) then
								Gibbs=compGibbs(tot1,bv(iVS,0:4),0)
								Sat(1,iCS)=Sat(1,iCS)-Gibbs*v_cloud(iCS,iVS)
							endif
						enddo
						Sat(1,iCS)=exp(-Sat(1,iCS))
				end select
				if(tot1.gt.maxT(iCS)) Sat(1,iCS)=0d0
				if(Sat(1,iCS).lt.1d-40) Sat(1,iCS)=1d-40
				do iVS=1,nVS
					if(v_include(iVS)) then
						Sat(1,iCS)=Sat(1,iCS)*(xv_bot(iVS)*CloudMMW(1)/muV(iVS))**v_cloud(iCS,iVS)
					endif
				enddo
				Sat(1,iCS)=(1d0/Sat(1,iCS))**(1d0/(v_H2(iCS)+sum(v_cloud(iCS,1:nVS))))			
			enddo
			write(20,form) tot1,Sat(1,1:nCS)
		enddo
		close(unit=20)
	endif


c Compute crystallinity
c values from Fabian et al. 2000
	nucryst=2d13
	Tcryst=41000d0

	do i=1,nnr
		do iCS=1,nCS
			select case(CSname(iCS))
				case("H2O")
					call PvapH2O(CloudT(i),Sat(i,iCS),liq)
					Sat(i,iCS)=CloudP(i)/Sat(i,iCS)
				case("NH3")
c	Gibbs energy as derived from Eq from GGChem paper does not work at high pressures. Use adjusted Pvap.
					call PvapNH3(CloudT(i),Sat(i,iCS),liq)
					Sat(i,iCS)=CloudP(i)/Sat(i,iCS)
				case("SiO")
c	Gibbs energy as derived from Eq from GGChem paper does not work at high pressures. Use adjusted Pvap.
					call PvapSiO(CloudT(i),Sat(i,iCS),liq)
					Sat(i,iCS)=CloudP(i)/Sat(i,iCS)
				case("Cr")
c	Gibbs energy as derived from Eq from GGChem paper does not work at high pressures. Use adjusted Pvap.
					call PvapCr(CloudT(i),Sat(i,iCS),liq)
					Sat(i,iCS)=CloudP(i)/Sat(i,iCS)
				case default
					Gibbs=compGibbs(CloudT(i),bc(iCS,0:4),ifit(iCS))
					Sat(i,iCS)=Gibbs-v_H2(iCS)*compGibbs(CloudT(i),bH2(0:4),0)
					do iVS=1,nVS
						if(v_cloud(iCS,iVS).gt.0d0.and.v_include(iVS)) then
							Gibbs=compGibbs(CloudT(i),bv(iVS,0:4),0)
							Sat(i,iCS)=Sat(i,iCS)-Gibbs*v_cloud(iCS,iVS)
						endif
					enddo
					Sat(i,iCS)=exp(-Sat(i,iCS))*CloudP(i)
			end select
			if(CloudT(i).gt.maxT(iCS)) Sat(i,iCS)=0d0
			if(Sat(i,iCS).lt.1d-40) Sat(i,iCS)=1d-40
		enddo
		CloudG(i)=Ggrav*Mplanet/CloudR(i)**2
		if(Cloud(ii)%hazetype.eq.'optEC') then
			Sn(i)=Clouddens(i)*CloudP(i)*CloudtauUV(i)
			if(i.eq.nnr) then
				tot=tot+abs(CloudR(i-1)-CloudR(i))*Sn(i)
			else if(i.eq.1) then
				tot=tot+abs(CloudR(i+1)-CloudR(i))*Sn(i)
			else
				tot=tot+abs(CloudR(i-1)-CloudR(i+1))*0.5*Sn(i)
			endif
		else
			Sn(i)=(Clouddens(i)*CloudG(i)*Sigmadot/(sigmastar*CloudP(i)*1d6*sqrt(2d0*pi)))*exp(-log(CloudP(i)/Pstar)**2/(2d0*sigmastar**2))
		endif
		tcrystinv(i)=nucryst*exp(-Tcryst/CloudT(i))
	enddo
	if(Cloud(ii)%hazetype.eq.'optEC') Sn=Sn*scaleUV*Sigmadot/tot

	Sat0=Sat


	if(Cloud(ii)%rainout) then
		allocate(c_rainout(nCS))
		c_rainout=.false.
		f=1d-4
		i=1
		do iter=1,10000
			xv(i,1:nVS)=0d0
			j=0
			do iCS=1,nCS
				tot1=1d200
				iVL(i,iCS)=1
				do iVS=1,nVS
					if(v_cloud(iCS,iVS).gt.0d0.and.v_include(iVS)) then
						tot2=xv_bot(iVS)/v_cloud(iCS,iVS)
						if(tot2.lt.tot1) then
							tot1=tot2
							iVL(i,iCS)=iVS
						endif
					endif
				enddo
				fSat(i,iCS)=CloudP(i)**v_H2(iCS)
				do iVS=1,nVS
					if(v_include(iVS).and.v_cloud(iCS,iVS).ne.0d0) then
						fSat(i,iCS)=fSat(i,iCS)*(CloudP(i)*xv_bot(iVS)*CloudMMW(i)/muV(iVS))**v_cloud(iCS,iVS)
					endif
				enddo
				tot1=Sat(i,iCS)*fSat(i,iCS)/Cloud(ii)%Srainout
				if(tot1.gt.1d0) then
					j=j+1
					tot1=xv_bot(iVL(i,iCS))*(1d0-1d0/tot1)
					c_rainout(iCS)=.true.
					do iVS=1,nVS
						if(v_cloud(iCS,iVS).gt.0d0.and.v_include(iVS)) then
							xv(i,iVS)=xv(i,iVS)+tot1*v_cloud(iCS,iVS)*muV(iVS)/(v_cloud(iCS,iVL(i,iCS))*muV(iVL(i,iCS)))
						endif
					enddo
				endif
			enddo
			if(j.eq.0) exit
			maxerr=1d200
			do iVS=1,nVS
				if(v_include(iVS).and.xv(i,iVS).gt.0d0) then
					err=xv_bot(iVS)/xv(i,iVS)
					if(err.lt.maxerr) maxerr=err
				endif
			enddo
			if(.not.maxerr.lt.1d100) exit
			xv_bot(1:nVS)=xv_bot(1:nVS)-maxerr*0.1*xv(i,1:nVS)
		enddo
c		do iCS=1,nCS
c			if(c_rainout(iCS)) then
c				print*,'Partial rainout of ',CSname(iCS)
c			endif
c		enddo
		deallocate(c_rainout)
	endif

	if(Cloud(ii)%computeJn) then
		fscale=1d-8
		pscale=(1d0/fscale)**(1d0/min(real(niter/2),50.0))
		f=0.98
		nfscale=0
		eps=1d-2
	else
		fscale=1d0
		pscale=1d0
		f=0.75
		nfscale=1000
		eps=1d-3

		fscale=1d-8
		pscale=(1d0/fscale)**(1d0/min(real(niter/2),50.0))
		nfscale=0
	endif

	allocate(vthv(nnr))
	allocate(Sc(nnr,nCS))
	allocate(IWORK(NN))
	allocate(AB(2*NStot+NStot+1,NN))

	iconv=0
	Jn_xv=0d0
	
	
	
	
c start the loop
	do iter=1,niter
	call tellertje(iter,niter)

	vsed=0d0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,tot1,tot2,cs,fsed,tot,iCS,iVS,lmfp,St,Dp,Jn_temp)
!$OMP& SHARED(nnr,nCS,nVS,iVL,v_cloud,CloudT,xv,muV,fSat,Sat0,Sat,v_include,CloudP,CloudMMW,v_H2,vth,ii,CloudR,Rplanet,sigmamol,
!$OMP&		vsed,Kd,CloudHp,rpart,rho_av,CloudG,Clouddens,Cloud,mpart,xv_bot,xn,xc,Km,vthv,Sc,muC,do_nuc,CSname,Nc_nuc,Jn_xv,m_nuc,
!$OMP&		iter,complexKzz,fstick,inuc,A_J,B_J,Nf_nuc,r0_nuc,sigma_nuc,fscale,nfscale,f)
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC, 1)
	do i=1,nnr
		do iCS=1,nCS
			tot1=1d200
			iVL(i,iCS)=1
			do iVS=1,nVS
				if(v_cloud(iCS,iVS).gt.0d0.and.v_include(iVS)) then
					tot2=xv(iVS,i)*(sqrt(8d0*kb*CloudT(i)/(pi*muV(iVS)*mp)))/(muV(iVS)*v_cloud(iCS,iVS))
					if(tot2.lt.tot1) then
						tot1=tot2
						iVL(i,iCS)=iVS
					endif
				endif
			enddo
			fSat(i,iCS)=CloudP(i)**v_H2(iCS)*(CloudMMW(i)/muV(iVL(i,iCS)))*
     &				(CloudP(i)*xv(iVL(i,iCS),i)*CloudMMW(i)/muV(iVL(i,iCS)))**(v_cloud(iCS,iVL(i,iCS))-1d0)
			do iVS=1,nVS
				if(v_include(iVS).and.iVS.ne.iVL(i,iCS)) then
					fSat(i,iCS)=fSat(i,iCS)*(CloudP(i)*xv(iVS,i)*CloudMMW(i)/muV(iVS))**v_cloud(iCS,iVS)
				endif
			enddo
			if(fSat(i,iCS).lt.1d-40) fSat(i,iCS)=1d-40
		enddo

		cs=sqrt(kb*CloudT(i)/(CloudMMW(i)*mp))
		vth(i)=sqrt(8d0*kb*CloudT(i)/(pi*CloudMMW(i)*mp))
		if(complexKzz) then
			St=rpart(i)*rho_av(i)*Km(i)/(vth(i)*Clouddens(i)*CloudHp(i)**2)
			Kd(i)=Km(i)/(1d0+St)
		endif
		if(Cloud(ii)%usefsed) then
			fsed=Cloud(ii)%fsed_alpha*exp((CloudR(i)-Rplanet)/(6d0*Cloud(ii)%fsed_beta*CloudHp(i)))+1d-2
			vsed(i)=-fsed*Kd(i)/CloudHp(i)
			lmfp=CloudMMW(i)*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
			rpart(i)=min( -vsed(i)*(Clouddens(i)*vth(i))/(rho_av(i)*CloudG(i)),
     &				 sqrt(-vsed(i)*(Clouddens(i)*vth(i))*(9d0*lmfp)/(4d0*rho_av(i)*CloudG(i))) )
			if(rpart(i).lt.Cloud(ii)%rnuc) then
				rpart(i)=Cloud(ii)%rnuc
				vsed(i)=-rpart(i)*rho_av(i)*CloudG(i)/(Clouddens(i)*vth(i))
				lmfp=CloudMMW(i)*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
				vsed(i)=vsed(i)*max(1d0,4d0*rpart(i)/(9d0*lmfp))
			endif
			mpart(i)=rho_av(i)*4d0*pi*(rpart(i)**3)/3d0
			if(iter.eq.1) then
				tot=0d0
				do iCS=1,nCS
					do iVS=1,nVS
						tot=tot+xv_bot(iVL(i,iCS))*v_cloud(iCS,iVS)*muV(iVS)/muV(iVL(i,iCS))
					enddo
				enddo
			else
				tot=sum(xc(1:nCS,i))
			endif
			xn(i)=xn(i)*f+(tot/mpart(i))*(1d0-f)
		else
			vsed(i)=-rpart(i)*rho_av(i)*CloudG(i)/(Clouddens(i)*vth(i))
			lmfp=CloudMMW(i)*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
			vsed(i)=vsed(i)*max(1d0,4d0*rpart(i)/(9d0*lmfp))
			mpart(i)=rho_av(i)*4d0*pi*rpart(i)**3/3d0
		endif
		Sat(i,1:nCS)=Sat0(i,1:nCS)*fSat(i,1:nCS)
		do iCS=1,nCS
			if(Sat(i,iCS).gt.1d0) Sat(i,iCS)=exp(log(Sat(i,iCS))*fscale)
		enddo

		do iCS=1,nCS
			vthv(i)=sqrt(8d0*kb*CloudT(i)/(pi*muV(iVL(i,iCS))*mp))

			Dp=kb*CloudT(i)*vthv(i)/(3d0*CloudP(i)*1d6*sigmamol)
			Sc(i,iCS)=fstick*Clouddens(i)**2*(muC(iCS)/(muV(iVL(i,iCS))*v_cloud(iCS,iVL(i,iCS))))*
     &						pi*rpart(i)*min(rpart(i)*vthv(i),4d0*Dp)
			if(do_nuc(iCS)) then
				tot1=Sat(i,iCS)*xv(iVL(i,iCS),i)
				if(tot1.gt.2d0) tot1=2d0*exp(log(tot1/2d0)*fscale)
				tot2=CloudP(i)*CloudMMW(i)/(muV(iVL(i,iCS))*kb*CloudT(i))
				select case(CSname(iCS))
					case('SiO','TiO2')
						call ComputeJ_xv(xv(iVL(i,iCS),i),tot2,CloudT(i),tot1,Jn_temp,A_J(iCS),B_J(iCS))
						Nc_nuc(i,iCS)=m_nuc/(muC(iCS)*mp)
					case default
						call ComputeJ(CloudT(i),tot1,tot2,xv(iVL(i,iCS),i),vthv(i),sigma_nuc(iCS),r0_nuc(iCS),Nf_nuc(iCS),Jn_temp,Nc_nuc(i,iCS))
				end select
				Jn_xv(i,iCS)=Jn_xv(i,iCS)*0.5+Jn_temp*(1d0-0.5)*fscale
			else
				Jn_xv(i,iCS)=0d0
			endif
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	fscale=fscale*pscale
	if(fscale.gt.1d0) then
		fscale=1d0
		nfscale=nfscale+1
	endif
	Nc_nuc=Nc_nuc*mp
	

C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================
c equations for material
C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================

	x=0d0

	AB=0d0

	KL=NStot
	KU=NStot

	j=0
	i=1

	if(.not.Cloud(ii)%usefsed) then
		j=j+1
		if(Cloud(ii)%freeflow_nuc) then
			x(j)=0d0
	
			dztot=(CloudR(i+1)-CloudR(i))
			dz=(CloudR(i)-CloudR(i+1))
	
			ik=KL+KU+1+j-ixn(i+1)
			AB(ik,ixn(i+1))=AB(ik,ixn(i+1))-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
	
			ik=KL+KU+1+j-ixn(i)
			AB(ik,ixn(i))=AB(ik,ixn(i))+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
			AB(ik,ixn(i))=AB(ik,ixn(i))+Clouddens(i)*vsed(i)/dztot
		else
			ik=KL+KU+1+j-ixn(i)
			AB(ik,ixn(i))=1d0
			x(j)=Cloud(ii)%xm_bot
		endif
	endif
	do iCS=1,nCS
		j=j+1
		if(Cloud(ii)%freeflow_con) then
c assume continuous flux at the bottom (dF/dz=Sc=0)
			x(j)=0d0

			dztot=(CloudR(i+1)-CloudR(i))
			dz=(CloudR(i)-CloudR(i+1))

			ik=KL+KU+1+j-ixc(iCS,i+1)
			AB(ik,ixc(iCS,i+1))=AB(ik,ixc(iCS,i+1))-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot

			ik=KL+KU+1+j-ixc(iCS,i)
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+Clouddens(i)*vsed(i)/dztot
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
		else
			ik=KL+KU+1+j-ixc(iCS,i)
			AB(ik,ixc(iCS,i))=1d0
		endif
	enddo
	do iVS=1,nVS
		if(v_include(iVS)) then
			j=j+1
			ik=KL+KU+1+j-ixv(iVS,i)
			AB(ik,ixv(iVS,i))=1d0
			x(j)=xv_bot(iVS)
		endif
	enddo
	do i=2,nnr-1
		if(.not.Cloud(ii)%usefsed) then
			j=j+1
	
			dztot=(CloudR(i+1)-CloudR(i-1))/2d0
			dz=(CloudR(i)-CloudR(i+1))
			ik=KL+KU+1+j-ixn(i+1)
			AB(ik,ixn(i+1))=AB(ik,ixn(i+1))-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
	
			ik=KL+KU+1+j-ixn(i)
			AB(ik,ixn(i))=AB(ik,ixn(i))+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot
	
			dz=(CloudR(i-1)-CloudR(i))
			ik=KL+KU+1+j-ixn(i-1)
			AB(ik,ixn(i-1))=AB(ik,ixn(i-1))-0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz/dztot
	
			ik=KL+KU+1+j-ixn(i)
			AB(ik,ixn(i))=AB(ik,ixn(i))+(0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz)/dztot
			AB(ik,ixn(i))=AB(ik,ixn(i))+Clouddens(i)*vsed(i)/dztot
	
			if(haze_nuc) then
				ik=KL+KU+1+j-ixv(ihaze,i)
				AB(ik,ixv(ihaze,i))=AB(ik,ixv(ihaze,i))+(CloudMMW(i)/muV(ihaze))*Sn(i)/m_nuc
			else
				x(j)=-Sn(i)/m_nuc
			endif
			
			do iCS=1,nCS
				if(do_nuc(iCS)) then
					ik=KL+KU+1+j-ixv(iVL(i,iCS),i)
					AB(ik,ixv(iVL(i,iCS),i))=AB(ik,ixv(iVL(i,iCS),i))+Jn_xv(i,iCS)
				endif
			enddo
c coagulation
			if(Cloud(ii)%coagulation) then
				npart=xn(i)*Clouddens(i)
				lmfp=CloudMMW(i)*mp/(sqrt(2d0)*Clouddens(i)*sigmamol)
				vmol=0.5d0*lmfp*vth(i)
				Dp=kb*CloudT(i)/(6d0*pi*rpart(i)*vmol*Clouddens(i))
				vBM=sqrt(16d0*kb*CloudT(i)/(pi*mpart(i)))
				if(Dp/rpart(i).lt.vBM) vBM=Dp/rpart(i)
				tcoaginv=npart*pi*rpart(i)**2*(abs(vsed(i))+2d0*vBM)
	
				if(.not.tcoaginv.gt.0d0) tcoaginv=0d0
				ffrag=1d0-2d0/(1d0+exp(-10d0*((abs(vsed(i))+2d0*vBM)-vfrag)/vfrag))
				tcoaginv=tcoaginv*ffrag
	
				tcinv(iter,i)=tcoaginv
				tcoaginv=sum(tcinv(max(1,iter-10):iter,i))/real(iter-max(1,iter-10)+1)
	
				ik=KL+KU+1+j-ixn(i)
				AB(ik,ixn(i))=AB(ik,ixn(i))-Clouddens(i)*tcoaginv
			endif
		endif
		do iCS=1,nCS
			j=j+1
	
			dztot=(CloudR(i+1)-CloudR(i-1))/2d0
			dz=(CloudR(i)-CloudR(i+1))
			ik=KL+KU+1+j-ixc(iCS,i+1)
			AB(ik,ixc(iCS,i+1))=AB(ik,ixc(iCS,i+1))-(Clouddens(i+1)*vsed(i+1)+0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot

			ik=KL+KU+1+j-ixc(iCS,i)
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+(0.5d0*(Kd(i+1)*Clouddens(i+1)+Kd(i)*Clouddens(i))/dz)/dztot

			dz=(CloudR(i-1)-CloudR(i))
			ik=KL+KU+1+j-ixc(iCS,i-1)
			AB(ik,ixc(iCS,i-1))=AB(ik,ixc(iCS,i-1))-0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz/dztot

			ik=KL+KU+1+j-ixc(iCS,i)
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+(0.5d0*(Kd(i-1)*Clouddens(i-1)+Kd(i)*Clouddens(i))/dz)/dztot
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+Clouddens(i)*vsed(i)/dztot
	
			if(iamorph(iCS).ne.0) then
				ik=KL+KU+1+j-ixc(iamorph(iCS),i)
				AB(ik,ixc(iamorph(iCS),i))=AB(ik,ixc(iamorph(iCS),i))+Clouddens(i)*tcrystinv(i)
			else
				ik=KL+KU+1+j-ixv(iVL(i,iCS),i)
				AB(ik,ixv(iVL(i,iCS),i))=AB(ik,ixv(iVL(i,iCS),i))+Sc(i,iCS)*xn(i)
			endif
			if(icryst(iCS).ne.0) then
				ik=KL+KU+1+j-ixc(iCS,i)
				AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))-Clouddens(i)*tcrystinv(i)
			endif

			if(do_nuc(iCS)) then
				ik=KL+KU+1+j-ixv(iVL(i,iCS),i)
				AB(ik,ixv(iVL(i,iCS),i))=AB(ik,ixv(iVL(i,iCS),i))+Jn_xv(i,iCS)*Nc_nuc(i,iCS)*muC(iCS)
			endif

			ik=KL+KU+1+j-ixc(iCS,i)
			AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))-Sc(i,iCS)/(Sat(i,iCS)*mpart(i))
		enddo

		do iVS=1,nVS
			if(v_include(iVS)) then
				j=j+1
	
				dztot=(CloudR(i+1)-CloudR(i-1))/2d0
				dz=(CloudR(i)-CloudR(i+1))
				ik=KL+KU+1+j-ixv(iVS,i+1)
				AB(ik,ixv(iVS,i+1))=AB(ik,ixv(iVS,i+1))-(0.5d0*(Kg(i+1)*Clouddens(i+1)+Kg(i)*Clouddens(i))/dz)/dztot

				ik=KL+KU+1+j-ixv(iVS,i)
				AB(ik,ixv(iVS,i))=AB(ik,ixv(iVS,i))+(0.5d0*(Kg(i+1)*Clouddens(i+1)+Kg(i)*Clouddens(i))/dz)/dztot

				dz=(CloudR(i-1)-CloudR(i))
				ik=KL+KU+1+j-ixv(iVS,i-1)
				AB(ik,ixv(iVS,i-1))=AB(ik,ixv(iVS,i-1))-0.5d0*(Kg(i-1)*Clouddens(i-1)+Kg(i)*Clouddens(i))/dz/dztot

				ik=KL+KU+1+j-ixv(iVS,i)
				AB(ik,ixv(iVS,i))=AB(ik,ixv(iVS,i))+(0.5d0*(Kg(i-1)*Clouddens(i-1)+Kg(i)*Clouddens(i))/dz)/dztot
	
				do iCS=1,nCS
					if(iamorph(iCS).eq.0) then
						ik=KL+KU+1+j-ixv(iVL(i,iCS),i)
						AB(ik,ixv(iVL(i,iCS),i))=AB(ik,ixv(iVL(i,iCS),i))-Sc(i,iCS)*xn(i)*v_cloud(iCS,iVS)*(muV(iVS)/muC(iCS))
					endif
					ik=KL+KU+1+j-ixc(iCS,i)
					AB(ik,ixc(iCS,i))=AB(ik,ixc(iCS,i))+Sc(i,iCS)*v_cloud(iCS,iVS)*(muV(iVS)/muC(iCS))/(Sat(i,iCS)*mpart(i))

					if(do_nuc(iCS)) then
						ik=KL+KU+1+j-ixv(iVL(i,iCS),i)
						AB(ik,ixv(iVL(i,iCS),i))=AB(ik,ixv(iVL(i,iCS),i))-Jn_xv(i,iCS)*Nc_nuc(i,iCS)*v_cloud(iCS,iVS)*muV(iVS)
					endif
				enddo

				if(haze_nuc.and.iVS.eq.ihaze) then
					ik=KL+KU+1+j-ixv(ihaze,i)
					AB(ik,ixv(ihaze,i))=AB(ik,ixv(ihaze,i))-(CloudMMW(i)/muV(ihaze))*Sn(i)
				endif

			endif
		enddo
	enddo
	i=nnr
	dz=CloudR(i)-CloudR(i-1)

	if(.not.Cloud(ii)%usefsed) then
		j=j+1
		ik=KL+KU+1+j-ixn(i)
		AB(ik,ixn(i))=Kd(i)/dz-vsed(i)
	
		ik=KL+KU+1+j-ixn(i-1)
		AB(ik,ixn(i-1))=-Kd(i)/dz
	
		x(j)=0d0
c		x(j)=F_IDP/((4d0*pi*(0.1e-4)**3*rho_av(i)/3d0)*Clouddens(i))
	endif
	do iCS=1,nCS
		j=j+1
		ik=KL+KU+1+j-ixc(iCS,i)
		AB(ik,ixc(iCS,i))=Kd(i)/dz-vsed(i)

		ik=KL+KU+1+j-ixc(iCS,i-1)
		AB(ik,ixc(iCS,i-1))=-Kd(i)/dz

		x(j)=0d0
c		if(CSname(iCS).eq.'MgSiO3') then
c			x(j)=F_IDP/(Clouddens(i))
c		endif
	enddo
	do iVS=1,nVS
		if(v_include(iVS)) then
			j=j+1
			ik=KL+KU+1+j-ixv(iVS,i)
			AB(ik,ixv(iVS,i))=Kg(i)/dz

			ik=KL+KU+1+j-ixv(iVS,i-1)
			AB(ik,ixv(iVS,i-1))=-Kg(i)/dz

			x(j)=0d0
		endif
	enddo
	
10	continue
	NRHS=1
	info=0

	j=2*KL+KU+1
	x0=x
	call DGBSV(NN,KL,KU,NRHS,AB,j,IWORK,x,NN,INFO)	

	do i=1,nnr
		if(.not.Cloud(ii)%usefsed) then
			if(.not.x(ixn(i)).gt.1d-50/m_nuc) x(ixn(i))=1d-50/m_nuc
		endif
		do iCS=1,nCS
			if(.not.x(ixc(iCS,i)).gt.0d0) x(ixc(iCS,i))=0d0
		enddo
		do iVS=1,nVS
			if(v_include(iVS)) then
				if(.not.x(ixv(iVS,i)).gt.1d-50) x(ixv(iVS,i))=1d-50
				if(.not.x(ixv(iVS,i)).lt.xv_bot(iVS)) x(ixv(iVS,i))=xv_bot(iVS)
			endif
		enddo
	enddo

	xc_iter(iter,1:nCS,1:nnr)=xc(1:nCS,1:nnr)
	xv_iter(iter,1:nVS,1:nnr)=xv(1:nVS,1:nnr)
	xn_iter(iter,1:nnr)=xn(1:nnr)

	maxerr=0d0
	do i=1,nnr
		if(i.gt.1) then
			if(.not.Cloud(ii)%usefsed) then
				err=abs(xn(i)-x(ixn(i)))/(xn(i)+x(ixn(i)))
				if(err.gt.maxerr.and.tot.gt.1d-20.and.(xn(i)*m_nuc.gt.1d-20.or.x(ixn(i))*m_nuc.gt.1d-20)) then
					maxerr=err
				endif
			endif
			do iCS=1,nCS
				err=abs(xc(iCS,i)-x(ixc(iCS,i)))/(xc(iCS,i)+x(ixc(iCS,i)))
				if(err.gt.maxerr.and.tot.gt.1d-20.and.(xc(iCS,i).gt.1d-20.or.x(ixc(iCS,i)).gt.1d-20)) then
					maxerr=err
				endif
			enddo
			do iVS=1,nVS
				if(v_include(iVS)) then
					err=abs(xv(iVS,i)-x(ixv(iVS,i)))/(xv(iVS,i)+x(ixv(iVS,i)))
					if(err.gt.maxerr.and.tot.gt.1d-20.and.(xv(iVS,i).gt.1d-20.or.x(ixv(iVS,i)).gt.1d-20)) then
						maxerr=err
					endif
				endif
			enddo
		endif
		if(.not.Cloud(ii)%usefsed) then
			if(iter.eq.1) then
				xn(i)=x(ixn(i))
			else
				xn(i)=xn(i)*f+x(ixn(i))*(1d0-f)
			endif
		endif
		do iCS=1,nCS
			if(iter.eq.1) then
				xc(iCS,i)=x(ixc(iCS,i))
			else
				xc(iCS,i)=xc(iCS,i)*f+x(ixc(iCS,i))*(1d0-f)
			endif
		enddo
		do iVS=1,nVS
			if(v_include(iVS)) then
				if(iter.eq.1) then
					xv(iVS,i)=x(ixv(iVS,i))
				else
					if(Cloud(ii)%computeJn.or.Cloud(ii)%usefsed) then
						xv(iVS,i)=xv(iVS,i)**f*x(ixv(iVS,i))**(1d0-f)
					else
						xv(iVS,i)=xv(iVS,i)*f+x(ixv(iVS,i))*(1d0-f)
					endif
				endif
			endif
		enddo
	enddo

C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================
C=========================================================================================


	do i=nnr,1,-1
		tot=0d0
		do iCS=1,nCS
			tot=tot+xc(iCS,i)/rhodust(iCS)
		enddo
		if(tot.gt.0d0) then
			rho_av(i)=sum(xc(1:nCS,i))/tot
		else if(i.lt.nnr) then
			rho_av(i)=rho_av(i+1)
		else
			rho_av(i)=sum(rhodust(1:nCS))/real(nCS)
		endif
		tot=sum(xc(1:nCS,i))
		if(xn(i).gt.0d0) then
			if(Cloud(ii)%computeJn) then
				rr=(3d0*(tot/xn(i))/(4d0*pi*rho_av(i)))**(1d0/3d0)
				if(.not.rr.ge.1e-8) rr=1e-8
			else
				rr=(3d0*(tot/xn(i))/(4d0*pi*rho_av(i))+Cloud(ii)%rnuc**3)**(1d0/3d0)
				if(.not.rr.ge.Cloud(ii)%rnuc) rr=Cloud(ii)%rnuc
			endif
		else
			rr=Cloud(ii)%rnuc
		endif
		err=abs(rr-rpart(i))/(rr+rpart(i))
		if(err.gt.maxerr.and.tot.gt.1d-20.and.i.gt.1.and..not.Cloud(ii)%usefsed) then
			maxerr=err
		endif
		if(.not.Cloud(ii)%usefsed) rpart(i)=rr!sqrt(rpart(i)*rr)
	enddo
	if(Cloud(ii)%computeJn.and.maxerr.lt.eps.and.nfscale.gt.100.and.iconv.ge.nconv.and.eps.gt.5d-3) then
		iconv=0
		eps=1d-3
		f=0.999
	endif
	if(maxerr.lt.eps.and.nfscale.gt.100) then
		iconv=iconv+1
		if(iconv.gt.nconv) exit
	else
		iconv=0
	endif
c	print*,iter,maxerr
20	continue
	enddo
c end the loop





	if(iter.gt.niter.and..not.retrieval) then
		if(iconv.eq.0.and.(nTiter.gt.4.or..not.computeT)) print*,'Cloud formation not converged: ',maxerr
		iter=niter
	endif
	xn=0d0
	xc=0d0
	xv=0d0
	do i=iter-nconv+1,iter
		xn(1:nnr)=xn(1:nnr)+xn_iter(i,1:nnr)/real(nconv)
		xc(1:nCS,1:nnr)=xc(1:nCS,1:nnr)+xc_iter(i,1:nCS,1:nnr)/real(nconv)
		xv(1:nVS,1:nnr)=xv(1:nVS,1:nnr)+xv_iter(i,1:nVS,1:nnr)/real(nconv)
	enddo
	do i=nnr,1,-1
		tot=0d0
		do iCS=1,nCS
			tot=tot+xc(iCS,i)/rhodust(iCS)
		enddo
		if(xn(i).gt.0d0) then
			if(Cloud(ii)%haze) then
				rr=(3d0*(tot/xn(i))/(4d0*pi))**(1d0/3d0)
			else
				rr=(3d0*(tot/xn(i))/(4d0*pi)+Cloud(ii)%rnuc**3)**(1d0/3d0)
			endif
			if(.not.rr.ge.Cloud(ii)%rnuc) then
				rr=Cloud(ii)%rnuc
				xn(i)=3d0*(tot/(rr**3))/(4d0*pi)
			endif
		else
			rr=Cloud(ii)%rnuc
			xn(i)=3d0*(tot/(rr**3))/(4d0*pi)
		endif
		if(tot.gt.0d0) then
			rho_av(i)=sum(xc(1:nCS,i))/tot
		else if(i.lt.nnr) then
			rho_av(i)=rho_av(i+1)
		else
			rho_av(i)=sum(rhodust(1:nCS))/real(nCS)
		endif
		if(.not.Cloud(ii)%usefsed) rpart(i)=rr
	enddo

	deallocate(vthv)
	deallocate(Sc)
	deallocate(IWORK)
	deallocate(AB)

	allocate(dx(nnr))
	logP(1:nr)=-log(P(1:nr))
	logCloudP(1:nnr)=-log(CloudP(1:nnr))

	x(1:nnr)=0d0
	do iCS=1,nCS
		x(1:nnr)=x(1:nnr)+xc(iCS,1:nnr)*Clouddens(1:nnr)
	enddo
	call regridarray(logCloudP,x,nnr,logP,cloud_dens(1:nr,ii),nr)

	x(1:nnr)=rpart(1:nnr)
	call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%rv(1:nr),nr)

	Cloud(ii)%frac(1:nr,1:Cloud(ii)%nmat)=0d0

	do iCS=1,nCS
		x(1:nnr)=xc(iCS,1:nnr)
		call regridarray(logCloudP,x,nnr,logP,Cloud(ii)%frac(1:nr,iCS),nr)
	enddo

	if(.not.retrieval) then
		do i=1,nnr
			do iCS=1,nCS
				if(do_nuc(iCS)) then
					Jn_xv(i,iCS)=Jn_xv(i,iCS)*xv(iVL(i,iCS),i)*Nc_nuc(i,iCS)*muC(iCS)/m_nuc
					Sn(i)=Sn(i)+Jn_xv(i,iCS)
				endif
			enddo
		enddo
		if(do3D) then
			open(unit=20,file=trim(outputdir) // '/cloudstructure' // trim(int2string(i3D,'(i0.4)')) // '.dat',
     &             FORM="FORMATTED",ACCESS="STREAM")
		else
			open(unit=20,file=trim(outputdir) // '/cloudstructure.dat',FORM="FORMATTED",ACCESS="STREAM")
		endif
		allocate(v_names_out(nVS),xv_out(nVS))
		j=0
		do iVS=1,nVS
			if(v_include(iVS)) then
				j=j+1
				v_names_out(j)=v_names(iVS)
			endif
		enddo
		form='("#",a18,a19,a19,' // trim(int2string(nCS+j+1,'(i4)')) // 'a23,a19,a19,a19,a19)'
		write(20,form) "P[bar]","dens[g/cm^3]","xn",(trim(CSname(i))//"[s]",i=1,nCS),
     &				(trim(v_names_out(i))//"[v]",i=1,j),"r[cm]","T[K]","Jstar","fsed"
		form='(es19.7E3,es19.7E3,es19.7E3,' // trim(int2string(nCS+j,'(i4)')) // 'es23.7E3,es19.7E3,es19.7E3,es19.7E3,es19.7E3)'
		do i=1,nnr
			j=0
			do iVS=1,nVS
				if(v_include(iVS)) then
					j=j+1
					xv_out(j)=xv(iVS,i)
				endif
			enddo
			write(20,form) CloudP(i),Clouddens(i),xn(i)*m_nuc,xc(1:nCS,i),xv_out(1:j),rpart(i),
     &							CloudT(i),Sn(i),-vsed(i)*CloudHp(i)/Kd(i)
		enddo
		close(unit=20)
	endif

	do i=1,nr
		tot=sum(Cloud(ii)%frac(i,1:Cloud(ii)%nmat))
		if(tot.gt.0d0) then
			Cloud(ii)%frac(i,1:Cloud(ii)%nmat)=Cloud(ii)%frac(i,1:Cloud(ii)%nmat)/tot
		else
			Cloud(ii)%frac(i,1:Cloud(ii)%nmat)=1d0/real(Cloud(ii)%nmat)
			cloud_dens(i,ii)=0d0
		endif
	enddo


c Elemental abundances
	allocate(at_ab(nr,N_atoms))
	at_ab=0d0
	do iVS=1,nVS
		if(v_include(iVS)) then
			x(1:nnr)=xv(iVS,1:nnr)
			call regridarray(logCloudP,x,nnr,logP,logx,nr)
			do j=1,N_atoms
				at_ab(1:nr,j)=at_ab(1:nr,j)+logx(1:nr)*mutot*v_atoms(iVS,j)/muV(iVS)
			enddo
		endif
	enddo

	call cpu_time(time)
	timecloud=timecloud+time
	call system_clock(itime)
	itimecloud=itimecloud+itime

c	open(unit=20,file=trim(outputdir) // '/atoms.dat',FORM="FORMATTED",ACCESS="STREAM")
	if(dochemistry) then
		dochemR=.false.
		dochemR(1)=.true.
		dochemR(nr)=.true.
		do i=1,nr,nrstepchem
			dochemR(i)=.true.
		enddo
		ini=.true.
		do i=1,nr
			call tellertje(i,nr)
			if(dochemR(i)) then
			molfracs_atoms(1:N_atoms)=at_ab(i,1:N_atoms)
			molfracs_atoms=molfracs_atoms+molfracs_atoms0
			molfracs_atoms(3)=molfracs_atoms(3)+COabun
			molfracs_atoms(5)=molfracs_atoms(5)+COabun
			do j=1,N_atoms
				if(.not.molfracs_atoms(j).gt.0d0) then
					molfracs_atoms(j)=0d0
				endif
			enddo
			tot=sum(molfracs_atoms(1:N_atoms))
			molfracs_atoms=molfracs_atoms/tot
			do j=1,N_atoms
				if(.not.molfracs_atoms(j).gt.1d-50) then
					molfracs_atoms(j)=1d-50
				endif
			enddo
			if(nPhotoReacts.gt.0) call doPhotoChemAtom(i)
			if((P(i).ge.mixP.or.i.eq.1).and.dochemistry) then
				call call_chemistry(T(i),P(i),mixrat_r(i,1:nmol),molname(1:nmol),nmol,ini,.false.,cloudspecies,
     &				XeqCloud(i,1:nclouds),nclouds,nabla_ad(i),MMW(i),didcondens(i),includemol,.false.)
			else
				mixrat_r(i,1:nmol)=mixrat_r(i-1,1:nmol)
				XeqCloud(i,1:nclouds)=XeqCloud(i-1,1:nclouds)
				nabla_ad(i)=nabla_ad(i-1)
				MMW(i)=MMW(i-1)
				didcondens(i)=didcondens(i-1)
			endif
c			write(20,*) P(i),molfracs_atoms(1:N_atoms)
			endif
		enddo
		if(nrstepchem.ne.1) then
			do i=1,nmol
				if(includemol(i).or.diseqmol(i)) then
					call fillblanks(P,mixrat_r(1:nr,i),nr,dochemR,.true.)
				endif
			enddo
			call fillblanks(P,MMW,nr,dochemR,.true.)
			call fillblanks(P,nabla_ad,nr,dochemR,.true.)
		endif
		if(fixMMW) MMW=MMW0
c		close(unit=20)
		if(disequilibrium) then
			call AddDiseqAtoms(Otot,Ctot,Ntot)
c			call disequilibrium code
c			input: 	R(1:nr+1) : These are the radial boundaries of the layers (bottom to top)
c			P(1:nr),T(1:nr) : These are the pressure and temperature inside the layers
c			molname(1:nmol) : names of the molecules included
c			Kzz_r(1:nr) : Diffusion coefficient
c			input/output:	mixrat_r(1:nr,1:nmol) : number densities inside each layer. Now set to equilibrium abundances.
			call output("==================================================================")
			call output("Computing disequilibrium chemistry")
			call diseq_calc(nr,R(1:nr+1),P(1:nr),T(1:nr),nmol,molname(1:nmol),mixrat_r(1:nr, 1:nmol),COratio,Kzz_g(1:nr))
			call CorrectDiseqAtoms(Otot,Ctot,Ntot)
		endif
	else
		do iVS=1,nVS
			if(v_include(iVS)) then
				do k=1,nmol
					if(includemol(k)) then
						if(v_names(iVS).eq.molname(k)) then
							x(1:nnr)=xv(iVS,1:nnr)
							call regridarray(logCloudP,x,nnr,logP,logx,nr)
							mixrat_r(1:nr,k)=logx(1:nr)*mutot/(muV(iVS))
						endif
					endif
				enddo
			endif
		enddo
	endif
	do i=1,nr
		tot=0d0
		do j=1,nmol
			if(mixrat_r(i,j).gt.0d0) tot=tot+mixrat_r(i,j)
		enddo
		if(tot.gt.0d0) mixrat_r(i,1:nmol)=mixrat_r(i,1:nmol)/tot
	enddo
	deallocate(at_ab)

	do i=1,nr
		do j=1,nmol
			if(.not.mixrat_r(i,j).gt.0d0) mixrat_r(i,j)=0d0
		enddo
	enddo
	if(nfixmol.gt.0) then
		do i=1,nfixmol
			mixrat_r(1:nr,ifixmol(i))=fixmol_abun(i)*exp(-P(1:nr)/fixmol_P(i))
		enddo
	endif
	call doPhotoChemMol()
	do j=1,nmol
		if(isotope(j).gt.0) then
			do i=1,nr
				mixrat_r(i,isotope(j))=mixrat_r(i,j)/f_isotope(j)
				mixrat_r(i,j)=mixrat_r(i,j)*(1d0-1d0/f_isotope(j))
			enddo
		endif
	enddo
	
	if(.not.retrieval) then
		if(complexKzz) then
			open(unit=50,file=trim(outputdir) // 'cloudKzz.dat',FORM="FORMATTED",ACCESS="STREAM")
			form='("#",a12,a13,a13,a13)'
			write(50,trim(form)) "Kzz [cm^2/s]","P [bar]","Kpart","Kgas"
			form='(es13.3E3,es13.3E3,es13.3E3,es13.3E3)'
			do i=1,nnr
				write(50,trim(form)) Km(i),CloudP(i),Kd(i),Kg(i)
			enddo
			close(unit=50)
		endif
	endif


	deallocate(Sat,docondense)
	deallocate(mpart)
	deallocate(rho_av)
	deallocate(y)
	deallocate(Sn)
	deallocate(vth)
	deallocate(tcinv,xn_iter,xc_iter,xv_iter)
	deallocate(vsed)
	deallocate(ixv)
	deallocate(ixc)
	deallocate(x,dx)
	deallocate(logCloudP)
	deallocate(Kd,Kg,Km)
	if(Cloud(ii)%hazetype.eq.'optEC') deallocate(CloudtauUV,CloudkappaUV)

	return
	end

	subroutine ComputeJ_xv(xv,scale,T,Sat,J,A,B)
	IMPLICIT NONE
	real*8 xv,scale,T,Sat,J,A,B
	
	if(Sat.lt.1d0) then
		J=0d0
	else
		J=xv*scale**2*exp(B-A/(T**3*(log(Sat))**2))
	endif
	if(.not.J.gt.0d0) J=0d0
	
	return
	end


	subroutine ComputeJ(T,Sat,scale,xv,vth,sigma,r0,Nf,J,Nc)
	IMPLICIT NONE
	real*8 T,nx,sigma,Sat,vth,J
	real*8 Nstar1,theta,r0,A0,logSat,Nstar
	real*8 Tmax,mu,tau,MMW,scale,xv
	real*8 pi,Z,dGRT,Nstar_inf,Nf,kb,mp,Nc,x0,x1,x2,x3,dgdn,thetaN,nst
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(kb=1.3806503d-16)
	parameter(mp=1.660539040d-24)	!atomic mass unit
	
	if(Sat.le.1d0) then
		J=0d0
		Nc=1d0
		return
	endif
	nx=scale*xv
	
	logSat=log(Sat)

	A0=4d0*pi*r0**2

	theta=A0*sigma/(kb*T)

	Nstar_inf=(2d0*theta/(3d0*logSat))**3
	Nstar1=(Nstar_inf/8d0)*(1d0+sqrt(1d0+2d0*(Nf/Nstar_inf)**(1d0/3d0))-2d0*(Nf/Nstar_inf)**(1d0/3d0))**3

	if(.not.Nstar1.gt.1d-6) Nstar1=1d-6

	Z=sqrt(theta*(2d0*Nstar1**(1d0/3d0)+3d0*Nf**(1d0/3d0))/(6d0*pi))/(Nstar1**(1d0/3d0)+Nf**(1d0/3d0))

	dGRT=theta*Nstar1/(Nstar1**(1d0/3d0)+Nf**(1d0/3d0))

	tau=nx*vth*(Nstar1+1d0)**(2./3.)*A0

	J=scale*tau*Z*exp(Nstar1*logSat-dGRT)
	Nc=Nstar1+1d0

      A0 = 4.d0*pi*r0**2
      theta = A0*sigma/kb
      x0     = 2.d0*theta/(3.d0*T*logSat)
      x1     = Nf**(1.d0/3.d0)
      x2     = x1/x0
      x3     = 0.5d0*(1.d0+DSQRT(1.d0+2.d0*x2)) - x2
      Nstar  = 1.d0 + (x0*x3)**3 
      if (Nstar<=1.d0) Nstar=1.000000001d0
 
 
      x0     = x0*logSat
      x2     = (Nstar-1.d0)**(1.d0/3.d0) + x1
      x3     = 1.d0/(Nstar-1.d0)**(2.d0/3.d0)
      dgdn   = x0*x3* ( 1.d0/x2**2 + x1/x2**3 ) / 3.d0
      Z = SQRT(dgdn/(2.d0*pi))
      thetaN = theta/(1.d0+(Nf/(Nstar-1.d0))**(1.d0/3.d0))
      x1     = (Nstar-1.d0)*logSat - (thetaN/T)
     &         *(Nstar-1.d0)**(2.d0/3.d0)
      nst    = scale*EXP(x1)
 
      J = tau*nst*Z
      Nc=Nstar

	if(.not.J.gt.0d0.or..not.Nc.gt.1d0) then
		J=0d0
		Nc=1d0
	endif

	return	
	end


	
	subroutine AddDiseqAtoms(Otot,Ctot,Ntot)	
	use GlobalSetup
	IMPLICIT NONE
	real*8 Otot(nr),Ctot(nr),Ntot(nr)
	integer iCO,iCO2,iH2O,iCH4,iN2,iNH3
	
	iH2O=1
	iCO2=2
	iCO=5
	iCH4=6
	iNH3=11
	iN2=22
	
	Otot(1:nr)=mixrat_r(1:nr,iH2O)+2.0*mixrat_r(1:nr,iCO2)+mixrat_r(1:nr,iCO)
	Ctot(1:nr)=mixrat_r(1:nr,iCO2)+mixrat_r(1:nr,iCO)+mixrat_r(1:nr,iCH4)
	Ntot(1:nr)=mixrat_r(1:nr,iNH3)+2.0*mixrat_r(1:nr,iN2)
	
	return
	end

	subroutine CorrectDiseqAtoms(Otot,Ctot,Ntot)	
	use GlobalSetup
	IMPLICIT NONE
	real*8 Otot(nr),Ctot(nr),Ntot(nr),f
	integer iCO,iCO2,iH2O,iCH4,iN2,iNH3
	integer i
	
	iH2O=1
	iCO2=2
	iCO=5
	iCH4=6
	iNH3=11
	iN2=22
	
	do i=1,nr
		f=(Otot(i)-mixrat_r(i,iCO))/(2.0*mixrat_r(i,iCO2)+mixrat_r(i,iH2O))
		if(f.lt.0d0) f=0d0
		mixrat_r(i,iCO2)=f*mixrat_r(i,iCO2)
		mixrat_r(i,iH2O)=f*mixrat_r(i,iH2O)
		f=(Ctot(i)-mixrat_r(i,iCO2))/(mixrat_r(i,iCH4)+mixrat_r(i,iCO))
		if(f.lt.0d0) f=0d0
		mixrat_r(i,iCO)=f*mixrat_r(i,iCO)
		mixrat_r(i,iCH4)=f*mixrat_r(i,iCH4)
		f=Ntot(i)/(mixrat_r(i,iNH3)+2.0*mixrat_r(i,iN2))
		mixrat_r(i,iNH3)=f*mixrat_r(i,iNH3)
		mixrat_r(i,iN2)=f*mixrat_r(i,iN2)
	enddo
	
	return
	end
	

	subroutine PvapNH3(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0,c1,c2
	parameter(c0=10.53,c1=-2161.0,c2=-86596.0)

	Pvap=exp(c0+c1/T+c2/T**2)
	liquid=.true.

	return
	end
	

	subroutine PvapSiO(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0,c1
	parameter(c0=-4.95200E+04,c1=3.25200E+01)

	Pvap=1d-6*exp(c0/T+c1)
	liquid=.true.

	return
	end
	
	

	subroutine PvapNaCl(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0_i,c1_i,c2_i,c3_i,c4_i
	real*8 c0_l,c1_l,c2_l,c3_l,c4_l
	parameter(	
     &	c0_i=-2.79146E+04,
     &	c1_i=3.46023E+01,
     &	c2_i=-3.11287E-03,
     &	c3_i=5.30965E-07,
     &	c4_i=-2.59584E-12,
     &	c0_l=-2.48880E+04,
     &	c1_l=3.18494E+01,
     &	c2_l=-3.08748E-03,
     &	c3_l=4.84990E-07,
     &	c4_l=-2.60359E-11)

	Pi=c0_i/T+c1_i+c2_i*T+c3_i*T**2+c4_i*T**3
	Pl=c0_l/T+c1_l+c2_l*T+c3_l*T**2+c4_l*T**3

	if(Pi.lt.Pl) then
		Pvap=Pi
		liquid=.false.
	else
		Pvap=Pl
		liquid=.true.
	endif
	Pvap=1d-6*exp(Pvap)
	
	return
	end
		
	
	subroutine PvapKCl(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0_i,c1_i,c2_i,c3_i,c4_i
	real*8 c0_l,c1_l,c2_l,c3_l,c4_l
	parameter(	
     &	c0_i=-2.69250E+04,
     &	c1_i=3.39574E+01,
     &	c2_i=-2.04903E-03,
     &	c3_i=-2.83957E-07,
     &	c4_i=1.82974E-10,
     &	c0_l=-2.50293E+04,
     &	c1_l=3.39453E+01,
     &	c2_l=-4.61815E-03,
     &	c3_l=7.36857E-07,
     &	c4_l=0.00000E+00)

	Pi=c0_i/T+c1_i+c2_i*T+c3_i*T**2+c4_i*T**3
	Pl=c0_l/T+c1_l+c2_l*T+c3_l*T**2+c4_l*T**3

	if(Pi.lt.Pl) then
		Pvap=Pi
		liquid=.false.
	else
		Pvap=Pl
		liquid=.true.
	endif
	Pvap=1d-6*exp(Pvap)
		
	return
	end

	subroutine PvapCr(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0_i,c1_i,c2_i,c3_i,c4_i
	real*8 c0_l,c1_l,c2_l,c3_l,c4_l
	parameter(	
     &	c0_i=-4.78455E+04,
     &	c1_i=3.22423E+01,
     &	c2_i=-5.28710E-04,
     &	c3_i=-6.17347E-08,
     &	c4_i=2.88469E-12,
     &	c0_l=-4.47712E+04,
     &	c1_l=3.09753E+01,
     &	c2_l=-7.84094E-04,
     &	c3_l=-5.92580E-10,
     &	c4_l=1.25866E-11)

	Pi=c0_i/T+c1_i+c2_i*T+c3_i*T**2+c4_i*T**3
	Pl=c0_l/T+c1_l+c2_l*T+c3_l*T**2+c4_l*T**3

	if(Pi.lt.Pl.and.T.lt.2981d0) then
		Pvap=Pi
		liquid=.false.
	else
		Pvap=Pl
		liquid=.true.
	endif
	Pvap=1d-6*exp(Pvap)
		
	return
	end

	subroutine PvapFeS(T,Pvap,liquid)
	IMPLICIT NONE
	real*8 T,Pvap,Pl,Pi
	logical liquid
	real*8 c0_i,c1_i,c2_i,c3_i,c4_i
	real*8 c0_l,c1_l,c2_l,c3_l,c4_l
	parameter(	
     &	c0_i=-5.69922E+04,
     &	c1_i=3.86753E+01,
     &	c2_i=-4.68301E-03,
     &	c3_i=1.03559E-06,
     &	c4_i=-8.42872E-11,
     &	c0_l=-5.26135E+04,
     &	c1_l=3.46138E+01,
     &	c2_l=-3.55056E-03,
     &	c3_l=7.59195E-07,
     &	c4_l=-6.94708E-11)

	Pi=c0_i/T+c1_i+c2_i*T+c3_i*T**2+c4_i*T**3
	Pl=c0_l/T+c1_l+c2_l*T+c3_l*T**2+c4_l*T**3

	if(Pi.lt.Pl) then
		Pvap=Pi
		liquid=.false.
	else
		Pvap=Pl
		liquid=.true.
	endif
	Pvap=1d-6*exp(Pvap)
		
	return
	end

	real*8 function compGibbs(T,b,ifit)
	IMPLICIT NONE
	real*8 b(0:4),T
	integer ifit
	
	select case(ifit)
		case(0)
			compGibbs=-(b(0)/T + b(1)*log(T) + b(2) + b(3)*T + b(4)*T**2)
		case(1)
			compGibbs=b(0)/T + b(1) + b(2)*T + b(3)* T**2 + b(4)*T**3
			compGibbs=compGibbs/(1.987*T)
		case(2)
			compGibbs=b(0)/T + b(1) + b(2)*T + b(3)* T**2 + b(4)*T**3
			compGibbs=compGibbs/(8.314*T)
		case default
			print*,'unknown ifit'
			stop
	end select
	
	return
	end
	
