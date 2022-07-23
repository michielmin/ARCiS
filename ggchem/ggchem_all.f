************************************************************************
      module PARAMETERS
************************************************************************
      character(len=200) :: elements,abund_file,struc_file,adapt_file
      integer :: abund_pick,model_dim,Npoints,model_struc,verbose
      integer :: method_eqcond,Nseq
      logical :: model_eqcond,model_pconst,pick_mfrac,initchem_info
      logical :: useDataBase,remove_condensates,phyllosilicates
      logical :: auto_atmos,adapt_cond,stop_after_init
      real*8  :: Tfast,Tmin,Tmax,pmin,pmax,nHmin,nHmax,Tseq(100)
      real*8  :: Mpl,Rpl,gamma,Tmin_atmos
      end

************************************************************************
      module DUST_DATA
************************************************************************
      character(len=200) :: DustChem_file
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: NELEM=41        ! number of elements (up to Zr + W)
      integer,parameter :: NDUSTmax=500    ! max number of condensed species
      integer :: NDUST                     ! number of condensed species
      integer :: NEPS                      ! number of affected elements
      
      character(len=2)  :: elnam(NELEM)       ! names of elements
      character(len=20) :: dust_nam(NDUSTmax) ! names of dust species
      integer :: elnr(NELEM),elcode(NELEM)    ! element cross-indices
      real(kind=qp) :: eps0(NELEM)            ! element abundances
      real*8  :: mass(NELEM)                  ! element masses
      real*8  :: dust_rho(NDUSTmax)           ! dust material densities
      real*8  :: dust_mass(NDUSTmax)          ! dust monomer volume
      real*8  :: dust_vol(NDUSTmax)           ! dust monomer volume
      real*8  :: Tmelt(NDUSTmax)              ! melting points
      real*8  :: Tcorr(NDUSTmax)
      logical :: is_liquid(NDUSTmax)
      integer :: dust_nel(NDUSTmax)           ! no of elements in dust
      integer :: dust_el(NDUSTmax,8)          ! indices of elements
      integer :: dust_nu(NDUSTmax,8)          ! stoichiometric coeffs
      
      integer :: fit(NDUSTmax)                ! fit-formular identifier
      real*8  :: cfit(NDUSTmax,0:4)           ! pvap fit coefficients
      integer :: Nfit(NDUSTmax)               ! number of BURCAT fits
      real*8  :: Tfit(NDUSTmax,9)             ! BURCAT Tfit
      real*8  :: Bfit(NDUSTmax,8,14)          ! BURCAT fit coefficients
      
      real(kind=qp) :: bk=1.380662Q-16        ! Boltzman constant
      real(kind=qp) :: bar=1.Q+6              ! 1 bar in dyn/cm2
      real(kind=qp) :: amu=1.66055Q-24        ! atomar mass unit
      real(kind=qp) :: atm=1.013Q+6           ! standard atmosphere pressure
      real(kind=qp) :: rgas=8.3144598Q+0      ! gas constant 
      real(kind=qp) :: mel=9.109389754Q-28    ! electron mass
      real*8        :: grav=6.6725985D-8      ! gravitational constant
      real*8        :: MEarth=5.9722D+27      ! Earth mass [g]
      real*8        :: REarth=6371.D+5        ! Earth radius [cm]
      real(kind=qp) :: muH                    ! rho/n<H>
      end

************************************************************************
      module CHEMISTRY
************************************************************************
      use DUST_DATA,ONLY: NELEM
      character(len=200) :: dispol_file(4)
      logical :: NewFullIt
      integer :: NewBackIt,NewFastLevel,NewPreMethod
      real*8  :: NewBackFac
      integer :: NMOLdim         ! max number of molecules
      integer :: NMOLE           ! number of molecules found
      integer :: NELM            ! number of elements found
      integer :: el=0,H=0,He=0,Li=0,Be=0,B=0,C=0,N=0,O=0,F=0,Ne=0
      integer :: Na=0,Mg=0,Al=0,Si=0,P=0,S=0,Cl=0,Ar=0,K=0,Ca=0
      integer :: Sc=0,Ti=0,V=0,Cr=0,Mn=0,Fe=0,Co=0,Ni=0,Cu=0,Zn=0
      integer :: Ga=0,Ge=0,As=0,Se=0,Br=0,Kr=0,Rb=0,Sr=0
      integer :: Y=0,Zr=0,W=0
      logical :: charge
      character(len=2) :: catm(NELEM)           ! names of elements
      character(len=20),allocatable :: cmol(:)  ! names of molecules
      integer :: elnum(NELEM)                   ! indices of found elements
      integer :: elion(NELEM)                   ! indices of ions
      integer,allocatable :: fit(:)             ! fit-formular identifier
      integer,allocatable :: natom(:)           ! no of atoms in molecule    
      integer,allocatable :: source(:)          ! no of source file
      integer,allocatable :: m_kind(:,:)        ! index of elements
      integer,allocatable :: m_anz(:,:)         ! stoichiometric coeffs
      real*8,allocatable  :: a(:,:)             ! kp fit-coeffs
      real*8,allocatable  :: error(:)           ! kp fit errors
      real*8 :: b_nasa(NELEM,0:13)              ! Added by Yui Kawashima
      integer :: i_nasa,c_nasa(NELEM)           ! Added by Yui Kawashima
      integer :: Natmax,Ncmax,STOImax(NELEM)
      real*8 :: th1,th2,th3,th4,TT1,TT2,TT3     
      end

************************************************************************
      module STRUCTURE
************************************************************************
      use DUST_DATA,ONLY: NELEM
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: Npmax=10000 
      real*8,dimension(Npmax) :: Tgas,press,pelec,dens,nHtot,zz
      real(kind=qp) :: estruc(Npmax,NELEM)
      end

************************************************************************
      module EXCHANGE
************************************************************************
      use CHEMISTRY,ONLY: NMOLE
      use DUST_DATA,ONLY: NELEM
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: nel,nat(NELEM),nion(NELEM)
      real(kind=qp),allocatable :: nmol(:),mmol(:)
      integer :: HII,HeII,CII,NII,OII,NaII,MgII,LiII,ClII
      integer :: AlII,KII,TiII,SII,SiII,FeII,CaII
      integer,parameter :: H=1,He=2,Li=3,Be=4,B=5,C=6,N=7,O=8,F=9
      integer,parameter :: Ne=10,Na=11,Mg=12,Al=13,Si=14,P=15,S=16
      integer,parameter :: Cl=17,Ar=18,K=19,Ca=20,Sc=21,Ti=22
      integer,parameter :: V=23,Cr=24,Mn=25,Fe=26,Co=27,Ni=28
      integer,parameter :: Cu=29,Zn=30,Ga=31,Ge=32,As=33,Se=34
      integer,parameter :: Br=35,Kr=36,Rb=37,Sr=38,Y=39,Zr=40,W=41
      integer*8 :: chemcall=0,chemiter=0,itransform=0,ieqcond=0
      integer*8 :: preIter=0,preEst=0,preUse=0,DUALcorr=0,HCOcorr=0
      integer*8 :: ieqconditer=0
      end
**********************************************************************
      MODULE DATABASE
**********************************************************************
      use dust_data,ONLY: NELEM,NDUSTmax,dust_nam
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,parameter :: DMAX = 2*10**5
      integer :: NDAT=0,NLAST=0,NMODI=0,NPICK1=0,NPICK2=0
      TYPE ENTRY
        real*8 :: ln
        real*8 :: lT
        real*8 :: eprod
        real(kind=qp) :: eps(NELEM)
        real(kind=qp) :: ddust(NDUSTmax)
      END TYPE ENTRY
      TYPE(ENTRY) :: dbase(DMAX)
      end MODULE DATABASE

**********************************************************************
      SUBROUTINE ERASE_DBASE
**********************************************************************
      use DATABASE,ONLY: NDAT,NLAST,NMODI
      implicit none
      NDAT=0
      NLAST=0
      NMODI=0
      end


**********************************************************************
      SUBROUTINE SAVE_DBASE
**********************************************************************
      use dust_data,ONLY: NELEM,NDUST,dust_nam
      use DATABASE,ONLY: NDAT,NLAST,dbase
      implicit none
      integer :: i
      character(len=80) :: filename="database.dat"
      if (NLAST==0) then
        open(unit=11,file=filename,form='unformatted',status='replace')
        write(11) NELEM,NDUST
        write(11) dust_nam
        do i=1,NDAT
          write(11) dbase(i)%ln 
          write(11) dbase(i)%lT
          write(11) dbase(i)%eprod
          write(11) dbase(i)%eps
          write(11) dbase(i)%ddust(1:NDUST)
        enddo 
        close(11)
      else if (NDAT>NLAST) then 
        open(unit=11,file=filename,form='unformatted',position='append')
        do i=NLAST+1,NDAT
          write(11) dbase(i)%ln 
          write(11) dbase(i)%lT
          write(11) dbase(i)%eprod
          write(11) dbase(i)%eps
          write(11) dbase(i)%ddust(1:NDUST)
        enddo 
        close(11)
      endif  
      NLAST = NDAT
      end

**********************************************************************
      SUBROUTINE LOAD_DBASE
**********************************************************************
      use dust_data,ONLY: NELEM,NDUST,dust_nam
      use DATABASE,ONLY: qp,NDAT,NLAST,dbase
      implicit none
      integer :: i,NELEM_read,NDUST_read
      logical :: ex
      character(len=20) :: dust_nam_read(NDUST)
      character(len=80) :: filename="database.dat"

      NDAT = 0
      NLAST = 0
      inquire(file=filename,exist=ex)
      if (.not.ex) goto 200
      open(unit=11,file=filename,form="unformatted",status="old")
      read(11) NELEM_read,NDUST_read
      if (NELEM_read.ne.NELEM) goto 200
      if (NDUST_read.ne.NDUST) goto 200
      read(11) dust_nam_read
      do i=1,NDUST
        if (dust_nam(i).ne.dust_nam_read(i)) goto 200
      enddo
      do i=1,999999
        read(11,end=100) dbase(i)%ln 
        read(11) dbase(i)%lT
        read(11) dbase(i)%eprod
        read(11) dbase(i)%eps
        read(11) dbase(i)%ddust(1:NDUST)
        NDAT = NDAT+1
        !print*,i,EXP(dbase(i)%ln),EXP(dbase(i)%lT)
      enddo 
 100  close(11)
      print*,"... having read ",NDAT," datasets." 
      NLAST = NDAT
      return
 200  close(11)
      print*,"... no / unsuitable database."
      end

**********************************************************************
      SUBROUTINE PUT_DATA(nH,T,eps,ddust,qbest,ibest,active)
**********************************************************************
      use PARAMETERS,ONLY: verbose
      use dust_data,ONLY: NELEM,NDUST,NEPS,dust_nam,eps0,elnr
      use DATABASE,ONLY: qp,NDAT,NLAST,NMODI,DMAX,dbase
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(in) :: nH,T,qbest
      integer,intent(in) :: ibest
      real*8 :: prod
      real(kind=qp),intent(in) :: eps(NELEM),ddust(NDUST)
      logical,intent(in) :: active(0:NDUST)
      integer :: i,j,e,el
      
      !if (qbest<1.d-8) then
      !  return 
      !else 
      if (qbest<1.d-3) then
        i = ibest
        if (verbose>=0) then
          write(*,'(" ... replacing database entry (",I6,") nH,T=",
     >      2(1pE15.7))') i,nH,T
        endif  
      else  
        NDAT = NDAT+1
        i = NDAT
        if (verbose>=0) then
          write(*,'(" ... adding database entry (",I6,") nH,T=",
     >      2(1pE15.7))') i,nH,T
        endif
        if (NDAT>DMAX) then
          print*,"*** NDAT>DMAX in PUT_DATA",NDAT,DMAX
          stop
        endif  
      endif  
      prod = 0.0
      do i=1,NELM
        if (i==iel) cycle
        el = elnum(i)
        prod = prod + LOG(eps0(el))
      enddo  
      dbase(i)%ln    = LOG(nH)
      dbase(i)%lT    = LOG(T)
      dbase(i)%eprod = prod
      dbase(i)%eps   = eps
      do j=1,NDUST
        !if (ddust(j)>0.Q0) print*,active(j),dust_nam(j),real(ddust(j))
        dbase(i)%ddust(j) = ddust(j)
        if (.not.active(j)) dbase(i)%ddust(j)=0.Q0
      enddo
      NMODI = i
      if (NDAT>NLAST+10) then
        call SAVE_DBASE
        print*,"... saved ",NDAT," datasets."
      endif  
      end


**********************************************************************
      subroutine GET_DATA(nH,T,eps,ddust,qbest,ibest,active,method)
**********************************************************************
      use PARAMETERS,ONLY: verbose,remove_condensates
      use dust_data,ONLY: NEPS,NELEM,NDUST,eps0,elnam,elnr,
     >                    dust_nel,dust_nu,dust_el,dust_nam
      use DATABASE,ONLY: qp,NDAT,NMODI,NPICK1,NPICK2,DMAX,dbase
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      implicit none
      real*8,intent(in) :: nH,T
      real*8,intent(out) :: qbest
      integer,intent(out) :: ibest,method
      real(kind=qp),intent(inout) :: eps(NELEM),ddust(NDUST)
      logical,intent(out) :: active(0:NDUST)
      real*8 :: prod,lp,ln,lT,lpread,lnread,lTread
      real*8 :: qual,qmodi,pot,pote,potn,rsort(NEPS)
      real(kind=qp) :: check(NELEM),error,errmax,sjk,sik
      real(kind=qp) :: stoich(NEPS,NEPS),xx(NEPS),rest(NEPS),tmp,val
      real(kind=qp) :: ecopy(NELEM),dcopy(NDUST),deps0(NELEM),echange
      integer :: i,j,k,el,elworst,b,bb,Nbuf,jj,ii,Nact,e
      integer :: isort(NEPS),jmain(NELEM),ibuf(NELEM)
      character(len=1) :: char
      character(len=80) :: frmt
      character(len=999) :: condensates
      logical :: IS_NAN,modified
      logical :: condensed(NELEM),used(NDUST)
      logical,save :: firstCall=.true.
      
      if (firstCall) then
        call LOAD_DBASE 
        firstCall = .false.
      endif  

      prod = 0.0
      do i=1,NELM
        if (i==iel) cycle
        el = elnum(i)
        prod = prod + LOG(eps0(el))
      enddo
      if (verbose>=0) then
        print'("looking for nH,T,eprod=",3(1pE13.5)," ...")',nH,T,prod
      endif  
      ln = LOG(nH)
      lT = LOG(T) 
      lp = prod
      qbest  = 9.d+99
      ibest  = 0
      pot    = -0.03
      pote   = 10.0
      potn   = 0.05
      !--- try last entry modified first ---
      if (NMODI>0) then
        i=NMODI
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        lpread = dbase(i)%eprod
        qual = potn*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + pote*ABS(lpread-lp)
        qbest = qual
        qmodi = qual
        ibest = i
        if (qbest<1.d-3) goto 100
      endif  
      !--- try around entry picked last time ---  
      do i=MAX(1,NPICK1-1),MIN(NDAT,NPICK1+1)
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        lpread = dbase(i)%eprod
        qual = potn*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + pote*ABS(lpread-lp)
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo
      do i=MAX(1,NPICK2-1),MIN(NDAT,NPICK2+1)
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        lpread = dbase(i)%eprod
        qual = potn*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + pote*ABS(lpread-lp)
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo
      if (verbose>=0) write(*,*) "entering full search ..."
      !--- check them all ---  
      do i=NDAT,1,-1
        lnread = dbase(i)%ln 
        lTread = dbase(i)%lT
        lpread = dbase(i)%eprod
        qual = potn*ABS(lnread-ln)+ABS((lTread-lT)+pot*(lnread-ln))
     >       + pote*ABS(lpread-lp)
        if (qual<qbest) then 
          qbest = qual
          ibest = i
          if (qbest<1.d-3) goto 100
        endif  
      enddo

 100  active = .false.
      condensates = ""
      if (ibest>0) then
        eps    = dbase(ibest)%eps
        ddust  = dbase(ibest)%ddust(1:NDUST)
        Nact = 0
        do i=1,NDUST
          if (ddust(i)>0.Q0) then
            Nact = Nact+1
            active(i)=.true.
            condensates = trim(condensates)//" "//trim(dust_nam(i))
          endif  
        enddo
        NPICK2 = NPICK1
        NPICK1 = ibest
        if (verbose>=0) then
          write(*,'(" ... found best dataset (",I6,")  
     >      nH,T,eprod,qual=",4(1pE13.5))')
     >      ibest,EXP(dbase(ibest)%ln),EXP(dbase(ibest)%lT),
     >      dbase(ibest)%eprod,qbest
          write(*,*) "active condensates: "//trim(condensates)
        endif  

        condensed(:) = .false.
        do i=1,NDUST
          do j=1,dust_nel(i)
            if (ddust(i)==0.Q0) cycle
            el = dust_el(i,j)
            condensed(el) = .true.
          enddo
        enddo
        modified = .false.
        do i=1,NELM
          if (i==iel) cycle
          el = elnum(i)
          if (.not.condensed(el)) then
            !print*,elnam(el)," has no condensates"
            if (ABS(eps(el)/eps0(el)-1.Q0).gt.1.E-4) modified=.true.
            eps(el) = eps0(el)
          endif
        enddo
        if (modified) qbest=MIN(qbest,0.2)
        ecopy = eps
        dcopy = ddust

        if (remove_condensates) then
          !--- check whether ddust=0 would result in perfect fit ---
          errmax = -1.Q0
          do e=1,NEPS
            el = elnr(e) 
            error = ABS(1.Q0-eps(el)/eps0(el))
            if (error.gt.errmax) then
              errmax = error
              elworst = el
            endif   
          enddo
          if (errmax<1.Q-25) then
            if (verbose>-1) print*,"return with ddust=0",errmax
            ddust = 0.Q0
            qbest = MIN(qbest,0.01)
            return
          endif
        endif
  
        check = eps
        do i=1,NDUST
          do j=1,dust_nel(i)
            if (ddust(i)==0.Q0) cycle
            el = dust_el(i,j)
            check(el) = check(el) + ddust(i)*dust_nu(i,j)    
          enddo
        enddo
        deps0 = eps0-check
        errmax = -1.Q0
        do e=1,NEPS
          el = elnr(e) 
          error = ABS(1.Q0-check(el)/eps0(el))
          if (error.gt.errmax) then
            errmax = error
            elworst = el
          endif   
        enddo
        if (verbose>=0) then
          print*,"need fitting? "//elnam(elworst),SNGL(errmax)
        endif  
        method = 0
        if (errmax<1.Q-25) return ! perfect fit - nothing to do

        !--------------------------------------------------------------
        ! ***  adapt eps and ddust to reach current eps0: method 1  ***
        !--------------------------------------------------------------
        if (verbose>=0) print*,"adjust with method 1 ..."
        method = 1

        !--- 1. sort elements from rare to abundant ---
        rsort = 9.d+99
        isort = 0
        do i=1,NEPS
          el = elnr(i)
          val = eps(el)
          if (.not.condensed(el)) val=val+1.E+10  ! non-condensed to bottom
          do j=NEPS+1,2,-1
            if (val>rsort(j-1)) exit
          enddo  
          isort(j+1:NEPS) = isort(j:NEPS-1)
          rsort(j+1:NEPS) = rsort(j:NEPS-1)
          isort(j) = el
          rsort(j) = val
        enddo
        jmain = 0                                 ! list of condensates
        used  = .false.
        Nbuf  = 0
        do j=1,NDUST                              
          if (ddust(j)==0.Q0) cycle
          Nbuf = Nbuf+1          
          el = isort(Nbuf)
          ibuf(Nbuf) = el
          jmain(el) = j
          used(j)=.true.
          !print*,elnam(el),dust_nam(j),REAL(eps(el))
        enddo  

        !--------------------------------------------------------------
        !*** 2. setup linear equation system                        ***
        !*** the idea is to find modified ddust, such that eps is   ***
        !*** as in the database:  ((stoich))*(ddust) = (eps0-eps)   ***
        !*** ibuf(1...Nbuf) : sorted list of elements               ***
        !*** jmain(el) : dust index of main reservoir               ***
        !-------------------------------------------------------------- 
        stoich = 0.Q0
        do b=1,Nbuf
          el = ibuf(b) 
          do bb=1,Nbuf
            j = jmain(ibuf(bb))
            do k=1,dust_nel(j)
              if (dust_el(j,k)==el) then
                 stoich(b,bb) = dust_nu(j,k) 
              endif
            enddo
          enddo
        enddo  
        check = eps0 - eps
        do j=1,NDUST
          if (ddust(j)==0.Q0) cycle
          if (used(j)) cycle 
          do k=1,dust_nel(j)
            el = dust_el(j,k) 
            check(el) = check(el) - ddust(j)*dust_nu(j,k)
          enddo
        enddo  
        do b=1,Nbuf
          el = ibuf(b) 
          j = jmain(el)
          rest(b) = check(el)
          if (verbose>=0) then
            write(frmt,'("(A2,2x,",I2,"(I2),A16,2(1pE13.6))")') Nbuf
            write(*,frmt) elnam(el),INT(stoich(b,1:Nbuf)),
     &           trim(dust_nam(j)),ddust(j),rest(b)
          endif  
        enddo  
        call GAUSS16( NEPS, Nbuf, stoich, xx, rest)
        do b=1,Nbuf
          el = ibuf(b)
          j = jmain(el)
          if (verbose>=0) then
            print'(A3,A16,1pE13.6," ->",1pE13.6)',
     &         elnam(el),trim(dust_nam(j)),ddust(j),xx(b)
          endif  
          ddust(j) = xx(b)
          if (xx(b)<0.Q0.or.IS_NAN(DBLE(xx(b)))) then
            if (verbose>=0) then
              print*,"*** negative dust abund. in database.f method 1"
            endif  
            goto 300
            !qbest = 9.d+99
            !return
          endif  
        enddo 
 
        !--- 3. correct remaining gas element abundances ---
        check = 0.Q0
        do i=1,NDUST
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            check(el) = check(el) + ddust(i)*dust_nu(i,j)    
          enddo
        enddo
        echange = 0.Q0
        do i=1,NEPS
          el = isort(i)
          if (jmain(el)==0) then
            tmp = eps0(el)-check(el)
            if (verbose>=0) then
              print'(A3,A16,1pE13.6," ->",1pE13.6)',
     &           elnam(el),"gas",eps(el),tmp
            endif  
            if (tmp<0.Q0.or.IS_NAN(DBLE(tmp))) then
              if (verbose>=0) then
                print*,"*** negative el. abund. in database.f method 1"
              endif  
              goto 300
              !qbest = 9.d+99
              !return
            endif  
            eps(el) = tmp
            echange = MAX(echange,ABS(LOG10(eps(el)/ecopy(el))))
          endif
        enddo
        if (verbose>=0) print*,'echange =',REAL(echange)
        if (echange<2.Q0) goto 500

 300    continue    
        !--------------------------------------------------------------
        ! ***  adapt eps and ddust to reach current eps0: method 2  ***
        !--------------------------------------------------------------
        if (verbose>=0) print*,"adjust with method 2 ..."
        method = 2
        !--------------------------------------------------------------
        ! ***  the idea is to minimise the change of eps            ***
        ! ***      Sum_k (deps_k/eps_k)**2 -> Min                   ***
        ! ***  where deps_k = eps_k - eps_k(table)                  ***
        ! ***  with eps_k        = eps0_k        - Sum_j sjk xj     ***
        ! ***  and  eps_k(table) = eps0_k(table) - Sum_j sjk xj(table)
        ! ***       ---------------------------------------------------
        ! ***       deps_k       = deps0_k       - Sum_j sjk xx     ***
        ! ***       ---------------------------------------------------
        ! ***  the differentiate with respect to each xx_i = 0      ***
        ! -------------------------------------------------------------
        ddust = dcopy
        eps   = ecopy
        Nact = 0
        stoich = 0.Q0
        rest = 0.Q0
        i = 0
        do ii=1,NDUST
          if (ddust(ii)==0.Q0) cycle
          Nact = Nact+1
          i = i+1
          do k=1,NEPS
            el = elnr(k)
            do e=1,dust_nel(ii)
              if (el==dust_el(ii,e)) exit
            enddo  
            sik = dust_nu(ii,e)
            !print*,dust_nam(ii),elnam(el),sik,deps0(el)
            rest(i) = rest(i) + sik*deps0(el)/eps(el)**2
            j = 0
            do jj=1,NDUST
              if (ddust(jj)==0.Q0) cycle
              j = j+1
              do e=1,dust_nel(jj)
                if (el==dust_el(jj,e)) exit
              enddo  
              sjk = dust_nu(jj,e)
              stoich(i,j) = stoich(i,j) + sik*sjk/eps(el)**2 
            enddo  
          enddo  
        enddo
        !print*,NEPS,Nact,rest
        call GAUSS16( NEPS, Nact, stoich, xx, rest)
        eps = eps + deps0
        j = 0
        do jj=1,NDUST
          if (ddust(jj)==0.Q0) cycle
          j = j+1
          if (verbose>=0) then
            print'(A12,2(1pE12.4))',dust_nam(jj),dcopy(jj),xx(j)
          endif  
          ddust(jj) = ddust(jj)+xx(j)
          do e=1,dust_nel(jj)
            el = dust_el(jj,e)
            eps(el) = eps(el)-dust_nu(jj,e)*xx(j)
          enddo
          if (ddust(jj)<0.Q0.or.IS_NAN(DBLE(ddust(jj)))) then
            method = 0
            if (verbose>=0) then
              print*,"*** negative dust abund. in database.f method 2"
            endif  
            qbest = 9.d+99
            return
          endif
        enddo 
        do e=1,NEPS
          el = elnr(e)
          if (verbose>=0) then
            print'(A3,2(1pE12.4))',elnam(el),ecopy(el),eps(el)-ecopy(el)
          endif  
          if (eps(el)<0.Q0.or.IS_NAN(DBLE(eps(el)))) then
            method = 0
            if (verbose>=0) then
              print*,"*** negative el. abund. in database.f method 2"
            endif  
            qbest = 9.d+99
            return
          endif
        enddo

 500    continue
        check = eps
        do i=1,NDUST
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            check(el) = check(el) + ddust(i)*dust_nu(i,j)    
          enddo
        enddo
        errmax = -1.Q0
        do i=1,NELM
          if (i==iel) cycle
          el = elnum(i)
          !print*,elnam(el),jmain(el),check(el),eps0(el)
          error = ABS(1.Q0-check(el)/eps0(el))
          if (error.gt.errmax) then
            errmax = error
            elworst = el
          endif   
        enddo  
        !print*,"check ",elnam(elworst),errmax
        if (errmax>1.Q-8) then
          print*,"*** element conservation violation in database.f"
          print*,elnam(elworst),errmax
          stop
        endif  
        !--- for remove_condensates models: force new entry --- 
        if (maxval(ddust)==0.Q0) qbest=MIN(qbest,0.01)
      endif

      end
*********************************************************************
      SUBROUTINE SUPERSAT(T,nat,nmol,Sat)
*********************************************************************
      use CHEMISTRY,ONLY: NMOLE,cmol
      use DUST_DATA,ONLY: NELEM,NDUST,bk,atm,rgas,bar,fit,cfit,
     &                    dust_nam,dust_nel,dust_el,dust_nu,dust_mass,
     &                    is_liquid,Tcorr,elnam,Nfit,Tfit,Bfit
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: T
      real(kind=qp),intent(in) :: nat(NELEM),nmol(NMOLE)
      real(kind=qp),intent(out):: Sat(NDUST)
      real(kind=qp),parameter :: cal=4.184Q+0       ! 1 cal in J
      real(kind=qp),parameter :: mmHg=1.3328Q+3     ! 1 mmHg in dyn/cm2
      real(kind=qp),parameter :: Joule=1.Q+7        ! 1 Joule in erg
      real(kind=qp),parameter :: mol=6.02214076Q+23 ! 1 mol
      real(kind=qp),parameter :: eV=1.60218Q-12     ! 1 eV in erg

      real(kind=qp) :: T1,T2,T3,TC,kT,RT,dG,lbruch,pst,psat,dGRT
      real(kind=qp) :: a(0:4),term,n1,natom,aa(0:6)
      integer :: i,j,l,STINDEX,el,imol,imol1,imol2,ifit
      character(len=20) :: search,upper,leer='                    '


      T1  = T
      T2  = T1**2
      T3  = T1**3
      kT  = bk*T1
      RT  = rgas*T1
      Sat = 0.Q0
      do i=1,NDUST
        a(:) = cfit(i,:)
        if (fit(i)==1) then
          !-------------------------------------
          !***  dG-fit Sharp & Huebner 1990  ***
          !-------------------------------------
          pst = atm
          dG = a(0)/T1 + a(1) + a(2)*T1 + a(3)*T2 + a(4)*T3
          dG = dG/(rgas/cal*T1)
          lbruch = 0.Q0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            term   = nat(el)*kT/pst
            lbruch = lbruch + dust_nu(i,j)*LOG(term)
          enddo
          Sat(i) = EXP(lbruch-dG)

        else if (fit(i)==2) then
          !-----------------------------------
          !***  dG polynom-fit NIST-Janaf  ***
          !-----------------------------------
          pst = bar
          dG = a(0)/T1 + a(1) + a(2)*T1 + a(3)*T2 + a(4)*T3
          dGRT = dG/RT
          lbruch = 0.Q0
          Natom = 0.0          
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            term   = nat(el)*kT/pst
            lbruch = lbruch + dust_nu(i,j)*LOG(term)
            Natom  = Natom + dust_nu(i,j)            
          enddo
          Sat(i) = EXP(lbruch-dGRT)
          !print*,dust_nam(i),T1,-dG*Joule/mol/eV/Natom,      ! eV/atom
     >    !                      -dG*Joule/mol/dust_mass(i)   ! erg/g          

        else if (fit(i)==3) then
          !-----------------------------------------
          !***  ln(pvap) polynom-fit NIST-Janaf  ***
          !-----------------------------------------
          psat = EXP(a(0)/T1 + a(1) + a(2)*T1 + a(3)*T2 + a(4)*T3)
          if (dust_nel(i)==1.and.dust_nu(i,1)==1) then
            n1 = nat(dust_el(i,1))
          else
            search = trim(dust_nam(i))
            l = index(search,'[')
            search = upper(search(1:l-1)//leer(l:20))
            imol = STINDEX(cmol,NMOLE,search)
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            n1 = nmol(imol)
          endif
          Sat(i) = n1*kT/psat

        else if (fit(i)==4) then
          !----------------------------------------
          !***  ln(pvap) 3-para-fit NIST-Janaf  ***
          !----------------------------------------
          psat = EXP(a(0) + a(1)/(T1 + a(2)))
          if (dust_nel(i)==1.and.dust_nu(i,1)==1) then
            n1 = nat(dust_el(i,1))
          else
            search = trim(dust_nam(i))
            l = index(search,'[')
            search = upper(search(1:l-1)//leer(l:20))
            imol = STINDEX(cmol,NMOLE,search)
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            n1 = nmol(imol)
          endif
          Sat(i) = n1*kT/psat

        else if (fit(i)==5) then
          !--------------------------
          !***  -dG/RT Stock-fit  ***
          !--------------------------
          pst = bar
          dGRT = a(0)/T1 + a(1)*LOG(T1) + a(2) + a(3)*T1 + a(4)*T2
          lbruch = 0.Q0
          Natom = 0.0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            term   = nat(el)*kT/pst
            lbruch = lbruch + dust_nu(i,j)*LOG(term)
            Natom  = Natom + dust_nu(i,j)
          enddo
          Sat(i) = EXP(lbruch+dGRT)
          !dG = dGRT*RT                                      ! Joule/mol
          !print*,dust_nam(i),T1,dG*Joule/mol/eV/Natom,      ! eV/atom
     >    !                      dG*Joule/mol/dust_mass(i)   ! erg/g

        else if (fit(i)==6) then
          !---------------------------------------------------------------
          !***  Yaws' Chemical Properties Handbook (McGraw-Hill 1999)  ***
          !---------------------------------------------------------------
          psat = a(0) + a(1)/T1 + a(2)*LOG10(T1) + a(3)*T1 + a(4)*T2
          psat = 10.Q0**psat * mmHg
          if (dust_nel(i)==1.and.dust_nu(i,1)==1) then
            n1 = nat(dust_el(i,1))
          else
            search = trim(dust_nam(i))
            l = index(search,'[')
            search = upper(search(1:l-1)//leer(l:20))
            imol = STINDEX(cmol,NMOLE,search)
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            n1 = nmol(imol)
          endif
          Sat(i) = n1*kT/psat

        else if (fit(i)==7) then
          !---------------------
          !***  BURCAT data  ***
          !---------------------
          do ifit=1,Nfit(i)
            if (T1<Tfit(i,ifit+1)) exit
            if (ifit==Nfit(i)) exit
          enddo
          if (T1>1000.0) then
            aa(0:6) = Bfit(i,ifit,1:7)
          else
            aa(0:6) = Bfit(i,ifit,8:14)
          endif  
          dGRT = aa(0) *(log(T1)-1.d0) 
     &         + aa(1) *T1   / 2.d0
     &         + aa(2) *T1**2/ 6.d0    
     &         + aa(3) *T1**3/12.d0
     &         + aa(4) *T1**4/20.d0    
     &         - aa(5) /T1  
     &         + aa(6)
          lbruch = 0.Q0
          Natom = 0.0
          do j=1,dust_nel(i)
            el     = dust_el(i,j)
            term   = nat(el)*kT/pst
            lbruch = lbruch + dust_nu(i,j)*LOG(term)
            Natom  = Natom + dust_nu(i,j)
          enddo
          Sat(i) = EXP(lbruch+dGRT)
          !print*,dust_nam(i),"BURCAT:",T1,ifit,dGRT,Sat(i)
          if (T1>1.2*Tfit(i,Nfit(i)+1)) then
            Sat(i) = 1.Q-999   ! BURCAT-fits go crasy beyond fit-regime
          endif  
          
        else if (fit(i)==99) then
          !-----------------------
          !***  special cases  ***
          !-----------------------
          if (dust_nam(i).eq.'H2O[l]') then
            !--- Ackerman & Marley 2001 ---
            TC   = MIN(2000.0,T1)-273.15         ! T[degree Celsius]
            psat = 6112.1*exp((18.729*TC - TC**2/227.3)/(TC + 257.87))
            imol = STINDEX(cmol,NMOLE,"H2O")
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'H2O[s]') then
            !--- Ackerman & Marley 2001 ---
            TC   = MIN(2000.0,T1)-273.15         ! T[degree Celsius]
            psat = 6111.5*exp((23.036*TC - TC**2/333.7)/(TC + 279.82))
            imol = STINDEX(cmol,NMOLE,"H2O")
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'NH3[s]') then
            !--- CRC Handbook of Chemistry and Physics (Weast 1971) ---
            psat = exp(10.53 - 2161.Q0/T1 - 86596.Q0/T2)*bar
            imol = STINDEX(cmol,NMOLE,"NH3")
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'CH4[s]') then
            !--- Prydz, R.; Goodwin, R.D., J. Chem. Thermodyn., 1972, 4,1 ---
            psat = 10.0**(3.9895 - 443.028/(T1-0.49))*bar
            imol = STINDEX(cmol,NMOLE,"CH4")
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'NH4SH[s]') then
            !--- G.Lee's fit to Walker & Lumsden (1897) ---
            psat  = 10.0**(7.8974 - 2409.4/T1)*bar
            imol1 = STINDEX(cmol,NMOLE,"NH3")
            imol2 = STINDEX(cmol,NMOLE,"H2S")
            if (imol1<=0.or.imol2<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            Sat(i) = SQRT(nmol(imol1)*kT/(psat*0.5))
     &             * SQRT(nmol(imol2)*kT/(psat*0.5))

          else if (dust_nam(i).eq.'H2S[s]') then
            !--- Stull (1947) ---
            if (T1 < 30.0) then ! Limiter for very cold T
              T1 = 30.0
            end if
            if (T1 < 212.8) then
              psat = 10.0**(4.43681 - 829.439/(T1-25.412))*bar
            else
              psat = 10.0**(4.52887 - 958.587/(T1-0.539))*bar
            end if
            imol = STINDEX(cmol,NMOLE,"H2S")
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'S2[s]') then
            !--- Zahnle et al. (2016) ---
            if (T1 < 413.0) then
              psat = exp(27.0 - 18500.0/T1)*bar
            else
              psat = exp(16.1 - 14000.0/T1)*bar
            end if
            !--- Lyons 2008 ---
            !write(50,'(F8.1,2(1pE13.4))')
     &      !         T1,psat,10.0**(7.024 - 6091.0/T1)*bar
            !psat = 10.0**(7.024 - 6091.0/T1)*bar
            imol = STINDEX(cmol,NMOLE,"S2")
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            Sat(i) = nmol(imol)*kT/psat

          else if (dust_nam(i).eq.'S8[s]') then
            !--- Zahnle et al. (2016) ---
            if (T1 < 413.0) then
              psat = exp(20.0 - 11800.0/T1)*bar
            else
              psat = exp(9.6 - 7510.0/T1)*bar
            end if
            !--- Lyons 2008 ---
            !psat = 10.0**(4.188 - 3269.0/T1)*bar
            imol = STINDEX(cmol,NMOLE,"S8")
            if (imol<=0) then
              print*,"*** supersat.f molecule not found ",dust_nam(i)
              stop
            endif
            Sat(i) = nmol(imol)*kT/psat

          else
            print*,"*** supersat.f fit=",fit(i)," ??? ",dust_nam(i)
            stop
          endif

        else
          print*,"*** supersat.f fit=",fit(i)," ??? ",dust_nam(i)
          stop
        endif

        if (Tcorr(i)>0.0) then
          if (is_liquid(i).and.T1<Tcorr(i)) then
            Sat(i) = Sat(i)/EXP(0.1*(Tcorr(i)-T1))
          endif
          if ((.not.is_liquid(i)).and.T1>Tcorr(i)) then
            Sat(i) = Sat(i)/EXP(0.1*(T1-Tcorr(i)))
          endif
        endif

      enddo

      RETURN
      end

************************************************************************
      integer function stindex(array,dim,string)
************************************************************************
      integer dim
      character*(*) array(dim),string
      integer i
      
      do i=1,dim
        if ( array(i) .eq. string ) then
          stindex = i
          return
        endif
      enddo
      write(6,*) 'not found: ' ,string
      stindex = 0
      end
************************************************************************
      SUBROUTINE smchem16 (anHges,Tg,eps,anmono,anmol,merk,verbose)
************************************************************************
*                                                                      *
*     small chemistry                                                  *
*     ---------------                                                  *
*     Diese Routine berechnet eine GG-Chemie                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     e i n g a b e :                                                  *
*     anHges : totale Dichte der Wasserstoffatome (~ nh+2*nh2)         *
*     tg     : Gastemperatur                                           *
*     eps    : Vektor mit linearen chemischen Haeufigkeiten ([H]=1)    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     a u s g a b e :                                                  *
*     anmono : vektor mit den dichten der monomere  falls relevant     *
*     anmol  : vektor mit den dichten der molekuele falls relevant     *
*                                                                      *
************************************************************************
      use CHEMISTRY,ONLY: NewBackIt,NewFullIt,NewBackFac,
     >                    NewFastLevel,NewPreMethod,
     >                    nml=>NMOLE,nel=>NELM,cmol,catm,
     >                    m_kind,m_anz,charge,elion,el,
     >                    Natom,Ncmax,STOImax,H,C,N,O,
     >                    th1,th2,th3,th4,TT1,TT2,TT3
      use EXCHANGE,ONLY: chemcall,chemiter,preUse,preIter,preEst,
     >                   DUALcorr,HCOcorr
      implicit none
*-----------------------------------------------------------------------
*  Dimensionierung fuer die Molekuel- und Atom Felder. Hier ist auf
*  Konsistenz mit dem aufrufenden Programm zu achten.
*-----------------------------------------------------------------------
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in)  :: anHges,Tg
      real(kind=qp),intent(in)  :: eps(nel)
      real(kind=qp),intent(out) :: anmono(nel),anmol(nml)
      integer,intent(inout) :: verbose
      real(kind=qp),parameter :: bk=1.380662Q-16
*-----------------------------------------------------------------------
*  Bei merk=.true. merkt sich die Routine die letzte konvergiert Loesung
*  und geht beim naechsten Mal von diesen Startwerten aus.
      logical,intent(INOUT) :: merk
*-----------------------------------------------------------------------
*  Die Variable "ngestst" entscheidet, ob die Elementerhaltung ueber-
*  prueft werden soll.
      logical,parameter :: ngestst=.false.
*-----------------------------------------------------------------------
      integer :: stindex,Nconv,switch,ido,iredo,nu,tmp
      integer :: Nact,all_to_act(nel),act_to_all(nel),switchoff(nel)
      integer :: e,i,j,j1,ii,jj,l,it,m1,m2,piter,ifatal,ipull,pullmax
      integer :: lmin,lmax,Nseq,iloop,imethod,enew,i1,i2,i3,s1,s2,s3,s4
      integer :: eseq(nel),imaj(nel),imaj2(nel)
      integer :: NpreLoop,NpreIt,Ntaken,Nestim,DUALco,HCOco
      integer,parameter :: itmax=200
      real(kind=qp) :: finish,qual,qual0,qual1,qual2,qual3,qmost,qq
      real(kind=qp) :: g(0:nml),limit
      real(kind=qp) :: pH2O,pCO2,pCH4,lpH2O,lpCO2,lpCH4,lpH,lpO,lpC
      real(kind=qp) :: kT,kT1,cc,nelek,ng,Sa,fak,lth,arg,term,f,fs
      real(kind=qp) :: pel,delta,pat,atfrac,atmax,lnp,lnp3
      real(kind=qp) :: nges(nel),coeff(-3:Ncmax),lnc(-3:Ncmax)
      real(kind=qp) :: DF(nel,nel),dp(nel),FF(nel),pmol,delp,nold
      real(kind=qp) :: DF0(nel,nel),FF0(nel),xscale(nel),fscale(nel)
      real(kind=qp) :: nsave(nel),null(nel),sortq(nel)
      real(kind=qp) :: conv(0:itmax,nel),converge(-1:itmax)
      real(kind=qp) :: soll,haben,abw,sum,maxs,mcharge,crit,crit1,crit2
      real(kind=qp) :: pbefore(nel),norm(nel),xx(nel)
      real(kind=qp) :: emax,pges,pwork,pwork2,pp,psc,ptest,aim
      logical :: from_merk,eact(nel),redo(nel),done(nel),affect,known
      logical :: relevant(nml)
      logical :: possible,ptake,HCOproblem,hasH,hasC,hasO,IS_NAN
      character(len=5000) :: mols
      character(len=100) :: txt
      character(len=1) :: char,bem
      integer,save :: TiC,H2O,CH4,CO2,ilauf=0
      real(kind=qp),allocatable,save :: amerk(:),ansave(:)
      real(kind=qp),allocatable,save :: badness(:),pcorr(:,:) 
      integer,allocatable,save :: pkey(:)
!$omp threadprivate(TiC,ilauf,amerk,ansave,badness,pcorr,pkey)
*-----------------------------------------------------------------------      

      ifatal = 0
      if (.not.allocated(amerk)) then
        allocate(badness(nel),pcorr(nel,nel),pkey(nel),
     >           amerk(nel),ansave(nel))
        TiC = stindex(cmol,nml,'TIC    ')
        H2O = stindex(cmol,nml,'H2O    ')
        CH4 = stindex(cmol,nml,'CH4    ')
        CO2 = stindex(cmol,nml,'CO2    ')
        badness = 1.Q0
        pcorr   = 1.Q0
        pkey    = 0
      endif

*-----------------------------------------------------------------------
*     ! zu niedrige Temperaturen abfangen und
*     ! Variable fuer die Dissoziationskonstanten berechnen
*     =====================================================
      TT1 = Tg
      TT2 = TT1*TT1
      TT3 = TT2*TT1
      th1 = 5040.Q0/TT1
      th2 = th1*th1
      th3 = th2*th1
      th4 = th3*th1
      kT  = bk*TT1
      kT1 = 1.Q0/kT
      
*-----------------------------------------------------------------------
*     ! init vectors
*     ==============
      anmono = 0.Q0
      anmol  = 0.Q0

* --------------------------------------------------------------------------
*     ! compute equilibrium constants
*     ===============================
      do i=1,nml
        if (i.ne.TiC) g(i)=gk(i)       ! compute all equil.constants lnk
      enddo  

*    TiC Gleichgewichtskonstante von Andreas Gauger ist anders
*        definiert als die Gleichgewichtskonstanten von Gail
*  Gauger: 
*  log Kp = 12.75293-5.4485*th1-1.56672*log(th1)+1.56041*(log(th1))**2
*           - 0.93275(log(th1))**3
*         = log ( p(A)p(B)/p(AB) )
*  Gail:
*   ln Kp = a0 + a1*th1 + a2*th1**2 + a3*th1**3 + a4*th1**4
*         =  ln ( p(AB)/p(A)p(B) )
*  Tsuji:
*  log Kp = a0 + a1*th1 + a2*th1**2 + a3*th1**3 + a4*th1**4
*         = log ( p(A)p(B)/p(AB) )
*
*  Umrechnung der Gauger-TiC-GG-Konstante in das Gail'sche System
*  -log(10)*log Kp(Gauger) = -2.30256*log Kp(Gauger) = ln Kp(Gail)

      lth = LOG10(th1)
      arg = 12.75293 - 5.44850*th1    - 1.56672*lth
     &               + 1.56041*lth**2 - 0.93275*lth**3
      g(TiC) = -2.30256*arg

*---------------------------------------------------------------------------
      !print'("smchem16 called ilauf,merk",i5,l2)',ilauf,merk
      !print'("T,n<H>,eps=",0pF10.2,99(1pE11.3))',Tg,anHges,eps
      NpreLoop = 0
      NpreIt = 0
      Ntaken = 0
      Nestim = 0
      DUALco = 0
      HCOco  = 0
      if ((ilauf.gt.10).and.merk) then
        do i=1,nel
          anmono(i) = amerk(i) * anhges
        enddo
        from_merk = .true.
        goto 200
      endif

*---------------------------------------------------------------------------
*     ! estimate electron density
*     =========================== 
 100  continue
      if (verbose>1) print'("SMCHEM_16:")'
      from_merk = .false.
      if (charge) then
        nelek = 0.Q0 
        do i=1,nel
          if (i==el) cycle 
          ng = anHges * eps(i)
          Sa = EXP(g(elion(i)))*kT1
          nelek = nelek + ng/(0.5Q0 + SQRT(0.25Q0 + ng/Sa))
        enddo
        anmono(el) = nelek
        pel = nelek*kT
        mcharge = pel
        !peest = pel
        !pel = pecorr*pel
        !anmono(el) = pecorr*anmono(el) 
        if (verbose>1) print'(" estimate pel=",1pE10.3)',pel
      endif  

*-----------------------------------------------------------------------
*     ! estimate atomic pressures: new method
*     =======================================
      sortq(:) = 1.Q0
 110  Nseq = nel
      done(:) = .false.                ! all elements to be estimated here
      !if (charge) done(el)=.true.     ! ... except for the electrons      
      !if (charge) Nseq=nel-1
      eseq(:) = 0                      ! hirachical sequence of elements
      ptake = .true.
      do ido=1,Nseq
        !---------------------------------------------------------
        ! search for the most abundant element not yet considered 
        !---------------------------------------------------------
        emax = 0.Q0
        enew = 0
        do e=1,nel
          if (done(e)) cycle
          norm(e) = eps(e)*sortq(e)
          if (e==el) norm(e)=anmono(el)*sortq(el)/anHges
          if (norm(e)<emax.or.(ido==1.and.e==el)) cycle
          emax = norm(e)
          enew = e
        enddo
        if (enew==0) then
          print*,catm 
          print*,eps
          print*,ido
          print*,done
          print*, "*** should not occur."
          return
        endif  
        if (verbose>1) print*,'estimate p'//trim(catm(enew))//' ...'
        if (enew.ne.pkey(ido)) ptake=.false.
        done(enew) = .true.
        eseq(ido) = enew               ! add to hirachical sequence 
        pges = eps(enew)*anHges*kT
        if (verbose>0) mols = ''
        !--------------------------------------------------------------------
        ! make rough estimate of atom pressure, considering that one molecule 
        ! takes all the element.  Store coeff for Sum_l coeff(l) p^l = pges 
        !--------------------------------------------------------------------
        coeff(:) = 0.Q0   
        lnc(:) = 0.Q0
        pwork = pges
        pwork2 = 0.Q0
        imaj(enew) = 0
        imaj2(enew) = 0
        lmin = +99
        lmax = -99
        do e=1,nel
          if (.not.done(e)) cycle
          if (e==enew) cycle
          xx(e) = LOG(anmono(e)*kT)
        enddo  
        do i=1,nml
          affect = .false.
          known  = .true.
          lnp = g(i)
          do j=1,m_kind(0,i)
            e  = m_kind(j,i)
            nu = m_anz(j,i)
            if (.not.done(e)) then
              known = .false.
              exit
            endif  
            if (e==enew) then
              l = nu
              affect = .true.
            else
              lnp = lnp + nu*xx(e)
            endif  
          enddo  
          if (.not.affect) cycle  
          if (.not.known) cycle
          if (verbose>0) mols = trim(mols)//" "//cmol(i)
          lmin = MIN(lmin,l)
          lmax = MAX(lmax,l)
          coeff(l) = coeff(l) + l*EXP(lnp)
          if (l>0.and.(enew/=el)) then
            lnc(l) = MAX(lnc(l),LOG(REAL(l,kind=qp))+lnp)
            sum = LOG(pges) - lnp - LOG(REAL(l,kind=qp))
            if (sum>300.q0*l) cycle
            ptest = EXP(sum/l)
            if (ptest<pwork) then
              pwork2 = pwork
              pwork = ptest
              imaj2(enew) = imaj(enew)
              imaj(enew) = i
              !print*,"first",cmol(i),pwork,lnc(l)
            else if (ptest<pwork2) then  
              pwork2 = ptest
              imaj2(enew) = i
              !print*,"second",cmol(i),pwork,lnc(l)
            endif
          endif
        enddo
        if (verbose>1) print*,trim(mols)
        if (enew==el) then
          !print*,coeff(-1),coeff(1)
          pel = SQRT(-coeff(-1)/(1.Q0+coeff(+1)))    ! 0 = pel - a/pel + b*pel
          if (verbose>1) print*,'pel=',REAL(anmono(el)*kT),REAL(pel)
          anmono(el) = pel*kT1
          mcharge = MAX(pel,MAX(coeff(-1)/pel,coeff(1)*pel))
        else   
          !----------------------------------------
          ! scale to make coeff(:) fit into real*16  
          !----------------------------------------
          aim = LOG(1.Q+10)
          psc = -LOG(pges)                           ! when atom dominates
          do l=lmin,lmax
            if (lnc(l)==0.Q0) cycle
            ptest = (aim-lnc(l))/l                   ! lnc(l)+l*ln(psc) = aim
            psc = MAX(psc,-ptest)
          enddo  
          psc = -psc
          !if (verbose>1) print*,"pwork=",REAL(pwork),"  psc=",REAL(psc)
          !------------------------------------------------
          ! store coeff for Sum_l coeff(l) (p*psc)^l = pges 
          !------------------------------------------------
          coeff(:) = 0.Q0
          do l=lmin,lmax
            if (lnc(l)==0.Q0) cycle
            coeff(l) = EXP(lnc(l)+l*psc)
          enddo
          psc = EXP(psc)
          !if (verbose>1) print'(" coeff(",I2,":",I2,")=",99(1pE10.2))',
     >    !               lmin,lmax,coeff(lmin:lmax)
          !----------------------------------------------
          ! solve 1d equation above with Newton's method 
          !----------------------------------------------
          pp = pwork/psc
          do piter=1,99                  
            f  = pp*psc-pges
            fs = psc
            do l=lmin,lmax
              if (coeff(l)==0.Q0) cycle
              f  = f  + coeff(l)*pp**l
              fs = fs + coeff(l)*l*pp**(l-1)
            enddo
            if (fs==0.Q0) print*, "*** fs=0 in smchem16 1d-pre-it."
            delta = f/fs
            pp = pp-delta
            if (verbose>1) print'(A2,I3,1pE25.15,1pE10.2)',
     >                     catm(enew),piter,pp*psc,delta/pp
            if (ABS(delta)<1.Q-4*ABS(pp)) exit 
          enddo  
          if (piter>=99) then
            write(*,*) "*** smchem16 no conv in 1D pre-it "//catm(enew)
            write(*,*) anHges,Tg
            do i=1,nel
              write(*,'(A2,1x,0pF30.26)') catm(i),eps(i)
            enddo
            write(*,*) "psc=",real(psc),lmin,lmax
            write(*,*) "  lnc:",lnc(lmin:lmax)
            write(*,*) "coeff:",coeff(lmin:lmax)
            if (ifatal<=1) then
              from_merk = .false.
              badness   = 1.Q0
              pcorr     = 1.Q0
              ifatal    = ifatal+1
              goto 100        ! try again from scratch before giving up
            endif  
            goto 1000
          endif  
          anmono(enew) = pp*psc*kT1
        endif  

        !-----------------------------------------------------------
        ! take into account feedback on elements considered before,
        ! unless they are much more abundant, with Newton-Raphson 
        !-----------------------------------------------------------
        nsave = anmono
 150    continue
        eact(:) = .false.
        Nact = 0
        !if (verbose>0) print*,"NewFastLevel,ptake=",NewFastLevel,ptake
        do iredo=MAX(1,ido-NewBackIt),ido
          e = eseq(iredo)
          if (norm(e)<NewBackFac*norm(enew)) then
            eact(e) = .true. 
            Nact = Nact+1
            all_to_act(e) = Nact
            act_to_all(Nact) = e
            pbefore(e) = anmono(e)
            if (NewFastLevel<3.and.ptake) then
              anmono(e) = anmono(e)*pcorr(enew,e)
              Ntaken = Ntaken+1
            else
              pcorr(enew,e) = 1.Q0
            endif
            Nestim = Nestim+1
          endif
        enddo
        if (verbose>1) then
          print*,catm(eseq(1:ido))
          print*,eact(eseq(1:ido))
        endif
        if (Nact==1) then
          pkey(ido) = enew
          e = eseq(ido)
          xx(e) = LOG(anmono(e)*kT)
          cycle
        endif  
        do i=1,nml
          affect = .false. 
          known  = .true.
          do j=1,m_kind(0,i)
            e = m_kind(j,i) 
            if (.not.done(e)) then
              known = .false.
              exit
            endif
            if (eact(e)) affect=.true.
          enddo  
          relevant(i) = (known.and.affect)
        enddo

        !---------------------------------------------------------
        !***  check if new important molecule has huge impact  ***
        !---------------------------------------------------------
        i = imaj(enew)
        if (i>0) then
          pmol = g(i)
          do j=1,m_kind(0,i)
            pmol = pmol + m_anz(j,i)*LOG(anmono(m_kind(j,i))*kT)
          enddo
          pmol = EXP(pmol)
          i1 = 0
          s1 = 0
          qmost = 0.Q0
          do j=1,m_kind(0,i)
            e = m_kind(j,i) 
            if (e==el) cycle
            if (e==enew) then
              i1 = enew
              s1 = m_anz(j,i) 
              cycle
            endif
            qq = pmol*m_anz(j,i)/(eps(e)*anHges*kT)
            !print*,catm(e),cmol(i),REAL(qq)
            if (qq>qmost) then
              qmost = qq
              i2 = e
              s2 = m_anz(j,i)
            endif  
          enddo
          possible = (qmost>0.95)
          if (possible) then
            imaj2(i2) = imaj(i2)
            imaj(i2) = i
            pmol = eps(i2)*anHges*kT/s2
            if (eps(i1)*anHges*kT-s1*pmol<0.Q0.and.imaj2(i2)>0) then
              tmp= i1
              i1 = i2
              i2 = tmp
              tmp= s1
              s1 = s2
              s2 = tmp
            endif
            pmol = eps(i2)*anHges*kT/s2
            if (eps(i1)*anHges*kT-s1*pmol<0.Q0) possible=.false.
            if (i2==el) possible=.false.
          endif  
          if (possible) then
            !print*,catm(i1),cmol(i),cmol(imaj2(i1)),
     >      !       REAL(eps(i1)*anHges*kT-s1*pmol)
            !print*,catm(i2),cmol(i),cmol(imaj2(i2)),
     >      !       REAL(eps(i2)*anHges*kT-s2*pmol)
            i3 = imaj2(i1)      ! second most important carrier of i1
            s3 = 0
            s4 = 0
            lnp3 = g(i3)
            do j=1,m_kind(0,i3)
              e = m_kind(j,i3) 
              if (e==i1) then
                s3 = m_anz(j,i3)
                cycle
              endif  
              if (e==i2) then
                s4 = m_anz(j,i3)
                cycle
              endif  
              lnp3 = lnp3 + m_anz(j,i3)*LOG(anmono(e)*kT)
            enddo
            lnp = g(i)
            do j=1,m_kind(0,i)
              e = m_kind(j,i)
              if (e==i1.or.e==i2) cycle
              lnp = lnp + m_anz(j,i)*LOG(anmono(e)*kT)
            enddo
            if (verbose>1) print*,"readjusting "//catm(i1)//" "
     >                   //catm(i2)
     >                   //" with molecules "//trim(cmol(i))//" "
     >                   //trim(cmol(i3))//" ..."
            !--- solve eps(i1)*anHges*kT = s1*pmol + s3*p_second ---
            !--- with  eps(i2)*anHges*kT = s2*pmol               ---
            FF(1) = LOG((eps(i1)*anHges*kT-s1*pmol)/s3)-lnp3
            FF(2) = LOG((eps(i2)*anHges*kT)/s2)-lnp
            DF(1,1) = s3
            DF(1,2) = s4
            DF(2,1) = s1
            DF(2,2) = s2
            call GAUSS16(nel,2,DF,dp,FF)
            possible = (dp(1)<1000.0).and.(dp(2)<1000.0)
          endif
          if (possible) then
            anmono(i1) = EXP(dp(1))*kT1
            anmono(i2) = EXP(dp(2))*kT1
            !print'(A2,9(1pE16.8))',catm(i1),eps(i1)*anHges*kT,
     >      !                       s1*EXP(lnp+s1*dp(1)+s2*dp(2))
     >      !                      +s3*EXP(lnp3+s3*dp(1)+s4*dp(2))
            !print'(A2,9(1pE16.8))',catm(i2),eps(i2)*anHges*kT,
     >      !                       s2*EXP(lnp+s1*dp(1)+s2*dp(2))
            pbefore(i1) = anmono(i1)
            pbefore(i2) = anmono(i2)
            if (NewFastLevel<3.and.ptake) then
              anmono(i1) = anmono(i1)*pcorr(enew,i1)
              anmono(i2) = anmono(i2)*pcorr(enew,i2)
            endif  
            DUALco = DUALco+1
          endif  
          !--- special case CH4+CO2+H2O ---
          HCOproblem = (Nact>2).and.(enew==H.or.enew==C.or.enew==O)
     >                 .and.(CH4>0).and.(H2O>0).and.(CO2>0)
          if (HCOproblem) then
            hasH = .false.
            hasC = .false.
            hasO = .false.
            do ii=1,Nact
              i = act_to_all(ii)
              if (i==H) hasH=.true.
              if (i==C) hasC=.true.
              if (i==O) hasO=.true.
            enddo  
            HCOproblem = hasH.and.hasO.and.hasC
          endif  
          if (HCOproblem) then
            pH2O = (2*eps(O)-4*eps(C)+eps(H))*anHges*kT/4
            pCO2 = (2*eps(O)+4*eps(C)-eps(H))*anHges*kT/8
            pCH4 = (eps(H)+4*eps(C)-2*eps(O))*anHges*kT/8
            HCOproblem = (pH2O>0).and.(pCO2>0).and.(pCH4>0)
          endif  
          if (HCOproblem) then
            lpH2O = LOG(pH2O)
            lpCO2 = LOG(pCO2)
            lpCH4 = LOG(pCH4)
            lpH = (lpCH4+lpH2O*2+g(CO2)-lpCO2-g(CH4)-g(H2O)*2)/8
            lpO = (lpCO2+lpH2O*2+g(CH4)-lpCH4-g(CO2)-g(H2O)*2)/4
            lpC = (lpCH4+lpCO2+g(H2O)*2-lpH2O*2-g(CH4)-g(CO2))/2
            if (verbose>1) print*,"applying H2O-CH4-CO2 correction ..."
            !print'(2(1pE16.8))',2*pH2O+4*pCH4,anHges*kT*eps(H)
            !print'(2(1pE16.8))',  pH2O+2*pCO2,anHges*kT*eps(O)
            !print'(2(1pE16.8))',  pCH4+  pCO2,anHges*kT*eps(C)
            !print'(2(1pE16.8))',pH2O,EXP(lpH*2+lpO+g(H2O))
            !print'(2(1pE16.8))',pCO2,EXP(lpC+lpO*2+g(CO2))
            !print'(2(1pE16.8))',pCH4,EXP(lpC+lpH*4+g(CH4))
            anmono(H) = EXP(lpH)*kT1
            anmono(O) = EXP(lpO)*kT1
            anmono(C) = EXP(lpC)*kT1
            pbefore(H) = anmono(H)
            pbefore(O) = anmono(O)
            pbefore(C) = anmono(C)
            if (NewFastLevel<3.and.ptake) then
              anmono(H) = anmono(H)*pcorr(enew,H)
              anmono(O) = anmono(O)*pcorr(enew,O)
              anmono(C) = anmono(C)*pcorr(enew,C)
            endif  
            HCOco = HCOco + 1
          endif  
        endif  
          
        !--- solve with different methods ---
        if (verbose>1) then 
          print'("corr",99(1pE11.2))',pcorr(enew,act_to_all(1:Nact))
        endif  
        qual = 9.Q+99
        do iloop=1,4
          if (iloop==1) then
            imethod = NewPreMethod                    ! first choice
          elseif (iloop==2) then
            imethod = NewPreMethod                    ! second try with pcorr=1.0
          elseif (iloop==3) then
            imethod = MOD(NewPreMethod,3)+1
          elseif (iloop==4) then
            imethod = MOD(NewPreMethod+1,3)+1
          endif
          if (verbose>1) print*,"imethod=",iloop,imethod
          !write(97,*) catm(eseq(1:ido))
          !write(97,*) eact(eseq(1:ido))
          !write(97,*) imethod
          !write(97,*) real(pcorr(enew,act_to_all(1:Nact)))          
          if (imethod==1) then
            !-------- method 1: xx=log(patm)-variables --------
            null = anmono
            qual0 = qual 
            qual1 = qual 
            qual2 = qual 
            do it=1,299
              do ii=1,Nact
                i = act_to_all(ii)
                xx(i) = LOG(anmono(i)*kT)
                xscale(i) = anmono(i)*kT
                if (i==el) then
                  fscale(ii) = 1.Q0/mcharge
                else  
                  fscale(ii) = 1.Q0/(anHges*eps(i)*kT)
                endif  
                FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
                DF(ii,:)  = 0.Q0
                DF(ii,ii) = -fscale(ii)*anmono(i)*kT
              enddo  
              do i=1,nml
                if (.not.relevant(i)) cycle 
                pmol = g(i)
                do j=1,m_kind(0,i)
                  pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
                enddo
                pmol = EXP(pmol)
                do j=1,m_kind(0,i)
                  m1 = m_kind(j,i)
                  if (.not.eact(m1)) cycle
                  ii = all_to_act(m1)
                  term = fscale(ii)*m_anz(j,i)*pmol
                  FF(ii) = FF(ii) - term
                  do l=1,m_kind(0,i)
                    m2 = m_kind(l,i)
                    if (.not.eact(m2)) cycle
                    jj = all_to_act(m2)
                    DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
                  enddo	    
                enddo
              enddo
              call GAUSS16(nel,Nact,DF,dp,FF)
              bem = " "
              qual3 = qual2
              qual2 = qual1
              qual1 = qual0
              qual0 = qual
              qual  = 0.Q0
              do ii=1,Nact
                qual = qual + ABS(dp(ii))
              enddo  
              maxs = 3.0
              if (it>30.and.(qual >qual0.or.qual0>qual1.or.
     >                       qual1>qual2.or.qual2>qual3)) then
                maxs = 3.0*exp(-MAX(0,it-30)/70.0)
                bem = "*"
              endif  
              do ii=1,Nact
                i = act_to_all(ii) 
                xx(i) = xx(i) - MAX(-maxs,MIN(maxs,dp(ii)))
                anmono(i) = exp(xx(i))*kT1
              enddo
              !write(97,'(I4,A2,99(1pE11.3E3))')
     >        !               it,bem,anmono(act_to_all(1:Nact))*kT,qual
              if (verbose>1) print'(I4,A2,99(1pE11.3E3))',
     >                       it,bem,anmono(act_to_all(1:Nact))*kT,qual
              NpreIt = NpreIt+1
              if (it>1.and.qual<1.Q-12) exit
            enddo
            NpreLoop = NpreLoop+1

          else if (imethod==2) then
            !-------- method 2: lin-variables with pullback --------
            dp(:) = 0.Q0
            qual = 9.Q+99
            fak = 1.Q0
            do it=1,199
              qual0 = qual 
              null  = anmono
              pullmax = 1
              if (it>100) pullmax=10
              do ipull=1,pullmax  ! pullback if quality gets worse
                !--- make a step ---
                do ii=1,Nact
                  i = act_to_all(ii)
                  anmono(i) = null(i)-fak*dp(ii)*kT1
                  xx(i) = LOG(anmono(i)*kT)
                enddo  
                !--- determine new FF and DF ---
                do ii=1,Nact
                  i = act_to_all(ii) 
                  xscale(i) = anmono(i)*kT
                  if (i==el) then
                    fscale(ii) = 1.Q0/mcharge
                  else  
                    fscale(ii) = 1.Q0/(anHges*eps(i)*kT)
                  endif  
                  FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
                  DF(ii,:) = 0.Q0
                  DF(ii,ii) = -fscale(ii)*anmono(i)*kT
                enddo
                do i=1,nml
                  if (.not.relevant(i)) cycle 
                  pmol = g(i)
                  do j=1,m_kind(0,i)
                    pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
                  enddo
                  pmol = EXP(pmol)
                  do j=1,m_kind(0,i)
                    m1 = m_kind(j,i)
                    if (.not.eact(m1)) cycle
                    ii = all_to_act(m1)
                    term = fscale(ii)*m_anz(j,i)*pmol
                    FF(ii) = FF(ii) - term
                    do l=1,m_kind(0,i)
                      m2 = m_kind(l,i)
                      if (.not.eact(m2)) cycle
                      jj = all_to_act(m2)
                      DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
                    enddo	    
                  enddo
                enddo
                !--- determine new quality ---
                qual = 0.Q0
                do ii=1,Nact
                  i = act_to_all(ii)
                  qual = qual + FF(ii)**2
                enddo  
                if (qual<qual0) exit
                if (ipull==pullmax) exit
                !write(97,'("pullback",3(1pE11.3))') fak,qual0,qual
                if (verbose>1) print'("pullback",3(1pE11.3))',
     >                         fak,qual0,qual
                fak = 0.5*fak   ! reduce NR-step
              enddo  
              !write(97,'(I4,99(1pE11.3))')
     >        !             it,anmono(act_to_all(1:Nact))*kT,qual
              if (verbose>1) print'(I4,99(1pE11.3))',
     >                     it,anmono(act_to_all(1:Nact))*kT,qual
              NpreIt = NpreIt+1
              if (it>1.and.qual<1.Q-12) exit
              !--- determine new NR-vector ---
              call GAUSS16(nel,Nact,DF,dp,FF)
              do ii=1,Nact
                i = act_to_all(ii) 
                dp(ii) = dp(ii)*xscale(i)
              enddo
              null = anmono
              !--- limit step physically, keep direction ---
              fak = 1.Q0
              do ii=1,Nact
                i = act_to_all(ii)
                if (null(i)*kT-fak*dp(ii)>5.Q0*null(i)*kT) then
                  fak=MIN(fak,-4.Q0*null(i)*kT/dp(ii))
                endif
                if (null(i)*kT-fak*dp(ii)<0.2Q0*null(i)*kT) then
                  fak=MIN(fak,0.8Q0*null(i)*kT/dp(ii))
                endif
              enddo
            enddo  
            NpreLoop = NpreLoop+1

          else if (imethod==3) then
            !-------- method 3: xx=log(patm)-variables --------
            null = anmono
            do it=1,299
              do ii=1,Nact
                i = act_to_all(ii)
                xx(i) = LOG(anmono(i)*kT)
                xscale(i) = anmono(i)*kT
                if (i==el) then
                  fscale(ii) = 1.Q0/mcharge
                else  
                  fscale(ii) = 1.Q0/(anHges*eps(i)*kT)
                endif  
                FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
                DF(ii,:)  = 0.Q0
                DF(ii,ii) = -fscale(ii)*anmono(i)*kT
              enddo  
              do i=1,nml
                if (.not.relevant(i)) cycle 
                pmol = g(i)
                do j=1,m_kind(0,i)
                  pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
                enddo
                pmol = EXP(pmol)
                do j=1,m_kind(0,i)
                  m1 = m_kind(j,i)
                  if (.not.eact(m1)) cycle
                  ii = all_to_act(m1)
                  term = fscale(ii)*m_anz(j,i)*pmol
                  FF(ii) = FF(ii) - term
                  do l=1,m_kind(0,i)
                    m2 = m_kind(l,i)
                    if (.not.eact(m2)) cycle
                    jj = all_to_act(m2)
                    DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
                  enddo	    
                enddo
              enddo
              call GAUSS16(nel,Nact,DF,dp,FF)
              qual = 0.Q0
              do ii=1,Nact
                qual = MAX(qual,ABS(dp(ii)))
              enddo  
              fak = MIN(1.Q0,3.Q0/qual)
              do ii=1,Nact
                i = act_to_all(ii) 
                xx(i) = xx(i) - fak*dp(ii)
                anmono(i) = exp(xx(i))*kT1
                if (IS_NAN(REAL(anmono(i)))) then
                  qual=2.Q+99
                  exit
                endif  
              enddo
              !write(97,'(I4,A2,99(1pE11.3E3))')
     >        !               it,bem,anmono(act_to_all(1:Nact))*kT,qual
              if (verbose>1) print'(I4,99(1pE11.3E3))',
     >                       it,anmono(act_to_all(1:Nact))*kT,fak,qual
              NpreIt = NpreIt+1
              if (qual<1.Q-12.or.qual>1.Q+99) exit
            enddo  
            NpreLoop = NpreLoop+1
          endif
          if (qual<1.Q-12) exit
          do ii=1,Nact
            i = act_to_all(ii)
            anmono(i) = pbefore(i)    ! do not use pre-it corrections again
          enddo 
          if (verbose>1) read(*,'(A1)') char
        enddo  
        if (qual>1.Q-4) then
          if (ptake) then
            anmono = nsave 
            ptake = .false.
            goto 150
          endif  
          write(*,*) "*** no convergence in NR pre-it "
          print*,"Tg=",Tg
          print*,catm(eseq(1:ido))
          print*,eact(eseq(1:ido))
          print*,pcorr(enew,act_to_all(1:Nact))
          if (ifatal<=1) then
            from_merk = .false.
            badness   = 1.Q0
            pcorr     = 1.Q0
            ifatal    = ifatal+1
            goto 100            ! try again from scratch before giving up
          endif  
          !sortq(enew)=3.0*sortq(enew)
          !if (sortq(enew)<10.0) goto 150
          goto 1000
        endif  
        !--- save ratio after/before for next run ---
        !print*,"ido,pkey(ido),enew=",ido,pkey(ido),enew
        pkey(ido) = enew
        pcorr(enew,:) = 1.Q0
        do ii=1,Nact
          i = act_to_all(ii)
          pcorr(enew,i) = anmono(i)/pbefore(i)    
        enddo 
        if (verbose>1) print'("corr",99(1pE11.2))',
     >                 pcorr(enew,act_to_all(1:Nact))
        if (verbose>1) read(*,'(A1)') char
      enddo  
*
*     ! redo electron density
*     =======================
      if (charge) then
        coeff(:) = 0.Q0
        relevant = .false.
        do i=1,nml
          do j=1,m_kind(0,i)
            if (m_kind(j,i)==el) relevant(i)=.true. 
          enddo  
          if (.not.relevant(i)) cycle
          pmol = g(i)
          l = 0
          do j=1,m_kind(0,i)
            pat = anmono(m_kind(j,i))*kT
            if (m_kind(j,i)==el) then
              l = m_anz(j,i)   
            else 
              pmol = pmol + m_anz(j,i)*LOG(pat)
            endif  
          enddo
          coeff(l) = coeff(l)+EXP(pmol)
        enddo
        pel = SQRT(coeff(-1)/(1.Q0+coeff(+1)))     ! 0 = pel - a/pel + b*pel
        mcharge = MAX(pel,MAX(coeff(-1)/pel,coeff(1)*pel))
        if (verbose>1) print'(" redo pel =",1pE17.10," ->",1pE17.10)',
     >                        REAL(anmono(el)*kT),REAL(pel)
        anmono(el) = pel*kT1
        !pecorr = pel/peest
      endif  

*     ! use memory of deviations between predicted atom pressures 
*     ! and converged atom pressures to improve the predictions
*     ============================================================
      ansave = anmono
      if (NewFastLevel<2.and.ptake) anmono = anmono*badness
*     
 200  continue
*-----------------------------------------------------------------------
      if (NewFullIt) then
*       ! Jacobi matrix and rhs vector for Newton-Raphson
*       =================================================
        it = 0
        eact(:) = .true.
        conv(:,:) = 1.Q+99
        converge(:) = 1.Q+99
        switchoff(:) = 0
        finish = 1.Q-24
        crit  = 1.Q+99
        crit1 = 1.Q+99
        crit2 = 1.Q+99
 300    continue
        !if (it>30) finish=10.Q0**(-22.0+19.0*(it-30.0)/(itmax-30.0))
        Nact = 0
        ii = 0
        do i=1,nel
          if (.not.eact(i)) cycle
          Nact = Nact+1
          ii = ii+1
          all_to_act(i) = ii
          act_to_all(ii) = i
          xx(i) = LOG(anmono(i)*kT)
          xscale(i) = anmono(i)*kT
          if (i==el) then
            fscale(ii) = 1.Q0/mcharge
          else  
            fscale(ii) = 1.Q0/(anHges*eps(i)*kT)
          endif  
          FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
          DF(ii,:) = 0.Q0
          DF(ii,ii) = -fscale(ii)*anmono(i)*kT
        enddo	
        do i=1,nml
          pmol = g(i)
          do j=1,m_kind(0,i)
            pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
          enddo
          pmol = EXP(pmol)
          anmol(i) = pmol*kT1
          do j=1,m_kind(0,i)
            m1 = m_kind(j,i)
            if (.not.eact(m1)) cycle
            ii = all_to_act(m1)
            term = m_anz(j,i)*pmol*fscale(ii)
            FF(ii) = FF(ii) - term
            do l=1,m_kind(0,i)
              m2 = m_kind(l,i)
              if (.not.eact(m2)) cycle
              jj = all_to_act(m2)
              DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
            enddo	    
          enddo
        enddo

*       ! compute the Newton-Naphson step	  
*       =================================
        FF0 = FF
        DF0 = DF
        call GAUSS16(nel,Nact,DF,dp,FF)
        !do i=1,Nact
        !  print'(99(1pE11.3E3))',DF0(i,1:Nact),dp(i),FF0(i)
        !enddo  

*       ! apply limited NR step and check convergence 	  
*       =============================================
        if (.true.) then        
          qual = 0.Q0
          do ii=1,Nact
            qual = MAX(qual,ABS(dp(ii)))
          enddo  
          fak = MIN(1.Q0,3.Q0/qual)
          converge(it) = 0.Q0
          Nconv = 0
          if (verbose>0) txt = ""
          do i=1,nel
            if (.not.eact(i)) then
              Nconv = Nconv+1 
              if (verbose>0) txt = trim(txt)//" "//catm(i)
            else 
              ii = all_to_act(i) 
              delp = -dp(ii)                             ! relative change dx/x
              conv(it,i) = delp
              converge(it) = MAX(converge(it),ABS(delp))
              if (ABS(delp)<finish) then
                Nconv = Nconv+1 
                if (verbose>0) txt = trim(txt)//" "//catm(i)
              endif  
            endif  
          enddo
          if (verbose>1.and.it<100) then
            print*,"fak=",fak
            print'(7x,A14,A14,A14)',"patom","dpatom","badness" 
            do ii=1,Nact
              i = act_to_all(ii) 
              print'(A7,3(1pE14.6))',catm(i),anmono(i)*kT,
     >             -dp(ii),badness(i)
            enddo
          endif
          do ii=1,Nact
            i = act_to_all(ii) 
            xx(i) = xx(i) - fak*dp(ii)
            anmono(i) = exp(xx(i))*kT1
          enddo
          limit = fak

        else
          !--- re-scale ---
          do ii=1,Nact
            i = act_to_all(ii) 
            dp(ii) = dp(ii)*xscale(i)
          enddo  
          fak = 5.Q0
          limit = 1.Q0                                   ! limit step, keep direction
          converge(it) = 0.Q0
          Nconv = 0
          if (verbose>0) txt = ""
          do i=1,nel
            if (.not.eact(i)) then
              Nconv = Nconv+1 
              if (verbose>0) txt = trim(txt)//" "//catm(i)
            else 
              ii = all_to_act(i) 
              delp = -dp(ii)/(anmono(i)*kT)              ! relative change dx/x
              conv(it,i) = delp
              converge(it) = MAX(converge(it),ABS(delp))
              if (ABS(delp)<finish) then
                Nconv = Nconv+1 
                if (verbose>0) txt = trim(txt)//" "//catm(i)
              endif  
              if (1.Q0+delp>fak) then
                limit = min(limit,(fak-1.Q0)/delp)       ! such that xnew=xold*fac 
              else if (1.Q0+delp<1.Q0/fak) then
                limit = min(limit,(1.Q0/fak-1.Q0)/delp)  ! such that xnew=xold/fac
              endif
            endif  
          enddo
          if (it<=10) then
            limit = 1.Q0
          else
            dp = dp*limit
          endif  
          if (verbose>1.and.it==0) then
            write(*,*) 
            print'(7x,A14,A14,A14)',"patom","dpatom","badness" 
            do ii=1,Nact
              i = act_to_all(ii) 
              print'(A7,3(1pE14.6))',catm(i),anmono(i)*kT,
     >             -dp(ii)/(anmono(i)*kT),badness(i)
            enddo
          endif

*         ! apply limited NR step
*         =======================
          !fak = 1.Q0+4.Q0*EXP(-(MAX(0,it-20))/13.Q0)
          do ii=1,nact
            i = act_to_all(ii)
            delp = -dp(ii)*kT1
            nold = anmono(i)
            anmono(i) = MAX(nold/fak,MIN(nold*fak,nold+delp))
          enddo
          if (it>itmax-10) then
c            verbose = 2
            do ii=1,Nact
              i = act_to_all(ii) 
              delp = -dp(ii)/(anmono(i)*kT)
              print'(A3,99(1pE12.4))',catm(i),anmono(i),delp
            enddo  
          endif  

        endif
        crit2 = crit1
        crit1 = crit
        crit  = converge(it)
        if (verbose>1) print'(i3,i3,2(1pE9.1)," converged(",i2,"):",
     >                    A50)',it,Nact,converge(it),limit,Nconv,txt
        if (it==itmax) then 
          write(*,*) '*** keine Konvergenz in SMCHEM16!'
          write(*,*) 'it, converge, ind =',it,converge(it),limit
          write(*,*) '  n<H>, T =',anhges,Tg
          write(*,*) 'from_merk,NewFastLevel=',from_merk,NewFastLevel
          if (ifatal<=1) then
            chemiter  = chemiter + it
            from_merk = .false.
            badness   = 1.Q0
            pcorr     = 1.Q0
            ifatal    = ifatal+1
c            if (ifatal==2) verbose=2
            goto 100        ! try again from scratch before giving up
          endif  
          goto 1000
        endif
        if (it>=5) then
          j = 0 
          do ii=1,Nact
            i = act_to_all(ii)
            if (MAXVAL(ABS(conv(it-5:it,i)))<finish) then
              switchoff(i) = it
              eact(i) = .false.
              j = j+1
              if (verbose>1) then
                print*,"switching off "//catm(i)//" ..."
              endif  
            endif
          enddo
          Nact = Nact-j
        endif  
        it = it + 1
        if (verbose.gt.1) read(*,'(a1)') char
        if (it==1) goto 300                                       ! at least two iterations
        if ((crit>finish*1.Q-5).and.(Nact>0).and.
     >      (MAX(crit,crit1,crit2)>finish*1.Q+4).and.
     >      (MAX(crit,crit1,crit2)>finish*1.Q+8.or.it<100)) goto 300   ! continue iterating
*
*       ! redo rare elements
*       ====================
        redo(:) = .false.
        do iredo=1,nel
          atmax = 0.Q0 
          e = 0
          do i=1,nel
            xx(i) = LOG(anmono(i)*kT)
            if (redo(i)) cycle   
            atfrac = anmono(i)/anHges
            if (atfrac>1.Q-30) cycle   
            if (atfrac<atmax) cycle
            atmax = atfrac
            e = i
          enddo  
          if (e==0) exit
          redo(e) = .true.
          psc = xx(e)
          coeff(:) = 0.Q0
          lmin = +99
          lmax = -99
          do i=1,nml
            pmol = g(i)
            affect = .false.
            do j=1,m_kind(0,i)
              nu = m_anz(j,i)   
              if (m_kind(j,i)==e) then
                l = nu
                affect = .true.
                pmol = pmol + l*psc
              else
                pmol = pmol + nu*xx(m_kind(j,i))
              endif  
            enddo
            if (.not.affect) cycle
            lmin = MIN(lmin,l)
            lmax = MAX(lmax,l)
            coeff(l) = coeff(l)+EXP(pmol)
          enddo
          psc = EXP(psc)
          if (verbose>1) then
            print'("redo rare element patom dp/patom ",A2)',catm(e)
            print*,"psc=",psc
            print'(" coeff(",I2,":",I2,")=",99(1pE10.2))',
     >                   lmin,lmax,coeff(lmin:lmax)
          endif  
          pat = anmono(e)*kT/psc
          do piter=1,99
            f  = pat*psc-eps(e)*anHges*kT
            fs = psc
            do l=lmin,lmax
              if (coeff(l)==0.Q0) cycle
              f  = f  + coeff(l)*l*pat**l
              fs = fs + coeff(l)*l**2*pat**(l-1)
            enddo
            delta = f/fs
            if (verbose>1) print'(1x,A2,1pE25.15,1pE10.2)',
     >                            catm(e),pat*psc,delta/pat
            pat = pat-delta
            if (ABS(delta)<finish*ABS(pat)) exit 
          enddo  
          if (piter>=99) then
            write(*,*) "*** no convergence in post-it "//catm(e)
            write(*,*) anmono(e),eps(e)
            write(*,*) psc,lmin,lmax
            write(*,*) anHges,Tg
            write(*,*) coeff(lmin:lmax)
            if (ifatal<=1) then
              from_merk = .false.
              badness   = 1.Q0
              pcorr     = 1.Q0
              ifatal    = ifatal+1
              goto 100        ! try again from scratch before giving up
            endif  
            goto 1000
          endif  
          anmono(e) = pat*psc*kT1
        enddo
*
*       ! how bad was the initial guess?
*       ================================
        if (.not.from_merk) then
          if (verbose>1) print'(7x,3(A14))',
     >                   "patom","conv.","init.guess"
          do i=1,nel
            badness(i) = anmono(i)/ansave(i)
            switch = switchoff(i)
            if (switch==0) switch=it-1
            if (verbose>1) then
              print'(A7,3(1pE14.6))',catm(i),anmono(i)*kT,
     >              conv(switch,i),badness(i)
            endif  
          enddo
!$omp critical(fort99)
          ilauf = ilauf+1
          !if (ilauf==1) write(99,'(A9,A10,A4,99(A10))') 
     >    !      'Tg','n<H>','it',catm(1:nel)
          !write(99,'(0pF9.3,1pE10.3,I4,99(1pE10.3))') 
     >    !      Tg,anHges,it,badness
!$omp end critical(fort99)
        endif  

      endif     ! NewFullIt

*     ! final anmol determination
*     ===========================
      do i=1,nel
        xx(i) = LOG(anmono(i)*kT)
      enddo  
      amerk = anmono/anHges
      anmol = 0.Q0
      do i=1,nml
        pmol = g(i)
        do j=1,m_kind(0,i)
          pat = anmono(m_kind(j,i))*kT
          pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
        enddo
        anmol(i) = EXP(pmol)*kT1
      enddo
      if (charge) pel=anmono(el)*kT

      if (ngestst) then
*       ! Test auf Elementerhaltung
*       ===========================
        do i=1,nel
          nges(i) = anmono(i)
        enddo
        do i=1,nml
          do j=1,m_kind(0,i)
            j1 = m_kind(j,i)
            nges(j1) = nges(j1) + m_anz(j,i)*anmol(i)
          enddo
        enddo
        do e=1,nel
          if (e==el) cycle 
          soll  = anHges * eps(e)
          haben = nges(e)
          abw   = ABS(soll-haben)/MAX(soll,haben)
          if (abw>1.Q-15) then
            if (verbose>1) then
              print'("*** element conservation error ",A2)',catm(e)
              print'(A12,1pE14.7)',catm(e),anmono(e)/(eps(e)*anHges)
            endif  
            sum = anmono(e)/(eps(e)*anHges)
            do i=1,nml
              do j=1,m_kind(0,i)
                j1 = m_kind(j,i)
                cc = m_anz(j,i)*anmol(i)/(eps(e)*anHges)
                if (j1==e.and.cc>1.Q-7) then
                  if (verbose>1) print'(A12,1pE14.7)',cmol(i),cc
                  sum = sum+cc
                endif  
              enddo
            enddo
            if (verbose>1) then
              print'(3(1pE14.7))',soll/anHges,haben/anHges,sum
            endif  
            from_merk = .false.
            ansave = anmono
c            verbose=2
            goto 200
          endif
        enddo
      endif
      
      if (verbose.gt.0) print '("  ==> smchem used it=",I3,
     &     " conv=",3(1pE9.2))',it,crit,crit1,crit2
      if (verbose.gt.0) print'("number of pre-iterations",I4,
     &     " -- used saved initial guesses",0pF5.1,"%")',
     &     NpreIt,REAL(Ntaken)/REAL(Nestim+1)*100
      if (verbose.gt.1) read(*,'(a1)') char

!$omp critical(counters)
      chemcall = chemcall + 1
      chemiter = chemiter + it
      preIter  = preIter  + NpreIt
      preUse   = preUse   + Ntaken
      preEst   = preEst   + Nestim
      DUALcorr = DUALcorr + DUALco
      HCOcorr  = HCOcorr  + HCOco
!$omp end critical(counters)

      return

 1000 continue
      open(unit=12,file='fatal.case')
      do i=1,nel
        write(12,'(A2,1x,0pF32.26)') catm(i),12.Q0+log10(eps(i))
      enddo  
      write(12,*) anhges,Tg
      close(12)
      print*, "***  giving up."
      return


      CONTAINS       ! internal functions - not visible to other units 
************************************************************************
      REAL(kind=qp) FUNCTION gk(i)
************************************************************************
*****  returns ln(kp) [cgs] for different fit formula              *****
************************************************************************
      use CHEMISTRY,ONLY: a,th1,th2,th3,th4,TT1,TT2,TT3,fit,natom,cmol,
     >                    NELEM,elnum,b_nasa,c_nasa
      implicit none
      integer,intent(in) :: i                 ! index of molecule
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp),parameter :: bar=1.Q+6, atm=1.013Q+6, Rcal=1.987Q+0
      real(kind=qp),parameter :: Rgas=8.3144598Q+0
      real(kind=qp),parameter :: ln10=LOG(10.Q0)
      real(kind=qp),parameter :: lnatm=LOG(atm), lnbar=LOG(bar)
      real(kind=qp) :: lnk,dG 
      real(kind=qp) :: h_rt,s_r               !Added by Yui Kawashima
      real(kind=qp) :: dG_rt_ref(NELEM),dG_rt !Added by Yui Kawashima
      integer:: k,j                           !Added by Yui Kawashima

      if (i.eq.0) then
        gk = -1000.Q0                  ! tiny kp for unassigned molecules
        return
      endif
      if (fit(i).eq.1) then
        !---------------------
        ! ***  Gail's fit  *** 
        !---------------------
        lnk = a(i,0) + a(i,1)*th1 + a(i,2)*th2 
     &               + a(i,3)*th3 + a(i,4)*th4 

      else if (fit(i).eq.2) then
        !---------------------------
        ! ***  Tsuji (1973) fit  *** 
        !---------------------------
        lnk = ln10*( - a(i,0) - a(i,1)*th1 - a(i,2)*th2
     &                        - a(i,3)*th3 - a(i,4)*th4 ) 

      else if (fit(i).eq.3) then  
        !---------------------------------
        ! ***  Sharp & Huebner (1990)  ***
        !---------------------------------
        dG  = a(i,0)/TT1 + a(i,1) + a(i,2)*TT1 + a(i,3)*TT2 + a(i,4)*TT3
        lnk = -dG/(Rcal*TT1) + (1-Natom(i))*lnatm

      else if (fit(i).eq.4) then
        !-----------------------------------
        ! ***  Stock (2008) & Kietzmann  ***
        !-----------------------------------
        dG  = a(i,0)/TT1+a(i,1)*LOG(TT1)+a(i,2)+a(i,3)*TT1+a(i,4)*TT2
        lnk = dG + (1-Natom(i))*lnbar

      else if (fit(i).eq.5) then
        !--------------------
        ! ***  dG(T)-fit  ***
        !--------------------
        dG  = a(i,0)/TT1 + a(i,1) + a(i,2)*TT1 + a(i,3)*TT2 + a(i,4)*TT3
        lnk = -dG/(Rgas*TT1) + (1-Natom(i))*lnbar
         
      else if (fit(i).eq.6) then
        !-------------------------------
        ! ***  Barklem & Collet fit  ***
        !-------------------------------
        lnk = a(i,0)/TT3 + a(i,1)/TT2 + a(i,2)/TT1 + a(i,3)/TT1**0.05d0
     &      + a(i,4)*LOG(TT1) + a(i,5) + a(i,6)*TT1 + a(i,7)*TT2

      else if (fit(i).eq.7) then
        !-----------------------------------------------------
        ! ***  NASA polynomial fit added by Yui Kawashima  ***
        !-----------------------------------------------------         
        if (TT1>1.Q3) then
          h_rt = a(i,0) + a(i,1)*TT1/2.Q0 
     &         + a(i,2)*TT1**2/3.Q0 + a(i,3)*TT1**3/4.Q0
     &         + a(i,4)*TT1**4/5.Q0 + a(i,5)/TT1
          s_r  = a(i,0)*log(TT1) + a(i,1)*TT1
     &         + a(i,2)*TT1**2/2.Q0 + a(i,3)*TT1**3/3.Q0
     &         + a(i,4)*TT1**4/4.Q0 + a(i,6)
        else
          h_rt = a(i,7) + a(i,8)*TT1/2.Q0
     &         + a(i,9) *TT1**2/3.Q0 + a(i,10)*TT1**3/4.Q0
     &         + a(i,11)*TT1**4/5.Q0 + a(i,12)/TT1
          s_r  = a(i,7)*log(TT1) + a(i,8)*TT1
     &         + a(i,9) *TT1**2/2.Q0 + a(i,10)*TT1**3/3.Q0
     &         + a(i,11)*TT1**4/4.Q0 + a(i,13)           
        endif
        dG_rt = h_rt - s_r
        do k=1,m_kind(0,i)
          j = elnum(m_kind(k,i))
          if (c_nasa(j)==0) then
            print*,"Provide the data in data/Burcat_ref-elements.dat"
     &            ," and edit nasa_polynomial.f for "
     &            ,trim(catm(m_kind(k,i)))
            stop
          else
            if (TT1>1.Q3) then
              h_rt = b_nasa(j,0) + b_nasa(j,1)*TT1/2.Q0 
     &             + b_nasa(j,2)*TT1**2/3.Q0
     &             + b_nasa(j,3)*TT1**3/4.Q0
     &             + b_nasa(j,4)*TT1**4/5.Q0 + b_nasa(j,5)/TT1
              s_r  = b_nasa(j,0)*log(TT1) + b_nasa(j,1)*TT1
     &             + b_nasa(j,2)*TT1**2/2.Q0
     &             + b_nasa(j,3)*TT1**3/3.Q0
     &             + b_nasa(j,4)*TT1**4/4.Q0 + b_nasa(j,6)
            else
              h_rt = b_nasa(j,7) + b_nasa(j,8)*TT1/2.Q0
     &             + b_nasa(j,9) *TT1**2/3.Q0
     &             + b_nasa(j,10)*TT1**3/4.Q0
     &             + b_nasa(j,11)*TT1**4/5.Q0
     &             + b_nasa(j,12)/TT1
              s_r  = b_nasa(j,7)*log(TT1) + b_nasa(j,8)*TT1
     &             + b_nasa(j,9) *TT1**2/2.Q0
     &             + b_nasa(j,10)*TT1**3/3.Q0
     &             + b_nasa(j,11)*TT1**4/4.Q0 + b_nasa(j,13)           
            endif
            dG_rt_ref(j) = h_rt - s_r
            dG_rt = dG_rt - m_anz(k,i)*dG_rt_ref(j)
          endif
        enddo
        dG = -dG_rt
        lnk = dG + (1-Natom(i))*lnbar

      else if (fit(i).eq.8) then
        !-----------------------------------------------
        ! ***  NASA 7-polynomial fit for BURCAT data ***
        !-----------------------------------------------
        if (TT1>1.Q3) then
          dG = a(i,0) *(log(TT1)-1.Q0) 
     &       + a(i,1) *TT1   / 2.Q0
     &       + a(i,2) *TT1**2/ 6.Q0    
     &       + a(i,3) *TT1**3/12.Q0
     &       + a(i,4) *TT1**4/20.Q0    
     &       - a(i,5) /TT1  
     &       + a(i,6)
        else
          dG = a(i,7) *(log(TT1)-1.Q0) 
     &       + a(i,8) *TT1   / 2.Q0
     &       + a(i,9) *TT1**2/ 6.Q0    
     &       + a(i,10)*TT1**3/12.Q0
     &       + a(i,11)*TT1**4/20.Q0    
     &       - a(i,12)/TT1  
     &       + a(i,13)
        endif
        lnk = dG + (1-Natom(i))*lnbar

      else
        print*,cmol(i),"i,fit=",i,fit(i)
        stop "???"
      endif  

      gk = lnk
      end FUNCTION gk

      end SUBROUTINE smchem16
************************************************************************
      SUBROUTINE smchem8 (anHges,Tg,eps,anmono,anmol,merk,verbose)
************************************************************************
*                                                                      *
*     small chemistry                                                  *
*     ---------------                                                  *
*     Diese Routine berechnet eine GG-Chemie                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     e i n g a b e :                                                  *
*     anHges : totale Dichte der Wasserstoffatome (~ nh+2*nh2)         *
*     tg     : Gastemperatur                                           *
*     eps    : Vektor mit linearen chemischen Haeufigkeiten ([H]=1)    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     a u s g a b e :                                                  *
*     anmono : vektor mit den dichten der monomere  falls relevant     *
*     anmol  : vektor mit den dichten der molekuele falls relevant     *
*                                                                      *
************************************************************************
      use CHEMISTRY,ONLY: NewBackIt,NewFullIt,NewBackFac,
     >                    NewFastLevel,NewPreMethod,
     >                    nml=>NMOLE,nel=>NELM,cmol,catm,
     >                    m_kind,m_anz,charge,elion,el,
     >                    Natom,Ncmax,STOImax,H,C,N,O,
     >                    th1,th2,th3,th4,TT1,TT2,TT3
      use EXCHANGE,ONLY: chemcall,chemiter,preUse,preIter,preEst,
     >                   DUALcorr,HCOcorr
      implicit none
*-----------------------------------------------------------------------
*  Dimensionierung fuer die Molekuel- und Atom Felder. Hier ist auf
*  Konsistenz mit dem aufrufenden Programm zu achten.
*-----------------------------------------------------------------------
      real*8,intent(in)  :: anHges,Tg
      real*8,intent(in)  :: eps(nel)
      real*8,intent(out) :: anmono(nel),anmol(nml)
      integer,intent(inout) :: verbose
      real*8,parameter :: bk=1.380662d-16
*-----------------------------------------------------------------------
*  Bei merk=.true. merkt sich die Routine die letzte konvergiert Loesung
*  und geht beim naechsten Mal von diesen Startwerten aus.
      logical,intent(INOUT) :: merk
*-----------------------------------------------------------------------
*  Die Variable "ngestst" entscheidet, ob die Elementerhaltung ueber-
*  prueft werden soll.
      logical,parameter :: ngestst=.false.
*-----------------------------------------------------------------------
      integer :: stindex,Nconv,switch,ido,iredo,nu,tmp
      integer :: Nact,all_to_act(nel),act_to_all(nel),switchoff(nel)
      integer :: e,i,j,j1,ii,jj,l,it,m1,m2,piter,ifatal,ipull,pullmax
      integer :: lmin,lmax,Nseq,iloop,imethod,enew,i1,i2,i3,s1,s2,s3,s4
      integer :: eseq(nel),imaj(nel),imaj2(nel)
      integer :: NpreLoop,NpreIt,Ntaken,Nestim,DUALco,HCOco
      integer :: ind,indx(nel),info,ipvt(nel)
      integer,parameter :: itmax=200
      real*8 :: finish,qual,qual0,qual1,qual2,qual3,qmost,qq
      real*8 :: g(0:nml),limit
      real*8 :: pH2O,pCO2,pCH4,lpH2O,lpCO2,lpCH4,lpH,lpO,lpC
      real*8 :: kT,kT1,cc,nelek,ng,Sa,fak,lth,arg,term,f,fs
      real*8 :: pel,delta,pat,atfrac,atmax,lnp,lnp3
      real*8 :: nges(nel),coeff(-3:Ncmax),lnc(-3:Ncmax)
      real*8 :: DF(nel,nel),dp(nel),FF(nel),pmol,crit,delp,nold
      real*8 :: DF0(nel,nel),FF0(nel),xscale(nel),fscale(nel)
      real*8 :: nsave(nel),null(nel),sortq(nel)
      real*8 :: conv(0:itmax,nel),converge(0:itmax)
      real*8 :: soll,haben,abw,sum,maxs,mcharge
      real*8 :: pbefore(nel),norm(nel),xx(nel)
      real*8 :: emax,pges,pwork,pwork2,pp,psc,ptest,aim
      real*8 :: work(nel*(nel+1))
      real*8 :: condnum1,work2(nel)
      logical :: from_merk,eact(nel),redo(nel),done(nel),affect,known
      logical :: relevant(nml)
      logical :: possible,ptake,HCOproblem,hasH,hasC,hasO,IS_NAN
      character(len=5000) :: mols
      character(len=100) :: txt
      character(len=1) :: char,bem
      integer,save :: TiC,H2O,CH4,CO2,ilauf=0
      real*8,allocatable,save :: amerk(:),ansave(:)
      real*8,allocatable,save :: badness(:),pcorr(:,:) 
      integer,allocatable,save :: pkey(:)
!$omp threadprivate(TiC,ilauf,amerk,ansave,badness,pcorr,pkey)
*-----------------------------------------------------------------------      

      ifatal = 0
      if (.not.allocated(amerk)) then
        allocate(badness(nel),pcorr(nel,nel),pkey(nel),
     >           amerk(nel),ansave(nel))
        TiC = stindex(cmol,nml,'TIC    ')
        H2O = stindex(cmol,nml,'H2O    ')
        CH4 = stindex(cmol,nml,'CH4    ')
        CO2 = stindex(cmol,nml,'CO2    ')
        badness = 1.d0
        pcorr   = 1.d0
        pkey    = 0
      endif

*-----------------------------------------------------------------------
*     ! zu niedrige Temperaturen abfangen und
*     ! Variable fuer die Dissoziationskonstanten berechnen
*     =====================================================
      TT1 = Tg
      TT2 = TT1*TT1
      TT3 = TT2*TT1
      th1 = 5040.d0/TT1
      th2 = th1*th1
      th3 = th2*th1
      th4 = th3*th1
      kT  = bk*TT1
      kT1 = 1.d0/kT
*      
*-----------------------------------------------------------------------
*     ! init vectors
*     ==============
      anmono = 0.d0
      anmol  = 0.d0

* --------------------------------------------------------------------------
*     ! compute equilibrium constants
*     ===============================
      do i=1,nml
        if (i.ne.TiC) g(i)=gk(i)       ! compute all equil.constants lnk
      enddo  

*    TiC Gleichgewichtskonstante von Andreas Gauger ist anders
*        definiert als die Gleichgewichtskonstanten von Gail
*  Gauger: 
*  log Kp = 12.75293-5.4485*th1-1.56672*log(th1)+1.56041*(log(th1))**2
*           - 0.93275(log(th1))**3
*         = log ( p(A)p(B)/p(AB) )
*  Gail:
*   ln Kp = a0 + a1*th1 + a2*th1**2 + a3*th1**3 + a4*th1**4
*         =  ln ( p(AB)/p(A)p(B) )
*  Tsuji:
*  log Kp = a0 + a1*th1 + a2*th1**2 + a3*th1**3 + a4*th1**4
*         = log ( p(A)p(B)/p(AB) )
*
*  Umrechnung der Gauger-TiC-GG-Konstante in das Gail'sche System
*  -log(10)*log Kp(Gauger) = -2.30256*log Kp(Gauger) = ln Kp(Gail)

      lth = LOG10(th1)
      arg = 12.75293 - 5.44850*th1    - 1.56672*lth
     &               + 1.56041*lth**2 - 0.93275*lth**3
      g(TiC) = -2.30256*arg

*---------------------------------------------------------------------------
      !print'("smchem8 called ilauf,merk",i5,l2)',ilauf,merk
      !print'("T,n<H>,eps=",0pF10.2,99(1pE11.3))',Tg,anHges,eps
      NpreLoop = 0
      NpreIt = 0
      Ntaken = 0
      Nestim = 0
      DUALco = 0
      HCOco  = 0
      if ((ilauf.gt.10).and.merk) then
        do i=1,nel
          anmono(i) = amerk(i) * anhges
        enddo
        from_merk = .true.
        goto 200
      endif

*---------------------------------------------------------------------------
*     ! estimate electron density
*     =========================== 
 100  continue
      if (verbose>1) print'("SMCHEM_8:")'
      from_merk = .false.
      if (charge) then
        nelek = 0.d0 
        do i=1,nel
          if (i==el) cycle 
          ng = anHges * eps(i)
          Sa = EXP(g(elion(i)))*kT1
          nelek = nelek + ng/(0.5d0 + SQRT(0.25d0 + ng/Sa))
        enddo
        anmono(el) = nelek
        pel = nelek*kT
        mcharge = pel
        !peest = pel
        !pel = pecorr*pel
        !anmono(el) = pecorr*anmono(el) 
        if (verbose>1) print'(" estimate pel=",1pE10.3)',pel
      endif  

*-----------------------------------------------------------------------
*     ! estimate atomic pressures: new method
*     =======================================
      sortq(:) = 1.Q0
 110  Nseq = nel
      done(:) = .false.                ! all elements to be estimated here
      !if (charge) done(el)=.true.     ! ... except for the electrons      
      !if (charge) Nseq=nel-1
      eseq(:) = 0                      ! hirachical sequence of elements
      ptake = .true.
      do ido=1,Nseq
        !---------------------------------------------------------
        ! search for the most abundant element not yet considered 
        !---------------------------------------------------------
        emax = 0.d0
        enew = 0
        do e=1,nel
          if (done(e)) cycle
          norm(e) = eps(e)*sortq(e)
          if (e==el) norm(e)=anmono(el)*sortq(el)/anHges
          if (norm(e)<emax.or.(ido==1.and.e==el)) cycle
          emax = norm(e)
          enew = e
        enddo
        if (enew==0) then
          print*,catm 
          print*,eps
          print*,ido
          print*,done
          print*, "*** should not occur."
          return
        endif  
        if (verbose>1) print*,'estimate p'//trim(catm(enew))//' ...'
        if (enew.ne.pkey(ido)) ptake=.false.
        done(enew) = .true.
        eseq(ido) = enew               ! add to hirachical sequence 
        pges = eps(enew)*anHges*kT
        if (verbose>0) mols = ''
        !--------------------------------------------------------------------
        ! make rough estimate of atom pressure, considering that one molecule 
        ! takes all the element.  Store coeff for Sum_l coeff(l) p^l = pges 
        !--------------------------------------------------------------------
        coeff(:) = 0.d0   
        lnc(:) = 0.d0
        pwork = pges
        pwork2 = 0.d0
        imaj(enew) = 0
        imaj2(enew) = 0
        lmin = +99
        lmax = -99
        do e=1,nel
          if (.not.done(e)) cycle
          if (e==enew) cycle
          xx(e) = LOG(anmono(e)*kT)
        enddo  
        do i=1,nml
          affect = .false.
          known  = .true.
          lnp = g(i)
          do j=1,m_kind(0,i)
            e  = m_kind(j,i)
            nu = m_anz(j,i)
            if (.not.done(e)) then
              known = .false.
              exit
            endif  
            if (e==enew) then
              l = nu
              affect = .true.
            else
              lnp = lnp + nu*xx(e)
            endif  
          enddo  
          if (.not.affect) cycle  
          if (.not.known) cycle
          if (verbose>0) mols = trim(mols)//" "//cmol(i)
          lmin = MIN(lmin,l)
          lmax = MAX(lmax,l)
          coeff(l) = coeff(l) + l*EXP(lnp)
          if (l>0.and.(enew/=el)) then
            lnc(l) = MAX(lnc(l),LOG(REAL(l))+lnp)
            sum = LOG(pges) - lnp - LOG(REAL(l))
            if (sum>300.0*l) cycle
            ptest = EXP(sum/l)
            if (ptest<pwork) then
              pwork2 = pwork
              pwork = ptest
              imaj2(enew) = imaj(enew)
              imaj(enew) = i
              !print*,"first",cmol(i),pwork,lnc(l)
            else if (ptest<pwork2) then  
              pwork2 = ptest
              imaj2(enew) = i
              !print*,"second",cmol(i),pwork,lnc(l)
            endif
          endif
        enddo
        if (verbose>1) print*,trim(mols)
        if (enew==el) then
          !print*,coeff(-1),coeff(1)
          pel = SQRT(-coeff(-1)/(1.d0+coeff(+1)))    ! 0 = pel - a/pel + b*pel
          if (verbose>1) print*,'pel=',anmono(el)*kT,pel
          anmono(el) = pel*kT1
          mcharge = MAX(pel,MAX(coeff(-1)/pel,coeff(1)*pel))
        else   
          !---------------------------------------
          ! scale to make coeff(:) fit into real*8  
          !---------------------------------------
          aim = LOG(1.d+10)
          psc = -LOG(pges)                           ! when atom dominates
          do l=lmin,lmax
            if (lnc(l)==0.d0) cycle
            ptest = (aim-lnc(l))/l                   ! lnc(l)+l*ln(psc) = aim
            psc = MAX(psc,-ptest)
          enddo  
          psc = -psc
          if (verbose>1) print*,"pwork=",pwork,"  psc=",psc
          !------------------------------------------------
          ! store coeff for Sum_l coeff(l) (p*psc)^l = pges 
          !------------------------------------------------
          coeff(:) = 0.d0
          do l=lmin,lmax
            if (lnc(l)==0.d0) cycle
            coeff(l) = EXP(lnc(l)+l*psc)
          enddo
          psc = EXP(psc)
          !if (verbose>1) print'(" coeff(",I2,":",I2,")=",99(1pE10.2))',
     >    !               lmin,lmax,coeff(lmin:lmax)
          !----------------------------------------------
          ! solve 1d equation above with Newton's method 
          !----------------------------------------------
          pp = pwork/psc
          do piter=1,99                  
            f  = pp*psc-pges
            fs = psc
            do l=lmin,lmax
              if (coeff(l)==0.d0) cycle
              f  = f  + coeff(l)*pp**l
              fs = fs + coeff(l)*l*pp**(l-1)
            enddo
            if (fs==0.d0) print*, "*** fs=0 in smchem8 1d-pre-it."
            delta = f/fs
            pp = pp-delta
            if (verbose>1) print'(A2,I3,1pE25.15,1pE10.2)',
     >                     catm(enew),piter,pp*psc,delta/pp
            if (ABS(delta)<1.d-4*ABS(pp)) exit 
          enddo  
          if (piter>=99) then
            write(*,*) "*** smchem8 no conv in 1D pre-it "//catm(enew)
            write(*,*) anHges,Tg
            do i=1,nel
              write(*,'(A2,1x,0pF30.26)') catm(i),eps(i)
            enddo
            write(*,*) "psc=",real(psc),lmin,lmax
            write(*,*) "  lnc:",lnc(lmin:lmax)
            write(*,*) "coeff:",coeff(lmin:lmax)
            goto 1000
          endif  
          anmono(enew) = pp*psc*kT1
        endif  

        !-----------------------------------------------------------
        ! take into account feedback on elements considered before,
        ! unless they are much more abundant, with Newton-Raphson 
        !-----------------------------------------------------------
        nsave = anmono
 150    continue
        eact(:) = .false.
        Nact = 0
        !if (verbose>0) print*,"NewFastLevel,ptake=",NewFastLevel,ptake
        do iredo=MAX(1,ido-NewBackIt),ido
          e = eseq(iredo)
          if (norm(e)<NewBackFac*norm(enew)) then
            eact(e) = .true. 
            Nact = Nact+1
            all_to_act(e) = Nact
            act_to_all(Nact) = e
            pbefore(e) = anmono(e)
            if (NewFastLevel<3.and.ptake) then
              anmono(e) = anmono(e)*pcorr(enew,e)
              Ntaken = Ntaken+1
            else
              pcorr(enew,e) = 1.d0
            endif
            Nestim = Nestim+1
          endif
        enddo
        if (verbose>1) then
          print*,catm(eseq(1:ido))
          print*,eact(eseq(1:ido))
        endif
        if (Nact==1) then
          pkey(ido) = enew
          e = eseq(ido)
          xx(e) = LOG(anmono(e)*kT)
          cycle
        endif  
        do i=1,nml
          affect = .false. 
          known  = .true.
          do j=1,m_kind(0,i)
            e = m_kind(j,i) 
            if (.not.done(e)) then
              known = .false.
              exit
            endif
            if (eact(e)) affect=.true.
          enddo  
          relevant(i) = (known.and.affect)
        enddo

        !---------------------------------------------------------
        !***  check if new important molecule has huge impact  ***
        !---------------------------------------------------------
        i = imaj(enew)
        if (i>0) then
          pmol = g(i)
          do j=1,m_kind(0,i)
            pmol = pmol + m_anz(j,i)*LOG(anmono(m_kind(j,i))*kT)
          enddo
          pmol = EXP(pmol)
          i1 = 0
          s1 = 0
          qmost = 0.d0
          do j=1,m_kind(0,i)
            e = m_kind(j,i)
            if (e==el) cycle
            if (e==enew) then
              i1 = enew
              s1 = m_anz(j,i) 
              cycle
            endif
            qq = pmol*m_anz(j,i)/(eps(e)*anHges*kT)
            if (qq>qmost) then
              qmost = qq
              i2 = e
              s2 = m_anz(j,i)
            endif  
          enddo
          possible = (qmost>0.95)
          if (possible) then
            imaj2(i2) = imaj(i2)
            imaj(i2) = i
            pmol = eps(i2)*anHges*kT/s2
            if (eps(i1)*anHges*kT-s1*pmol<0.d0.and.imaj2(i2)>0) then
              tmp= i1
              i1 = i2
              i2 = tmp
              tmp= s1
              s1 = s2
              s2 = tmp
              pmol = eps(i2)*anHges*kT/s2
            endif
            if (eps(i1)*anHges*kT-s1*pmol<0.d0) possible=.false.
            if (i2==el) possible=.false.
          endif  
          if (possible) then
            !print*,catm(i1),cmol(i),cmol(imaj2(i1)),
     >      !       REAL(eps(i1)*anHges*kT-s1*pmol)
            !print*,catm(i2),cmol(i),cmol(imaj2(i2)),
     >      !       REAL(eps(i2)*anHges*kT-s2*pmol)
            i3 = imaj2(i1)      ! second most important carrier of i1
            s3 = 0
            s4 = 0
            lnp3 = g(i3)
            do j=1,m_kind(0,i3)
              e = m_kind(j,i3) 
              if (e==i1) then
                s3 = m_anz(j,i3)
                cycle
              endif  
              if (e==i2) then
                s4 = m_anz(j,i3)
                cycle
              endif  
              lnp3 = lnp3 + m_anz(j,i3)*LOG(anmono(e)*kT)
            enddo
            lnp = g(i)
            do j=1,m_kind(0,i)
              e = m_kind(j,i)
              if (e==i1.or.e==i2) cycle
              lnp = lnp + m_anz(j,i)*LOG(anmono(e)*kT)
            enddo
            if (verbose>1) print*,"readjusting "//catm(i1)//" "
     >                   //catm(i2)
     >                   //" with molecules "//trim(cmol(i))//" "
     >                   //trim(cmol(i3))//" ..."
            !--- solve eps(i1)*anHges*kT = s1*pmol + s3*p_second ---
            !--- with  eps(i2)*anHges*kT = s2*pmol               ---
            FF(1) = LOG((eps(i1)*anHges*kT-s1*pmol)/s3)-lnp3
            FF(2) = LOG((eps(i2)*anHges*kT)/s2)-lnp
            DF(1,1) = s3
            DF(1,2) = s4
            DF(2,1) = s1
            DF(2,2) = s2
            call GAUSS8(nel,2,DF,dp,FF)
            possible = (dp(1)<1000.0).and.(dp(2)<1000.0)
          endif
          if (possible) then
            anmono(i1) = EXP(dp(1))*kT1
            anmono(i2) = EXP(dp(2))*kT1
            !print'(A2,9(1pE16.8))',catm(i1),eps(i1)*anHges*kT,
     >      !                       s1*EXP(lnp+s1*dp(1)+s2*dp(2))
     >      !                      +s3*EXP(lnp3+s3*dp(1)+s4*dp(2))
            !print'(A2,9(1pE16.8))',catm(i2),eps(i2)*anHges*kT,
     >      !                       s2*EXP(lnp+s1*dp(1)+s2*dp(2))
            pbefore(i1) = anmono(i1)
            pbefore(i2) = anmono(i2)
            if (NewFastLevel<3.and.ptake) then
              anmono(i1) = anmono(i1)*pcorr(enew,i1)
              anmono(i2) = anmono(i2)*pcorr(enew,i2)
            endif  
            DUALco = DUALco+1
          endif  
          !--- special case CH4+CO2+H2O ---
          HCOproblem = (Nact>2).and.(enew==H.or.enew==C.or.enew==O)
     >                 .and.(CH4>0).and.(H2O>0).and.(CO2>0)
          if (HCOproblem) then
            hasH = .false.
            hasC = .false.
            hasO = .false.
            do ii=1,Nact
              i = act_to_all(ii)
              if (i==H) hasH=.true.
              if (i==C) hasC=.true.
              if (i==O) hasO=.true.
            enddo  
            HCOproblem = hasH.and.hasO.and.hasC
          endif  
          if (HCOproblem) then
            pH2O = (2*eps(O)-4*eps(C)+eps(H))*anHges*kT/4
            pCO2 = (2*eps(O)+4*eps(C)-eps(H))*anHges*kT/8
            pCH4 = (eps(H)+4*eps(C)-2*eps(O))*anHges*kT/8
            HCOproblem = (pH2O>0).and.(pCO2>0).and.(pCH4>0)
          endif  
          if (HCOproblem) then
            lpH2O = LOG(pH2O)
            lpCO2 = LOG(pCO2)
            lpCH4 = LOG(pCH4)
            lpH = (lpCH4+lpH2O*2+g(CO2)-lpCO2-g(CH4)-g(H2O)*2)/8
            lpO = (lpCO2+lpH2O*2+g(CH4)-lpCH4-g(CO2)-g(H2O)*2)/4
            lpC = (lpCH4+lpCO2+g(H2O)*2-lpH2O*2-g(CH4)-g(CO2))/2
            if (verbose>1) print*,"applying H2O-CH4-CO2 correction ..."
            !print'(2(1pE16.8))',2*pH2O+4*pCH4,anHges*kT*eps(H)
            !print'(2(1pE16.8))',  pH2O+2*pCO2,anHges*kT*eps(O)
            !print'(2(1pE16.8))',  pCH4+  pCO2,anHges*kT*eps(C)
            !print'(2(1pE16.8))',pH2O,EXP(lpH*2+lpO+g(H2O))
            !print'(2(1pE16.8))',pCO2,EXP(lpC+lpO*2+g(CO2))
            !print'(2(1pE16.8))',pCH4,EXP(lpC+lpH*4+g(CH4))
            anmono(H) = EXP(lpH)*kT1
            anmono(O) = EXP(lpO)*kT1
            anmono(C) = EXP(lpC)*kT1
            pbefore(H) = anmono(H)
            pbefore(O) = anmono(O)
            pbefore(C) = anmono(C)
            if (NewFastLevel<3.and.ptake) then
              anmono(H) = anmono(H)*pcorr(enew,H)
              anmono(O) = anmono(O)*pcorr(enew,O)
              anmono(C) = anmono(C)*pcorr(enew,C)
            endif  
            HCOco = HCOco + 1
          endif  
        endif
  
        !--- solve with different methods ---
        if (verbose>1) then 
          print'("corr",99(1pE11.2))',pcorr(enew,act_to_all(1:Nact))
        endif  
        qual  = 9.d+99
        do iloop=1,4
          if (iloop==1) then
            imethod = NewPreMethod                    ! first choice
          elseif (iloop==2) then
            imethod = NewPreMethod                    ! second try with pcorr=1.0
          elseif (iloop==3) then
            imethod = MOD(NewPreMethod,3)+1
          elseif (iloop==4) then
            imethod = MOD(NewPreMethod+1,3)+1
          endif
          if (verbose>1) print*,"imethod=",iloop,imethod
          !write(97,*) catm(eseq(1:ido))
          !write(97,*) eact(eseq(1:ido))
          !write(97,*) imethod
          !write(97,*) real(pcorr(enew,act_to_all(1:Nact)))          
          if (imethod==1) then
            !-------- method 1: xx=log(patm)-variables --------
            null = anmono
            qual0 = qual 
            qual1 = qual 
            qual2 = qual 
            do it=1,299
              do ii=1,Nact
                i = act_to_all(ii)
                xx(i) = LOG(anmono(i)*kT)
                xscale(i) = anmono(i)*kT
                if (i==el) then
                  fscale(ii) = 1.d0/mcharge
                else  
                  fscale(ii) = 1.d0/(anHges*eps(i)*kT)
                endif  
                FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
                DF(ii,:)  = 0.d0
                DF(ii,ii) = -fscale(ii)*anmono(i)*kT
              enddo  
              do i=1,nml
                if (.not.relevant(i)) cycle 
                pmol = g(i)
                do j=1,m_kind(0,i)
                  pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
                enddo
                pmol = EXP(pmol)
                do j=1,m_kind(0,i)
                  m1 = m_kind(j,i)
                  if (.not.eact(m1)) cycle
                  ii = all_to_act(m1)
                  term = fscale(ii)*m_anz(j,i)*pmol
                  FF(ii) = FF(ii) - term
                  do l=1,m_kind(0,i)
                    m2 = m_kind(l,i)
                    if (.not.eact(m2)) cycle
                    jj = all_to_act(m2)
                    DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
                  enddo	    
                enddo
              enddo
              call GAUSS8(nel,Nact,DF,dp,FF)
              bem = " "
              qual3 = qual2
              qual2 = qual1
              qual1 = qual0
              qual0 = qual
              qual  = 0.d0
              do ii=1,Nact
                qual = qual + ABS(dp(ii))
              enddo  
              maxs = 3.0
              if (it>30.and.(qual >qual0.or.qual0>qual1.or.
     >                       qual1>qual2.or.qual2>qual3)) then
                maxs = 3.0*exp(-MAX(0,it-30)/70.0)
                bem = "*"
              endif  
              do ii=1,Nact
                i = act_to_all(ii) 
                xx(i) = xx(i) - MAX(-maxs,MIN(maxs,dp(ii)))
                anmono(i) = exp(xx(i))*kT1
              enddo
              !write(97,'(I4,A2,99(1pE11.3E3))')
     >        !               it,bem,anmono(act_to_all(1:Nact))*kT,qual
              if (verbose>1) print'(I4,A2,99(1pE11.3E3))',
     >                       it,bem,anmono(act_to_all(1:Nact))*kT,qual
              NpreIt = NpreIt+1
              if (it>1.and.qual<1.d-10) exit
            enddo
            NpreLoop = NpreLoop+1

          else if (imethod==2) then
            !-------- method 2: lin-variables with pullback --------
            dp(:) = 0.d0
            qual  = 9.d+99
            fak   = 1.d0
            do it=1,199
              qual0 = qual 
              null  = anmono
              pullmax = 1
              if (it>100) pullmax=10
              do ipull=1,pullmax  ! pullback if quality gets worse
                !--- make a step ---
                do ii=1,Nact
                  i = act_to_all(ii)
                  anmono(i) = null(i)-fak*dp(ii)*kT1
                  xx(i) = LOG(anmono(i)*kT)
                enddo  
                !--- determine new FF and DF ---
                do ii=1,Nact
                  i = act_to_all(ii) 
                  xscale(i) = anmono(i)*kT
                  if (i==el) then
                    fscale(ii) = 1.d0/mcharge
                  else  
                    fscale(ii) = 1.d0/(anHges*eps(i)*kT)
                  endif  
                  FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
                  DF(ii,:) = 0.d0
                  DF(ii,ii) = -fscale(ii)*anmono(i)*kT
                enddo
                do i=1,nml
                  if (.not.relevant(i)) cycle 
                  pmol = g(i)
                  do j=1,m_kind(0,i)
                    pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
                  enddo
                  pmol = EXP(pmol)
                  do j=1,m_kind(0,i)
                    m1 = m_kind(j,i)
                    if (.not.eact(m1)) cycle
                    ii = all_to_act(m1)
                    term = fscale(ii)*m_anz(j,i)*pmol
                    FF(ii) = FF(ii) - term
                    do l=1,m_kind(0,i)
                      m2 = m_kind(l,i)
                      if (.not.eact(m2)) cycle
                      jj = all_to_act(m2)
                      DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
                    enddo	    
                  enddo
                enddo
                !--- determine new quality ---
                qual = 0.d0
                do ii=1,Nact
                  i = act_to_all(ii)
                  qual = qual + FF(ii)**2
                enddo  
                if (qual<qual0) exit
                if (ipull==pullmax) exit
                !write(97,'("pullback",3(1pE11.3))') fak,qual0,qual
                if (verbose>1) print'("pullback",3(1pE11.3))',
     >                         fak,qual0,qual
                fak = 0.5*fak   ! reduce NR-step
              enddo  
              !write(97,'(I4,99(1pE11.3))')
     >        !             it,anmono(act_to_all(1:Nact))*kT,qual
              if (verbose>1) print'(I4,99(1pE11.3))',
     >                     it,anmono(act_to_all(1:Nact))*kT,qual
              NpreIt = NpreIt+1
              if (it>1.and.qual<1.d-10) exit
              !--- determine new NR-vector ---
              call GAUSS8(nel,Nact,DF,dp,FF)
              do ii=1,Nact
                i = act_to_all(ii) 
                dp(ii) = dp(ii)*xscale(i)
              enddo
              null = anmono
              !--- limit step physically, keep direction ---
              fak = 1.d0
              do ii=1,Nact
                i = act_to_all(ii)
                if (null(i)*kT-fak*dp(ii)>5.d0*null(i)*kT) then
                  fak=MIN(fak,-4.d0*null(i)*kT/dp(ii))
                endif
                if (null(i)*kT-fak*dp(ii)<0.2d0*null(i)*kT) then
                  fak=MIN(fak,0.8d0*null(i)*kT/dp(ii))
                endif
              enddo
            enddo  
            NpreLoop = NpreLoop+1

          else if (imethod==3) then
            !-------- method 3: xx=log(patm)-variables --------
            null = anmono
            do it=1,299
              do ii=1,Nact
                i = act_to_all(ii)
                xx(i) = LOG(anmono(i)*kT)
                xscale(i) = anmono(i)*kT
                if (i==el) then
                  fscale(ii) = 1.d0/mcharge
                else  
                  fscale(ii) = 1.d0/(anHges*eps(i)*kT)
                endif  
                FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
                DF(ii,:)  = 0.d0
                DF(ii,ii) = -fscale(ii)*anmono(i)*kT
              enddo  
              do i=1,nml
                if (.not.relevant(i)) cycle 
                pmol = g(i)
                do j=1,m_kind(0,i)
                  pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
                enddo
                pmol = EXP(pmol)
                do j=1,m_kind(0,i)
                  m1 = m_kind(j,i)
                  if (.not.eact(m1)) cycle
                  ii = all_to_act(m1)
                  term = fscale(ii)*m_anz(j,i)*pmol
                  FF(ii) = FF(ii) - term
                  do l=1,m_kind(0,i)
                    m2 = m_kind(l,i)
                    if (.not.eact(m2)) cycle
                    jj = all_to_act(m2)
                    DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
                  enddo	    
                enddo
              enddo
              call GAUSS8(nel,Nact,DF,dp,FF)
              qual = 0.d0
              do ii=1,Nact
                qual = MAX(qual,ABS(dp(ii)))
              enddo  
              fak = MIN(1.0,3.0/qual)
              do ii=1,Nact
                i = act_to_all(ii) 
                xx(i) = xx(i) - fak*dp(ii)
                anmono(i) = exp(xx(i))*kT1
                if (IS_NAN(anmono(i))) then
                  qual=2.d+99
                  exit
                endif  
              enddo
              !write(97,'(I4,A2,99(1pE11.3E3))')
     >        !               it,bem,anmono(act_to_all(1:Nact))*kT,qual
              if (verbose>1) print'(I4,99(1pE11.3E3))',
     >                       it,anmono(act_to_all(1:Nact))*kT,fak,qual
              NpreIt = NpreIt+1
              if (qual<1.d-10.or.qual>1.d+99) exit
            enddo  
            NpreLoop = NpreLoop+1
          endif
          if (qual<1.d-5) exit
          do ii=1,Nact
            i = act_to_all(ii)
            anmono(i) = pbefore(i)    ! do not use pre-it corrections again
          enddo 
          if (verbose>1) read(*,'(A1)') char
        enddo  
        if (qual>1.d-4) then
          if (ptake) then
            anmono = nsave 
            ptake = .false.
            goto 150
          endif  
          write(*,*) "*** no convergence in NR pre-it "
          print*,"Tg=",Tg
          print*,catm(eseq(1:ido))
          print*,eact(eseq(1:ido))
          print*,pcorr(enew,act_to_all(1:Nact))
          goto 1000
        endif  
        !--- save ratio after/before for next run ---
        !print*,"ido,pkey(ido),enew=",ido,pkey(ido),enew
        pkey(ido) = enew
        pcorr(enew,:) = 1.d0
        do ii=1,Nact
          i = act_to_all(ii)
          pcorr(enew,i) = anmono(i)/pbefore(i)    
        enddo 
        if (verbose>1) print'("corr",99(1pE11.2))',
     >                 pcorr(enew,act_to_all(1:Nact))
        if (verbose>1) read(*,'(A1)') char
      enddo  
*
*     ! redo electron density
*     =======================
      if (charge) then
        coeff(:) = 0.d0
        relevant = .false.
        do i=1,nml
          do j=1,m_kind(0,i)
            if (m_kind(j,i)==el) relevant(i)=.true. 
          enddo  
          if (.not.relevant(i)) cycle
          pmol = g(i)
          l = 0
          do j=1,m_kind(0,i)
            pat = anmono(m_kind(j,i))*kT
            if (m_kind(j,i)==el) then
              l = m_anz(j,i)   
            else 
              pmol = pmol + m_anz(j,i)*LOG(pat)
            endif  
          enddo
          coeff(l) = coeff(l)+EXP(pmol)
        enddo
        pel = SQRT(coeff(-1)/(1.d0+coeff(+1)))     ! 0 = pel - a/pel + b*pel
        mcharge = MAX(pel,MAX(coeff(-1)/pel,coeff(1)*pel))
        if (verbose>1) print'(" redo pel =",1pE17.10," ->",1pE17.10)',
     >                        anmono(el)*kT,pel
        anmono(el) = pel*kT1
        !pecorr = pel/peest
      endif  

*     ! use memory of deviations between predicted atom pressures 
*     ! and converged atom pressures to improve the predictions
*     ============================================================
      ansave = anmono
      if (NewFastLevel<2.and.ptake) anmono = anmono*badness
*     
 200  continue
*-----------------------------------------------------------------------
      if (NewFullIt) then
*       ! Jacobi matrix and rhs vector for Newton-Raphson
*       =================================================
        it = 0
        eact(:) = .true.
        conv(:,:) = 9.d+99
        switchoff(:) = 0
        finish = 1.d-12
 300    continue
        if (it>30) finish=10.d0**(-12.0+3.0*(it-30.0)/(itmax-30.0))
        Nact = 0
        ii = 0
        do i=1,nel
          if (.not.eact(i)) cycle
          Nact = Nact+1
          ii = ii+1
          all_to_act(i) = ii
          act_to_all(ii) = i
          xx(i) = LOG(anmono(i)*kT)
          xscale(i) = anmono(i)*kT
          if (i==el) then
            fscale(ii) = 1.d0/mcharge
          else  
            fscale(ii) = 1.d0/(anHges*eps(i)*kT)
          endif  
          FF(ii) = fscale(ii)*(anHges*eps(i)*kT - anmono(i)*kT)
          DF(ii,:)  = 0.d0
          DF(ii,ii) = -fscale(ii)*anmono(i)*kT
        enddo	
        do i=1,nml
          pmol = g(i)
          do j=1,m_kind(0,i)
            pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
          enddo
          pmol = EXP(pmol)
          anmol(i) = pmol*kT1
          do j=1,m_kind(0,i)
            m1 = m_kind(j,i)
            if (.not.eact(m1)) cycle
            ii = all_to_act(m1)
            term = m_anz(j,i)*pmol*fscale(ii)
            FF(ii) = FF(ii) - term
            do l=1,m_kind(0,i)
              m2 = m_kind(l,i)
              if (.not.eact(m2)) cycle
              jj = all_to_act(m2)
              DF(ii,jj) = DF(ii,jj) - m_anz(l,i)*term
            enddo	    
          enddo
        enddo

*       ! compute the Newton-Naphson step	  
*       =================================
        FF0 = FF
        DF0 = DF
        !call QGEFA ( DF, nel, nel, ipvt, info )
        !call QGECO ( DF, nel, nel, ipvt, condnum1, work2 )
        !call QGESL ( DF, nel, nel, ipvt, FF, 0 )
        !dp  = FF
        !print'("condnum1 = ",1pE12.2E3)',condnum1
        call SGEIR(DF,nel,Nact,FF,1,ind,work,indx)
        dp = FF
        if (ind<0) then
          FF = FF0
          DF = DF0
          call GAUSS8(nel,Nact,DF,dp,FF)
        endif  
        !do i=1,Nact
        !  print'(99(1pE11.3E3))',DF0(i,1:Nact),dp(i),FF0(i)
        !enddo  
        !print*,ind

*       ! apply limited NR step and check convergence 	  
*       =============================================
        if (.true.) then        
          qual = 0.d0
          do ii=1,Nact
            qual = MAX(qual,ABS(dp(ii)))
          enddo  
          fak = MIN(1.0,3.0/qual)
          converge(it) = 0.d0
          Nconv = 0
          if (verbose>0) txt = ""
          do i=1,nel
            if (.not.eact(i)) then
              Nconv = Nconv+1 
              if (verbose>0) txt = trim(txt)//" "//catm(i)
            else 
              ii = all_to_act(i) 
              delp = -dp(ii)                             ! relative change dx/x
              conv(it,i) = delp
              converge(it) = MAX(converge(it),ABS(delp))
              if (ABS(delp)<finish) then
                Nconv = Nconv+1 
                if (verbose>0) txt = trim(txt)//" "//catm(i)
              endif  
            endif  
          enddo
          if (verbose>1.and.it<100) then
            print*,"ind,fak=",ind,fak
            print'(7x,A14,A14,A14)',"patom","dpatom","badness" 
            do ii=1,Nact
              i = act_to_all(ii) 
              print'(A7,3(1pE14.6))',catm(i),anmono(i)*kT,
     >             -dp(ii),badness(i)
            enddo
          endif
          do ii=1,Nact
            i = act_to_all(ii) 
            xx(i) = xx(i) - fak*dp(ii)
            anmono(i) = exp(xx(i))*kT1
          enddo
          limit = fak

        else
          !--- re-scale ---
          do ii=1,Nact
            i = act_to_all(ii) 
            dp(ii) = dp(ii)*xscale(i)
          enddo  
          fak = 5.d0
          limit = 1.d0                                   ! limit step, keep direction
          converge(it) = 0.d0
          Nconv = 0
          if (verbose>0) txt = ""
          do i=1,nel
            if (.not.eact(i)) then
              Nconv = Nconv+1 
              if (verbose>0) txt = trim(txt)//" "//catm(i)
            else 
              ii = all_to_act(i) 
              delp = -dp(ii)/(anmono(i)*kT)              ! relative change dx/x
              conv(it,i) = delp
              converge(it) = MAX(converge(it),ABS(delp))
              if (ABS(delp)<finish) then
                Nconv = Nconv+1 
                if (verbose>0) txt = trim(txt)//" "//catm(i)
              endif  
              if (1.d0+delp>fak) then
                limit = min(limit,(fak-1.d0)/delp)       ! such that xnew=xold*fac 
              else if (1.d0+delp<1.d0/fak) then
                limit = min(limit,(1.d0/fak-1.d0)/delp)  ! such that xnew=xold/fac
              endif
            endif  
          enddo
          if (it<=10) then
            limit = 1.d0
          else
            dp = dp*limit
          endif  
          if (verbose>1.and.it==0) then
            write(*,*) 
            print'(7x,A14,A14,A14)',"patom","dpatom","badness" 
            do ii=1,Nact
              i = act_to_all(ii) 
              print'(A7,3(1pE14.6))',catm(i),anmono(i)*kT,
     >             -dp(ii)/(anmono(i)*kT),badness(i)
            enddo
          endif

*         ! apply limited NR step
*         =======================
          !fak = 1.d0+4.d0*EXP(-(MAX(0,it-20))/13.d0)
          do ii=1,nact
            i = act_to_all(ii)
            delp = -dp(ii)*kT1
            nold = anmono(i)
            anmono(i) = MAX(nold/fak,MIN(nold*fak,nold+delp))
          enddo
          if (it>itmax-10) then
c            verbose = 2
            do ii=1,Nact
              i = act_to_all(ii) 
              delp = -dp(ii)/(anmono(i)*kT)
              print'(A3,99(1pE12.4))',catm(i),anmono(i),delp
            enddo  
          endif  

        endif
        crit = MAXVAL(converge(MAX(0,it-1):it))
        if (verbose>1) print'(3(i3),2(1pE9.1)," converged(",i2,"):",
     >                 A50)',it,Nact,ind,converge(it),limit,Nconv,txt
        if (it==itmax) then 
          write(*,*) '*** keine Konvergenz in SMCHEM8!'
          write(*,*) 'it, converge, ind =',it,converge(it),limit
          write(*,*) '  n<H>, T =',anhges,Tg
          write(*,*) 'from_merk,NewFastLevel=',from_merk,NewFastLevel
          if (ifatal<=1) then
            chemiter  = chemiter + it
            from_merk = .false.
            badness   = 1.Q0
            pcorr     = 1.Q0
            ifatal    = ifatal+1
c            if (ifatal==2) verbose=2
            goto 100        ! try again from scratch before giving up
          endif  
          goto 1000
        endif
        if (it>=5) then
          j = 0 
          do ii=1,Nact
            i = act_to_all(ii)
            if (MAXVAL(ABS(conv(it-5:it,i)))<finish) then
              switchoff(i) = it
              eact(i) = .false.
              j = j+1
              if (verbose>1) then
                print*,"switching off "//catm(i)//" ..."
              endif  
            endif
          enddo
          Nact = Nact-j
        endif  
        it = it + 1
        if (verbose.gt.1) read(*,'(a1)') char
        if (crit>finish.and.Nact>0) goto 300       ! continue iterating
*
*       ! redo rare elements
*       ====================
        redo(:) = .false.
        do iredo=1,nel
          atmax = 0.d0 
          e = 0
          do i=1,nel
            xx(i) = LOG(anmono(i)*kT)
            if (redo(i)) cycle   
            atfrac = anmono(i)/anHges
            if (atfrac>1.d-30) cycle   
            if (atfrac<atmax) cycle
            atmax = atfrac
            e = i
          enddo  
          if (e==0) exit
          redo(e) = .true.
          psc = xx(e)
          coeff(:) = 0.d0
          lmin = +99
          lmax = -99
          do i=1,nml
            pmol = g(i)
            affect = .false.
            do j=1,m_kind(0,i)
              nu = m_anz(j,i)   
              if (m_kind(j,i)==e) then
                l = nu
                affect = .true.
                pmol = pmol + l*psc
              else
                pmol = pmol + nu*xx(m_kind(j,i))
              endif  
            enddo
            if (.not.affect) cycle
            lmin = MIN(lmin,l)
            lmax = MAX(lmax,l)
            coeff(l) = coeff(l)+EXP(pmol)
          enddo
          psc = EXP(psc)
          if (verbose>1) then
            print'("redo rare element patom dp/patom ",A2)',catm(e)
            print*,"psc=",psc
            print'(" coeff(",I2,":",I2,")=",99(1pE10.2))',
     >                   lmin,lmax,coeff(lmin:lmax)
          endif  
          pat = anmono(e)*kT/psc
          do piter=1,99
            f  = pat*psc-eps(e)*anHges*kT
            fs = psc
            do l=lmin,lmax
              if (coeff(l)==0.d0) cycle
              f  = f  + coeff(l)*l*pat**l
              fs = fs + coeff(l)*l**2*pat**(l-1)
            enddo
            delta = f/fs
            if (verbose>1) print'(1x,A2,1pE25.15,1pE10.2)',
     >                            catm(e),pat*psc,delta/pat
            pat = pat-delta
            if (ABS(delta)<finish*ABS(pat)) exit 
          enddo  
          if (piter>=99) then
            write(*,*) "*** no convergence in post-it "//catm(e)
            write(*,*) anmono(e),eps(e)
            write(*,*) psc,lmin,lmax
            write(*,*) anHges,Tg
            write(*,*) coeff(lmin:lmax)
            goto 1000
          endif  
          anmono(e) = pat*psc*kT1
        enddo
*
*       ! how bad was the initial guess?
*       ================================
        if (.not.from_merk) then
          if (verbose>1) print'(7x,3(A14))',
     >                   "patom","conv.","init.guess"
          do i=1,nel
            badness(i) = anmono(i)/ansave(i)
            switch = switchoff(i)
            if (switch==0) switch=it-1
            if (verbose>1) then
              print'(A7,3(1pE14.6))',catm(i),anmono(i)*kT,
     >              conv(switch,i),badness(i)
            endif  
          enddo
!$omp critical(fort99)
          ilauf = ilauf+1
          !if (ilauf==1) write(99,'(A9,A10,A4,99(A10))') 
     >    !      'Tg','n<H>','it',catm(1:nel)
          !write(99,'(0pF9.3,1pE10.3,I4,99(1pE10.3))') 
     >    !      Tg,anHges,it,badness
!$omp end critical(fort99)
        endif  

      endif     ! NewFullIt

*     ! final anmol determination
*     ===========================
      do i=1,nel
        xx(i) = LOG(anmono(i)*kT)
      enddo  
      amerk = anmono/anHges
      anmol = 0.d0
      do i=1,nml
        pmol = g(i)
        do j=1,m_kind(0,i)
          pat = anmono(m_kind(j,i))*kT
          pmol = pmol + m_anz(j,i)*xx(m_kind(j,i))
        enddo
        anmol(i) = EXP(pmol)*kT1
      enddo
      if (charge) pel=anmono(el)*kT

      if (ngestst) then
*       ! Test auf Elementerhaltung
*       ===========================
        do i=1,nel
          nges(i) = anmono(i)
        enddo
        do i=1,nml
          do j=1,m_kind(0,i)
            j1 = m_kind(j,i)
            nges(j1) = nges(j1) + m_anz(j,i)*anmol(i)
          enddo
        enddo
        do e=1,nel
          if (e==el) cycle 
          soll  = anHges * eps(e)
          haben = nges(e)
          abw   = ABS(soll-haben)/MAX(soll,haben)
          if (abw>1.d-5) then
            if (verbose>1) then
              print'("*** element conservation error ",A2)',catm(e)
              print'(A12,1pE14.7)',catm(e),anmono(e)/(eps(e)*anHges)
            endif  
            sum = anmono(e)/(eps(e)*anHges)
            do i=1,nml
              do j=1,m_kind(0,i)
                j1 = m_kind(j,i)
                cc = m_anz(j,i)*anmol(i)/(eps(e)*anHges)
                if (j1==e.and.cc>1.d-7) then
                  if (verbose>1) print'(A12,1pE14.7)',cmol(i),cc
                  sum = sum+cc
                endif  
              enddo
            enddo
            if (verbose>1) then
              print'(3(1pE14.7))',soll/anHges,haben/anHges,sum
            endif  
            from_merk = .false.
            ansave = anmono
c            verbose=2
            goto 200
          endif
        enddo
      endif
      
      if (verbose.gt.0) print '("  ==> smchem used it=",I3,
     &                          " conv=",1pE9.2)',it,crit
      if (verbose.gt.0) print'("number of pre-iterations",I4,
     &     " -- used saved initial guesses",0pF5.1,"%")',
     &     NpreIt,REAL(Ntaken)/REAL(Nestim+1)*100
      if (verbose.gt.1) read(*,'(a1)') char

!$omp critical(counters)
      chemcall = chemcall + 1
      chemiter = chemiter + it
      preIter  = preIter  + NpreIt
      preUse   = preUse   + Ntaken
      preEst   = preEst   + Nestim
      DUALcorr = DUALcorr + DUALco
      HCOcorr  = HCOcorr  + HCOco
!$omp end critical(counters)

      return

 1000 continue
      open(unit=12,file='fatal.case')
      do i=1,nel
        write(12,'(A2,1x,0pF22.16)') catm(i),12.d0+log10(eps(i))
      enddo  
      write(12,*) anhges,Tg
      close(12)
      print*, "***  giving up."
      return


      CONTAINS       ! internal functions - not visible to other units 
************************************************************************
      REAL*8 FUNCTION gk(i)
************************************************************************
*****  returns  ln(kp) [cgs] for different fit formula             *****
************************************************************************
      use CHEMISTRY,ONLY: a,th1,th2,th3,th4,TT1,TT2,TT3,fit,natom,cmol,
     >                    NELEM,elnum,b_nasa,c_nasa
      implicit none
      integer,intent(in) :: i          ! index of molecule
      real*8,parameter :: bar=1.d+6, atm=1.013d+6, Rcal=1.987d+0
      real*8,parameter :: Rgas=8.3144598d+0
      real*8,parameter :: ln10=LOG(10.d0)
      real*8,parameter :: lnatm=LOG(atm), lnbar=LOG(bar)
      real*8 :: lnk,dG
      real*8 :: h_rt,s_r               !Added by Yui Kawashima
      real*8 :: dG_rt_ref(NELEM),dG_rt !Added by Yui Kawashima
      integer:: k,j                    !Added by Yui Kawashima

      if (i.eq.0) then
        gk = -1000.d0                  ! tiny kp for unassigned molecules
        return
      endif
      if (fit(i).eq.1) then
        !---------------------
        ! ***  Gail's fit  *** 
        !---------------------
        lnk = a(i,0) + a(i,1)*th1 + a(i,2)*th2 
     &               + a(i,3)*th3 + a(i,4)*th4 

      else if (fit(i).eq.2) then
        !---------------------------
        ! ***  Tsuji (1973) fit  *** 
        !---------------------------
        lnk = ln10*( - a(i,0) - a(i,1)*th1 - a(i,2)*th2
     &                        - a(i,3)*th3 - a(i,4)*th4 ) 

      else if (fit(i).eq.3) then  
        !---------------------------------
        ! ***  Sharp & Huebner (1990)  ***
        !---------------------------------
        dG  = a(i,0)/TT1 + a(i,1) + a(i,2)*TT1 + a(i,3)*TT2 + a(i,4)*TT3
        lnk = -dG/(Rcal*TT1) + (1-Natom(i))*lnatm

      else if (fit(i).eq.4) then
        !-----------------------------------
        ! ***  Stock (2008) & Kietzmann  ***
        !-----------------------------------
        dG  = a(i,0)/TT1+a(i,1)*LOG(TT1)+a(i,2)+a(i,3)*TT1+a(i,4)*TT2
        lnk = dG + (1-Natom(i))*lnbar

      else if (fit(i).eq.5) then
        !--------------------
        ! ***  dG(T)-fit  ***
        !--------------------
        dG  = a(i,0)/TT1 + a(i,1) + a(i,2)*TT1 + a(i,3)*TT2 + a(i,4)*TT3
        lnk = -dG/(Rgas*TT1) + (1-Natom(i))*lnbar
         
      else if (fit(i).eq.6) then
        !-------------------------------
        ! ***  Barklem & Collet fit  ***
        !-------------------------------
        lnk = a(i,0)/TT3 + a(i,1)/TT2 + a(i,2)/TT1 + a(i,3)/TT1**0.05d0
     &      + a(i,4)*LOG(TT1) + a(i,5) + a(i,6)*TT1 + a(i,7)*TT2

      else if (fit(i).eq.7) then
        !-----------------------------------------------------
        ! ***  NASA polynomial fit added by Yui Kawashima  ***
        !-----------------------------------------------------         
        if (TT1>1.d3) then
          h_rt = a(i,0) + a(i,1)*TT1/2.d0 
     &         + a(i,2)*TT1**2/3.d0 + a(i,3)*TT1**3/4.d0
     &         + a(i,4)*TT1**4/5.d0 + a(i,5)/TT1
          s_r  = a(i,0)*log(TT1) + a(i,1)*TT1
     &         + a(i,2)*TT1**2/2.d0 + a(i,3)*TT1**3/3.d0
     &         + a(i,4)*TT1**4/4.d0 + a(i,6)
        else
          h_rt = a(i,7) + a(i,8)*TT1/2.d0
     &         + a(i,9) *TT1**2/3.d0 + a(i,10)*TT1**3/4.d0
     &         + a(i,11)*TT1**4/5.d0 + a(i,12)/TT1
          s_r  = a(i,7)*log(TT1) + a(i,8)*TT1
     &         + a(i,9) *TT1**2/2.d0 + a(i,10)*TT1**3/3.d0
     &         + a(i,11)*TT1**4/4.d0 + a(i,13)           
        endif
        dG_rt = h_rt - s_r
        do k=1,m_kind(0,i)
          j = elnum(m_kind(k,i))
          if (c_nasa(j)==0) then
            print*,"Provide the data in data/Burcat_ref-elements.dat"
     &            ," and edit nasa_polynomial.f for "
     &            ,trim(catm(m_kind(k,i)))
            stop
          else
            if (TT1>1.d3) then
              h_rt = b_nasa(j,0) + b_nasa(j,1)*TT1/2.d0 
     &             + b_nasa(j,2)*TT1**2/3.d0
     &             + b_nasa(j,3)*TT1**3/4.d0
     &             + b_nasa(j,4)*TT1**4/5.d0 + b_nasa(j,5)/TT1
              s_r  = b_nasa(j,0)*log(TT1) + b_nasa(j,1)*TT1
     &             + b_nasa(j,2)*TT1**2/2.d0
     &             + b_nasa(j,3)*TT1**3/3.d0
     &             + b_nasa(j,4)*TT1**4/4.d0 + b_nasa(j,6)
            else
              h_rt = b_nasa(j,7) + b_nasa(j,8)*TT1/2.d0
     &             + b_nasa(j,9) *TT1**2/3.d0
     &             + b_nasa(j,10)*TT1**3/4.d0
     &             + b_nasa(j,11)*TT1**4/5.d0
     &             + b_nasa(j,12)/TT1
              s_r  = b_nasa(j,7)*log(TT1) + b_nasa(j,8)*TT1
     &             + b_nasa(j,9) *TT1**2/2.d0
     &             + b_nasa(j,10)*TT1**3/3.d0
     &             + b_nasa(j,11)*TT1**4/4.d0 + b_nasa(j,13)           
            endif
            dG_rt_ref(j) = h_rt - s_r
            dG_rt = dG_rt - m_anz(k,i)*dG_rt_ref(j)
          endif
        enddo
        dG = -dG_rt
        lnk = dG + (1-Natom(i))*lnbar

      else if (fit(i).eq.8) then
        !-----------------------------------------------
        ! ***  NASA 7-polynomial fit for BURCAT data ***
        !-----------------------------------------------
        if (TT1>1.d3) then
          dG = a(i,0) *(log(TT1)-1.d0) 
     &       + a(i,1) *TT1   / 2.d0
     &       + a(i,2) *TT1**2/ 6.d0    
     &       + a(i,3) *TT1**3/12.d0
     &       + a(i,4) *TT1**4/20.d0    
     &       - a(i,5) /TT1  
     &       + a(i,6)
        else
          dG = a(i,7) *(log(TT1)-1.d0) 
     &       + a(i,8) *TT1   / 2.d0
     &       + a(i,9) *TT1**2/ 6.d0    
     &       + a(i,10)*TT1**3/12.d0
     &       + a(i,11)*TT1**4/20.d0    
     &       - a(i,12)/TT1  
     &       + a(i,13)
        endif
        lnk = dG + (1-Natom(i))*lnbar

      else
        print*,cmol(i),"i,fit=",i,fit(i)
        stop "???"
      endif  

      gk = lnk
      end FUNCTION gk

      end SUBROUTINE smchem8
*********************************************************************
      SUBROUTINE NUCLEATION(species,T,V0,n1in,SSin,Jstar,Nstar)
*********************************************************************
*****                                                           *****
*****  computes nucleation rate according to                    *****
*****  classical nucleation theory (Gail et al 1984)            *****
*****                                                           *****
*****  INPUT:  species = name of nucleating species             *****
*****          T  = Gastemperatur [K]                           *****
*****          V0 = monomer volume [cm-3]                       *****
*****          n1 = monomer particle density [cm-3]             *****
*****          SS = supersaturation ratio                       ***** 
*****                                                           *****
*****  OUTPUT: Jstar = nucleation rate [cm^-3 s^-1]             *****
*****          Nstar = critical cluster size                    *****
*****                                                           *****
*********************************************************************
      use DUST_DATA,only: qp,mass,elnam
      use EXCHANGE,ONLY: NMOLE,nat,nmol,C,W
      use CHEMISTRY,ONLY: m_kind,m_anz,cmol,elnum
      implicit none
      character(len=*),intent(in) :: species
      real*8,intent(in) :: T,V0
      real(kind=qp),intent(in) :: n1in,SSin
      real*8,intent(out) :: Jstar,Nstar
      real*8 :: pi,bk,amu,a0,f0,sigma,Nf,alpha,slog
      real*8 :: thetun,thetaN,x0,x1,x2,x3,dgdn,nst,n1,SS
      real*8 :: zeldof,vth,beta,fNst,m,stoich
      integer :: i,j,iel,el
      logical :: found
      data pi/3.14159265358979D+0/, bk/1.38066D-16/, amu/1.66055D-24/

      n1 = n1in
      SS = SSin
      if (trim(species)=='C') then
        Nf    = 5.0                        ! fit for sigma=sigma(N)
        sigma = 1400.                      ! erg/cm2 (Gail+1984)
        iel   = C
      else if (trim(species)=='W') then
        Nf    = 10.0 
        sigma = 3340.                      ! erg/cm2 (R.Tran+2016)
        iel   = W
      else
        print*,"*** surface tension not known in NUCLEATION.f"
        stop
      endif  

*     -----------------------------------------
*     ***  monomer radius and surface area  *** 
*     -----------------------------------------
      a0 = (3.d0*V0/(4.d0*pi))**(1.d0/3.d0)
      f0 = 4.d0*pi*a0**2
      !print*,T,V0,a0,f0,n1,SS

*     -------------------------------
*     ***  supersaturation ratio  ***
*     -------------------------------
      if (SS.le.1.d0) then
        Jstar = 0.d+0
        Nstar = 9.d+99
        goto 500
      end if  
      slog = LOG(SS)

*     -------------------------------------------------------------
*     ***  size of critical cluster according to droplet model  ***
*     -------------------------------------------------------------
      thetun = f0*sigma/bk
      x0     = 2.d0*thetun/(3.d0*T*slog)
      x1     = Nf**(1.d0/3.d0)
      x2     = x1/x0
      x3     = 0.5d0*(1.d0+DSQRT(1.d0+2.d0*x2)) - x2
      Nstar  = 1.d0 + (x0*x3)**3 
      if (Nstar<=1.d0) Nstar=1.000000001d0
*
*     --------------------------------------------
*     ***  number density of critical cluster  ***
*     --------------------------------------------
      x0     = x0*slog
      x2     = (Nstar-1.d0)**(1.d0/3.d0) + x1
      x3     = 1.d0/(Nstar-1.d0)**(2.d0/3.d0)
      dgdn   = x0*x3* ( 1.d0/x2**2 + x1/x2**3 ) / 3.d0
      zeldof = SQRT(dgdn/(2.d0*pi))
      thetaN = thetun/(1.d0+(Nf/(Nstar-1.d0))**(1.d0/3.d0))
      x1     = (Nstar-1.d0)*slog - (thetaN/T)
     &         *(Nstar-1.d0)**(2.d0/3.d0)
      nst    = n1*EXP(x1)
      fNst   = f0*Nstar**(2.d0/3.d0)
*
*     -------------------------
*     ***  growth velocity  ***
*     -------------------------
      alpha = 1.d0
      vth   = SQRT(bk*T/(2.d0*pi*mass(iel)))
      beta  = alpha*nat(iel)*vth
      !print*,elnam(iel),nat(iel),mass(iel)/amu
      do i=1,NMOLE
        m = 0.d0
        stoich = 0.0
        found = .false.
        do j=1,m_kind(0,i)
          el = elnum(m_kind(j,i))
          if (el>0) then
            !print*,cmol(i),m_anz(j,i),elnam(el)
            m = m + m_anz(j,i)*mass(el)
          endif  
          if (iel==el) then
            stoich = m_anz(j,i)
            found = .true.
          endif
        enddo  
        if (found) then
          vth  = SQRT(bk*T/(2.d0*pi*m))
          beta = beta + alpha*nmol(i)*vth*stoich
          !print*,cmol(i),nmol(i),stoich,m/amu
        endif
      enddo  
*
*     -------------------------
*     ***  nucleation rate  ***
*     -------------------------
      Jstar = beta*nst*fNst*zeldof

      !print*,T,SS,Nstar,n1,nst/n1,Jstar,vth
      !if (Nstar<20) stop

 500  continue
      RETURN
      end
***********************************************************************
      SUBROUTINE GGCHEM(nHges,Tg,eps,merk,verbose)
***********************************************************************
*****                                                             *****
*****  Ruft lediglich SMCHEM auf (mit kompatiber Datenstruktur)   *****
*****                                                             *****
***********************************************************************
      use PARAMETERS,ONLY: Tfast
      use CHEMISTRY,ONLY: NMOLE,NELEM,NELM,elnum,elion,el,charge,catm
      use EXCHANGE,ONLY: nel,nat,nion,nmol
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: nHges,Tg
      real*8 :: epsi8(NELM),anmono8(NELM),nmol8(NMOLE)
      real(kind=qp),intent(in) :: eps(NELEM)
      real(kind=qp) :: epsi(NELM),anmono(NELM)
      logical,intent(in) :: merk
      integer,intent(in) :: verbose
      integer :: i,verb
      logical :: mk

      if (charge) epsi(el)=0.Q0
      do i=1,NELM
        if (i==el) cycle
        epsi(i) = eps(elnum(i))
        !print*,i,catm(i),elnum(i),epsi(i)
      enddo  
      verb = verbose
      mk = merk

      if (Tg<Tfast) then
        call SMCHEM16(nHges, Tg, epsi, anmono, nmol, mk, verb)
      else
        epsi8 = epsi
        call SMCHEM8 (nHges, Tg,epsi8,anmono8,nmol8, mk, verb)
        anmono = anmono8
        nmol   = nmol8
      endif  

      if (charge) nel=anmono(el)
      do i=1,NELM
        if (i==el) cycle                   ! slot for electrons
        nat(elnum(i))  = anmono(i)         ! neutral atoms
        if (charge) then
          nion(elnum(i)) = nmol(elion(i))  ! positive ions
        endif  
      enddo  
             
      RETURN
      end 




**********************************************************************
      SUBROUTINE INIT_DUSTCHEM
**********************************************************************
      use PARAMETERS,ONLY: model_eqcond,phyllosilicates
      use CHEMISTRY,ONLY: NMOLE,NELM,catm
      use DUST_DATA,ONLY: NDUSTmax,NEPS,NELEM,NDUST,eps0,amu,
     &                    dust_nam,dust_rho,dust_vol,dust_mass,
     &                    dust_nel,dust_nu,dust_el,fit,cfit,
     &                    Nfit,Tfit,Bfit,
     &                    elnr,elcode,elnam,mass,Tmelt,Tcorr,
     &                    DustChem_file
      use EXCHANGE,ONLY: H,Si,Al,Ca
      implicit none
      integer :: i,imax,j,k,el,j1,j2
      real*8 :: dmass,prec(NDUSTmax)
      character(len=10000) :: allcond
      character(len=200):: zeile,lastzeile
      character(len=100) :: trivial(NDUSTmax),tmp
      character(len=2)  :: name
      logical :: found,allfound,hasH,hasSi,hasAl,hasCa

      write(*,*) 
      write(*,*) "reading "//trim(DustChem_file)//" ..."
      write(*,*) "========================"
      trivial(:)=' '

c      open(12, file='data/'//trim(DustChem_file), status='old')
      open(12, file=trim(DustChem_file), status='old')
 
      write(*,*) '--- dust species ---'
      read(12,1000) zeile
      read(12,1000) zeile
      read(12,*) imax
      read(12,1000) zeile
      allcond = " "
      NDUST = 1
      do i=1,imax
        read(12,1000) zeile
        read(zeile,*) dust_nam(NDUST)
        j1 = index(zeile,' ')
        read(zeile(j1+1:),*) trivial(NDUST)
        if (index(zeile,'[l]')>0) then
          j2 = index(zeile,trim(trivial(NDUST)))
     &       + len(trim(trivial(NDUST)))
          read(zeile(j2+1:),*) Tmelt(NDUST)
          trivial(NDUST)=' '
        endif
        read(12,*) dust_rho(NDUST)
        read(12,*) dust_nel(NDUST)
        dmass = 0.d0
        allfound = .true.
        hasH  = .false.
        hasSi = .false.
        hasCa = .false.
        hasAl = .false.
        do j=1,dust_nel(NDUST)
          read(12,1030) dust_nu(NDUST,j),name
          found = .false. 
          do k=1,NELEM
            if (elnam(k).eq.name) then
              dust_el(NDUST,j) = k
              dmass = dmass + dust_nu(NDUST,j)*mass(k)
              found = .true.
              if (k==H)  hasH =.true.
              if (k==Si) hasSi=.true.
              if (k==Al) hasAl=.true.
              if (k==Ca) hasCa=.true.
            endif
          enddo
          if (.not.found) then
            print*,trim(dust_nam(NDUST)),name
            print*,elnam(1:NELEM)
            stop 'Element in dust species not found'
          endif  
          found = .false.
          do k=1,NELM
            if (catm(k).eq.name) then
              found = .true.
              exit
            endif
          enddo
          if (.not.found) allfound=.false.
        enddo
        found = .false.
        do 
          lastzeile = zeile 
          read(12,1000) zeile
          if (trim(zeile)=='') exit
          if (zeile(1:1)=='#') cycle
          read(zeile,*) fit(NDUST)
          prec(NDUST) = 0.0
          if (fit(NDUST).ne.7) then
            read(zeile,*) fit(NDUST),cfit(NDUST,0:4)
            j1 = index(lastzeile,'+/-')
            j2 = index(lastzeile,':')
            if (j1>0) then
              tmp = lastzeile(j1+3:)
              if (j2>j1) tmp=lastzeile(j1+3:j2-1)            
              read(tmp,*) prec(NDUST)
            endif  
            !print*,trim(tmp),prec(NDUST)
          else
            read(zeile,*) fit(NDUST),Nfit(NDUST),
     >                    Tfit(NDUST,1:Nfit(NDUST)+1)
            do k=1,Nfit(NDUST)
              read(12,*) Bfit(NDUST,k,1:14)
            enddo  
          endif  
          found = .true.
        enddo
        if (.not.found) then
          print*,"*** syntax error in DustChem.dat, condensate=",
     &         dust_nam(NDUST)
          stop
        endif  
        j1 = index(allcond," "//trim(dust_nam(NDUST)))
        if (j1>0) then
          print*,"*** double condensate in DustChem.dat"
          print*,dust_nam(NDUST)
          stop
        endif  
        if ((.not.phyllosilicates).and.hasH
     &      .and.(hasSi.or.hasAl.or.hasCa)) allfound=.false.
        if (allfound) then
          dust_mass(NDUST) = dmass
          dust_vol(NDUST) = dmass/dust_rho(NDUST)
          write(*,1060) NDUST,dust_nam(NDUST),dust_rho(NDUST),
     &                  dust_vol(NDUST), (dust_nu(NDUST,j),
     &                  elnam(dust_el(NDUST,j)),j=1,dust_nel(NDUST))
          allcond = " "//trim(allcond)//" "//trim(dust_nam(NDUST))
          NDUST = NDUST+1
        endif
      enddo
      NDUST=NDUST-1
      write(*,*) NDUST," condensed species"
      write(*,*)
      write(*,*) '--- involved elements ---'
      NEPS=0
      elcode(:)=0
      do i=1,NDUST
        do j=1,dust_nel(i)
          name = elnam(dust_el(i,j)) 
          do k=1,NELEM
            if (elnam(k).eq.name) then
              el = k
              exit
            endif
          enddo
          found = .false.           
          do k=1,NEPS
            if (el==elnr(k)) found=.true.
          enddo
          if (.not.found) then
            NEPS = NEPS+1 
            elnr(NEPS) = el
            elcode(el) = NEPS
            write(*,*) elcode(elnr(NEPS)),' ',name,el
          endif
        enddo
      enddo

      Tcorr(:) = -1.d0
      if (model_eqcond) call CHECK_MELTING
      write(*,*)
      
      !open(unit=1,file='condensates.tex')
      !do i=1,NDUST
      !  limit = ' '
      !  j = index(dust_nam(i),"[l]")
      !  if (Tcorr(i)>0.and.j>0) then
      !    write(limit,'("$>$",I4,"K")') int(Tcorr(i))
      !  else if (Tcorr(i)>0) then
      !    write(limit,'("$<$",I4,"K")') int(Tcorr(i))
      !  endif  
      !  if (prec(i)>0.0) then
      !    write(1,3000)
     &!      i,dust_nam(i),trivial(i),dust_rho(i),
     &!      fit(i),limit,cfit(i,0:4),prec(i)
      !  else  
      !    write(1,3001)
     &!      i,dust_nam(i),trivial(i),dust_rho(i),
     &!      fit(i),limit,cfit(i,0:4)
      !  endif  
      !enddo  
      !close(1)

      RETURN 
 1000 format(a200)
 1010 format(a2)
 1020 format(2(l1,1x),i1,1x,a10)
 1030 format(i2,1x,a2)
 1040 format(i2,1x,a10)
 1050 format(1x,a10,i4,' mass=',0pf7.3," amu")
 1060 format(I4,1x,a20," rhod=",0pf6.3," V0=",1pe11.3,2x,
     &       99(i2,"x",A2,1x))
 1070 format(1x,a10,99(i1,"x",i2,1x))
 2011 format(1(I2,1x,a8),22x,'->',I2,1x,a10,99(I2,1x,a8))
 2021 format(2(I2,1x,a8),11x,'->',I2,1x,a10,99(I2,1x,a8))
 2031 format(3(I2,1x,a8)    ,'->',I2,1x,a10,99(I2,1x,a8))
 3000 format(I3," & ",A20," & ",A25," & ",0pF5.2," & ",I2," & ",A8,
     &       " & ",5(1pE12.5," & "),"$\pm$",0pF4.2,"\\")
 3001 format(I3," & ",A20," & ",A25," & ",0pF5.2," & ",I2," & ",A8,
     &       " & ",5(1pE12.5," & "),9x,"\\")
      end 

***********************************************************************
      SUBROUTINE CHECK_MELTING
***********************************************************************
      use CHEMISTRY,ONLY: NMOLE,NELM,catm
      use DUST_DATA,ONLY: qp,NELEM,NDUST,dust_nam,Tmelt,Tcorr,is_liquid
      implicit none
      real*8 :: T
      real(kind=qp) :: nat(NELEM),nmol(NMOLE),Sat(NDUST)
      real(kind=qp) :: old,new,S(NDUST,10000)
      integer :: i,j,k,iT,Ncheck,il,is
      integer :: iliq(NDUST),isol(NDUST)
      character(len=15) :: search

      !--------------------------------------
      ! ***  identify solid/liquid pairs  ***
      !--------------------------------------
      is_liquid(:) = .false.
      Ncheck = 0
      do i=1,NDUST
        k = index(dust_nam(i),'[l]')
        if (k>0) then
          is_liquid(i) = .true. 
          Ncheck = Ncheck+1 
          iliq(Ncheck) = i
          isol(Ncheck) = 0
          search = dust_nam(i)
          search = search(1:k-1)//'[s]'
          do j=1,NDUST
            if (search==dust_nam(j)) then
              isol(Ncheck) = j
            endif
          enddo
          if (isol(Ncheck)==0) then
            print*,"*** liquid without solid "//trim(dust_nam(i))
            stop
          endif  
        endif
      enddo  
      if (Ncheck==0) return

      !-------------------------------
      ! ***  check melting points  ***
      !-------------------------------
      print*
      print*,'auto-correction for spurious liquid <-> solid '//
     &       'phase transitions ...'
      nat = 1.Q-100
      nmol = 1.Q-100
      do iT=100,10000
        T = DBLE(iT) 
        call SUPERSAT(T,nat,nmol,Sat)
        S(:,iT) = Sat(:)
      enddo  
      do i=1,Ncheck
        il = iliq(i)
        is = isol(i)
        do iT=101,10000
          T = DBLE(iT) 
          old = S(is,iT-1)/S(il,iT-1)
          new = S(is,iT)/S(il,iT)
          if (old>1.Q0.and.new<1.Q0) then
            !print'(A15,"-> ",A15,":",2(0pF8.1))',
     &      !     dust_nam(is),dust_nam(il),T,Tmelt(il)
          else if (old<1.Q0.and.new>1.Q0) then
            !print'(A15,"<- ",A15,":",0pF8.1,
     &      !     " false intersection point")',
     &      !     dust_nam(is),dust_nam(il),T
            if (T<Tmelt(il)) then
              Tcorr(il) = 0.5*(T+Tmelt(il))  
              print'(" ... correct ",A15," T <",0pF7.1)',
     &             dust_nam(il),Tcorr(il) 
            else  
              Tcorr(is) = 0.5*(T+Tmelt(il))  !correct solid
              print'(" ... correct ",A15," T >",0pF7.1)',
     &             dust_nam(is),Tcorr(is) 
            endif  
          endif  
        enddo   
      enddo
      do iT=100,10000
        T = DBLE(iT) 
        call SUPERSAT(T,nat,nmol,Sat)
        S(:,iT) = Sat(:)
      enddo  
      print'(26x,"melting point[K]  should be")'
      do i=1,Ncheck
        il = iliq(i)
        is = isol(i)
        do iT=101,10000
          T = DBLE(iT) 
          old = S(is,iT-1)/S(il,iT-1)
          new = S(is,iT)/S(il,iT)
          if (old>1.Q0.and.new<1.Q0) then
            print'(A15,"-> ",A15,":",2(0pF8.1))',
     &           dust_nam(is),dust_nam(il),T,Tmelt(il)
          else if (old<1.Q0.and.new>1.Q0) then
            print'(A15,"<- ",A15,":",0pF8.1,
     &           " false intersection point")',
     &           dust_nam(is),dust_nam(il),T
            stop
          endif  
        enddo   
      enddo
      end
************************************************************************
      subroutine INIT_CHEMISTRY
************************************************************************
      use PARAMETERS,ONLY: elements,initchem_info
      use CHEMISTRY,ONLY: NMOLdim,NMOLE,NELM,catm,cmol,el,
     &    dispol_file,source,fit,natom,a,error,i_nasa,
     &    m_kind,m_anz,elnum,elion,charge,Natmax,Ncmax,STOImax,
     &    el,H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,
     &    Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,W
      use DUST_DATA,ONLY: mass,mel,amu
      use EXCHANGE,ONLY: nmol,mmol
      implicit none
      integer :: loop,i,ii,j,iel,e,smax,ret
      character(len=2) :: cel(40),elnam
      character(len=20) :: molname,upper,leer='                    '
      character(len=200) :: filename
      character(len=300) :: line
      character(len=1) :: char
      logical :: found,charged
      real*8 :: fiterr

      cel(:) = '.'
      read(elements,*,end=100) cel
 100  NELM = 0
      charge = .false.
      do i=1,99
        if (cel(i)=='.') exit
        elnam = cel(i)
        found = .false.
        do e=1,NELM
          if (elnam.eq.catm(e)) then
            found=.true.
            exit
          endif  
        enddo
        NELM = NELM+1
        catm(NELM) = elnam
        if     (elnam=='el') then; el=NELM ; charge=.true.
        elseif (elnam=='H')  then;  H=NELM ; elnum(NELM)=1 
        elseif (elnam=='He') then; He=NELM ; elnum(NELM)=2 
        elseif (elnam=='Li') then; Li=NELM ; elnum(NELM)=3
        elseif (elnam=='Be') then; Be=NELM ; elnum(NELM)=4 
        elseif (elnam=='B')  then;  B=NELM ; elnum(NELM)=5 
        elseif (elnam=='C')  then;  C=NELM ; elnum(NELM)=6 
        elseif (elnam=='N')  then;  N=NELM ; elnum(NELM)=7 
        elseif (elnam=='O')  then;  O=NELM ; elnum(NELM)=8 
        elseif (elnam=='F')  then;  F=NELM ; elnum(NELM)=9 
        elseif (elnam=='Ne') then; Ne=NELM ; elnum(NELM)=10
        elseif (elnam=='Na') then; Na=NELM ; elnum(NELM)=11 
        elseif (elnam=='Mg') then; Mg=NELM ; elnum(NELM)=12
        elseif (elnam=='Al') then; Al=NELM ; elnum(NELM)=13
        elseif (elnam=='Si') then; Si=NELM ; elnum(NELM)=14
        elseif (elnam=='P')  then;  P=NELM ; elnum(NELM)=15 
        elseif (elnam=='S')  then;  S=NELM ; elnum(NELM)=16
        elseif (elnam=='Cl') then; Cl=NELM ; elnum(NELM)=17
        elseif (elnam=='Ar') then; Ar=NELM ; elnum(NELM)=18
        elseif (elnam=='K')  then;  K=NELM ; elnum(NELM)=19
        elseif (elnam=='Ca') then; Ca=NELM ; elnum(NELM)=20
        elseif (elnam=='Sc') then; Sc=NELM ; elnum(NELM)=21
        elseif (elnam=='Ti') then; Ti=NELM ; elnum(NELM)=22
        elseif (elnam=='V')  then;  V=NELM ; elnum(NELM)=23
        elseif (elnam=='Cr') then; Cr=NELM ; elnum(NELM)=24
        elseif (elnam=='Mn') then; Mn=NELM ; elnum(NELM)=25
        elseif (elnam=='Fe') then; Fe=NELM ; elnum(NELM)=26
        elseif (elnam=='Co') then; Co=NELM ; elnum(NELM)=27
        elseif (elnam=='Ni') then; Ni=NELM ; elnum(NELM)=28
        elseif (elnam=='Cu') then; Cu=NELM ; elnum(NELM)=29
        elseif (elnam=='Zn') then; Zn=NELM ; elnum(NELM)=30
        elseif (elnam=='Ga') then; Ga=NELM ; elnum(NELM)=31
        elseif (elnam=='Ge') then; Ge=NELM ; elnum(NELM)=32 
        elseif (elnam=='As') then; As=NELM ; elnum(NELM)=33 
        elseif (elnam=='Se') then; Se=NELM ; elnum(NELM)=34 
        elseif (elnam=='Br') then; Br=NELM ; elnum(NELM)=35 
        elseif (elnam=='Kr') then; Kr=NELM ; elnum(NELM)=36 
        elseif (elnam=='Rb') then; Rb=NELM ; elnum(NELM)=37 
        elseif (elnam=='Sr') then; Sr=NELM ; elnum(NELM)=38 
        elseif (elnam=='Y')  then;  Y=NELM ; elnum(NELM)=39 
        elseif (elnam=='Zr') then; Zr=NELM ; elnum(NELM)=40
        elseif (elnam=='W')  then;  W=NELM ; elnum(NELM)=41
        else
          stop "*** unknown element "
        endif
        if (initchem_info) then
          print*,'element '//elnam,elnum(NELM)
        endif
      enddo

      NMOLdim = 10000
      allocate(cmol(NMOLdim),fit(NMOLdim),natom(NMOLdim))
      allocate(error(NMOLdim),a(NMOLdim,0:13))
      allocate(source(NMOLdim),m_kind(0:6,NMOLdim),m_anz(6,NMOLdim))
      i_nasa = 0
      i = 1
      do loop=1,4
        filename = trim(dispol_file(loop))
        if (filename=='') exit
c        filename = 'data/'//trim(filename)
        if (initchem_info) write(*,*)
        write(*,*) 'reading kp-data from '
     &             //trim(filename)//" ..."
        open(unit=12, file=filename, status='old')
        read(12,*) NMOLdim
        do ii=1,NMOLdim
          read(12,'(A300)') line
          read(line,*) molname,iel,cel(1:iel),m_anz(1:iel,i)
          molname=trim(molname)
          fiterr = 0.0
          j = index(line,"+/-")
          if (j>0) read(line(j+3:),*) fiterr
          error(i) = fiterr
          read(12,'(A300)') line
          read(line,*) fit(i)
          if (fit(i)==6) then
            read(line,*) fit(i),(a(i,j),j=0,7)
          else if (fit(i)==7) then
            i_nasa = 1
            read(line,*) fit(i),(a(i,j),j=0,13)
          else if (fit(i)==8) then
            read(line,*) fit(i),(a(i,j),j=0,13)
          else   
            read(line,*) fit(i),(a(i,j),j=0,4)
          endif  
          m_kind(0,i) = iel
          natom(i) = 0
          found = .true.
          smax  = 0
          do j=1,m_kind(0,i)
            natom(i) = natom(i)+m_anz(j,i)
            if (index(elements,cel(j))<=0) found=.false. 
            smax = MAX(smax,ABS(m_anz(j,i)))
          enddo  
          if (.not.found) cycle    ! molecule has unselected element 
          if (smax>Natmax) cycle   ! stoichiometric coefficient > Natmax
          if (m_kind(0,i)==1.and.natom(i)==1) cycle  ! pure atom
          j = index(molname,"_")
          if (j>1) then
            cmol(i) = upper(molname(j+1:)//leer(1:j))
          else
            cmol(i) = upper(molname)
          endif
          charged = .false.
          do j=1,m_kind(0,i)
            elnam = cel(j)
            found = .false.
            do e=1,NELM
              if (elnam.eq.catm(e)) then
                found=.true.
                exit
              endif  
            enddo
            if (.not.found) stop "*** should not occur"
            m_kind(j,i) = e
            if (e==el) charged=.true. 
          enddo  
          if (fit(i)==6.and.charged) cycle ! old charged BarklemCollet
          source(i) = loop
          call CHECK_DOUBLE(cmol(i),m_kind(:,i),m_anz(:,i),i,loop,ret)
          if (ret>0) then
            source(ret) = loop
            cmol(ret) = cmol(i)
            fit(ret)  = fit(i)
            a(ret,:)  = a(i,:)
            error(ret)= error(i)
            if (initchem_info) then
              write(line,'(I4,A20,1x,99(I3,1x,A2,1x))')
     &           ret,trim(cmol(ret)),(m_anz(j,ret),cel(j),j=1,iel)
              print*,trim(line)//"    OVERWRITE"
            endif  
          else  
            if (initchem_info) then
              write(line,'(I4,A20,1x,99(I3,1x,A2,1x))')
     &            i,trim(cmol(i)),(m_anz(j,i),catm(m_kind(j,i)),j=1,iel)
              if (loop==1) then
                print*,trim(line)
              else
                print*,trim(line)//"    --> NEW" 
              endif
            endif  
            if (iel==2.and.
     >       ((m_kind(1,i)==el.and.m_anz(1,i)==-1.and.m_anz(2,i)==1).or.
     >        (m_kind(2,i)==el.and.m_anz(2,i)==-1.and.m_anz(1,i)==1))
     >        ) then
              e = m_kind(1,i)
              if (e==el) e=m_kind(2,i)
              elion(e) = i
            endif
            i = i+1
          endif  
        enddo
 200    close(12)
      enddo  
      NMOLE = i-1
      allocate(nmol(NMOLE),mmol(NMOLE))

      if (i_nasa==1) call NASA_POLYNOMIAL !Added by Yui Kawashima

      if (loop>1.and.initchem_info) then
        print* 
        do i=1,NMOLE
          print*,i,cmol(i),' ->  '//trim(dispol_file(source(i)))
        enddo
      endif  
  
      !open(unit=1,file='chemicals.tex')
      !write(1,*) NMOLE
      !do i=1,NMOLE
      !  if (error(i)>0.0) then 
      !    write(1,3000)
     &!      i,cmol(i),source(i),fit(i),a(i,0:4),error(i)
      !  else  
      !    write(1,3010)
     &!      i,cmol(i),source(i),fit(i),a(i,0:4)
      !  endif  
      !enddo  
      !close(1)
      !stop

      STOImax(:) = 0
      do i=1,NMOLE
        mmol(i) = 0.d0
        do j=1,m_kind(0,i)
          e = m_kind(j,i)
          if (e==el) then
            mmol(i) = mmol(i) + m_anz(j,i)*mel
          else
            mmol(i) = mmol(i) + m_anz(j,i)*mass(elnum(e))
          endif
          STOImax(e) = MAX(STOImax(e),m_anz(j,i))
        enddo
        !print*,cmol(i),mmol(i)/amu
      enddo  
      Ncmax = MAXVAL(STOImax)   ! maximum stoichiometric factor

      print* 
      print*,NMOLE,' species'
      print*,NELM,' elements'
      print'(99(A4))',(trim(catm(j)),j=1,NELM)
      print'(99(I4))',elnum(1:NELM)
      print'(99(I4))',STOImax(1:NELM)
      if (charge) then
        print'(1x,99(A4))',(trim(cmol(elion(j))),j=1,el-1),'  ',
     >                     (trim(cmol(elion(j))),j=el+1,NELM)
      endif  
      
 3000 format(I4," & ",A12," & (",I1,") & ",I1," & ",
     &       5(1pE12.5," & "),"$\pm$",0pF4.2,"\\")
 3010 format(I4," & ",A12," & (",I1,") & ",I1," & ",
     &       5(1pE12.5," & "),"\\")
      end

************************************************************************
      subroutine CHECK_DOUBLE(molname,kind,anz,N,loop,ret)
************************************************************************
      use PARAMETERS,ONLY: initchem_info
      use CHEMISTRY,ONLY: cmol,m_kind,m_anz,dispol_file,source
      implicit none
      character(len=20) :: molname
      integer,intent(IN) :: kind(0:6),anz(6),N,loop
      integer,intent(OUT) :: ret
      integer :: i,j,jj,el,ambi
      logical :: found,allfound,eqname,eqsource
      character(len=1) :: char1

      ret  = 0
      ambi = 0
      do i=1,N-1
        if (kind(0).ne.m_kind(0,i)) cycle   ! different no(elements)
        allfound=.true.
        do j=1,kind(0)
          el = kind(j)
          found = .false.
          do jj=1,m_kind(0,i)
            if (el.ne.m_kind(jj,i)) cycle
            found = .true.
            exit
          enddo
          if (.not.found) then
            allfound = .false.
            exit                            ! different elements
          else if (anz(j).ne.m_anz(jj,i)) then
            allfound = .false.
            exit                            ! different stoich.fac.
          endif
        enddo
        if (.not.allfound) cycle            
        eqname = (trim(molname)==trim(cmol(i)))
        eqsource = (loop==source(i))
        if (eqname.and.eqsource) then
          print*,"*** double molecule in "//dispol_file(loop)
          print*,trim(molname)//", "//trim(cmol(i))
          stop
        else if ((.not.eqname).and.eqsource.and.loop==1) then
          if (initchem_info) then
            print*,trim(molname)//", "//trim(cmol(i))//
     &           " different isomere in first source is OK"
          endif
          return  
        else if (eqname.and.(.not.eqsource)) then  
          ret = i
          return
        else
          ambi = i 
        endif
      enddo
      if (ambi>0) then
        if (source(ambi)==loop) then 
          if (initchem_info) then
            print*,trim(molname)//", "//trim(cmol(ambi))//
     &           " different isomere in subsequent source is OK"
          endif  
          ret = 0
          return
        else if (trim(dispol_file(loop))=="dispol_BURCAT.dat") then
          print*,"*** "//trim(molname)//", "//trim(cmol(ambi))//
     &         " ambiguous names in ..."
          print*,trim(dispol_file(loop))//
     &         ", "//trim(dispol_file(source(ambi)))
          print*,"will overwrite ..."
          ret = ambi
          !read(*,'(A1)') char1
        else  
          print*,"*** "//trim(molname)//", "//trim(cmol(ambi))//
     &         " ambiguous names in ..."
          print*,trim(dispol_file(loop))//
     &         ", "//trim(dispol_file(source(ambi)))
          print*,"please equalise in both data files."
          stop 
        endif  
      endif  
      end
**********************************************************************
      SUBROUTINE GAUSS_NM(Ndim,Mdim,N,M,A,x,b,info)
**********************************************************************
*****  tries to solve a system of N equations for M unknowns     *****
*****                   A x = b                                  *****
*****  info = 0  means that a solution was found                 *****
*****  info = 1  no solution in case N>M                         *****
*****  info = 2  NaNs produced (which normally occurs if N<M)    *****
**********************************************************************
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: Ndim,Mdim,N,M
      real(kind=qp),intent(inout) :: A(Ndim,Mdim),b(Ndim)
      real(kind=qp),intent(out) :: x(Mdim)
      integer,intent(out) :: info
      integer :: i,j,k,kmax,D
      real(kind=qp) :: c,Amax !,crit1,crit2

      !write(89,*) "GAUS_NM ..."

      D = min(N-1,M)

      do i=1,D

        !do k=1,N
        !  print'(99(1pE10.3))',A(k,1:M),b(k) 
        !enddo
        !print*

        !-------------------------------------
        !***  MAX-exchange of i-th column  ***      
        !-------------------------------------
        kmax = i
        Amax = ABS(A(i,i))
        do k=i+1,N
          if (ABS(A(k,i))>Amax) then
            kmax = k
            Amax = ABS(A(k,i))
          endif
        enddo  
        if (kmax.ne.i) then
          do j=1,M
            c         = A(i,j)
            A(i,j)    = A(kmax,j)
            A(kmax,j) = c 
          enddo  
          c       = b(i)
          b(i)    = b(kmax)
          b(kmax) = c 
        endif
        !-------------------------------
        !***  make triangular shape  ***
        !-------------------------------
        do k=i+1,N
          if (A(i,i)==0.Q0) then
            info = 2
            return
          endif  
          c = A(k,i) / A(i,i)
          A(k,i) = 0.Q0
          do j=i+1,M
            A(k,j) = A(k,j) - c * A(i,j)
          enddo  
          b(k) = b(k) - c * b(i)
        enddo  
      enddo  

      !do i=1,N
      !  print'(99(1pE10.3))',A(i,1:M),b(i) 
      !enddo
      !print*

      !-------------------
      !***  resolve x  ***
      !-------------------
      info = 0
      x(:) = 0.Q0
      !crit1 = 0.Q0
      !crit2 = 9.Q+99
      do i=M,1,-1
        if (A(i,i)==0.Q0) then
          info = 2
          return
        endif  
        c = 0.Q0
        do j=i+1,M
          c = c + A(i,j) * x(j)
        enddo  
        x(i) = (b(i) - c) / A(i,i)
        !crit1 = MAX(crit1,ABS(x(i)))
        if (ABS(x(i))>1.Q+25) info=2
      enddo  

      !print*,x(1:M)
      !print*

      if (N>M) then
        !--- more equations than unknowns --- 
        do i=M+1,N
          c = 0.Q0         
          do j=1,M
            c = c + A(i,j)*x(j)
          enddo
          !crit2 = MIN(crit2,ABS(c-b(i)))
          if (ABS(c-b(i))>1.Q-25) info=1 
        enddo 
      endif   

      !write(89,*) "GAUS_NM: crit1,crit2=",REAL(crit1),REAL(crit2)

      !if (info==0) then
      !  print*,'solution found.'
      !else if (info==1) then
      !  print*,'more eqs than unknowns, no linear combi'
      !else if (info==2) then
      !  print*,'more unknowns than eqs, no linear combi'
      !endif  

      end
**********************************************************************
           SUBROUTINE GAUSS8 (Ndim,N,a,x,b)
**********************************************************************
*****                                                            *****
*****   Diese Routine loesst ein lineares Gleichungssystem       *****
*****   der Form    (( a )) * ( x ) = ( b )     nach x auf.      *****
*****   Der Algorithmus funktioniert, indem die Matrix a         *****
*****   auf Dreiecksform gebracht wird.                          *****
*****                                                            *****
*****   EINGABE:  n = Dimension der Vektoren, der Matrix         *****
*****             a = (N x N)-Matrix                             *****
*****             b = (N)-Vektor                                 *****
*****   AUSGABE:  x = (N)-Vektor                                 *****
*****                                                            *****
**********************************************************************
      implicit none
      integer Ndim,N,i,j,k,kmax
      real*8  a(Ndim,Ndim),x(Ndim),b(Ndim),c,amax
*
c      integer ipiv(N),info      
c      call DGESV( N, 1, a(1,1), N, ipiv, b, N, info )
c      if (info.eq.0) then
c        x(:) = b(:)
c        RETURN
c      else  
c        write(*,*) 'linear equation system not solvable with DGESV'
c      endif

      do 500 i=1,N-1
*       ------------------------------------------
*       ***  MAX-Zeilentausch der i-ten Zeile  ***      
*       ------------------------------------------
        kmax = i
        amax = DABS(a(i,i))
        do 200 k=i+1,N
          if ( DABS(a(k,i)) .gt. amax ) then
            kmax = k
            amax = DABS(a(k,i))
          endif
  200   continue
        if (kmax.ne.i) then
          do 210 j=1,N
            c         = a(i,j)
            a(i,j)    = a(kmax,j)
            a(kmax,j) = c 
  210     continue
          c       = b(i)
          b(i)    = b(kmax)
          b(kmax) = c 
        endif
*
*       ---------------------------------
*       ***  bringe auf Dreiecksform  ***
*       ---------------------------------
        do 310 k=i+1,N
          c = a(k,i) / a(i,i)
          a(k,i) = 0.d0
          do 300 j=i+1,N
            a(k,j) = a(k,j) - c * a(i,j)
  300     continue        
          b(k) = b(k) - c * b(i)
  310   continue
*
  500 continue
*
*     --------------------------
*     ***  loese nach x auf  ***
*     --------------------------
      do 610 i=N,1,-1
        c = 0.d0
        if (i.lt.N) then
          do 600 j=i+1,N
            c = c + a(i,j) * x(j)
  600     continue
        end if
        x(i) = (b(i) - c) / a(i,i)
  610 continue
      RETURN
      end
**********************************************************************
           SUBROUTINE GAUSS16 (Ndim,N,a,x,b)
**********************************************************************
*****                                                            *****
*****   Diese Routine loesst ein lineares Gleichungssystem       *****
*****   der Form    (( a )) * ( x ) = ( b )     nach x auf.      *****
*****   Der Algorithmus funktioniert, indem die Matrix a         *****
*****   auf Dreiecksform gebracht wird.                          *****
*****                                                            *****
*****   EINGABE:  n = Dimension der Vektoren, der Matrix         *****
*****             a = (N x N)-Matrix                             *****
*****             b = (N)-Vektor                                 *****
*****   AUSGABE:  x = (N)-Vektor                                 *****
*****                                                            *****
**********************************************************************
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer :: Ndim,N,i,j,k,kmax
      real(kind=qp) :: a(Ndim,Ndim),x(Ndim),b(Ndim),c,amax
*
      do 500 i=1,N-1
*       ------------------------------------------
*       ***  MAX-Zeilentausch der i-ten Zeile  ***      
*       ------------------------------------------
        kmax = i
        amax = ABS(a(i,i))
        do 200 k=i+1,N
          if ( ABS(a(k,i)) .gt. amax ) then
            kmax = k
            amax = ABS(a(k,i))
          endif
  200   continue
        if (kmax.ne.i) then
          do 210 j=1,N
            c         = a(i,j)
            a(i,j)    = a(kmax,j)
            a(kmax,j) = c 
  210     continue
          c       = b(i)
          b(i)    = b(kmax)
          b(kmax) = c 
        endif
*
*       ---------------------------------
*       ***  bringe auf Dreiecksform  ***
*       ---------------------------------
        do 310 k=i+1,N
          c = a(k,i) / a(i,i)
          a(k,i) = 0.Q0
          do 300 j=i+1,N
            a(k,j) = a(k,j) - c * a(i,j)
  300     continue        
          b(k) = b(k) - c * b(i)
  310   continue
*
  500 continue
*
*     --------------------------
*     ***  loese nach x auf  ***
*     --------------------------
      do 610 i=N,1,-1
        c = 0.Q0
        if (i.lt.N) then
          do 600 j=i+1,N
            c = c + a(i,j) * x(j)
  600     continue
        end if
        x(i) = (b(i) - c) / a(i,i)
  610 continue
      RETURN
      end
      module CONVERSION
      use DUST_DATA,ONLY: NELEM,NDUSTmax
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer :: Nind,Ndep,Iindex(NELEM),Dindex(NDUSTmax+NELEM)
      logical :: is_dust(NDUSTmax+NELEM)
      real(kind=qp) :: conv(NDUSTmax+NELEM,NELEM)
!$omp threadprivate(Nind,Ndep,Iindex,Dindex,is_dust,conv)
      end

!-------------------------------------------------------------------------
      SUBROUTINE EQUIL_COND(nHtot,T,eps,Sat,ddust,verbose)
!-------------------------------------------------------------------------
! ***  computes the gas element abundances eps, the saturation ratios  ***
! ***  Sat and the dust number densities ddust after equilibrium       ***
! ***  condensation. The routine determines the state in which each    ***
! ***  solid is either undersaturated or saturated S<=1 (all solids).  ***
! ***  The book-keeping is special: Although in principle redundant,   ***
! ***  eps(NELEM) and ddust(NDUST) are changed by dx(NELEM) separately ***
! ***  to avoid numerical problems close to complete condensation      ***
!-------------------------------------------------------------------------
      use PARAMETERS,ONLY: Tfast,useDatabase,method_eqcond
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nam,dust_nel,dust_nu,dust_el,
     >                    eps0,elnam,elcode,NEPS,elnr
      use CHEMISTRY,ONLY: NewFastLevel,NELM,elnum,iel=>el
      use CONVERSION,ONLY: Nind,Ndep,Iindex,Dindex,is_dust,conv
      use EXCHANGE,ONLY: Fe,Mg,Si,Al,Ca,Ti,C,O,S,Na,Cl,H,Li,Mn,W,Ni,Cr,
     >                   Fluor=>F,Kalium=>K,Zr,V,
     >                   itransform,ieqcond,ieqconditer
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in) :: nHtot                ! H nuclei density [cm-3]
      real*8,intent(in) :: T                    ! temperature [K]
      real(kind=qp),intent(out) :: eps(NELEM)   ! gas element abundances
      real(kind=qp),intent(out) :: Sat(NDUST)   ! saturation ratio
      real(kind=qp),intent(out) :: ddust(NDUST) ! condensed units per <H> [-]
      integer,intent(inout) :: verbose
      real(kind=qp),dimension(NELEM) :: eps00,epsread,check,FF,Fsav,dx
      real(kind=qp),dimension(NELEM) :: eps_save,vec,xstep,Iabund,work
      real(kind=qp),dimension(NELEM) :: null,scale,bvec,LastRow
      real(kind=qp),dimension(NDUST) :: ddustread,dscale,pot,dust_save
      real(kind=qp),dimension(NDUST) :: Sat0,Sat1,Sat2,xvec,slin,tmp2
      real(kind=qp),dimension(NDUST) :: dstep
      real(kind=qp),dimension(NELEM,NELEM) :: DF,DFsav,emat,vecs,base_el
      real(kind=qp),dimension(NDUST,NELEM) :: mat
      real(kind=qp),dimension(NELEM,NDUST) :: AA,base
      real(kind=qp) :: worst,xmin,Smax,Smin,qual,qold,SQUAL,del
      real(kind=qp) :: abun,tmp1,amax
      real(kind=qp) :: turnon,maxon,minoff,fac,fac2,amount,Nt,dscalemax
      real(kind=qp) :: deps,esum,emax,NRstep,dev,LastDev,val1,val2,test
      real(kind=qp) :: det(2),converge(5000,NELEM),crit,cbest,qualDF
      real(kind=qp) :: deplete1,deplete2,small=1.Q-30
      real(kind=qp) :: dterm,dtmax,dlim,varfac,target,Qfinish,Sfinish
      integer,parameter :: itmax=5000
      integer,dimension(NELEM) :: elem,Nslot,eind
      integer,dimension(NDUST) :: dind,dlin,switchedON,switchedOFF
      integer,dimension(NELEM,NDUST) :: dustkind,stoich
      integer :: it,i,j,k,el,el2,Nact,Nact_read,Neq,slots,sl,dk,eq
      integer :: itry,knowns,unknowns,unknown,ii,jj,lastit,laston=0
      integer :: dtry,dtry_last=0,imaxon,iminoff,info,ipvt(NELEM)
      integer :: e_num(NELEM),e_num_save(NELEM)
      integer :: Nunsolved,unsolved(NELEM),Nvar1,Nvar2,var(NELEM)
      integer :: Nsolve,ebest,dbest,nonzero,itrivial,iread,ioff,method
      integer :: ifail,Nall,imax,swap,irow,erow,Eact,Nlin,iback,e1,e2
      integer :: act_to_elem(NELEM),act_to_dust(NELEM)
      integer :: Nzero,Ntrivial,etrivial(NELEM),dtrivial(NELEM)
      integer :: no_action,itJac,last_tri,method_failed
      logical,dimension(NELEM) :: e_resolved,e_act,e_taken,is_esolved
      logical,dimension(NELEM) :: e_eliminated,eblocked
      logical,dimension(0:NDUST) :: active,act_read,act_old
      logical,dimension(NDUST) :: is_dsolved,d_resolved,d_eliminated
      logical,dimension(NDUST) :: itried
      logical :: action,changed,solved,limited,ok,found
      logical :: limdust,dtry_break,IS_NAN
      character(len=1) :: char1,txt0
      character(len=2) :: rem,tnum
      character(len=6) :: dum6
      character(len=500) :: txt,txt1,txt2,text
      logical,save :: firstCall=.true.
      integer,save :: iAl2O3=0,iFe=0,iFeS=0,iNa2SiO3=0,iMgSiO3=0
      integer,save :: iMg2SiO4=0,iTi4O7=0,iCaSiO3=0,iCaMgSi2O6=0
      integer,save :: iNaAlSi3O8=0,iMgAl2O4=0,iCaTiO3=0,iSiO=0,iSiO2=0
      integer,save :: iTiO2=0,iMgTi2O5=0,iSiC=0,iCaS=0,iFe2SiO4=0,iFeO=0
      integer,save :: iNaCl=0,iKCl=0,iKAlSi3O8=0,iFe_l=0,iH2O=0,iH2O_l=0
      integer,save :: iFeS_l=0,iNaCl_l=0,iTiO2_l=0,iSiO2_l=0,iCaSO4=0
      integer,save :: iNa2SiO3_l=0,iMgAl2O4_l=0,iMg2SiO4_l=0,iMgSiO3_l=0
      integer,save :: iAl2O3_l=0,iCaAl2Si2O8=0,iC=0,iTiC=0,iFe2O3=0
      integer,save :: iMgO=0,iNa=0,iS=0,iMgS=0,iLiCl=0,iSiS2=0,iFeS2=0
      integer,save :: iH2SO4_l=0,iNa2S=0,iAlCl3=0,iNH3=0,iCaO=0,iNa_l=0
      integer,save :: iKCl_l=0,iCaCl2_l=0,iLiCl_l=0,iTi4O7_l=0,iFeO_l=0
      integer,save :: iCH4=0,iCaCl2=0,iLiH=0,iLiH_l=0,iMgTi2O5_l=0
      integer,save :: iMgO_l=0,iAlCl3_l=0,iMgTiO3=0,iMgTiO3_l=0,iCaO_l=0
      integer,save :: iS_l=0,iK2SiO3=0,iK2SiO3_l=0,iTiC_l=0,iTi=0
      integer,save :: iTi_l=0,iTiO=0,iTiO_l=0,iSiS2_l=0,iLiOH=0
      integer,save :: iLiOH_l=0,iMnS=0,iW=0,iW_l=0,iZrO2=0,iZrSiO4=0
      integer,save :: iVO=0,iV2O3=0,iCr=0,iH2SO4=0,iS2=0
      integer,save :: iNi=0,iNi_l,iNi3S2=0,iFe3O4=0,iKMg3AlSi3O12H2=0
      integer,save :: iKFe3AlSi3O12H2=0,iMg3Si2O9H4=0,iFe3Si2O9H4=0
      integer,save :: iMgCr2O4=0,iCr2O3=0,iMn3Al2Si3O12=0,iMn2SiO4=0
      integer,save :: iKAlSi2O6=0,iCa3Al2Si3O12=0,iFeAl2SiO7H2=0
      integer,save :: iNaMg3AlSi3O12H2=0,iNaAlSiO4=0,iCa2MgSi2O7=0
      integer,save :: iCa2Al2SiO7=0
      integer,save :: iCaTiSiO5=0,iNaAlSi2O6=0,iKAlSiO4=0,iMg3Si4O12H2=0
      integer,save :: iCa3MgSi2O8=0,iCaMgSiO4=0
      real*8 :: time0,time1,qread

      if (firstCall) then
        do i=1,NDUST
          if (dust_nam(i).eq.'Al2O3[s]')      iAl2O3=i 
          if (dust_nam(i).eq.'Al2O3[l]')      iAl2O3_l=i 
          if (dust_nam(i).eq.'Fe2O3[s]')      iFe2O3=i 
          if (dust_nam(i).eq.'SiO[s]')        iSiO=i 
          if (dust_nam(i).eq.'SiO2[s]')       iSiO2=i
          if (dust_nam(i).eq.'SiO2[l]')       iSiO2_l=i 
          if (dust_nam(i).eq.'SiS2[s]')       iSiS2=i
          if (dust_nam(i).eq.'SiS2[l]')       iSiS2_l=i
          if (dust_nam(i).eq.'Fe[s]')         iFe=i 
          if (dust_nam(i).eq.'Fe[l]')         iFe_l=i 
          if (dust_nam(i).eq.'FeS[s]')        iFeS=i 
          if (dust_nam(i).eq.'FeS[l]')        iFeS_l=i 
          if (dust_nam(i).eq.'FeS2[s]')       iFeS2=i 
          if (dust_nam(i).eq.'Na2SiO3[s]')    iNa2SiO3=i
          if (dust_nam(i).eq.'Na2SiO3[l]')    iNa2SiO3_l=i 
          if (dust_nam(i).eq.'MgSiO3[s]')     iMgSiO3=i 
          if (dust_nam(i).eq.'MgSiO3[l]')     iMgSiO3_l=i 
          if (dust_nam(i).eq.'Mg2SiO4[s]')    iMg2SiO4=i 
          if (dust_nam(i).eq.'Mg2SiO4[l]')    iMg2SiO4_l=i 
          if (dust_nam(i).eq.'CaSiO3[s]')     iCaSiO3=i 
          if (dust_nam(i).eq.'CaTiO3[s]')     iCaTiO3=i 
          if (dust_nam(i).eq.'CaMgSi2O6[s]')  iCaMgSi2O6=i
          if (dust_nam(i).eq.'CaAl2Si2O8[s]') iCaAl2Si2O8=i
          if (dust_nam(i).eq.'Ca2Al2SiO7[s]') iCa2Al2SiO7=i
          if (dust_nam(i).eq.'NaAlSi3O8[s]')  iNaAlSi3O8=i
          if (dust_nam(i).eq.'MgAl2O4[s]')    iMgAl2O4=i
          if (dust_nam(i).eq.'MgAl2O4[l]')    iMgAl2O4_l=i
          if (dust_nam(i).eq.'Ti4O7[s]')      iTi4O7=i
          if (dust_nam(i).eq.'Ti4O7[l]')      iTi4O7_l=i
          if (dust_nam(i).eq.'TiO2[s]')       iTiO2=i
          if (dust_nam(i).eq.'TiO2[l]')       iTiO2_l=i
          if (dust_nam(i).eq.'MgTi2O5[s]')    iMgTi2O5=i
          if (dust_nam(i).eq.'MgTi2O5[l]')    iMgTi2O5_l=i
          if (dust_nam(i).eq.'MgTiO3[s]')     iMgTiO3=i
          if (dust_nam(i).eq.'MgTiO3[l]')     iMgTiO3_l=i
          if (dust_nam(i).eq.'SiC[s]')        iSiC=i
          if (dust_nam(i).eq.'CaS[s]')        iCaS=i
          if (dust_nam(i).eq.'CaSO4[s]')      iCaSO4=i
          if (dust_nam(i).eq.'Fe2SiO4[s]')    iFe2SiO4=i
          if (dust_nam(i).eq.'Fe3O4[s]')      iFe3O4=i
          if (dust_nam(i).eq.'FeO[s]')        iFeO=i
          if (dust_nam(i).eq.'FeO[l]')        iFeO_l=i
          if (dust_nam(i).eq.'NaCl[s]')       iNaCl=i
          if (dust_nam(i).eq.'NaCl[l]')       iNaCl_l=i
          if (dust_nam(i).eq.'LiCl[s]')       iLiCl=i
          if (dust_nam(i).eq.'LiCl[l]')       iLiCl_l=i
          if (dust_nam(i).eq.'KCl[s]')        iKCl=i
          if (dust_nam(i).eq.'KCl[l]')        iKCl_l=i
          if (dust_nam(i).eq.'KAlSi3O8[s]')   iKAlSi3O8=i
          if (dust_nam(i).eq.'KAlSi2O6[s]')   iKAlSi2O6=i
          if (dust_nam(i).eq.'K2SiO3[s]')     iK2SiO3=i
          if (dust_nam(i).eq.'K2SiO3[l]')     iK2SiO3_l=i
          if (dust_nam(i).eq.'H2O[s]')        iH2O=i
          if (dust_nam(i).eq.'H2O[l]')        iH2O_l=i
          if (dust_nam(i).eq.'H2SO4[l]')      iH2SO4_l=i
          if (dust_nam(i).eq.'C[s]')          iC=i
          if (dust_nam(i).eq.'S[s]')          iS=i
          if (dust_nam(i).eq.'S[l]')          iS_l=i
          if (dust_nam(i).eq.'Ti[s]')         iTi=i
          if (dust_nam(i).eq.'Ti[l]')         iTi_l=i
          if (dust_nam(i).eq.'TiC[s]')        iTiC=i
          if (dust_nam(i).eq.'TiC[l]')        iTiC_l=i
          if (dust_nam(i).eq.'MgO[s]')        iMgO=i
          if (dust_nam(i).eq.'MgO[l]')        iMgO_l=i
          if (dust_nam(i).eq.'CaO[s]')        iCaO=i
          if (dust_nam(i).eq.'CaO[l]')        iCaO_l=i
          if (dust_nam(i).eq.'MgS[s]')        iMgS=i
          if (dust_nam(i).eq.'Na[s]')         iNa=i 
          if (dust_nam(i).eq.'Na[l]')         iNa_l=i 
          if (dust_nam(i).eq.'Na2S[s]')       iNa2S=i
          if (dust_nam(i).eq.'AlCl3[s]')      iAlCl3=i 
          if (dust_nam(i).eq.'AlCl3[l]')      iAlCl3_l=i 
          if (dust_nam(i).eq.'NH3[s]')        iNH3=i
          if (dust_nam(i).eq.'CaCl2[s]')      iCaCl2=i
          if (dust_nam(i).eq.'CaCl2[l]')      iCaCl2_l=i
          if (dust_nam(i).eq.'CH4[s]')        iCH4=i
          if (dust_nam(i).eq.'LiH[s]')        iLiH=i
          if (dust_nam(i).eq.'LiH[l]')        iLiH_l=i
          if (dust_nam(i).eq.'TiO[s]')        iTiO=i
          if (dust_nam(i).eq.'TiO[l]')        iTiO_l=i
          if (dust_nam(i).eq.'LiOH[s]')       iLiOH=i
          if (dust_nam(i).eq.'LiOH[l]')       iLiOH_l=i
          if (dust_nam(i).eq.'MnS[s]')        iMnS=i 
          if (dust_nam(i).eq.'W[s]')          iW=i 
          if (dust_nam(i).eq.'W[l]')          iW_l=i 
          if (dust_nam(i).eq.'Ni[s]')         iNi=i 
          if (dust_nam(i).eq.'Ni[l]')         iNi_l=i 
          if (dust_nam(i).eq.'ZrO2[s]')       iZrO2=i 
          if (dust_nam(i).eq.'ZrSiO4[s]')     iZrSiO4=i 
          if (dust_nam(i).eq.'KAlSiO4[s]')    iKAlSiO4=i
          if (dust_nam(i).eq.'Ni3S2[s]')      iNi3S2=i
          if (dust_nam(i).eq.'Mg3Si2O9H4[s]') iMg3Si2O9H4=i
          if (dust_nam(i).eq.'Mg3Si4O12H2[s]') iMg3Si4O12H2=i
          if (dust_nam(i).eq.'Fe3Si2O9H4[s]') iFe3Si2O9H4=i
          if (dust_nam(i).eq.'MgCr2O4[s]')    iMgCr2O4=i
          if (dust_nam(i).eq.'Cr[s]')         iCr=i
          if (dust_nam(i).eq.'Cr2O3[s]')      iCr2O3=i
          if (dust_nam(i).eq.'Mn2SiO4[s]')    iMn2SiO4=i
          if (dust_nam(i).eq.'NaAlSi2O6[s]')  iNaAlSi2O6=i
          if (dust_nam(i).eq.'NaAlSiO4[s]')   iNaAlSiO4=i
          if (dust_nam(i).eq.'CaMgSiO4[s]')   iCaMgSiO4=i
          if (dust_nam(i).eq.'CaTiSiO5[s]')   iCaTiSiO5=i
          if (dust_nam(i).eq.'Ca2MgSi2O7[s]') iCa2MgSi2O7=i
          if (dust_nam(i).eq.'Ca3MgSi2O8[s]') iCa3MgSi2O8=i
          if (dust_nam(i).eq.'Ca3Al2Si3O12[s]') iCa3Al2Si3O12=i
          if (dust_nam(i).eq.'Mn3Al2Si3O12[s]') iMn3Al2Si3O12=i
          if (dust_nam(i).eq.'FeAl2SiO7H2[s]')  iFeAl2SiO7H2=i
          if (dust_nam(i).eq.'KMg3AlSi3O12H2[s]') iKMg3AlSi3O12H2=i
          if (dust_nam(i).eq.'KFe3AlSi3O12H2[s]') iKFe3AlSi3O12H2=i
          if (dust_nam(i).eq.'NaMg3AlSi3O12H2[s]') iNaMg3AlSi3O12H2=i
          if (dust_nam(i).eq.'VO[s]')         iVO=i
          if (dust_nam(i).eq.'V2O3[s]')       iV2O3=i
          if (dust_nam(i).eq.'H2SO4[s]')      iH2SO4=i
          if (dust_nam(i).eq.'S2[s]')         iS2=i
        enddo
        firstCall = .false. 
      endif

      open(unit=1,file='Last_abund.in',status='replace')
      do i=1,NEPS
        el = elnr(i)
        write(1,'(A2,2x,0pF25.20)') elnam(el),12.0+LOG10(eps0(el))
      enddo
      write(1,*) nHtot,T
      close(1)

      if (verbose>=0) then
        write(*,*)
        write(*,'("EQUIL_COND started")')
      endif  
      call CPU_TIME(time0)

      !------------------------
      ! ***  initial state  ***
      !------------------------
      ddust  = 0.Q0                 ! initial state dust-free
      eps    = eps0                 ! initial gas abundances 
      active = .false.              ! no solid condensing

      !--------------------------------------------
      ! ***  load initial state from database?  ***
      !--------------------------------------------
      Nact = 0
      Nact_read = 0
      act_read = .false.
      if (useDatabase) then
        call GET_DATA(nHtot,T,epsread,ddustread,qread,iread,
     >                act_read,method)
        if (qread.lt.1.0) then
          eps    = epsread
          ddust  = ddustread
          active = act_read
          text = "active solids:"
          Nact_read = 0
          do i=1,NDUST
            if (.not.act_read(i)) cycle
            Nact_read = Nact_read + 1
            text = trim(text)//" "//trim(dust_nam(i))
          enddo
          Nact = Nact_read
          !verbose = 0
          !if (qread>1.Q-3.and.Nact>0) verbose=2
          !if (qread>1.Q-3.and.iread==207) verbose=2
          !if (method==2) verbose=2
          if (verbose>0) then
            write(*,'(" ... using database entry (",I6,
     >            ") qual=",1pE15.7)') iread,qread
            write(*,*) trim(text)
          endif  
        endif  
      endif
  
      !----------------------------------------------------
      ! ***  recompute eps00 from initial state,        ***
      ! ***  because these can very slightly drift      ***
      ! ***  eps00 = total element abundances gas+dust  ***
      !----------------------------------------------------
      check = eps
      do i=1,NDUST
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          check(el) = check(el) + ddust(i)*dust_nu(i,j)    
        enddo
      enddo
      worst = 0.Q0
      do i=1,NEPS
        el = elnr(i)
        worst = MAX(worst,ABS(1.Q0-check(el)/eps0(el)))
      enddo
      !eps00 = check
      eps00 = eps0
      if (verbose>0) then
        write(*,*) "element conservation error:",worst
        write(*,*) "initial gas fractions ..."
        do i=1,NELEM
          if (elcode(i)==0) cycle
          print'(3x,A2,2(1pE15.6))',elnam(i),eps(i),eps(i)/eps0(i)
        enddo
      endif  
      if (worst>1.Q-8) print*, "*** worst>1.Q-8 in equil_cond"

      !----------------------------------------------------------
      ! ***  compute maximum possible dust abundances dscale  ***
      ! ***  if all elements turn into selected dust species  ***
      !----------------------------------------------------------
      do i=1,NDUST
        xmin = 9.Q+99 
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          xmin = min(xmin,eps00(el)/REAL(dust_nu(i,j),kind=qp))    
        enddo
        dscale(i) = xmin                         ! max dust abundances
      enddo   
      dscalemax = MAXVAL(dscale(1:NDUST))

      null(:) = 0.Q0             
      call SUPER(nHtot,T,null,eps,Sat0,.false.) ! from scratch
      qual = SQUAL(Sat0,active)
      Smax = maxval(Sat0)
      if (verbose>=0) then
        print'("it =",I4,"  qual =",1pE11.4,"  Smax-1 =",1pE11.2E4)',
     >          0,qual,Smax-1.Q0
      endif  
      if (T>Tfast) then
        Qfinish = 1.Q-11
        Sfinish = 1.Q-09
      else
        Qfinish = 1.Q-18
        Sfinish = 1.Q-15
      endif
      switchedON(:) = 0
      switchedOFF(:) = 0
      if (verbose>-1) write(97,*)
      act_old = active
      lastit = -99
      iminoff = 0
      limited = .false.
      ifail = 0
      no_action = 0
      
!---------------------------  start of main iteration loop  ----------------------------

      do it=1,itmax

        !---------------------------------------
        ! ***  selection of solids to solve  ***
        !---------------------------------------
        changed = .false.
        Smax = maxval(Sat0)
        ioff = 0
        if ((qread<0.5).and.(it<=3).and.(Nact_read>0)
     >                             .and.(qual>0.Q0)) then
          active = act_read
          Nact = Nact_read
          pot = 0.0
        else if (it>lastit+3) then
          maxon   = 0.Q0 
          minoff  = 0.Q0 
          imaxon  = 0
          if (qual==0.Q0) limited=.false.
          do i=1,NDUST
            xmin = 9.Q+99 
            esum = 0.Q0
            emax = 0.Q0
            do j=1,dust_nel(i)
              el = dust_el(i,j)
              esum = esum + dust_nu(i,j)
              if (eps(el)/REAL(dust_nu(i,j),kind=qp).lt.xmin) then    
                xmin = eps(el)/REAL(dust_nu(i,j),kind=qp)
                emax = REAL(dust_nu(i,j),kind=qp)
              endif
            enddo
            esum = esum**1.65            ! simple condensates first
            abun = dscalemax/dscale(i)   ! potentially abundant condensates first
            abun = MIN(abun,1.Q+10)
            pot(i)  = 1.Q0/(esum+0.05*abun+1.0*switchedOFF(i))
            Sat1(i) = Sat0(i)**pot(i)
            if (verbose>1.and.(.not.active(i)).and.Sat0(i)>1.Q0) then
              print'(A20," simplicity=",1pE10.3," abun=",1pE10.3,
     >                   " pot=",3(1pE10.3))',
     >           dust_nam(i),1.Q0/esum,1.Q0/abun,pot(i),Sat0(i),Sat1(i)
            endif  
          enddo 
          Smax = 0.Q0
          imax = 0
          do i=1,NDUST
            if (Sat1(i)>Smax) then
              Smax = Sat1(i)
              imax = i
            endif  
            if (Sat1(i)>1.Q0.and.(.not.active(i))) then
              turnon = Sat1(i)-1.Q0 
              if (turnon>maxon.and.(.not.limited)) then
                maxon  = turnon
                imaxon = i
              endif  
            endif  
          enddo
          if (qual>maxon.and.no_action<50) then  ! keep on iterating without switching on yet
            maxon = 0.0
            imaxon = 0
          endif  
          if (verbose>0) print'("limited=",L1,
     >                   "  Smax=",1pE10.3,2x,A18)',
     >                   limited,Smax,dust_nam(imax)
          if (verbose>0.and.maxon>0.Q0) print'("  maxon =",
     >                   1pE10.2,2x,A18)',maxon,dust_nam(imaxon)

          if (maxon>0.0*MAX(Smax-1.Q0,0.Q0)) then
            if (imaxon.ne.iminoff) then
              txt = dust_nam(imaxon)
              i = index(txt,'[l]')
              if (i>0) then
                txt = trim(txt(1:i-1))//'[s]'
                do i=1,NDUST
                  if (active(i).and.trim(dust_nam(i)).eq.trim(txt)) then
                    !print*,dust_nam(imaxon),Sat1(imaxon)
                    !print*,dust_nam(i),Sat1(i)
                    if (Sat1(imaxon)<Sat1(i)) imaxon=0
                  endif
                enddo  
              endif  
            endif
            if (imaxon>0) then
              active(imaxon) = .true.
              if (verbose>=0) then
                print*,"switch on ",trim(dust_nam(imaxon))
              endif  
            endif  
          endif  
          Nact = 0
          do i=1,NDUST
            active(i) = active(i).or.(ddust(i)>0.Q0)
            if (active(i).neqv.act_old(i)) changed=.true.
            if (active(i)) Nact=Nact+1
          enddo
          
          if (imaxon>0.and.Nact>1) then
            !----------------------------------------
            ! ***  eliminate linear combinations  ***
            !----------------------------------------
            eps_save = eps
            dust_save = ddust
            ioff = 0
            e_act(:) = .false.
            do i=1,NDUST
              if (.not.active(i)) cycle 
              do j=1,dust_nel(i)
                el = dust_el(i,j)
                e_act(el) = .true.
              enddo
            enddo  
            eind(:) = 0
            erow = 0
            do el=1,NELEM
              if (.not.e_act(el)) cycle
              erow = erow+1
              eind(el) = erow
            enddo  
            Eact = erow
            AA(:,:) = 0.Q0
            bvec(:) = 0.Q0
            irow = 0
            do i=1,NDUST
              if (.not.active(i)) cycle 
              if (i.ne.imaxon) then
                irow = irow+1
                dind(irow) = i
                do j=1,dust_nel(i)
                  el = dust_el(i,j)
                  erow = eind(el)
                  AA(erow,irow) = REAL(dust_nu(i,j),kind=qp)  
                enddo
              else
                dind(Nact) = i
                do j=1,dust_nel(i)
                  el = dust_el(i,j)
                  erow = eind(el)
                  bvec(erow) = REAL(dust_nu(i,j),kind=qp)  
                enddo
              endif  
            enddo  
            if (verbose>0) then
              print*,"searching for linear combination ..." 
              print'(2x,99(A10))',(trim(dust_nam(dind(i))),i=1,Nact) 
              do el=1,NELEM
                if (.not.e_act(el)) cycle
                erow = eind(el)
                print'(A2,99(F10.1))',trim(elnam(el)),AA(erow,1:Nact-1),
     >                                bvec(erow)
              enddo
            endif  
            call GAUSS_NM(NELEM,NDUST,Eact,Nact-1,AA,xvec,bvec,info)
            if (verbose>0) print'(" GAUSS_NM info =",I2)',info
            if (info.eq.0) then
              Nlin = 1
              dlin(1) = imaxon
              slin(imaxon) = -1.Q0
              do i=1,Nact-1
                if (ABS(xvec(i))<1.Q-25) cycle
                Nlin = Nlin+1
                dlin(Nlin) = dind(i)
                slin(dind(i)) = xvec(i)
              enddo  
              txt = trim(dust_nam(dlin(1)))//" <-> "
              do i=2,Nlin
                dk = dlin(i) 
                write(dum6,'(F6.3)') slin(dk)
                txt = trim(txt)//dum6//" "//trim(dust_nam(dk))
              enddo
              if (verbose>=0) then
                print*,"linear combination found: "//trim(txt)
              endif  
              itried(:) = .false.
              do
                Smin = 9.Q+99
                ioff = 0
                do i=1,Nlin
                  dk = dlin(i) 
                  if (dk==imaxon) cycle
                  if (Sat0(dk)<Smin.and.(.not.itried(dk))) then
                    Smin = Sat0(dk) 
                    ioff = dk
                  endif  
                enddo
                if (ioff==0) print*, "*** ioff=0 should not occur"
                amount = ddust(ioff)
                ok = .true.
                do i=1,Nlin
                  dk = dlin(i)
                  if (dk==ioff) cycle
                  if (ddust(dk)-slin(dk)/slin(ioff)*amount<0.Q0) then
                    ok=.false.
                  endif
                enddo
                if (ok) exit
                itried(ioff) = .true.
              enddo
              changed = .true.
              active(ioff) = .false.  
              Nt = REAL(Nlin-1,kind=qp)
              amount = ddust(ioff)/Nt
              do i=1,Nlin
                dk = dlin(i)
                if (dk==ioff) cycle
                call TRANSFORM(ioff,dk,amount,-slin(dk)/slin(ioff)*Nt,
     >                         ddust,eps,dscale,ok)
              enddo  
              ddust(ioff) = 0.Q0
              eps = eps_save
              if (verbose>=0) then
                print*,"switching off ",dust_nam(ioff)
              endif  
            endif  
          endif
        endif  
        if (ioff>0) itransform=itransform+1

        action = .false.
        do i=1,NDUST
          if (active(i).and.(.not.act_old(i))) then
            action = .true.
            laston = i
            switchedON(i) = switchedON(i)+1
            if (verbose>-1) then
              write(97,*) it,"on  ",dust_nam(i),switchedON(i)
            endif  
          else if (act_old(i).and.(.not.active(i))) then
            action = .true.
            switchedOFF(i) = switchedOFF(i)+1
            if (verbose>-1) then
              write(97,*) it,"off ",dust_nam(i),switchedOFF(i)
            endif  
          endif
        enddo
        no_action = no_action+1
        if (action) no_action=0
  
        if (changed) then
          Nact = 0 
          do i=1,NDUST
            if (active(i)) Nact=Nact+1
          enddo
          call SUPER(nHtot,T,null,eps,Sat0,NewFastLevel<1)
          qual = SQUAL(Sat0,active)
          lastit = it
        endif
        if (verbose>0) then
          do i=1,NDUST
            rem = "  "
            if (active(i)) rem=" *"
            if (verbose>=0.and.(active(i).or.Sat0(i)>0.1)) then
              write(*,'(3x,A18,2(1pE11.3)1pE19.10,I3,1pE11.3,A2)') 
     >          dust_nam(i),ddust(i),ddust(i)/dscale(i),Sat0(i),
     >          switchedOFF(i),pot(i),rem
            endif  
          enddo
          do i=1,NDUST
            if (active(i).and.(.not.act_old(i))) then
              print*,"... switching on "//trim(dust_nam(i)) 
            else if (.not.active(i).and.act_old(i)) then
              print*,"... switching off "//trim(dust_nam(i)) 
            endif
          enddo   
        endif
        act_old = active
        if (Nact==0.and.qual<1.Q-30) exit   ! no solid supersaturated 

        method = method_eqcond
        method_failed = 0
 50     continue
        
        !----------------------------  use method 1  ------------------------------------
        if (method==1) then
        !--------------------------------------
        ! ***  select independent elements  ***
        !--------------------------------------
        if (it>1.and.(.not.changed)) goto 100
        Nind = 1
        Iabund(1) = 9.D+99
        e_act(:) = .false.
        do i=1,NDUST
          if (.not.active(i)) cycle
          if (dust_nel(i)>1) cycle          ! include pure metals first
          el = dust_el(i,1)
          if (e_act(el)) cycle
          e_act(el) = .true.
          Iindex(2:Nind+1) = Iindex(1:Nind)
          Iabund(2:Nind+1) = Iabund(1:Nind)
          Iindex(1) = el
          Iabund(1) = 0.Q0
          !print*,1,elnam(Iindex(1:Nind))
          Nind = Nind+1
        enddo  
        do i=1,NDUST
          if (.not.active(i)) cycle         ! include elements in compounds
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            if (e_act(el)) cycle
            e_act(el) = .true.
            do ii=1,Nind                  
              if (eps(el)<Iabund(ii)) exit  ! sort by element abundance 
            enddo
            Iindex(ii+1:Nind+1) = Iindex(ii:Nind)
            Iabund(ii+1:Nind+1) = Iabund(ii:Nind)
            Iindex(ii) = el
            Iabund(ii) = eps(el) 
            !print*,2,elnam(Iindex(1:Nind))
            Nind = Nind+1
          enddo  
        enddo  
        e_act(:) = .false.
        e_num(:) = 0
        do i=1,Nact
          el = Iindex(i)
          e_act(el) = .true.
        enddo 
        do i=1,NDUST 
          if (.not.active(i)) cycle
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            e_num(el) = e_num(el)+1
          enddo  
        enddo  
        if (Nind-1<Nact) print*, "*** Nind<Nact in equil_cond."
        Nall = Nind-1
        Nind = Nact                         ! truncate at number of condensates
        if (verbose>1) print'(99(A3))',(elnam(Iindex(j)),j=1,Nall)
        if (verbose>1) print'(99(I3))',e_num(Iindex(1:Nall)) 

        !-----------------------------------------------
        ! ***  check and correct choice of elements  ***
        !-----------------------------------------------
        eblocked(:) = .false.
        do dk=1,NDUST
          if (.not.active(dk)) cycle
          e1 = 0
          e2 = 0
          do j=1,dust_nel(dk)
            el=dust_el(dk,j)
            if (e_num(el)==1) then
              e2 = e1
              e1 = el
            endif
          enddo
          if (e1>0.and.e2>0) then
            eblocked(e1) = .true.
            eblocked(e2) = .true.
            if (verbose>1) then
              print*,"*** blocked ",dust_nam(dk),elnam(e1),elnam(e2)
            endif  
          endif  
        enddo  
        e_num_save(:) = e_num(:)
        dtry = 0
 200    continue
        if (Nact<Nall) then
          e_num(:) = e_num_save
          e_eliminated(:) = .false.
          d_eliminated(:) = .false.
          do 
            if (verbose==2) then
              txt  = ''
              txt1 = ''
              txt2 = ''
              do dk=1,NDUST
                if (active(dk).and.(.not.d_eliminated(dk))) then
                  txt = trim(txt)//" "//trim(dust_nam(dk))
                endif  
              enddo
              do i=1,Nall
                el = Iindex(i)
                if (.not.e_eliminated(el)) then
                  txt1 = trim(txt1)//" "//elnam(el)
                  write(tnum,'(I2)') e_num(el) 
                  txt2 = trim(txt2)//" "//tnum
                endif  
              enddo
              print*,trim(txt)
              print*," "//trim(txt1)
              print*,trim(txt2)
            endif  
            found = .false. 
            do i=1,Nact
              el = Iindex(i)
              if (e_eliminated(el)) cycle
              if (e_num(el)==1) then
                do dk=1,NDUST
                  if (.not.active(dk)) cycle
                  if (d_eliminated(dk)) cycle
                  do j=1,dust_nel(dk)
                    if (el.eq.dust_el(dk,j)) then
                      if (verbose==2) then
                        print*,elnam(el)//" "//trim(dust_nam(dk)) 
                      endif  
                      found = .true.
                      exit
                    endif   
                  enddo  
                  if (found) exit
                enddo   
              endif   
              if (found) exit
            enddo  
            if (.not.found) exit
            e_eliminated(el) = .true.
            d_eliminated(dk) = .true.
            do j=1,dust_nel(dk)
              el = dust_el(dk,j)
              e_num(el) = e_num(el)-1
            enddo   
          enddo   
          !--- is there still a selected element with e_num=0?  ---
          found = .false.
          do i=1,Nact
            el = Iindex(i)
            if (e_eliminated(el)) cycle
            if (e_num(el)==0) then
              found = .true. 
              exit
            endif   
          enddo
          if (found) then
            found = .false. 
            do j=Nact+1,Nall 
              el = Iindex(j)
              if (e_num(el)>0) then
                found=.true.
                exit
              endif
            enddo  
            if (found) then
              !do i=1,Nact
              !  el = Iindex(i)
              !  print*,elnam(el),e_eliminated(el),e_num(el)
              !enddo  
              !do i=Nact,1,-1
              !  el = Iindex(i)
              !  if (e_eliminated(el)) cycle
              !  if (e_num(el)>0) exit
              !enddo
              !if (i<1) stop "*** should not occur"
              if (verbose>=0) then
                print*,"... exchanging1 "//elnam(Iindex(i))//
     >               " for "//elnam(Iindex(j))
                if (verbose>1) read(*,'(a1)') char1
              endif  
              swap = Iindex(i)   
              Iindex(i) = Iindex(j)
              Iindex(j) = swap
              e_act(Iindex(i)) = .true.
              e_act(Iindex(j)) = .false.
            endif
            if (.not.found) then
              print*,"*** no alternative element selection found."
              return 
            endif   
            goto 200 
          endif 
          !--- is there an unselected element with e_num=1?
          found = .false.
          !print*,"Nall,Nact=",Nall,Nact
          !if (Nall>Nact+1) then
          if (Nall>Nact) then
            do i=Nact+1,Nall
              el = Iindex(i)
              if (e_num(el)==1) then
                found = .true. 
                exit
              endif   
            enddo          
          endif
          if (found) then
            found = .false. 
            do j=Nact,1,-1
              el = Iindex(j)
              if (e_num(el)>0) then 
                found = .true.
                exit
              endif 
            enddo
            deplete1 = eps(Iindex(i))**2/eps0(Iindex(i))
            deplete2 = eps(Iindex(j))**2/eps0(Iindex(j))
            if (verbose>1) then
              print*,"deplete",elnam(Iindex(i)),elnam(Iindex(j))
              print*,deplete1,deplete2
            endif  
            if (dtry==0.and.deplete1>1.Q+2*deplete2) found=.false.
            if (found) then
              if (verbose>=0) then
                print*,"... exchanging2 "//elnam(Iindex(j))//
     >               " for "//elnam(Iindex(i))
                if (verbose>1) read(*,'(a1)') char1
              endif  
              swap = Iindex(i)   
              Iindex(i) = Iindex(j)
              Iindex(j) = swap
              e_act(Iindex(i)) = .true.
              e_act(Iindex(j)) = .false.
              goto 200 
            endif  
          endif
          !--- special cases ---
          found = (Nall==Nact+1.and.Iindex(Nall)==H.and.Iindex(Nact)==O)
     >       .and.active(iFe3O4).and.active(iFe2SiO4)
     >       .and.active(iFeAl2SiO7H2).and.active(iMg3Si2O9H4)
     >       .and.active(iCaMgSi2O6).and.active(iCa3Al2Si3O12)     
          if (found) then
            ! there is a linear-combination disregarding hydrogen
            i = Nall
            j = Nact
            if (verbose>=0) then
              print*,"... exchanging3 "//elnam(Iindex(j))//
     >             " for "//elnam(Iindex(i))
              if (verbose>1) read(*,'(a1)') char1
            endif  
            swap = Iindex(i)   
            Iindex(i) = Iindex(j)
            Iindex(j) = swap
            e_act(Iindex(i)) = .true.
            e_act(Iindex(j)) = .false.
            goto 200 
          endif
          !print*,Nall,Nact,elnam(Iindex(Nall)),e_num(Iindex(Nall))
          found = (Nall>=Nact+1.and.Iindex(Nall)==O
     >         .and.e_num(Iindex(Nall))==1)
     >         .and.active(iCaS).and.active(iCaSO4)
          if (found) then
            found = .false.
            do j=Nact,1,-1
              !print*,j,elnam(Iindex(j)),S,e_num(Iindex(j))
              if (Iindex(j)==S.and.e_num(Iindex(j))==2) then
                found = .true.
                exit
              endif 
            enddo
          endif  
          if (found) then
            ! can't have O as independent variable in this case
            i = Nall
            if (verbose>=0) then
              print*,"... exchanging4 "//elnam(Iindex(j))//
     >           " for "//elnam(Iindex(i))
              if (verbose>1) read(*,'(a1)') char1
            endif  
            swap = Iindex(i)   
            Iindex(i) = Iindex(j)
            Iindex(j) = swap
            e_act(Iindex(i)) = .true.
            e_act(Iindex(j)) = .false.
            goto 200 
          endif  
        endif   

 210    continue
        dtry_break = .false.
        if (verbose>1) then
          print*,"dependent elements ... ",
     >         (elnam(Iindex(i))//" ",i=Nind+1,Nall)
          print*,"solving for ... ",
     >         (elnam(Iindex(i))//" ",i=1,Nind)
          print'(99(1pE11.3))',(eps(Iindex(i)),i=1,Nind)
        endif

        !------------------------------------------------
        ! ***  determine dependent dust and elements  ***
        !------------------------------------------------
        d_resolved(1:NDUST) = .false.
        e_resolved(1:NELEM) = .false.
        e_act(:) = .false.
        e_taken(:) = .false.
        do i=1,Nind
          el = Iindex(i)
          e_act(el) = .true.
          e_resolved(el) = .true.
        enddo  
        Ndep = 0
        do i=1,NDUST
          if (.not.active(i)) cycle
          Ndep = Ndep+1
          Dindex(Ndep) = i
          is_dust(Ndep) =.true.
        enddo
        do i=1,NDUST
          if (.not.active(i)) cycle
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            if (e_resolved(el)) cycle
            if (e_taken(el)) cycle
            Ndep = Ndep+1
            Dindex(Ndep) = el
            is_dust(Ndep) = .false.
            e_taken(el) = .true.
          enddo  
        enddo

        !----------------------------
        ! ***  conversion matrix  ***
        !----------------------------
        Neq = 1
        do el=1,NELEM
          slots = 0
          do i=1,NDUST
            if (.not.active(i)) cycle
            do j=1,dust_nel(i)
              if (el.ne.dust_el(i,j)) cycle
              elem(Neq) = el
              slots = slots+1
              dustkind(Neq,slots) = i
              stoich(Neq,slots) = dust_nu(i,j)
            enddo
          enddo
          if (slots>0) then
            Nslot(Neq) = slots
            Neq=Neq+1
          endif  
        enddo  
        Neq = Neq-1
        do itry=1,99
          action = .false. 
          solved = .true.
          Nunsolved = 0
          do eq=1,Neq
            el = elem(eq)
            slots = Nslot(eq)
            knowns = 0
            unknown = 0
            if (e_resolved(el)) knowns=1
            do sl=1,slots
              dk = dustkind(eq,sl) 
              if (d_resolved(dk)) then
                knowns = knowns + 1
              else
                unknown = sl
              endif  
            enddo 
            unknowns = slots+1-knowns
            text = ""
            do sl=1,slots
              dk = dustkind(eq,sl)  
              write(txt1,'(I2)') stoich(eq,sl)
              txt0 = " "
              if (d_resolved(dk)) txt0="*"
              text = trim(text)//" "//trim(txt1)//" "
     >             //trim(dust_nam(dk))//txt0
            enddo
            write(txt1,'(I2,": ")') eq 
            write(txt2,'(" (",I2,")")') unknowns 
            txt0 = " "
            if (e_resolved(el)) txt0="*"
            if (verbose>1) print*,trim(txt1)//trim(text)//" + "
     >                    //trim(elnam(el))//txt0//trim(txt2)
            if (unknowns>1) solved=.false.
            if (unknowns==1) then
              action = .true. 
              call GET_KNOWNS(eq,mat,emat,e_act,d_resolved,e_resolved,
     >                        elem,Nslot,dustkind,stoich,vec,verbose)
              if (unknown>0) then 
                dk = dustkind(eq,unknown)
                d_resolved(dk) = .true.
                mat(dk,:) = -vec(:)/REAL(stoich(eq,unknown),kind=qp)
                do j=1,Nind
                  el2 = Iindex(j)
                  if (mat(dk,el2).eq.0.Q0) cycle
                  if (verbose>1) print*,"in1 ",trim(dust_nam(dk))
     >                      //" "//elnam(el2),REAL(mat(dk,el2))
                enddo  
              else
                e_resolved(el) = .true.
                emat(el,:) = -vec(:)
                do j=1,Nind
                  el2 = Iindex(j)
                  if (emat(el,el2).eq.0.Q0) cycle
                  if (verbose>1) print*,"in2 ",trim(elnam(el))//" "
     >                      //elnam(el2),REAL(emat(el,el2))
                enddo  
              endif
              exit
            else if (unknowns>1) then 
              Nunsolved = Nunsolved+1 
              unsolved(Nunsolved) = eq
            endif  
          enddo
          if (.not.action.and.Nunsolved==0) exit
          if (.not.action.and.Nunsolved>1) then
            !--------------------------------------------
            ! ***  solve N equations with N unknowns?  *** 
            !--------------------------------------------
            Nvar1 = 0
            do i=1,NDUST
              if (.not.active(i)) cycle
              if (d_resolved(i)) cycle
              Nvar1 = Nvar1+1
              var(Nvar1) = i
            enddo  
            Nvar2 = 0            
            do el=1,NELEM
              if (.not.e_taken(el)) cycle
              if (e_resolved(el)) cycle
              Nvar2 = Nvar2+1
              var(Nvar1+Nvar2) = el
            enddo  
            !write(*,*) Nunsolved,Nvar1,Nvar2
            !write(*,*) unsolved(1:Nunsolved)
            !write(*,*) dust_nam(var(1:Nvar1)), 
     >      !           elnam(var(Nvar1+1:Nvar1+Nvar2))
            if (verbose>1) print'("solving",I2," equations with",I2,
     >                 " unknowns ",99(A18))',Nunsolved,Nvar1+Nvar2,
     >                 dust_nam(var(1:Nvar1)), 
     >                 elnam(var(Nvar1+1:Nvar1+Nvar2))
            if (Nunsolved/=Nvar1+Nvar2) then
              print*,"... is impossible"
              dtry = dtry+1
              dtry_break = .true.
              if (dtry<Nind) exit
              return
            endif  
            if (Nunsolved==Nvar1+Nvar2) then
              DF = 0.Q0 
              FF = 0.Q0
              do i=1,Nunsolved
                eq = unsolved(i)
                call GET_KNOWNS(eq,mat,emat,e_act,d_resolved,e_resolved,
     >                     elem,Nslot,dustkind,stoich,vecs(i,:),verbose)
                slots = Nslot(eq)
                do j=1,Nvar1
                  dk = var(j)
                  do sl=1,slots
                    if (dk==dustkind(eq,sl)) then  
                      DF(i,j) = stoich(eq,sl)
                    endif  
                  enddo
                enddo
                el = elem(eq)
                do j=Nvar1+1,Nvar1+Nvar2
                  el2 = var(j)
                  if (el2==el) then
                    DF(i,j) = 1.Q0
                  endif  
                enddo
              enddo
              !--- compute inverse matrix ---
              DFsav = DF
              call QGEFA ( DF, NELEM, Nunsolved, ipvt, info )
              call QGEDI ( DF, NELEM, Nunsolved, ipvt, det, work, 1 )
              if (info.ne.0) then
                print*,"*** singular matrix in QGEFA: info=",info
                dtry = dtry+1
                dtry_break = .true.
                if (dtry<Nind) exit
                do i=1,Nunsolved
                  print'(99(1pE12.3))',(DFsav(i,j),j=1,Nunsolved)
                enddo
                return
              endif   
              do i=1,Nvar1
                dk = var(i)
                vec(:) = 0.Q0
                do j=1,Nunsolved
                  vec = vec + DF(i,j)*vecs(j,:)
                enddo
                mat(dk,:) = -vec(:)
                d_resolved(dk) = .true.
                do j=1,Nind
                  el2 = Iindex(j)
                  if (mat(dk,el2).eq.0.Q0) cycle
                  if (verbose>1) print*,"in3 ",trim(dust_nam(dk))//" "
     >                                //elnam(el2),REAL(mat(dk,el2))
                enddo  
              enddo
              do i=Nvar1+1,Nvar1+Nvar2
                el = var(i)
                vec(:) = 0.Q0
                do j=1,Nunsolved
                  vec = vec + DF(i,j)*vecs(j,:)
                enddo
                e_resolved(el) = .true.
                emat(el,:) = -vec(:)
                do j=1,Nind
                  el2 = Iindex(j)
                  if (emat(el,el2).eq.0.Q0) cycle
                  if (verbose>1) print*,"in4 ",trim(elnam(el))//" "
     >                            //elnam(el2),REAL(emat(el,el2))
                enddo  
              enddo  
            endif
            if (dtry>0.and.dtry_break) exit
          endif    
          if (itry==100) print*, "*** itry==100"
        enddo
        if (.not.solved) then
          if (dtry==1) goto 200   ! may work by relaxing the depletion-criterium
          dtry_last = dtry_last+1
          !print*,Nind,Ndep,Nall
          !do i=1,Nall
          !  el2 = Iindex(i)
          !  print*,elnam(el2),e_num(el2)
          !enddo  
          do i=Nind+1-dtry_last,1,-1
            el2 = Iindex(i)
            !print*,"from ",elnam(el2),e_num(el2)
            if (e_num(el2)>0) exit
          enddo  
          do j=Nind+1,Nall
            el2 = Iindex(j)
            !print*,"to   ",elnam(el2),e_num(el2),eblocked(el2)
            if (e_num(el2)>0.and.(.not.eblocked(el2))) exit
          enddo
          if (i>=1.and.j<=Nall) then
            if (verbose>=0) then
              print*,"... dtry_last=",dtry_last
              print*,"... exchanging5 "//elnam(Iindex(i))//
     >               " for "//elnam(Iindex(j))
              if (verbose>1) read(*,'(a1)') char1
            endif  
            swap = Iindex(i)   
            Iindex(i) = Iindex(j)
            Iindex(j) = swap
            e_act(Iindex(i)) = .true.
            e_act(Iindex(j)) = .false.
            goto 210
          endif
          write(*,*) "*** couldn't resolve the conversion matrix."
          return
        endif   
        do i=1,Nind
          el = Iindex(i) 
          do j=1,Ndep 
            if (is_dust(j)) then         
              dk  = Dindex(j)
              conv(j,i) = mat(dk,el)
            else
              el2 = Dindex(j) 
              conv(j,i) = emat(el2,el)
            endif
          enddo
        enddo  
        if (verbose>0) then
          print'(A24,99(A7))',"solving for ...",elnam(Iindex(1:Nind))
          do j=1,Ndep 
            if (is_dust(j)) then
              dk  = Dindex(j)
              txt = dust_nam(dk)
            else
              el  = Dindex(j) 
              txt = elnam(el)
            endif   
            print'(" conv.mat ",A14,99(0pF7.3))',trim(txt),
     >             (conv(j,i),i=1,Nind)
          enddo  
          !print'("max,min(conv.mat)=",2(1pE20.22))',
     >    !       MAXVAL(conv(1:Nind,1:Nind)),
     >    !       MINVAL(conv(1:Nind,1:Nind))
          !write(89,*) "---"
          !do i=1,Nind
          !  write(89,'(99(1pE30.22))'),(conv(i,j),j=1,Nind)
          !enddo
          !write(89,*)
        endif  

 100    continue
        !-----------------------------------------------
        ! ***  stop iteration of parts of solution?  ***
        !-----------------------------------------------
        if (changed.or.it==1) then
          converge(:,:) = 9.Q+99
          is_esolved(:) = .false.
          is_dsolved(:) = .false.
        else if (it>3) then
          Ntrivial = 0
          dtrivial = 0
          etrivial = 0
          do ii=1,Nsolve
            Nzero = 0
            nonzero = 0
            do jj=1,Nsolve
              if (ABS(DFsav(jj,ii))<1.Q-3) then
                Nzero = Nzero+1
              else
                nonzero = jj
              endif
            enddo
            if (Nzero==Nsolve-1) then
              Ntrivial = Ntrivial+1
              dtrivial(Ntrivial) = nonzero
              etrivial(Ntrivial) = ii
            endif  
          enddo  
          ebest = 0
          dbest = 0
          cbest = 9.Q+99
          do itrivial=1,Ntrivial
            ii = etrivial(itrivial) 
            jj = dtrivial(itrivial)
            i  = act_to_elem(ii)
            j  = act_to_dust(jj) 
            el = Iindex(i)
            dk = Dindex(j)
            crit = MAXVAL(ABS(converge(it-3:it-1,i)))  
            !print'(A2,A13,3(1pE12.3))',elnam(el),trim(dust_nam(dk)),
     >      !     converge(it-1,i),crit
            if (crit<1.Q-15.and.crit<cbest) then
              cbest = crit
              ebest = i
              dbest = j
            endif
          enddo
          if (ebest.ne.0) then
            el = Iindex(ebest) 
            dk = Dindex(dbest)
            is_esolved(ebest) = .true.
            is_dsolved(dbest) = .true.
            print*,elnam(el)//" (->"//trim(dust_nam(dk))//
     >             ") has converged."
          endif  
        endif  
        !-----------------------------------------------
        ! ***  fill in r.h.s. vector FF and          ***
        ! ***  determine current quality of solution ***
        !-----------------------------------------------
        !!! not necessary: xstep(:)= 0.Q0             
        !!! call SUPER(nHtot,T,xstep,eps,Sat0,NewFastLevel<1)
        qual = SQUAL(Sat0,active)
        ii = 0
        do i=1,Nind
          if (is_dsolved(i)) cycle
          ii = ii+1
          act_to_dust(ii) = i
          dk = Dindex(i)
          !FF(ii) = Sat0(dk)-1.Q0/Sat0(dk)
          FF(ii) = LOG(Sat0(dk))
        enddo 
        Nsolve = ii

        !------------------------------------------
        ! ***  compute numerical derivative DF  ***
        !------------------------------------------
        if (T>Tfast) then
          varfac = 1.Q-7
          target = 1.Q-6
          dlim   = 1.Q-12
        else
          varfac = 1.Q-20
          target = 1.Q-19
          dlim   = 1.Q-30
        endif  
        jj = 0
        do j=1,Nind
          if (is_esolved(j)) cycle
          jj = jj+1
          act_to_elem(jj) = j
          el = Iindex(j) 
          deps = -varfac*eps(el)                 ! limited by el abundance
          scale(j) = eps(el)                     ! consider decrease of el
          do i=1,Ndep                            
            if (is_dust(i)) cycle
            if (conv(i,j)==0.Q0) cycle
            el2 = Dindex(i)                       
            del = -varfac*eps(el2)/conv(i,j)        
            !print'(A3," <-",A3,99(1pE10.2))',elnam(el),elnam(el2),
     >      !     eps(el),deps,del,conv(i,j)
            if (ABS(del)<ABS(deps)) then         ! limited by dep. element?
              deps = del                         ! consider decrease of el2
              scale(j) = eps(el2)
            endif  
          enddo
          !print*,elnam(el),REAL(eps(el)),REAL(deps)
          LastRow = 0.Q0
          LastDev = 9.Q+99
          do itJac=1,99
            xstep(:) = 0.Q0
            xstep(j) = deps
            call SUPER(nHtot,T,xstep,eps,Sat2,NewFastLevel<1)
            dtmax = 0.Q0
            do ii=1,Nsolve
              i  = act_to_dust(ii) 
              dk = Dindex(i)
              !dterm = FF(ii)-Sat2(dk)+1.Q0/Sat2(dk)
              dterm = LOG(Sat0(dk)/Sat2(dk))
              dtmax = MAX(dtmax,ABS(dterm))
              DF(ii,jj) = dterm/deps*scale(j)
            enddo  
            !print'("JAC:",A3,99(1pE11.3))',elnam(el),
     >      !    eps(el),deps/eps(el),dtmax,DF(1:Nsolve,jj)
            if (dtmax<target) exit
            dev = 0.Q0
            do ii=1,Nsolve
              val1 = LastRow(ii)+SIGN(1.Q-8,LastRow(ii))
              val2 = DF(ii,jj)  +SIGN(1.Q-8,DF(ii,jj))
              dev = dev+(val1/val2-1.Q0)**2
            enddo
            !print*,"dev=",dev
            if (itJac>1.and.dev<1.Q-10.and.dev<LastDev) exit
            if (itJac>1.and.LastDev<1.E-2.and.dev>LastDev) then
              DF(1:Nsolve,jj) = LastRow(1:Nsolve) 
              exit
            endif  
            if (ABS(deps)<dlim*scale(j)) exit
            !--- decrease deps to get more precice DF-entries ---
            deps = deps * 2.Q-1
            LastRow(1:Nsolve) = DF(1:Nsolve,jj)
            LastDev = dev
          enddo  
        enddo
        !---- check functionality of DF ----
        !do jj=1,Nsolve
        !  j  = act_to_elem(jj)
        !  el = Iindex(j)
        !  deps = 1.Q+5*varfac*scale(j)
        !  xstep(:) = 0.Q0
        !  xstep(j) = deps
        !  call SUPER(nHtot,T,xstep,eps,Sat2,NewFastLevel<1)
        !  qualDF = 0.Q0
        !  do ii=1,Nsolve
        !    i  = act_to_dust(ii) 
        !    dk = Dindex(i)
        !    val1 = deps/scale(j)*DF(ii,jj)
        !    val2 = FF(ii)-LOG(Sat2(dk))
        !    qualDF = MAX(qualDF,ABS(val1-val2))
        !    print'(A14,A3,99(1pE18.11))',dust_nam(dk),elnam(el),
     >  !           LOG(Sat0(dk)),val1,val2,val1-val2
        !  enddo
        !enddo  
        !print*,"quality of DF():",REAL(qualDF)

        !---- print DF and FF ----
        !if (verbose>1) then
        !  print'(12x,99(A11))',elnam(Iindex(act_to_elem(1:Nsolve)))
        !  do ii=1,Nsolve
        !    i  = act_to_dust(ii) 
        !    dk = Dindex(i)
        !    print'(A18,99(1pE11.3))',dust_nam(dk),DF(ii,1:Nsolve),FF(ii)
        !  enddo  
        !endif

        !---------------------------------
        ! ***  Newton-Raphson step dx  ***
        !---------------------------------
        Fsav  = FF
        DFsav = DF
        !call GAUSS16( NELEM, Nsolve, DF, dx, FF)
        call QGEFA( DF, NELEM, Nsolve, ipvt, info )
        call QGESL( DF, NELEM, Nsolve, ipvt, FF, 0 )
        dx = FF
        if (verbose>1) print*,"QGESL info=",info
        if (info.ne.0) then
          print*,"*** singular matrix in QGEFA NR-step: info=",info
          if (method_failed==0) then
            method_failed = method_failed+1
            method = 2
            Nact = 0 
            do i=1,NDUST
              if (active(i)) Nact=Nact+1
            enddo
            changed = .true.
            print*,"trying eqcond_method 2 ..."
            goto 50
          endif  
          return
        endif
  
        !--- re-scale ---
        if (it>1) converge(it,:) = converge(it-1,:)
        do ii=1,Nsolve
          i = act_to_elem(ii) 
          el = Iindex(i) 
          dx(ii) = dx(ii)*scale(i)
          converge(it,i) = dx(ii)/eps(el)
        enddo  

        !-----------------------------------
        ! ***  limit NR step physically  ***
        !-----------------------------------
        fac = 1.Q0
        iminoff = 0
        limdust = .false.
        do ii=1,Nsolve
          i  = act_to_elem(ii) 
          el = Iindex(i)
          if (eps(el)+dx(ii)<0.02*eps(el)) then
            fac2 = (-0.98*eps(el))/dx(ii)               ! eps+fac*dx = 0.02*eps
            if (verbose>0) print'(" *** limiting element1 ",A2,
     >        " eps=",1pE9.2,"  fac=",1pE9.2)',elnam(el),eps(el),fac2
            if (fac2<fac) then
              fac = fac2 
            endif
          endif  
        enddo
        do j=Ndep,1,-1
          del = 0.Q0 
          do ii=1,Nsolve
            i = act_to_elem(ii)
            del = del + conv(j,i)*dx(ii)
          enddo 
          if (is_dust(j)) then
            dk = Dindex(j)
            if (laston>0.and.ddust(dk)+del<0.Q0) 
     >        print'("dk,laston,it,lastit=",2(A20),2(i4))',
     >             trim(dust_nam(dk)),trim(dust_nam(laston)),it,lastit
            if (del<0.Q0.and.ddust(dk)>0.1*dscale(dk)) then
              fac2 = (-ddust(dk)+0.05*dscale(dk))/del   ! ddust+fac*del = 0.05*dscale
              if (fac2<1.0.and.verbose>0) print*,"*** limiting dust 1 "
     >                              //dust_nam(dk),REAL(fac2)
              if (fac2<fac) then
                fac = fac2
                iminoff = 0
                limdust = .true.
              endif  
            else if (ddust(dk)+del<0.Q0.and.
     >               dk==laston.and.it<lastit+5) then
              fac2 = -0.9*ddust(dk)/del                 ! ddust+fac*del = ddust/10
              if (fac2<1.0.and.verbose>0) print*,"*** limiting dust 2 "
     >                              //dust_nam(dk),REAL(fac2)
              if (fac2<fac) then
                fac = fac2
                iminoff = 0
                limdust = .true.
              endif  
            else if (ddust(dk)+del<0.Q0) then
              fac2 = (-ddust(dk)-small*dscale(dk))/del  ! ddust+fac*del = -small*dscale
              if (fac2<1.0.and.verbose>0) print*,"*** limiting dust 3 "
     >                              //dust_nam(dk),REAL(fac2)
              if (fac2<fac) then
                fac = fac2 
                iminoff = dk
                limdust = .true.
              endif
            endif  
          else  
            el = Dindex(j)
            if (eps(el)+del<0.3*eps(el)) then
              fac2 = (-0.7*eps(el))/del                 ! eps+fac*del = 0.3*eps
              if (verbose>0) print'(" *** limiting element2 ",A2,
     >          " eps=",1pE9.2,"  fac=",1pE9.2)',elnam(el),eps(el),fac2
              if (fac2<fac) then
                fac = fac2 
                iminoff = 0
                limdust = .false.
              endif  
            endif
            !if (eps(el)+del>3.0*eps(el)) then
            !  fac2 = (2.0*eps(el))/del                ! eps+fac*dx = 3.0*eps
            !  if (verbose>0) print'(" *** limiting element2a ",A2,
     >      !    " eps=",1pE9.2,"  fac=",1pE9.2)',elnam(el),eps(el),fac2
            !  if (fac2<fac) then
            !    fac = fac2 
            !    iminoff = 0
            !    limdust = .false.
            !  endif  
            !endif
          endif  
        enddo  
        dx = dx*fac
        limited = (fac<1.Q0)
        !if (iminoff>0.and.(iminoff.ne.laston)) then
        if (iminoff>0) then
          if (verbose>=0) print*,"will switch off ",dust_nam(iminoff) 
          active(iminoff) = .false.
          lastit = -99
          !if (iminoff.eq.laston) then
          !  print*,"=> fall back"
          !  active = save_active
          !  eps = save_eps
          !  ddust = save_ddust
          !endif  
        endif

        !------------------------------------
        ! ***  apply dx to ddust and eps  ***
        !------------------------------------
        eps_save = eps
        dust_save = ddust
        NRstep = 1.Q0
        qold = qual
        do iback=1,7
          eps = eps_save
          ddust = dust_save
          do ii=1,Nsolve
            i = act_to_elem(ii) 
            el = Iindex(i)
            eps(el) = eps(el) + NRstep*dx(ii)      ! direct effect
          enddo  
          do j=1,Ndep
            del = 0.Q0 
            do ii=1,Nsolve
              i = act_to_elem(ii)
              del = del + conv(j,i)*NRstep*dx(ii)  ! effect of indep.element i
            enddo 
            if (is_dust(j)) then
              dk = Dindex(j)
              ddust(dk) = ddust(dk) + del
            else  
              el = Dindex(j)
              eps(el) = eps(el) + del
            endif  
          enddo
          call SUPER(nHtot,T,null,eps,Sat0,NewFastLevel<1)
          qual = SQUAL(Sat0,active)
          if (verbose>0) print'("--> pullback",i3,0pF6.3," Q=",1pE11.3,
     >                          " ->",1pE11.3)',iback,NRstep,qold,qual
          if (qual<qold*1.5) exit
          if (qual<1.0) exit
          if (limdust) exit
          if (iback==1) then
            NRstep = 0.9*NRstep
          else  
            NRstep = 0.5*NRstep
          endif  
        enddo
        endif ! method 1


        !--------------------------  use method 2  ------------------------------------
        if (method==2) then
        !------------------------
        ! ***  sort elements  ***
        !------------------------
        if (it>1.and.(.not.changed).and.(it<last_tri+10)) goto 300
        last_tri = it
        Nact = 0
        Iabund(1) = 9.D+99
        e_act(:) = .false.
        ii = 0
        do i=1,NDUST
          if (.not.active(i)) cycle         
          ii = ii+1
          act_to_dust(ii) = i               ! index of active condensate 
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            if (e_act(el)) cycle
            e_act(el) = .true.
            Nact = Nact + 1
            do k=1,Nact               
              if (eps(el)<Iabund(k)) exit   ! sort by element abundance 
            enddo
            Iindex(k+1:Nact+1) = Iindex(k:Nact)
            Iabund(k+1:Nact+1) = Iabund(k:Nact)
            Iindex(k) = el
            Iabund(k) = eps(el) 
          enddo  
        enddo  
        Nsolve = ii
        do i=1,Nact
          el = Iindex(i)
          e_num(el) = i
        enddo  

        !---------------------------------------------------------------
        ! ***  find useful linear combination of solids to vary ...  ***
        !---------------------------------------------------------------
        if (verbose>0) then
          print'(99(A8))',elnam(Iindex(1:Nact))
          print'(99(1pE8.1))',Iabund(1:Nact)
          print'(20x,99(A4))',elnam(Iindex(1:Nact))
        endif  
        DF(:,:) = 0.Q0
        do ii=1,Nsolve
          i = act_to_dust(ii)
          base(ii,:) = 0.Q0
          base(ii,i) = 1.Q0
          do j=1,dust_nel(i)
            el = dust_el(i,j)
            jj = e_num(el)
            DF(ii,jj) = dust_nu(i,j)
          enddo
          if (verbose>1) print'(A20,99(0pF4.0),A20)',
     >                   dust_nam(i),DF(ii,1:Nact)
        enddo
        !--------------------------------------------------------------
        ! ***  ... by triangulation of the stoichiometry matrix DF  ***
        !--------------------------------------------------------------
        do ii=1,Nsolve-1
          imax = ii
          amax = ABS(DF(ii,ii))
          do i=ii+1,Nsolve
            if (ABS(DF(i,ii))>amax) then
              imax = i
              amax = ABS(DF(i,ii))
            endif
          enddo  
          if (imax>ii) then
            do jj=1,Nact
              tmp1        = DF(ii,jj)
              DF(ii,jj)   = DF(imax,jj)
              DF(imax,jj) = tmp1 
            enddo
            tmp2(:)      = base(ii,:)
            base(ii,:)   = base(imax,:)
            base(imax,:) = tmp2(:) 
          endif
          if (DF(ii,ii)==0.Q0) then
            print*,"*** WARNING: triangulation impossible in equil_cond"
            cycle  
          endif  
          do k=ii+1,Nsolve
            fac = DF(k,ii)/DF(ii,ii)
            do jj=1,Nact
              DF(k,jj) = DF(k,jj) - fac*DF(ii,jj)
            enddo
            DF(k,ii) = 0.Q0                            ! exact by construction
            base(k,:) = base(k,:) - fac*base(ii,:)
          enddo
        enddo  
        base_el = 0.Q0                                 ! effect of 1 base-vec on elements
        do ii=1,Nsolve
          do k=1,Nsolve
            i = act_to_dust(k)
            if (base(ii,i)==0.Q0) cycle
            do j=1,dust_nel(i)
              el = dust_el(i,j)
              base_el(ii,el) = base_el(ii,el) + base(ii,i)*dust_nu(i,j)
            enddo
          enddo
          do i=1,NELM
            if (i==iel) cycle
            el = elnum(i)
            if (ABS(base_el(ii,el))<1.Q-10) base_el(ii,el)=0.Q0   ! avoid round-off errors
          enddo  
        enddo  
        if (verbose>0) then
          print'(99(A6))',elnam(Iindex(1:Nact))
          do ii=1,Nsolve
            print'(99(0pF6.2))',DF(ii,1:Nact)
          enddo  
          do ii=1,Nsolve
            i = act_to_dust(ii)
            text = ''
            do jj=1,Nsolve
              j = act_to_dust(jj)
              if (base(ii,j)==0.Q0) cycle
              write(txt1,'(0pF7.4)') base(ii,j)
              text = trim(text)//" "//trim(txt1)//" "
     >               //trim(dust_nam(j))
            enddo
            print'(" base-vec ",I2,": ",A)',ii,trim(text)
            text = ''
            do i=1,NELM
              if (i==iel) cycle
              el = elnum(i)
              if (base_el(ii,el)==0.Q0) cycle
              text = trim(text)//" "//trim(elnam(el))
            enddo
            print'("          modifies: ",A)',trim(text)
          enddo
        endif
  
 300    continue
        !-----------------------------------------------
        ! ***  fill in r.h.s. vector FF and          ***
        ! ***  determine current quality of solution ***
        !-----------------------------------------------
        Nind = 0  
        Ndep = 0  
        !not necessary call SUPER(nHtot,T,null,eps,Sat0,NewFastLevel<1)
        qual = SQUAL(Sat0,active)
        ii = 0
        do ii=1,Nsolve
          i = act_to_dust(ii)                ! index of active condensate 
          FF(ii) = LOG(Sat0(i))              ! the function to be nullified
        enddo  

        !------------------------------------------
        ! ***  compute numerical derivative DF  ***
        !------------------------------------------
        if (T>Tfast) then
          varfac = 1.Q-7
          target = 1.Q-5
          dlim   = 1.Q-12
        else
          varfac = 1.Q-20
          target = 1.Q-18
          dlim   = 1.Q-30
        endif  
        DF(:,:) = 0.Q0
        LastRow = 0.Q0
        LastDev = 9.Q+99
        do jj=1,Nsolve
          scale(jj) = 9.Q+99
          do i=1,NELM
            if (i==iel) cycle
            el = elnum(i)
            if (base_el(jj,el)==0.Q0) cycle
            scale(jj) = MIN(scale(jj),eps(el))
          enddo  
          del = varfac*scale(jj)             ! tiny amount of lin.comb.of.cond. to evaporate
          do itJac=1,99
            check(:) = eps(:) + del*base_el(jj,:)
            call SUPER(nHtot,T,null,check,Sat2,NewFastLevel<1)
            dtmax = 0.Q0            
            do ii=1,Nsolve
              i = act_to_dust(ii) 
              dterm = LOG(Sat0(i)/Sat2(i))
              dtmax = MAX(dtmax,ABS(dterm))
              DF(ii,jj) = dterm/del*scale(jj)
            enddo  
            !print'("JAC: base-vec ",I2,99(1pE10.2))',
     >      !    jj,del/scale(jj),dtmax,DF(1:Nsolve,jj)
            if (dtmax<target) exit
            dev = 0.Q0
            do ii=1,Nsolve
              val1 = LastRow(ii)+SIGN(1.Q-8,LastRow(ii))
              val2 = DF(ii,jj)  +SIGN(1.Q-8,DF(ii,jj))
              dev = dev+(val1/val2-1.Q0)**2
            enddo
            !print*,"dev=",dev
            if (itJac>1.and.dev<1.Q-10.and.dev<LastDev) exit
            if (itJac>1.and.LastDev<1.E-2.and.dev>LastDev) then
              DF(1:Nsolve,jj) = LastRow(1:Nsolve) 
              exit
            endif  
            if (ABS(del)<dlim*scale(jj)) exit
            !--- decrease del to get more precice DF-entries ---
            del = del * 2.Q-1
            LastRow(1:Nsolve) = DF(1:Nsolve,jj)
            LastDev = dev
          enddo
          !j = act_to_dust(jj)
          !print'("JAC:",A18,99(1pE10.2))',
     >    !     trim(dust_nam(j)),DF(1:Nsolve,jj)
        enddo
  
        !---- check functionality of DF ----
        !do jj=1,Nsolve
        !  j = act_to_dust(jj)
        !  del = 1.Q-5*scale(jj)
        !  check = eps             
        !  do i=1,dust_nel(j)
        !    el = dust_el(j,i)
        !    check(el) = check(el)+del*dust_nu(j,i)
        !  enddo
        !  call SUPER(nHtot,T,null,check,Sat2,NewFastLevel<1)
        !  qualDF = 0.Q0
        !  do ii=1,Nsolve
        !    i  = act_to_dust(ii) 
        !    val1 = del/scale(jj)*DF(ii,jj)
        !    val2 = FF(ii)-LOG(Sat2(i))
        !    qualDF = MAX(qualDF,ABS(val1-val2))
        !    !print'(A14,A14,99(1pE18.11))',dust_nam(i),dust_nam(j),
     >  !    !       LOG(Sat0(i)),val1,val2,val1-val2
        !  enddo
        !enddo  
        !if (verbose>1) print*,"quality of DF():",REAL(qualDF)

        !---------------------------------
        ! ***  Newton-Raphson step dx  ***
        !---------------------------------
        Fsav  = FF
        DFsav = DF
        call QGEFA( DF, NELEM, Nsolve, ipvt, info )
        call QGESL( DF, NELEM, Nsolve, ipvt, FF, 0 )
        dx = FF
        if (verbose>1) print*,"QGESL info=",info
        if (info.ne.0) then
          print*,"*** singular matrix in QGEFA NR-step: info=",info
          if (method_failed==0) then
            method_failed = method_failed+1
            method = 1
            Nact = 0 
            do i=1,NDUST
              if (active(i)) Nact=Nact+1
            enddo
            changed = .true.
            print*,"trying eqcond_method 1 ..."
            goto 50
          endif  
          return
        endif
        dstep = 0.Q0                 ! the NR-step in condensate abundances 
        do ii=1,Nsolve
          dx(ii) = dx(ii)*scale(ii)  ! de-normalisation
          do k=1,Nsolve
            i = act_to_dust(k)
            if (base(ii,i)==0.Q0) cycle
            dstep(i) = dstep(i) + dx(ii)*base(ii,i)
          enddo
          !test = 0.Q0
          !do jj=1,Nsolve
          !  test = test + DFsav(ii,jj)*(dx(jj)/scale(jj))
          !enddo
          !print'(I3,2(1pE18.10),1pE9.2)',ii,Fsav(ii),test, 
     >    !                                  Fsav(ii)/test-1.Q0
          if (IS_NAN(REAL(dx(ii)))) then
            print*,NELEM,Nsolve
            do i=1,Nsolve
              print'(99(1pE13.6))',DFsav(i,1:Nsolve),Fsav(i)
            enddo  
            print'(99(1pE13.6))',dx(1:Nsolve)
        	print*, "*** dx=NaN"
            return
          endif  
        enddo  
        xstep = 0.Q0                 ! the corresponding NR-step in element abundances
        do ii=1,Nsolve
          xstep(:) = xstep(:) + dx(ii)*base_el(ii,:)
        enddo  

        !-----------------------------------
        ! ***  limit NR step physically  ***
        !-----------------------------------
        fac = 1.Q0
        iminoff = 0
        limdust = .false.
        do ii=1,Nsolve
          i = act_to_dust(ii)
          del = -dstep(i)
          if (del<0.Q0.and.ddust(i)>0.1*dscale(i)) then ! try small ddust before switching off
            fac2 = (-ddust(i)+0.05*dscale(i))/del       ! ddust+fac*del = 0.05*dscale
            if (fac2<1.Q0.and.verbose>0) print*,"*** limiting dust 1 "
     >                                   //dust_nam(i),REAL(fac2)
            if (fac2<fac) then
              fac = fac2
              iminoff = 0
              !limdust = .true.
            endif  
          else if (ddust(i)+del<0.Q0.and.
     >             i==laston.and.it<lastit+5) then      ! just switched on: keep on trying
            fac2 = -0.9*ddust(i)/del                    ! ddust+fac*del = ddust/10
            if (fac2<1.Q0.and.verbose>0) print*,"*** limiting dust 2 "
     >                                  //dust_nam(i),REAL(fac2)
            if (fac2<fac) then
              fac = fac2
              iminoff = 0
              limdust = .true.
            endif  
          else if (ddust(i)+del<0.Q0) then              ! preparation to switch off next
            fac2 = (-ddust(i)-small*dscale(i))/del      ! ddust+fac*del = -small*dscale
            if (fac2<1.Q0.and.verbose>0) print*,"*** limiting dust 3 "
     >                                   //dust_nam(i),REAL(fac2)
            if (fac2<fac) then
              fac = fac2 
              iminoff = i
              limdust = .true.
            endif
          endif  
        enddo
        do i=1,NELM
          if (i==iel) cycle
          el = elnum(i)
          del = xstep(el)
          if (eps(el)+del<0.02*eps(el)) then
            fac2 = (-0.98*eps(el))/del                  ! eps+fac*del = 0.02*eps
            if (verbose>0) print'(" *** limiting element ",A2,
     >        " eps=",1pE9.2,"  fac=",1pE9.2)',elnam(el),eps(el),fac2
            if (fac2<fac) then
              fac = fac2 
              iminoff = 0
              limdust = .false.
            endif
          endif
        enddo  
        dstep = dstep*fac
        xstep = xstep*fac
        limited = (fac<1.Q0)
        if (iminoff>0) then
          if (verbose>=0) print*,"will switch off ",dust_nam(iminoff) 
          active(iminoff) = .false.
          lastit = -99
        endif

        !if (verbose>0) then
        !  do ii=1,Nsolve
        !    i = act_to_dust(ii) 
        !    print'(A16,3(1pE18.10))',dust_nam(i),ddust(i),-dstep(i)
        !  enddo  
        !  do i=1,NELM
        !    if (i==iel) cycle
        !    el = elnum(i)
        !    print'(A3,2(1pE18.10))',elnam(el),eps(el),xstep(el)
        !  enddo
        !endif  

        !------------------------------------------------
        ! ***  apply dstep to ddust and xstep to eps  ***
        !------------------------------------------------
        eps_save = eps
        dust_save = ddust
        NRstep = 1.Q0
        qold = qual
        do iback=1,7
          eps = eps_save
          ddust = dust_save
          do ii=1,Nsolve
            i = act_to_dust(ii)
            ddust(i) = ddust(i) - NRstep*dstep(i)
          enddo
          do i=1,NELM
            if (i==iel) cycle
            el = elnum(i)
            eps(el) = eps(el) + NRstep*xstep(el)
          enddo  
          call SUPER(nHtot,T,null,eps,Sat0,NewFastLevel<1)
          qual = SQUAL(Sat0,active)
          if (verbose>0) print'("--> pullback",i3,0pF6.3," Q=",1pE11.3,
     >                          " ->",1pE11.3)',iback,NRstep,qold,qual
          if (qual<qold*1.5) exit
          if (qual<1.0) exit
          if (limdust) exit
          if (iback==1) then
            NRstep = 0.9*NRstep
          else  
            NRstep = 0.5*NRstep
          endif  
        enddo
        endif ! method 2

        !-------------------------------------
        ! ***  check element conservation  ***
        !-------------------------------------
        !check = eps
        !do i=1,NDUST
        !  do j=1,dust_nel(i)
        !    el = dust_el(i,j)
        !    check(el) = check(el) + ddust(i)*dust_nu(i,j)    
        !  enddo
        !enddo
        !worst = 0.Q0
        !do i=1,NEPS
        !  el = elnr(i)
        !  worst = MAX(worst,ABS(1.Q0-check(el)/eps00(el)))
        !enddo
        !print*,"element conservation error 2:",worst

        !--- not necessary to call again ---
        !xstep(:) = 0.Q0
        !call SUPER(nHtot,T,xstep,eps,Sat0,NewFastLevel<1)
        !qual = SQUAL(Sat0,active)
        Smax = maxval(Sat0)
        print'("it =",I4,"  qual =",1pE11.4,"  Smax-1 =",1pE11.2E4)',
     >          it,qual,Smax-1.Q0
        if ((Smax<1.Q0+Sfinish).and.(qual<Qfinish)) exit
        if (verbose>0) read(*,'(a1)') char1
        if (verbose>0.and.char1=='1') method=1;changed=.true.
        if (verbose>0.and.char1=='2') method=2

      enddo  
      Sat = Sat0

      call CPU_TIME(time1)
      if (it.lt.itmax) then
        if (verbose>=0) then
          write(*,'("EQUIL_COND converged after ",I3," iter, time =",
     >          0pF7.3," CPU sec.")') it,time1-time0
        endif  
      else
        write(*,'("*** EQUIL_COND failed after ",I3," iter,  time =",
     >            0pF9.4," CPU sec.")') it,time1-time0 
        return
      endif

      !-------------------------
      ! ***  check solution  ***
      !-------------------------
      do i=1,NDUST
        if (ddust(i)>0.Q0.and.Sat(i)<0.9999) then
          print*,"*** error: ddust>0 but S<1"
          print*,dust_nam(i),REAL(ddust(i)),REAL(Sat(i))
          return
        endif
        if (Sat(i)>1.00001) then
          print*,"*** error: S>1"
          print*,dust_nam(i),REAL(ddust(i)),REAL(Sat(i))
          return
        endif
        if (ddust(i)<-10*small*dscale(i)) then
          print*,"*** error: ddust<0"
          print*,dust_nam(i),REAL(ddust(i)),REAL(Sat(i))
          return
        endif  
      enddo  

      !-------------------------------------
      ! ***  check element conservation  ***
      !-------------------------------------
      check = eps
      do i=1,NDUST
        do j=1,dust_nel(i)
          el = dust_el(i,j)
          check(el) = check(el) + ddust(i)*dust_nu(i,j)
        enddo
      enddo
      do i=1,NEPS
        el = elnr(i)
        if (ABS(1.Q0-check(el)/eps00(el))>1.Q-8) then
          print*,"*** element conservation error 1"
          print*,elnam(el),check(el),eps0(el),eps00(el)
          return
        endif  
      enddo

      !----------------------------------
      ! ***  save result to database  ***
      !----------------------------------
      !if (qual<1.Q-10.and.useDatabase) then
      if (useDatabase) then
        call PUT_DATA(nHtot,T,eps,ddust,qread,iread,active)
      endif  
      ieqcond = ieqcond + 1
      ieqconditer = ieqconditer + it

      end
            

!-------------------------------------------------------------------------
      subroutine SUPER(nHtot,T,xx,eps,Sat,merk)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nel,dust_nu,dust_el,
     >                    dust_nam,elnam
      use EXCHANGE,ONLY: nel,nat,nion,nmol
      use CHEMISTRY,ONLY: NELM,elnum,iel=>el
      use CONVERSION,ONLY: Ndep,Nind,Dindex,Iindex,is_dust,conv
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real*8,intent(in)  :: nHtot,T
      real(kind=qp),intent(in) :: xx(NELEM),eps(NELEM)
      real(kind=qp),intent(out) :: Sat(NDUST)
      logical,intent(in) :: merk
      real(kind=qp) :: eps1(NELEM),dx
      integer :: i,j,el

      !-------------------------------------------
      ! ***  compute remaining gas abundances  ***
      !-------------------------------------------
      eps1 = eps
      do i=1,Nind
        el = Iindex(i)
        eps1(el) = eps1(el) + xx(i)  ! direct effect
      enddo  
      do j=1,Ndep
        if (is_dust(j)) cycle
        dx = 0.Q0 
        do i=1,Nind
          dx = dx + conv(j,i)*xx(i)  ! effect of indep.elements i on el j
        enddo 
        el = Dindex(j)
        eps1(el) = eps1(el) + dx
      enddo

      do i=1,NELM
        if (i==iel) cycle
        el = elnum(i)
        if (eps1(el).le.0.Q0) then
          write(*,*) "*** negative el.abund. SUPER",elnam(el),eps1(el)
          return
        endif  
      enddo
      
      !----------------------------------------------
      ! ***  compute chemistry & supersaturation  ***
      !----------------------------------------------
      call GGCHEM(nHtot,T,eps1,merk,0)
      call SUPERSAT(T,nat,nmol,Sat)
      
      end
      

!-------------------------------------------------------------------------
      function SQUAL(Sat,active)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NDUST
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp),intent(in) :: Sat(NDUST)
      real(kind=qp) :: SQUAL,qual
      logical,intent(in) :: active(0:NDUST)
      integer :: i

      qual = 0.d0
      do i=1,NDUST
        if (active(i)) then
         !qual = qual + (1.Q0-Sat(i))**2
         !qual = qual + (Sat(i)-1.Q0/Sat(i))**2
          qual = qual + LOG(Sat(i))**2
        else if (Sat(i).gt.1.Q0) then
         !qual = qual + MIN(Sat(i)-1.Q0,1.Q0)
        endif  
      enddo
      SQUAL = qual
      end

!-------------------------------------------------------------------------
      subroutine VAPORIZE(i,ddust,eps)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nel,dust_nu,dust_el,dust_nam
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: i
      real(kind=qp),intent(inout) :: ddust(NDUST),eps(NELEM)
      real(kind=qp) :: del
      integer :: j,el
      
      del = ddust(i)
      print*," ==>  vaporize "//trim(dust_nam(i)),REAL(del)
      ddust(i) = 0.Q0
      do j=1,dust_nel(i)
        el = dust_el(i,j)
        eps(el) = eps(el) + del*dust_nu(i,j)    
      enddo
      end

!-------------------------------------------------------------------------
      subroutine TRANSFORM(i1,i2,del,fac,ddust,eps,dscale,ok)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nel,dust_nu,dust_el,dust_nam,
     >                    eps0
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: i1,i2
      real(kind=qp),parameter :: dsmall=1.Q-30
      real(kind=qp),intent(inout) :: ddust(NDUST),eps(NELEM)
      real(kind=qp),intent(in) :: del,fac,dscale(NDUST)
      logical,intent(inout) :: ok
      integer :: j,el
      
      print*," ==>  transform "//trim(dust_nam(i1))//" -> "
     &       //trim(dust_nam(i2)),REAL(fac*del/dscale(i1))
      ddust(i1) = ddust(i1)-del
      ddust(i2) = ddust(i2)+fac*del
      do j=1,dust_nel(i1)
        el = dust_el(i1,j)
        eps(el) = eps(el) + del*dust_nu(i1,j)    
      enddo
      do j=1,dust_nel(i2)
        el = dust_el(i2,j)
        eps(el) = eps(el) - fac*del*dust_nu(i2,j)    
      enddo

      if (ddust(i1)<-dsmall.or.ddust(i2)<-dsmall) ok=.false.
      !if (ddust(i1)<-dsmall) then
      !  call VAPORIZE(i1,ddust,eps)
      !  active(i1) = .false.
      !  vap = .true.
      !endif  
      !if (ddust(i2)<-dsmall) then
      !  call VAPORIZE(i2,ddust,eps)
      !  active(i2) = .false.
      !  vap = .true.
      !endif
  
      !-------------------------------------
      ! ***  check element conservation  ***
      !-------------------------------------
      !check = eps
      !do i=1,NDUST
      !  do j=1,dust_nel(i)
      !    el = dust_el(i,j)
      !    check(el) = check(el) + ddust(i)*dust_nu(i,j)    
      !  enddo
      !enddo
      !worst = 0.d0
      !do i=1,NEPS
      !  el = elnr(i)
      !  worst = MAX(worst,ABS(1.Q0-check(el)/eps00(el)))
      !enddo
      !write(*,*) "element conservation error 3:",worst
      !if (worst>1.Q-8) stop

      end

!-------------------------------------------------------------------------
      subroutine GET_KNOWNS(eq,mat,emat,e_act,d_resolved,e_resolved,
     >                      elem,Nslot,dustkind,stoich,vec,verbose)
!-------------------------------------------------------------------------
      use DUST_DATA,ONLY: NELEM,NDUST,dust_nam,dust_nel,dust_nu,dust_el,
     >                    elnam 
      use CONVERSION,ONLY: Nind,Ndep,Iindex,Dindex,is_dust
      implicit none
      integer,parameter :: qp = selected_real_kind ( 33, 4931 )
      integer,intent(in) :: eq
      integer,intent(in),dimension(NELEM) :: elem,Nslot
      integer,intent(in),dimension(NELEM,NDUST) :: dustkind,stoich
      logical,intent(in),dimension(NDUST) :: d_resolved
      logical,intent(in),dimension(NELEM) :: e_resolved,e_act
      integer,intent(in) :: verbose
      real(kind=qp),intent(in) :: mat(NDUST,NELEM)
      real(kind=qp),intent(in) :: emat(NELEM,NELEM)
      real(kind=qp),intent(out):: vec(NELEM)
      integer :: i,j,el,el2,slots,sl

      vec(:) = 0.Q0
      slots = Nslot(eq)
      do sl=1,slots
        i = dustkind(eq,sl)
        if (.not.d_resolved(i)) cycle
        do j=1,Nind
          el2 = Iindex(j)
          if (mat(i,el2).eq.0.Q0) cycle
          if (verbose>1) print*,"out1 ",trim(dust_nam(i))//" "
     >           //elnam(el2),REAL(mat(i,el2)),stoich(eq,sl)
          vec(el2) = vec(el2)+ mat(i,el2)*stoich(eq,sl)
        enddo
      enddo  
      el = elem(eq)
      if (e_act(el)) then
        vec(el) = vec(el) + 1.Q0 
      else if (e_resolved(el)) then
        do j=1,Nind
          el2 = Iindex(j)
          if (emat(el,el2).eq.0.Q0) cycle
          vec(el) = vec(el) + emat(el,el2)
          if (verbose>1) print*,"out2 ",trim(elnam(el))//" "
     >           //elnam(el2),REAL(emat(el,el2))
        enddo  
      endif
      end
      function upper(strIn) result(strOut)
      implicit none
      character(len=*),intent(in) :: strIn
      character(len=len(strIn)) :: strOut
      integer :: i,j,l
      logical :: change

      change = .true.
      l = len(strIn)
      do i=1,l
        if (i<l-4) then 
          if (strIn(i:i+4)=='trans') change=.false. 
        endif  
        if (i<l-2) then 
          if (strIn(i:i+2)=='cis'  ) change=.false. 
        endif  
        j = iachar(strIn(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") .and. change) then
          strOut(i:i) = achar(iachar(strIn(i:i))-32)
        else
          strOut(i:i) = strIn(i:i)
        end if
      enddo

      end
