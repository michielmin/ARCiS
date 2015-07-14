	subroutine SetupOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE

	if(compute_opac) then
		call ComputeOpacities()
	else
		call ReadOpacities()
	endif

	return
	end


	subroutine ReadOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer imol
	real*8 kappa(ng),nu1,nu2,tanscale
	real*8 x1,x2,rr,gasdev,random,dnu,Saver,starttime,stoptime
	real*8,allocatable :: k_line(:),nu_line(:),dnu_line(:)
	real*8,allocatable :: opac_tot(:,:),cont_tot(:),kaver(:),kappa_mol(:,:,:)
	integer n_nu_line,iT
	integer i,j,ir,k,nl,ig
	integer,allocatable :: inu1(:),inu2(:)
	character*500 filename

	allocate(cont_tot(nlam))
	allocate(kaver(nlam))
	allocate(opac_tot(nlam,ng))
	allocate(kappa_mol(nmol,nlam,ng))

	n_nu_line=ng*nmol
	allocate(nu_line(n_nu_line))
	
	opac_tot=0d0

c	call cpu_time(starttime)
	do ir=nr,1,-1
		call tellertje(nr-ir+1,nr)
c		call output("Opacities for layer: " // 
c     &		trim(int2string(ir,'(i4)')) // " of " // trim(int2string(nr,'(i4)')))
c		call output("T = " // trim(dbl2string(T(ir),'(f8.2)')) // " K")
c		call output("P = " // trim(dbl2string(P(ir),'(es8.2)')) // " Ba")
		cont_tot=0d0
		do i=1,ncia
			if(T(ir).lt.CIA(i)%T(1)) then
				iT=1
			else if(T(iR).gt.CIA(i)%T(CIA(i)%nT)) then
				iT=CIA(i)%nT
			else
				do iT=1,CIA(i)%nT-1
					if(T(ir).ge.CIA(i)%T(iT).and.T(ir).le.CIA(i)%T(iT+1)) exit
				enddo
			endif
			cont_tot(1:nlam)=cont_tot(1:nlam)+CIA(i)%Cabs(iT,1:nlam)*Ndens(ir)*cia_mixrat(CIA(i)%imol1)*cia_mixrat(CIA(i)%imol2)
		enddo
		do imol=1,nmol
			kappa_mol(imol,1:nlam,1:ng)=0d0
			if(mixrat_r(ir,imol).gt.0d0) call ReadOpacityFITS(kappa_mol,imol,ir)
		enddo
c		call output("Compute k-tables")
		nu1=0d0
		nu2=1d0
		do i=1,n_nu_line
			nu_line(i)=1d0-real(i-1)/real(n_nu_line-1)
		enddo
		allocate(k_line(n_nu_line))
		do i=1,nlam-1
			nu1=-1d0
			nu2=2d0
			do j=1,n_nu_line
				k_line(j)=0d0
				do imol=1,nmol
					if(mixrat_r(ir,imol).gt.0d0) then
						ig=random(idum)*real(ng)+1
						k_line(j)=k_line(j)+kappa_mol(imol,i,ig)*mixrat_r(ir,imol)
					endif
				enddo
			enddo
			call ComputeKtable(ir,nu1,nu2,nu_line,k_line,n_nu_line,kappa,cont_tot(i))
			Cabs(ir,i,1:ng)=kappa(1:ng)
			Csca(ir,i)=8.4909d-45*(freq(i+1)*freq(i))**2
			do j=1,ng
				if(Cabs(ir,i,j).lt.Csca(ir,i)*1d-4) Cabs(ir,i,j)=Csca(ir,i)*1d-4 
			enddo
			opac_tot(i,1:ng)=opac_tot(i,1:ng)+(Cabs(ir,i,1:ng)+Csca(ir,i))*Ndens(ir)*(R(ir+1)-R(ir))
		enddo
		deallocate(k_line)
		if(outputopacity) then
			call WriteOpacity(ir,"ktab",freq,Cabs(ir,1:nlam-1,1:ng),nlam-1,ng)
			do i=1,nlam-1
				kaver(i)=0d0
				do j=1,ng
					kaver(i)=kaver(i)+Cabs(ir,i,j)/real(ng)
				enddo
			enddo
			call WriteOpacity(ir,"aver",freq,kaver(1:nlam-1),nlam-1,1)
		endif
	
c		call cpu_time(stoptime)

c		open(unit=30,file=trim(outputdir) // "opticaldepth" // trim(int2string(ir,'(i0.3)')) // ".dat",RECL=6000)
c		write(30,'("#",a13,a19)') "lambda [mu]","total average tau"
c		do i=1,nlam-1
c			write(30,'(f12.6,e19.7)') sqrt(lam(i)*lam(i+1))/micron,sum(opac_tot(i,1:ng))/real(ng)
c		enddo
c		close(unit=30)
	enddo

	open(unit=30,file=trim(outputdir) // "opticaldepth.dat",RECL=6000)
	write(30,'("#",a13,a19)') "lambda [mu]","total average tau"
	do i=1,nlam-1
		write(30,'(f12.6,e19.7)') sqrt(lam(i)*lam(i+1))/micron,sum(opac_tot(i,1:ng))/real(ng)
	enddo
	close(unit=30)

	return
	end




	subroutine ComputeOpacities()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer imol
	real*8 kappa(ng),nu1,nu2,tanscale
	real*8 x1,x2,rr,gasdev,random,dnu,Saver,starttime,stoptime
	real*8,allocatable :: k_line(:),nu_line(:),dnu_line(:)
	real*8,allocatable :: opac_tot(:,:),cont_tot(:),kaver(:)
	integer n_nu_line,iT
	integer i,j,ir,k,nl
	integer,allocatable :: inu1(:),inu2(:)
	character*500 filename

	allocate(cont_tot(nlam))
	allocate(kaver(nlam))
	allocate(opac_tot(nlam,ng))
	allocate(inu1(nlam))
	allocate(inu2(nlam))

	n_voigt=1d6
	allocate(a_therm(n_voigt))
	allocate(a_press(n_voigt))
	tanscale=atan(cutoff_lor)
	do i=1,n_voigt
		x1=gasdev(idum)/sqrt(2d0)
		x2=tan((random(idum)-0.5d0)*2d0*tanscale)
		rr=random(idum)
		if(rr.lt.0.25) then
			x1=-x1
		else if(rr.lt.0.5) then
			x2=-x2
		else if(rr.lt.0.75) then
			x1=-x1
			x2=-x2
		endif
		a_therm(i)=x1
		a_press(i)=x2
	enddo

	opac_tot=0d0

	call cpu_time(starttime)
	do ir=nr,1,-1
		call output("Opacities for layer: " // 
     &		trim(int2string(ir,'(i4)')) // " of " // trim(int2string(nr,'(i4)')))
		call output("T = " // trim(dbl2string(T(ir),'(f8.2)')) // " K")
		call output("P = " // trim(dbl2string(P(ir),'(es8.2)')) // " Ba")
		call LineStrengthWidth(ir,dnu,freq(nlam),freq(1))
		dnu=dnu/2d0
		n_nu_line=abs(freq(1)/freq(nlam))/dnu
		nu1=freq(1)
		n_nu_line=1
		do while(nu1.gt.freq(nlam))
			nu1=nu1-nu1*dnu
			n_nu_line=n_nu_line+1
		enddo
		allocate(k_line(n_nu_line))
		allocate(nu_line(n_nu_line))
		allocate(dnu_line(n_nu_line))
		do i=1,n_nu_line
			nu_line(i)=exp(log(freq(1))+log(freq(nlam)/freq(1))*real(i-1)/real(n_nu_line-1))
			dnu_line(i)=nu_line(i)*dnu
		enddo
		cont_tot=0d0
		do i=1,ncia
			if(T(ir).lt.CIA(i)%T(1)) then
				iT=1
			else if(T(iR).gt.CIA(i)%T(CIA(i)%nT)) then
				iT=CIA(i)%nT
			else
				do iT=1,CIA(i)%nT-1
					if(T(ir).ge.CIA(i)%T(iT).and.T(ir).le.CIA(i)%T(iT+1)) exit
				enddo
			endif
			cont_tot(1:nlam)=cont_tot(1:nlam)+CIA(i)%Cabs(iT,1:nlam)*Ndens(ir)*cia_mixrat(CIA(i)%imol1)*cia_mixrat(CIA(i)%imol2)
		enddo
		call output("Compute lines")
		call ComputeKline(ir,nu_line,k_line,n_nu_line,dnu_line)
		if(outputopacity) call WriteOpacity(ir,"line",nu_line,k_line,n_nu_line,1)
		call output("Compute k-tables")
		call tellertje(1,nlam-1)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,nu1,nu2,kappa,k)
!$OMP& SHARED(nlam,freq,ir,nu_line,k_line,n_nu_line,cont_tot,Cabs,Csca,opac_tot,Ndens,R,ng,nclouds,cloud_dens,Cloud)
!$OMP DO
		do i=1,nlam-1
			call tellertje(i+1,nlam+1)
			nu1=freq(i+1)
			nu2=freq(i)
			call ComputeKtable(ir,nu1,nu2,nu_line,k_line,n_nu_line,kappa,cont_tot(i))
			Cabs(ir,i,1:ng)=kappa(1:ng)
			Csca(ir,i)=8.4909d-45*(nu1*nu2)**2
			do j=1,ng
				if(Cabs(ir,i,j).lt.Csca(ir,i)*1d-4) Cabs(ir,i,j)=Csca(ir,i)*1d-4 
			enddo
			opac_tot(i,1:ng)=opac_tot(i,1:ng)+(Cabs(ir,i,1:ng)+Csca(ir,i))*Ndens(ir)*(R(ir+1)-R(ir))
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		if(outputopacity) then
			call WriteOpacity(ir,"ktab",freq,Cabs(ir,1:nlam-1,1:ng),nlam-1,ng)
			do i=1,nlam-1
				kaver(i)=0d0
				do j=1,ng
					kaver(i)=kaver(i)+Cabs(ir,i,j)/real(ng)
				enddo
			enddo
			call WriteOpacity(ir,"aver",freq,kaver(1:nlam-1),nlam-1,1)
		endif
		call tellertje(nlam,nlam)
	
		open(unit=30,file=trim(outputdir) // "opticaldepth.dat",RECL=6000)
		write(30,'("#",a13,a19)') "lambda [mu]","total average tau"
		do i=1,nlam-1
			write(30,'(f12.6,e19.7)') sqrt(lam(i)*lam(i+1))/micron,sum(opac_tot(i,1:ng))/real(ng)
		enddo
		close(unit=30)

		deallocate(k_line)
		deallocate(nu_line)
		deallocate(dnu_line)
		call cpu_time(stoptime)
		call output("Time for this layer: " // trim(dbl2string((stoptime-starttime),'(f10.2)')) // " s")
		call output("==================================================================")
		starttime=stoptime
		if(minval(opac_tot(1:nlam-1,1:ng)).gt.maxtau.and.maxtau.gt.0d0) then
			call output("Maximum optical depth reached at all wavelengths")
			call output("ignoring layers: 1 to " // trim(int2string(ir-1,'(i4)')))
			if(ir.gt.1) then
				do i=1,ir-1
					Cabs(i,1:nlam-1,1:ng)=Cabs(ir,1:nlam-1,1:ng)
					Csca(i,1:nlam-1)=Csca(ir,1:nlam-1)
				enddo
			endif
			exit
		endif
	enddo

	deallocate(a_therm)
	deallocate(a_press)

	return
	end
	
	subroutine LineStrengthWidth(ir,minw,nu1,nu2)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 w,x1,x2,minw,nu1,nu2,Saver(nlam),gamma,Saver0(nlam)
	integer imol,iT,i,ir,iiso,nl,nl0,ilam,nlines_lam(nlam)

	call output("Line strengths and widths")
	
	call hunt(TZ,nTZ,T(ir),iT)
	if(iT.lt.1) iT=1
	if(iT.gt.nTZ) iT=nTZ

	call tellertje(1,nlines)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,imol,iiso,x1,x2)
!$OMP& SHARED(nlines,nu1,nu2,minw,iT,Mmol,P,ir,ZZ,mixrat_r,T,
!$OMP&     L_do,L_S,L_S0,L_freq,L_imol,L_iiso,L_Elow)
!$OMP DO
	do i=1,nlines
		call tellertje(i+1,nlines+2)
		L_do(i)=.false.
		if((L_freq(i)).gt.nu1.and.(L_freq(i)).lt.nu2) then
			L_do(i)=.true.
			imol=L_imol(i)
			iiso=L_iiso(i)
c line strength
			x1=exp(-hplanck*clight*L_Elow(i)/(kb*T(ir)))
			x2=exp(-hplanck*clight*L_freq(i)/(kb*T(ir)))

			L_S(i)=L_S0(i)*(x1*(1d0-x2))*mixrat_r(ir,imol)/ZZ(imol,iiso,iT)
		endif
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	minw=sqrt(minw)
	call tellertje(nlines,nlines)

	Saver=0d0
	nl=nlines*2

1	continue
	nl0=nl
	Saver0=Saver*eps_lines
	Saver=0d0
	nl=0
	nlines_lam=0
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i)
!$OMP& SHARED(L_do,nlines,L_ilam,L_S,Saver0,Saver,nlines_lam)
!$OMP DO
	do i=1,nlines
		if(L_do(i)) then
			if(L_S(i).gt.Saver0(L_ilam(i))) then
				Saver(L_ilam(i))=Saver(L_ilam(i))+L_S(i)
				nlines_lam(L_ilam(i))=nlines_lam(L_ilam(i))+1
			else
				L_do(i)=.false.
			endif
		endif
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	nl=sum(nlines_lam(1:nlam))

	Saver=Saver/real(nlines_lam)
	if(real(nl0).gt.(real(nl)*1.1)) goto 1

	minw=0.01d0
	minw=(1d0/specres)**2
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,imol,w,gamma,x1,x2)
!$OMP& SHARED(nlines,minw,iT,Mmol,P,ir,T,Saver,L_Saver,L_ilam,nlines_lam,L_nclose,mixrat_r,
!$OMP&     L_do,L_gamma_air,L_gamma_self,L_a_therm,L_a_press,L_freq,L_imol,L_n)
!$OMP DO
	do i=1,nlines
		if(L_do(i)) then
			imol=L_imol(i)
c thermal broadening
			w=(sqrt(2d0*kb*T(ir)/(Mmol(imol)*mp)))
			L_a_therm(i)=w*L_freq(i)/clight
c pressure broadening
			x1=L_gamma_air(i)*(1d0-mixrat_r(ir,imol))
			x2=L_gamma_self(i)*mixrat_r(ir,imol)
			L_a_press(i)=((x1+x2)*P(ir)/atm)*(296d0/T(ir))**L_n(i)

			gamma=((L_a_press(i)*4d0)**2+L_a_therm(i)**2)/L_freq(i)**2
			if(gamma.lt.minw) then
				minw=gamma
			endif
			L_nclose(i)=nlines_lam(L_ilam(i))
			L_Saver(i)=Saver(L_ilam(i))
		endif
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	minw=sqrt(minw)

	call output("number of lines: " // trim(dbl2string(dble(nl),'(es7.1)')))
	
	return
	end
	
	
	subroutine ComputeKtable(ir,nu1,nu2,nu,kline,nnu,kappa,Ccont)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nnu,nkdis
	parameter(nkdis=100)
	real*8 nu1,nu2,kappa(ng),g(ng),w,dnu,gamma,fact,eps
	real*8 nu(nnu),kline(nnu),Ccont,kap
	real*8 kdis(nkdis),dis(nkdis),gb(ng+1)
	real*8 Eu,El,A,x,kmax,kmin,V,scale,x1,x2,gasdev,random,rr,gu,gl
	integer inu,iT,imol,i,ju,jl,j,NV,nl,k,iiso,ir,NV0,iter
	integer i_therm,i_press,inu1,inu2,ig(ng+1)
	logical converged
	real*8 f,a_t,a_p
	
	do i=1,ng
		g(i)=(real(i)-0.5)/real(ng)
	enddo
	do i=1,ng+1
		gb(i)=real(i-1)/real(ng)
	enddo

	call hunt(nu,nnu,nu1,inu2)
	call hunt(nu,nnu,nu2,inu1)
	inu1=inu1+1
	if(inu1.gt.nnu) inu1=nnu
	if(inu2.gt.nnu) inu2=nnu
	if(inu1.lt.2) inu1=2
	if(inu2.lt.2) inu2=2

	kmin=1d200
	kmax=Ccont
	do inu=inu1,inu2
		kap=kline(inu)+Ccont
		if(kap.gt.kmax) kmax=kap
		if(kap.lt.kmin) kmin=kap
	enddo
	if(kmax.eq.0d0) then
		kappa=0d0
		return
	endif
	if(kmin.lt.kmax/1d10) kmin=kmax/1d10
	kmin=log10(kmin)
	kmax=log10(kmax)

	do i=1,nkdis
		kdis(i)=10d0**(kmin+real(i-1)*(kmax-kmin)/real(nkdis-1))
	enddo
	dis=0d0
	do inu=inu1,inu2
		do i=nkdis,1,-1
			if((kline(inu)+Ccont).le.kdis(i)) then
				dis(i)=dis(i)+1d0
			else
				exit
			endif
		enddo
	enddo

	dis(1:nkdis-1)=dis(1:nkdis-1)/dis(nkdis)
	dis(nkdis)=1d0

	do i=1,ng+1
		call hunt(dis,nkdis,gb(i),ig(i))
		if(ig(i).lt.1) ig(i)=1
		if(ig(i).gt.nkdis) ig(i)=nkdis
	enddo
	do i=1,ng
		kappa(i)=sqrt(kdis(ig(i))*kdis(ig(i+1)))
	enddo

	return
	end



	subroutine ComputeKline(ir,nu,kline,nnu,dnu)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nnu
	real*8 w,gamma,fact
	real*8 nu(nnu),dnu(nnu)
	real*8 kline(nnu)
	real*8 Eu,El,A,x,kmax,kmin,V,scale,x1,x2,gasdev,random,rr,gu,gl
	integer iT,imol,i,ju,jl,j,nkdis,NV,k,iiso,ir,NV0,iter
	integer i_therm,i_press,il,idnu,inu1,inu2,inu
	real*8 f,a_t,a_p
	
	fact=50d0
	kline=0d0

	scale=real(nnu-1)/log(nu(nnu)/nu(1))

c	NV0=real(nnu)*100d0/real(nl+1)+250d0
c	if(NV0.gt.25000) NV0=25000

	call hunt(TZ,nTZ,T(ir),iT)

	call tellertje(1,nlines)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,imol,gamma,A,a_t,a_p,f,x1,x2,rr,x,inu,i_press,i_therm,idnu,inu1,inu2,NV,NV0)
!$OMP& SHARED(fact,mixrat_r,scale,kline,nnu,nu,nlines,a_therm,a_press,n_voigt,P,ir,cutoff_abs,
!$OMP&     L_do,L_S,L_a_therm,L_a_press,L_freq,L_imol,L_nclose,nlam,L_Saver)
!$OMP DO SCHEDULE (STATIC,1)
	do i=1,nlines
		call tellertje(i+1,nlines+2)
		if(L_do(i)) then
			imol=L_imol(i)
			A=L_S(i)
			a_t=L_a_therm(i)
			a_p=L_a_press(i)
			gamma=sqrt(a_t**2+a_p**2)
			f=L_freq(i)
c	Random sampling of the Voigt profile
			NV0=real(nnu)*50d0/real(nlam*L_nclose(i)+1)+100d0
			if(NV0.gt.20000) NV0=20000
			NV=real(NV0)*A/L_Saver(i)
			if(NV.gt.100*NV0) NV=100*NV0
			if(NV.lt.25) NV=25

			A=A/real(NV)
			i_therm=random(idum)*real(n_voigt)
			i_press=random(idum)*real(n_voigt)
			do j=1,NV
				i_therm=i_therm+1
				i_press=i_press+1
				if(i_therm.gt.n_voigt) i_therm=1
				if(i_press.gt.n_voigt) i_press=1
				x1=a_therm(i_therm)
				x1=x1*a_t
				x2=a_press(i_press)
				x2=x2*a_p
				x=x1+x2
				if(abs(x).lt.cutoff_abs) then
					idnu=abs(x/gamma)
					x=x+f
					x=log(x/nu(1))*scale
					inu=int(x)+1
					inu1=inu-idnu
					inu2=inu+idnu
					if(inu1.le.nnu.and.inu2.ge.1) then
						if(inu1.lt.1) inu1=1
						if(inu2.gt.nnu) inu2=nnu
						kline(inu1:inu2)=kline(inu1:inu2)+A/real(inu2-inu1+1)
					endif
				endif
			enddo
		endif
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(nlines,nlines)

	kline=kline/dnu

	return
	end



