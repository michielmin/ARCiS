	subroutine SetupOpacities()
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
		call LineStrengthWidth(ir,dnu,Saver,nl,freq(nlam),freq(1))
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
		call ComputeKline(ir,nu_line,k_line,n_nu_line,dnu_line,Saver,nl)
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
		if(minval(opac_tot(1:nlam-1,1:ng)).gt.maxtau) then
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

	return
	end
	
	subroutine LineStrengthWidth(ir,minw,Saver,nl,nu1,nu2)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 w,x1,x2,x3,x4,minw,nu1,nu2,Saver,gamma,Saver0
	integer imol,iT,i,ir,iiso,nl,nl0

	call output("Line strengths and widths")
	
	minw=0.01d0
	minw=(1d0/specres)**2
	call hunt(TZ,nTZ,T(ir),iT)
	if(iT.lt.1) iT=1
	if(iT.gt.nTZ) iT=nTZ

	Saver=0d0
	nl=0
	call tellertje(1,nlines)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,imol,iiso,w,x1,x2,x3,x4,gamma)
!$OMP& SHARED(nlines,nu1,nu2,minw,iT,Mmol,P,ir,ZZ,mixrat_r,T,Saver,nl,
!$OMP&     L_do,L_S,L_S0,L_gamma_air,L_gamma_self,L_a_therm,L_a_press,L_freq,L_imol,L_iiso,L_Elow,L_n)
!$OMP DO
	do i=1,nlines
		call tellertje(i+1,nlines+2)
		L_do(i)=.false.
		if((L_freq(i)).gt.nu1.and.(L_freq(i)).lt.nu2) then
			L_do(i)=.true.
			imol=L_imol(i)
			iiso=L_iiso(i)
c thermal broadening
			w=(sqrt(2d0*kb*T(ir)/(Mmol(imol)*mp)))
			L_a_therm(i)=w*L_freq(i)/clight
c pressure broadening
			L_a_press(i)=L_gamma_air(i)*P(ir)*(1d0-mixrat_r(ir,imol))/atm
			L_a_press(i)=L_a_press(i)+L_gamma_self(i)*P(ir)*mixrat_r(ir,imol)/atm
			L_a_press(i)=L_a_press(i)*(296d0/T(ir))**L_n(i)
c line strength
			x1=exp(-hplanck*clight*L_Elow(i)/(kb*T(ir)))
			x2=exp(-hplanck*clight*L_freq(i)/(kb*T(ir)))
c			x3=exp(-hplanck*clight*L%Elow/(kb*296d0))
c			x4=exp(-hplanck*clight*L%freq/(kb*296d0))

c			L%S=L%S0*(x1*(1d0-x2))/(x3*ZZ(imol,iiso,iT)*(1d0-x4))
			L_S(i)=L_S0(i)*(x1*(1d0-x2))/ZZ(imol,iiso,iT)
			L_S(i)=L_S(i)*mixrat_r(ir,imol)

			gamma=((L_a_press(i)*4d0)**2+L_a_therm(i)**2)/L_freq(i)**2
			if(gamma.lt.minw) then
				minw=gamma
			endif
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
	do i=1,nlines
		if(L_do(i)) then
			imol=L_imol(i)
			if(L_S(i).gt.Saver0) then
				nl=nl+1
				Saver=Saver+L_S(i)
			else
				L_do(i)=.false.
			endif
		endif
	enddo

	Saver=Saver/real(nl)
	if(real(nl0).gt.(real(nl)*1.1)) goto 1

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
	integer inu,iT,imol,i,ju,jl,j,NV,nl,k,iiso,ir,NV0,iter,maxiter
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

	kmin=1d200
	kmax=Ccont
	do inu=inu1+1,inu2-1
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



	subroutine ComputeKline(ir,nu,kline,nnu,dnu,Saver,nl)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nnu
	real*8 w,gamma,fact
	real*8 nu(nnu),dnu(nnu)
	real*8,target :: kline(nnu)
	real*8 Eu,El,A,x,kmax,kmin,V,scale,x1,x2,gasdev,random,rr,gu,gl,Saver
	integer iT,imol,i,ju,jl,j,nkdis,NV,nl,k,iiso,ir,NV0,iter,maxiter
	integer i_therm,i_press,il,idnu,inu1,inu2,inu
	real*8 f,a_t,a_p
	
	fact=50d0
	kline=0d0

	scale=real(nnu-1)/log(nu(nnu)/nu(1))

	NV0=real(nnu)*100d0/real(nl+1)+250d0
	
	call hunt(TZ,nTZ,T(ir),iT)

	call tellertje(1,nlines)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,imol,gamma,A,a_t,a_p,f,x1,x2,rr,x,inu,i_press,i_therm,idnu,inu1,inu2,NV)
!$OMP& SHARED(fact,mixrat_r,scale,NV0,kline,nnu,nu,nlines,a_therm,a_press,n_voigt,P,ir,cutoff_abs,Saver,
!$OMP&     L_do,L_S,L_a_therm,L_a_press,L_freq,L_imol)
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
			NV=real(NV0)*A/Saver
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



