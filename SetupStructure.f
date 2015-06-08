	subroutine SetupStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 g,dp,dz,dlogp,RgasBar,sh
	parameter(RgasBar=82.05736*1.01325)
	integer i,imol

	R(1)=Rplanet
	mu=1d0
	do imol=1,nmol
		if(mixrat(imol).gt.0d0) mu=mu-mixrat(imol)
	enddo
	mu=mu*2.0
	do imol=1,nmol
		if(mixrat(imol).gt.0d0) mu=mu+mixrat(imol)*Mmol(imol)
	enddo
	call output("Mean molecular weight: " // dbl2string(mu,'(f8.3)'))
	do i=1,nr
		mu=1d0
		do imol=1,nmol
			if(mixrat_r(i,imol).gt.0d0) mu=mu-mixrat_r(i,imol)
		enddo
		mu=mu*2.0
		do imol=1,nmol
			if(mixrat_r(i,imol).gt.0d0) mu=mu+mixrat_r(i,imol)*Mmol(imol)
		enddo
c		call output("Layer: " // 
c     &		trim(int2string(i,'(i4)')) // " of " // trim(int2string(nr,'(i4)')))
c		call output("Mean molecular weight: " // dbl2string(mu,'(f8.3)'))

		g=Ggrav*Mplanet/R(i)**2
		
		if(i.eq.nr) then
			dp=P(nr-1)-P(nr)
			dlogp=log(P(nr-1)/P(nr))
		else
			dp=(P(i)-P(i+1))
			dlogp=log(P(i)/P(i+1))
		endif
		Ndens(i)=P(i)*Avogadro/(RgasBar*T(i))
		dens(i)=Ndens(i)*mp*mu
		sh=(T(i)*kb)/(g*mp*mu)
		dz=dp/(dens(i)*g)
		dz=dlogp*sh
		R(i+1)=R(i)+dz
	enddo

	do i=1,nclouds
		call SetupCloud(i)
	enddo

	open(unit=50,file=trim(outputdir) // 'densityprofile.dat',RECL=100)
	write(50,'("#",a14,a15,a15,a13,a10,a10)') "radius [cm]","height [cm]","dens [g/cm^3]","N [1/cm^3]","T [K]","P [Ba]"
	do i=1,nr
		write(50,'(es15.7,es15.4,es15.4,es13.4,f10.3,es10.3)') sqrt(R(i)*R(i+1)),sqrt(R(i)*R(i+1))-Rplanet
     &			,dens(i),Ndens(i),T(i),P(i)
	enddo
	close(unit=20)

	return
	end
	

	subroutine SetupCloud(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 h,dh,findpressure,column,rr,tot
	integer ii,i,j,nsubr
	
	if(Cloud(ii)%H.lt.0d0) then
		if(Cloud(ii)%P.gt.0d0) then
			h=findpressure(Cloud(ii)%P)
			if(Cloud(ii)%dP.gt.0d0) then
				dh=abs(findpressure(Cloud(ii)%P*(1d0+Cloud(ii)%dP))-h)
			else if(Cloud(ii)%dH.gt.0d0) then
				dh=Cloud(ii)%dH
			else
				call output("Cloud thickness not set.")
				stop
			endif
		else
			call output("Cloud height not set.")
			stop
		endif
	else
		h=Cloud(ii)%H
		dh=Cloud(ii)%dH
	endif
	
	column=0d0
	nsubr=100
	do i=1,nr
		do j=1,nsubr
			rr=10d0**(log10(R(i))+log10(R(i+1)/R(i))*real(j)/real(nsubr+1))
			cloud_dens(i,ii)=exp(-((rr-h)/dh)**2/2d0)/real(nsubr)
			column=column+cloud_dens(i,ii)*(R(i+1)-R(i))
		enddo
	enddo

	cloud_dens(1:nr,ii)=cloud_dens(1:nr,ii)*Cloud(ii)%column/column

	tot=0d0
	do i=1,Cloud(ii)%nsize
		Cloud(ii)%w(i)=(1q4*Cloud(ii)%rv(i))**(1q0+(1q0-3q0*Cloud(ii)%veff)/Cloud(ii)%veff)*
     &					qexp(-1q4*Cloud(ii)%rv(i)/(Cloud(ii)%reff*Cloud(ii)%veff))
		Cloud(ii)%w(i)=Cloud(ii)%w(i)*Cloud(ii)%rv(i)**3
		tot=tot+Cloud(ii)%w(i)
	enddo
	Cloud(ii)%w=Cloud(ii)%w/tot

	return
	end	


	real*8 function findpressure(P0)
	use GlobalSetup
	IMPLICIT NONE
	real*8 P0,r1,r2
	integer i

	do i=1,nr-1
		if(P(i).gt.P0.and.P(i+1).lt.P0) then
			r1=log10(sqrt(R(i)*R(i+1)))
			r2=log10(sqrt(R(i+1)*R(i+2)))
			findpressure=10d0**(r1+(r2-r1)*(P(i)-P0)/(P(i)-P(i+1)))
			return
		endif
	enddo
	call output("Error in finding pressure...")
	findpressure=R(nr+1)
	
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine SetupPartCloud(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,is,ilam,j
	real*8 phi,thet,tot,tot2,fact
	
	allocate(Cloud(ii)%rv(Cloud(ii)%nsize))
	allocate(Cloud(ii)%Kabs(Cloud(ii)%nsize,nlam))
	allocate(Cloud(ii)%Ksca(Cloud(ii)%nsize,nlam))
	allocate(Cloud(ii)%Kext(Cloud(ii)%nsize,nlam))
	allocate(Cloud(ii)%F(Cloud(ii)%nsize,nlam))

	select case(Cloud(ii)%ptype)
		case("COMPUTE")
			do is=1,Cloud(ii)%nsize
				call ComputePart(Cloud(ii),ii,is)
			enddo
c		case("PARTFILE")
c			call ReadParticle(Cloud(ii),ii)
		case default
			call output("I did not understand what I was trying to do. Sorry!")
	end select
	
	if(nspike.gt.0.and.nspike.le.180) call output("making the first " // trim(int2string(nspike,'(i3)')) // " degrees isotropic")
	do ilam=1,nlam
		do is=1,Cloud(ii)%nsize
			tot=0d0
			tot2=0d0
			do j=1,180
				tot=tot+Cloud(ii)%F(is,ilam)%F11(j)*sin(pi*(real(j)-0.5)/180d0)
				tot2=tot2+sin(pi*(real(j)-0.5)/180d0)
			enddo
			do j=1,180
				Cloud(ii)%F(is,ilam)%F11(j)=tot2*Cloud(ii)%F(is,ilam)%F11(j)/tot
				Cloud(ii)%F(is,ilam)%F12(j)=tot2*Cloud(ii)%F(is,ilam)%F12(j)/tot
				Cloud(ii)%F(is,ilam)%F22(j)=tot2*Cloud(ii)%F(is,ilam)%F22(j)/tot
				Cloud(ii)%F(is,ilam)%F33(j)=tot2*Cloud(ii)%F(is,ilam)%F33(j)/tot
				Cloud(ii)%F(is,ilam)%F34(j)=tot2*Cloud(ii)%F(is,ilam)%F34(j)/tot
				Cloud(ii)%F(is,ilam)%F44(j)=tot2*Cloud(ii)%F(is,ilam)%F44(j)/tot
			enddo

			if(nspike.gt.0.and.nspike.le.180) then
c the nspike parameter removes the n degree spike in the forward direction.
				do j=1,nspike
					fact=Cloud(ii)%F(is,ilam)%F11(nspike+1)/Cloud(ii)%F(is,ilam)%F11(j)
					Cloud(ii)%F(is,ilam)%F12(j)=Cloud(ii)%F(is,ilam)%F12(j)*fact
					Cloud(ii)%F(is,ilam)%F22(j)=Cloud(ii)%F(is,ilam)%F22(j)*fact
					Cloud(ii)%F(is,ilam)%F33(j)=Cloud(ii)%F(is,ilam)%F33(j)*fact
					Cloud(ii)%F(is,ilam)%F34(j)=Cloud(ii)%F(is,ilam)%F34(j)*fact
					Cloud(ii)%F(is,ilam)%F44(j)=Cloud(ii)%F(is,ilam)%F44(j)*fact
					Cloud(ii)%F(is,ilam)%F11(j)=Cloud(ii)%F(is,ilam)%F11(nspike+1)
				enddo

				tot=0d0
				tot2=0d0
				do j=1,180
					tot=tot+Cloud(ii)%F(is,ilam)%F11(j)*sin(pi*(real(j)-0.5)/180d0)
					tot2=tot2+sin(pi*(real(j)-0.5)/180d0)
				enddo
				Cloud(ii)%Ksca(is,ilam)=Cloud(ii)%Ksca(is,ilam)*tot/tot2
				Cloud(ii)%Kext(is,ilam)=Cloud(ii)%Kabs(is,ilam)+Cloud(ii)%Ksca(is,ilam)
				do j=1,180
					Cloud(ii)%F(is,ilam)%F11(j)=tot2*Cloud(ii)%F(is,ilam)%F11(j)/tot
					Cloud(ii)%F(is,ilam)%F12(j)=tot2*Cloud(ii)%F(is,ilam)%F12(j)/tot
					Cloud(ii)%F(is,ilam)%F22(j)=tot2*Cloud(ii)%F(is,ilam)%F22(j)/tot
					Cloud(ii)%F(is,ilam)%F33(j)=tot2*Cloud(ii)%F(is,ilam)%F33(j)/tot
					Cloud(ii)%F(is,ilam)%F34(j)=tot2*Cloud(ii)%F(is,ilam)%F34(j)/tot
					Cloud(ii)%F(is,ilam)%F44(j)=tot2*Cloud(ii)%F(is,ilam)%F44(j)/tot
				enddo
			endif
		enddo
	enddo

	do ilam=1,nlam
		do is=1,Cloud(ii)%nsize
			Cloud(ii)%F(is,ilam)%IF11=0d0
			Cloud(ii)%F(is,ilam)%IF12=0d0
			do j=1,180
				thet=pi*(real(j)-0.5d0)/180d0
				Cloud(ii)%F(is,ilam)%IF11=Cloud(ii)%F(is,ilam)%IF11+pi*sin(thet)
     &			*Cloud(ii)%F(is,ilam)%F11(j)/180d0
				Cloud(ii)%F(is,ilam)%IF12=Cloud(ii)%F(is,ilam)%IF12+pi*sin(thet)
     &			*Cloud(ii)%F(is,ilam)%F12(j)/180d0
			enddo
		enddo
	enddo
c	do j=0,360
c		phi=pi*real(j-1)/179.5d0
c		cos2phi(j)=cos(2d0*phi)
c		sin2phi(j)=sin(2d0*phi)
c	enddo
	
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	