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

	open(unit=50,file=trim(outputdir) // 'densityprofile.dat',RECL=100)
	write(50,'("#",a14,a15,a15,a13,a10,a10)') "radius [cm]","height [cm]","dens [g/cm^3]","N [1/cm^3]","T [K]","P [Ba]"
	do i=1,nr
		write(50,'(es15.7,es15.4,es15.4,es13.4,f10.3,es10.3)') sqrt(R(i)*R(i+1)),sqrt(R(i)*R(i+1))-Rplanet
     &			,dens(i),Ndens(i),T(i),P(i)
	enddo
	close(unit=20)
	
	
	return
	end
	
	

