	subroutine SetupStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 g,dp,dz,dlogp,RgasBar,sh
	parameter(RgasBar=82.05736*1.01325)
	integer i

	g=Ggrav*Mplanet/Rplanet**2
	R(1)=Rplanet
	do i=1,nr
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
		write(50,'(e15.7,e15.4,e15.4,e13.4,f10.3,e10.3)') sqrt(R(i)*R(i+1)),sqrt(R(i)*R(i+1))-Rplanet
     &			,dens(i),Ndens(i),T(i),P(i)
	enddo
	close(unit=20)
	
	
	return
	end
	
	
