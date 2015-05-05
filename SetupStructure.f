	subroutine SetupStructure()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 g,dp,dz
	integer i

	g=Ggrav*Mplanet/Rplanet**2
	R(1)=Rplanet
	do i=1,nr
		if(i.eq.1) then
			dp=P(1)-P(2)
		else if(i.eq.nr) then
			dp=P(nr-1)-P(nr)
		else
			dp=(P(i-1)-P(i+1))/2d0
		endif
		Ndens(i)=P(i)/(kb*T(i))
		dens(i)=Ndens(i)*mp*mu
		dz=dp/(dens(i)*g)
		R(i+1)=R(i)+dz
	enddo

	open(unit=50,file=trim(outputdir) // 'densityprofile.dat',RECL=1000)
	do i=1,nr
		write(50,*) sqrt(R(i)*R(i+1)),dens(i),Ndens(i),T(i),P(i)
	enddo
	close(unit=20)
	
	
	return
	end
	
	
