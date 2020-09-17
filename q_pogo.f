	module NeuralOpacity
	use mod_network, only: network_type
	use mod_kinds, only: ik, rk
	IMPLICIT NONE
	
	type(network_type) :: net_abs,net_sca
	real*8 minabs,maxabs,minsca,maxsca
	
	end module NeuralOpacity
	
	subroutine InitNeural(dirabs,dirsca)
	use NeuralOpacity
	IMPLICIT NONE
	character*100 dirabs,dirsca
	
	call net_abs%load(trim(dirabs) // '/network.txt')
	open(unit=20,file=trim(dirabs) // '/minmax.txt')
	read(20,*) minabs,maxabs
	close(unit=20)
	
	call net_sca%load(trim(dirsca) // '/network.txt')
	open(unit=20,file=trim(dirsca) // '/minmax.txt')
	read(20,*) minsca,maxsca
	close(unit=20)
	
	return
	end

	
	subroutine q_pogo(e1,e2,lam,rad,cext,csca,g,maxf)
	use NeuralOpacity
	IMPLICIT NONE
	real*8 e1,e2,lam,rad,cext,csca,g,maxf,cabs,x,xx,pi
	parameter(pi=3.1415926536)
	real(rk) in(4),out(1)

	g=0d0

	in(1)=log10(e1)
	in(2)=log10(e2)
	x=2d0*pi*rad/lam
	in(3)=log10(x)
	in(4)=maxf

	if(in(3).gt.3.0) in(3)=3.0
	if(in(3).lt.-2.0) in(3)=-2.0

	xx=10d0**in(3)
	
	out=net_abs%output(in)
	cabs=out(1)
	cabs=cabs*(maxabs-minabs)+minabs
	cabs=10d0**cabs
	cabs=cabs*xx**2/(1d0+1d0/xx)
	if(x.gt.1d3) cabs=cabs*((x**2)/1d6)
	if(x.lt.1d-2) cabs=cabs*((x**3)/1d-6)

	out=net_sca%output(in)
	csca=out(1)
	csca=csca*(maxsca-minsca)+minsca
	csca=10d0**csca
	csca=csca*(xx**2)/(1d0+1d0/xx**4)
	if(x.gt.1d3) csca=csca*((x**2)/1d6)
	if(x.lt.1d-2) csca=csca*((x**6)/1d-12)

	csca=csca*(lam/(2d0*pi))**2
	cabs=cabs*(lam/(2d0*pi))**2

	cext=csca+cabs

	if(cabs.lt.cext*1d-4) then
		cabs=cext*1d-4
		cext=cabs+csca
	endif


	return
	end
	