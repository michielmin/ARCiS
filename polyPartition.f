	subroutine polyPartition(imol,T,Q)
	IMPLICIT NONE
	real*8 a(0:5,70),lnT,T,Q
	integer imol,i

	if(imol.lt.56.or.imol.gt.58.and.imol.ne.48) then
		print*,'Species not tabulated'
		return
	endif

	i=48	! He
	a(0,i)=-3.76575219d-1
	a(1,i)=2.33951687d-1
	a(2,i)=-5.79755525d-2
	a(3,i)=7.16333160d-3
	a(4,i)=-4.41302573d-4
	a(5,i)=1.08442997d-5
	
	i=56	! Na
	a(0,i)=-2.60507178d3
	a(1,i)=1.71419244d3
	a(2,i)=-4.50632658d2
	a(3,i)=5.91751503d1
	a(4,i)=-3.88164070d0
	a(5,i)=1.01752936d-1
	
	i=57	! K
	a(0,i)=-1.05889824d3
	a(1,i)=7.37734975d2
	a(2,i)=-2.05032627d2
	a(3,i)=2.84353380d1
	a(4,i)=-1.96809457d0
	a(5,i)=5.43892001d-2
	
	i=58	! TiO
	a(0,i)=7.93847646d1
	a(1,i)=-3.32551562d1
	a(2,i)=4.67690671d0
	a(3,i)=-6.24510956d-2
	a(4,i)=-3.08185663d-2
	a(5,i)=1.72675950d-3
	
	Q=0d0
	lnT=log(T)
	do i=0,5
		Q=Q+a(i,imol)*lnT**real(i)
	enddo
	Q=exp(Q)
	
	return
	end
	