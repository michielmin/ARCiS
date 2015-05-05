	subroutine Raytrace(iobs)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs
	real*8 rr,xx1,xx2,si,exp_tau,A,d,s,fluxg,Planck,fact,tau,freq0
	real*8,allocatable :: rtrace(:)
	integer nrtrace,ndisk,i,ir,ir_next,ilam,ig
	logical in
	
	ndisk=10
	
	nrtrace=nr+ndisk-1
	allocate(rtrace(nrtrace))
	
	do i=1,ndisk
		rtrace(i)=Rplanet*real(i-1)/real(ndisk-1)
	enddo
	do i=1,nr
		rtrace(i+ndisk-1)=R(i+1)
	enddo

	do ilam=1,nlam-1
		freq0=sqrt(freq(ilam)*freq(ilam+1))
		obs(iobs)%lam(ilam)=sqrt(lam(ilam)*lam(ilam+1))
		obs(iobs)%flux(ilam)=0d0
		do ig=1,ng
			fluxg=0d0
			do i=1,nrtrace-1
				fact=1d0
				A=pi*(rtrace(i+1)**2-rtrace(i)**2)
				rr=sqrt(rtrace(i)*rtrace(i+1))
				ir=nr
				si=-1d0
				xx1=si*sqrt(R(ir+1)**2-rr**2)
				in=.true.
1				continue
				if(in) then
					xx2=(R(ir)**2-rr**2)
					if(xx2.gt.0d0) then
						xx2=si*sqrt(xx2)
						d=abs(xx1-xx2)
						ir_next=ir-1
						goto 2
					else
						si=-si
						xx2=si*sqrt(R(ir+1)**2-rr**2)
						d=abs(xx1-xx2)
						ir_next=ir+1
						in=.false.
						goto 2
					endif
				else
					xx2=si*sqrt(R(ir+1)**2-rr**2)
					d=abs(xx1-xx2)
					ir_next=ir+1
					goto 2
				endif
2				continue
				tau=d*Ndens(ir)*opac(ir,ilam,ig)
				exp_tau=exp(-tau)
				fluxg=fluxg+A*Planck(T(ir),freq0)*(1d0-exp_tau)*fact
				fact=fact*exp_tau
				if(ir_next.gt.0.and.ir_next.le.nr) then
					ir=ir_next
					xx1=xx2
					goto 1
				endif
			enddo
		enddo
		obs(iobs)%flux(ilam)=obs(iobs)%flux(ilam)+fluxg/real(ng)
	enddo
	
	deallocate(rtrace)
	
	return
	end
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	real*8 function Planck(T,nu)
	use Constants
	IMPLICIT NONE
	real*8 T,nu,x

	x=hplanck*nu*clight/(kb*T)
	Planck=(2d0*hplanck*nu**3/clight**2)/(exp(x)-1d0)

	return
	end


