	subroutine Genetic(geneticfun,var0,dvar0,nvars,nobs,npop,ngen,idum,gene_cross)
	IMPLICIT NONE
	integer i,j,k,l,nvars,nobs,npop,ngen,idum,ii,jj,jj1,jj2,ig,idupl,i1,i2,nmate,icreep,ncomp
	integer mut_dir,imut,nclones
	real*8,allocatable :: var(:,:,:),fit(:,:),fitobs(:,:,:),w(:,:),mate(:,:),x(:),y(:)
	real*8,allocatable :: ym(:),yf(:),male(:),female(:)
	real*8 ran1,r,var0(nvars),dvar0(2,nvars)
	real*8 av,dav1,dav2,tot,max,mutate,mut_min,mut_max,dCor,gasdev
	real*8 geneticfun,dfit,dfit2
	real*8 r1(nvars),r2(nvars),error(nvars)
	external geneticfun
	character*8,allocatable :: gene(:,:,:)
	character*1 v
	real*8 ra1,ra2,v1,v2,xx
	integer,allocatable :: same_pop(:),same_gen(:),ipar(:)
	logical,allocatable :: same(:)
	logical precomputed,creep,gene_cross
	logical,allocatable :: finished(:)

	allocate(var(ngen+1,0:npop,nvars))
	allocate(fit(ngen,0:npop))
	allocate(fitobs(ngen,0:npop,nobs))
	allocate(w(nobs,nvars))
	allocate(mate(0:npop,0:nobs))
	allocate(male(0:npop))
	allocate(female(0:npop))
	allocate(ipar(0:npop))
	allocate(x((ngen+1)*npop))
	allocate(y((ngen+1)*npop))
	allocate(ym((ngen+1)*npop))
	allocate(yf((ngen+1)*npop))
	allocate(gene(ngen+1,0:npop,nvars))
	allocate(same(npop))
	allocate(same_gen(npop))
	allocate(same_pop(npop))
	allocate(finished(npop+1))

	error=0d0
	mut_min=0.0005
	mut_max=0.1
	mutate=mut_min
	mut_dir=1
	imut=0
	dfit2=0d0

	same=.false.
	
	var(1,1,1:nvars)=var0(1:nvars)
	do j=2,npop
		do k=1,nvars
			var(1,j,k)=-1d0
			do while(var(1,j,k).lt.0d0.or.var(1,j,k).gt.1d0)
				var(1,j,k)=gasdev(idum)*dvar0(1,k)+var0(k)
			enddo
		enddo
	enddo

	do j=1,npop
		do k=1,nvars
			write(gene(1,j,k),'(f7.5)') var(1,j,k)
		enddo
	enddo
	
	max=0d0
	ncomp=0
	same_gen=0
	same_pop=npop+1
	do i=1,ngen
		print*,'========================================='
		print*,'Computing generation',i,ncomp
		print*,'========================================='
		if(i.eq.1) then
			var(i,0,1:nvars)=var(i,1,1:nvars)
			fit(i,0)=fit(i,1)
			fitobs(i,0,1:nobs)=fitobs(i,1,1:nobs)
			gene(i,0,1:nvars)(1:7)=gene(i,1,1:nvars)(1:7)
		else
			var(i,0,1:nvars)=var(i-1,0,1:nvars)
			fit(i,0)=fit(i-1,0)
			fitobs(i,0,1:nobs)=fitobs(i-1,0,1:nobs)
			gene(i,0,1:nvars)(1:7)=gene(i-1,0,1:nvars)(1:7)
		endif

		finished=.false.

		do j=1,npop
			if(same(j)) then
				print*,'individual ',j+(i-1)*npop,' is a clone of ',same_pop(j)+(same_gen(j)-1)*npop
				fit(i,j)=fit(same_gen(j),same_pop(j))
				fitobs(i,j,1:nobs)=fitobs(same_gen(j),same_pop(j),1:nobs)
			else
				fit(i,j)=geneticfun(j+(i-1)*npop,nvars,var(i,j,1:nvars),nobs,fitobs(i,j,1:nobs),.true.,error)
				ncomp=ncomp+1
				finished(j)=.true.
				if(fit(i,j).gt.max) then
					print*,'fittest individual: ',j+(i-1)*npop
					print*,' (generation ',i,' individual ',j,')'
					max=fit(i,j)
					var(i,0,1:nvars)=var(i,j,1:nvars)
					fit(i,0)=fit(i,j)
					fitobs(i,0,1:nobs)=fitobs(i,j,1:nobs)
					gene(i,0,1:nvars)(1:7)=gene(i,j,1:nvars)(1:7)
					call WriteStructure()
					call WriteOutput()
					call WriteRetrieval(j+(i-1)*npop,1d0/fit(i,j),var(i,j,1:nvars))
				endif
			endif
		enddo

		do k=1,nvars
			do l=1,nobs
				jj=0
				i1=i-10
				i2=i
				if(i1.lt.1) i1=1
				do ii=i1,i2
					do j=1,npop
						jj=jj+1
						x(jj)=var(ii,j,k)
						y(jj)=fitobs(ii,j,l)
					enddo
				enddo
				call ranking(y,ym,yf,jj,jj1,jj2)
				w(l,k)=dCor(x,yf,jj)
			enddo
		enddo

		do l=1,nobs
			mate(0:npop,l)=fitobs(i,0:npop,l)/sum(fitobs(i,0:npop,l))
			call ranking(mate(0:npop,l),male,female,npop+1,jj1,jj2)
		enddo
		
		mate(0:npop,0)=fit(i,0:npop)/sum(fit(i,0:npop))
		call ranking(mate(0:npop,0),male,female,npop+1,jj1,jj2)
		jj1=jj1-1
		jj2=jj2-1

		dfit=0d0
		do j=1,nvars
			dfit=dfit+(var(i,jj1,j)-var(i,jj2,j))**2
		enddo
		dfit=sqrt(dfit/real(nvars+1))

		if(i.gt.1) dfit=abs(fit(i,0)-fit(i-1,0))/(fit(i,0)+fit(i-1,0))

		if(mut_dir.lt.0) then
			mutate=10d0**(log10(mutate*mut_max)/2d0)
		else
			mutate=10d0**(log10(mutate*mut_min)/2d0)
		endif

		if(dfit.le.dfit2/4d0) imut=imut+1
		if(imut.gt.3) then
			mut_dir=-mut_dir
			imut=0
		endif
		dfit2=dfit

		print*,'========================================='
		print*,'Evolution sets in'
		print*,'Mutation rate: ',mutate
		print*,'========================================='
		ipar=0
		nclones=0
		do j=1,npop
			if(j.lt.npop/5) then
c random freakshow
			do k=1,nvars
				var(i+1,j,k)=ran1(idum)
				write(gene(i+1,j,k),'(f7.5)') var(i+1,j,k)
			enddo

			else
c normal evolution
1			continue

			r=ran1(idum)
			if(r.lt.0.70d0) then
				call select_mate(male,npop,jj1,idum,mate(0:npop,1:nobs),nobs,-1)
				call select_mate(male,npop,jj2,idum,mate(0:npop,1:nobs),nobs,jj1)
			else if(r.lt.0.80d0) then
				call select_mate(female,npop,jj1,idum,mate(0:npop,1:nobs),nobs,-1)
				call select_mate(female,npop,jj2,idum,mate(0:npop,1:nobs),nobs,jj1)
			else if(r.lt.0.90d0) then
				call select_mate(male,npop,jj1,idum,mate(0:npop,1:nobs),nobs,-1)
				call select_mate(female,npop,jj2,idum,mate(0:npop,1:nobs),nobs,jj1)
			else
				call select_mate(female,npop,jj1,idum,mate(0:npop,1:nobs),nobs,-1)
				call select_mate(male,npop,jj2,idum,mate(0:npop,1:nobs),nobs,jj1)
			endif

			if(jj1.eq.jj2) goto 1
			if(ipar(jj1).gt.(npop/4).or.ipar(jj2).gt.(npop/4)) goto 1
			if(jj1.eq.0.or.jj2.eq.0.and.ipar(0).gt.2) goto 1

			idupl=0
2			continue
			idupl=idupl+1
			do k=1,nvars
				r1(k)=0d0
				r2(k)=0d0
				do l=1,nobs
					r1(k)=r1(k)+mate(jj1,l)*w(l,k)
					r2(k)=r2(k)+mate(jj2,l)*w(l,k)
				enddo
			enddo
			r1(1:nvars)=r1(1:nvars)/sum(r1(1:nvars))
			r2(1:nvars)=r2(1:nvars)/sum(r2(1:nvars))

			creep=.false.
			icreep=int(ran1(idum)*5.0)+3
			if(ran1(idum).lt.0.5) then
				creep=.true.
				icreep=int(ran1(idum)*5.0)+1
			endif
			do k=1,nvars
				write(gene(i+1,j  ,k)(1:2),'("0.")')
			if(ran1(idum).gt.0.5d0) then	!gene_cross) then
c	use complex cross-over method with genes
				do ig=3,7
					r=ran1(idum)
					if(r.lt.0.5d0) then		!(r1(k)/(r1(k)+r2(k)))) then
						gene(i+1,j  ,k)(ig:ig)=gene(i,jj1,k)(ig:ig)
					else
						gene(i+1,j  ,k)(ig:ig)=gene(i,jj2,k)(ig:ig)
					endif
				enddo
				read(gene(i+1,j  ,k)(1:7),*,err=100) var(i+1,j  ,k)
100				continue
			else
c	use simple cross-over method
4				continue
c				xx=gasdev(idum)*(mutate+abs(var(i,jj1,k)-var(i,jj2,k)))+(r1(k)/(r1(k)+r2(k)))
				xx=gasdev(idum)+0.5d0
				var(i+1,j,k)=var(i,jj1,k)*xx+var(i,jj2,k)*(1d0-xx)
				if(var(i+1,j,k).lt.0d0.or.var(i+1,j,k).ge.1d0) goto 4
				write(gene(i+1,j,k),'(f7.5)') var(i+1,j,k)
			endif
				do l=1,5
					if(ran1(idum).lt.mutate) then
						if(.not.creep) then
c				one point mutation
							ig=int(ran1(idum)*5.0)+3
							ig=icreep
							write(gene(i+1,j,k)(ig:ig),'(i1)') int(ran1(idum)*10d0)
							read(gene(i+1,j,k)(1:7),*,err=1) var(i+1,j,k)
						else
c				creep mutation
3							continue
							ig=int(ran1(idum)*5.0)+1
							ig=icreep
							if(ran1(idum).lt.0.5) then
								var(i+1,j,k)=var(i+1,j,k)+10d0**(-real(ig))
							else
								var(i+1,j,k)=var(i+1,j,k)-10d0**(-real(ig))
							endif
							if(var(i+1,j,k).gt.1d0.or.var(i+1,j,k).lt.0d0) then
								read(gene(i+1,j,k)(1:7),*,err=1) var(i+1,j,k)
								goto 3
							endif
							write(gene(i+1,j,k),'(f7.5)') var(i+1,j,k)
						endif
					endif
				enddo
			enddo
			if(idupl.lt.npop*ngen) then
				do ii=1,i
					do jj=1,npop
						same(j)=.true.
						do k=1,nvars
							do ig=1,7
								if(gene(i+1,j,k)(ig:ig).ne.gene(ii,jj,k)(ig:ig)) same(j)=.false.
							enddo
						enddo
						if(same(j)) then
							if(nclones.ge.npop/3) goto 1
							nclones=nclones+1
							same_pop(j)=jj
							same_gen(j)=ii
							goto 50
						endif
					enddo
				enddo
				do jj=1,j-1
					same(j)=.true.
					do k=1,nvars
						do ig=1,7
							if(gene(i+1,j,k)(ig:ig).ne.gene(i+1,jj,k)(ig:ig)) same(j)=.false.
						enddo
					enddo
					if(same(j)) then
						if(nclones.ge.npop/4) goto 1
						nclones=nclones+1
						same_pop(j)=jj
						same_gen(j)=i+1
						goto 50
					endif
				enddo				
			endif
50			continue
			ipar(jj1)=ipar(jj1)+1
			ipar(jj2)=ipar(jj2)+1

			endif
		enddo
	enddo

	var0(1:nvars)=var(ngen,0,1:nvars)
	
	deallocate(var)
	deallocate(fit)
	deallocate(fitobs)
	deallocate(w)
	deallocate(mate)
	deallocate(x)
	deallocate(y)
	deallocate(gene)

	return
	end
	
	subroutine select_mate(mate,npop,jj,idum,mate_obs,nobs,jpartner)
	IMPLICIT NONE
	integer npop,jj,idum,j,jpartner,nobs,i
	real*8 mate(0:npop),r,ran1,mate_obs(0:npop,1:nobs),selection(0:npop),f

	selection(0:npop)=mate(0:npop)
	if(jpartner.gt.0.and.ran1(idum).gt.0.5) then
		do j=0,npop
			f=0d0
			do i=1,nobs
				f=f+(mate_obs(j,i)/sum(mate_obs(j,1:nobs)
     &				-mate_obs(jpartner,i)/sum(mate_obs(jpartner,1:nobs))))**2
			enddo
			selection(j)=selection(j)*(sqrt(f/real(nobs))+real(npop))
		enddo
	endif
	f=sum(selection(0:npop))
	selection(0:npop)=selection(0:npop)/f

	r=ran1(idum)
	do jj=0,npop
		r=r-selection(jj)
		if(r.lt.0d0) return
	enddo

	jj=0
	
	return
	end


	subroutine ranking(mate,male,female,npop,jbest,jmedian)
	IMPLICIT NONE
	integer npop,i,j,rank(npop),k,jbest,jmedian
	real*8 mate(npop),max,male(npop),female(npop),p

	female(1:npop)=mate(1:npop)

	rank=1
	do i=1,npop
		rank(i)=i
	enddo
	do i=1,npop
		max=-1d100
		k=i
		do j=1,npop
			if(mate(j).ge.max) then
				max=mate(j)
				k=j
			endif
		enddo
		rank(k)=i
		mate(k)=-1d200
	enddo
	
	do i=1,npop
		mate(i)=real(npop+1-rank(i))
		if(rank(i).eq.2) jbest=i
		if(rank(i)-npop/2.ne.0) jmedian=i
	enddo
	mate(1:npop)=mate(1:npop)/sum(mate(1:npop))

	male(1:npop)=mate(1:npop)

	female(1:npop)=female(1:npop)*male(1:npop)
	female(1:npop)=female(1:npop)/sum(female(1:npop))

	return
	end
	
	





	real*8 function dCor(x,y,n)
	IMPLICIT NONE
	integer n,i,j
        real*8,allocatable :: Ax(:,:), Ay(:,:)
	real*8 x(n),y(n)
	real*8 avx1(n),avx2(n),avX
	real*8 avy1(n),avy2(n),avY
	real*8 dCov2,dVarX2,dVarY2

	allocate(Ax(n,n))
	allocate(Ay(n,n))
	
	avX=0d0
	avY=0d0
	do i=1,n
		do j=1,n
			Ax(i,j)=abs(x(i)-x(j))
			avX=avX+Ax(i,j)
			Ay(i,j)=abs(y(i)-y(j))
			avY=avY+Ay(i,j)
		enddo
	enddo
	avX=avX/real(n*n)
	avY=avY/real(n*n)
	do i=1,n
		avx1(i)=0d0
		avx2(i)=0d0
		avy1(i)=0d0
		avy2(i)=0d0
		do j=1,n
			avx1(i)=avx1(i)+Ax(i,j)
			avx2(i)=avx2(i)+Ax(j,i)
			avy1(i)=avy1(i)+Ay(i,j)
			avy2(i)=avy2(i)+Ay(j,i)
		enddo
		avx1(i)=avx1(i)/real(n)
		avx2(i)=avx2(i)/real(n)
		avy1(i)=avy1(i)/real(n)
		avy2(i)=avy2(i)/real(n)
	enddo
	dCov2=0d0
	dVarX2=0d0
	dVarY2=0d0
	do i=1,n
		do j=1,n
			Ax(i,j)=Ax(i,j)-avx1(i)-avx2(j)+avX
			Ay(i,j)=Ay(i,j)-avy1(i)-avy2(j)+avY
			dCov2=dCov2+Ax(i,j)*Ay(i,j)
			dVarX2=dVarX2+Ax(i,j)*Ax(i,j)
			dVarY2=dVarY2+Ay(i,j)*Ay(i,j)
		enddo
	enddo

	dCor=sqrt(dCov2)/(sqrt(sqrt(dVarX2)*sqrt(dVarY2)))

	deallocate(Ax)
	deallocate(Ay)
	
	return
	end
	




