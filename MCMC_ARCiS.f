      subroutine MCMC(likelihood,x0,ny,NDIM,NBURN,N_UNIQUE,epsinit)
      implicit none

      integer NDIM, N_UNIQUE, NBURN, ny
c      parameter (NDIM=3)
c      parameter (N_UNIQUE=10000)
c      parameter (NBURN=300)

      double precision SD, EPS, epsinit
      parameter (EPS=1.0D-6)

      double precision x(NDIM), xnew(NDIM), mean(NDIM), z(NDIM), xbest(NDIM)
      double precision cov(NDIM,NDIM), chol(NDIM,NDIM), x0(NDIM)
      double precision logp, logp_new, alpha, u, likelihood, lambda, logp_best
      double precision target_acc,start_acc,end_acc,p_acc,curr_acc,eta_curr,eta_end
      double precision samples(NDIM, N_UNIQUE+NBURN), x_mapped(NDIM),eta_start
      double precision beta_curr,beta_start
      double precision samples_mapped(NDIM, N_UNIQUE+NBURN), cov0(NDIM,NDIM)
      integer weights(N_UNIQUE+NBURN)
      integer i, j, k, accept, step, nacc
      external likelihood
	integer*4 counts, count_rate, count_max
	real*8 starttime,stoptime,remaining
	logical burnstart
	
      SD=(2.38D0**2)/real(NDIM)
      call random_seed_f77()

      accept = 0
      step = 1
      nacc = 1
      weights = 0
      lambda = 1.0

	end_acc=0.234
	start_acc=0.233
	p_acc=2.0
	eta_start=(1d0-0.01**(10d0*start_acc/real(NBURN))) !0.01
	eta_end=(1d0-0.01**(2d0*end_acc/real(NBURN))) !0.001
	eta_curr=eta_start
	
	beta_start=0.1
	beta_curr=beta_start
	
	curr_acc=end_acc

      do j = 1, NDIM
         x(j) = x0(j)
         mean(j) = x(j)
         do k = 1, NDIM
            if (j .eq. k) then
               cov(j,k) = 1.0D0
            else
               cov(j,k) = 0.0D0
            endif
         enddo
      enddo
      cov=cov*epsinit**2

      logp = (likelihood(x, ny, x_mapped))
    	xbest=xnew
        logp_best=logp_new
      do j = 1, NDIM
         samples(j,1) = x(j)
         samples_mapped(j,1) = x_mapped(j)
      enddo
      weights(1) = 1
	call write_pew_output(samples_mapped(1:NDIM,1:NBURN),weights(1:NBURN),NDIM,1,0)
	call SYSTEM_CLOCK(counts, count_rate, count_max)
	starttime = DBLE(counts)/DBLE(count_rate)
	burnstart=.true.

10    continue
         step = step + 1

c         if (nacc > 4 .and. nacc <= NBURN) then
         if (nacc > NBURN/2 .and. nacc <= NBURN) then
            call compute_mean(samples(1:NDIM,NBURN/2:nacc), weights(NBURN/2:nacc), mean, NDIM, nacc+1-NBURN/2)
            call compute_cov(samples(1:NDIM,NBURN/2:nacc), weights(NBURN/2:nacc), mean, cov, NDIM, nacc+1-NBURN/2, EPS)
		else if(nacc <= 4) then
			cov=cov/2d0
         endif
		if(nacc.lt.NBURN/2) then
			target_acc=start_acc+(end_acc-start_acc)*(real(nacc)*2d0/real(NBURN/2)-1d0)**p_acc
		else
			target_acc=end_acc
		endif
         if (nacc <= NBURN) then
			lambda=lambda*exp(((curr_acc)-target_acc)/(real(step)**0.6))
         endif

         call cholesky(cov, chol, NDIM)
         call random_normal_vec(z, NDIM)
		if(burnstart.and.nacc.ge.NBURN) then
			x=xbest
			logp=logp_best
			write(*,*) '===================================='
			write(*,*) 'Burn-in fase done, starting sampling'
			write(*,*) '     starting from best model so far'
			write(*,*) '===================================='
			burnstart=.false.
			curr_acc=end_acc
		endif

         do j = 1, NDIM
            xnew(j) = x(j)
            do k = 1, NDIM
               xnew(j) = xnew(j) + lambda**2 * sqrt(SD) * chol(j,k) * z(k)
            enddo
			if(nacc.gt.NBURN) then
				if(xnew(j).gt.1d0.or.xnew(j).lt.0d0) then
					weights(nacc)=weights(nacc) + 1
					goto 10
				endif
			endif
         enddo
		if(nacc.le.NBURN) call fold_MCMC(xnew,NDIM)
         logp_new = (likelihood(xnew, ny, x_mapped))

		if(nacc.le.NBURN) then
			eta_curr=10d0**(log10(eta_start)+log10(eta_end/eta_start)*real(nacc)/real(NBURN))
		else
			eta_curr=eta_end
		endif
		if(nacc.le.NBURN/2) then
			beta_curr=10d0**(log10(beta_start)+log10(1d0/beta_start)*real(nacc)/real(NBURN/2))
		else
			beta_curr=1d0
		endif

         if(logp_new.gt.logp_best) then
         	xbest=xnew
         	logp_best=logp_new
         endif
         call random_number(u)
         alpha = min(1.0D0, dexp((logp_new - logp)*beta_curr))

		curr_acc=curr_acc*(1d0-eta_curr)
         if (u < alpha) then
			curr_acc=curr_acc+eta_curr
			if(nacc.ge.NBURN) then
		      call write_pew_output(samples_mapped(1:NDIM,nacc),weights(nacc),NDIM,1,1)
		    endif
            do j = 1, NDIM
               x(j) = xnew(j)
            enddo
            logp = logp_new
            accept = accept + 1
            nacc = nacc + 1
            if (nacc > N_UNIQUE+NBURN) goto 99
            do j = 1, NDIM
               samples(j,nacc) = x(j)
               samples_mapped(j,nacc) = x_mapped(j)
            enddo
            weights(nacc) = 1

			if(nacc.lt.NBURN) then
				write(*,'(a,f6.1,a)') "Burn-in at: ",100d0*real(nacc)/real(NBURN),"%"
			else
				write(*,'(a,f6.1,a)') "Sampling at: ",100d0*real(nacc-NBURN)/real(N_UNIQUE),"%"
			endif
			call SYSTEM_CLOCK(counts, count_rate, count_max)
			stoptime = DBLE(counts)/DBLE(count_rate)

			if(nacc.ge.NBURN) then
				remaining=real(N_UNIQUE-(nacc-NBURN))/end_acc
			else if(nacc.ge.NBURN/2) then
				remaining=real(N_UNIQUE)/end_acc
				remaining=remaining+real(NBURN-nacc)/end_acc
			else
				remaining=real(NBURN/2)*(-atan(-sqrt((end_acc-start_acc)/end_acc))
     &		+atan(sqrt((end_acc-start_acc)/end_acc)*real(2*nacc-NBURN/2)/real(NBURN/2)))/(2d0*sqrt(start_acc*(end_acc-start_acc)))
				remaining=(real(step)/remaining)*real(NBURN/2)*(atan(sqrt((end_acc-start_acc)/end_acc))
     &		-atan(sqrt((end_acc-start_acc)/end_acc)*real(2*nacc-NBURN/2)/real(NBURN/2)))/(2d0*sqrt(start_acc*(end_acc-start_acc)))
				remaining=remaining+real(N_UNIQUE)/end_acc
				remaining=remaining+real(NBURN/2)/end_acc
    		endif
			remaining=(stoptime-starttime)*remaining/real(step)
			if(remaining.gt.3600d0*24d0) then
				write(*,'(a,f6.1,a)') "time remaining: " ,remaining/3600d0/24d0, " days"
			else if(remaining.gt.3600d0) then
				write(*,'(a,f6.1,a)') "time remaining: " ,remaining/3600d0, " hours"
			else if(remaining.gt.60d0) then
				write(*,'(a,f6.1,a)') "time remaining: " ,remaining/60d0, " minutes"
			else
				write(*,'(a,f6.1,a)') "time remaining: " ,remaining, " seconds"
			endif
			write(*,'(a,f6.1,a)') "acceptance: ",100d0*curr_acc,"%"
         else
            weights(nacc) = weights(nacc) + 1
         endif

         goto 10

99    continue
      close(10)
      print *, 'Final unique samples:', N_UNIQUE
      print *, 'Total steps taken:  ', step
      print *, 'Acceptance rate:    ', dble(accept)/dble(step)
      call write_pew_output(samples_mapped(1:NDIM,NBURN+1:N_UNIQUE+NBURN),weights(NBURN+1:N_UNIQUE+NBURN),
     &		NDIM,N_UNIQUE,2)

c      call write_pew_output(samples_mapped(1:NDIM,NBURN+1:N_UNIQUE+NBURN),weights(NBURN+1:N_UNIQUE+NBURN),
c     &		NDIM,N_UNIQUE,0)
c      call write_pew_output(samples_mapped(1:NDIM,NBURN+1:N_UNIQUE+NBURN),weights(NBURN+1:N_UNIQUE+NBURN),
c     &		NDIM,N_UNIQUE,1)
c      call write_pew_output(samples_mapped(1:NDIM,NBURN+1:N_UNIQUE+NBURN),weights(NBURN+1:N_UNIQUE+NBURN),
c     &		NDIM,N_UNIQUE,2)
      
      end

	subroutine write_pew_output(samples,weights,ndim,nsamples,init)
	use GlobalSetup
	IMPLICIT NONE
	integer init,i,j,ndim,nsamples
	integer weights(nsamples)
	real*8 samples(ndim,nsamples)
	character*2000 line
	
	if(init.eq.0) then
		open(unit=83,file=trim(outputdir) // "posteriorMCMC.dat",FORM="FORMATTED",ACCESS="STREAM")
		line="# "
		do i=1,ndim
			if(i.eq.1) then
				write(line(3:12),'(a10)') RetPar(i)%keyword
			else
				write(line(i*12-11:i*12+1),'(" ",a11)') RetPar(i)%keyword
			endif
		enddo
		write(83,'(a)') trim(line)
		flush(83)
	else if(init.eq.1) then
		do i=1,nsamples
			do j=1,weights(i)
				write(line,'("(",i0.4,"es12.4)")') ndim
				write(83,line) samples(1:ndim,i)
			enddo
		enddo
		flush(83)
	else
		close(unit=83)
	endif
	return
	end

            
      subroutine compute_mean(samples, weights, mean, ndim, nsamp)
      implicit none
      integer ndim, nsamp, i, j, tot, weights(nsamp)
      double precision samples(ndim,nsamp), mean(ndim)
	tot=0
	do j=1,nsamp
		tot=tot+weights(j)
	enddo
      do i = 1, ndim
         mean(i) = 0.0D0
         do j = 1, nsamp
            mean(i) = mean(i) + samples(i,j)*dble(weights(j))
         enddo
         mean(i) = mean(i) / dble(tot)
      enddo
      return
      end

      subroutine compute_cov(samples, weights, mean, cov, ndim, nsamp, eps)
      implicit none
      integer ndim, nsamp, i, j, k, tot, weights(nsamp)
      double precision samples(ndim,nsamp), mean(ndim)
      double precision cov(ndim,ndim), eps, dx(ndim)
	tot=0
	do j=1,nsamp
		tot=tot+weights(j)
	enddo

      do i = 1, ndim
         do j = 1, ndim
            cov(i,j) = 0.0D0
         enddo
      enddo      
      do k = 1, nsamp
         do i = 1, ndim
            dx(i) = samples(i,k) - mean(i)
         enddo
         do i = 1, ndim
            do j = 1, ndim
               cov(i,j) = cov(i,j) + dx(i)*dx(j)*dble(weights(k))
            enddo
         enddo
      enddo
      do i = 1, ndim
         do j = 1, ndim
            cov(i,j) = cov(i,j) / dble(tot - 1)
            if (i .eq. j) cov(i,j) = max(1d-8,cov(i,j)*(1d0 + eps))
         enddo
      enddo
      return
      end

      subroutine cholesky_lapack(A, L, n)
      implicit none
      integer n, info
      double precision A(n,n), L(n,n)

C     Local copy of A because DPOTRF works in-place
      double precision Acopy(n,n)
      integer i, j

C     Copy A into Acopy
      do i = 1, n
         do j = 1, n
            Acopy(i,j) = A(i,j)
         enddo
      enddo

C     Call LAPACK Cholesky decomposition: Acopy -> lower triangle of L
      call DPOTRF('L', n, Acopy, n, info)

      if (info /= 0) then
         print *, 'DPOTRF failed, INFO =', info
         stop 'Cholesky decomposition error'
      endif

C     Extract lower triangle into L, set upper triangle to zero
      do i = 1, n
         do j = 1, n
            if (i >= j) then
               L(i,j) = Acopy(i,j)
            else
               L(i,j) = 0.0D0
            endif
         enddo
      enddo

      return
      end

      subroutine cholesky(A, L, n)
      implicit none
      integer n, i, j, k
      double precision A(n,n), L(n,n), sum
      do i = 1, n
         do j = 1, n
            L(i,j) = 0.0D0
         enddo
      enddo
      do i = 1, n
         do j = 1, i
            sum = A(i,j)
            do k = 1, j-1
               sum = sum - L(i,k)*L(j,k)
            enddo
            if (i .eq. j) then
               if (sum > 0.0D0) then
                  L(i,j) = sqrt(sum)
               else
                  L(i,j) = 0.0D0
               endif
            else
               if (L(j,j) .ne. 0.0D0) then
                  L(i,j) = sum / L(j,j)
               else
                  L(i,j) = 0.0D0
               endif
            endif
         enddo
      enddo
      return
      end

      subroutine random_normal_vec(z, ndim)
      implicit none
      integer ndim, i
      double precision z(ndim), u1, u2
      do i = 1, ndim, 2
         call random_number(u1)
         call random_number(u2)
         z(i) = sqrt(-2.0D0*log(u1)) * cos(6.283185307179586D0*u2)
         if (i+1 <= ndim) then
            z(i+1) = sqrt(-2.0D0*log(u1)) * sin(6.283185307179586D0*u2)
         endif
      enddo
      return
      end

      subroutine random_seed_f77()
      implicit none
      double precision dummy
      call random_number(dummy)
      return
      end
 
      subroutine identity_matrix(A, n)
      implicit none
      integer n, i, j
      double precision A(n,n)
      do i = 1, n
         do j = 1, n
            if (i .eq. j) then
               A(i,j) = 1.0D0
            else
               A(i,j) = 0.0D0
            endif
         enddo
      enddo
      return
      end
     
	subroutine fold_MCMC(var,nvar)
	IMPLICIT NONE
	integer nvar,i
	real*8 var(nvar)
	
	do i=1,nvar
1		continue
		if(var(i).gt.1d0) then
			var(i)=2d0-var(i)
			goto 1
		endif
		if(var(i).lt.0d0) then
			var(i)=-var(i)
			goto 1
		endif
	enddo

	return
	end







*> \brief \b DISNAN tests input for NaN.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DISNAN + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       LOGICAL FUNCTION DISNAN( DIN )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION, INTENT(IN) :: DIN
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*> future.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIN
*> \verbatim
*>          DIN is DOUBLE PRECISION
*>          Input to test for NaN.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup isnan
*
*  =====================================================================
      LOGICAL FUNCTION DISNAN( DIN )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN
*     ..
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
*  ..
*  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END
*> \brief \b DLAISNAN tests input for NaN by comparing two arguments for inequality.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAISNAN + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaisnan.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaisnan.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaisnan.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> This routine is not for general use.  It exists solely to avoid
*> over-optimization in DISNAN.
*>
*> DLAISNAN checks for NaNs by comparing its two arguments for
*> inequality.  NaN is the only floating-point value where NaN != NaN
*> returns .TRUE.  To check for NaNs, pass the same variable as both
*> arguments.
*>
*> A compiler must assume that the two arguments are
*> not the same variable, and the test will not be optimized away.
*> Interprocedural or whole-program optimization may delete this
*> test.  The ISNAN functions will be replaced by the correct
*> Fortran 03 intrinsic once the intrinsic is widely available.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIN1
*> \verbatim
*>          DIN1 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] DIN2
*> \verbatim
*>          DIN2 is DOUBLE PRECISION
*>          Two numbers to compare for inequality.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup laisnan
*
*  =====================================================================
      LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
*     ..
*
*  =====================================================================
*
*  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END
*> \brief \b DPOTRF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DPOTRF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpotrf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpotrf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpotrf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DPOTRF computes the Cholesky factorization of a real symmetric
*> positive definite matrix A.
*>
*> The factorization has the form
*>    A = U**T * U,  if UPLO = 'U', or
*>    A = L  * L**T,  if UPLO = 'L',
*> where U is an upper triangular matrix and L is lower triangular.
*>
*> This is the block version of the algorithm, calling Level 3 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          N-by-N upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky
*>          factorization A = U**T*U or A = L*L**T.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the leading principal minor of order i
*>                is not positive, and the factorization could not be
*>                completed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup potrf
*
*  =====================================================================
      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DPOTRF2, DSYRK, DTRSM,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         CALL DPOTRF2( UPLO, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code.
*
         IF( UPPER ) THEN
*
*           Compute the Cholesky factorization A = U**T*U.
*
            DO 10 J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Upper', 'Transpose', JB, J-1, -ONE,
     $                     A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTRF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block row.
*
                  CALL DGEMM( 'Transpose', 'No transpose', JB,
     $                        N-J-JB+1,
     $                        J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ),
     $                        LDA, ONE, A( J, J+JB ), LDA )
                  CALL DTRSM( 'Left', 'Upper', 'Transpose',
     $                        'Non-unit',
     $                        JB, N-J-JB+1, ONE, A( J, J ), LDA,
     $                        A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
*
         ELSE
*
*           Compute the Cholesky factorization A = L*L**T.
*
            DO 20 J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Lower', 'No transpose', JB, J-1, -ONE,
     $                     A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTRF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block column.
*
                  CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1,
     $                        JB,
     $                        J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ),
     $                        LDA, ONE, A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'Transpose',
     $                        'Non-unit',
     $                        N-J-JB+1, JB, ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = INFO + J - 1
*
   40 CONTINUE
      RETURN
*
*     End of DPOTRF
*
      END
*> \brief \b DPOTRF2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE DPOTRF2( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DPOTRF2 computes the Cholesky factorization of a real symmetric
*> positive definite matrix A using the recursive algorithm.
*>
*> The factorization has the form
*>    A = U**T * U,  if UPLO = 'U', or
*>    A = L  * L**T,  if UPLO = 'L',
*> where U is an upper triangular matrix and L is lower triangular.
*>
*> This is the recursive version of the algorithm. It divides
*> the matrix into four submatrices:
*>
*>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
*>    A = [ -----|----- ]  with n1 = n/2
*>        [  A21 | A22  ]       n2 = n-n1
*>
*> The subroutine calls itself to factor A11. Update and scale A21
*> or A12, update A22 then calls itself to factor A22.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          N-by-N upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky
*>          factorization A = U**T*U or A = L*L**T.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the leading principal minor of order i
*>                is not positive, and the factorization could not be
*>                completed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup potrf2
*
*  =====================================================================
      RECURSIVE SUBROUTINE DPOTRF2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            N1, N2, IINFO
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSYRK, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     N=1 case
*
      IF( N.EQ.1 ) THEN
*
*        Test for non-positive-definiteness
*
         IF( A( 1, 1 ).LE.ZERO.OR.DISNAN( A( 1, 1 ) ) ) THEN
            INFO = 1
            RETURN
         END IF
*
*        Factor
*
         A( 1, 1 ) = SQRT( A( 1, 1 ) )
*
*     Use recursive code
*
      ELSE
         N1 = N/2
         N2 = N-N1
*
*        Factor A11
*
         CALL DPOTRF2( UPLO, N1, A( 1, 1 ), LDA, IINFO )
         IF ( IINFO.NE.0 ) THEN
            INFO = IINFO
            RETURN
         END IF
*
*        Compute the Cholesky factorization A = U**T*U
*
         IF( UPPER ) THEN
*
*           Update and scale A12
*
            CALL DTRSM( 'L', 'U', 'T', 'N', N1, N2, ONE,
     $                  A( 1, 1 ), LDA, A( 1, N1+1 ), LDA )
*
*           Update and factor A22
*
            CALL DSYRK( UPLO, 'T', N2, N1, -ONE, A( 1, N1+1 ), LDA,
     $                  ONE, A( N1+1, N1+1 ), LDA )
            CALL DPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
*
*        Compute the Cholesky factorization A = L*L**T
*
         ELSE
*
*           Update and scale A21
*
            CALL DTRSM( 'R', 'L', 'T', 'N', N2, N1, ONE,
     $                  A( 1, 1 ), LDA, A( N1+1, 1 ), LDA )
*
*           Update and factor A22
*
            CALL DSYRK( UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA,
     $                  ONE, A( N1+1, N1+1 ), LDA )
            CALL DPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
         END IF
      END IF
      RETURN
*
*     End of DPOTRF2
*
      END
*> \brief \b DSYRK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDC,N
*       CHARACTER TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),C(LDC,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYRK  performs one of the symmetric rank k operations
*>
*>    C := alpha*A*A**T + beta*C,
*>
*> or
*>
*>    C := alpha*A**T*A + beta*C,
*>
*> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
*> in the second case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*>           triangular  part  of the  array  C  is to be  referenced  as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*>                                  is to be referenced.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry,  TRANS  specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
*>
*>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
*>
*>              TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N specifies the order of the matrix C.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*>           of  columns   of  the   matrix   A,   and  on   entry   with
*>           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*>           of rows of the matrix  A.  K must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by n  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*>           be at least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension ( LDC, N )
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*>           upper triangular part of the array C must contain the upper
*>           triangular part  of the  symmetric matrix  and the strictly
*>           lower triangular part of C is not referenced.  On exit, the
*>           upper triangular part of the array  C is overwritten by the
*>           upper triangular part of the updated matrix.
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*>           lower triangular part of the array C must contain the lower
*>           triangular part  of the  symmetric matrix  and the strictly
*>           upper triangular part of C is not referenced.  On exit, the
*>           lower triangular part of the array  C is overwritten by the
*>           lower triangular part of the updated matrix.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, n ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup herk
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Test the input parameters.
*
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T')) .AND.
     +         (.NOT.LSAME(TRANS,'C'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYRK ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  C := alpha*A*A**T + beta*C.
*
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*A**T*A + beta*C.
*
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP = ZERO
                      DO 190 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP = ZERO
                      DO 220 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYRK
*
      END
