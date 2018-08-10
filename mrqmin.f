      SUBROUTINE mrqmin(y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *funcs,alamda,beta)
      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL*8 alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),
     *sig(ndata),y(ndata)
      PARAMETER (MMAX=100)
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,mfit
      REAL*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,da,mfit!,beta
      if(alamda.lt.0.)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001
		ochisq=1d200
        call mrqcof(y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,ochisq,funcs)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      do 14 j=1,mfit
        do 13 k=1,mfit
          covar(j,k)=alpha(j,k)
13      continue
        covar(j,j)=alpha(j,j)*(1.+alamda)
        da(j)=beta(j)
14    continue
      call gaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.)then
        call covsrt(covar,nca,ma,ia,mfit)
        call covsrt(alpha,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,ochisq,funcs)
      if(chisq.lt.ochisq)then
        alamda=0.1*alamda
        ochisq=chisq
        do 17 j=1,mfit
          do 16 k=1,mfit
            alpha(j,k)=covar(j,k)
16        continue
          beta(j)=da(j)
17      continue
        do 18 l=1,ma
          a(l)=atry(l)
18      continue
      else
        alamda=10.*alamda
        chisq=ochisq
      endif
      return
      END

      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
      INTEGER ma,mfit,npc,ia(ma)
      REAL*8 covar(npc,npc)
      INTEGER i,j,k
      REAL*8 swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.
          covar(j,i)=0.
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=100)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL*8 big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END

      SUBROUTINE mrqcof(y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,ochisq,
     *funcs)
      INTEGER ma,nalp,ndata,ia(ma)
      REAL*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),
     *y(ndata),ochisq
      EXTERNAL funcs
      INTEGER mfit,i,j,k,l,m,idum
      REAL*8 dy,sig2i,wt,ymod(ndata),dyda(ndata,ma),ran2
      mfit=0
      do 11 j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
11    continue
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.
12      continue
        beta(j)=0.
13    continue
      chisq=0.
      i=1
		call funcs(a,ymod,dyda,ma,ndata,i)
	do i=1,ndata
        sig2i=1./(sig(i)*sig(i))
        dy=y(i)-ymod(i)
        chisq=chisq+dy*dy*sig2i
	enddo
	if(chisq.gt.ochisq.and.ran2(idum).gt.0.1d0) return
		i=2
		call funcs(a,ymod,dyda,ma,ndata,i)
		chisq=0.
      do 16 i=1,ndata
        sig2i=1./(sig(i)*sig(i))
        dy=y(i)-ymod(i)
        j=0
        do 15 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            wt=dyda(i,l)*sig2i
            k=0
            do 14 m=1,l
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(i,m)
              endif
14          continue
            beta(j)=beta(j)+dy*wt
          endif
15      continue
        chisq=chisq+dy*dy*sig2i
16    continue
      do 18 j=2,mfit
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
17      continue
18    continue
      return
      END
