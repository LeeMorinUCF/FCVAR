      subroutine olsqc(nobs,nomax,nvar,nvmax,i1,ssr,beta,
     & xpy,xpxi,xmat,yvect,resid)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
c nomax is dimension of things that logically have length nobs (so that
c   nomax => nobs); nvmax is dimension of things that logically have
c length nvar (so that nvmax => nvar).
c This version returns beta, ssr, resid.
c It also returns xpxi, and can be told to use an xpxi given to it.
c It does not compute the covariance matrix; multiply xpxi by the
c square of sigest to do that
c
      real*8 xmat(nomax,nvmax), yvect(nomax), beta(nvmax), resid(nomax)
      real*8 xpy(nvmax), xpxi(nvmax,nvmax)
c
c  form and invert xpx matrix unless already done
c
      if(i1.ne.0) go to 250
c only one triangle of xpxi computed.
      do i=1,nvar
      do j=i,nvar
        sum = 0.d0
        do k=1,nobs
          sum = sum + xmat(k,i)*xmat(k,j)
        end do
        xpxi(j,i) = sum
        xpxi(i,j) = sum
      end do
      end do
c
c  now invert xpx
c
      call chol2(xpxi,nvmax,nvar,kf2)
  250 continue
c
c preliminary work has now been done.
c
c form xpy
c
      do i=1,nvar
        sum = 0.d0
        do k=1,nobs
          sum = sum + xmat(k,i)*yvect(k)
        end do
        xpy(i) = sum
      end do
c
c  now form estimates of beta.
c
      do i=1,nvar
        beta(i) = 0.d0
        do j=1,nvar
          beta(i) = beta(i) + xpxi(i,j)*xpy(j)
        end do
      end do
c
      ssr = 0.d0
      do k=1,nobs
        sum = yvect(k)
        do i=1,nvar
          sum = sum - xmat(k,i)*beta(i)
        end do
        ssr = ssr + sum**2
        resid(k) = sum
      end do
      return
      end
      subroutine chol2(a,m,n,kf2)
c
c Copyright (c) James G. MacKinnon 2011
c
c This routine uses the cholesky decomposition to invert a real
c symmetric matrix.
c
      implicit real*8 (a-h,o-z)
      real*8 a(m,m)
      kf2 = 0
      do 8 i=1,n
        kl = i - 1
        do 7 j=i,n
          if (i.gt.1) then
            do k=1,kl
              a(i,j) = a(i,j) - a(k,i)*a(k,j)
            end do
          else
            if (a(i,i).le.0.d0) then
              kf2 = i
              go to 20
            end if
          end if
          if (i.eq.j) then
            a(i,i) = dsqrt(a(i,i))
          else
            if (j.eq.i+1) ooa = 1.d0/a(i,i)
            a(i,j) = a(i,j)*ooa
          end if
 7      continue
 8    continue
      do 13 i=1,n
        do j=i,n
          ooa = 1.d0/a(j,j)
          if (i.ge.j) then
            t = 1.d0
            go to 12
          end if
          kl = j - 1
          t = 0.d0
          do k=i,kl
            t = t - a(i,k)*a(k,j)
          end do
 12       a(i,j) = t*ooa
        end do
 13   continue
      do 16 i=1,n
        do 15 j=i,n
          t = 0.d0
          do k=j,n
            t = t + a(i,k)*a(j,k)
          end do
          a(i,j) = t
          a(j,i) = t
 15     continue
 16   continue
 20   return
      end
