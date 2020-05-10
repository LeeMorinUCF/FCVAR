      subroutine fpval(np,iq,stat,probs,bedf,ginv,pval)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
c This routine calculates P values.
c
c stat is test statistic.
c np is number of points for local approximation (probably 9).
c bedf contains quantiles of numerical distribution for specified.
c value of b or values of b and d.
c ginv contains quantiles of approximating chi-squared distribution.
c
      real*8 probs(221), bedf(221), ginv(221), gamma(3)
      real*8 xmat(25,3), yvect(25), xpy(3), xpxi(3,3), resid(25)
      nomax = 25
      nvmax = 3
      ndf = iq**2
c
c deal with extreme cases
c
      btiny = 0.5d0*bedf(1)
      bbig = 2.d0*bedf(221)
      if (stat.lt.btiny) then
        pval = 1.d0
        return
      end if
      if (stat.gt.bbig) then
        pval = 0.d0
        return
      end if
c
c find critical value closest to test statistic
c
      diffm = 1000.d0
      imin = 0
      do i=1,221
        diff = abs(stat - bedf(i))
        if (diff.lt.diffm) then
          diffm = diff
          imin = i
        end if
      end do
      nph = np/2
      nptop = 221 - nph
      if (imin.gt.nph.and.imin.lt.nptop) then
c
c imin is not too close to the end. Use np points around stat.
c
        do i=1,np
          ic = imin - nph - 1 + i
          yvect(i) = ginv(ic)
          xmat(i,1) = 1.d0
          xmat(i,2) = bedf(ic)
          xmat(i,3) = xmat(i,2)*bedf(ic)
        end do
        nvar = 3
        i1 = 0
        call olsqc(np,nomax,nvar,nvmax,i1,ssr,gamma,xpy,xpxi,xmat,
     &    yvect,resid)
        crfit = gamma(1) + gamma(2)*stat + gamma(3)*stat**2
        if (crfit.lt.1.d-6) crfit = 1.d-6
        call chicdf(crfit,pval,ndf)
        pval = 1.d0 - pval
      else
c
c imin is close to one of the ends. Use points from imin +/- nph to end.
c
        if (imin.le.nph) then
          np1 = imin + nph
          if (np1.lt.5) np1 = 5
          do i=1,np1
            yvect(i) = ginv(i)
            xmat(i,1) = 1.d0
            xmat(i,2) = bedf(i)
            xmat(i,3) = xmat(i,2)*bedf(i)
          end do
        else
          np1 = 222 - imin + nph
          if (np1.lt.5) np1 = 5
          do i=1,np1
            ic = 222 - i
            yvect(i) = ginv(ic)
            xmat(i,1) = 1.d0
            xmat(i,2) = bedf(ic)
            xmat(i,3) = xmat(i,2)*bedf(ic)
          end do
        end if
        nvar = 3
        i1 = 0
        call olsqc(np1,nomax,nvar,nvmax,i1,ssr,gamma,xpy,xpxi,xmat,
     &    yvect,resid)
        crfit = gamma(1) + gamma(2)*stat + gamma(3)*stat**2
        if (crfit.lt.1.d-6) crfit = 1.d-6
        call chicdf(crfit,pval,ndf)
        pval = 1.d0 - pval
      end if
      return
      end
      subroutine chicdf(arg,value,ndf)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
      arg1 = ndf/2.d0
      arg2 = arg/2.d0
      value = gammp(arg1,arg2)
      return
      end
      function gammp(a,x)
      implicit real*8 (a-h,o-z)
      if(x.lt.0.d0.or.a.le.0.d0) then
        write(6,*) 'Trouble in gammp.'
        stop
      end if
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammp = gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp = 1.d0 - gammcf
      endif
      return
      end
      subroutine gser(gamser,a,x,gln)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
      parameter (itmax=200,eps=1.d-7)
      gln = gammln(a)
      if (x.le.0.d0) then
        if (x.lt.0.d0) then
          write(6,*) 'Trouble in gser.'
          stop
        end if
        gamser = 0.d0
        return
      endif
      ap = a
      sum = 1.d0/a
      del = sum
      do n=1,itmax
        ap = ap + 1.d0
        del = del*x/ap
        sum = sum + del
        if(abs(del).lt.abs(sum)*eps) go to 1
      end do
      write(6,*) 'Warning! a too large or itmax too small in gser.'
 1    gamser = sum*exp(-x + a*log(x) - gln)
      return
      end
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
      function gammln(xx)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
      real*8 cof(6),stp,half,one,fpf,x,tmp,ser
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     &    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x = xx - one
      tmp = x + fpf
      tmp = (x + half)*log(tmp) - tmp
      ser = one
      do j=1,6
        x = x + one
        ser = ser + cof(j)/x
      end do
      gammln = tmp + log(stp*ser)
      return
      end
      subroutine gcf(gammcf,a,x,gln)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
      parameter (itmax=200,eps=1.d-7)
      gln = gammln(a)
      gold = 0.d0
      a0 = 1.d0
      a1 = x
      b0 = 0.d0
      b1 = 1.d0
      fac = 1.d0
      do n=1,itmax
        an = dble(n)
        ana = an - a
        a0 = (a1 + a0*ana)*fac
        b0 = (b1 + b0*ana)*fac
        anf = an*fac
        a1 = x*a0 + anf*a1
        b1 = x*b0 + anf*b1
        if(a1.ne.0.d0)then
          fac = 1.d0/a1
          g = b1*fac
          if(abs((g-gold)/g).lt.eps)go to 1
          gold = g
        endif
      end do
      write(6,*) 'Warning! a too large or itmax too small in gcf.'
   1  gammcf = dexp(-x + a*dlog(x) - gln)*g
      return
      end
