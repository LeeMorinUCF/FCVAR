
      implicit real*8 (a-h,o-z)
c
c read_test_fracdist.f reads test cases for fracdist.f, which is a
c program to compute P values or critical values for tests of
c fractional unit roots and cointegration.
c Version 1 (May, 14 2020) that only tests file IO.
c Copyright (c) Lealand Morin, 2020
c
c Need to put declaration statements before other commands.
      character*16 pval_header
      character*18 cval_header
c
c Set parameters for number of test cases.
      npval = 2400
      ncval = 720
c
c Open file to read test cases for p-values.
c
      open(unit=12,file='test_fpval.txt',status='OLD',err=333)
      go to 334
 333  write(6,*) 'File for p-value test cases not opened.'
      stop
 334  continue
c
c Open file to read test cases for p-values.
c
      open(unit=14,file='test_fcval.txt',status='OLD',err=335)
      go to 336
 335  write(6,*) 'File for critical value test cases not opened.'
      stop
 336  continue
c
c Open files to save results.
c
      open (unit=13,file='test_fpval.out')
      open (unit=15,file='test_fcval.out')
c
c Extract parameters for test case in each line.
c Evaluate p-values.
c Write result to output file.
c
      write(6,100) 'Evaluating test cases for p-values'
 100  format(/,A34,/)
      read(12,200) pval_header
 200  format(a16)
      write(6,*) 'File of test cases contains the following variables:'
      write(6,*) pval_header
      write(13,101) pval_header, 'pval'
 101  format(a16,1x,a4)
      isave = 0
      ipc = 1
      do ip=1,npval
        read(12,201) iscon, iq, bb, stat
 201    format(i1,1x,i2,1x,f5.3,1x,f8.4)
        call frval(iq,isave,ipc,iscon,bb,stat,pval,clevel,ccrit)
        write(13,102) iscon, iq, bb, stat, pval
 102    format(i1,1x,i2,1x,f5.3,1x,f8.4,1x,f6.4)
      end do
      write(6,103) 'See output in file test_fpval.out'
 103  format(/,A33,/)
c
c Extract parameters for test case in each line.
c Evaluate critical values.
c Write result to output file.
c
      write(6,104) 'Evaluating test cases for critical values'
  104 format(/,A41,/)
      read(14,202) cval_header
  202 format(a18)
      write(6,*) 'File of test cases contains the following variables:'
      write(6,*) cval_header
      write(15,105) cval_header, 'cval'
  105 format(a18,1x,a4)
      isave = 0
      ipc = 2
      do ip=1,ncval
        read(14,203) iscon, iq, bb, clevel
  203   format(i1,1x,i2,1x,f5.3,1x,f4.2)
c        ccrit = 123.12345678
        call frval(iq,isave,ipc,iscon,bb,stat,pval,clevel,ccrit)
        write(15,106) iscon, iq, bb, clevel, ccrit
  106   format(i1,1x,i2,1x,f5.3,1x,f4.2,1x,f8.4)
      end do
      write(6,107) 'See output in file test_fcval.out'
  107 format(/,A33,/)
      stop
      end
c
c Subroutines from fracdist.f:
c
      subroutine frval(iq,isave,ipc,iscon,bb,stat,pval,clevel,ccrit)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
c This is the main routine that is called for both P values and critical
c values.
c
c isave = 0 or 1. If isave=1, this call to frvalnc has the same values of
c iq and bb as the last one, so data files do not need to be read again.
c iscon = 0 means no constant term; any other value means there is one.
c
      real*8 xndf(221,31), bval(31), probs(221), ginv(221)
      real*8 bedf(221), estcrit(31)
      common/saveallnc/ probs, bval, xndf, bedf, estcrit, ginv
      nb = 31
      np = 221
      if (isave.eq.0) then
        call readdata(iq,iscon,probs,bval,xndf)
        do i=1,np
          do j=1,nb
            estcrit(j) = xndf(i,j)
          end do
          call blocal(nb,bb,estcrit,bval,bedf(i))
        end do
        call gcinv(iq,np,probs,ginv)
      end if
      np = 9
      if (ipc.eq.1) then
        call fpval(np,iq,stat,probs,bedf,ginv,pval)
      else
        call fpcrit(np,iq,clevel,probs,bedf,ginv,ccrit)
      end if
      return
      end
      subroutine gcinv(iq,np,probs,ginv)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2010
c
c This routine calculates inverse of chi-squared distribution for each
c probability
c
      real*8 probs(221), ginv(221)
      ndf = iq**2
      do i=1,np
        call ccdfinv(probs(i),ginv(i),ndf)
      end do
      return
      end
      subroutine readdata(iq,iscon,probs,bbb,xndf)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
c This routine reads data from whichever input file is appropriate
c for the no-constant case.
c If iscon=0, files are frmapp01.txt through frmapp12.txt
c If iscon.ne.0, files are frcapp01.txt through frcapp12.txt
c
      real*8 xndf(221,31), bbb(31), probs(221)
      character*6 dfirst
      character*4 dlast
      character*2 dq
      character*12 dshort
cunix character*31 dname
cunix character*19 dirname
cunix dirname = '/usr/local/urcdist/'
cwin  character*24 dname
cwin  character*12 dirname
cwin  dirname = 'C:\fracdist\'
      character*24 dname
      character*12 dirname
      dirname = 'C:\fracdist\'
c
      if (iscon.eq.0) then
        dfirst = 'frmapp'
      else
        dfirst = 'frcapp'
      end if
      dlast = '.txt'
      if (iq.ge.10) then
        write(dq,101) iq
 101    format(i2)
      else
        write(dq,102) iq
 102    format('0',i1)
      end if
c
c First try to read input file from current directory.
c
      dshort = dfirst//dq//dlast
      open(unit=2,file=dshort,status='OLD',err=333)
      go to 334
 333  continue
c
c Input file not in current directory.
c
      dname = dirname//dfirst//dq//dlast
      open(unit=2,file=dname)
 334  continue
      do ib=1,31
        read(2,200) bbb(ib)
 200    format(4x,f5.2)
        do j=1,221
          read(2,201) probs(j), xndf(j,ib)
 201      format(f6.4,3x,f16.12)
        end do
      end do
      return
      end
      subroutine blocal(nb,bb,estcrit,bval,bfit)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
c This routine calculates a fitted value for any specified b (bb)
c between 0.51 and 2.0, by fitting a locally quadratic approximation
c using nearby values for which quantiles have been estimated.
c
      real*8 estcrit(31), bval(31), weight(31)
      real*8 yvect(9), xmat(9,3), resid(9), xpy(3), xpxi(3,3), beta(3)
      nomax = 9
      nvmax = 3
      if (bb.lt.0.51d0.or.bb.gt.2.0d0) then
        write(6,*) 'Specified value of b out of range 0.51 to 2.0'
        write(6,*) 'Value of b is ', bb
        return
      end if
      jbot = 0
c
c jbot will be index of lowest b that gets positive weight
c
      nobs = 0
      nvar = 3
      do j=1,nb
        weight(j) = 1.d0 - 5.d0*abs(bval(j) - bb)
        if (weight(j).lt.1.d-12) then
          weight(j) = 0.d0
        else
          if (jbot.eq.0) jbot = j
          nobs = nobs + 1
        end if
      end do
      do i=1,nobs
        yvect(i) = weight(jbot + i - 1)*estcrit(jbot+i-1)
        xmat(i,1) = weight(jbot + i - 1)
        xmat(i,2) = bval(jbot + i - 1)*xmat(i,1)
        xmat(i,3) = bval(jbot + i - 1)*xmat(i,2)
      end do
      i1 = 0
      call olsqc(nobs,nomax,nvar,nvmax,i1,ssr,beta,xpy,xpxi,xmat,
     &  yvect,resid)
      bfit = beta(1) + beta(2)*bb + beta(3)*bb**2
      return
      end
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
      subroutine fpcrit(np,iq,clevel,probs,bedf,ginv,ccrit)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
c This routine calculates critical values.
c
c clevel is level for test.
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
c Handle extreme cases.
c
      ptiny = 0.0001d0
      pbig = 0.9999d0
      if (clevel.lt.ptiny) then
        ccrit = bedf(221)
        return
      end if
      if (clevel.gt.pbig) then
        ccrit = bedf(1)
        return
      end if
c
c find probability closest to test level
c
      diffm = 2.d0
      imin = 0
      cquant = 1.d0 - clevel
      call ccdfinv(cquant,gcq,ndf)
      do i=1,221
        diff = abs(cquant - probs(i))
        if (diff.lt.diffm) then
          diffm = diff
          imin = i
        end if
      end do
c
      nph = np/2
      nptop = 221 - nph
      if (imin.gt.nph.and.imin.lt.nptop) then
c
c imin is not too close to the end. Use np points around cquant.
c
        do i=1,np
          ic = imin - nph - 1 + i
          yvect(i) = bedf(ic)
          xmat(i,1) = 1.d0
          xmat(i,2) = ginv(ic)
          xmat(i,3) = xmat(i,2)*ginv(ic)
        end do
        nvar = 3
        i1 = 0
        call olsqc(np,nomax,nvar,nvmax,i1,ssr,gamma,xpy,xpxi,xmat,
     &    yvect,resid)
        ccrit = gamma(1) + gamma(2)*gcq + gamma(3)*gcq**2
      else
c
c imin is close to one of the ends. Use points from imin +/- nph to end.
c
        if (imin.le.nph) then
          np1 = imin + nph
          if (np1.lt.5) np1 = 5
          do i=1,np1
            yvect(i) = bedf(i)
            xmat(i,1) = 1.d0
            xmat(i,2) = ginv(i)
            xmat(i,3) = xmat(i,2)*ginv(i)
          end do
        else
          np1 = 222 - imin + nph
          if (np1.lt.5) np1 = 5
          do i=1,np1
            ic = 222 - i
            yvect(i) = bedf(ic)
            xmat(i,1) = 1.d0
            xmat(i,2) = ginv(ic)
            xmat(i,3) = xmat(i,2)*ginv(ic)
          end do
        end if
        nvar = 3
        i1 = 0
        call olsqc(np1,nomax,nvar,nvmax,i1,ssr,gamma,xpy,xpxi,xmat,
     &    yvect,resid)
        ccrit = gamma(1) + gamma(2)*gcq + gamma(3)*gcq**2
      end if
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
      subroutine ccdfinv(prob,quantile,ndf)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
c This routine inverts a chi-squared distribution.
c It is intended to be more accurate than chisqd for ndf > 2
c It uses chisqd for an initial estimate.
c
      if (ndf.le.2) then
        call chisqd(ndf,prob,quantile)
        return
      end if
      toler = 1.d-10
      call chisqd(ndf,prob,q1)
      call chicdf(q1,p1,ndf)
      if (p1.gt.prob) then
        qlow = .998d0*q1
        qhi = q1
 100    continue
        call chicdf(qlow,plow,ndf)
        qnew = .5d0*(qlow + qhi)
        call chicdf(qnew,pnew,ndf)
        if (pnew.gt.prob) then
          qhi = qnew
        else
          qlow = qnew
        end if
        diff = (qhi - qlow)/q1
        if (diff.lt.toler) then
          quantile = qnew
          return
        end if
        go to 100
      else
        qhi = 1.002d0*q1
        qlow = q1
 101  continue
      call chicdf(qlow,plow,ndf)
      qnew = .5d0*(qlow + qhi)
      call chicdf(qnew,pnew,ndf)
      if (pnew.gt.prob) then
        qhi = qnew
      else
        qlow = qnew
      end if
      diff = (qhi - qlow)/q1
      if (diff.lt.toler) then
        quantile = qnew
        return
      end if
      go to 101
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
      subroutine chisqd(n,p,q)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon 2011
c
c This provides a not very accurate inverse of the chi-squared distribution.
c Based on TOMS451.
c
      dimension c(21), a(19)
      data c(1)/1.565326d-3/, c(2)/1.060438d-3/,
     &  c(3)/-6.950356d-3/, c(4)/-1.323293d-2/,
     &  c(5)/2.277679d-2/, c(6)/-8.986007d-3/,
     &  c(7)/-1.513904d-2/, c(8)/2.530010d-3/,
     &  c(9)/-1.450117d-3/, c(10)/5.169654d-3/,
     &  c(11)/-1.153761d-2/, c(12)/1.128186d-2/,
     &  c(13)/2.607083d-2/, c(14)/-0.2237368/,
     &  c(15)/9.780499d-5/, c(16)/-8.426812d-4/,
     &  c(17)/3.125580d-3/, c(18)/-8.553069d-3/,
     &  c(19)/1.348028d-4/, c(20)/0.4713941d0/, c(21)/1.0000886d0/
      data a(1)/1.264616d-2/, a(2)/-1.425296d-2/,
     &  a(3)/1.400483d-2/, a(4)/-5.886090d-3/,
     &  a(5)/-1.091214d-2/, a(6)/-2.304527d-2/,
     &  a(7)/3.135411d-3/, a(8)/-2.728484d-4/,
     &  a(9)/-9.699681d-3/, a(10)/1.316872d-2/,
     &  a(11)/2.618914d-2/, a(12)/-0.2222222/,
     &  a(13)/5.406674d-5/, a(14)/3.483789d-5/,
     &  a(15)/-7.274761d-4/, a(16)/3.292181d-3/,
     &  a(17)/-8.729713d-3/, a(18)/0.4714045d0/, a(19)/1.d0/
      if (n.eq.1) then
        prob = 0.5d0*(1.d0 - p)
        call innorm(prob,q)
        q = q*q
        return
      else if (n.eq.2) then
        q = -2.0d0*dlog(1.d0 - p)
        return
      end if
      f = n
      f1 = 1.d0 / f
      call innorm(p,t)
      f2 = sqrt(f1) * t
      if (n.lt.(2+int(4.d0*abs(t)))) then
        q = (((((((c(1)*f2+c(2))*f2+c(3))*f2+c(4))*f2
     &    +c(5))*f2+c(6))*f2+c(7))*f1+((((((c(8)+c(9)*f2)*f2
     &    +c(10))*f2+c(11))*f2+c(12))*f2+c(13))*f2+c(14)))*f1 +
     &    (((((c(15)*f2+c(16))*f2+c(17))*f2+c(18))*f2
     &    +c(19))*f2+c(20))*f2+c(21)
      else
40      q = (((a(1)+a(2)*f2)*f1+(((a(3)+a(4)*f2)*f2
     &  +a(5))*f2+a(6)))*f1+(((((a(7)+a(8)*f2)*f2+a(9))*f2
     &  +a(10))*f2+a(11))*f2+a(12)))*f1 + (((((a(13)*f2
     &  +a(14))*f2+a(15))*f2+a(16))*f2+a(17))*f2*f2
     &  +a(18))*f2+a(19)
      end if
50    q = q*q*q*f
      return
      end
      subroutine innorm(prob,anorm)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) 1994 James G. MacKinnon
c
c Excellent inverse normal routine that adjusts crude result twice
c It seems to be accurate to about 14 digits
c It calls ddnor, which is included in this file.
c
c Simple formula is taken from Abramowitz & Stegun (1968)
c It should have abs. error < 4.5 * 10^-4
c 1994-7-11
c
      data c0/2.515517d0/, d1/1.432788d0/, c1/0.802853d0/
      data c2/0.010328d0/, d3/0.001308d0/, d2/0.189269d0/
      data const/.398942280401432678d0/
      if (prob.lt.0.d0.or.prob.gt.1.d0) then
         write(6,*) 'Attempt to find inverse normal of ', prob
         stop
      end if
      pr = prob
      if (prob.gt.0.5d0) pr = 1.d0 - prob
      arg = 1.d0/pr**2
      t = sqrt(log(arg))
      anorm = t - (c0 + c1*t + c2*t**2)/
     &  (1.d0 + d1*t + d2*t**2 + d3*t**3)
c
c Now correct crude result by direct method
c
      call ddnor(anorm,prob2)
      pr2 = 1.d0 - prob2
      arg = 1/pr2**2
      t = sqrt(log(arg))
      anorm2 = t - (c0 + c1*t + c2*t**2)/
     & (1.d0 + d1*t + d2*t**2 + d3*t**3)
      anorm = anorm + anorm - anorm2
      if (prob.lt.0.5d0) anorm = -anorm
c
c Now correct better result, using Taylor series approximation
c
      call ddnor(anorm,prob2)
      error = prob2 - prob
      dens = const*dexp(-.5d0*anorm**2)
      anorm = anorm - error/dens
      return
      end
      subroutine ddnor(ystar,gauss)
c
c Copyright (c) 1993 James G. MacKinnon
c
c This subroutine uses Cody's method to evaluate the cumulative
c normal distribution. it is probably accurate to 19 or 20
c significant digits. It was written by James MacKinnon late in
c 1977, based on the Cody article referred to in the documentation
c for IMSL subroutine mdnor.
c
c Modified 1993 to avoid changing the argument
c
      implicit real*8(a-h,o-z)
      real*8 p(6), q(5), a(9), b(8), c(5), d(4)
      data p(1)/-6.58749161529837803157d-04/,
     1     p(2)/-1.60837851487422766278d-02/,
     2     p(3)/-1.25781726111229246204d-01/,
     3     p(4)/-3.60344899949804439429d-01/,
     4     p(5)/-3.05326634961232344035d-01/,
     5     p(6)/-1.63153871373020978498d-02/
      data q(1)/2.33520497626869185443d-03/,
     1     q(2)/6.05183413124413191178d-02/,
     2     q(3)/5.27905102951428412248d-01/,
     3     q(4)/1.87295284992346047209d00/,
     4     q(5)/2.56852019228982242072d00/
      data a(1)/1.23033935479799725272d03/,
     1     a(2)/2.05107837782607146532d03/,
     2     a(3)/1.71204761263407058314d03/,
     3     a(4)/8.81952221241769090411d02/,
     4     a(5)/2.98635138197400131132d02/,
     5     a(6)/6.61191906371416294775d01/,
     6     a(7)/8.88314979438837594118d00/,
     7     a(8)/5.64188496988670089180d-01/,
     8     a(9)/2.15311535474403846343d-08/
      data b(1)/1.23033935480374942043d03/,
     1     b(2)/3.43936767414372163696d03/,
     2     b(3)/4.36261909014324715820d03/,
     3     b(4)/3.29079923573345962678d03/,
     4     b(5)/1.62138957456669018874d03/,
     5     b(6)/5.37181101862009857509d02/,
     6     b(7)/1.17693950891312499305d02/,
     7     b(8)/1.57449261107098347253d01/
      data c(1)/3.209377589138469472562d03/,
     1     c(2)/3.774852376853020208137d02/,
     2     c(3)/1.138641541510501556495d02/,
     3     c(4)/3.161123743870565596947d00/,
     4     c(5)/1.857777061846031526730d-01/
      data d(1)/2.844236833439170622273d03/,
     1     d(2)/1.282616526077372275645d03/,
     2     d(3)/2.440246379344441733056d02/,
     3     d(4)/2.360129095234412093499d01/
      data orpi/.5641895835477562869483d0/,
     1   root2/.70710678118654752440083d0/
c
      isw = 1
      y = ystar
      if (ystar.lt.-16.d0) y = -16.d0
      if (ystar.gt.16.d0) y = 16.d0
      x = -y*root2
      if(x.gt.0.d0) go to 1
      if(x.lt.0.d0) go to 2
      gauss = .5d0
      return
    2 continue
      x = - x
      isw = -1
    1 continue
      if(x.lt..477d0) go to 10
      if(x.le.4.d0) go to 20
c  evaluate erfc for x.gt.4.0
      x2 = x*x
      xm2 = 1.d0/x2
      xm4 = xm2*xm2
      xm6 = xm4*xm2
      xm8 = xm4*xm4
      xm10 = xm6*xm4
      top = p(1) + p(2)*xm2 + p(3)*xm4 + p(4)*xm6 + p(5)*xm8 + p(6)*xm10
      bot = q(1) + q(2)*xm2 + q(3)*xm4 + q(4)*xm6 + q(5)*xm8 + xm10
      crap = orpi + top/(bot*x2)
      erfc = dexp(-x2)*crap/x
c
      if(isw.eq.-1) erfc = 2.d0 - erfc
      gauss = erfc*.5d0
      return
   20 continue
c  evaluate erfc for .477.lt.x.le.4.0
      x2 = x*x
      x3 = x2*x
      x4 = x2*x2
      x5 = x3*x2
      x6 = x3*x3
      x7 = x3*x4
      x8 = x4*x4
      top = a(1) + a(2)*x + a(3)*x2 + a(4)*x3 + a(5)*x4 + a(6)*x5 +
     1 a(7)*x6 + a(8)*x7 + a(9)*x8
      bot = b(1) + b(2)*x + b(3)*x2 + b(4)*x3 + b(5)*x4 + b(6)*x5 +
     1 b(7)*x6 + b(8)*x7 + x8
      erfc = dexp(-x2)*top/bot
c
      if(isw.eq.-1) erfc = 2.d0 - erfc
      gauss = erfc*.5d0
      return
   10 continue
c  evaluate erf for x.lt..477
      x2 = x*x
      x4 = x2*x2
      x6 = x4*x2
      x8 = x4*x4
      top = c(1) + c(2)*x2 + c(3)*x4 + c(4)*x6 + c(5)*x8
      bot = d(1) + d(2)*x2 + d(3)*x4 + d(4)*x6 + x8
      erf = x*top/bot
c
      erf = erf*isw
      erfc = 1.d0 - erf
      gauss = erfc*.5d0
      return
      end
