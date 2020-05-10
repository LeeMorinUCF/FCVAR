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
