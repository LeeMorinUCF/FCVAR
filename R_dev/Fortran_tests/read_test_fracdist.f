
      implicit real*8 (a-h,o-z)
c
c read_test_fracdist.f reads test cases for fracdist.f, which is a
c program to compute P values or critical values for tests of
c fractional unit roots and cointegration.
c Version 1 (May, 14 2020) that only tests file IO.
c Copyright (c) Lealand Morin, 2020
c
c Not sure why these are 'unexpected data declaration statements'
      character*16 pval_header
c      character*18 cval_header
c
c Set parameters for number of test cases.
c      npval = 2400
      npval = 24
      ncval = 720
c
c Open file to read test cases for p-values.
c
      open(unit=2,file='test_fpval.txt',status='OLD',err=333)
      go to 334
 333  write(6,*) 'File of test cases not opened.'
      stop
 334  continue
c
c Open file to save results.
c
      open (unit=13,file='test_fpval.out')
c
c Extract parameters for test case in each line.
c Evaluate p-values.
c Write result to output file.
c
      read(2,200) pval_header
 200    format(a16)
      write(6,*) 'File of test cases contains the following variables:'
      write(6,*) pval_header
      write(13,101) pval_header, 'pval'
 101  format(a16,1x,a4)
      do ip=1,npval
        read(2,201) iscon, iq, bb, stat
 201    format(i1,1x,i2,1x,f5.3,1x,f8.4)
        isave = 0
        ipc = 0
        clevel = 0
        ccrit = 0
c        pval = 0.12345678
        call frval(iq,isave,ipc,iscon,bb,stat,pval,clevel,ccrit)
        write(13,102) iscon, iq, bb, stat, pval
 102    format(i1,1x,i2,1x,f5.3,1x,f8.4,1x,f6.4)
      end do
      write(6,*) 'See output in file test_fpval.out'
      stop
      end
c
c Subroutines from fracdist.f:
c

c
c End
c
