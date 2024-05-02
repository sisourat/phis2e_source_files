module f11mod
use general, only : pi
implicit none

contains

      subroutine computef11(n,x,f11)
! compute nth derivative of M(0.5,1.5,-x**2) = cte * x * erf(x)

      integer n
      double complex x, z, chg, i1cerf
      double complex f11(0:n), mabx(0:n)
      double precision a, b

      integer i

      if(abs(x)<50d0) then
         do i=0,n
          call cchg (0.5d0+i, 1.5d0+i, x, chg)
          f11(i) = (0.5d0/(0.5d0+i))*chg
         enddo
      else
         z = sqrt(-x)
         mabx(0) = dsqrt(pi)*i1cerf(0,z)/(2d0*z)
         f11(0) = dsqrt(pi)*i1cerf(0,z)/(2d0*z)
        if(n>0) then
         mabx(1) = 1.5d0*(exp(x)-mabx(0))/x
         f11(1) = mabx(1)/3d0
          do i = 2, n
            a = 0.5d0+(i-1)
            b = 1.5d0+(i-1)
            mabx(i) = ( b*(1d0-b+x)*mabx(i-1) + b*(b-1)*mabx(i-2) )/(a*x)
            f11(i) = (0.5d0/(0.5d0+i)) * mabx(i)
          enddo
        endif
      endif

      end subroutine computef11

      end module f11mod
