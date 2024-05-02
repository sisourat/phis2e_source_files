module nuclear_mod
use general, only : pi
use f11mod
implicit none

  double complex, dimension(:,:,:,:), save, allocatable :: nuc
  integer, dimension(:,:,:,:), save, allocatable :: done

!$OMP THREADPRIVATE(nuc, done)


  private
  public nuclear_driver

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nuclear_driver(lmax,t,gam,vx,vy,vz,bet,cx,cy,cz,tmax,potmat)

  integer, intent(in) :: lmax, t
  double precision, intent(in) :: gam, vx, vy, vz
  double precision, intent(in) :: bet, cx, cy, cz

  integer, intent(in) :: tmax

  double complex, dimension(0:lmax,0:lmax,0:lmax), intent(inout) :: potmat

  double complex :: x, expx
  double complex, dimension(:), allocatable :: f11, ft
  double precision :: rc2, rc, alp

  integer :: l, m, n, it

!!$OMP THREADPRIVATE(x, expx, f11, ft,  rc2, rc, alp, l, m, n, it )


  allocate(f11(0:tmax),ft(0:tmax))
  allocate(nuc(0:lmax,0:lmax,0:lmax,0:tmax),done(0:lmax,0:lmax,0:lmax,0:tmax))

   rc2 = cx**2+cy**2+cz**2 ; rc = dsqrt(rc2)
   alp = Gam + bet
 ! x = -B2/(4*alp)
   x = ( 0.25d0*(vx**2+vy**2+vz**2) - Gam**2*rc2 + dcmplx(0d0,Gam*(vx*cx+vy*cy+vz*cz)) )/alp
   call computef11(tmax,x,f11)

 ! x = A2/(4*alp)
   x = ( -0.25d0*(vx**2+vy**2+vz**2) + Bet**2*rc2 + dcmplx(0d0,Bet*(vx*cx+vy*cy+vz*cz)) )/alp
!   write(*,*)x
   expx = exp(-bet*rc2+x)

   do it = 0, tmax
      ft(it) = expx*f11(it)
!      write(101,'(i5,100(f20.16,1X))')it,f11(it),expx,gam,x
!      write(101,*)x
   enddo

  done(:,:,:,:) = -1
! nuc and done are used to store intermediate integrals (see end of subroutine nuclear)
! avoid recomputing them => huge speed-up


  do l = 0, lmax
    do m = 0, lmax
      do n = 0, lmax

         potmat(n,m,l) = nuclear(t,l,m,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft)
!         write(100,'(3(i3,1X),100(f20.16,1X))')n,m,l,potmat(n,m,l),Gam,bet

      enddo
    enddo
  enddo

  deallocate(f11,ft)
  deallocate(nuc,done)

end subroutine nuclear_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


recursive function nuclear(t,l,m,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) result(pot)
use general
implicit none

integer, intent(in) :: t, l, m, n
double precision, intent(in) :: gam, vx, vy, vz
double precision, intent(in) :: bet, cx, cy, cz

integer, intent(in) :: tmax
double complex, dimension(0:tmax), intent(in) :: ft

double complex :: mx, my, mz, px, py, pz

!!$OMP  THREADPRIVATE( mx, my, mz, px, py, pz )

double precision :: expo, alp
double complex :: pot

!!$OMP  THREADPRIVATE(expo, alp)
 

! check that the matrix elements is not already computed
if(done(l,m,n,t) == 0) then
  pot = nuc(l,m,n,t)
  return
endif 

alp = gam + bet

if(l>1) then
 
 pot = -0.5d0/alp*( (l-1)*(-nuclear(t,l-2,m,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) + nuclear(t+1,l-2,m,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) ) &
  &- dcmplx(2d0*bet*cx,vx) * nuclear(t,l-1,m,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) - dcmplx(2d0*gam*cx,-vx) * nuclear(t+1,l-1,m,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) ) 
      
elseif(l==1) then

    if(m>1) then
   
        pot = -0.5d0/alp*( (m-1)*(-nuclear(t,l,m-2,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) + nuclear(t+1,l,m-2,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) ) &
          &- dcmplx(2d0*bet*cy,vy) * nuclear(t,l,m-1,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) - dcmplx(2d0*gam*cy,-vy) * nuclear(t+1,l,m-1,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) )

    elseif(m==1) then

         if(n>1) then
  
          pot = -0.5d0/alp*( (n-1)*(-nuclear(t,l,m,n-2,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) + nuclear(t+1,l,m,n-2,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) ) &
          &- dcmplx(2d0*bet*cz,vz) * nuclear(t,l,m,n-1,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) - dcmplx(2d0*gam*cz,-vz) * nuclear(t+1,l,m,n-1,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) )

         elseif(n==1) then
!           Vt(1,1,1)
            mx = dcmplx(-vx,2d0*bet*cx)
            my = dcmplx(-vy,2d0*bet*cy)
            mz = dcmplx(-vz,2d0*bet*cz)

            px = dcmplx(vx,2d0*gam*cx)
            py = dcmplx(vy,2d0*gam*cy)
            pz = dcmplx(vz,2d0*gam*cz)

            pot = dcmplx(0d0,0.25d0*pi/alp**4)*( mx*my*mz*ft(t) + ( mz*(mx*py+px*my) + pz*mx*my )*ft(t+1) + ( pz*(mx*py+px*my) + mz*px*py )*ft(t+2) &
            &+ px*py*pz*ft(t+3) )
         else ! n = 0
!           Vt(1,1,0)
            mx = dcmplx(-vx,2d0*bet*cx)
            my = dcmplx(-vy,2d0*bet*cy)

            px = dcmplx(vx,2d0*gam*cx)
            py = dcmplx(vy,2d0*gam*cy)

            pot = -0.5d0*pi/alp**3*( mx*my*ft(t) + (mx*py+px*my)*ft(t+1) + px*py*ft(t+2) ) 
         endif

    else !m==0

         if(n>1) then

          pot = -0.5d0/alp*( (n-1)*(-nuclear(t,l,m,n-2,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) + nuclear(t+1,l,m,n-2,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) ) &
          &- dcmplx(2d0*bet*cz,vz) * nuclear(t,l,m,n-1,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) - dcmplx(2d0*gam*cz,-vz) * nuclear(t+1,l,m,n-1,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) )

         elseif(n==1) then
!           Vt(1,0,1)
            mx = dcmplx(-vx,2d0*bet*cx)
            mz = dcmplx(-vz,2d0*bet*cz)

            px = dcmplx(vx,2d0*gam*cx)
            pz = dcmplx(vz,2d0*gam*cz)

            pot = -0.5d0*pi/alp**3*( mx*mz*ft(t) + (mx*pz+px*mz)*ft(t+1) + px*pz*ft(t+2) ) 
         else ! n==0
!           Vt(1,0,0)
            pot = dcmplx(0d0,-pi/alp**2)*( dcmplx(-vx,2d0*bet*cx)*ft(t) + dcmplx(vx,2d0*gam*cx)*ft(t+1) )
         endif

    endif

else ! l==0

    if(m>1) then

     pot = -0.5d0/alp*( (m-1)*(-nuclear(t,l,m-2,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) + nuclear(t+1,l,m-2,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) ) &
        &- dcmplx(2d0*bet*cy,vy) * nuclear(t,l,m-1,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) - dcmplx(2d0*gam*cy,-vy) * nuclear(t+1,l,m-1,n,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) )

    elseif(m==1) then

         if(n>1) then

          pot = -0.5d0/alp*( (n-1)*(-nuclear(t,l,m,n-2,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) + nuclear(t+1,l,m,n-2,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) ) &
          &- dcmplx(2d0*bet*cz,vz) * nuclear(t,l,m,n-1,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) - dcmplx(2d0*gam*cz,-vz) * nuclear(t+1,l,m,n-1,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) )

         elseif(n==1) then
!           Vt(0,1,1)
            my = dcmplx(-vy,2d0*bet*cy)
            mz = dcmplx(-vz,2d0*bet*cz)

            py = dcmplx(vy,2d0*gam*cy)
            pz = dcmplx(vz,2d0*gam*cz)

            pot = -0.5d0*pi/alp**3*( my*mz*ft(t) + (my*pz+py*mz)*ft(t+1) + py*pz*ft(t+2) ) 
         else
!           Vt(0,1,0)
            pot = dcmplx(0d0,-pi/alp**2)*( dcmplx(-vy,2d0*bet*cy)*ft(t) + dcmplx(vy,2d0*gam*cy)*ft(t+1) )
         endif

    else !m==0

         if(n>1) then

          pot = -0.5d0/alp*( (n-1)*(-nuclear(t,l,m,n-2,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) + nuclear(t+1,l,m,n-2,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) ) &
          &- dcmplx(2d0*bet*cz,vz) * nuclear(t,l,m,n-1,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) - dcmplx(2d0*gam*cz,-vz) * nuclear(t+1,l,m,n-1,gam,vx,vy,vz,bet,cx,cy,cz,tmax,ft) )

         elseif(n==1) then
!           Vt(0,0,1)
            pot = dcmplx(0d0,-pi/alp**2)*( dcmplx(-vz,2d0*bet*cz)*ft(t) + dcmplx(vz,2d0*gam*cz)*ft(t+1) )
         else
!           Vt(0,0,0)
            pot = 2d0*pi/alp*ft(t)
         endif
    endif
endif

         nuc(l,m,n,t) = pot
         done(l,m,n,t) = 0

end function nuclear

end module nuclear_mod
