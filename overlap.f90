module overlap_mod
use general, only : pi
use OMP_LIB
implicit none

  double complex, dimension(:,:,:), save, allocatable :: ovt
  integer, dimension(:,:,:), save, allocatable :: done

!$OMP THREADPRIVATE(done, ovt)

  private
  public overlap_driver

contains


subroutine overlap_driver(lmax,alp,vx,vy,vz,ovmat) 

  integer, intent(in) :: lmax
  double precision, intent(in) :: alp, vx, vy, vz

  double complex, dimension(0:lmax,0:lmax,0:lmax), intent(inout) :: ovmat

  integer :: l, m, n
!!$OMP THREADPRIVATE(l, m, n)


  allocate(ovt(0:lmax,0:lmax,0:lmax),done(0:lmax,0:lmax,0:lmax))

  done(:,:,:) = -1

  !WRITE(*,*)"OV inp",alp,vx,vy,vz
  do l = 0, lmax
    do m = 0, lmax
      do n = 0, lmax
         ovmat(n,m,l) = overlap(l,m,n,alp,vx,vy,vz)
         !write(*,*) "l,m,n,alp,vx,vy,vz",n,m,l,alp,vx,vy,vz
         !WRITE(*,*) "OV out",n,m,l,ovmat(n,m,l)
      enddo
    enddo
  enddo

  deallocate(ovt,done)

end subroutine overlap_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive function overlap(l,m,n,alp,vx,vy,vz) result(ov)

integer, intent(in) :: l, m, n
double precision, intent(in) :: alp, vx, vy, vz
double precision :: expo
double complex :: ov

!!$OMP THREADPRIVATE(expo)

! check that the matrix elements is not already computed
if(done(l,m,n) == 0) then
  ov = ovt(l,m,n)
  return
endif 

if(l>1) then

 ov = ( (l-1)*overlap(l-2,m,n,alp,vx,vy,vz) + dcmplx(0d0,vx)*overlap(l-1,m,n,alp,vx,vy,vz) )/(2d0*alp)

elseif(l==1) then

    if(m>1) then
         ov = ( (m-1)*overlap(l,m-2,n,alp,vx,vy,vz) + dcmplx(0d0,vy)*overlap(l,m-1,n,alp,vx,vy,vz)) /(2d0*alp)
   
    elseif(m==1) then

         if(n>1) then
              ov = ( (n-1)*overlap(l,m,n-2,alp,vx,vy,vz) + dcmplx(0d0,vz)*overlap(l,m,n-1,alp,vx,vy,vz)) /(2d0*alp)
         elseif(n==1) then
              !s(1,1,1)
              expo = -(vx**2+vy**2+vz**2)/(4d0*alp)
              ov = -dcmplx(0d0,0.125d0*vx*vy*vz*(pi/alp)**1.5d0*dexp(expo)/(alp**3))
         else
              !s(1,1,0)
              expo = -(vx**2+vy**2+vz**2)/(4d0*alp)
              ov = -0.25d0*vx*vy*(pi/alp)**1.5d0*dexp(expo)/(alp**2)
         endif

    else !m==0

         if(n>1) then
              ov = ( (n-1)*overlap(l,m,n-2,alp,vx,vy,vz) + dcmplx(0d0,vz)*overlap(l,m,n-1,alp,vx,vy,vz)) /(2d0*alp)
         elseif(n==1) then
              !s(1,0,1)
              expo = -(vx**2+vy**2+vz**2)/(4d0*alp)
              ov = -0.25d0*vx*vz*(pi/alp)**1.5d0*dexp(expo)/(alp**2)
         else ! n==0
              !s(1,0,0)
              expo = -(vx**2+vy**2+vz**2)/(4d0*alp)
              ov = dcmplx(0d0,0.5d0*vx*(pi/alp)**1.5d0*dexp(expo)/alp)
         endif

    endif

else ! l==0
    if(m>1) then
         ov = ( (m-1)*overlap(l,m-2,n,alp,vx,vy,vz) + dcmplx(0d0,vy)*overlap(l,m-1,n,alp,vx,vy,vz)) /(2d0*alp)
    elseif(m==1) then

         if(n>1) then
              ov = ( (n-1)*overlap(l,m,n-2,alp,vx,vy,vz) + dcmplx(0d0,vz)*overlap(l,m,n-1,alp,vx,vy,vz)) /(2d0*alp)
         elseif(n==1) then
              !s(0,1,1)
              expo = -(vx**2+vy**2+vz**2)/(4d0*alp)
              ov = -0.25d0*vy*vz*(pi/alp)**1.5d0*dexp(expo)/(alp**2)
         else
              !s(0,1,0)
              expo = -(vx**2+vy**2+vz**2)/(4d0*alp)
              ov = dcmplx(0d0,0.5d0*vy*(pi/alp)**1.5d0*dexp(expo)/alp)
         endif

    else !m==0
         if(n>1) then
              ov = ( (n-1)*overlap(l,m,n-2,alp,vx,vy,vz) + dcmplx(0d0,vz)*overlap(l,m,n-1,alp,vx,vy,vz)) /(2d0*alp)
         elseif(n==1) then
              !s(0,0,1)
              expo = -(vx**2+vy**2+vz**2)/(4d0*alp)
              ov = dcmplx(0d0,0.5d0*vz*(pi/alp)**1.5d0*dexp(expo)/alp)
         else
              !s(0,0,0)
              expo = -(vx**2+vy**2+vz**2)/(4d0*alp)
              ov = (pi/alp)**1.5d0*dexp(expo)
         endif
    endif
endif
         ovt(l,m,n) = ov
         done(l,m,n) = 0
end function overlap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module overlap_mod
