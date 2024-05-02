!------------------------Created by J.W.GAO, 4th Jan. 2017------------------

module repulsion_mod
use general, only : pi
use f11mod
use OMP_lib
implicit none

  double complex, dimension(:,:,:,:,:,:,:),save, allocatable :: repul
  integer, dimension(:,:,:,:,:,:,:),save, allocatable :: done

!$OMP THREADPRIVATE( repul,done )


  private
  public repulsion_driver

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine repulsion_driver(lmax,t,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,repulmat)
!use general, only : pi
!use omp_lib
  integer, intent(in) :: lmax, t
  double precision, intent(in) :: alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz
  double precision, intent(in) :: cx, cy, cz

  integer, intent(in) :: tmax

  double complex, dimension(0:lmax, 0:lmax, 0:lmax, 0:lmax, 0:lmax, 0:lmax), intent(inout) :: repulmat

  double complex :: x, expx
  double complex, dimension(:),allocatable :: f11, ft
  double precision :: rc2, rc, alp, bet
  double precision :: a2, b2, v2


!!$OMP THREADPRIVATE( f11, ft)

!!$OMP THREADPRIVATE(  x, expx,rc2, rc, alp, bet,a2, b2, v2 )

  integer :: la, ma, na, lb, nb, mb, it

!!$OMP THREADPRIVATE( la, ma, na, lb, nb, mb, it)


  allocate(f11(0:tmax),ft(0:tmax))
  allocate( repul( 0:lmax, 0:lmax, 0:lmax, 0:lmax, 0:lmax, 0:lmax, 0:tmax) )
  allocate( done( 0:lmax, 0:lmax, 0:lmax, 0:lmax, 0:lmax, 0:lmax, 0:tmax)  )

!  repul=0.d0

   a2 = ax**2+ay**2+az**2 
   b2 = bx**2+by**2+bz**2 
   v2 = vx**2+vy**2+vz**2 
   rc2 = cx**2+cy**2+cz**2 ; rc = dsqrt(rc2)

!----x = -bet2/alp------------------------------------------------

   x=( alp_b*a2/(4.d0*alp_a*(alp_a+alp_b)) )+( alp_a*b2/(4.d0*alp_b*(alp_a+alp_b)) )-((ax*bx+ay*by+az*bz)/(2.d0*(alp_a+alp_b))  )&
     &+v2*(alp_a+alp_b)/(4.d0*alp_a*alp_b) + (ax*vx+ay*vy+az*vz)/(2.d0*alp_a) - (bx*vx+by*vy+bz*vz)/(2.d0*alp_b) &
     &-(alp_a*alp_b*rc2)/(alp_a+alp_b) +  dcmplx(0d0,alp_b*(ax*cx+ay*cy+az*cz)/(alp_a+alp_b) )&
     &- dcmplx(0d0,alp_a*(bx*cx+by*cy+bz*cz)/(alp_a+alp_b) ) + dcmplx(0d0,(cx*vx+cy*vy+cz*vz) ) 
    
    
   call computef11(tmax,x,f11)


!---alp=1/4alp_a + 1/4alp_b------------------------------------------------------------------------------

   alp=1.d0/(4.d0*alp_a)+1.d0/(4.d0*alp_b)

!---x = B2/(4*alp_a)+D2/(4*alp_b)!----- -ibR-------------------------------------------------------------

   x=-0.25d0*(a2+2.d0*(ax*vx+ay*vy+az*vz)+v2)/alp_a +0.25d0*(4.d0*rc2*alp_b**2.d0 -b2 - v2 + 2.d0*(bx*vx+by*vy+bz*vz)&
     & + dcmplx(0d0,(4.d0*alp_b*(cx*(bx-vx)+cy*(by-vy)+cz*(bz-vz)))) ) /alp_b! - dcmplx(0d0,(bx*cx+by*cy+bz*cz) )

!--------------------------------------------------------------------------------exp(iar1+ibr2)----------
                   

   expx =( pi**2.5d0/( alp**0.5d0 * (alp_a*alp_b)**1.5d0) ) * exp(-alp_b*rc2+x)
!  write(*,*) expx


  do it = 0, tmax
     ft(it) = expx*f11(it)
  enddo


  done(:,:,:,:,:,:, 0:max( 0, 2*(lmax-1)) ) = -1  ! The sum of x,y,z exponents is L.....



  do la = 0,lmax
    do ma = 0, lmax
      do na = 0, lmax
        do lb = 0, lmax
          do mb = 0, lmax
            do nb = 0, lmax

             repulmat(nb,mb,lb,na,ma,la) = repulsion(t,la,ma,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)
        
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo


  deallocate(f11,ft)
  deallocate(repul,done)

end subroutine repulsion_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive function repulsion(t,la,ma,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) result(rep)

use general
implicit none

integer, intent(in) :: t, la, ma, na, lb, mb, nb

double precision, intent(in) :: alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz
double precision, intent(in) :: cx, cy, cz


integer, intent(in) :: tmax
double complex, dimension(0:tmax), intent(in) :: ft

double complex :: mx, my, mz, px, py, pz

!!$OMP THREADPRIVATE( mx, my, mz, px, py, pz )


double precision :: Acx,Acy,Acz,Ccx,Ccy,Ccz
double complex :: Bcx,Bcy,Bcz,Dcx,Dcy,Dcz
 
!!$OMP THREADPRIVATE( Acx,Acy,Acz,Ccx,Ccy,Ccz, Bcx,Bcy,Bcz,Dcx,Dcy,Dcz)


double complex :: rep

!--For any term involving "-1, -2, ..." should be 0. "-" comes from cross term in the recursion relationship--------------

if( (la.eq.-1).or.(ma.eq.-1).or.(na.eq.-1).or.(lb.eq.-1).or.(mb.eq.-1).or.(nb.eq.-1) ) then   

   rep = 0.d0
   return

endif


!--check that the matrix elements is not already computed---

if(done(la,ma,na,lb,mb,nb,t) == 0) then

  rep = repul(la,ma,na,lb,mb,nb,t)
  return

endif 

Acx= (ax+vx)
Acy= (ay+vy)
Acz= (az+vz)
Ccx=-vx+bx
Ccy=-vy+by
Ccz=-vz+bz

Bcx=  dcmplx( (alp_b*ax-alp_a*bx+vx*(alp_a+alp_b)), 2.d0*alp_a*alp_b*cx )
Bcy=  dcmplx( (alp_b*ay-alp_a*by+vy*(alp_a+alp_b)), 2.d0*alp_a*alp_b*cy )
Bcz=  dcmplx( (alp_b*az-alp_a*bz+vz*(alp_a+alp_b)), 2.d0*alp_a*alp_b*cz )
Dcx=  -Bcx   !-dcmplx( (alp_b*ax-alp_a*bx+vx*(alp_a+alp_b)), 2.d0*alp_a*alp_b*cx )
Dcy=  -Bcy   !-dcmplx( (alp_b*ay-alp_a*by+vy*(alp_a+alp_b)), 2.d0*alp_a*alp_b*cy )
Dcz=  -Bcz   !-dcmplx( (alp_b*az-alp_a*bz+vz*(alp_a+alp_b)), 2.d0*alp_a*alp_b*cz )





if(la>1) then

   rep = 0.5d0*(la-1)/alp_a * ( repulsion(t,la-2,ma,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
       &alp_b/(alp_a+alp_b) * (repulsion(t+1,la-2,ma,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
           &+( (alp_a/(alp_a+alp_b))*(lb/(la-1.d0)) *repulsion(t+1,la-1,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
       &-1.d0/(2.d0*alp_a) * ( dcmplx(0.d0, -Acx )* repulsion(t,la-1,ma,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
       &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Bcx* repulsion(t+1,la-1,ma,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )


elseif(la==1) then
 
     if(ma>1) then

   rep = 0.5d0*(ma-1)/alp_a * ( repulsion(t,la,ma-2,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
       &alp_b/(alp_a+alp_b) * (repulsion(t+1,la,ma-2,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
           &+( (alp_a/(alp_a+alp_b))*(mb/(ma-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
       &-1.d0/(2.d0*alp_a) * ( dcmplx(0.d0, -Acy )* repulsion(t,la,ma-1,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
       &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Bcy* repulsion(t+1,la,ma-1,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
        
     elseif(ma==1) then

          if(na>1) then
        rep = 0.5d0*(na-1)/alp_a * ( repulsion(t,la,ma,na-2,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_b/(alp_a+alp_b) * (repulsion(t+1,la,ma,na-2,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_a/(alp_a+alp_b))*(nb/(na-1.d0)) *repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_a) * ( dcmplx(0.d0, -Acz )* repulsion(t,la,ma,na-1,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Bcz* repulsion(t+1,la,ma,na-1,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

          elseif(na==1) then

                if(lb>1)then

        rep = 0.5d0*(lb-1)/alp_b * ( repulsion(t,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(la/(lb-1.d0)) *repulsion(t+1,la-1,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccx )* repulsion(t,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcx* repulsion(t+1,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                elseif(lb==1) then

                      if(mb>1)then

        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                            elseif(nb==1) then 

        rep = -1.d0/(64.d0*(alp_a*alp_b)**3.d0) *( Acx*Acy*Acz*Ccx*Ccy*Ccz*ft(t) &
            &-ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Acz*(Ccx*Ccy*Dcz+Ccy*Ccz*Dcx+Ccz*Ccx*Dcy )+Ccx*Ccy*Ccz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy)&
                 & + 2.d0*alp_a*alp_b*(Acx*Acy*Ccx*Ccy+Acx*Acz*Ccx*Ccz +Acz*Acy*Ccz*Ccy  ) ) &
            &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Acx*Acy*Acz*(Ccx*Dcy*Dcz+Ccy*Dcz*Dcx+Ccz*Dcx*Dcy) &
                &+(Ccx*Ccy*Dcz+Ccy*Ccz*Dcx+Ccz*Ccx*Dcy )*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) &
                &+Ccx*Ccy*Ccz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy) &
                 &+2.d0*alp_a*alp_b*(Acx*Ccx*(Acy*Dcy+Acz*Dcz+Bcy*Ccy+Bcz*Ccz)+Acy*Ccy*(Acx*Dcx+Acz*Dcz+Bcx*Ccx+Bcz*Ccz)&
                 &+Acz*Ccz*(Acy*Dcy+Acx*Dcx+Bcy*Ccy+Bcx*Ccx)  )+(2.d0*alp_a*alp_b)**2.d0 *(Acx*Ccx+Acy*Ccy+Acz*Ccz)      ) &
            &-ft(t+3)/(alp_a+alp_b)**3.d0 * ((Ccx*Dcy*Dcz+Ccy*Dcz*Dcx+Ccz*Dcx*Dcy)*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy)&
                &+(Ccx*Ccy*Dcz+Ccy*Ccz*Dcx+Ccz*Ccx*Dcy )*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy) +Acx*Acy*Acz*Dcx*Dcy*Dcz+ Ccx*Ccy*Ccz*Bcx*Bcy*Bcz &
                 &+2.d0*alp_a*alp_b*((Bcz*Ccz+Acz*Dcz)*(Acy*Dcy+Acx*Dcx+Bcy*Ccy+Bcx*Ccx) + (Acx*Dcx+Bcx*Ccx)*(Acy*Dcy+Bcy*Ccy) &
                 &+ Bcx*Dcx*(Acy*Ccy+Acz*Ccz) + Bcy*Dcy*(Acx*Ccx+Acz*Ccz) +Bcz*Dcz*(Acy*Ccy+Acx*Ccx)  ) &
                 &+ (2.d0*alp_a*alp_b)**2.d0 *(Acx*Dcx+Bcx*Ccx+Acy*Dcy+Bcy*Ccy+Acz*Dcz+Bcz*Ccz ) +(2.d0*alp_a*alp_b)**3.d0   )&
            &+ft(t+4)/(alp_a+alp_b)**4.d0 * (Bcx*Bcy*Bcz*(Ccx*Ccy*Dcz+Ccy*Ccz*Dcx+Ccz*Ccx*Dcy ) &
                &+ (Ccx*Dcy*Dcz+Ccy*Dcz*Dcx+Ccz*Dcx*Dcy)*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)&
                &+Dcx*Dcy*Dcz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy)  &
                  &+2.d0*alp_a*alp_b*(Bcx*Dcx*(Acy*Dcy+Acz*Dcz+Bcy*Ccy+Bcz*Ccz)+Bcy*Dcy*(Acx*Dcx+Acz*Dcz+Bcx*Ccx+Bcz*Ccz)&
                 &+Bcz*Dcz*(Acy*Dcy+Acx*Dcx+Bcy*Ccy+Bcx*Ccx)  )+(2.d0*alp_a*alp_b)**2.d0 *(Bcx*Dcx+Bcy*Dcy+Bcz*Dcz)   ) &
            &-ft(t+5)/(alp_a+alp_b)**5.d0 *(Dcx*Dcy*Dcz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)+Bcx*Bcy*Bcz*(Ccx*Dcy*Dcz+Ccy*Dcz*Dcx+Ccz*Dcx*Dcy) &
               & + 2.d0*alp_a*alp_b*(Bcx*Bcy*Dcx*Dcy+Bcx*Bcz*Dcx*Dcz +Bcz*Bcy*Dcz*Dcy  )  )&
            &+ft(t+6)/(alp_a+alp_b)**6.d0 * Bcx*Bcy*Bcz*Dcx*Dcy*Dcz  )        ! Vt(1,1,1,1,1,1)

                            else ! nb = 0

        rep = -dcmplx(0.d0, 1.d0/(32.d0*alp_a*alp_a*alp_a*alp_b*alp_b))*(-Acx*Acy*Acz*Ccx*Ccy*ft(t) &
            &+ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Acz*Ccy*Dcx+Acx*Acy*Acz*Ccx*Dcy+Ccx*Ccy*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) &
               &+2.d0*alp_a*alp_b*Acz*(Acx*Ccx+Acy*Ccy)    ) &
            &-ft(t+2)/(alp_a+alp_b)**2.d0 * (Acx*Acy*Acz*Dcx*Dcy +Ccx*Dcy*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy)&
              &+Ccx*Ccy*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)+Dcx*Ccy*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy)  &
               &+2.d0*alp_a*alp_b*Bcz*(Acx*Ccx+Acy*Ccy) +Acz*(2.d0*alp_a*alp_b*(Acx*Dcx+Acy*Dcy+Bcx*Ccx+Bcy*Ccy)+4.d0*(alp_a*alp_b)**2.d0 ) ) &
            &+ft(t+3)/(alp_a+alp_b)**3.d0 * (Ccx*Dcy*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)+Dcx*Dcy*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy)&
              &+Ccx*Ccy*Bcx*Bcy*Bcz+ Dcx*Ccy*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)  &
               &+2.d0*alp_a*alp_b*Acz*(Bcx*Dcx+Bcy*Dcy) +Bcz*(2.d0*alp_a*alp_b*(Acx*Dcx+Acy*Dcy+Bcx*Ccx+Bcy*Ccy)+4.d0*(alp_a*alp_b)**2.d0 ) ) &
            &-ft(t+4)/(alp_a+alp_b)**4.d0 * ( Ccx*Dcy*Bcx*Bcy*Bcz + Dcx*Dcy*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy) + Dcx*Ccy*Bcx*Bcy*Bcz &
                            &+2.d0*alp_a*alp_b*Bcz*(Bcx*Dcx+Bcy*Dcy)  )&
            &+ft(t+5)/(alp_a+alp_b)**5.d0 * Dcx*Dcy*Bcx*Bcy*Bcz )             ! Vt(1,1,1,1,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 

        rep = -dcmplx(0.d0, 1.d0/(32.d0*alp_a*alp_a*alp_a*alp_b*alp_b))*(-Acx*Acy*Acz*Ccx*Ccz*ft(t) &
            &+ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Acz*Ccz*Dcx+Acx*Acy*Acz*Ccx*Dcz+Ccx*Ccz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) &
                    &+2.d0*alp_a*alp_b*Acy*(Acx*Ccx+Acz*Ccz)  ) &
            &-ft(t+2)/(alp_a+alp_b)**2.d0 * (Acx*Acy*Acz*Dcx*Dcz +Ccx*Dcz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy)&
              &+Ccx*Ccz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)+Dcx*Ccz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) &
                 &+2.d0*alp_a*alp_b*Bcy*(Acx*Ccx+Acz*Ccz) +Acy*(2.d0*alp_a*alp_b*(Acx*Dcx+Acz*Dcz+Bcx*Ccx+Bcz*Ccz)+4.d0*(alp_a*alp_b)**2.d0 )  ) &
            &+ft(t+3)/(alp_a+alp_b)**3.d0 * (Ccx*Dcz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)+Dcx*Dcz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy)&
              &+Ccx*Ccz*Bcx*Bcy*Bcz+ Dcx*Ccz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)  &
                 &+2.d0*alp_a*alp_b*Acy*(Bcx*Dcx+Bcz*Dcz) +Bcy*(2.d0*alp_a*alp_b*(Acx*Dcx+Acz*Dcz+Bcx*Ccx+Bcz*Ccz)+4.d0*(alp_a*alp_b)**2.d0 ) ) &
            &-ft(t+4)/(alp_a+alp_b)**4.d0 * ( Ccx*Dcz*Bcx*Bcy*Bcz + Dcx*Dcz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy) + Dcx*Ccz*Bcx*Bcy*Bcz &
                            &+2.d0*alp_a*alp_b*Bcy*(Bcx*Dcx+Bcz*Dcz) )&
            &+ft(t+5)/(alp_a+alp_b)**5.d0 * Dcx*Dcz*Bcx*Bcy*Bcz )                             !Vt(1,1,1,1,0,1)

                            else ! nb = 0

        rep = 1.d0/(16.d0*alp_a*alp_a*alp_a*alp_b) * ( Acx*Acy*Acz*Ccx*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Acz*Dcx + Ccx*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) + 2.d0*alp_a*alp_b*Acy*Acz) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccx*(Bcx*Bcy*Acz+Bcy*Bcz*Acx+Bcz*Bcx*Acy)+Dcx*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) &
                & + 2.d0*alp_a*alp_b*(Acy*Bcz+Acz*Bcy) ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcx*Bcy*Bcz*Ccx + Dcx*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)  + 2.d0*alp_a*alp_b*Bcy*Bcz) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcx*Bcy*Bcz*Dcx  )                            !Vt(1,1,1,1,0,0)

                            endif
                      endif
                else   !lb==0

                      if(mb>1)then
        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
             &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                 &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
             &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
             &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 

        rep = -dcmplx(0.d0, 1.d0/(32.d0*alp_a*alp_a*alp_a*alp_b*alp_b))*(-Acx*Acy*Acz*Ccy*Ccz*ft(t) &
            &+ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Acz*Ccz*Dcy+Acx*Acy*Acz*Ccy*Dcz+Ccy*Ccz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) &
               &+2.d0*alp_a*alp_b*Acx*(Acz*Ccz+Acy*Ccy)  ) &
            &-ft(t+2)/(alp_a+alp_b)**2.d0 * (Acx*Acy*Acz*Dcy*Dcz +Ccy*Dcz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy)&
              &+Ccy*Ccz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)+Dcy*Ccz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) &
                &+2.d0*alp_a*alp_b*Bcx*(Acz*Ccz+Acy*Ccy) +Acx*(2.d0*alp_a*alp_b*(Acz*Dcz+Acy*Dcy+Bcz*Ccz+Bcy*Ccy)+4.d0*(alp_a*alp_b)**2.d0 ) ) &
            &+ft(t+3)/(alp_a+alp_b)**3.d0 * (Ccy*Dcz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)+Dcy*Dcz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy)&
              &+Ccy*Ccz*Bcx*Bcy*Bcz+ Dcy*Ccz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy) &
               &+2.d0*alp_a*alp_b*Acx*(Bcz*Dcz+Bcy*Dcy) +Bcx*(2.d0*alp_a*alp_b*(Acz*Dcz+Acy*Dcy+Bcz*Ccz+Bcy*Ccy)+4.d0*(alp_a*alp_b)**2.d0 )   ) &
            &-ft(t+4)/(alp_a+alp_b)**4.d0 * ( Ccy*Dcz*Bcx*Bcy*Bcz + Dcy*Dcz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy) + Dcy*Ccz*Bcx*Bcy*Bcz &
                             &+2.d0*alp_a*alp_b*Bcx*(Bcz*Dcz+Bcy*Dcy)  )&
            &+ft(t+5)/(alp_a+alp_b)**5.d0 * Dcy*Dcz*Bcx*Bcy*Bcz )                               ! Vt(1,1,1,0,1,1)
                            else ! nb = 0
        rep = 1.d0/(16.d0*alp_a*alp_a*alp_a*alp_b) * ( Acx*Acy*Acz*Ccy*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Acz*Dcy + Ccy*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) +2.d0*alp_a*alp_b*Acx*Acz ) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccy*(Bcx*Bcy*Acz+Bcy*Bcz*Acx+Bcz*Bcx*Acy)+Dcy*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) &
               &+ 2.d0*alp_a*alp_b*(Acx*Bcz+Acz*Bcx) ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcx*Bcy*Bcz*Ccy + Dcy*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy)  + 2.d0*alp_a*alp_b*Bcx*Bcz ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcx*Bcy*Bcz*Dcy  )                               ! Vt(1,1,1,0,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then
        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = 1.d0/(16.d0*alp_a*alp_a*alp_a*alp_b) * ( Acx*Acy*Acz*Ccz*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Acz*Dcz + Ccz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) +2.d0*alp_a*alp_b*Acy*Acx  ) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccz*(Bcx*Bcy*Acz+Bcy*Bcz*Acx+Bcz*Bcx*Acy)+Dcz*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) &
               &+2.d0*alp_a*alp_b*(Acx*Bcy+Acy*Bcx)  ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcx*Bcy*Bcz*Ccz + Dcz*(Acx*Bcy*Bcz+Acy*Bcz*Bcx+Acz*Bcx*Bcy) +2.d0*alp_a*alp_b*Bcy*Bcx  ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcx*Bcy*Bcz*Dcz  )                               !Vt(1,1,1,0,0,1)
                            else ! nn = 0
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_a*alp_a) ) * ( -Acx*Acy*Acz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acx*Acy*Bcz+Acy*Acz*Bcx+Acz*Acx*Bcy) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcx*Bcy*Acz+Bcy*Bcz*Acx+Bcz*Bcx*Acy) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcx*Bcy*Bcz )                          !Vt(1,1,1,0,0,0)
          
                            endif
                      endif
                endif
          else !na==0

                if(lb>1)then

         rep = 0.5d0*(lb-1)/alp_b * ( repulsion(t,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(la/(lb-1.d0)) *repulsion(t+1,la-1,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccx )* repulsion(t,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcx* repulsion(t+1,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                elseif(lb==1) then

                      if(mb>1)then

      rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = -dcmplx(0.d0, 1.d0/(32.d0*alp_a*alp_a*alp_b*alp_b*alp_b))*(-Acx*Acy*Ccz*Ccx*Ccy*ft(t) &
            &+ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Ccz*Ccy*Dcx+Acx*Acy*Ccz*Ccx*Dcy+Ccx*Ccy*(Acx*Acy*Dcz+Acy*Ccz*Bcx+Ccz*Acx*Bcy) &
               &+2.d0*alp_a*alp_b*Ccz*(Acx*Ccx+Acy*Ccy)  ) &
            &-ft(t+2)/(alp_a+alp_b)**2.d0 * (Acx*Acy*Ccz*Dcx*Dcy +Ccx*Dcy*(Acx*Acy*Dcz+Acy*Ccz*Bcx+Ccz*Acx*Bcy)&
              &+Ccx*Ccy*(Acx*Bcy*Dcz+Acy*Dcz*Bcx+Ccz*Bcx*Bcy)+Dcx*Ccy*(Acx*Acy*Dcz+Acy*Ccz*Bcx+Ccz*Acx*Bcy) &
              &+2.d0*alp_a*alp_b*Dcz*(Acx*Ccx+Acy*Ccy) +Ccz*(2.d0*alp_a*alp_b*(Acx*Dcx+Acy*Dcy+Bcx*Ccx+Bcy*Ccy)+4.d0*(alp_a*alp_b)**2.d0 ) ) &
            &+ft(t+3)/(alp_a+alp_b)**3.d0 * (Ccx*Dcy*(Acx*Bcy*Dcz+Acy*Dcz*Bcx+Ccz*Bcx*Bcy)+Dcx*Dcy*(Acx*Acy*Dcz+Acy*Ccz*Bcx+Ccz*Acx*Bcy)&
              &+Ccx*Ccy*Bcx*Bcy*Dcz+ Dcx*Ccy*(Acx*Bcy*Dcz+Acy*Dcz*Bcx+Ccz*Bcx*Bcy) &
               &+2.d0*alp_a*alp_b*Ccz*(Bcx*Dcx+Bcy*Dcy) +Dcz*(2.d0*alp_a*alp_b*(Acx*Dcx+Acy*Dcy+Bcx*Ccx+Bcy*Ccy)+4.d0*(alp_a*alp_b)**2.d0 )  ) &
            &-ft(t+4)/(alp_a+alp_b)**4.d0 * ( Ccx*Dcy*Bcx*Bcy*Dcz + Dcx*Dcy*(Acx*Bcy*Dcz+Acy*Dcz*Bcx+Ccz*Bcx*Bcy) + Dcx*Ccy*Bcx*Bcy*Dcz &
                           &+2.d0*alp_a*alp_b*Dcz*(Bcx*Dcx+Bcy*Dcy)  )&
            &+ft(t+5)/(alp_a+alp_b)**5.d0 * Dcx*Dcy*Bcx*Bcy*Dcz )                              ! Vt(1,1,0,1,1,1)

                            else ! nb = 0
        rep = 1.d0/(16.d0*alp_a*alp_a*alp_b*alp_b) * ( Acx*Acy*Ccy*Ccx*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Ccy*Dcx + Ccx*(Acx*Acy*Dcy+Acy*Ccy*Bcx+Ccy*Acx*Bcy) +2.d0*alp_a*alp_b*(Acx*Ccx +Acy*Ccy) ) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccx*(Bcx*Bcy*Ccy+Bcy*Dcy*Acx+Dcy*Bcx*Acy)+Dcx*(Acx*Acy*Dcy+Acy*Ccy*Bcx+Ccy*Acx*Bcy) &
              &+ 2.d0*alp_a*alp_b*(Acx*Dcx +Acy*Dcy+ Bcx*Ccx+Bcy*Ccy) +(2.d0*alp_a*alp_b)**2.d0 ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcx*Bcy*Dcy*Ccx + Dcx*(Acx*Bcy*Dcy+Acy*Dcy*Bcx+Ccy*Bcx*Bcy) + 2.d0*alp_a*alp_b*(Bcx*Dcx +Bcy*Dcy) ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcx*Bcy*Dcy*Dcx  )                                ! Vt(1,1,0,1,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = 1.d0/(16.d0*alp_a*alp_a*alp_b*alp_b) * ( Acx*Acy*Ccz*Ccx*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Ccz*Dcx + Ccx*(Acx*Acy*Dcz+Acy*Ccy*Bcx+Ccz*Acx*Bcy) +2.d0*alp_a*alp_b*Acy*Ccz) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccx*(Bcx*Bcy*Ccz+Bcy*Dcz*Acx+Dcz*Bcx*Acy)+Dcx*(Acx*Acy*Dcz+Acy*Ccz*Bcx+Ccz*Acx*Bcy) &
              &+2.d0*alp_a*alp_b*(Acy*Dcz+Ccz*Bcy) ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcx*Bcy*Dcz*Ccx + Dcx*(Acx*Bcy*Dcz+Acy*Dcz*Bcx+Ccz*Bcx*Bcy) +2.d0*alp_a*alp_b*Bcy*Dcz ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcx*Bcy*Dcz*Dcx  )                               !Vt(1,1,0,1,0,1)
                            else ! nn = 0
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_a*alp_b) ) * ( -Acx*Acy*Ccx*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acx*Acy*Dcx+Acy*Ccx*Bcx+Ccx*Acx*Bcy+2.d0*alp_a*alp_b*Acy) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcx*Bcy*Ccx+Bcy*Dcx*Acx+Dcx*Bcx*Acy&
               &+2.d0*alp_a*alp_b*Bcy) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcx*Bcy*Dcx )                   !Vt(1,1,0,1,0,0)
                            endif
                       endif
                else   !lb==0

                      if(mb>1)then

        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                            elseif(nb==1) then 
        rep = 1.d0/(16.d0*alp_a*alp_a*alp_b*alp_b) * ( Acx*Acy*Ccz*Ccy*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acx*Acy*Ccz*Dcy + Ccy*(Acx*Acy*Dcz+Acy*Ccy*Bcx+Ccz*Acx*Bcy) +2.d0*alp_a*alp_b*Acx*Ccz) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccy*(Bcx*Bcy*Ccz+Bcy*Dcz*Acx+Dcz*Bcx*Acy)+Dcy*(Acx*Acy*Dcz+Acy*Ccz*Bcx+Ccz*Acx*Bcy)&
              &+2.d0*alp_a*alp_b*(Acx*Dcz+Bcx*Ccz) ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcx*Bcy*Dcz*Ccy + Dcy*(Acx*Bcy*Dcz+Acy*Dcz*Bcx+Ccz*Bcx*Bcy)+2.d0*alp_a*alp_b*Bcx*Dcz ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcx*Bcy*Dcz*Dcy  )                                ! Vt(1,1,0,0,1,1)
                            else ! nb = 0
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_a*alp_b) ) * ( -Acx*Acy*Ccy*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acx*Acy*Dcy+Acy*Ccy*Bcx+Ccy*Acx*Bcy+2.d0*alp_a*alp_b*Acx) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcx*Bcy*Ccy+Bcy*Dcy*Acx+Dcy*Bcx*Acy &
               &+2.d0*alp_a*alp_b*Bcx) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcx*Bcy*Dcy )                            ! Vt(1,1,0,0,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                            elseif(nb==1) then 
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_a*alp_b) ) * ( -Acx*Acy*Ccz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acx*Acy*Dcz+Acy*Ccz*Bcx+Ccz*Acx*Bcy) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcx*Bcy*Ccz+Bcy*Dcz*Acx+Dcz*Bcx*Acy) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcx*Bcy*Dcz )                            !Vt(1,1,0,0,0,1)
                            else ! nb = 0
        rep = -1.d0/(4.d0*alp_a*alp_a) * ( Acx*Acy*ft(t)-(Acx*Bcy+Acy*Bcx)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Bcx*Bcy/(alp_a+alp_b)**2.d0 ) !Vt(1,1,0,0,0,0)
                            endif
                      endif
                endif

          endif
!
     else !ma==0

          if(na>1) then

        rep = 0.5d0*(na-1)/alp_a * ( repulsion(t,la,ma,na-2,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_b/(alp_a+alp_b) * (repulsion(t+1,la,ma,na-2,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_a/(alp_a+alp_b))*(nb/(na-1.d0)) *repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_a) * ( dcmplx(0.d0, -Acz )* repulsion(t,la,ma,na-1,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Bcz* repulsion(t+1,la,ma,na-1,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
         
          elseif(na==1) then

                if(lb>1)then
        rep = 0.5d0*(lb-1)/alp_b * ( repulsion(t,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(la/(lb-1.d0)) *repulsion(t+1,la-1,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccx )* repulsion(t,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcx* repulsion(t+1,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                elseif(lb==1) then

                      if(mb>1)then

        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                 &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
 
                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                            elseif(nb==1) then 
        rep = -dcmplx(0.d0, 1.d0/(32.d0*alp_a*alp_a*alp_b*alp_b*alp_b))*(-Acx*Acz*Ccz*Ccx*Ccy*ft(t) &
            &+ft(t+1)/(alp_a+alp_b) * (Acx*Acz*Ccz*Ccy*Dcx+Acx*Acz*Ccz*Ccx*Dcy+Ccx*Ccy*(Acx*Acz*Dcz+Acz*Ccz*Bcx+Ccz*Acx*Bcz)&
                   &+2.d0*alp_a*alp_b*Ccy*(Acx*Ccx+Acz*Ccz) ) &
            &-ft(t+2)/(alp_a+alp_b)**2.d0 * (Acx*Acz*Ccz*Dcx*Dcy +Ccx*Dcy*(Acx*Acz*Dcz+Acz*Ccz*Bcx+Ccz*Acx*Bcz)&
              &+Ccx*Ccy*(Acx*Bcz*Dcz+Acz*Dcz*Bcx+Ccz*Bcx*Bcz)+Dcx*Ccy*(Acx*Acz*Dcz+Acz*Ccz*Bcx+Ccz*Acx*Bcz ) &
                 &+2.d0*alp_a*alp_b*Dcy*(Acx*Ccx+Acz*Ccz) +Ccy*(2.d0*alp_a*alp_b*(Acx*Dcx+Acz*Dcz+Bcx*Ccx+Bcz*Ccz)+4.d0*(alp_a*alp_b)**2.d0 ) ) &
            &+ft(t+3)/(alp_a+alp_b)**3.d0 * (Ccx*Dcy*(Acx*Bcz*Dcz+Acz*Dcz*Bcx+Ccz*Bcx*Bcz)+Dcx*Dcy*(Acx*Acz*Dcz+Acz*Ccz*Bcx+Ccz*Acx*Bcz)&
              &+Ccx*Ccy*Bcx*Bcz*Dcz+ Dcx*Ccy*(Acx*Bcz*Dcz+Acz*Dcz*Bcx+Ccz*Bcx*Bcz) &
                 &+2.d0*alp_a*alp_b*Ccy*(Bcx*Dcx+Bcz*Dcz) +Dcy*(2.d0*alp_a*alp_b*(Acx*Dcx+Acz*Dcz+Bcx*Ccx+Bcz*Ccz)+4.d0*(alp_a*alp_b)**2.d0 )    ) &
            &-ft(t+4)/(alp_a+alp_b)**4.d0 * ( Ccx*Dcy*Bcx*Bcz*Dcz + Dcx*Dcy*(Acx*Bcz*Dcz+Acz*Dcz*Bcx+Ccz*Bcx*Bcz) + Dcx*Ccy*Bcx*Bcz*Dcz  &
                            &+2.d0*alp_a*alp_b*Dcy*(Bcx*Dcx+Bcz*Dcz))&
            &+ft(t+5)/(alp_a+alp_b)**5.d0 * Dcx*Dcy*Bcx*Bcz*Dcz )                              ! Vt(1,0,1,1,1,1)


                            else ! nn = 0
        rep = 1.d0/(16.d0*alp_a*alp_a*alp_b*alp_b) * ( Acx*Acz*Ccy*Ccx*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acx*Acz*Ccy*Dcx + Ccx*(Acx*Acz*Dcy+Acz*Ccy*Bcx+Ccy*Acx*Bcz) +2.d0*alp_a*alp_b*Acz*Ccy) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccx*(Bcx*Bcz*Ccy+Bcz*Dcy*Acx+Dcy*Bcx*Acz)+Dcx*(Acx*Acz*Dcy+Acz*Ccy*Bcx+Ccy*Acx*Bcz) &
              &+2.d0*alp_a*alp_b*(Acz*Dcy+ Ccy*Bcz) ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcx*Bcz*Dcy*Ccx + Dcx*(Acx*Bcz*Dcy+Acz*Dcy*Bcx+Ccy*Bcx*Bcz) +2.d0*alp_a*alp_b*Bcz*Dcy ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcx*Bcz*Dcy*Dcx  )                               ! Vt(1,0,1,1,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                            elseif(nb==1) then 
        rep = 1.d0/(16.d0*alp_a*alp_a*alp_b*alp_b) * ( Acx*Acz*Ccz*Ccx*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acx*Acz*Ccz*Dcx + Ccx*(Acx*Acz*Dcz+Acz*Ccz*Bcx+Ccz*Acx*Bcz) +2.d0*alp_a*alp_b*(Acx*Ccx +Acz*Ccz)) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccx*(Bcx*Bcz*Ccz+Bcz*Dcz*Acx+Dcz*Bcx*Acz)+Dcx*(Acx*Acz*Dcz+Acz*Ccz*Bcx+Ccz*Acx*Bcz)  &
                     &+ 2.d0*alp_a*alp_b*(Acx*Dcx +Acz*Dcz+ Bcx*Ccx+Bcz*Ccz) +(2.d0*alp_a*alp_b)**2.d0)  &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcx*Bcz*Dcz*Ccx + Dcx*(Acx*Bcz*Dcz+Acz*Dcz*Bcx+Ccz*Bcx*Bcz)+2.d0*alp_a*alp_b*(Bcx*Dcx +Bcz*Dcz) ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcx*Bcz*Dcz*Dcx  )                               !Vt(1,0,1,1,0,1)

                            else ! nb = 0
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_a*alp_b) ) * ( -Acx*Ccx*Acz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acx*Ccx*Bcz+Ccx*Acz*Bcx+Acz*Acx*Dcx+2.d0*alp_a*alp_b*Acz) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcx*Dcx*Acz+Dcx*Bcz*Acx+Bcz*Bcx*Ccx&
                  &+2.d0*alp_a*alp_b*Bcz) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcx*Dcx*Bcz )                               !Vt(1,0,1,1,0,0)
                            endif
                      endif
                else   !lb==0

                      if(mb>1)then

        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 

        rep = 1.d0/(16.d0*alp_a*alp_a*alp_b*alp_b) * ( Acx*Acz*Ccy*Ccz*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acx*Acz*Ccy*Dcz + Ccz*(Acx*Acz*Dcy+Acz*Ccy*Bcx+Ccy*Acx*Bcz) +2.d0*alp_a*alp_b*Acx*Ccy) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccz*(Bcx*Bcz*Ccy+Bcz*Dcy*Acx+Dcy*Bcx*Acz)+Dcz*(Acx*Acz*Dcy+Acz*Ccy*Bcx+Ccy*Acx*Bcz) &
              &+2.d0*alp_a*alp_b*(Acx*Dcy+ Ccy*Bcx) ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcx*Bcz*Dcy*Ccz + Dcz*(Acx*Bcz*Dcy+Acz*Dcy*Bcx+Ccy*Bcx*Bcz) + 2.d0*alp_a*alp_b*Bcx*Dcy ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcx*Bcz*Dcy*Dcz  )                             ! Vt(1,0,1,0,1,1)

                            else ! nb = 0
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_a*alp_b) ) * ( -Acx*Ccy*Acz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acx*Ccy*Bcz+Ccy*Acz*Bcx+Acz*Acx*Dcy) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcx*Dcy*Acz+Dcy*Bcz*Acx+Bcz*Bcx*Ccy) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcx*Dcy*Bcz )           ! Vt(1,0,1,0,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_a*alp_b) ) * ( -Acx*Ccz*Acz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acx*Ccz*Bcz+Ccz*Acz*Bcx+Acz*Acx*Dcz+2.d0*alp_a*alp_b*Acx) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcx*Dcz*Acz+Dcz*Bcz*Acx+Bcz*Bcx*Ccz&
                &+2.d0*alp_a*alp_b*Bcx) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcx*Dcz*Bcz )                             !Vt(1,0,1,0,0,1)
                            else ! nb = 0
        rep = -1.d0/(4.d0*alp_a*alp_a) * ( Acx*Acz*ft(t)-(Acx*Bcz+Acz*Bcx)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Bcx*Bcz/(alp_a+alp_b)**2.d0 ) !Vt(1,0,1,0,0,0)
                            endif
                      endif
                endif
          else !na==0

                if(lb>1)then

        rep = 0.5d0*(lb-1)/alp_b * ( repulsion(t,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                 &+( (alp_b/(alp_a+alp_b))*(la/(lb-1.d0)) *repulsion(t+1,la-1,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccx )* repulsion(t,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcx* repulsion(t+1,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
 
                elseif(lb==1) then

                      if(mb>1)then     
        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                            elseif(nb==1) then 
        rep = 1.d0/(16.d0*alp_a*alp_b*alp_b*alp_b) * ( Acx*Ccx*Ccz*Ccy*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acx*Ccx*Ccz*Dcy + Ccy*(Acx*Ccx*Dcz+Ccx*Ccz*Bcx+Ccz*Acx*Dcx) +2.d0*alp_a*alp_b*Ccz*Ccy) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccy*(Bcx*Dcx*Ccz+Dcx*Dcz*Acx+Dcz*Bcx*Ccx)+Dcy*(Acx*Ccx*Dcz+Ccx*Ccz*Bcx+Ccz*Acx*Dcx) &
               &+2.d0*alp_a*alp_b*(Dcz*Ccy+Dcy*Ccz) ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcx*Dcx*Dcz*Ccy + Dcy*(Acx*Dcx*Dcz+Ccx*Dcz*Bcx+Ccz*Bcx*Dcx) +2.d0*alp_a*alp_b*Dcz*Dcy ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcx*Dcx*Dcz*Dcy  )                                ! Vt(1,0,0,1,1,1)
                            else ! nb = 0
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_b*alp_b) ) * ( -Acx*Ccx*Ccy*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acx*Ccx*Dcy+Ccx*Ccy*Bcx+Ccy*Acx*Dcx+2.d0*alp_a*alp_b*Ccy) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcx*Dcx*Ccy+Dcx*Dcy*Acx+Dcy*Bcx*Ccx&
                 +2.d0*alp_a*alp_b*Dcy) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcx*Dcx*Dcy )      ! Vt(1,0,0,1,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                            elseif(nb==1) then 
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_b*alp_b) ) * ( -Acx*Ccx*Ccz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acx*Ccx*Dcz+Ccx*Ccz*Bcx+Ccz*Acx*Dcx+2.d0*alp_a*alp_b*Ccz) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcx*Dcx*Ccz+Dcx*Dcz*Acx+Dcz*Bcx*Ccx &
                &+2.d0*alp_a*alp_b*Dcz) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcx*Dcx*Dcz )         !Vt(1,0,0,1,0,1)
                            else ! nb = 0
        rep = -1.d0/(4.d0*alp_a*alp_b) * ( Acx*Ccx*ft(t)-(Acx*Dcx+Ccx*Bcx+2.d0*alp_a*alp_b)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Bcx*Dcx/(alp_a+alp_b)**2.d0 )
                                                                   !Vt(1,0,0,1,0,0)
                            endif
                      endif

                else   !lb==0

                      if(mb>1)then

        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                            elseif(nb==1) then 
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_b*alp_b) ) * ( -Acx*Ccy*Ccz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acx*Ccy*Dcz+Ccy*Ccz*Bcx+Ccz*Acx*Dcy) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcx*Dcy*Ccz+Dcy*Dcz*Acx+Dcz*Bcx*Ccy) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcx*Dcy*Dcz )                                ! Vt(1,0,0,0,1,1)
                            else ! nb = 0
        rep = -1.d0/(4.d0*alp_a*alp_b) * ( Acx*Ccy*ft(t)-(Acx*Dcy+Ccy*Bcx)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Bcx*Dcy/(alp_a+alp_b)**2.d0 ) ! Vt(1,0,0,0,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = -1.d0/(4.d0*alp_a*alp_b) * ( Acx*Ccz*ft(t)-(Acx*Dcz+Ccz*Bcx)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Bcx*Dcz/(alp_a+alp_b)**2.d0 ) !Vt(1,0,0,0,0,1)
                            else ! nb = 0
        rep= -dcmplx( 0.d0, 1.d0/(2.d0*alp_a) ) * (-Acx*ft(t) + Bcx/(alp_a+alp_b) * ft(t+1) )   !Vt(1,0,0,0,0,0)
                            endif
                      endif
                endif
          endif
     endif

else !la==0

     if(ma>1) then

        rep = 0.5d0*(ma-1)/alp_a * ( repulsion(t,la,ma-2,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_b/(alp_a+alp_b) * (repulsion(t+1,la,ma-2,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_a/(alp_a+alp_b))*(mb/(ma-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_a) * ( dcmplx(0.d0, -Acy )* repulsion(t,la,ma-1,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Bcy* repulsion(t+1,la,ma-1,na,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

     elseif(ma==1) then

          if(na>1) then

        rep = 0.5d0*(na-1)/alp_a * ( repulsion(t,la,ma,na-2,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
             &alp_b/(alp_a+alp_b) * (repulsion(t+1,la,ma,na-2,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_a/(alp_a+alp_b))*(nb/(na-1.d0)) *repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
             &-1.d0/(2.d0*alp_a) * ( dcmplx(0.d0, -Acz )* repulsion(t,la,ma,na-1,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
             &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Bcz* repulsion(t+1,la,ma,na-1,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

          elseif(na==1) then

                if(lb>1)then

        rep = 0.5d0*(lb-1)/alp_b * ( repulsion(t,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(la/(lb-1.d0)) *repulsion(t+1,la-1,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccx )* repulsion(t,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcx* repulsion(t+1,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                elseif(lb==1) then

                      if(mb>1)then
 
        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                      elseif(mb==1) then 
                                    
                            if(nb>1)then
        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = -dcmplx(0.d0, 1.d0/(32.d0*alp_a*alp_a*alp_b*alp_b*alp_b))*(-Acy*Acz*Ccz*Ccx*Ccy*ft(t) &
            &+ft(t+1)/(alp_a+alp_b) * (Acy*Acz*Ccz*Ccy*Dcx+Acy*Acz*Ccz*Ccx*Dcy+Ccx*Ccy*(Acy*Acz*Dcz+Acz*Ccz*Bcy+Ccz*Acy*Bcz) &
               &+2.d0*alp_a*alp_b*Ccx*(Acz*Ccz+Acy*Ccy) ) &
            &-ft(t+2)/(alp_a+alp_b)**2.d0 * (Acy*Acz*Ccz*Dcx*Dcy +Ccx*Dcy*(Acy*Acz*Dcz+Acz*Ccz*Bcy+Ccz*Acy*Bcz)&
              &+Ccx*Ccy*(Acy*Bcz*Dcz+Acz*Dcz*Bcy+Ccz*Bcy*Bcz)+Dcx*Ccy*(Acy*Acz*Dcz+Acz*Ccz*Bcy+Ccz*Acy*Bcz)  &
                &+2.d0*alp_a*alp_b*Dcx*(Acz*Ccz+Acy*Ccy) +Ccx*(2.d0*alp_a*alp_b*(Acz*Dcz+Acy*Dcy+Bcz*Ccz+Bcy*Ccy)+4.d0*(alp_a*alp_b)**2.d0 )) &
            &+ft(t+3)/(alp_a+alp_b)**3.d0 * (Ccx*Dcy*(Acy*Bcz*Dcz+Acz*Dcz*Bcy+Ccz*Bcy*Bcz)+Dcx*Dcy*(Acy*Acz*Dcz+Acz*Ccz*Bcy+Ccz*Acy*Bcz)&
              &+Ccx*Ccy*Bcy*Bcz*Dcz+ Dcx*Ccy*(Acy*Bcz*Dcz+Acz*Dcz*Bcy+Ccz*Bcy*Bcz)  &
              &+2.d0*alp_a*alp_b*Ccx*(Bcz*Dcz+Bcy*Dcy) +Dcx*(2.d0*alp_a*alp_b*(Acz*Dcz+Acy*Dcy+Bcz*Ccz+Bcy*Ccy)+4.d0*(alp_a*alp_b)**2.d0 )   ) &
            &-ft(t+4)/(alp_a+alp_b)**4.d0 * ( Ccx*Dcy*Bcy*Bcz*Dcz + Dcx*Dcy*(Acy*Bcz*Dcz+Acz*Dcz*Bcy+Ccz*Bcy*Bcz) + Dcx*Ccy*Bcy*Bcz*Dcz  &
                             &+2.d0*alp_a*alp_b*Dcx*(Bcz*Dcz+Bcy*Dcy) )&
            &+ft(t+5)/(alp_a+alp_b)**5.d0 * Dcx*Dcy*Bcy*Bcz*Dcz )                              ! Vt(0,1,1,1,1,1)
                            else ! nb = 0
        rep = 1.d0/(16.d0*alp_a*alp_a*alp_b*alp_b) * ( Ccy*Acy*Acz*Ccx*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Ccy*Acy*Acz*Dcx + Ccx*(Ccy*Acy*Bcz+Acy*Acz*Dcy+Acz*Ccy*Bcy) +2.d0*alp_a*alp_b*Ccx*Acz ) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccx*(Dcy*Bcy*Acz+Bcy*Bcz*Ccy+Bcz*Dcy*Acy)+Dcx*(Ccy*Acy*Bcz+Acy*Acz*Dcy+Acz*Ccy*Bcy) &
                 &+ 2.d0*alp_a*alp_b*(Dcx*Acz+Bcz*Ccx) ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Dcy*Bcy*Bcz*Ccx + Dcx*(Ccy*Bcy*Bcz+Acy*Bcz*Dcy+Acz*Dcy*Bcy)+ 2.d0*alp_a*alp_b*Dcx*Bcz ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Dcy*Bcy*Bcz*Dcx  )                               ! Vt(0,1,1,1,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then
        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                            elseif(nb==1) then 
        rep = 1.d0/(16.d0*alp_a*alp_a*alp_b*alp_b) * ( Ccz*Acy*Acz*Ccx*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Ccz*Acy*Acz*Dcx + Ccx*(Ccz*Acy*Bcz+Acy*Acz*Dcz+Acz*Ccz*Bcy)+  2.d0*alp_a*alp_b*Acy*Ccx ) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccx*(Dcz*Bcy*Acz+Bcy*Bcz*Ccz+Bcz*Dcz*Acy)+Dcx*(Ccz*Acy*Bcz+Acy*Acz*Dcz+Acz*Ccz*Bcy) &
              &+ 2.d0*alp_a*alp_b*(Acy*Dcx+Ccx*Bcy) ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Dcz*Bcy*Bcz*Ccx + Dcx*(Ccz*Bcy*Bcz+Acy*Bcz*Dcz+Acz*Dcz*Bcy) + 2.d0*alp_a*alp_b*Bcy*Dcx ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Dcz*Bcy*Bcz*Dcx  )                              !Vt(0,1,1,1,0,1)
                            else ! nb = 0
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_a*alp_b) ) * ( -Ccx*Acy*Acz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Ccx*Acy*Bcz+Acy*Acz*Dcx+Acz*Ccx*Bcy) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Dcx*Bcy*Acz+Bcy*Bcz*Ccx+Bcz*Dcx*Acy) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Dcx*Bcy*Bcz )                       !Vt(0,1,1,1,0,0)
                            endif
                      endif
                else   !lb==0

                      if(mb>1)then

        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                      elseif(mb==1) then 
                                    
                            if(nb>1)then


        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0)) *repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = 1.d0/(16.d0*alp_a*alp_a*alp_b*alp_b) * ( Ccy*Acy*Acz*Ccz*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Ccy*Acy*Acz*Dcz + Ccz*(Ccy*Acy*Bcz+Acy*Acz*Dcy+Acz*Ccy*Bcy) +2.d0*alp_a*alp_b*(Acz*Ccz +Acy*Ccy)) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccz*(Dcy*Bcy*Acz+Bcy*Bcz*Ccy+Bcz*Dcy*Acy)+Dcz*(Ccy*Acy*Bcz+Acy*Acz*Dcy+Acz*Ccy*Bcy) &
               &+ 2.d0*alp_a*alp_b*(Acz*Dcz +Acy*Dcy+ Bcz*Ccz+Bcy*Ccy) +(2.d0*alp_a*alp_b)**2.d0  ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Dcy*Bcy*Bcz*Ccz + Dcz*(Ccy*Bcy*Bcz+Acy*Bcz*Dcy+Acz*Dcy*Bcy) +2.d0*alp_a*alp_b*(Bcz*Dcz +Bcy*Dcy) ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Dcy*Bcy*Bcz*Dcz  )                                ! Vt(0,1,1,0,1,1)

                            else ! nb = 0
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_a*alp_b) ) * ( -Ccy*Acy*Acz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Ccy*Acy*Bcz+Acy*Acz*Dcy+Acz*Ccy*Bcy+2.d0*alp_a*alp_b*Acz) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Dcy*Bcy*Acz+Bcy*Bcz*Ccy+Bcz*Dcy*Acy &
                 &+2.d0*alp_a*alp_b*Bcz) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Dcy*Bcy*Bcz )                               ! Vt(0,1,1,0,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
             &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
             &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
             &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_a*alp_b) ) * ( -Ccz*Acy*Acz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Ccz*Acy*Bcz+Acy*Acz*Dcz+Acz*Ccz*Bcy+2.d0*alp_a*alp_b*Acy) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Dcz*Bcy*Acz+Bcy*Bcz*Ccz+Bcz*Dcz*Acy&
                 &+2.d0*alp_a*alp_b*Bcy) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Dcz*Bcy*Bcz )                               !Vt(0,1,1,0,0,1)
                            else ! nn = 0
        rep = -1.d0/(4.d0*alp_a*alp_a) * ( Acy*Acz*ft(t)-(Acy*Bcz+Acz*Bcy)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Bcy*Bcz/(alp_a+alp_b)**2.d0 ) !Vt(0,1,1,0,0,0)
                            endif
                      endif
                endif
          else !na==0

                if(lb>1)then
        rep = 0.5d0*(lb-1)/alp_b * ( repulsion(t,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
             &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(la/(lb-1.d0)) *repulsion(t+1,la-1,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
             &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccx )* repulsion(t,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
             &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcx* repulsion(t+1,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )


             
                elseif(lb==1) then

                      if(mb>1)then

        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                            elseif(nb==1) then 
        rep = 1.d0/(16.d0*alp_a*alp_b*alp_b*alp_b) * ( Acy*Ccx*Ccz*Ccy*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acy*Ccx*Ccz*Dcy + Ccy*(Acy*Ccx*Dcz+Ccx*Ccz*Bcy+Ccz*Acy*Dcx) +2.d0*alp_a*alp_b*Ccx*Ccz ) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccy*(Bcy*Dcx*Ccz+Dcx*Dcz*Acy+Dcz*Bcy*Ccx)+Dcy*(Acy*Ccx*Dcz+Ccx*Ccz*Bcy+Ccz*Acy*Dcx)  &
               &+2.d0*alp_a*alp_b*(Ccx*Dcz+Dcx*Ccz) ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcy*Dcx*Dcz*Ccy + Dcy*(Acy*Dcx*Dcz+Ccx*Dcz*Bcy+Ccz*Bcy*Dcx)+ 2.d0*alp_a*alp_b*Dcx*Dcz  ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcy*Dcx*Dcz*Dcy  )                              ! Vt(0,1,0,1,1,1)

                            else ! nn = 0
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_b*alp_b) ) * ( -Acy*Ccx*Ccy*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acy*Ccx*Dcy+Ccx*Ccy*Bcy+Ccy*Acy*Dcx+2.d0*alp_a*alp_b*Ccx) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcy*Dcx*Ccy+Dcx*Dcy*Acy+Dcy*Bcy*Ccx &
               &+2.d0*alp_a*alp_b*Dcx) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcy*Dcx*Dcy )         ! Vt(0,1,0,1,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                            elseif(nb==1) then 
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_b*alp_b) ) * ( -Acy*Ccz*Ccx*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acy*Ccz*Dcx+Ccz*Ccx*Bcy+Ccx*Acy*Dcz) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcy*Dcz*Ccx+Dcz*Dcx*Acy+Dcx*Bcy*Ccz) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcy*Dcz*Dcx )                              !Vt(0,1,0,1,0,1)
                            else ! nn = 0
        rep = -1.d0/(4.d0*alp_a*alp_b) * ( Acy*Ccx*ft(t)-(Acy*Dcx+Ccx*Bcy)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Bcy*Dcx/(alp_a+alp_b)**2.d0 ) !Vt(0,1,0,1,0,0)
                            endif
                       endif
                else   !lb==0

                      if(mb>1)then
        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_b*alp_b) ) * ( -Acy*Ccz*Ccy*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acy*Ccz*Dcy+Ccz*Ccy*Bcy+Ccy*Acy*Dcz+2.d0*alp_a*alp_b*Ccz) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcy*Dcz*Ccy+Dcz*Dcy*Acy+Dcy*Bcy*Ccz &
                 &+2.d0*alp_a*alp_b*Dcz) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcy*Dcz*Dcy )                                ! Vt(0,1,0,0,1,1)
                            else ! nb = 0
        rep = -1.d0/(4.d0*alp_a*alp_b) * ( Acy*Ccy*ft(t)-(Acy*Dcy+Ccy*Bcy+2.d0*alp_a*alp_b)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Bcy*Dcy/(alp_a+alp_b)**2.d0 )
                                                                                          ! Vt(0,1,0,0,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                            elseif(nb==1) then 
        rep = -1.d0/(4.d0*alp_a*alp_b) * ( Acy*Ccz*ft(t)-(Acy*Dcz+Ccz*Bcy)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Bcy*Dcz/(alp_a+alp_b)**2.d0 ) !Vt(0,1,0,0,0,1)
                            else ! nb = 0
        rep= -dcmplx( 0.d0, 1.d0/(2.d0*alp_a) ) * (-Acy*ft(t) + Bcy/(alp_a+alp_b) * ft(t+1) )  !Vt(0,1,0,0,0,0)
                            endif
                      endif
                endif
          endif
     else !ma==0

          if(na>1) then

        rep = 0.5d0*(na-1)/alp_a * ( repulsion(t,la,ma,na-2,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_b/(alp_a+alp_b) * (repulsion(t+1,la,ma,na-2,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_a/(alp_a+alp_b))*(nb/(na-1.d0)) *repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_a) * ( dcmplx(0.d0, -Acz )* repulsion(t,la,ma,na-1,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Bcz* repulsion(t+1,la,ma,na-1,lb,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
         
          elseif(na==1) then

                if(lb>1)then

         rep = 0.5d0*(lb-1)/alp_b * ( repulsion(t,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(la/(lb-1.d0)) *repulsion(t+1,la-1,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccx )* repulsion(t,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcx* repulsion(t+1,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )


             
                elseif(lb==1) then

                      if(mb>1)then

        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                            elseif(nb==1) then 
        rep = 1.d0/(16.d0*alp_a*alp_b*alp_b*alp_b) * ( Acz*Ccx*Ccz*Ccy*ft(t) &
             &-ft(t+1)/(alp_a+alp_b) * (Acz*Ccx*Ccz*Dcy + Ccy*(Acz*Ccx*Dcz+Ccx*Ccz*Bcz+Ccz*Acz*Dcx) + 2.d0*alp_a*alp_b*Ccx*Ccy ) &
             &+ft(t+2)/(alp_a+alp_b)**2.d0 * (Ccy*(Bcz*Dcx*Ccz+Dcx*Dcz*Acz+Dcz*Bcz*Ccx)+Dcy*(Acz*Ccx*Dcz+Ccx*Ccz*Bcz+Ccz*Acz*Dcx) &
              &+2.d0*alp_a*alp_b*(Ccx*Dcy+Dcx*Ccy)  ) &
             &-ft(t+3)/(alp_a+alp_b)**3.d0 * (Bcz*Dcx*Dcz*Ccy + Dcy*(Acz*Dcx*Dcz+Ccx*Dcz*Bcz+Ccz*Bcz*Dcx) +2.d0*alp_a*alp_b*Dcx*Dcy  ) &
             &+ft(t+4)/(alp_a+alp_b)**4.d0 *  Bcz*Dcx*Dcz*Dcy  )                               ! Vt(0,0,1,1,1,1)
                            else ! nn = 0
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_b*alp_b) ) * ( -Acz*Ccx*Ccy*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acz*Ccx*Dcy+Ccx*Ccy*Bcz+Ccy*Acz*Dcx) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcz*Dcx*Ccy+Dcx*Dcy*Acz+Dcy*Bcz*Ccx) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcz*Dcx*Dcy )                              ! Vt(0,0,1,1,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                            elseif(nb==1) then 
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_b*alp_b) ) * ( -Acz*Ccx*Ccz*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acz*Ccx*Dcz+Ccx*Ccz*Bcz+Ccz*Acz*Dcx+2.d0*alp_a*alp_b*Ccx) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcz*Dcx*Ccz+Dcx*Dcz*Acz+Dcz*Bcz*Ccx &
               &+2.d0*alp_a*alp_b*Dcx) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcz*Dcx*Dcz )                                !Vt(0,0,1,1,0,1)
                            else ! nb = 0
        rep = -1.d0/(4.d0*alp_a*alp_b) * ( Ccx*Acz*ft(t)-(Ccx*Bcz+Acz*Dcx)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Dcx*Bcz/(alp_a+alp_b)**2.d0 ) !Vt(0,0,1,1,0,0)
                            endif
                      endif
                else   !lb==0

                      if(mb>1)then

        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

             
                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                            elseif(nb==1) then 
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_a*alp_b*alp_b) ) * ( -Acz*Ccz*Ccy*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Acz*Ccz*Dcy+Ccz*Ccy*Bcz+Ccy*Acz*Dcz+2.d0*alp_a*alp_b*Ccy) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Bcz*Dcz*Ccy+Dcz*Dcy*Acz+Dcy*Bcz*Ccz &
              &+2.d0*alp_a*alp_b*Dcy) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Bcz*Dcz*Dcy )                                 ! Vt(0,0,1,0,1,1)
                            else ! nb = 0
        rep = -1.d0/(4.d0*alp_a*alp_b) * ( Ccy*Acz*ft(t)-(Ccy*Bcz+Acz*Dcy)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Dcy*Bcz/(alp_a+alp_b)**2.d0 ) ! Vt(0,0,1,0,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = -1.d0/(4.d0*alp_a*alp_b) * ( Ccz*Acz*ft(t)-(Ccz*Bcz+Acz*Dcz+2.d0*alp_a*alp_b)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Dcz*Bcz/(alp_a+alp_b)**2.d0 )
                                                                                                 !Vt(0,0,1,0,0,1)
                                  else ! nb = 0
        rep= -dcmplx( 0.d0, 1.d0/(2.d0*alp_a) ) * (-Acz*ft(t) + Bcz/(alp_a+alp_b) * ft(t+1) )  !Vt(0,0,1,0,0,0)
                            endif
                      endif
                endif
          else !na==0

                if(lb>1)then

         rep = 0.5d0*(lb-1)/alp_b * ( repulsion(t,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb-2,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(la/(lb-1.d0)) *repulsion(t+1,la-1,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccx )* repulsion(t,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcx* repulsion(t+1,la,ma,na,lb-1,mb,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                elseif(lb==1) then

                      if(mb>1)then
        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep = 0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                            elseif(nb==1) then 
        rep = dcmplx(0.d0, 1.d0/(8.d0*alp_b*alp_b*alp_b) ) * ( -Ccx*Ccz*Ccy*ft(t) +ft(t+1)/(alp_a+alp_b)&
             &*(Ccx*Ccz*Dcy+Ccz*Ccy*Dcx+Ccy*Ccx*Dcz) - ft(t+2)/(alp_a+alp_b)**2.d0 *(Dcx*Dcz*Ccy+Dcz*Dcy*Ccx+Dcy*Dcx*Ccz) &
             &+ ft(t+3)/(alp_a+alp_b)**3.d0 *Dcx*Dcz*Dcy )                                ! Vt(0,0,0,1,1,1)

                            else ! nb = 0
        rep = -1.d0/(4.d0*alp_b*alp_b) * ( Ccx*Ccy*ft(t)-(Ccx*Dcy+Ccy*Dcx)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Dcx*Dcy/(alp_a+alp_b)**2.d0 ) ! Vt(0,0,0,1,1,0)
                            endif
                      else !mb==0
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )


             
                            elseif(nb==1) then 
        rep = -1.d0/(4.d0*alp_b*alp_b) * ( Ccx*Ccz*ft(t)-(Ccx*Dcz+Ccz*Dcx)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Dcx*Dcz/(alp_a+alp_b)**2.d0 )  !Vt(0,0,0,1,0,1)
                            else ! nb = 0
        rep= -dcmplx( 0.d0, 1.d0/(2.d0*alp_b) ) * (-Ccx*ft(t) + Dcx/(alp_a+alp_b) * ft(t+1) )     !Vt(0,0,0,1,0,0)
                            endif
                      endif

                else   !lb==0

                      if(mb>1)then

        rep = 0.5d0*(mb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
            &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb-2,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(ma/(mb-1.d0)) *repulsion(t+1,la,ma-1,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
            &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccy )* repulsion(t,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
            &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcy* repulsion(t+1,la,ma,na,lb,mb-1,nb,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
             
                      elseif(mb==1) then 
                                    
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )

                            elseif(nb==1) then 
        rep = -1.d0/(4.d0*alp_b*alp_b) * ( Ccy*Ccz*ft(t)-(Ccy*Dcz+Ccz*Dcy)/(alp_a+alp_b)*ft(t+1)+ft(t+2)*Dcy*Dcz/(alp_a+alp_b)**2.d0 ) ! Vt(0,0,0,0,1,1)

                            else ! nb = 0
        rep= -dcmplx( 0.d0, 1.d0/(2.d0*alp_b) ) * (-Ccy*ft(t) + Dcy/(alp_a+alp_b) * ft(t+1) )    ! Vt(0,0,0,0,1,0)


                            endif
                      else !mb==0
                            if(nb>1)then

        rep =  0.5d0*(nb-1)/alp_b * ( repulsion(t,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) - &
           &alp_a/(alp_a+alp_b) * (repulsion(t+1,la,ma,na,lb,mb,nb-2,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft))&
                &+( (alp_b/(alp_a+alp_b))*(na/(nb-1.d0))*repulsion(t+1,la,ma,na-1,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)) )&
           &-1.d0/(2.d0*alp_b) * ( dcmplx(0.d0, -Ccz )* repulsion(t,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft)&
           &+1.d0/(alp_a+alp_b)*dcmplx(0.d0,1.d0)*Dcz* repulsion(t+1,la,ma,na,lb,mb,nb-1,alp_a,ax,ay,az,alp_b,bx,by,bz,vx,vy,vz,cx,cy,cz,tmax,ft) )
      
                            elseif(nb==1) then 

        rep= -dcmplx( 0.d0, 1.d0/(2.d0*alp_b) ) * (-Ccz*ft(t) + Dcz/(alp_a+alp_b) * ft(t+1) )     !Vt(0,0,0,0,0,1)

                            else ! nb = 0
                            rep = ft(t)      !Vt(0,0,0,0,0,0)
                            endif
                      endif
                endif
          endif
     endif

endif !la



         repul(la,ma,na,lb,mb,nb,t) = rep
         done(la,ma,na,lb,mb,nb,t)  = 0

!17 continue

end function repulsion

end module repulsion_mod
