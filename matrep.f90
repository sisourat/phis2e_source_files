!--------------22.Feb.2017--by J.W.GAO.--follow-Nico's-one-electron-code--------------
module MatRep
use newtypes
use repulsion_mod
use fdn_mod
use cgto, only : normgto, genPolgto
use centerlib
use matrices
use omp_lib
implicit none
 

  contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine StaMatRep(geom)

  integer, dimension(*) :: geom

  integer :: i, j, k, l, iy, jy, ky, ly, li, lj, lk, ll, nli, nlj, nlk, nll !, k2ij,k2kl,k2cij,k2ckl
  integer :: igto, jgto, kgto, lgto, icgto, jcgto, kcgto, lcgto, im, jm, km, lm, kicgto, kjcgto, kkcgto, klcgto
  integer :: lmax,lmax1,lmax2, tmax, nsizeGto
  integer :: itrans, jtrans, ktrans, ltrans, mi, mj, mk, ml, kirep,  kjrep,  kkrep,  klrep



 ! integer :: kcirep1,  kcjrep1,  kckrep1,  kclrep1, kcirep2,  kcjrep2,  kckrep2,  kclrep2
 ! integer :: ich,jch,kch,lch



  double precision :: xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl, xpik, ypik, zpik, xpjl, ypjl, zpjl
  double precision :: alpi, alpj, alpk, alpl, Gamik, Gamjl, rabik2, rabjl2, kfactorik, kfactorjl 
  double precision :: cx2c, cy2c, cz2c

  double precision :: vx, vy, vz

  double complex, dimension(:,:,:,:,:,:), save, allocatable :: repmat

  double complex, dimension(:,:,:,:), save, allocatable ::  repGto
  double complex, dimension(:,:,:,:), allocatable ::  repGto1
  double precision, dimension(:,:,:), save, allocatable ::  fkx, fky, fkz, hkx, hky, hkz

  double precision, dimension(:,:,:), save, allocatable :: normGtoI, normGtoJ, normGtoK, normGtoL

  integer :: la, ma, na, lb, nb, mb

  !$OMP THREADPRIVATE( repmat,repGto, fkx, fky, fkz, hkx, hky, hkz,normGtoI, normGtoJ, normGtoK, normGtoL)



   nsizeGto = sum(Centers(:)%nsizegtoBlocks)

   allocate( repGto1(nsizeGto,nsizeGto,nsizeGto,nsizeGto )  ) 

!   repGto1=0.d0
               

  do i = 1, nCenters
   do iy = 1, nYlmMax


    if(Centers(i)%block(iy)%ncgto==0) cycle

    li = iy - 1
    nli = (li+1)*(li+2)/2

     do j = 1, nCenters
      do jy = 1, nYlmMax

        if(Centers(j)%block(jy)%ncgto==0) cycle

        lj = jy - 1
        nlj = (lj+1)*(lj+2)/2

       
        do k = 1, nCenters
         do ky = 1, nYlmMax

          if(Centers(k)%block(ky)%ncgto==0) cycle

          lk = ky - 1
          nlk = (lk+1)*(lk+2)/2

           do l = 1, nCenters
            do ly = 1, nYlmMax

             if(Centers(l)%block(ly)%ncgto==0) cycle

             ll = ly - 1
             nll = (ll+1)*(ll+2)/2


!**********************************************************************************************************
!**********************************************************************************************************

!$OMP PARALLEL PRIVATE(igto,jgto,kgto,lgto,alpi,alpj,alpk,alpl,Gamik,Gamjl,kfactorik,kfactorjl,xpik,xpjl,ypik,ypjl,zpik,zpjl,cx2c,cy2c,cz2c,kirep,kjrep,kkrep,klrep,mi,mj,mk,ml) 
!$OMP DO COLLAPSE(4) SCHEDULE(DYNAMIC)



                do igto = 1, Centers(i)%block(iy)%ngto
                 do jgto = 1, Centers(j)%block(jy)%ngto

                  do kgto =1, Centers(k)%block(ky)%ngto
                   do lgto = 1, Centers(l)%block(ly)%ngto


                      xi = Centers(i)%x(geom(i)); yi = Centers(i)%y(geom(i)); zi = Centers(i)%z(geom(i),1)
                      xj = Centers(j)%x(geom(j)); yj = Centers(j)%y(geom(j)); zj = Centers(j)%z(geom(j),1)

                      xk = Centers(k)%x(geom(k)); yk = Centers(k)%y(geom(k)); zk = Centers(k)%z(geom(k),1)
                      xl = Centers(l)%x(geom(l)); yl = Centers(l)%y(geom(l)); zl = Centers(l)%z(geom(l),1)


                      alpi = Centers(i)%block(iy)%gto(igto); alpj = Centers(j)%block(jy)%gto(jgto)
                      alpk = Centers(k)%block(ky)%gto(kgto); alpl = Centers(l)%block(ly)%gto(lgto)

                      Gamik =  alpi + alpk
                      Gamjl =  alpj + alpl

                      rabik2 = (xi-xk)**2 + (yi-yk)**2 + (zi-zk)**2
                      rabjl2 = (xj-xl)**2 + (yj-yl)**2 + (zj-zl)**2

                      kfactorik = exp(-alpi*alpk*rabik2/Gamik)
                      kfactorjl = exp(-alpj*alpl*rabjl2/Gamjl)
            
                      xpik = (alpi*xi + alpk*xk)/Gamik ; ypik = (alpi*yi + alpk*yk)/Gamik ; zpik = (alpi*zi + alpk*zk)/Gamik
                      xpjl = (alpj*xj + alpl*xl)/Gamjl ; ypjl = (alpj*yj + alpl*yl)/Gamjl ; zpjl = (alpj*zj + alpl*zl)/Gamjl

                      cx2c = xpik - xpjl ; cy2c = ypik - ypjl ; cz2c = zpik - zpjl  ! the distance between the centers of two gaussians
!nico China 2019 change ci in -ci to be consistent with fdn
!      cx2c = -cx2c; cy2c = - cy2c; cz2c = - cz2c  

 
                      lmax1 = (li+lk) ; lmax2=(lj+ll) 

                      lmax = max(lmax1 , lmax2)   ; tmax = 6*lmax ! How to choose the Lmax and Tmax??
  

!-----repmat-repulsion-integrals------------------------------------------------------------------------------


                      allocate( repmat(0:lmax,0:lmax,0:lmax,0:lmax,0:lmax,0:lmax) )
                      allocate(normGtoI(0:lmax,0:lmax,0:lmax),normGtoJ(0:lmax,0:lmax,0:lmax))
                      allocate(normGtoK(0:lmax,0:lmax,0:lmax),normGtoL(0:lmax,0:lmax,0:lmax))

                      call normgto(lmax,alpi,normGtoI)
                      call normgto(lmax,alpj,normGtoJ)
                      call normgto(lmax,alpk,normGtoK)
                      call normgto(lmax,alpl,normGtoL)



!---------------------repulsion-integrals, the most time-consuming part of this subroutine---here, set a=0.d0, b=0.d0, v=0.d0-------------------


                      call repulsion_driver(lmax,0,Gamik,0.d0,0.d0,0.d0,Gamjl,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,cx2c,cy2c,cz2c,tmax,repmat)


              
                                                                                                       
                      repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) !&
!                           & * EXP( dcmplx(0d0,2.d0*0.0*xpik) ) !---exp(2ivRs) comes from the shift of ETF, no need for now ...--
                                         
                     
! fkx, fky, fkz and hkx, hky, hkz contain the coeff. from Gaussian product rule---------------

                      allocate(fkx(0:lmax,0:lmax,0:lmax),fky(0:lmax,0:lmax,0:lmax),fkz(0:lmax,0:lmax,0:lmax))
                      allocate(hkx(0:lmax,0:lmax,0:lmax),hky(0:lmax,0:lmax,0:lmax),hkz(0:lmax,0:lmax,0:lmax))
          

                      call fdn(xi-xpik,xk-xpik,lmax,fkx) 
                      call fdn(yi-ypik,yk-ypik,lmax,fky) 
                      call fdn(zi-zpik,zk-zpik,lmax,fkz) 
                      call fdn(xj-xpjl,xl-xpjl,lmax,hkx) 
                      call fdn(yj-ypjl,yl-ypjl,lmax,hky) 
                      call fdn(zj-zpjl,zl-zpjl,lmax,hkz) 

! now the matrices are computed according to the Gaussian product rule for all
! (nli+1)*(nli+2)/2 and (nlj+1)*(nlj+2)/2 cartesian gaussian functions

                      allocate(  repGto(nll,nlk,nlj,nli)  ) 

                      call sumfkhk(lmax,li,lj,lk,ll,nli,nlj,nlk,nll,fkx,fky,fkz,hkx,hky,hkz,normGtoI,normGtoJ,normGtoK,normGtoL,repmat,repGto)

                    !  write(*,*) igto,jgto,kgto,lgto
                    !  write(*,*) real(repGto)
                    !  read(*,*)

!----------save repulsion integrals-------------------------------------------------------------------

                      kirep=0;kjrep=0;kkrep=0;klrep=0

                      do mi=1,nli
                        do mj=1,nlj
                          do mk=1,nlk
                            do ml=1,nll

                                kirep = sum(Centers(1:i-1)%nsizegtoBlocks) + sum(Centers(i)%nsizegtoPerBlock(1:iy-1)) + (igto-1)*(nli) + mi
                                kjrep = sum(Centers(1:j-1)%nsizegtoBlocks) + sum(Centers(j)%nsizegtoPerBlock(1:jy-1)) + (jgto-1)*(nlj) + mj
                                kkrep = sum(Centers(1:k-1)%nsizegtoBlocks) + sum(Centers(k)%nsizegtoPerBlock(1:ky-1)) + (kgto-1)*(nlk) + mk
                                klrep = sum(Centers(1:l-1)%nsizegtoBlocks) + sum(Centers(l)%nsizegtoPerBlock(1:ly-1)) + (lgto-1)*(nll) + ml      


                                repGto1 (klrep,kkrep,kjrep,kirep)  =  repGto(ml,mk,mj,mi)  
!-----------------------------------------interchange basis on electron 1------------------------------------------------------------------  
                      
!                               if (  repGto1 (klrep,kirep,kjrep,kkrep) .eq. 0.0)then

                            !    repGto1 (klrep,kirep,kjrep,kkrep)  =  repGto(ml,mk,mj,mi)  
!                               endif
!-----------------------------------------interchange basis on electron 2------------------------------------------------------------------  

!                               if ( repGto1 (kjrep,kkrep,klrep,kirep).eq. 0.0)then
                            !    repGto1 (kjrep,kkrep,klrep,kirep)  =  repGto(ml,mk,mj,mi)  
!                               endif

!-----------------------------------------interchange basis on both electron 1 & 2------------------------------------------------------------------  

!                               if ( repGto1 (kjrep,kirep,klrep,kkrep) .eq. 0.0)then
                            !    repGto1 (kjrep,kirep,klrep,kkrep)  =  repGto(ml,mk,mj,mi)  

!                               endif

                              

                            enddo         
                          enddo   
                        enddo   
                      enddo  

                      deallocate(repmat,repGto,fkx,fky,fkz,hkx,hky,hkz,normGtoI,normGtoJ,normGtoK,normGtoL)

                   enddo  !lgto
                  enddo   !kgto
                 enddo    !jgto
                enddo     !igto

!$OMP END DO 
!$OMP END PARALLEL       
 



            enddo   !ly
           enddo    !l
 
         enddo      !ky
        enddo       !k
 
      enddo         !jy
     enddo          !j

   enddo            !iy
  enddo             !i




!-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-



 do i = 1, nCenters
   do iy = 1, nYlmMax

    if(Centers(i)%block(iy)%ncgto==0) cycle

    li = iy - 1
    nli = (li+1)*(li+2)/2

     do j = 1, nCenters
      do jy = 1, nYlmMax

        if(Centers(j)%block(jy)%ncgto==0) cycle

        lj = jy - 1
        nlj = (lj+1)*(lj+2)/2

       
        do k = 1, nCenters
         do ky = 1, nYlmMax

          if(Centers(k)%block(ky)%ncgto==0) cycle

          lk = ky - 1
          nlk = (lk+1)*(lk+2)/2

           do l = 1, nCenters
            do ly = 1, nYlmMax

             if(Centers(l)%block(ly)%ncgto==0) cycle

             ll = ly - 1
             nll = (ll+1)*(ll+2)/2


!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
               
! finally the matrices are computed in the CGto basis-

                do igto = 1, Centers(i)%block(iy)%ngto
                 do jgto = 1, Centers(j)%block(jy)%ngto

                  do kgto = 1, Centers(k)%block(ky)%ngto
                   do lgto = 1, Centers(l)%block(ly)%ngto

                      kicgto = 0; kjcgto = 0; kkcgto = 0; klcgto = 0
                      kirep = 0 ; kjrep = 0; kkrep = 0; klrep = 0


                      do icgto = 1, Centers(i)%block(iy)%ncgtoreal(igto)                           
                       do jcgto = 1, Centers(j)%block(jy)%ncgtoreal(jgto)

                        do kcgto = 1, Centers(k)%block(ky)%ncgtoreal(kgto)
                         do lcgto = 1, Centers(l)%block(ly)%ncgtoreal(lgto)


                           itrans = Centers(i)%block(iy)%nincgto(icgto,igto)
                           jtrans = Centers(j)%block(jy)%nincgto(jcgto,jgto) 
                           ktrans = Centers(k)%block(ky)%nincgto(kcgto,kgto) 
                           ltrans = Centers(l)%block(ly)%nincgto(lcgto,lgto) 

                            do im = 1, nli
                             do jm = 1, nlj
                              do km = 1, nlk
                               do lm = 1, nll
 

           
                                kicgto = sum(Centers(1:i-1)%nsizeBlocks) + sum(Centers(i)%nsizePerBlock(1:iy-1)) + (itrans-1)*(nli) + im
                                kjcgto = sum(Centers(1:j-1)%nsizeBlocks) + sum(Centers(j)%nsizePerBlock(1:jy-1)) + (jtrans-1)*(nlj) + jm
                                kkcgto = sum(Centers(1:k-1)%nsizeBlocks) + sum(Centers(k)%nsizePerBlock(1:ky-1)) + (ktrans-1)*(nlk) + km
                                klcgto = sum(Centers(1:l-1)%nsizeBlocks) + sum(Centers(l)%nsizePerBlock(1:ly-1)) + (ltrans-1)*(nll) + lm


                                kirep = sum(Centers(1:i-1)%nsizegtoBlocks) + sum(Centers(i)%nsizegtoPerBlock(1:iy-1)) + (igto-1)*(nli) + im
                                kjrep = sum(Centers(1:j-1)%nsizegtoBlocks) + sum(Centers(j)%nsizegtoPerBlock(1:jy-1)) + (jgto-1)*(nlj) + jm
                                kkrep = sum(Centers(1:k-1)%nsizegtoBlocks) + sum(Centers(k)%nsizegtoPerBlock(1:ky-1)) + (kgto-1)*(nlk) + km
                                klrep = sum(Centers(1:l-1)%nsizegtoBlocks) + sum(Centers(l)%nsizegtoPerBlock(1:ly-1)) + (lgto-1)*(nll) + lm



!-----------save-repulsion-matrices------------------------

                                SrepCGTO(klcgto,kkcgto,kjcgto,kicgto) =  SrepCGTO(klcgto,kkcgto,kjcgto,kicgto) &
                                   & + repGto1(klrep,kkrep,kjrep,kirep)&
                                   & *  Centers(i)%block(iy)%cgto(itrans,igto)*Centers(j)%block(jy)%cgto(jtrans,jgto)&
                                   & *  Centers(k)%block(ky)%cgto(ktrans,kgto)*Centers(l)%block(ly)%cgto(ltrans,lgto)


                            
                             
                               enddo
                              enddo
                             enddo
                            enddo

                         enddo
                        enddo
                       enddo
                      enddo
                      

                    !  deallocate(repmat,repGto,fkx,fky,fkz,hkx,hky,hkz,normGtoI,normGtoJ,normGtoK,normGtoL)


                   enddo  !lgto
                  enddo   !kgto
                 enddo    !jgto
                enddo     !igto
               
 
            enddo   !ly
           enddo    !l
 
         enddo      !ky
        enddo       !k
 
      enddo         !jy
     enddo          !j

   enddo            !iy
  enddo             !i


!write(*,'(16f10.5)') real(repGto1)

  deallocate (repGto1)






  end subroutine StaMatRep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine CollMatRep(geom,izgrid,ax,ay,az,bx,by,bz,vx,vy,vz,ACenters,nACenters,BCenters,nBCenters,CCenters,nCCenters,DCenters,nDCenters,SrepCGtoTP,Krs)

  character(len=1), intent(in) :: Krs

  type(center), dimension(:), intent(inout) :: ACenters
  type(center), dimension(:), intent(inout) :: BCenters
  type(center), dimension(:), intent(inout) :: CCenters
  type(center), dimension(:), intent(inout) :: DCenters

  double complex, dimension(:,:,:,:), intent(inout) :: SrepCGtoTP

  integer, intent(in) :: nACenters, nBCenters, nCCenters, nDCenters


  double precision, intent(in) :: ax,ay,az,bx,by,bz,vx,vy,vz
  integer, intent(in) :: izgrid

  integer, dimension(*) :: geom

  integer :: i, j, k, l, iy, jy, ky, ly, li, lj, lk, ll, nli, nlj, nlk, nll  !, k2ij,k2kl,k2cij,k2ckl
  integer :: igto, jgto, kgto, lgto, icgto, jcgto, kcgto, lcgto, im, jm, km, lm, kicgto, kjcgto, kkcgto, klcgto
  integer :: lmax,lmax1,lmax2, tmax, nsizeGtoA, nsizeGtoB, nsizeGtoC, nsizeGtoD
  integer :: itrans, jtrans, ktrans, ltrans  !, mi, mj, mk, ml, kirep,  kjrep,  kkrep,  klrep


!!$OMP THREADPRIVATE ( i, j, k, l, iy, jy, ky, ly, li, lj, lk, ll, nli, nlj, nlk, nll )
!!$OMP THREADPRIVATE ( igto, jgto, kgto, lgto, icgto, jcgto, kcgto, lcgto, im, jm, km, lm, kicgto, kjcgto, kkcgto, klcgto)
!!$OMP THREADPRIVATE ( lmax,lmax1,lmax2, tmax, nsizeGtoA, nsizeGtoB, nsizeGtoC, nsizeGtoD )
!!$OMP THREADPRIVATE ( itrans, jtrans, ktrans, ltrans )




  double precision :: xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl, xpik, ypik, zpik, xpjl, ypjl, zpjl
  double precision :: alpi, alpj, alpk, alpl, Gamik, Gamjl, rabik2, rabjl2, kfactorik, kfactorjl 
  double precision :: cx2c, cy2c, cz2c

!!$OMP THREADPRIVATE (xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl, xpik, ypik, zpik, xpjl, ypjl, zpjl) 
!!$OMP THREADPRIVATE (alpi, alpj, alpk, alpl, Gamik, Gamjl, rabik2, rabjl2, kfactorik, kfactorjl, cx2c, cy2c, cz2c )

!  double precision :: vx, vy, vz

  double complex, dimension(:,:,:,:,:,:), allocatable :: repmat

  double complex, dimension(:,:,:,:), allocatable ::  repGto  !, repGto1
  double precision, dimension(:,:,:), allocatable ::  fkx, fky, fkz, hkx, hky, hkz
  double precision, dimension(:,:,:), allocatable :: normGtoI, normGtoJ, normGtoK, normGtoL


!!$OMP THREADPRIVATE (repmat, repGto,fkx, fky, fkz, hkx, hky, hkz,normGtoI, normGtoJ, normGtoK, normGtoL)

 

           

  do i = 1, nACenters
   do iy = 1, nYlmMax

    if(ACenters(i)%block(iy)%ncgto==0) cycle

    li = iy - 1
    nli = (li+1)*(li+2)/2

     do j = 1, nBCenters
      do jy = 1, nYlmMax

        if(BCenters(j)%block(jy)%ncgto==0) cycle

        lj = jy - 1
        nlj = (lj+1)*(lj+2)/2

       
        do k = 1, nCCenters
         do ky = 1, nYlmMax

          if(CCenters(k)%block(ky)%ncgto==0) cycle

          lk = ky - 1
          nlk = (lk+1)*(lk+2)/2

           do l = 1, nDCenters
            do ly = 1, nYlmMax

             if(DCenters(l)%block(ly)%ncgto==0) cycle

             ll = ly - 1
             nll = (ll+1)*(ll+2)/2


                do igto = 1, ACenters(i)%block(iy)%ngto
                 do jgto = 1, BCenters(j)%block(jy)%ngto

                  do kgto = 1, CCenters(k)%block(ky)%ngto
                   do lgto = 1, DCenters(l)%block(ly)%ngto



                      xi = ACenters(i)%x(geom(i)); yi = ACenters(i)%y(geom(i)); zi = ACenters(i)%z(geom(i),izgrid)
                      xj = BCenters(j)%x(geom(j)); yj = BCenters(j)%y(geom(j)); zj = BCenters(j)%z(geom(j),izgrid)

                      xk = CCenters(k)%x(geom(k)); yk = CCenters(k)%y(geom(k)); zk = CCenters(k)%z(geom(k),izgrid)
                      xl = DCenters(l)%x(geom(l)); yl = DCenters(l)%y(geom(l)); zl = DCenters(l)%z(geom(l),izgrid)


                      alpi = ACenters(i)%block(iy)%gto(igto); alpj = BCenters(j)%block(jy)%gto(jgto)
                      alpk = CCenters(k)%block(ky)%gto(kgto); alpl = DCenters(l)%block(ly)%gto(lgto)

                      Gamik =  alpi + alpk
                      Gamjl =  alpj + alpl

                      rabik2 = (xi-xk)**2 + (yi-yk)**2 + (zi-zk)**2
                      rabjl2 = (xj-xl)**2 + (yj-yl)**2 + (zj-zl)**2

                      kfactorik = exp(-alpi*alpk*rabik2/Gamik)
                      kfactorjl = exp(-alpj*alpl*rabjl2/Gamjl)


      
            
                      xpik = (alpi*xi + alpk*xk)/Gamik ; ypik = (alpi*yi + alpk*yk)/Gamik ; zpik = (alpi*zi + alpk*zk)/Gamik
                      xpjl = (alpj*xj + alpl*xl)/Gamjl ; ypjl = (alpj*yj + alpl*yl)/Gamjl ; zpjl = (alpj*zj + alpl*zl)/Gamjl

                      cx2c = xpik - xpjl ; cy2c = ypik - ypjl ; cz2c = zpik - zpjl  ! the distance between the centers of two gaussians
!nico China 2019 change ci in -ci to be consistent with fdn
!      cx2c = -cx2c; cy2c = - cy2c; cz2c = - cz2c  
 
                      lmax1 = (li+lk) ; lmax2=(lj+ll) 

                      lmax = max(lmax1 , lmax2)   ; tmax = 6*lmax ! How to choose the Lmax and Tmax

!-----repmat-repulsion-integrals------------------------------------------------------------------------------


                      allocate( repmat(0:lmax,0:lmax,0:lmax,0:lmax,0:lmax,0:lmax) )
                      allocate(normGtoI(0:lmax,0:lmax,0:lmax),normGtoJ(0:lmax,0:lmax,0:lmax))
                      allocate(normGtoK(0:lmax,0:lmax,0:lmax),normGtoL(0:lmax,0:lmax,0:lmax))

                      call normgto(lmax,alpi,normGtoI)
                      call normgto(lmax,alpj,normGtoJ)
                      call normgto(lmax,alpk,normGtoK)
                      call normgto(lmax,alpl,normGtoL)


                      if ( (kfactorik .lt. 1d-8) .or. (kfactorjl .lt. 1d-8) )then
                              repmat=0.d0
                              goto  9
                      endif  
!---------------------repulsion-integrals-------------a,b are used to determine the types of ETF----------------------
                       call repulsion_driver(lmax,0,Gamik,ax,ay,az,Gamjl,bx,by,bz,vx,vy,vz,cx2c,cy2c,cz2c,tmax,repmat)



                      ! if(Krs == 'a' .or.Krs == 'A'  ) then 
                     !     !write(*,*) 
                     !    repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) 
                     !  endif  

                       SELECT CASE(Krs)

                              CASE('a')

                               repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) 

                              CASE('b')
                               repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) &
                                                      & * EXP( dcmplx(0d0,-1.d0*vz*zpik) )
                                                      !!!---exp(2ivRs) comes from the shift of ETF, the velocity always on z axis!--
                              CASE('c')
                               repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) &
                                                      & * EXP( dcmplx(0d0,-2.d0*vz*zpik) ) 
                              CASE('e')
                               repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) &
                                                      & * EXP( dcmplx(0d0,+1.d0*vz*zpik) )
                       ENDSELECT

9                    continue

                     !  if(Krs == 'b' .or.Krs == 'B'  ) then 

                    !      repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) &
                     !       & * EXP( dcmplx(0d0,-1.d0*vz*zpik) ) !---exp(2ivRs) comes from the shift of ETF, the velocity always on z axis!--

                     !  endif 

                    !   if(Krs == 'c' .or.Krs == 'C'  ) then 
                   !       repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) &
                   !          & * EXP( dcmplx(0d0,-2.d0*vz*zpik) ) !---exp(2ivRs) comes from the shift of ETF, the velocity always on z axis!--


                    !   endif 

                     !  if(Krs == 'd' .or.Krs == 'D'  ) then 
                      !    repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) &
                      !       & * EXP( dcmplx(0d0,+2.d0*vz*zpik) ) !---exp(2ivRs) comes from the shift of ETF, the velocity always on z axis!--


                      ! endif 

                    !   if(Krs == 'e' .or.Krs == 'E'  ) then 
                    !      repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) &
                    !         & * EXP( dcmplx(0d0,+1.d0*vz*zpik) ) !---exp(2ivRs) comes from the shift of ETF, the velocity always on z axis!--


                    !   endif 


 !                     if(Krs == 'd' .or.Krs == 'D'  ) then 
!                         repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) 
!                      endif                                                                           
                         
!                      repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:)! &
                           ! & * EXP( dcmplx(0d0,2.d0*vz*xpik) ) !---exp(2ivRs) comes from the shift of ETF, the velocity always on z axis!--
                                         
                     
! fkx, fky, fkz and hkx, hky, hkz contain the coeff. from Gaussian product rule---------------


                      allocate(fkx(0:lmax,0:lmax,0:lmax),fky(0:lmax,0:lmax,0:lmax),fkz(0:lmax,0:lmax,0:lmax))
                      allocate(hkx(0:lmax,0:lmax,0:lmax),hky(0:lmax,0:lmax,0:lmax),hkz(0:lmax,0:lmax,0:lmax))
          

                      call fdn(xi-xpik,xk-xpik,lmax,fkx) 
                      call fdn(yi-ypik,yk-ypik,lmax,fky) 
                      call fdn(zi-zpik,zk-zpik,lmax,fkz) 
                      call fdn(xj-xpjl,xl-xpjl,lmax,hkx) 
                      call fdn(yj-ypjl,yl-ypjl,lmax,hky) 
                      call fdn(zj-zpjl,zl-zpjl,lmax,hkz) 

! now the matrices are computed according to the Gaussian product rule for all
! (nli+1)*(nli+2)/2 and (nlj+1)*(nlj+2)/2 cartesian gaussian functions

                      allocate(  repGto(nll,nlk,nlj,nli)  ) 

                      call sumfkhk(lmax,li,lj,lk,ll,nli,nlj,nlk,nll,fkx,fky,fkz,hkx,hky,hkz,normGtoI,normGtoJ,normGtoK,normGtoL,repmat,repGto)


 


! finally the matrices are computed in the CGto basis-               

                      kicgto = 0; kjcgto = 0; kkcgto = 0; klcgto = 0

                      do icgto = 1, ACenters(i)%block(iy)%ncgtoreal(igto)                           
                       do jcgto = 1, BCenters(j)%block(jy)%ncgtoreal(jgto)

                        do kcgto = 1, CCenters(k)%block(ky)%ncgtoreal(kgto)
                         do lcgto = 1, DCenters(l)%block(ly)%ncgtoreal(lgto)


                           itrans = ACenters(i)%block(iy)%nincgto(icgto,igto)
                           jtrans = BCenters(j)%block(jy)%nincgto(jcgto,jgto) 
                           ktrans = CCenters(k)%block(ky)%nincgto(kcgto,kgto) 
                           ltrans = DCenters(l)%block(ly)%nincgto(lcgto,lgto) 


                            do im = 1, nli
                             do jm = 1, nlj
                              do km = 1, nlk
                               do lm = 1, nll

                                kicgto = sum(ACenters(1:i-1)%nsizeBlocks) + sum(ACenters(i)%nsizePerBlock(1:iy-1)) + (itrans-1)*(nli) + im
                                kjcgto = sum(BCenters(1:j-1)%nsizeBlocks) + sum(BCenters(j)%nsizePerBlock(1:jy-1)) + (jtrans-1)*(nlj) + jm
                                kkcgto = sum(CCenters(1:k-1)%nsizeBlocks) + sum(CCenters(k)%nsizePerBlock(1:ky-1)) + (ktrans-1)*(nlk) + km
                                klcgto = sum(DCenters(1:l-1)%nsizeBlocks) + sum(DCenters(l)%nsizePerBlock(1:ly-1)) + (ltrans-1)*(nll) + lm



!-----------save-repulsion-matrices------------------------


                                SrepCGtoTP(klcgto,kkcgto,kjcgto,kicgto) =  SrepCGtoTP(klcgto,kkcgto,kjcgto,kicgto) &
                                   & +  dconjg(repGto(lm,km,jm,im) )&
                                   & *  ACenters(i)%block(iy)%cgto(itrans,igto)*BCenters(j)%block(jy)%cgto(jtrans,jgto)&
                                   & *  CCenters(k)%block(ky)%cgto(ktrans,kgto)*DCenters(l)%block(ly)%cgto(ltrans,lgto)



                            
                             
                               enddo
                              enddo
                             enddo
                            enddo

                         enddo
                        enddo
                       enddo
                      enddo

                     deallocate(repmat,repGto,fkx,fky,fkz,hkx,hky,hkz,normGtoI,normGtoJ,normGtoK,normGtoL)


                   enddo  !lgto
                  enddo   !kgto
                 enddo    !jgto
                enddo     !igto

 
            enddo   !ly
           enddo    !l
 
         enddo      !ky
        enddo       !k
 
      enddo         !jy
     enddo          !j

   enddo            !iy
  enddo             !i


  end subroutine CollMatRep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine CollMatRepTypeA(geom,izgrid,ax,ay,az,bx,by,bz,vx,vy,vz,ACenters,nACenters,BCenters,nBCenters,CCenters,nCCenters,DCenters,nDCenters,SrepCGtoTP)



  type(center), dimension(:), intent(inout) :: ACenters
  type(center), dimension(:), intent(inout) :: BCenters
  type(center), dimension(:), intent(inout) :: CCenters
  type(center), dimension(:), intent(inout) :: DCenters

  double complex, dimension(:,:,:,:), intent(inout) :: SrepCGtoTP

  integer,intent(inout) :: nACenters, nBCenters, nCCenters, nDCenters
  double precision, intent(in) :: ax,ay,az,bx,by,bz,vx,vy,vz
  integer, intent(in) :: izgrid

  integer, dimension(*) :: geom

  integer :: i, j, k, l, iy, jy, ky, ly, li, lj, lk, ll, nli, nlj, nlk, nll!, k2ij,k2kl,k2cij,k2ckl
  integer :: igto, jgto, kgto, lgto, icgto, jcgto, kcgto, lcgto, im, jm, km, lm, kicgto, kjcgto, kkcgto, klcgto
  integer :: lmax,lmax1,lmax2, tmax, nsizeGtoA, nsizeGtoB, nsizeGtoC, nsizeGtoD
  integer :: itrans, jtrans, ktrans, ltrans, mi, mj, mk, ml, kirep,  kjrep,  kkrep,  klrep


!!$OMP THREADPRIVATE ( i, j, k, l, iy, jy, ky, ly, li, lj, lk, ll, nli, nlj, nlk, nll )
!!$OMP THREADPRIVATE ( igto, jgto, kgto, lgto, icgto, jcgto, kcgto, lcgto, im, jm, km, lm, kicgto, kjcgto, kkcgto, klcgto)
!!$OMP THREADPRIVATE ( lmax,lmax1,lmax2, tmax, nsizeGtoA, nsizeGtoB, nsizeGtoC, nsizeGtoD )
!!$OMP THREADPRIVATE ( itrans, jtrans, ktrans, ltrans, mi, mj, mk, ml, kirep,  kjrep,  kkrep,  klrep )



  double precision :: xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl, xpik, ypik, zpik, xpjl, ypjl, zpjl
  double precision :: alpi, alpj, alpk, alpl, Gamik, Gamjl, rabik2, rabjl2, kfactorik, kfactorjl 
  double precision :: cx2c, cy2c, cz2c


!!$OMP THREADPRIVATE (xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl, xpik, ypik, zpik, xpjl, ypjl, zpjl) 
!!$OMP THREADPRIVATE (alpi, alpj, alpk, alpl, Gamik, Gamjl, rabik2, rabjl2, kfactorik, kfactorjl, cx2c, cy2c, cz2c )


!  double precision :: vx, vy, vz

  double complex, dimension(:,:,:,:,:,:), allocatable :: repmat

  double complex, dimension(:,:,:,:), allocatable ::  repGto, repGto1
  double precision, dimension(:,:,:), allocatable ::  fkx, fky, fkz, hkx, hky, hkz

  double precision, dimension(:,:,:), allocatable :: normGtoI, normGtoJ, normGtoK, normGtoL


!!$OMP THREADPRIVATE (repmat, repGto1, repGto,fkx, fky, fkz, hkx, hky, hkz,normGtoI, normGtoJ, normGtoK, normGtoL)



   nsizeGtoA = sum(ACenters(:)%nsizegtoBlocks)
   nsizeGtoB = sum(BCenters(:)%nsizegtoBlocks)
   nsizeGtoC = sum(CCenters(:)%nsizegtoBlocks)
   nsizeGtoD = sum(DCenters(:)%nsizegtoBlocks)

   allocate( repGto1(nsizeGtoD,nsizeGtoC,nsizeGtoB,nsizeGtoA )  ) 

!   repGto1=0.d0

            

  do i = 1, nACenters
   do iy = 1, nYlmMax

    if(ACenters(i)%block(iy)%ncgto==0) cycle

    li = iy - 1
    nli = (li+1)*(li+2)/2

     do j = 1, nBCenters
      do jy = 1, nYlmMax

        if(BCenters(j)%block(jy)%ncgto==0) cycle

        lj = jy - 1
        nlj = (lj+1)*(lj+2)/2

       
        do k = 1, nCCenters
         do ky = 1, nYlmMax

          if(CCenters(k)%block(ky)%ncgto==0) cycle

          lk = ky - 1
          nlk = (lk+1)*(lk+2)/2

           do l = 1, nDCenters
            do ly = 1, nYlmMax

             if(DCenters(l)%block(ly)%ncgto==0) cycle

             ll = ly - 1
             nll = (ll+1)*(ll+2)/2


                do igto = 1, ACenters(i)%block(iy)%ngto
                 do jgto = 1, BCenters(j)%block(jy)%ngto

                  do kgto = igto, CCenters(k)%block(ky)%ngto
                   do lgto = jgto, DCenters(l)%block(ly)%ngto
  



                      xi = ACenters(i)%x(geom(i)); yi = ACenters(i)%y(geom(i)); zi = ACenters(i)%z(geom(i),izgrid)
                      xj = BCenters(j)%x(geom(j)); yj = BCenters(j)%y(geom(j)); zj = BCenters(j)%z(geom(j),izgrid)

                      xk = CCenters(k)%x(geom(k)); yk = CCenters(k)%y(geom(k)); zk = CCenters(k)%z(geom(k),izgrid)
                      xl = DCenters(l)%x(geom(l)); yl = DCenters(l)%y(geom(l)); zl = DCenters(l)%z(geom(l),izgrid)


                      alpi = ACenters(i)%block(iy)%gto(igto); alpj = BCenters(j)%block(jy)%gto(jgto)
                      alpk = CCenters(k)%block(ky)%gto(kgto); alpl = DCenters(l)%block(ly)%gto(lgto)

                      Gamik =  alpi + alpk
                      Gamjl =  alpj + alpl

                      rabik2 = (xi-xk)**2 + (yi-yk)**2 + (zi-zk)**2
                      rabjl2 = (xj-xl)**2 + (yj-yl)**2 + (zj-zl)**2

                      kfactorik = exp(-alpi*alpk*rabik2/Gamik)
                      kfactorjl = exp(-alpj*alpl*rabjl2/Gamjl)
            
                      xpik = (alpi*xi + alpk*xk)/Gamik ; ypik = (alpi*yi + alpk*yk)/Gamik ; zpik = (alpi*zi + alpk*zk)/Gamik
                      xpjl = (alpj*xj + alpl*xl)/Gamjl ; ypjl = (alpj*yj + alpl*yl)/Gamjl ; zpjl = (alpj*zj + alpl*zl)/Gamjl

                      cx2c = xpik - xpjl ; cy2c = ypik - ypjl ; cz2c = zpik - zpjl  ! the distance between the centers of two gaussians
!nico China 2019 change ci in -ci to be consistent with fdn
!      cx2c = -cx2c; cy2c = - cy2c; cz2c = - cz2c  
 
                      lmax1 = (li+lk) ; lmax2=(lj+ll) 

                      lmax = max(lmax1 , lmax2)   ; tmax = 6*lmax ! How to choose the Lmax and Tmax

!-----repmat-repulsion-integrals------------------------------------------------------------------------------


                      allocate( repmat(0:lmax,0:lmax,0:lmax,0:lmax,0:lmax,0:lmax) )
                      allocate(normGtoI(0:lmax,0:lmax,0:lmax),normGtoJ(0:lmax,0:lmax,0:lmax))
                      allocate(normGtoK(0:lmax,0:lmax,0:lmax),normGtoL(0:lmax,0:lmax,0:lmax))

                      call normgto(lmax,alpi,normGtoI)
                      call normgto(lmax,alpj,normGtoJ)
                      call normgto(lmax,alpk,normGtoK)
                      call normgto(lmax,alpl,normGtoL)


                      if ( (kfactorik .lt. 1d-8) .or. (kfactorjl .lt. 1d-8) )then
                              repmat=0.d0
                              goto  9
                      endif  
  
!---------------------repulsion-integrals-------------a,b are used to determine the types of ETF----------------------
                      call repulsion_driver(lmax,0,Gamik,ax,ay,az,Gamjl,bx,by,bz,vx,vy,vz,cx2c,cy2c,cz2c,tmax,repmat)

                          repmat(:,:,:,:,:,:) = kfactorik*kfactorjl*repmat(:,:,:,:,:,:) 


9                     continue
                                         
                     
! fkx, fky, fkz and hkx, hky, hkz contain the coeff. from Gaussian product rule---------------
                      allocate(fkx(0:lmax,0:lmax,0:lmax),fky(0:lmax,0:lmax,0:lmax),fkz(0:lmax,0:lmax,0:lmax))
                      allocate(hkx(0:lmax,0:lmax,0:lmax),hky(0:lmax,0:lmax,0:lmax),hkz(0:lmax,0:lmax,0:lmax))
          

                      call fdn(xi-xpik,xk-xpik,lmax,fkx) 
                      call fdn(yi-ypik,yk-ypik,lmax,fky) 
                      call fdn(zi-zpik,zk-zpik,lmax,fkz) 
                      call fdn(xj-xpjl,xl-xpjl,lmax,hkx) 
                      call fdn(yj-ypjl,yl-ypjl,lmax,hky) 
                      call fdn(zj-zpjl,zl-zpjl,lmax,hkz) 

! now the matrices are computed according to the Gaussian product rule for all
! (nli+1)*(nli+2)/2 and (nlj+1)*(nlj+2)/2 cartesian gaussian functions

                      allocate(  repGto(nll,nlk,nlj,nli)  ) 

                      call sumfkhk(lmax,li,lj,lk,ll,nli,nlj,nlk,nll,fkx,fky,fkz,hkx,hky,hkz,normGtoI,normGtoJ,normGtoK,normGtoL,repmat,repGto)


                      kirep=0;kjrep=0;kkrep=0;klrep=0

                      do mi=1,nli
                        do mj=1,nlj
                          do mk=1,nlk
                            do ml=1,nll

                                kirep = sum(ACenters(1:i-1)%nsizegtoBlocks) + sum(ACenters(i)%nsizegtoPerBlock(1:iy-1)) + (igto-1)*(nli) + mi
                                kjrep = sum(BCenters(1:j-1)%nsizegtoBlocks) + sum(BCenters(j)%nsizegtoPerBlock(1:jy-1)) + (jgto-1)*(nlj) + mj
                                kkrep = sum(CCenters(1:k-1)%nsizegtoBlocks) + sum(CCenters(k)%nsizegtoPerBlock(1:ky-1)) + (kgto-1)*(nlk) + mk
                                klrep = sum(DCenters(1:l-1)%nsizegtoBlocks) + sum(DCenters(l)%nsizegtoPerBlock(1:ly-1)) + (lgto-1)*(nll) + ml

                                repGto1 (klrep,kkrep,kjrep,kirep)  =  repGto(ml,mk,mj,mi)  

!-----------------------------------------interchange basis on electron 1------------------------------------------------------------------  


                                repGto1 (klrep,kirep,kjrep,kkrep)  =   repGto(ml,mk,mj,mi)  
!-----------------------------------------interchange basis on electron 2------------------------------------------------------------------  

                                repGto1 (kjrep,kkrep,klrep,kirep)  =   repGto(ml,mk,mj,mi) 

!-----------------------------------------interchange basis on both electron 1 & 2------------------------------------------------------------------  

                                repGto1 (kjrep,kirep,klrep,kkrep)  =   repGto(ml,mk,mj,mi)  

                              

                            enddo         
                          enddo   
                        enddo   
                      enddo  

                      deallocate(repmat,repGto,fkx,fky,fkz,hkx,hky,hkz,normGtoI,normGtoJ,normGtoK,normGtoL)

                   enddo  !lgto
                  enddo   !kgto
                 enddo    !jgto
                enddo     !igto


            enddo   !ly
           enddo    !l
 
         enddo      !ky
        enddo       !k
 
      enddo         !jy
     enddo          !j

   enddo            !iy
  enddo             !i




!-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

! finally the matrices are computed in the CGto basis-

  do i = 1, nACenters
   do iy = 1, nYlmMax

    if(ACenters(i)%block(iy)%ncgto==0) cycle

    li = iy - 1
    nli = (li+1)*(li+2)/2

     do j = 1, nBCenters
      do jy = 1, nYlmMax

        if(BCenters(j)%block(jy)%ncgto==0) cycle

        lj = jy - 1
        nlj = (lj+1)*(lj+2)/2

       
        do k = 1, nCCenters
         do ky = 1, nYlmMax

          if(CCenters(k)%block(ky)%ncgto==0) cycle

          lk = ky - 1
          nlk = (lk+1)*(lk+2)/2

           do l = 1, nDCenters
            do ly = 1, nYlmMax

             if(DCenters(l)%block(ly)%ncgto==0) cycle

             ll = ly - 1
             nll = (ll+1)*(ll+2)/2

                do igto = 1, ACenters(i)%block(iy)%ngto
                 do jgto = 1, BCenters(j)%block(jy)%ngto

                  do kgto = 1, CCenters(k)%block(ky)%ngto
                   do lgto = 1, DCenters(l)%block(ly)%ngto


               

                      kicgto = 0; kjcgto = 0; kkcgto = 0; klcgto = 0

                      do icgto = 1, ACenters(i)%block(iy)%ncgtoreal(igto)                           
                       do jcgto = 1, BCenters(j)%block(jy)%ncgtoreal(jgto)

                        do kcgto = 1, CCenters(k)%block(ky)%ncgtoreal(kgto)
                         do lcgto = 1, DCenters(l)%block(ly)%ncgtoreal(lgto)


                           itrans = ACenters(i)%block(iy)%nincgto(icgto,igto)
                           jtrans = BCenters(j)%block(jy)%nincgto(jcgto,jgto) 
                           ktrans = CCenters(k)%block(ky)%nincgto(kcgto,kgto) 
                           ltrans = DCenters(l)%block(ly)%nincgto(lcgto,lgto) 


                            do im = 1, nli
                             do jm = 1, nlj
                              do km = 1, nlk
                               do lm = 1, nll

                                kicgto = sum(ACenters(1:i-1)%nsizeBlocks) + sum(ACenters(i)%nsizePerBlock(1:iy-1)) + (itrans-1)*(nli) + im
                                kjcgto = sum(BCenters(1:j-1)%nsizeBlocks) + sum(BCenters(j)%nsizePerBlock(1:jy-1)) + (jtrans-1)*(nlj) + jm
                                kkcgto = sum(CCenters(1:k-1)%nsizeBlocks) + sum(CCenters(k)%nsizePerBlock(1:ky-1)) + (ktrans-1)*(nlk) + km
                                klcgto = sum(DCenters(1:l-1)%nsizeBlocks) + sum(DCenters(l)%nsizePerBlock(1:ly-1)) + (ltrans-1)*(nll) + lm


                                kirep = sum(ACenters(1:i-1)%nsizegtoBlocks) + sum(ACenters(i)%nsizegtoPerBlock(1:iy-1)) + (igto-1)*(nli) + im
                                kjrep = sum(BCenters(1:j-1)%nsizegtoBlocks) + sum(BCenters(j)%nsizegtoPerBlock(1:jy-1)) + (jgto-1)*(nlj) + jm
                                kkrep = sum(CCenters(1:k-1)%nsizegtoBlocks) + sum(CCenters(k)%nsizegtoPerBlock(1:ky-1)) + (kgto-1)*(nlk) + km
                                klrep = sum(DCenters(1:l-1)%nsizegtoBlocks) + sum(DCenters(l)%nsizegtoPerBlock(1:ly-1)) + (lgto-1)*(nll) + lm


!-----------save-repulsion-matrices------------------------


                                SrepCGtoTP(klcgto,kkcgto,kjcgto,kicgto) =  SrepCGtoTP(klcgto,kkcgto,kjcgto,kicgto) &
                                   & +  dconjg(repGto1(klrep,kkrep,kjrep,kirep) )&
                                   & *  ACenters(i)%block(iy)%cgto(itrans,igto)*BCenters(j)%block(jy)%cgto(jtrans,jgto)&
                                   & *  CCenters(k)%block(ky)%cgto(ktrans,kgto)*DCenters(l)%block(ly)%cgto(ltrans,lgto)



                            
                             
                               enddo
                              enddo
                             enddo
                            enddo

                         enddo
                        enddo
                       enddo
                      enddo



                   enddo  !lgto
                  enddo   !kgto
                 enddo    !jgto
                enddo     !igto

 
            enddo   !ly
           enddo    !l
 
         enddo      !ky
        enddo       !k
 
      enddo         !jy
     enddo          !j

   enddo            !iy
  enddo             !i

  deallocate (repGto1)


  end subroutine CollMatRepTypeA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  subroutine sumfkhk(lmax,li,lj,lk,ll,nli,nlj,nlk,nll,fkx,fky,fkz,hkx,hky,hkz,normGtoI,normGtoJ,normGtoK,normGtoL,mat,mgto)
  implicit none
  integer, intent(in) :: lmax, li,lj,lk,ll,nli,nlj,nlk,nll
  double complex, dimension(0:lmax,0:lmax,0:lmax,0:lmax,0:lmax,0:lmax)  :: mat
  double precision, dimension(0:lmax,0:lmax,0:lmax) :: fkx, fky, fkz, hkx, hky, hkz, normGtoI, normGtoJ, normGtoK, normGtoL
  double complex, dimension(nll,nlk,nlj,nli) :: mgto


  
  integer, dimension(3,nli) :: indi
  integer, dimension(3,nlj) :: indj
  integer, dimension(3,nlk) :: indk
  integer, dimension(3,nll) :: indl



  integer :: ic, jc, kc, lc
  integer :: ix, iy, iz, jx, jy, jz, kx, ky, kz, lx ,ly, lz
  integer :: ikx, iky, ikz, jlx, jly, jlz


!!$OMP THREADPRIVATE (ic, jc, kc, lc, ix, iy, iz, jx, jy, jz, kx, ky, kz, lx ,ly, lz, ikx, iky, ikz, jlx, jly, jlz)


!! order i : x,y,z |  xx, xy, xz, yy, yz, zz | xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz | ....
  
  mgto(:,:,:,:) = 0d0

! creates all possible (l+1)*(l+2)/2  polynomials for a given angular momentum li, lj, lk, ll
  call genPolgto(li,nli,indi)
  call genPolgto(lj,nlj,indj)
  call genPolgto(lk,nlk,indk)
  call genPolgto(ll,nll,indl)
  
  do ic = 1, nli
    ix = indi(1,ic)
    iy = indi(2,ic)
    iz = indi(3,ic)
    do jc = 1, nlj
      jx = indj(1,jc)
      jy = indj(2,jc)
      jz = indj(3,jc)
      do kc = 1, nlk
        kx = indk(1,kc)
        ky = indk(2,kc)
        kz = indk(3,kc)
        do lc = 1, nll
          lx = indl(1,lc)
          ly = indl(2,lc)
          lz = indl(3,lc)

!       WRITE(*,*)ix,iy,iz,"----",kx,ky,kz

          do ikx = 0, ix + kx 
           do iky = 0, iy + ky 
            do ikz = 0, iz + kz 
               do jlx = 0, jx + lx 
                do jly = 0, jy + ly 
                 do jlz = 0, jz + lz 
                 
                   mgto(lc,kc,jc,ic) = mgto(lc,kc,jc,ic) + fkx(ikx,kx,ix)*fky(iky,ky,iy)*fkz(ikz,kz,iz)*hkx(jlx,lx,jx)*hky(jly,ly,jy)*hkz(jlz,lz,jz)&
                                      &*mat(jlz,jly,jlx,ikz,iky,ikx)*normgtoI(ix,iy,iz)*normgtoJ(jx,jy,jz)*normgtoK(kx,ky,kz)*normgtoL(lx,ly,lz)


                 enddo
                enddo
               enddo
            enddo
           enddo
          enddo


        enddo
      enddo
    enddo
  enddo


  end subroutine sumfkhk


 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine printmat(l,m,n,mat)
  integer, intent(in) :: l, m, n
  double complex, dimension(0:n,0:m,0:l), intent(in) :: mat
  integer :: iprint, jprint, kprint 

    do iprint = 0, l
      do jprint = 0, m 
       do kprint = 0, n 
           write(*,'(3i3,100(f20.12,1X))')iprint, jprint, kprint, mat(kprint,jprint,iprint)
       enddo
     enddo
   enddo

  end subroutine printmat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module MatRep
