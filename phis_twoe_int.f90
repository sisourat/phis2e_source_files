program phis_twoe_int
use general
use colllib_2e
use setup
use centerlib
use inputlib
use cgto
use matrices
!use colldyn
!use collint
!use collint_2e
use twoe_int_cgto
use onee_int_cgto
use OMP_LIB
implicit none

  character(len=lenmax) :: input, option, suff, finputxml, finput
  character(len=lenmax) :: fout, ffout, famp
  double precision :: spac
  double precision :: readtmp
!!!!!!!!!!!
  integer :: ilen,jlen
  integer :: iunitone, iunittwo
  character(len=lenmax) :: filemv
  character(len=lenmax) :: creat_impab, impab
  character(len=lenmax) :: impactb_value

!!!!!!!!!!!!

!  double complex :: temp_tp
!   double complex, dimension(4) :: ttemp_tpvec

  double precision :: t1, t2, tmec, tdyn, ttrans, tsm, ttot
 
  integer :: ib, is, iunit, iunit2,iunit3, izgrid
  integer :: i, j, k   !,ii,jj,kk

  double precision, dimension(:), allocatable :: sta2e,sta1et,sta1ep
  integer, dimension(:), allocatable :: vecordert,vecorderp

  call init_constants

! set up the calculations : directories, input, allocation, etc...

  call getarg(1,input)
  call getarg(2,option)

  suff = "_int_2e "
  call setdir(input,suff,option)
  
  finputxml='input.xml'
  call inputcoll_2e(finputxml) ! read and set up input.xml file

  call phis_cpfiles_2e() ! cp target and projectile eigenvalues, eigenvectors in workdir

! read target eigenvalues/eigenvectors
  write(*,'(/,a)')" # Target data"
  finputxml='tinput.xml'
  call inputsta(finputxml)
  ntCenters = nCenters
  allocate(tCenters(ntCenters))
  tCenters = Centers
  deallocate(Centers)

  do i = 1, ntCenters
    call read_cgto(tCenters(i))
  enddo
  
  ntcgto = sum(tCenters(:)%nsizeBlocks)

  open(unit=10,file='tmo.txt')
  read(10,*)ntmo, ntao
  write(1234,*) "ntmo, ntao=,",ntmo, ntao
  allocate(tmo(ntao,ntmo))

   do i = 1, ntmo
     do j = 1, ntao
       read(10,*) readtmp 
       tmo(j,i)=dcmplx(readtmp, 0.0d0)
       write(1234,*) "ntao=j,ntmo=i,tmo=,",j,i,tmo(j,i)
     enddo
   enddo
  close(10)



!-------------p-----------------------------------
  write(*,'(/,a)')" # Projectile data"
  finputxml='pinput.xml'
  call inputsta(finputxml)
  npCenters = nCenters
  allocate(pCenters(npCenters))
  pCenters = Centers
  deallocate(Centers)

  do i = 1, npCenters
    call read_cgto(pCenters(i))
  enddo

  npcgto = sum(pCenters(:)%nsizeBlocks)

  open(unit=10,file='pmo.txt')
   read(10,*)npmo, npao
  allocate(pmo(npao,npmo))
   do i = 1, npmo
     do j = 1, npao
       read(10,*) readtmp
       pmo(j,i)=dcmplx(readtmp, 0.0d0)
     enddo
   enddo
  close(10)


  ntotmo = ntmo + npmo
! here comes the real stuffs


  vproj = v%vz
  if(v%vx/=0d0 .or. v%vy/=0d0) then
   write(*,*)"Vproj should be along z axis"
   stop 
  endif
  tgrid%na = zgrid%na +1
  allocate(zgrid%a(zgrid%na +1),tgrid%a(tgrid%na))
  if(zgrid%typ=='linear' .or. zgrid%typ=='Linear' .or. zgrid%typ=='LINEAR') then

    do i = 1, zgrid%na +1
      zgrid%a(i) = zgrid%amin + (i-1)*(zgrid%amax-zgrid%amin)/zgrid%na
      tgrid%a(i) = zgrid%a(i)/vproj
    enddo

  elseif(zgrid%typ=='exp' .or. zgrid%typ=='Exp' .or. zgrid%typ=='EXP') then

    spac=0d0 

    do i = 1, int(DBLE(zgrid%na)/DBLE(2)) !zgrid%na/2
     spac=spac+1.1d0**DBLE(i)
    enddo
     spac=(zgrid%amax-zgrid%amin)/(2d0*spac)

     !zgrid%a(zgrid%na/2 +1)=0d0
     !zgrid%a(int(DBLE(zgrid%na)/DBLE(2))+1) = 0.0d0
     zgrid%a = 0.0d0
     j = 0
    do i = int(DBLE(zgrid%na)/DBLE(2))+2, zgrid%na
      j = j + 1
      zgrid%a(i) = zgrid%a(i-1)+1.1d0**DBLE(j)*spac
      zgrid%a(zgrid%na +2 -i) = -zgrid%a(i)
    enddo
      zgrid%a(1) = zgrid%amin
      zgrid%a(zgrid%na+1) = zgrid%amax

    do i = 1, zgrid%na +1
      tgrid%a(i) = zgrid%a(i)/vproj
      write(888,*) i, zgrid%a(i), tgrid%a(i)
    enddo

   else
     write(*,*)"No other z grids implemented"
     stop
   endif

  nbproj = b%n
  write(1234,*) "ccjia, code-2e: nbproj = b%n,",nbproj, b%n
  allocate(bproj(1:nbproj))
  do i = 1, nbproj
   if(b%typ=='linear' .or. b%typ=='Linear' .or. b%typ=='LINEAR') then
    bproj(i) = b%bmin + (i-1)*(b%bmax-b%bmin)/b%n
    write(1234,*) "ccjia, code-2e: i = 1, nbproj. bproj(i) ",i,bproj(i)

!----add impact parameters input by hand---26 Oct. 2017--Junwen-------

   elseif(b%typ=='input')then
    bproj(i)=b%binp(i)
   else
     write(*,*)"No other b grids implemented"
     stop
   endif
  enddo

!  ntotcgto = ntcgto*(ntcgto+1)/2 + ntcgto*npcgto + npcgto*(npcgto+1)/2
  ntotcgto =  ntcgto+npcgto

  write(*,'(/,a)')"# Running twoe int"
!-----------------------------------------------------------------------
  allocate(mcoup(zgrid%na +1,ntmo+npmo,ntmo+npmo),movl(zgrid%na +1,ntmo+npmo,ntmo+npmo))
  allocate(mcgtocoup(zgrid%na +1,ntcgto+npcgto,ntcgto+npcgto),mcgtoovl(zgrid%na +1,ntcgto+npcgto,ntcgto+npcgto))

  do ib = 1, nbproj
!!!!!!!
    write(impactb_value,"(f10.2)") bproj(ib)
    impab="impactb_"
    creat_impab= Trim(AdjustL(impab))//Trim(AdjustL(impactb_value))
    ilen=index(creat_impab,' ')
    call system('mkdir '// creat_impab(1:ilen-1))
!!!!!!!
 
    pCenters(1)%x(1) = bproj(ib) 
    pCenters(1)%y(1) = 0d0

   write(*,*)'## Compute 2e ints'

!   call cpu_time(t1)

   do izgrid = 1, zgrid%na +1
     pCenters(1)%z(1,izgrid) = zgrid%a(izgrid)
   enddo

!call omp_set_num_threads(20)

   call system("date")
   !   t1=OMP_get_wtime()

!-------!OMP time loops are inside the subroutine now.---Junwen, DEC. 12, 2017.--------------

   call Compute_1eint_cgto()
   call inttransmo
   open(unit=10,file='onee_int_tt.txt',status='new')
    write(10,*)bproj(ib)
    do izgrid = 1, zgrid%na +1
     write(10,*)zgrid%a(izgrid)
       do i = 1, ntmo
         do j = 1, ntmo
             write(10,'(2(i5,1X),4(f30.16,3X))')j,i,mcoup(izgrid,j,i),movl(izgrid,j,i)
         enddo
       enddo
    enddo
   close(10)

   open(unit=10,file='onee_int_pp.txt',status='new')
    write(10,*)bproj(ib)
    do izgrid = 1, zgrid%na +1
     write(10,*)zgrid%a(izgrid)
       do i = 1, npmo
         do j = 1, npmo
             write(10,'(2(i5,1X),4(f30.16,3X))')j,i,mcoup(izgrid,j+ntmo,i+ntmo),movl(izgrid,j+ntmo,i+ntmo)
         enddo
       enddo
    enddo
   close(10)
!!! !!! means MAT(j,i) to (Tj,Ti) (Pj,Pi)
!!! ccjia change: tp to pt, and pt to tp
   open(unit=10,file='onee_int_tp.txt',status='new')
   open(unit=11,file='onee_int_pt.txt',status='new')
    write(10,*)bproj(ib)
    write(11,*)bproj(ib)
    do izgrid = 1, zgrid%na +1
     write(10,*)zgrid%a(izgrid)
     write(11,*)zgrid%a(izgrid)
       do i = 1, npmo
         do j = 1, ntmo
             write(10,'(2(i5,1X),4(f30.16,1X))')j,i,mcoup(izgrid,j,i+ntmo),movl(izgrid,j,i+ntmo)
             !write(11,'(2(i5,1X),4(f25.16,1X))')j,i,mcoup(izgrid,j+ntmo,i),movl(izgrid,j+ntmo,i)
         enddo
       enddo
       do i = 1, ntmo
         do j = 1, npmo
             !write(10,'(2(i5,1X),4(f25.16,1X))')i,j,mcoup(izgrid,i,j+ntmo),movl(izgrid,i,j+ntmo)
             write(11,'(2(i5,1X),4(f30.16,1X))')j,i,mcoup(izgrid,j+ntmo,i),movl(izgrid,j+ntmo,i)
         enddo
       enddo
    enddo
   close(10)
   close(11)
!!!
   filemv= "onee_int_*"
   jlen=index(filemv,' ')
   call system('mv '//filemv(1:jlen-1)//' '//creat_impab(1:ilen-1))
   !call integrals_filemv(filemv,creat_impab)
!!!
   call Compute_2eint_cgto(ib,bproj(ib))
!!!
   filemv= "twoe_int_*"
   jlen=index(filemv,' ')
   call system('mv '//filemv(1:jlen-1)//' '//creat_impab(1:ilen-1))
   !call integrals_filemv(filemv,creat_impab)
!!!
   call system("date")
   !t2=OMP_get_wtime()

   tmec = t2-t1

   write(*,*)'Computing 2e ints takes',tmec

  enddo ! ib
  
  !do is = 1, ninitsta
  ! iunit = 10 + is
  ! close(iunit)
  !enddo

  deallocate(bproj)
  deallocate(tmo,pmo,mcoup,movl)

!  deallocate(ttcgtocoup,ttcgtoovl,tpcgtocoup,tpcgtoovl,ptcgtocoup,ptcgtoovl,ppcgtocoup,ppcgtoovl)
  deallocate( b%binp)
  call freeCenters()
  call freeColl()

  write(*,*)
  call system("date")

end program phis_twoe_int
