module colllib_2e
use general, only : lenmax
use newtypes, only : impactparameter, velocity, defgrid
use flib_dom
implicit none

 type(impactparameter) :: b
 type(velocity) :: v

 integer :: ninitsta
 integer, dimension(:), allocatable :: initsta
 
 character(len=lenmax) :: tinput, pinput

 integer :: ntsta, npsta, ntpsta, ntgto, npgto, ntpgto, ntcgto, npcgto
 integer :: tstai, tstaf, pstai, pstaf, tpstai, tpstaf

 double complex, dimension(:), allocatable :: testa, pesta, tpesta
 double precision, dimension(:,:), allocatable :: teigvec, peigvec, tpeigvec

 double complex, dimension(:), allocatable :: ttesta, ppesta
 double precision, dimension(:,:), allocatable :: tteigvec, ppeigvec

 integer :: nbproj
 double precision :: vproj
 double precision, dimension(:), allocatable :: bproj

 type(defgrid) :: zgrid, tgrid

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inputcoll_2e(finput)
  character(len=lenmax) :: k, finput

  character(len=lenmax),dimension(:),allocatable :: inistat_typ


  type(fnode), pointer     :: myInput, Node, InitNode, TargetNode,ProjNode,Tar_ProjNode, GridNode
  type(fnodeList), pointer :: List, ImpactVel, ImpactParam, InitList, TargetList, ProjList, Tar_ProjList, GridList

  integer :: i

  myInput => parsefile(finput,verbose=.false.)

! Get the Collision parameters

  ImpactParam => getElementsByTagName(myInput,"ImpactParam")
  ImpactVel => getElementsByTagName(myInput,"ImpactVel")
  GridList => getElementsByTagName(myInput,"Zgrid")
  InitList => getElementsByTagName(myInput,"InitState")

! print the input info

  Node => item(ImpactParam,0)
  b%typ = getAttribute(Node,"type")

  k = getAttribute(Node,"bmin"); read (k,*) b%bmin
  k = getAttribute(Node,"bmax"); read (k,*) b%bmax
  k = getAttribute(Node,"nb"); read (k,*) b%n

!----add impact parameters input by hand---26 Oct. 2017--Junwen--------

  allocate( b%binp(b%n))  
  if(b%typ.eq."input") then
 
  if(b%n.le.20) then
    k = getAttribute(Node,"binp1to20"); read (k,*) b%binp(1:b%n)
  elseif((b%n.GT.20).and.(b%n.LE.40))then
    k = getAttribute(Node,"binp1to20");  read (k,*) b%binp(1:20)
    k = getAttribute(Node,"binp21to40"); read (k,*) b%binp(21:b%n)
  elseif((b%n.GT.40).and.(b%n.LE.60))then
    k = getAttribute(Node,"binp1to20");  read (k,*) b%binp(1:20)
    k = getAttribute(Node,"binp21to40"); read (k,*) b%binp(21:40)
    k = getAttribute(Node,"binp41to60"); read (k,*) b%binp(41:b%n)
  elseif((b%n.GT.60).and.(b%n.LE.80))then
    k = getAttribute(Node,"binp1to20");  read (k,*) b%binp(1:20)
    k = getAttribute(Node,"binp21to40"); read (k,*) b%binp(21:40)
    k = getAttribute(Node,"binp41to60"); read (k,*) b%binp(41:60)
    k = getAttribute(Node,"binp61to80"); read (k,*) b%binp(61:b%n)

  elseif((b%n.GT.80).and.(b%n.LE.100))then
    k = getAttribute(Node,"binp1to20");  read (k,*) b%binp(1:20)
    k = getAttribute(Node,"binp21to40"); read (k,*) b%binp(21:40)
    k = getAttribute(Node,"binp41to60"); read (k,*) b%binp(41:60)
    k = getAttribute(Node,"binp61to80"); read (k,*) b%binp(61:80)
    k = getAttribute(Node,"binp81to100"); read (k,*) b%binp(81:b%n)
  elseif((b%n.GT.100).and.(b%n.LE.120))then
    k = getAttribute(Node,"binp1to20");  read (k,*) b%binp(1:20)
    k = getAttribute(Node,"binp21to40"); read (k,*) b%binp(21:40)
    k = getAttribute(Node,"binp41to60"); read (k,*) b%binp(41:60)
    k = getAttribute(Node,"binp61to80"); read (k,*) b%binp(61:80)
    k = getAttribute(Node,"binp81to100"); read (k,*) b%binp(81:100)
    k = getAttribute(Node,"binp101to120"); read (k,*) b%binp(101:b%n)
  else
    write(*,*) "Impact parameter input is wrong "
    write(*,*) "Or the code does not support to input too many b by hand now." 
  endif

  write(*,*) b%n
  write(*,*)  b%binp(1:b%n)

  endif
   
  if(b%typ=='linear' .or. b%typ=='Linear' .or. b%typ=='LINEAR') then
  write(*,'(a,a10,2(f12.6,1X),i5)')'b grid ',  b%typ, b%bmin, b%bmax, b%n
  endif

  Node => item(ImpactVel,0)
  k = getAttribute(Node,"vx"); read (k,*) v%vx
  k = getAttribute(Node,"vy"); read (k,*) v%vy
  k = getAttribute(Node,"vz"); read (k,*) v%vz
  write(*,'(a,3(f12.6,1X))') 'vx, vy, vz', v

  GridNode => item(GridList,0)
  zgrid%typ = getAttribute(GridNode,"type")
  k = getAttribute(GridNode,"zmin"); read (k,*) zgrid%amin
  k = getAttribute(GridNode,"zmax"); read (k,*) zgrid%amax
  k = getAttribute(GridNode,"nzgrid"); read (k,*) zgrid%na

  write(*,'(a,2(f12.6,1X),i5)') 'zgrid', zgrid%amin, zgrid%amax, zgrid%na

  ninitsta = getLength(InitList)
  write(*,*)ninitsta
  allocate(initsta(ninitsta))
 
  allocate(inistat_typ(ninitsta))
  do i = 0, ninitsta - 1
 
     InitNode => item(InitList,i)
     inistat_typ(i+1) = getAttribute(InitNode,"type")

     k = getAttribute(InitNode,"state"); read (k,*) initsta(i+1)

  enddo

! Target info---TT---------
   TargetList => getElementsByTagName(myInput,"Target")
   TargetNode => item(TargetList,0)
   tinput = getAttribute(TargetNode,"file")
   k = getAttribute(TargetNode,"stamin"); read (k,*)tstai
   k = getAttribute(TargetNode,"stamax"); read (k,*)tstaf
   ntsta = tstaf - tstai + 1
   if(tstaf<=0) ntsta=0

! Target-Proj.info---TP----
   Tar_ProjList => getElementsByTagName(myInput,"Tar_Proj")
   Tar_ProjNode => item(Tar_ProjList,0)
!   pinput = getAttribute(Tar_ProjNode,"file")

   k = getAttribute(Tar_ProjNode,"stamin"); read (k,*)tpstai


   k = getAttribute(Tar_ProjNode,"stamax"); read (k,*)tpstaf
   ntpsta = tpstaf - tpstai + 1

   if(tpstaf<=0) ntpsta=0

! Proj info---PP------------
   ProjList => getElementsByTagName(myInput,"Projectile")
   ProjNode => item(ProjList,0)
   pinput = getAttribute(ProjNode,"file")
   k = getAttribute(ProjNode,"stamin"); read (k,*)pstai
   k = getAttribute(ProjNode,"stamax"); read (k,*)pstaf
   npsta = pstaf - pstai + 1
   if(pstaf<=0) npsta=0


!----reset initial state by Target, Projectile or Tar_Proj------------

  do i = 0, ninitsta - 1

     if(inistat_typ(i+1) == "Tar_Proj")then
     initsta(i+1)=initsta(i+1) + ntsta
     elseif(inistat_typ(i+1) == "Projectile")then
     initsta(i+1)=initsta(i+1) + ntsta+ntpsta
     endif

     write(*,*)"Check the initial state position in the TT, TP, PP: "

     write(*,*) "NO.",initsta(i+1) 

  enddo

  deallocate(inistat_typ)

  end subroutine inputcoll_2e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cpfiles_2e()
! cp target and projectile eigenvalues, eigenvectors in workdir

  character(lenmax) ::  finput
  character(lenmax) :: workdir

  logical :: dir_e, file_e
  integer :: ilen, jlen

! target input xml
  jlen=index(tinput,' ')
  finput='../'//tinput(1:jlen-1)//'/input.xml'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/tinput.xml')

! target eigenvalues
!-------t------------------
  jlen=index(tinput,' ')
  finput='../'//tinput(1:jlen-1)//'/eigenvalues-t'
 
  jlen=index(finput,' ')

  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/teigenvalues')
!-------t----coupling matrices <T|T1+Vt1|T>, used in collint_2e-----
  jlen=index(tinput,' ')

  finput='../'//tinput(1:jlen-1)//'/mcoup-t'
 
  jlen=index(finput,' ')

  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/mcoup-t')

!-------tt-----------------
 jlen=index(tinput,' ')
 finput='../'//tinput(1:jlen-1)//'/eigenvalues-tt'

  jlen=index(finput,' ')

  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/tteigenvalues')

! target eigenvectors
!--------t------------------
  jlen=index(tinput,' ')
  finput='../'//tinput(1:jlen-1)//'/eigenvectors-t'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/teigenvectors')
!----------tt--------------------
  jlen=index(tinput,' ')
  finput='../'//tinput(1:jlen-1)//'/eigenvectors-tt'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/tteigenvectors')

! projectile input xml
  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/input.xml'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/pinput.xml')

! projectile eigenvalues
!---------p----------------------------------
  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/eigenvalues-t'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/peigenvalues')
!-------p------------coupling matrices <P|T1+Vp1|P>, used in collint_2e ----
  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/mcoup-t'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/mcoup-p')


!----------pp--------------------------------------
  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/eigenvalues-tt'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/ppeigenvalues')

! projectile eigenvectors
!--------------p--------------------------------
  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/eigenvectors-t'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/peigenvectors')
!--------------pp------------------------------
  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/eigenvectors-tt'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/ppeigenvectors')


  end subroutine cpfiles_2e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_eig(nsta,ncgto,esta,eigvec,ni,nf,finput,typ)

  integer, intent(in) :: nsta, ncgto
  double complex, dimension(nsta), intent(inout) :: esta
  double precision, dimension(ncgto,nsta), intent(inout) :: eigvec
  integer, intent(in) :: ni, nf

  character(len=lenmax), intent(in) :: finput
  character(len=lenmax) :: a

  double precision :: etmp

  integer :: nv, ne
  integer :: iunit

  integer :: i, j, k, ista

  double precision :: d, ecap

  character(2) :: typ

  d = 1d20
  ecap = 1d20

  iunit = 10
!write(*,*) nsta
!  write(*,*)"# List of states"
  open(unit=iunit,file=finput)

  open(621,file='2e_states')

   read(iunit,*)a

   read(iunit,*)nv, ne

   if(ne/=ncgto) then
     write(*,*)"error in get_eig input"
     write(*,*)"wrong numbers of ncgto in ", finput
     stop
   endif

   if(nv<nsta) then
     write(*,*)"error in get_eig input"
     write(*,*)"wrong numbers of nstates in ", finput
     stop
   endif
   
  
   k = 0
   do i = 1, nv
     read(iunit,*)ista, etmp
     if(ista>=ni .and. ista<=nf) then
       k = k + 1
       if(etmp<ecap) then 
         esta(k) = etmp
       else
! heuristic model to include ionization (see Saalfrank JCP 2009)
!         esta(k) = dcmplx(etmp,-0.5d0*dsqrt(etmp)/d)
         esta(k) = dcmplx(etmp,-0.5d0/d)
       endif

!       if(typ == 't1' .or. typ == 'p1') then

       write(*,'(i5,2(f26.16,4X))')k,esta(k)
    
!      read(*,*) 
!------used to print states inorder. junwen Oct 31, 2017-------------------------------------
       if(typ == 't2' .or. typ == 'p2') then
       write(621,'((f26.16,4X))') real(esta(k))
       endif
!-------------------------------------------


       do j = 1, ne
         read(iunit,*)eigvec(j,k) 
       enddo
     else
       do j = 1, ne
         read(iunit,*)
       enddo
     endif
   enddo
   
  close(iunit)
!  close(621)

end subroutine get_eig



! subroutine get_tpeig(nsta,estat,eigvect,estap,eigvecp,ncgtot,ncgtop,estatp,eigvectp,ni,nf)


! integer, intent(in) :: nsta, ncgtot,ncgtop
! double complex, dimension(ncgtot), intent(inout) :: estat
! double complex, dimension(ncgtop), intent(inout) :: estap
! double complex, dimension(nsta), intent(inout) :: estatp

! double precision, dimension(ncgtot,ncgtot), intent(inout) :: eigvect
! double precision, dimension(ncgtop,ncgtop), intent(inout) :: eigvecp
! double precision, dimension(ncgtot*ncgtop,nsta), intent(inout) :: eigvectp

! integer, dimension(ncgtot*ncgtop,2) :: ordertp
! integer, dimension(2) :: ordertp_temp
!
! integer, intent(in) :: ni, nf
! double complex :: temp_tp
! double complex, dimension(:), allocatable :: temp_tpvec

!
! integer :: nv, ne

! integer :: i, j, k, ista
! integer :: ii,jj,kk

! write(*,*)"# List of states"


!  if(ncgtot*ncgtop < nf) then
!    write(*,*)"error in get_tpeig"
!    write(*,*)"wrong numbers of nstates in for TP "!, finput
!    stop
!  endif


! allocate(tpesta1(ncgtot*ncgtop),tpeigvec1(ncgtot*ncgtop,ncgtot*ncgtop),tpeigvec2(ncgtot*ncgtop,ncgtot*ncgtop),tpeigvec3(ncgtot*ncgtop,ncgtot*ncgtop))
! k=0
! do i=1, ncgtot
!   do ii=1, ncgtop
!      k=k+1;kk=0
!      do j=1,ncgtot
!       do jj=1,ncgtop
!         kk=kk+1
!         tpesta1(k)=estat(i)+estap(ii)
!         ordertp(k,1)=i
!         ordertp(k,2)=ii
!         tpeigvec2(kk,k)=eigvecp(jj,ii)*eigvect(j,i)
!       enddo
!      enddo
!   enddo
! enddo

! k=0
! do i=1, ncgtop
!   do ii=1, ncgtot
!      k=k+1;kk=0
!      do j=1,ncgtop
!       do jj=1,ncgtot
!         kk=kk+1
!
!         tpeigvec3(kk,k)=eigvecp(j,i)*eigvect(jj,ii)
!       enddo
!      enddo
!   enddo
! enddo
! tpeigvec1(:,1)=tpeigvec2(:,1)
! tpeigvec1(:,4)=tpeigvec2(:,4)
! k=0

! do i=1, ncgtot
!   do j=1, ncgtop
!      k=k+1;kk=0

!      do ii=1,ncgtot
!       do jj=1,ncgtop
!         kk=kk+1

!          tpeigvec1(kk,k)=(eigvect(ii,i)*eigvecp(jj,j) + eigvect(jj,j)*eigvecp(ii,i))/sqrt(4.d0)

!         enddo
!      enddo
!   enddo
! enddo

!     write(*,'(4f10.5)') real(tpeigvec1  )


!
! tpeigvec1(:,2:3)= (tpeigvec2(:,2:3) + tpeigvec2(:,2:3))/sqrt(2.d0)

!write(*,'(2f10.5)')eigvect
!write(*,*) 

!write(*,'(2f10.5)')eigvecp

! write(*,*) 'ss'

! write(*,'(4f10.5)') real(tpeigvec2)
!write(*,*) 

!write(*,'(4f10.5)') real(tpeigvec3)
!write(*,*) 
!write(*,'(4f10.5)') (real(tpeigvec3)+ real(tpeigvec2))!/sqrt(2.d0)

!stop
! tpeigvec1(2,:)= (tpeigvec3(3,:))
! tpeigvec3(3,:)= (tpeigvec3(2,:))
! tpeigvec3(2,:)= (tpeigvec1(2,:))

!write(*,'(4f10.5)') real(tpeigvec3)
!stop


! tpeigvec1(3,:)=  ((tpeigvec2(3,:)) +(tpeigvec3(2,:)))/sqrt(2.d0)
! tpeigvec2(:,4)=  ((tpeigvec2(:,4)) +(tpeigvec3(:,4)))/2.d0!/sqrt(2.d0)


!  tpeigvec2(2,:)=tpeigvec1(2,:)
!  tpeigvec2(3,:)=tpeigvec1(3,:)

!write(*,*)sqrt(2.d0)
!tpeigvec2(:,1) =2.d0*tpeigvec2(:,1)/sqrt(2.d0)
!tpeigvec2(:,4) =2.d0*tpeigvec2(:,4)/sqrt(2.d0)

!  tpeigvec2(:,1) = ( tpeigvec2(:,1) +tpeigvec3(:,1) )/sqrt(2.d0)
!  tpeigvec1(:,2)=  ((tpeigvec2(:,2)) +(tpeigvec3(:,2)))/sqrt(2.d0)
!  tpeigvec2(:,3)=  ((tpeigvec2(:,3)) +(tpeigvec3(:,3)))/sqrt(2.d0)
!  tpeigvec2(:,4)=  ((tpeigvec2(:,4)) +(tpeigvec3(:,4)))/sqrt(2.d0)


! tpeigvec2(:,2)=tpeigvec1(:,2)
! tpeigvec2(:,3)=tpeigvec1(:,2)

!  tpeigvec2(:,3)=tpeigvec1(:,3)

!stop

!tpeigvec1=tpeigvec2!*2.d0/sqrt(2.d0)!+tpeigvec3)!/sqrt(2.d0)
!  !(real(tpeigvec3)+ real(tpeigvec2))/sqrt(2.d0)

!write(*,'(4f10.5)') tpeigvec1

!stop
! !tpeigvec3!/sqrt(2.d0)
!stop


! allocate(temp_tpvec(ncgtot*ncgtop))
!---------rerank the order--------------------------
!do i=1, ncgtot*ncgtop-1
!  do j=i+1, ncgtot*ncgtop

!      if( real( tpesta1(i)).gt.real(tpesta1(j)))then
!        temp_tp=tpesta1(j)
!        tpesta1(j)=tpesta1(i)
!        tpesta1(i)=temp_tp
!        ordertp_temp(:)=ordertp(j,:)
!        ordertp(j,:)=ordertp(i,:)
!        ordertp(i,:)=ordertp_temp(:)

!        temp_tpvec=tpeigvec1(:,j)
!        tpeigvec1(:,j)=tpeigvec1(:,i)
!        tpeigvec1(:,i)=temp_tpvec
!      endif

!        if( real( tpesta1(i)).eq.real(tpesta1(j)))then
!        
!          tpeigvec1(i,:)= tpeigvec1(j,:)
!        endif



!   enddo
! enddo


!write(*,*) ordertp(:,2)!, ordertp(:,2)
! deallocate(temp_tpvec)

! estatp(:)=tpesta1(ni:nf)
! eigvectp(:,:)=tpeigvec1(:,ni:nf)

! deallocate(tpesta1,tpeigvec1) 
! do i=1,nsta
!  write(*,'(i4,a,(f20.16,1X),i4,a,(f20.16,1X),i4,a,(f20.16,1X))') i,' TP: ',real(estatp(i)),&
!      & ordertp(i,1),' T:',real(estat(ordertp(i,1))),ordertp(i,2),' P:',real(estap(ordertp(i,2)))!!,ordertp(i,1)
! enddo

!end subroutine get_tpeig


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine freeColl()
   if(allocated(initsta)) deallocate(initsta)
  end subroutine freeColl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine phis_cpfiles_2e()
! cp target and projectile eigenvalues, eigenvectors in workdir

  character(lenmax) ::  finput
  character(lenmax) :: workdir

  logical :: dir_e, file_e
  integer :: ilen, jlen

! target input xml
  jlen=index(tinput,' ')
  finput='../'//tinput(1:jlen-1)//'/tinput.xml'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/tinput.xml')

  jlen=index(tinput,' ')
  finput='../'//tinput(1:jlen-1)//'/taobasis*.bas'
  jlen=index(finput,' ')
!  inquire( file="./"//finput(1:jlen-1), exist=file_e )
!  if ( file_e .eqv. .false. ) then
!    write(*,*) finput(1:jlen-1), " does not exist"
!    stop
!  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/')

  jlen=index(tinput,' ')
  finput='../'//tinput(1:jlen-1)//'/tmobasis.txt'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/tmo.txt')


! projectile input xml
  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/pinput.xml'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/pinput.xml')

  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/pmobasis.txt'
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/pmo.txt')

  jlen=index(pinput,' ')
  finput='../'//pinput(1:jlen-1)//'/paobasis*.bas'
  jlen=index(finput,' ')
!  inquire( file="./"//finput(1:jlen-1), exist=file_e )
!  if ( file_e .eqv. .false. ) then
!    write(*,*) finput(1:jlen-1), " does not exist"
!    stop
!  endif
  call getcwd( workdir)
  ilen=index(workdir,' ')
  call system('cp '//finput(1:jlen-1)//' '//workdir(:ilen-1)//'/')



  end subroutine phis_cpfiles_2e


end module colllib_2e

