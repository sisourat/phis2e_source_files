module matrices
use omp_lib

  integer :: nsizeCGto
  double complex, dimension(:,:),  allocatable :: ovCGto, potCGto, kinCGto, tCGto,Ov2eCGTO_2d,pot2eCGTO_2d,repCGto_2d,kin2eCGTO_2d,t2eCGto_2d
  double complex, dimension(:,:,:,:), allocatable :: SrepCGto, repCGto , t2eCGTO, Ov2eCGTO, pot2eCGTO, kin2eCGTO, Coup2eCGTO   !, tptpvec

!!$OMP  THREADPRIVATE( SrepCGto, repCGto , t2eCGTO, Ov2eCGTO, pot2eCGTO, kin2eCGTO, Coup2eCGTO )

!!$OMP  THREADPRIVATE( SrepCGto, repCGto , t2eCGTO, Ov2eCGTO, pot2eCGTO, kin2eCGTO, Coup2eCGTO)

  integer :: ntotsta, ntotcgto

!  double complex, dimension(:,:,:), save, allocatable :: ttcgtoovl
!  double complex, dimension(:,:,:), save,  allocatable :: ppcgtoovl


  double precision, dimension(:,:), allocatable :: tcoup,pcoup



  double complex, dimension(:,:,:),  allocatable :: mcgtocoup, mcgtoovl, m2ecgtocoup, m2ecgtoovl
  double complex, dimension(:,:,:),  allocatable :: mcoup, movl, m2ecoup, m2eovl
!!$OMP  THREADPRIVATE( m2ecgtocoup, m2ecgtoovl)


!!$OMP  THREADPRIVATE( ttcgtocoup, ttcgtoovl, tpcgtocoup, tpcgtoovl, ppcgtocoup, ppcgtoovl, ptcgtocoup, ptcgtoovl )

!!$OMP  THREADPRIVATE( ttcgtocoup,tpcgtocoup,ptcgtocoup,ppcgtocoup )



! double complex, dimension(:,:), allocatable :: ttttcgtocoup,ttttcgtoovl
! double complex, dimension(:,:), allocatable :: ppppcgtocoup,ppppcgtoovl
! double complex, dimension(:,:), allocatable :: tptpcgtocoup,tptpcgtoovl
! double complex, dimension(:,:), allocatable :: tttpcgtocoup,tttpcgtoovl
! double complex, dimension(:,:), allocatable :: ttppcgtocoup,ttppcgtoovl
! double complex, dimension(:,:), allocatable :: tpppcgtocoup,tpppcgtoovl
! double complex, dimension(:,:), allocatable :: ppttcgtocoup,ppttcgtoovl
! double complex, dimension(:,:), allocatable :: pptpcgtocoup,pptpcgtoovl
! double complex, dimension(:,:), allocatable :: tpttcgtocoup,tpttcgtoovl





  double complex, dimension(:,:,:,:), allocatable :: SrepCGto1,SrepCGto2,SrepCGto3,SrepCGto4

!!$OMP  THREADPRIVATE(SrepCGto1,SrepCGto2,SrepCGto3,SrepCGto4)

end module matrices
