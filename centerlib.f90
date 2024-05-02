module centerlib
use newtypes, only : center
implicit none

 integer :: nCenters
 type(center), dimension(:),allocatable :: Centers

 integer :: ntCenters
 type(center), dimension(:), allocatable :: tCenters

 integer :: npCenters
 type(center), dimension(:), allocatable :: pCenters


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine freeCenters
    if(allocated(Centers)) deallocate(Centers)
    if(allocated(tCenters)) deallocate(tCenters)
    if(allocated(pCenters)) deallocate(pCenters)
  end subroutine freeCenters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module centerlib
