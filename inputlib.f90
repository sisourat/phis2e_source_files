!! created by Nico 23.07.2014
!! contains getstaInput, getstaEnd

module inputlib
use general, only : lenmax
use newtypes
use flib_dom
use centerlib
implicit none


 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inputsta(finput)
  character(len=lenmax) :: k, finput

  type(fnode), pointer     :: myInput, CenterNode, PositionNode, FileNode, PotNode
  type(fnodeList), pointer :: myList, CenterList, PositionList, FileList, PotList

  integer :: i, j, ipos, ipara
  
  write(*,*)finput

  myInput => parsefile(finput,verbose=.false.)
 
  CenterList => getElementsByTagName(myInput,"Center")
  nCenters = getLength(CenterList)

  allocate(Centers(nCenters))

  do i = 0, nCenters-1
   CenterNode => item(CenterList,i)

! Position Center i
   PositionList => getElementsByTagName(CenterNode,"Position")
   Centers(i+1)%npos = getLength(PositionList)

   allocate(Centers(i+1)%x(Centers(i+1)%npos),Centers(i+1)%y(Centers(i+1)%npos) )

   allocate(  Centers(i+1)%z( Centers(i+1)%npos, 1000))
   do ipos = 0, Centers(i+1)%npos-1
     PositionNode => item(PositionList,ipos)
     k = getAttribute(PositionNode,"x"); read (k,*) Centers(i+1)%x(ipos+1)
     k = getAttribute(PositionNode,"y"); read (k,*) Centers(i+1)%y(ipos+1)
     k = getAttribute(PositionNode,"z"); read (k,*) Centers(i+1)%z(ipos+1,1)
!---------------------------------
     do ipara=2,1000
      Centers(i+1)%z(ipos+1,ipara) = Centers(i+1)%z(ipos+1,1)
     enddo
!---------------------------------

   enddo

! Basis set file on Center i
   FileList => getElementsByTagName(CenterNode,"Basisset")
   FileNode => item(FileList,0)
   Centers(i+1)%fbasis = getAttribute(FileNode,"file")

! Potential on Center i
   PotList => getElementsByTagName(CenterNode,"Potential")
   Centers(i+1)%pot%n = getLength(PotList)

   allocate(Centers(i+1)%pot%c(Centers(i+1)%pot%n), Centers(i+1)%pot%alp(Centers(i+1)%pot%n))

   do j = 0, Centers(i+1)%pot%n-1
      PotNode => item(PotList,j)
      k = getAttribute(PotNode,"charge"); read (k,*) Centers(i+1)%pot%c(j+1)

      !Centers(i+1)%pot%c(j+1)=
      k = getAttribute(PotNode,"exponent"); read (k,*) Centers(i+1)%pot%alp(j+1)
      
   enddo
  enddo

  end subroutine inputsta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module inputlib
