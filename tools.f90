module tools
use general, only : lenmax
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mysearch(sword,channel,found)
!
!     search for sword in file unit = channel
!
character*(*) sword
integer channel
integer found

integer length
character(len=lenmax) test

 test=' '
 found = 0
 length=index(sword,'&')-1
 rewind(unit=channel)

 do while(test(1:length) .ne. sword(1:length))
   read(channel,'(a)',end=71) test(1:length)
   if (test(1:length) .eq. sword(1:length)) return  
 enddo

 71   continue
 found = -1
 
end subroutine mysearch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tools
