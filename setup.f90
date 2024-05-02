module setup
use general, only : lenmax
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setdir(namedir,suffdir,option)
implicit none

character(lenmax), intent(in) ::  namedir, suffdir, option
character(lenmax) ::  finput
character(lenmax) :: creadir, workdir
integer :: ilen, jlen
logical :: dir_e, file_e
integer :: i, j

jlen=index(namedir,' ')
finput=namedir(1:jlen-1)//'.xml'
jlen=index(finput,' ')

inquire( file="./"//finput(1:jlen-1), exist=file_e )
if ( file_e .eqv. .false. ) then
 write(*,*) finput(1:jlen-1), " does not exist"
 stop
endif

ilen=index(namedir,' ')
creadir=namedir(1:ilen-1)//suffdir
ilen=index(creadir,' ')

write(*,*)"# Output files in ", creadir(1:ilen-1)

inquire( file="./"//creadir(1:ilen-1)//"/.", exist=dir_e )
if ( dir_e ) then
  if(option=="-rep") then
    write(*,*),"# ", creadir(:ilen-1), " replaced!"
    call system('rm -rf '// creadir(:ilen-1))
    call system('mkdir '// creadir(:ilen-1))
    call system('cp '//finput(1:jlen-1)//' '//creadir(:ilen-1)//'/input.xml')
    call chdir(creadir(:ilen-1))
    call getcwd( workdir)
    write(*,*)"# Move to ", workdir
  else
    write(*,*), creadir(:ilen-1), " exists!"
    write(*,*), "if you want to replace it, use option -rep"
    stop
  endif
else
  call system('mkdir '// creadir(:ilen-1))
  call system('cp '//finput(1:jlen-1)//' '//creadir(:ilen-1)//'/input.xml')
  call chdir(creadir(:ilen-1))
  call getcwd( workdir)
  write(*,*)"# Move to ", workdir
end if

end subroutine setdir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module setup
