module general
implicit none

integer, parameter :: nYlmMax = 7
character(14), dimension(nYlmMax), parameter :: YlmBlock = ['! s functions&', '! p functions&', '! d functions&', '! f functions&', '! g functions&', '! h functions&', '! i functions&' ] 
integer, parameter :: lenmax = 200

double precision, dimension(-1:4*nYlmMax+3) :: fact, fact2
double precision, dimension(-1:2*nYlmMax,0:2*nYlmMax) ::  cij

double precision, parameter :: pi = dacos(-1d0)
double complex, parameter :: imag = dcmplx(0d0,1d0)

 contains

 subroutine init_constants

 integer :: i, j

   fact(:)=1d0
  do i=1,2*nYlmMax
   fact(i)=fact(i-1)*DBLE(i)
  enddo

    fact2(:)=1d0
  do i=2,2*nYlmMax
    fact2(i)=DBLE(i)    
    do j=i-2,2,-2
     fact2(i)=fact2(i)*DBLE(j)
    enddo
  end do

  do i = 0, 2*nYlmMax
    do j = 0, i! 2*nYlmMax
          cij(i,j) = fact(i)/(fact(j)*fact(i-j))
    enddo
  enddo

 end subroutine init_constants

 !subroutine integrals_filemv(filemv,creat_impab)
 
  !integer :: ilen,jlen
  !character(len=lenmax),intent(in) :: filemv
  !character(len=lenmax),intent(in) :: creat_impab

   !ilen=index(creat_impab,' ')
   !jlen=index(filemv,' ')
   !call system('mv '//filemv(1:jlen-1)//' '//creat_impab(1:ilen-1)) 

 !end subroutine integrals_filemv

end module general
