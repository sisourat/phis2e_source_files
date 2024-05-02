module SlaterDeterminant
implicit none

type Sdeterminant
 integer :: nalpha, nbeta
 integer, dimension(:), allocatable :: alpha, beta
end type Sdeterminant

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function compute_det(n,mat)
implicit none

integer :: n, lda, info, nperm
double complex, dimension(n,n) :: mat
double complex :: compute_det, det
integer, dimension(n) :: ipiv
integer :: i

call ZGETRF(n,n,mat,lda,ipiv,info)

! number of perm
nperm=0
det=1d0
do i = 1, n
 if(ipiv(i)/=i) nperm=nperm+1
 det=det*mat(i,i)
enddo
compute_det = (-1d0)**nperm*(det)

end function compute_det

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lowdin(ne,nmo,ovmo,h1emo,r12mo,det1,det2,ov,h1e,r12)
implicit none

!! compute the matrix elements between two Slater determinants (det1 and
!det2) formed with non-orthogonal MOs
!! for efficiency, the matrix elements in the MO basis (ovmo, h1emo and
!r12mo) are global variables
!! ne is the number of electrons
!! ov, h1e, r12 are the matrix elements between the two determinants
!det1 and det2

!! FOR EFFICIENCY, r12mo SHOULD BE CHANGED AS GLOBAL VARIABLES (AVOID
!COPYING LARGE ARRAYS)

integer, intent(in) :: ne, nmo
double complex, dimension(nmo,nmo), intent(in) :: ovmo, h1emo
double complex, dimension(nmo,nmo,nmo,nmo), intent(in) :: r12mo
type(Sdeterminant), intent(in) :: det1, det2

double complex, intent(out) :: ov, h1e, r12

double complex, dimension(ne) :: vectemp
double complex, dimension(ne,ne) :: ovmat, ovstore
double complex, dimension(ne-1,ne-1) :: comat
double complex, dimension(ne-2,ne-2) :: comat2

integer :: i, j , k, l, ia, ib, ja, jb, ka, kb, la, lb

!! compute and store overlap matrix between the determinants
ovmat(:,:) = 0d0
do i = 1, det1%nalpha
 do j = 1, det2%nalpha
  ja =  det2%alpha(j)
  ia =  det1%alpha(i)
  ovmat(j,i) = ovmo(ja,ia)
 enddo
enddo

do i = 1, det1%nbeta
 do j = 1, det2%nbeta
  jb =  det2%beta(j)
  ib =  det1%beta(i)
  ovmat(j+det2%nalpha,i+det1%nalpha) = ovmo(jb,ib)
 enddo
enddo

ovstore(:,:) = ovmat(:,:)

!! compute the overlap between the two determinants
ov = compute_det(ne,ovmat)

!! compute the 1e hamiltonian between the two determinants
h1e = 0d0
do i = 1, det1%nalpha
 do j = 1, det2%nalpha
  ovmat(:,:) = ovstore(:,:)
  ovmat(j,:) = ovmat(ne,:)
  ovmat(:,i) = ovmat(:,ne)
  comat(:,:) = ovmat(1:ne-1,1:ne-1)
 
!  comat(1:j-1,1:i-1) = ovstore(1:j-1,1:i-1)
!  comat(j:ne,i:ne) = ovstore(j+1:ne,i+1:ne)
  ja =  det2%alpha(j)
  ia =  det1%alpha(i)
  h1e = h1e + (-1d0)**(i+j)*h1emo(ja,ia)*compute_det(ne-1,comat)
 enddo
enddo

do i = 1, det1%nbeta
 do j = 1, det2%nbeta
  ovmat(:,:) = ovstore(:,:)
  ovmat(j+det2%nalpha,:) = ovmat(ne,:)
  ovmat(:,i+det1%nalpha) = ovmat(:,ne)
  comat(:,:) = ovmat(1:ne-1,1:ne-1)

!  comat(1:det2%nalpha,1:det1%nalpha) =
!  ovstore(1:det2%nalpha,1:det1%nalpha)
!  comat(det2%nalpha+1:j-1,det1%nalpha+1:i-1) =
!  ovstore(det2%nalpha+1:j-1,det1%nalpha+1:i-1)
!  comat(det2%nalpha+j:ne,det1%nalpha+i:ne) =
!  ovstore(det2%nalpha+j+1:ne,det1%nalpha+i+1:ne)

  jb =  det2%beta(j)
  ib =  det1%beta(i)
  h1e = h1e + (-1d0)**(det1%nalpha+i+det2%nalpha+j)*h1emo(jb,ib)*compute_det(ne-1,comat)
 enddo
enddo

!! compute the 2e repulsion between the two determinants
r12 = 0d0

!! alpha el for e1 and e2
do i = 1, det1%nalpha  !e1
 do j = 1, det2%nalpha  !e1
  ja =  det2%alpha(j)
  ia =  det1%alpha(i)

   do k = 1, det1%nalpha  !e2
     do l = 1, det2%nalpha  !e2
       la =  det2%alpha(l)
       ka =  det1%alpha(k)

       ovmat(:,:) = ovstore(:,:)
       ovmat(:,j) = ovmat(:,ne)
       ovmat(:,i) = ovmat(:,ne-1)
       ovmat(l,:) = ovmat(ne,:)
       ovmat(k,:) = ovmat(ne-1,:)
       comat2(:,:) = ovmat(1:ne-2,1:ne-2)

       r12 = r12 + r12mo(la,ka,ja,ia)*(-1d0)**(i+j+k+l)*compute_det(ne-2,comat2)

     enddo
   enddo

 enddo
enddo

!! alpha el for e1 and beta el for e2
do i = 1, det1%nalpha  !e1
 do j = 1, det2%nalpha  !e1
  ja =  det2%alpha(j)
  ia =  det1%alpha(i)

   do k = 1, det1%nbeta  !e2
     do l = 1, det2%nbeta  !e2
       lb =  det2%beta(l)
       kb =  det1%beta(k)
       ovmat(:,:) = ovstore(:,:)
       ovmat(:,j) = ovmat(:,ne)
       ovmat(:,i) = ovmat(:,ne-1)
       ovmat(l+det2%nalpha,:) = ovmat(ne,:)
       ovmat(k+det1%nalpha,:) = ovmat(ne-1,:)
       comat2(:,:) = ovmat(1:ne-2,1:ne-2)

       r12 = r12 + r12mo(lb,kb,ja,ia)*(-1d0)**(i+j+k+l+det1%nalpha+det2%nalpha)*compute_det(ne-2,comat2)

     enddo
   enddo

 enddo
enddo

!! beta el for e1 and alpha el for e2
do i = 1, det1%nbeta  !e1
 do j = 1, det2%nbeta  !e1
  jb =  det2%beta(j)
  ib =  det1%beta(i)

   do k = 1, det1%nalpha  !e2
     do l = 1, det2%nalpha  !e2
       la =  det2%alpha(l)
       ka =  det1%alpha(k)
       ovmat(:,:) = ovstore(:,:)
       ovmat(:,j+det2%nalpha) = ovmat(:,ne)
       ovmat(:,i+det1%nalpha) = ovmat(:,ne-1)
       ovmat(l,:) = ovmat(ne,:)
       ovmat(k,:) = ovmat(ne-1,:)
       comat2(:,:) = ovmat(1:ne-2,1:ne-2)

       r12 = r12 + r12mo(la,ka,jb,ib)*(-1d0)**(i+j+k+l+det1%nalpha+det2%nalpha)*compute_det(ne-2,comat2)

     enddo
   enddo

 enddo
enddo

!! beta el for e1 and alpha el for e2
do i = 1, det1%nbeta  !e1
 do j = 1, det2%nbeta  !e1
  jb =  det2%beta(j)
  ib =  det1%beta(i)

   do k = 1, det1%nbeta  !e2
     do l = 1, det2%nbeta  !e2
       lb =  det2%beta(j)
       kb =  det1%beta(i)
       ovmat(:,:) = ovstore(:,:)
       ovmat(:,j+det2%nalpha) = ovmat(:,ne)
       ovmat(:,i+det1%nalpha) = ovmat(:,ne-1)
       ovmat(l+det2%nalpha,:) = ovmat(ne,:)
       ovmat(k+det1%nalpha,:) = ovmat(ne-1,:)
       comat2(:,:) = ovmat(1:ne-2,1:ne-2)

       r12 = r12 + r12mo(lb,kb,jb,ib)*(-1d0)**(i+j+k+l+2*det1%nalpha+2*det2%nalpha)*compute_det(ne-2,comat2)

     enddo
   enddo

 enddo
enddo

end subroutine lowdin


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module SlaterDeterminant
