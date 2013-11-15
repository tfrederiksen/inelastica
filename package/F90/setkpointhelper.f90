!
!     Fortran version of setkpointhelper.

subroutine f90setkpointhelper(nnzs, Sparse, kpoint, no_u, &
     numh, rcell, xij, listhptr, listh, Full)
  
  implicit none

! INPUT:
! Dimension of sparse matrix
  integer, intent(in)  :: nnzs
! Sparse matrix
  real*8, intent(in) :: Sparse(nnzs)
! K-point
  real*8, intent(in) :: kpoint(3)
! Orbitals in unitcell
  integer, intent(in)  :: no_u
! Number of nonzero column elements
  integer, intent(in)  :: numh(no_u)
! Reciprocal cell
  real*8, intent(in) :: rcell(3,3) 
! Vector connecting unitcells 
  real*8, intent(in) :: xij(nnzs,3) 
! Start of row in Sparse matrix
  integer, intent(in)  :: listhptr(no_u) 
! Column number
  integer, intent(in)  :: listh(nnzs) 

! OUTPUT:
! Output full matrix
  complex*16, intent(out) :: Full(no_u,no_u) 

! Local variables
! Temporary array to calculate phases 
  real*8 :: tmp(3) 
! Loop indecies
  integer :: i, j
  integer :: iuo, juo, ind
! Temp variable
  real*8 :: kxij

  ! convert the reciprocal cell together with the Pi-prefactor
  do i = 1 , 3
     tmp(i) = 3.141592653589793238_8*2._8*sum(kpoint(:)*rcell(i,:))
  end do

  ! initialize
  Full(:,:) = 0._8
  
  do iuo = 1 , no_u
     do j = 1 , numh(iuo)
        ind = listhptr(iuo) + j
        juo = idxuo(listh(ind),no_u)
        
        kxij = &
             tmp(1) * xij(ind,1) + &
             tmp(2) * xij(ind,2) + &
             tmp(3) * xij(ind,3)
        
        Full(iuo,juo) = Full(iuo,juo) + Sparse(ind)*exp(dcmplx(0._8,kxij))
        
     end do
  end do
  
contains

  pure function idxuo(io,no_u) result(uo)
    integer, intent(in) :: io,no_u
    integer :: uo
    uo = mod(io,no_u)
    if ( uo == 0 ) uo = no_u
  end function idxuo

end subroutine f90setkpointhelper
