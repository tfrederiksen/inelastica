!
!     Fortran version of removeUnitCellXij.
!
!    def removeUnitCellXij(self):
!        """
!        Remove displacements within unitcell from xij
!        NOTE: We remove the in cell difference so xij corresponds to 
!              lattice vectors to the relevant part of the supercell.
!        NOTE: xij = Rj-Ri where Ri,j corresponds to positions of the orbitals H_{i,j} 
!        TODO: Check why some orbitals in sparse matrix reported within cell but have xij!
!        """
!        for iuo in range(self.nuo):
!            for jnz in range(self.numh[iuo]):
!                juo = self.indxuo[self.listh[self.listhptr[iuo]+jnz]-1]-1
!                ia,ja = self.atomindx[iuo]-1, self.atomindx[juo]-1
!                self.xij[self.listhptr[iuo]+jnz,:] = self.xij[self.listhptr[iuo]+jnz,:]-\
!                    (self.xa[ja,:]-self.xa[ia,:])

subroutine f90removeunitcellxij(nnzs, no_u, na_u, &
     numh, xij, xa, listhptr, listh, atomindx, xijo)

  implicit none

! INPUT:
! Dimension of sparse matrix
  integer, intent(in)  :: nnzs                    
! Orbitals in unitcell
  integer, intent(in)  :: no_u                     
! Atoms in unitcell
  integer, intent(in)  :: na_u                      
! Number of nonzero column elements
  integer, intent(in)  :: numh(no_u) 
! Vector connecting unitcells 
  real*8, intent(in) :: xij(nnzs,3)         
! Position of atom in unitcell 
  real*8, intent(in) :: xa(na_u,3)         
! Start of row in Sparse matrix
  integer, intent(in)  :: listhptr(no_u)        
! Column number
  integer, intent(in)  :: listh(nnzs)         
! Corresponding orbital in unitcell
  integer, intent(in)  :: atomindx(no_u)

! OUTPUT:
! Output full matrix
  real*8, intent(out) :: xijo(nnzs,3)         

! Loop indecies
  integer :: iuo, j, ind, juo
  
  do iuo = 1 , no_u
     do j = 1 , numh(iuo)
        ind = listhptr(iuo) + j
        juo = idxuo(listh(ind),no_u)
        xijo(ind,:) = xij(ind,:) - &
             (xa(atomindx(juo),:)-xa(atomindx(iuo),:))
     enddo
  enddo
      
contains
  
  pure function idxuo(io,no_u) result(uo)
    integer, intent(in) :: io,no_u
    integer :: uo
    uo = mod(io,no_u)
    if ( uo == 0 ) uo = no_u
  end function idxuo

end subroutine f90removeunitcellxij
