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

      subroutine f90removeunitcellxij( maxnh, no, nuo, nua,
     +     numh, xij, xa, listhptr, listh, indxuo, atomindx, xijo)

      implicit none

!     INPUT:
!     Dimension of sparse matrix
      integer, intent(in) ::  maxnh                    
!     Total number of orbitals in supercell
      integer, intent(in) ::   no                       
!     Orbitals in unitcell
      integer, intent(in) ::   nuo                      
!     Atoms in unitcell
      integer, intent(in) ::   nua                      
!     Number of nonzero column elements
      integer, intent(in) ::   numh(0:nuo-1)            
!     Vector connecting unitcells 
      real*8, intent(in) ::   xij(0:maxnh-1,0:2)         
!     Position of atom in unitcell 
      real*8, intent(in) ::   xa(0:nua-1,0:2)         
!     Start of row in Sparse matrix
      integer, intent(in) ::   listhptr(0:nuo-1)        
!     Column number
      integer, intent(in) ::   listh(0:maxnh-1)         
!     Corresponding orbital in unitcell
      integer, intent(in) ::   indxuo(0:no-1)           
      integer, intent(in) ::   atomindx(0:nuo-1)

!     OUTPUT:
!     Output full matrix
      real*8, intent(out) ::   xijo(0:maxnh-1,0:2)         

!     =======================================================================
!     Local variables 

!     Loop indecies
      integer      iuo,jnz       
      integer      juo,ia,ja 

!        for iuo in range(self.nuo):
!            for jnz in range(self.numh[iuo]):
      do iuo=0,nuo-1
         do jnz=0,numh(iuo)-1
!                juo = self.indxuo[self.listh[self.listhptr[iuo]+jnz]-1]-1
            juo=indxuo(listh(listhptr(iuo)+jnz)-1)-1
!                ia,ja = self.atomindx[iuo]-1, self.atomindx[juo]-1
            ia=atomindx(iuo)-1
            ja=atomindx(juo)-1
!                self.xij[self.listhptr[iuo]+jnz,:] = self.xij[self.listhptr[iuo]+jnz,:]-\
!                    (self.xa[ja,:]-self.xa[ia,:])
            xijo(listhptr(iuo)+jnz,:)=xij(listhptr(iuo)+jnz,:)-
     +           (xa(ja,:)-xa(ia,:))
         enddo
      enddo

      return 
      end 
