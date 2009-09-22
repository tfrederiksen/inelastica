!
!     Fortran version of setkpointhelper.

      subroutine f90setkpointhelper(Sparse, maxnh, kpoint, no, nuo, 
     +     numh, rcell, xij, listhptr, listh, indxuo, Full)

      implicit none

!     INPUT:
!     Sparse matrix
      real*8, intent(in) ::  Sparse(0:maxnh-1)        
!     Dimension of sparse matrix
      integer, intent(in) ::  maxnh                    
!     K-point
      real*8, intent(in) ::   kpoint(0:2)              
!     Total number of orbitals in supercell
      integer, intent(in) ::   no                       
!     Orbitals in unitcell
      integer, intent(in) ::   nuo                      
!     Number of nonzero column elements
      integer, intent(in) ::   numh(0:nuo-1)            
!     Reciprocal cell
      real*8, intent(in) ::   rcell(0:2,0:2)           
!     Vector connecting unitcells 
      real*8, intent(in) ::   xij(0:maxnh-1,0:2)         
!     Start of row in Sparse matrix
      integer, intent(in) ::   listhptr(0:nuo-1)        
!     Column number
      integer, intent(in) ::   listh(0:maxnh-1)         
!     Corresponding orbital in unitcell
      integer, intent(in) ::   indxuo(0:no-1)           

!     OUTPUT:
!     Output full matrix
      complex*16, intent(out) :: Full(0:nuo-1,0:nuo-1)   

!     =======================================================================
!     Local variables 

      complex*16   phase(0:maxnh-1)   
!     Temporary array to calculate phases   
      real*8      tmp(0:2)           
!     Loop indecies
      integer      ii,jj,iuo,jz       
!     Temporary
      integer      si, juo            
!     Temp variable 
      real*8      rtmp               
!     
!    def setkpointhelper(self, Sparse, kpoint):
!        """
!        Make full matrices from sparse for specific k-point
!        NOTE: Assumption for Fourier transform                     
!        Psi(i) =           sum_R exp(i k.R) Psi_k(i)               NOTE sign!
!        (Full waveunction)                   (Unit cell part of wavefunction)
!              which gives:
!        i,j part of unitcell gives from the rows of the full H: 
!            H_k(i,j) = H_(i,j) + sum_R H(i,j+R) exp(i*k*R)         NOTE sign!
!            where R corresponds to R_j-R_0 which is Xij            NOTE sign!
!
!        TODO: Efficiency. This is the main timewaster of pyTBT
!        """        

!        # Phase factor 
!        tmp=mm(kpoint,N.transpose(self.rcell),N.transpose(self.xij))

!     kpoint*transpose(rcell)
      tmp=0.0d0
      do ii=0,2
         do jj=0,2
            tmp(ii)=tmp(ii)+kpoint(jj)*rcell(ii,jj)
         enddo
      enddo
!     tmp*transpose(xij)
!     phase=N.exp(2.0j*N.pi*tmp)    # exp(2 pi i k*(Rj-Ri)) where i,j from Hij
      do ii=0,maxnh-1
         rtmp=0.0d0
         do jj=0,2
            rtmp=rtmp+tmp(jj)*xij(ii,jj)
         enddo
         phase(ii)=cdexp(dcmplx(0d0,6.2831853071795864769d0*rtmp))
      enddo

!        Full = N.zeros((self.nuo,self.nuo),N.Complex)

      Full = 0.0d0

!        for iuo in range(self.nuo):
!            for jz in range(self.numh[iuo]):

      do iuo = 0,nuo-1
         do jz = 0, numh(iuo)-1
!                si = self.listhptr[iuo]+jz
            si=listhptr(iuo)+jz
!                juo = self.indxuo[self.listh[si]-1]-1
            juo = indxuo(listh(si)-1)-1
!                Full[iuo,juo]=Full[iuo,juo]+Sparse[si]*phase[si]

            Full(iuo,juo) = Full(iuo,juo) + Sparse(si)*phase(si)
         enddo
      enddo

!
!        return Full
      return 
      end 
