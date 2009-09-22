!
!     Fortran version of readNewTSHS
!
!def ReadNewTSHSFile(filename):
!    """
!    Read new version of TSHS file.
!    For return see code:
!      Note that onlyS -> does not return xij or Hsparse
!                gamma -> sets indxuo by hand to 1:nou, -1 for nou+1:nos! and sets xij=0.0
!    xa[atomnr,xyz] : Atom positions
!    ucell[nr,xyz]  : Unitcell
!    xij[nr,xyz]
!   """
      module readnewtshs

      implicit none

      PUBLIC
      real*8, dimension(:), allocatable :: S          
      real*8, dimension(:,:), allocatable :: H, xij, xa          
      integer, dimension(:), allocatable :: &
           indxuo, listhptr, numh, listh, lasto, isa

      CONTAINS

      subroutine read(fname, gamma, onlyS, no_u, &
           no_s, Enspin, &
           maxnh, qtot, &
           temp, na_u, ef, &  
           ucell, ts_kscell, ts_kdispl, &      
           ts_gamma_scf, istep, ia1)


!     INPUT:
!     filename
      character(*), intent(in) ::   fname

!     OUTPUT:
!     Output full matrix
      logical, intent(out) :: gamma, onlyS, ts_gamma_scf
      integer, intent(out) :: no_u, Enspin, maxnh
      integer, intent(out) :: istep, ia1, na_u, no_s
      real*8, intent(out) :: qtot, ef, temp, ucell(3,3)          
      integer, intent(out) :: ts_kscell(3,3)
      real*8, intent(out) :: ts_kdispl(3)

!     Local variables:
      integer :: iu, ih, is, k

      iu=12334
      open( iu, file=fname, form='unformatted', status='old' )      

      read(iu) na_u, no_u, no_s, Enspin, maxnh

!     Deallocate memmory
      if (allocated(S)) then 
        deallocate(S, H, xij, xa)
	deallocate(indxuo, listhptr, numh, listh, lasto, isa)        
      end if 
!     Reallocate it
      allocate(xa(3,na_u))
      allocate(isa(na_u)) 
      allocate(lasto(0:na_u))
      allocate(indxuo(1:no_s))
      allocate(numh(no_u))
      allocate(listhptr(no_u))
      allocate(listh(maxnh))
      allocate(S(maxnh))
      allocate(H(maxnh,Enspin))
      allocate(xij(3,maxnh))

      read(iu) xa
      read(iu) isa   
      read(iu) ucell  
      read(iu) gamma  
      read(iu) onlyS
      read(iu) ts_gamma_scf       
      read(iu) ts_kscell
      read(iu) ts_kdispl 
      read(iu) istep, ia1
      read(iu) lasto
      
      indxuo(:)=-1
      do ih=1,no_u
         indxuo(ih)=ih
      enddo
      if (.not. gamma) then
         read(iu) (indxuo(ih),ih=1,no_s)
      endif

      read(iu) numh(1:no_u)
      read(iu) qtot,temp
      read(iu) ef
      
      listhptr(1) = 0
      do ih = 2,no_u
         listhptr(ih) = listhptr(ih-1) + numh(ih-1)
      enddo
      
      do ih = 1,no_u
         read(iu) listh(listhptr(ih)+1:listhptr(ih)+numh(ih))
      enddo
      
      do ih = 1,no_u
         read(iu) S(listhptr(ih)+1:listhptr(ih)+numh(ih))
      enddo
      
      if (.not.onlyS) then
         do is = 1,Enspin
            do ih = 1,no_u
               read(iu) H(listhptr(ih)+1:listhptr(ih)+numh(ih),is)
            enddo
         enddo
      endif

      xij=0.0
      if (.not.gamma) then
         do ih = 1,no_u
            read(iu) (xij(k,listhptr(ih)+1:listhptr(ih)+numh(ih)),k=1,3)
         enddo
      end if

      close(iu)

      end subroutine read

      end module readnewtshs
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------


