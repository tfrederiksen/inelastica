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
module readTSHS

  implicit none

  public

  real(kind=8), dimension(:), allocatable :: S
  real(kind=8), dimension(:,:), allocatable :: H, xij, xa
  integer, dimension(:), allocatable :: listhptr, numh, listh, lasto

CONTAINS

  subroutine deallocate()
    if (allocated(S)) then
       deallocate(S, xij, xa)
       deallocate(listhptr, numh, listh, lasto) 
    end if
    if (allocated(H)) deallocate(H)
  end subroutine deallocate

  subroutine read(fname, gamma, onlyS, no_u, &
       no_s, nspin, &
       maxnh, qtot, &
       Temp, na_u, Ef, & 
       ucell, ts_kscell, ts_kdispl, & 
       ts_gamma_scf, istep, ia1)

  ! We do direct conversion instead of doing it in Inelastica
    real(kind=8), parameter :: Bohr2Ang = 0.529177_8
    real(kind=8), parameter :: Ry2eV = 13.6058_8

! INPUT:
! filename
    character(len=*), intent(in) :: fname
    
! OUTPUT:
! Output full matrix
    logical, intent(out) :: gamma, onlyS, ts_gamma_scf
    integer, intent(out) :: no_u, nspin, maxnh
    integer, intent(out) :: istep, ia1, na_u, no_s
    real(kind=8), intent(out) :: qtot, ef, temp, ucell(3,3)
    integer, intent(out) :: ts_kscell(3,3)
    real(kind=8), intent(out) :: ts_kdispl(3)

!     Local variables:
    integer :: iu, is, ind, io, i, j, ia, ja
    integer :: version, tm(3)
    logical :: exist
    integer :: n_s
    integer, allocatable :: offsets(:,:), indxs(:)

    ! Get file version
    version = tshs_version(fname)

    select case ( version )
    case ( 0 , 1 )
       ! do nothing
    case default
       stop 'Unsupported TSHS version file'
    end select
    inquire(file=fname,exist=exist)
    if ( .not. exist ) then
       stop 'Error could not read file given'
    end if

    ! Open file
    iu = 12334
    open( iu, file=fname, form='unformatted', status='old' )
       
    ! Read Dimensions Information
    if ( version /= 0 ) then
       read(iu) version
    end if
    read(iu) na_u, no_u, no_s, nspin, maxnh

    call deallocate()
    
    ! Read Geometry information
    allocate(xa(3,na_u)) 
    if ( version == 0 ) then
       read(iu) xa
       read(iu) ! iza
       read(iu) ucell
    else if ( version == 1 ) then
       read(iu) ucell, xa
    end if
    xa = xa * Bohr2Ang
    ucell = ucell * Bohr2Ang

    ! Read k-point sampling information
    if ( version == 0 ) then
       read(iu) Gamma
       read(iu) onlyS
       read(iu) ts_gamma_scf
       read(iu) ts_kscell
       read(iu) ts_kdispl
    else if ( version == 1 ) then
       read(iu) Gamma, ts_gamma_scf, onlyS
       read(iu) ts_kscell, ts_kdispl
       read(iu) Ef, Qtot, Temp
    end if
    read(iu) istep, ia1

    allocate(lasto(0:na_u))
    read(iu) lasto

    if ( .not. Gamma ) then
       read(iu) ! indxuo
    endif

    allocate(numh(no_u))
    read(iu) numh

    allocate(listhptr(no_u))
    listhptr(1) = 0
    do i = 2 , no_u
       listhptr(i) = listhptr(i-1) + numh(i-1)
    end do
    
    if ( version == 0 ) then
       ! Read Electronic Structure Information
       read(iu) Qtot,Temp
       read(iu) Ef
    end if

    Ef   = Ef   * Ry2eV
    Temp = Temp * Ry2eV

    ! Read listh
    allocate(listh(maxnh))
    do i = 1 , no_u
       read(iu) listh(listhptr(i)+1:listhptr(i)+numh(i))
    end do
    
    ! Read Overlap matrix
    allocate(S(maxnh))
    do i = 1 , no_u
       read(iu) S(listhptr(i)+1:listhptr(i)+numh(i))
    end do
    
    if ( .not. onlyS ) then
       allocate(H(maxnh,nspin))
       do is = 1 , nspin 
          do i = 1 , no_u
             read(iu) H(listhptr(i)+1:listhptr(i)+numh(i),is)
             H(listhptr(i)+1:listhptr(i)+numh(i),is) = &
                  H(listhptr(i)+1:listhptr(i)+numh(i),is) * Ry2eV
          end do
       end do
    end if  ! onlyS

    if ( .not. Gamma ) then
       
       allocate(xij(3,maxnh))

       if ( version == 0 ) then
          do i = 1 , no_u
             read(iu) (xij(j,listhptr(i)+1:listhptr(i)+numh(i)),j=1,3)
             xij(:,listhptr(i)+1:listhptr(i)+numh(i)) = &
                  xij(:,listhptr(i)+1:listhptr(i)+numh(i)) * Bohr2Ang
          end do
       else if ( version == 1 ) then

          allocate(indxs(maxnh))
          read(iu) n_s, indxs
          allocate(offsets(3,n_s))
          read(iu) offsets
          
          do io = 1 , no_u
             ia = iaorb(io,lasto)
             do j = 1 , numh(io)
                ind = listhptr(io) + j
                ja = iaorb(ucorb(listh(ind),no_u),lasto)
                
                tm(:) = offsets(:,indxs(ind))
                xij(:,ind) = ucell(:,1) * tm(1) &
                     + ucell(:,2) * tm(2) &
                     + ucell(:,3) * tm(3) &
                     + xa(:,ja) - xa(:,ia)
                
             end do
          end do
          
          ! clean-up
          deallocate(offsets,indxs)
          
       end if
    end if
       
    close(iu)

  contains
    
    function TSHS_version(fname) result(version)
      character(len=*), intent(in) :: fname
      integer :: version
      integer :: iu
      integer :: na_u, no_u, no_s, nspin, maxnh, err

      iu = 21356
      open( iu, file=fname, form='unformatted', status='unknown' )

      read(iu,iostat=err) na_u, no_u, no_s, nspin, maxnh
      if ( err == 0 ) then
         ! we can successfully read 5 integers
         version = 0
      else
         backspace(iu)
         read(iu,iostat=err) version
      end if

      close(iu)

    end function TSHS_version


  function iaorb(iorb,lasto) result(ia)
! *********************
! * INPUT variables   *
! *********************
    integer, intent(in) :: iorb
    integer, intent(in) :: lasto(:) ! actually lasto(0:na_u) !
! *********************
! * OUTPUT variables  *
! *********************
    integer :: ia, iiorb, na

    ! retrieve number of atoms
    na = ubound(lasto,dim=1)-1

    ! We need to transfer to the unit-cell
    iiorb = ucorb(iorb,lasto(na+1))

    ! Start searching in the "best" guess of the placement
    ! Our best guess is the number of orbitals on the first atom
    ! 1+int(iorb / orb(1)) == atom (if all atoms have the same orb. count, this is exact)
    ! This is the most qualified quess we can make
    ! It makes no sense to take the mean of the orbital count...!
    ! If iiorb and lasto are both > 0 then this will always
    ! return in the range [1;na]
    ia = max(1,min(int(real(iiorb,8) / lasto(2)),na))

    ia_loop: do
       if ( iiorb < lasto(ia) ) then
          ! 1. We have overestimated the orbital position
          ia = ia - 1
          cycle ia_loop
       else if ( iiorb > lasto(ia+1) ) then
          ! 2. We have overestimated the orbital position
          ia = ia + 1
          cycle ia_loop
       else if ( iiorb == lasto(ia) ) then
          ! 3. it is on the former atom
          ia = ia - 1
       end if
       ! We have found it!!!
       return
    end do ia_loop

    stop 'SOMETHING WENT TERRIBLY WRONG IN IAORB'
  end function iaorb

  elemental function ucorb(a,p) result(MODP)
    integer, intent(in) :: a,p
    integer :: MODP
    if ( a > p ) then
       MODP = MOD(a,p)
       if ( MODP == 0 ) MODP = p
    else
       MODP = a
    end if
  end function UCORB

end subroutine read
  
end module readTSHS
