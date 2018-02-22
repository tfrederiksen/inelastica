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

  integer :: version = 0
  real(kind=8), allocatable :: S(:)
  real(kind=8), allocatable :: H(:,:)
  real(kind=8), allocatable :: xij(:,:)
  real(kind=8), allocatable :: xa(:,:)
  integer, allocatable :: listhptr(:)
  integer, allocatable :: numh(:)
  integer, allocatable :: listh(:)
  integer, allocatable :: lasto(:)

CONTAINS

  subroutine deallocate()
    if ( allocated(S) ) then
       deallocate(S, xij, xa)
       deallocate(listhptr, numh, listh, lasto) 
       if ( allocated(H) ) then
          deallocate(H)
       end if
    end if
  end subroutine deallocate

  subroutine read_size(fname, no_u, no_s, nspin, na_u,maxnh)
    character(len=*), intent(in) :: fname
    integer, intent(out) :: no_u, nspin
    integer, intent(out) ::  na_u, no_s, maxnh
    integer :: iu, version
    logical :: exist

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

    close(iu)

  end subroutine read_size


  subroutine read(fname, gamma, onlyS, no_u, &
       no_s, nspin, &
       maxnh, qtot, &
       Temp, na_u, Ef, & 
       ucell, ts_kscell, ts_kdispl, & 
       ts_gamma_scf, istep, ia1)

  ! We do direct conversion instead of doing it in Inelastica
    real(kind=8), parameter :: Bohr2Ang = 0.529177d0
    real(kind=8), parameter :: Ry2eV = 13.6058d0

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
    integer :: tm(3)
    logical :: exist
    integer :: n_s
    integer, allocatable :: offsets(:,:)

    ! Get file version (global variable)
    version = tshs_version(fname)

    select case ( version )
    case ( 0 , 1 )
       ! do nothing
    case default
       stop 'Unsupported TSHS version file, convert to 0,1'
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
       read(iu) ! nsc
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

    if ( version == 0 .and. .not. Gamma ) then
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

    allocate(xij(3,maxnh))

    if ( Gamma ) then

       xij(:,:) = 0._8

    else

       if ( version == 0 ) then
          do i = 1 , no_u
             read(iu) (xij(j,listhptr(i)+1:listhptr(i)+numh(i)),j=1,3)
             xij(:,listhptr(i)+1:listhptr(i)+numh(i)) = &
                  xij(:,listhptr(i)+1:listhptr(i)+numh(i)) * Bohr2Ang
          end do
       else if ( version == 1 ) then
          n_s = no_s / no_u
          allocate(offsets(3,0:n_s-1))
          read(iu) offsets
          do io = 1 , no_u
             ia = iaorb(io,lasto)
             do j = 1 , numh(io)
                ind = listhptr(io) + j
                ja = iaorb(ucorb(listh(ind),no_u),lasto)
                
                ! supercell index
                is = (listh(ind)-1)/no_u
                tm(:) = offsets(:,is)
                xij(:,ind) = ucell(:,1) * tm(1) &
                     + ucell(:,2) * tm(2) &
                     + ucell(:,3) * tm(3)
                
             end do
          end do
          
          ! clean-up
          deallocate(offsets)
          
       end if
    end if
       
    close(iu)

  contains

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
    
  end subroutine read

  subroutine remove_Atoms(list,na_u,no_u,maxnh)
    integer, intent(in) :: list(:)
    integer, intent(out) :: na_u, no_u, maxnh

    real(kind=8), allocatable :: nS(:)
    real(kind=8), allocatable :: nH(:,:)
    real(kind=8), allocatable :: nxij(:,:)
    real(kind=8), allocatable :: nxa(:,:)
    integer, allocatable :: nlisthptr(:)
    integer, allocatable :: nnumh(:)
    integer, allocatable :: nlisth(:)
    integer, allocatable :: nlasto(:)
    
    integer :: ona, ono, na
    integer :: i, ind, j, io, no, ia, jo, a, nind, is, diffo
    logical :: has_H
    
    if ( .not. allocated(S) ) stop 'Your arrays have not been populated!'

    na = size(list)
    ! get sizes
    ona = size(lasto) - 1
    ono = lasto(ona)
    na_u = ona - na
    if ( ona == na_u ) then
       na_u = ona
       no_u = ono
       maxnh = size(S)
       return
    end if
    if ( any(list > ona) .or. any( 1 > list) ) then
       stop 'Hello!!!! Your are removing a non existing atom?'
    end if

    allocate(nxa(3,na_u),nlasto(0:na_u))
    nlasto(0) = 0 

    j = 0
    do i = 1 , ona
       if ( any(i==list) ) cycle
       j = j + 1
       nxa(:,j) = xa(:,i)
       nlasto(j) = nlasto(j-1) + lasto(i) - lasto(i-1)
    end do
    if ( j /= na_u ) stop 'Hello!!! Error in number of atoms'
    no_u = nlasto(na_u)

    ! Calculate difference in orbitals
    diffo = ono - no_u

    allocate(nnumh(no_u),nlisthptr(no_u))
    nlisthptr(1) = 0

    j = 0
    do io = 1 , ono
       ! Count number of orbitals not in the buffer atoms
       if ( remove(io,ona,lasto,na,list) ) cycle
       no = 0
       do ind = listhptr(io) + 1 , listhptr(io) + numh(io)
          i = ucorb(listh(ind),ono)
          if ( remove(i,ona,lasto,na,list) ) cycle
          no = no + 1
       end do
       ! we have now counted number of orbitals in the
       ! io'th row
       ! Save the data
       j = j + 1
       nnumh(j) = no
       if ( j < no_u ) then
          nlisthptr(j+1) = nlisthptr(j) + no
       end if

    end do
    if ( j /= no_u ) stop 'Hello!!! Error in number of rows'

    ! Number of non-zero elements
    maxnh = nlisthptr(no_u) + nnumh(no_u)

    ! allocate sparsity patterns
    allocate(nS(maxnh),nxij(3,maxnh),nlisth(maxnh))
    has_H = allocated(H)
    if ( has_H ) then
       allocate(nH(maxnh,size(H,2)))
    end if

    ! we have now counted all elements in the new reduced sparsity pattern
    nind = 0
    jo = 0
    do io = 1 , ono

       ! Count number of orbitals not in the buffer atoms
       if ( remove(io,ona,lasto,na,list) ) cycle
       jo = jo + 1
       
       ! copy over data
       do ind = listhptr(io) + 1 , listhptr(io) + numh(io)
          i = ucorb(listh(ind),ono)
          if ( remove(i,ona,lasto,na,list) ) cycle

          ! Update index
          nind = nind + 1 

          ! supercell index, we need to correct for the supercell
          ! offset
          is = (listh(ind)-1)/ono

          nlisth(nind) = listh(ind) - diffo * is ! copy over column
          do a = 1 , na
             ia = list(a)
             if ( i > lasto(ia) ) then
                ! We need to correct for this atom
                nlisth(nind) = nlisth(nind) - lasto(ia) + lasto(ia-1)
             end if
          end do

          nxij(:,nind) = xij(:,ind)
          nS(nind) = S(ind)
          if ( has_H ) then
             nH(nind,:) = H(ind,:)
          end if

       end do

       if ( nind /= nlisthptr(jo) + nnumh(jo) ) then
          stop 'Error in sparsity creation'
       end if

    end do

    if ( nind /= maxnh ) then
       stop 'Error in sparsity creation, not all traversed.'
    end if

    ! deallocate everything
    call deallocate()

    allocate(S(maxnh),xij(3,maxnh),xa(3,na_u))
    allocate(listhptr(no_u),numh(no_u),listh(maxnh),lasto(0:na_u))
    
    ! point to new data
    S = nS
    xij = nxij
    xa = nxa
    listhptr = nlisthptr
    numh = nnumh
    listh = nlisth
    lasto = nlasto
    if ( has_H ) then
       allocate(H(maxnh,size(nH,2)))
       H = nH
       deallocate(nH)
    end if
    deallocate(nS,nxij,nxa,nlisthptr,nnumh,nlisth,nlasto)

  contains
    
    function remove(io,na_u,lasto,na,list) result(rem)
      integer, intent(in) :: io,na_u, lasto(0:na_u),na,list(na)
      logical :: rem
      integer :: i, ia
      do i = 1 , na
         ia = list(i)
         if ( lasto(ia-1) < io .and. io <= lasto(ia) ) then
            rem = .true. 
            return
         end if
      end do
      rem = .false.
    end function remove
    
  end subroutine remove_Atoms

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

  function ucorb(a,p) result(MODP)
    integer, intent(in) :: a,p
    integer :: MODP
    if ( a > p ) then
       MODP = MOD(a,p)
       if ( MODP == 0 ) MODP = p
    else
       MODP = a
    end if
  end function UCORB
  
end module readTSHS
