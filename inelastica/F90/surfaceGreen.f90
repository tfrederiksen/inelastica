  ! Calculates the surface Green's function for the electrodes
  ! Handles both the left and right one
  subroutine surfaceGreen(no,ZE,H00,S00,H01,S01,accur,is_left,GS)

! ***************** INPUT **********************************************
! integer     no      : Number of orbitals in the electrode
! complex(dp) ZE      : The energy of the Green's function evaluation
! complex(dp) H00     : Hamiltonian within the first unit cell (discarding T-direction)
! complex(dp) S00     : Overlap matrix within the first unit cell (discarding T-direction)
! complex(dp) H01     : Transfer matrix from H00 to the neighbouring cell (in T-direction)
! complex(dp) S01     : Transfer matrix from S00 to the neighbouring cell (in T-direction)
! ***************** OUTPUT *********************************************
! complex(dp) GS      : Surface Green's function of the electrode
! **********************************************************************
    implicit none

    complex(kind=8), parameter :: z_1  = dcmplx(1._8,0._8)
    complex(kind=8), parameter :: z_m1 = dcmplx(-1._8,0._8)
    complex(kind=8), parameter :: z_0  = dcmplx(0._8,0._8)

! ***********************
! * INPUT variables     *
! ***********************
    integer,     intent(in) :: no
    complex(kind=8), intent(in) :: ZE 
    complex(kind=8), intent(in) :: H00(no*no),S00(no*no)
    complex(kind=8), intent(in) :: H01(no*no),S01(no*no)

    real(kind=8)   , intent(in) :: accur
    logical,     intent(in) :: is_left

! ***********************
! * OUTPUT variables    *
! ***********************
    complex(kind=8), intent(out) :: GS(no*no)

! ***********************
! * LOCAL variables     *
! ***********************
    complex(kind=8), allocatable, target :: zwork(:)
    integer :: nom1, no2, nosq
    integer :: ierr             !error in inversion
    integer :: i,j,ic,ic2
    logical :: LEFT, as_first

    real(kind=8) :: ro

    ! on the stack...
    integer :: ipvt(no)
    complex(kind=8), dimension(:), pointer :: rh,rh1,w,alpha,beta,gb

    ! The inelastica way is opposite, hence we have to reverse it here.
    LEFT = .not. is_left

    nom1 = no - 1
    no2  = 2 * no
    nosq = no * no

    allocate(zwork(8*nosq))

    i = 0
    rh  => zwork(i+1:i+2*nosq) 
    i = i + 2*nosq
    rh1 => zwork(i+1:i+2*nosq) 
    i = i + 2*nosq
    alpha => zwork(i+1:i+nosq) 
    i = i + nosq
    beta => zwork(i+1:i+nosq) 
    i = i + nosq
    w => zwork(i+1:i+nosq)
    i = i + nosq
    gb => zwork(i+1:i+nosq) 

! gb    =   Z*S00-H00
! alpha = -(Z*S01-H01)
! gs  = Z*S00-H00
    do i = 1 , nosq
       gb(i)    = ZE * S00(i) - H00(i)
       GS(i)    = gb(i)
       alpha(i) = H01(i) - ZE * S01(i)
    end do

! beta = -(Z*S10-H10)
    do j = 1 , no
       ic = no * (j-1)
       do i = 1 , no
          ic2 = j + no*(i-1)
          beta(ic+i) = dconjg(H01(ic2)) - ZE * dconjg(S01(ic2))
       end do
    end do

    ! Initialize loop
    ro = accur + 1._8
    as_first = .false.
    do while ( ro > accur ) 

! rh = -(Z*S01-H01) ,j<no
! rh = -(Z*S10-H10) ,j>no
       do i = 1, nosq
          rh(i)       = alpha(i)
          rh(nosq+i)  = beta(i)
       end do

! w = Z*S00-H00
       w(:) = gb(:)

! rh =  rh1^(-1)*rh
! rh =  t0
       call zgesv(no, no2, w, no, ipvt, rh, no, ierr)

       if ( ierr /= 0 ) then
          write(*,*) 'ERROR: calc_green 1 MATRIX INVERSION FAILED'
          write(*,*) 'ERROR: LAPACK INFO = ',ierr
       end if

       ! switch pointers instead of copying elements
       call switch_alpha_beta_rh1(as_first)

! alpha = -(Z*S01-H01)*t0
       call zgemm('N','N',no,no,no,z_1,rh1(1),no,rh(1),no,z_0,alpha,no)
! beta  = -(Z*S10-H10)*t0 ??
       call zgemm('N','N',no,no,no,z_1,rh1(nosq+1),no,rh(nosq+1),no,z_0,beta,no)

       if ( LEFT ) then
! ba    = (Z*S10-H10)*t0b
          call zgemm('N','N',no,no,no,z_m1,rh1(nosq+1),no,rh(1),no,z_0,w,no)
          gb(:) = gb(:) + w(:)
          gs(:) = gs(:) + w(:)
       else
! gb = gb + [ba    = (Z*S10-H10)*t0b]
          call zgemm('N','N',no,no,no,z_m1,rh1(nosq+1),no,rh(1),no,z_1,gb,no)
       end if

! ab    = (Z*S01-H01)*t0
       call zgemm('N','N',no,no,no,z_m1,rh1(1),no,rh(nosq+1),no,z_0,w,no)
       
       gb(:) = gb(:) + w(:)
       if ( .not. LEFT ) then
          gs(:) = gs(:) + w(:)
       end if
       
       ro = -1._8
       do i = 1 , nosq
          ro = max(ro,abs(w(i)))
       end do

    end do

    call zgetrf(no, no, GS, no, ipvt, ierr )
    if ( ierr /= 0 ) print *,"WARNING: Error in LU-decomposition"
    call zgetri(no, GS, no, ipvt, zwork, nosq, ierr)
    if ( ierr /= 0 ) then
       write(*,*) 'ERROR: surfaceGreen GS MATRIX INVERSION FAILED'
       write(*,*) 'ERROR: LAPACK INFO = ',ierr
    end if

    deallocate(zwork)

  contains

    ! We supply a routine to switch the pointer position of alpha,beta / rh1
    subroutine switch_alpha_beta_rh1(as_first)
      logical, intent(inout) :: as_first
      integer :: i 
      ! start
      i = 2 * nosq

      if ( as_first ) then
         rh1 => zwork(i+1:i+2*nosq) 
         i = i + 2*nosq
         alpha => zwork(i+1:i+nosq) 
         i = i + nosq
         beta => zwork(i+1:i+nosq) 
      else
         alpha => zwork(i+1:i+nosq) 
         i = i + nosq
         beta => zwork(i+1:i+nosq) 
         i = i + nosq
         rh1 => zwork(i+1:i+2*nosq) 
      end if
      as_first = .not. as_first

    end subroutine switch_alpha_beta_rh1

  end subroutine surfaceGreen
