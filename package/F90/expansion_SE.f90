subroutine expansion_SE(no_u,no_s,NA1,NA2,kpt, &
     na_u,lasto,ESH,G0,oESH,SGF)

  implicit none
  
  real(kind=8), parameter :: Pi = 3.141592653589793238d0
! ********************
! * INPUT variables  *
! ********************
  integer,  intent(in) :: no_u, no_s, NA1, NA2
  real(kind=8), intent(in) :: kpt(2)
  integer,  intent(in) :: na_u, lasto(0:na_u)
  complex(kind=8), dimension(no_u,no_u,NA1*NA2), intent(in) :: ESH, G0
! ********************
! * OUTPUT variables *
! ********************
  complex(kind=8), dimension(no_s,no_s), intent(out) :: oESH, SGF
! ********************
! * LOCAL variables  *
! ********************
  integer :: iq
  integer :: iow,iau,ia2,ia1,iuo
  integer :: jow,jau,ja2,ja1,juo
  complex(kind=8) :: ph
  real(kind=8) :: qmPi(2,NA1*NA2), wq

  ! We know that NA1*NA2 > 1

  iq = 0
  do ia2 = 0 , NA2 - 1
  do ia1 = 0 , NA1 - 1
     iq = iq + 1
     qmPi(1,iq) = - 2._8 * Pi * (kpt(1)+ia1)/NA1
     qmPi(2,iq) = - 2._8 * Pi * (kpt(2)+ia2)/NA2
  end do
  end do

  wq = 1._8 / real(NA1*NA2,8)

! This is the crucial calcuation.
! If we use bulk values in the electrodes
! we need not add the expanded H and S values to get the 
! electrode \Sigma. Hence, we need only expand
! surface Green's function
  iq = 1
  iow = 0
  do iau = 1 , na_u
   do ia2 = 1 , NA2
   do ia1 = 1 , NA1
     do iuo = 1 + lasto(iau-1) , lasto(iau)
      iow = iow + 1
      jow = 0
      do jau = 1 , na_u
       do ja2 = 1 , NA2
       do ja1 = 1 , NA1
         ph = wq * cdexp(dcmplx(0._8, &
              (ia1-ja1)*qmPi(1,iq) + (ia2-ja2)*qmPi(2,iq) ) )
         do juo = 1 + lasto(jau-1) , lasto(jau)
           jow = jow + 1
           
           SGF(jow,iow)  = ph * G0(juo,iuo,iq)
           
           oESH(jow,iow) = ph * ESH(juo,iuo,iq)
           
         end do !juo
       end do !ja1
       end do !ja2
      end do !jau
     end do !iuo
   end do !ia1
   end do !ia2
  end do !iau
  do iq = 2 , NA1*NA2
   iow = 0
   do iau = 1 , na_u
    do ia2 = 1 , NA2
    do ia1 = 1 , NA1
      do iuo = 1 + lasto(iau-1) , lasto(iau)
       iow = iow + 1
       jow = 0
       do jau = 1 , na_u
        do ja2 = 1 , NA2
        do ja1 = 1 , NA1
          ph = wq * cdexp(dcmplx(0._8, &
              (ia1-ja1)*qmPi(1,iq) + (ia2-ja2)*qmPi(2,iq) ) )
          do juo = 1 + lasto(jau-1) , lasto(jau)
             jow = jow + 1
             
             SGF(jow,iow)  = SGF(jow,iow) + ph * G0(juo,iuo,iq)
   
             oESH(jow,iow) = oESH(jow,iow) + ph * ESH(juo,iuo,iq)
   
          end do !juo
        end do !ja1
        end do !ja2
       end do !jau
      end do !iuo
    end do !ia1
    end do !ia2
   end do !iau
  end do !q-points

end subroutine expansion_SE
