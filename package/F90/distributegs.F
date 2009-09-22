!
!     Fortran version of part of getSig
!
!      ESmH, SGF = f90distributegs(loop, nuo, nua, NA1, NA2, kpoint, matEsmH, g0, ESmH, SGF)
!

      subroutine f90distributegs(loop, nuo, nua, NA1, NA2, kpoint, 
     +     matEsmH, g0, ESmH, SGF, outESmH, outSGF)

      implicit none

!     INPUT:
!     Loop indexes
      integer, intent(in) ::  loop(0:nua*NA1*NA2-1,0:6)        
!     Dimensions
      integer, intent(in) ::  nuo, nua, NA1, NA2                    
!     K-point
      real*8, intent(in) ::   kpoint(0:1)              
!     matESmH 
      complex*16, intent(in) ::   matESmH(0:nuo-1,0:nuo-1)                      
!     g0 
      complex*16, intent(in) ::   g0(0:nuo-1,0:nuo-1)                       
!     Matrices:
      complex*16, intent(in) :: ESmH(0:nuo*NA1*NA2-1,0:nuo*NA1*NA2-1)   
      complex*16, intent(in) :: SGF(0:nuo*NA1*NA2-1,0:nuo*NA1*NA2-1)   

!     OUTPUT:
      complex*16,intent(out) :: outESmH(0:nuo*NA1*NA2-1,0:nuo*NA1*NA2-1)   
      complex*16,intent(out) :: outSGF(0:nuo*NA1*NA2-1,0:nuo*NA1*NA2-1)   

!     =======================================================================
!     Local variables 

      complex*16   phase
!     Loop indecies
      integer      ii, ia, i1, i2, iSGFs, iSGFe, ig0s, ig0e
      integer      jj, ja, j1, j2, jSGFs, jSGFe, jg0s, jg0e

!      for ia, i1, i2, iSGFs, iSGFe, ig0s, ig0e in loop:
!          for ja, j1, j2, jSGFs, jSGFe, jg0s, jg0e in loop:


      outESmH = ESmH
      outSGF = SGF

      do ii=0,nua*NA1*NA2-1
         ia=   loop(ii,0)
         i1=   loop(ii,1)
         i2=   loop(ii,2)
         iSGFs=loop(ii,3)
         iSGFe=loop(ii,4)
         ig0s= loop(ii,5)
         ig0e= loop(ii,6)

         do jj=0,nua*NA1*NA2-1
            ja=   loop(jj,0)
            j1=   loop(jj,1)
            j2=   loop(jj,2)
            jSGFs=loop(jj,3)
            jSGFe=loop(jj,4)
            jg0s= loop(jj,5)
            jg0e= loop(jj,6)

!              # Same convention as for H_ij above: exp(2 pi i k * (jatom-iatom))
!              # The phases etc have been checked by comparing the self-energy from 
!              # 1x1 and 3x3 electrode calculations
!              phase = 1.0/(NA1*NA2)*N.exp(2.0j*N.pi*((j1-i1)*kpoint[0]+(j2-i2)*kpoint[1])) 
!              ESmH[iSGFs:iSGFe+1,jSGFs:jSGFe+1]=ESmH[iSGFs:iSGFe+1,jSGFs:jSGFe+1]+\
!                N.conjugate(phase)*matESmH[ig0s:ig0e+1,jg0s:jg0e+1]       
!              SGF[iSGFs:iSGFe+1,jSGFs:jSGFe+1]=SGF[iSGFs:iSGFe+1,jSGFs:jSGFe+1]+\
!                N.conjugate(phase)*g0[ig0s:ig0e+1,jg0s:jg0e+1]
            phase = 1.0d0/(NA1*NA2)*cdexp(
     +           dcmplx(0d0,-6.283185307179586477d0)*
     +           ((j1-i1)*kpoint(0)+(j2-i2)*kpoint(1)))           ! Note changed sign
            ! Here fortran is different in the second index by one 
            outESmH(iSGFs:iSGFe,jSGFs:jSGFe) = 
     +           outESmH(iSGFs:iSGFe,jSGFs:jSGFe)+
     +           phase*matESmH(ig0s:ig0e,jg0s:jg0e)        
            outSGF(iSGFs:iSGFe,jSGFs:jSGFe)=
     +           outSGF(iSGFs:iSGFe,jSGFs:jSGFe)+
     +           phase*g0(ig0s:ig0e,jg0s:jg0e)
            
         enddo
      enddo

      return 
      end 
