! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2021
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with Buckingham [1] and Matsui [2] potential
! [1] https://en.wikipedia.org/wiki/Buckingham_potential
! [2] Matsui, Phys Chem Minerals 16, 234 (1988) 


MODULE MD_Pot_Buck
!use Universal_constants
use Little_subroutines, only: fast_pow

implicit none

real(8) :: m_Matsui_f


parameter(m_Matsui_f = 0.04336410066)    ! "a standard force" [2] [eV/A*atom]


 contains



function Matsui_potential(A, B, C, r) result(Pot) ! Matsui potential [2]
   real(8) :: Pot  ! [eV]
   real(8), intent(in) :: A, B, C
   real(8), intent(in) :: r ! [A] interatomic distance
   real(8) :: r6
   r6 = fast_pow(r,6)   ! function from module "Little_subroutines"
   Pot = m_Matsui_f*B*exp((A-r)/B) - C/r6
end function Matsui_potential


function d_Matsui_potential(A, B, C, r) result(dPot)  ! derivative by r of Matsui potential
   real(8) :: dPot    ! [eV/A]
   real(8), intent(in) :: A, B, C
   real(8), intent(in) :: r ! [A] interatomic distance
   real(8) :: r6, r7
   r6 = fast_pow(r,6) ! function from module "Little_subroutines"
   r7 = r6*r
   dPot = -m_Matsui_f*exp((A-r)/B) + 6.0d0*C/r7
end function d_Matsui_potential
 


function Buck_potential(A, B, C, r) result(Pot) ! Buckingham potential [1]
   real(8) :: Pot  ! [eV]
   real(8), intent(in) :: A, B, C
   real(8), intent(in) :: r ! [A] interatomic distance
   real(8) :: r6
   r6 = fast_pow(r,6)   ! function from module "Little_subroutines"
   Pot = A*exp(-B*r) - C/r6
end function Buck_potential


function d_Buck_potential(A, B, C, r) result(dPot)  ! derivative by r of Buckingham potential
   real(8) :: dPot    ! [eV/A]
   real(8), intent(in) :: A, B, C
   real(8), intent(in) :: r ! [A] interatomic distance
   real(8) :: r6, r7
   r6 = fast_pow(r,6) ! function from module "Little_subroutines"
   r7 = r6*r
   dPot = -A*B*exp(-B*r) + 6.0d0*C/r7
end function d_Buck_potential


END MODULE MD_Pot_Buck
