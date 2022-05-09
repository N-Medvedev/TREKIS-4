! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2021
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with simple potentials (fitting into one subroutine),
! such as power-law potential (convenient to use for combined potentials like C/r^n),
! Lennard-Jones, exponential, etc.
! [1] https://en.wikipedia.org/wiki/Lennard-Jones_potential

MODULE MD_Pot_Simple
!use Universal_constants
use Little_subroutines, only: fast_pow


implicit none


 contains



function exp_potential(C, k, r) result(Pot) ! exponential potential
   real(8) Pot              ! [eV]
   real(8), intent(in) :: C ! [eV]
   real(8), intent(in) :: k ! [1/A]
   real(8), intent(in) :: r ! [A] interatomic distance

   Pot = C * exp(-k*r)
end function exp_potential


function d_exp_potential(C, k, r) result (dPot) ! derivative of exponential potential
   real(8) dPot             ! [eV/A]
   real(8), intent(in) :: C ! [eV]
   real(8), intent(in) :: k ! [1/A]
   real(8), intent(in) :: r ! [A] interatomic distance

   dPot = -C*k * exp(-k*r)
end function d_exp_potential



function LJ_potential(C12, C6, r) ! Lennard-Jones potential (in AB form [1])
   real(8) :: LJ_potential  ! [eV]
   real(8), intent(in) :: C12, C6   ! [eV*r^6], [eV*r^12] LJ coefficients
   real(8), intent(in) :: r ! [A] interatomic distance
   real(8) :: a_r6, a_r12
   a_r6 = fast_pow(r,6) ! function from module "Little_subroutines"
   a_r12 = a_r6*a_r6
   LJ_potential = C12/a_r12 - C6/a_r6
end function LJ_potential


function d_LJ_potential(C12, C6, r) ! derivative of Lennard-Jones potential
   real(8) :: d_LJ_potential    ! [eV/A]
   real(8), intent(in) :: C12, C6   ! [eV*r^6], [eV*r^12] LJ coefficients
   real(8), intent(in) :: r ! [A] interatomic distance
   real(8) :: a_r6, a_r7, a_r13
   a_r6 = fast_pow(r,6) ! function from module "Little_subroutines"
   a_r7 = a_r6*r
   a_r13 = a_r6*a_r7
   d_LJ_potential = -12.0d0*C12/a_r13 + 6.0d0*C6/a_r7
end function d_LJ_potential



function power_potential(C, n, r) result(Pot) ! Potential C*r^n
   real(8) Pot              ! [eV]
   real(8), intent(in) :: C ! [eV*r^(n-1)]
   real(8), intent(in) :: n ! power
   real(8), intent(in) :: r ! [A] interatomic distance

   Pot = C * r**n
end function power_potential


function d_power_potential(C, n, r) result (dPot) ! derivative of power potential C*r^n
   real(8) dPot             ! [eV/A]
   real(8), intent(in) :: C ! [eV*r^(n-1)]
   real(8), intent(in) :: n ! power
   real(8), intent(in) :: r ! [A] interatomic distance

   dPot = C*n * r**(n-1.0d0)
end function d_power_potential



END MODULE MD_Pot_Simple
