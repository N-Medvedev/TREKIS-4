! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by F. Akhmetov
! in 2021
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with Mermin-type CDF
! [1] N. Medvedev, et al., J. Phys. D: Appl. Phys. 48 (2015) 355303
! [2] R. Rymzhanov et al., Nucl. Instrum. and Methods B 388 (2016) 41
! [3] J. Lindhard, Dan. Vid. Selsk Mat.-Fys. Medd., 1954


module CDF_Mermin

use Universal_constants
use CDF_Ritchi, only: Single_diel_func, m_one_third

implicit none

real(8) :: m_over4Pi, m_3Pi2

parameter (m_over4Pi = 1.0d0 / (4.0d0 * g_Pi) )
parameter (m_3Pi2 = 3.0d0*g_Pi*g_Pi )




contains


pure function Mermin_Diel_func(A, E, Gamm, hw, hq) result(f) ! full CDF as sum of single scillators
   real(8) :: f     ! full CDF as a sum of individual functions
   real(8), dimension(:), intent(in) :: A, E, Gamm     ! CDF parameters
   real(8), intent(in) :: hw, hq        ! parameters and variable: transferred energy [eV], recoil momentum [eV]
   !--------------------------------------
   real(8) :: CDF, Mass1, qlim
   integer :: i, N
   
   ! Sum up contributions from all CDF oscillators:
   N = size(A)      ! number of CDF fitting functions
   CDF = 0.0d0      ! to start with
   do i = 1, N      ! for all functions
      CDF = CDF + Single_Mermin_func(A(i), E(i), Gamm(i), hw, hq) ! see below
   enddo
   f = CDF
end function Mermin_Diel_func


! Single Mermin oscillator for CDF, Mermin file:
pure function Single_Mermin_func(A, E, Gamm, W, Q) result(Imepsilon) ! (we use atomic untis below)
   real(8) Imepsilon ! function itself
   real(8), intent(in) :: A, E, Gamm, W, Q ! CDF parameters, W=hw and Q=hq^2/2me [eV]
   !----------------------------------------
   real(8) :: Eat, Gat, Wat, Qat, Qat2  ! CDF parameters in atomic units
   real(8) :: n_e, q_TF2, eps1, eps2, numer, denom, denom1, denom2, eps11, Gateps11, Qat2q_TF2, Qat2q_TF2_1    ! artificial parameters of oscillator
   complex(8) :: eps
   ! Transfer all CDF parameters to atomic units
   Eat = E*g_ev2au
   Gat = Gamm*g_ev2au
   Wat = W*g_ev2au
   Qat2 = 2.0d0*Q*g_ev2au
   !Qat = sqrt(2.0d0*Q*g_ev2au) ! Qat: hq^2/2me -> Qat^2/2
   Qat = sqrt(Qat2)
   !n_e = Eat**2/4.0d0/g_Pi ! artifical density of oscillator
   n_e = Eat*Eat * m_over4Pi
   q_TF2 = (4.0d0*(3.0d0*n_e/g_Pi)**m_one_third) ! square of artifical screening parameter        
   if (Qat .le. 1.0e-3) then    ! near q=0, use Drude-like (Ritchie) oscillator
      Imepsilon = Single_diel_func(A, E, Gamm, W, Q)    ! module "CDF_Ritchi"
   else ! for q>0, calculate full Mermin:
      eps = Lindhard_func(Wat, Qat, Gat, Eat) ! Lidnhard function, below
      !eps1 = real(eps)
      eps1 = dble(eps)
      eps2 = aimag(eps)
      eps11 = eps1 - 1.0d0
      Gateps11 = Gat*eps11
      Qat2q_TF2 = Qat2/q_TF2
      Qat2q_TF2_1 = 1.0d0 + Qat2q_TF2
      ! Im(-1/eps_M) = numerator/denominator, eq.(16) in Mermin file. chi1 ~ 2/(pi*q^3) * (eps1 - 1), chi2 ~ 2/(pi*q^3) * eps2
      !numer = (Wat*(Gat*(eps1 - 1.0d0) + Wat*eps2) - Gat*Wat*Qat**2.0d0/q_TF2*((eps1 - 1.0d0)**2.0d0 + eps2**2.0d0))
      numer = (Wat*(Gateps11 + Wat*eps2) - Gat*Wat*Qat2q_TF2*(eps11*eps11 + eps2*eps2))
      !denom = (Wat*eps1 - eps2*Gat*(1.0d0 + Qat**2.0d0/q_TF2))**2.0d0 &
      !          + (Gat*(eps1 - 1.0d0)*(1.0d0 + Qat**2.0d0/q_TF2) + Wat*eps2)**2.0d0
      denom1 = Wat*eps1 - eps2*Gat*Qat2q_TF2_1
      denom2 = Gateps11*Qat2q_TF2_1 + Wat*eps2
      denom = denom1*denom1 + denom2*denom2
      Imepsilon = numer/denom*A/E**2.0d0
    end if
end function Single_Mermin_func



! Lindhard dielectric function, see [3] Eq.(3.6):
pure function Lindhard_func(W, Q, Gamm, W_pl) result(Lindhard_Epsilon) ! (we use atomic units below)
   complex(8) Lindhard_Epsilon ! function itself
   real(8), intent(in) :: W, Q, Gamm, W_pl ! transferred energy, transferred momentum, damping coefficient, plasma frequency in atomic units
   !----------------------------------------
   real(8) :: n_e, k_f, z, W_pl2, z8
   complex(8) :: u, f, zmu, zpu

   !n_e = W_pl**2/4.0d0/g_Pi ! artifical density of oscillator
   W_pl2 = W_pl*W_pl
   n_e = W_pl2 * m_over4Pi ! artifical density of oscillator
   !k_f = (3.0d0*g_Pi**2*n_e)**m_one_third ! artifical Fermi momentum
   k_f = (m_3Pi2*n_e)**m_one_third ! artifical Fermi momentum
   z = 0.5d0*Q/k_f
   u = cmplx(W, Gamm)/Q/k_f
   zmu = cmplx(z,0.0d0) - u
   zpu = cmplx(z,0.0d0) + u
   z8 = 1.0d0/(8.0d0*z)
   !f = dcmplx(0.5d0, 0.0d0) + (1.0d0 - (z - u)**2.0d0)*log((z - u + 1.0d0)/(z - u - 1.0d0))/(8.0d0*z) + &
   !         (1.0d0 - (z + u)**2.0d0)*log((z + u + 1.0d0)/(z + u - 1.0d0))/(8.0d0*z)
   f = cmplx(0.5d0, 0.0d0) + (1.0d0 - zmu*zmu)*log((zmu + 1.0d0)/(zmu - 1.0d0)) * z8 + &
            (1.0d0 - zpu*zpu)*log((zpu + 1.0d0)/(zpu - 1.0d0)) * z8
   Lindhard_Epsilon = cmplx(1.0d0, 0.0d0) + 3.0d0*W_pl2/(Q*Q*k_f*k_f)*f
end function Lindhard_func


end module CDF_Mermin
