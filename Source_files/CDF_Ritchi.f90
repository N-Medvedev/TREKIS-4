! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with optical data and convert them into CDF
! [1] N. Medvedev, et al., J. Phys. D: Appl. Phys. 48 (2015) 355303
! [2] R. Rymzhanov et al., Nucl. Instrum. and Methods B 388 (2016) 41
! [3] M. Medvedev, A.E. Volkov, J. Phys. D: Appl. Phys. 53 (2020) 235302

MODULE CDF_Ritchi
use Universal_constants
use Little_subroutines, only: find_in_array_monoton
use Relativity, only: kinetic_energy_from_momentum

implicit none

real(8) :: m_one_third, m_two_third

parameter (m_one_third = 1.0d0/3.0d0)
parameter (m_two_third = 2.0d0/3.0d0)

 contains

 

pure function Diel_func(A, E, Gamma, hw, hq, v_f, E0_model, Mass) result(f) ! full CDF as sum of single scillators
   real(8) :: f     ! full CDF as a sum of individual functions
   real(8), dimension(:), intent(in) :: A, E, Gamma     ! CDF parameters
   real(8), intent(in) :: hw, hq        ! parameters and variable: transferred energy [eV], recoil momentum [eV]
   real(8), intent(in), optional ::  v_f    ! target fermi velosity [sqrt(eV/kg)]
   integer, intent(in), optional :: E0_model     ! whilch model to use for the extension of the dielectric functions
   real(8), intent(in), optional ::  Mass           ! mass [kg]
   !--------------------------------------
   real(8) :: CDF, Mass1, qlim
   integer :: i, N, j
   
   ! Sum up contributions from all CDF oscillators:
   N = size(A)      ! number of CDF fitting functions
   CDF = 0.0d0     ! to start with
   ! Select among various models of extension of the optical dielectric function into finite-q space:
   if (present(v_f) .and. present(E0_model) .and. present(Mass)) then
      do i = 1, N   ! for all functions
         CDF = CDF + Single_diel_func(A(i), E(i), Gamma(i), hw, hq, v_f=v_f, E0_model=E0_model, Mass=Mass)  ! see below
      enddo
   elseif (present(E0_model)) then
      do i = 1, N   ! for all functions
         CDF = CDF + Single_diel_func(A(i), E(i), Gamma(i), hw, hq, E0_model=E0_model)   ! see below
      enddo
   else ! no extra parameters
      do i = 1, N   ! for all functions
         CDF = CDF + Single_diel_func(A(i), E(i), Gamma(i), hw, hq) ! see below
      enddo
   endif
   f = CDF
end function Diel_func


pure function Single_diel_func(A, E, Gamma, W, Q, v_f, E0_model, Mass) result(Diel_func) ! single oscillator for CDF, [1] Eq.(6)
   real(8) Diel_func    ! function itself
   real(8), intent(in) :: A, E, Gamma, W, Q     ! CDF parameters, W=hw and Q=hq^2/2me [eV]
   real(8), intent(in), optional ::  v_f                ! target effective mass; fermi velosity in units of [sqrt(eV/kg)]
   integer, intent(in), optional :: E0_model      ! whilch model to use for the extension of the dielectric function
   real(8), intent(in), optional ::  Mass            ! mass [kg]
   !----------------------------------------
   real(8) :: E0, hw2, E02, hw2E02, Gamma1, hq2

   Gamma1 = Gamma   ! unless changed later in some of the models
   ! Model for E0 and Gamma parameters extension into finite dq-space:
   if (present(E0_model)) then
      select case (E0_model)    ! target dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie
      case (1)  ! Plasmon-pole approximation
         if (present(v_f) .and. present(Mass)) then
!             pause 'test 0'
            hq2 = Q*Mass    ! [kg*eV]
            E0 = sqrt(E*E + v_f*v_f*hq2*m_one_third + Q*Q)  ! [2] Eq.(3)
         else   ! Free electron
!             pause 'test 1'
            E0 = E + Q  ! relativistic version of [2] below Eq.(2)
         endif
      case (2)  ! Extended dielectric model of Ritchie
!          pause 'test 2'
         E0 = (E**m_two_third + (Q)**m_two_third)**1.5d0    ! [2] Eq.(4)
         Gamma1 = sqrt(Gamma*Gamma + Q*Q)
      case default  ! Free electron
!          pause 'test 3'
         E0 = E + Q ! relativistic version of [2] below Eq.(2)
      end select
   else ! by default, it's a free-electron dispersion relation:
      E0 = E + Q    ! relativistic version of [2] below Eq.(2)
!       print*, 'test 4'
   endif

   hw2 = W*W
   E02 = E0*E0
   hw2E02 = hw2 - E02
   Diel_func = A*Gamma1*W/(hw2E02*hw2E02 + Gamma1*Gamma1*hw2)   ! [1] Eq.(6)
!    write(*,'(a,f,f,f,f,f)') 'E0', W, Q, Diel_func
end function Single_diel_func


pure function Find_Ritchie_peak_Q(E0i,W) result(peak_Q)
   real(8) peak_Q   ! position of the peak of the Ritchie (Drude) oscillator in Q
   real(8), intent(in) :: E0i, W    ! oscillators position [eV] and transferred eenrgy [eV]
   peak_Q = W-E0i   ! use the positive solution; the negative solutions are: -E0i; -W-E0i
end function Find_Ritchie_peak_Q


pure function Find_Ritchie_oscillators_center(E0i,Gamm,Q) result(peak_E)
   real(8) peak_E   ! position of the peak of the Ritchie (Drude) oscillator in E
   real(8), intent(in) :: E0i, Q, Gamm    ! oscillators position [eV] and transferred momentum [eV], and oscillator width [eV]
   real(8) :: E0iQ2, G2
   E0iQ2 = (E0i + Q)**2
   G2 = Gamm**2
   peak_E = 1.0d0/6.0d0 * sqrt( 12.0d0*E0iQ2 - 6.0d0*G2 + 6.0d0*sqrt(16.0d0*E0iQ2**2 - 4.0d0*E0iQ2*G2 + G2**2) ) ! Eq.(22) [3]
end function Find_Ritchie_oscillators_center


pure function estimate_A(yi, xi, E0i, Gammai) result(A)	! for given Eo and Gamma, estimate A as the coefficient that makes fit coincide with the peak height
   real(8) :: A
   real(8), intent(in) :: yi, xi, E0i, Gammai
   real(8) :: x2, x_E0
   x2 = xi*xi
   x_E0 = x2 - E0i*E0i
   A = yi * ( x_E0*x_E0 + Gammai*Gammai*x2 )/(Gammai * xi)
end function estimate_A


pure function dCDF_dA(A,E,Gamma,x) result(df_dA)	! derivative of CDF by parameter A
   real(8) :: df_dA	! derivative by parameter A
   real(8), intent(in) :: A, E, Gamma, x		! parameters of CDF
   real(8) :: ksi, br, x2
   x2 = x*x
   br = x2 - E*E
   ksi = br*br + Gamma*Gamma*x2
   df_dA = Gamma*x / ksi
end function dCDF_dA


pure function dCDF_dGamma(A,E,Gamma,x) result(df_dg)	! derivative of CDF by parameter gamma
   real(8) :: df_dg	! derivative by parameter gamma
   real(8), intent(in) :: A, E, Gamma, x		! parameters of CDF
   real(8) :: ksi, br, x2, g2
   x2 = x*x
   br = x2 - E*E
   g2 = Gamma*Gamma
   ksi = br*br + g2*x2
   df_dg = A*x / ksi * (1.0d0 - 2.0d0*g2*x2 / ksi)
end function dCDF_dGamma


pure function dCDF_dE0(A,E,Gamma,x) result(df_dE0)	! derivative of CDF by parameter E0
   real(8) :: df_dE0	! derivative by parameter E0
   real(8), intent(in) :: A, E, Gamma, x		! parameters of CDF
   real(8) :: ksi, br, x2, g2, E2
   x2 = x*x
   E2 = E*E
   br = x2 - E2
   g2 = Gamma*Gamma
   ksi = br*br + g2*x2
   df_dE0 = 4.0d0*A*Gamma*E*x*br / (ksi*ksi)
end function dCDF_dE0


pure function w_plasma(At_dens, Mass)	! plamsa frequency as defined per atom, see [1] under Eq.(8)
   real(8) w_plasma ! Squared atomic plasma frequency [1/s]^2
   real(8), intent(in) :: At_dens	! [atom/m^3]
   real(8), intent(in), optional :: Mass    ! [kg]
!    w_plasma = (4.0d0*g_Pi*At_dens*g_e*g_e/(4.0d0*g_Pi*g_e0*g_me))
   if (present(Mass)) then  ! use user-provided mass:
      w_plasma = At_dens*g_e*g_e/(g_e0*Mass)
   else ! assume free electron
      w_plasma = At_dens*g_e*g_e/(g_e0*g_me)
   endif
end function w_plasma


subroutine sum_rules(Afit, Efit, Gfit, x_min, Omega, ksum, fsum)	! function to calculated both sum rules for given shell and CDF parameters [1], Eqs.(8-9)
    real(8), dimension(:),  intent(in), target :: Afit ! A coefficient
    real(8), dimension(:),  intent(in), target :: Efit ! E0 coefficient
    real(8), dimension(:),  intent(in), target :: Gfit ! Gamma coefficient
    real(8), intent(in) ::  x_min, Omega ! ionization potneital [eV] and plasma frequency [1/s^2]
    real(8), intent (out), optional :: ksum, fsum
    !-----------------------
    real(8), pointer ::  A1, E1, Gamma1
    real(8) ne, f
    integer j
    ne = 0.0d0
    f = 0.0d0
    do j = 1,size(Afit)
        A1 => Afit(j)
        E1 => Efit(j)
        Gamma1 => Gfit(j)
        ! Integrals from Ip to infinity (infinity = 1d20):
        ne = ne + Int_Ritchi_x(A1,E1,Gamma1,1d20) - Int_Ritchi_x(A1,E1,Gamma1,x_min)	! functions below
        f = f + Int_Ritchi_p_x(A1,E1,Gamma1,1d20) - Int_Ritchi_p_x(A1,E1,Gamma1,x_min)	! functions below
    enddo
    if (present(ksum)) ksum = 2.0d0*g_e*g_e/(g_Pi*Omega*g_h*g_h)*ne ! factors are converting Omega from [1/s]^2 -> [eV]^2
    if (present(fsum)) fsum = f*2.0d0/g_Pi
    nullify(A1, E1, Gamma1)
end subroutine sum_rules


pure function Int_Q_Ritchi_by_W(A, E, G, Q, Wmin, Wmax) result(Int_eps) ! analytical integral of the (1/Q*Ritchi) by W for Q/=0
    real(8) Int_eps
    real(8), intent(in) :: A, E, G  ! parameters and variable
    real(8), intent(in) :: Q    ! [eV] transfered momentum
    real(8), intent(in) :: Wmin, Wmax   ! [eV] min and max transferred energy
    real(8) :: E2, g2, EQ, denom_sqrt, atan_r_max, atan_r_min

    E2 = E*E
    g2 = G*G
    EQ = E*Q
    ! Factor to be reused:
    denom_sqrt = sqrt(4.0d0*(E2 + 2.0d0*EQ + Q**2) - g2)
    ! Upper integration limit:
    atan_r_max = atan_for_int_Rithie(E, G, Q, Wmax, denom_sqrt) ! below
    ! Lower integration limit:
    atan_r_min = atan_for_int_Rithie(E, G, Q, Wmin, denom_sqrt) ! below
    ! Combine all together:
    Int_eps = A/(Q*denom_sqrt)*(atan_r_min - atan_r_max)
end function Int_Q_Ritchi_by_W


pure function atan_for_int_Rithie(E, G, Q, W, denom_sqrt) result(atan_ritchie)
   real(8) atan_ritchie
   real(8), intent(in) :: E, G   ! CDF parameters
   real(8), intent(in) :: Q             ! [eV] transfered momentum
   real(8), intent(in) :: W             ! transferred energy
   real(8), intent(in) :: denom_sqrt    ! precalculated denominator
   real(8) :: E2
   E2 = E*E
   atan_ritchie = ATAN( ( 2.0d0*(-W**2 + E2 + 2.0d0*E*Q + Q**2) - G**2)/(G*denom_sqrt) )
end function atan_for_int_Rithie



pure function Int_Ritchi(A,E,Gamma,x) ! analytical integral of the Ritchi function for q=0
    real(8), intent(in) :: A, E, Gamma, x  ! parameters and variable
    real(8) Int_Ritchi ! function itself
    real(8) S, E2, g2
    complex(8) Sc
    E2 = E*E
    g2 = Gamma*Gamma
    !Sc = (sqrt((2.0d0*E)*(2.0d0*E) - Gamma*Gamma))
    Sc = (sqrt(4.0d0*E2 - g2))
    S = dble(Sc)
    Int_Ritchi = A/S*atan( (2.0d0*(x*x-E2) + g2)/(Gamma*S) )
end function Int_Ritchi


pure function Int_Ritchi_x(A,E,Gamma,x) ! analytical integral of the Ritchi*x (k-sum rule), analytical integral of Eq.(8) [1]
   ! Note that this sum rule is for an electronic system; for atomic one, there must be accounting of the atomic mass in the prefactor!
    real(8), intent(in) :: A, E, Gamma, x  ! parameters and variable
    real(8) Int_Ritchi_x ! function itself
    real(8) S, sq2, G, s_plus, s_minus,B, g2
    complex(8) Sc, Gc, s_plus_c, s_minus_c, Bc, Ic, arg,arg2, oneI, GammaSc, Ic1
    sq2 = sqrt(2.0d0)
!     oneI = cmplx(0.0d0,1.0d0)
    oneI = g_CI
    g2 = Gamma*Gamma
    if ((-(2.0d0*E)*(2.0d0*E) + g2) .GE. 0.0d0) then
        Sc = cmplx(sqrt(-(2.0d0*E)*(2.0d0*E) + g2),0.0d0)
    else
        Sc = cmplx(0.0d0, sqrt(ABS(-(2.0d0*E)*(2.0d0*E) + g2)))
    endif
    Gc = cmplx(g2 - 2.0d0*E*E,0.0d0)
    GammaSc = Gamma*Sc
    s_plus_c = sq2*sqrt(Gc + GammaSc)
    s_minus_c = sq2*sqrt(Gc - GammaSc)
    Bc=Gc/(GammaSc)
    !Int_Ritchi_x = real(A*Gamma*( ATAN(2.0e0*x/s_minus)/s_minus*(1.0e0-B) + ATAN(2.0e0*x/s_plus)/s_plus*(1.0e0+B) ))
    if (s_minus_c == cmplx(0.0d0,0.0d0)) then
       Ic1 = cmplx(0.0d0,0.0d0)
    else
       arg = 2.0d0*x/(s_minus_c)
       Ic1 = 0.5d0*oneI*(log(1.0d0-oneI*arg)-log(1.0d0+oneI*arg))/s_minus_c*(1.0d0-Bc)
    endif
    arg2 = 2.0d0*x/(s_plus_c)
    !Ic = (A*Gamma*( ATAN2(aimag(arg),real(arg))/s_minus_c*(1.0e0-Bc) + ATAN2(aimag(arg2),real(arg2))/s_plus_c*(1.0e0+Bc) ))
    Ic = Ic1 + (0.5d0*oneI*(log(1.0d0-oneI*arg2)-log(1.0d0+oneI*arg2)))/s_plus_c*(1.0d0+Bc)
    Ic = A*Gamma*Ic
    Int_Ritchi_x = dble(Ic)
end function Int_Ritchi_x


pure function Int_Ritchi_p_x(A,E,Gamma,x) ! analytical integral of the Ritchi/x (ff-sum rule), analytical integral of Eq.(9) [1]
    real(8), intent(in) :: A, E, Gamma, x  ! parameters and variable
    real(8) Int_Ritchi_p_x ! function itself
    real(8) S, sq2, G, s_plus, s_minus
    complex(8) Sc, Gc, s_plus_c, s_minus_c, oneI, In_c, arg, arg2, ln_c1
    sq2 = sqrt(2.0d0)
    oneI = cmplx(0.0d0,1.0d0)
    if ((-(2.0d0*E)*(2.0d0*E) + Gamma*Gamma) .GE. 0.0d0) then
        Sc = cmplx(sqrt(-(2.0d0*E)*(2.0d0*E) + Gamma*Gamma),0.0d0)
    else
        Sc = cmplx(0.0d0, sqrt(ABS(-(2.0d0*E)*(2.0d0*E) + Gamma*Gamma)))
    endif
    Gc = cmplx(Gamma*Gamma - 2.0d0*E*E,0.0d0)
    s_plus_c = sqrt(Gc + Gamma*Sc)
    s_minus_c = sqrt(Gc - Gamma*Sc)
    !Int_Ritchi_p_x = sq2*A/S*( ATAN(sq2*x/s_minus)/s_minus - ATAN(sq2*x/s_plus)/s_plus )
    if (s_minus_c == cmplx(0.0d0,0.0d0)) then
       ln_c1 = cmplx(0.0d0,0.0d0)
    else
       arg = sq2*x/s_minus_c
       ln_c1 = (0.5d0*oneI*(log(1.0d0-oneI*arg)-log(1.0d0+oneI*arg)))/s_minus_c
    endif
    arg2 = sq2*x/s_plus_c
    !In_c = sq2*A/Sc*( (0.5d0*oneI*(log(1.0d0-oneI*arg)-log(1.0d0+oneI*arg)))/s_minus_c - (0.5d0*oneI*(log(1.0d0-oneI*arg2)-log(1.0d0+oneI*arg2)))/s_plus_c )
    In_c = ln_c1 - (0.5d0*oneI*(log(1.0d0-oneI*arg2)-log(1.0d0+oneI*arg2)))/s_plus_c
    In_c = sq2*A/Sc*In_c
    Int_Ritchi_p_x = dble(In_c)
end function Int_Ritchi_p_x

 
 
END MODULE CDF_Ritchi
