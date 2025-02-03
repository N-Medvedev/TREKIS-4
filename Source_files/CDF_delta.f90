! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019-2020
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains subroutines for delta-function (an extension to Liljequist) model of CDF
! References used:
! [1] N.Medvedev, A.E.Volkov, J. PHys. D (2020)    DOI: https://doi.org/10.1088/1361-6463/ab7c09
! [2] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2008: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2009)


MODULE CDF_delta
use Universal_constants
use Objects
use CDF_Ritchi, only : Int_Ritchi_x
use CS_general_tools, only : W_from_impact_parameter
use CS_integration_limits, only : W_min, W_max, find_Wmax_equal_Wmin, minimal_sufficient_E
use Relativity, only : rest_energy, beta_factor, velosity_from_kinetic_energy
use Little_subroutines, only: find_order_of_number

implicit none
 contains



!-----------------------------------------------------------------------------
! Relativistic CDF cross section from [1], using Delta-CDF
function CDF_total_CS_delta(El_inelast, Ekin, Mass, Zeff, Ip, n_target, CDF, Mt, identical, used_SHI, hw_phonon, Emax_in) result(sigma)
   real(8) :: sigma    ! [A^2]
   integer, intent(in) :: El_inelast    ! type of cross section to be used
   real(8), intent(in) :: Ekin, Zeff ! [eV] kinetic energy of incident particle, effective charge 
   real(8), intent(in) :: Mass, Mt    ! [kg] indident and target (scatterrer) particles masses
   real(8), intent(in) :: Ip, n_target      ! [eV] ionization potential;  [1/cm^3] atomic density
   type(Ritchi_CDF), intent(in) :: CDF	! CDF coefficients
   logical, intent(in) :: identical   ! identical particles, true or not
   integer, intent(in), optional :: used_SHI  ! Z atomic number of SHI, if SHI is used (0 if it is not an SHI)
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8), intent(in), optional :: Emax_in ! [eV] integration limit, used for defining transfered energy
   !-------------------------------
   integer :: used_SHI_local

   if (present(used_SHI)) then ! it is an SHI:
      used_SHI_local = used_SHI
   else ! it is not an SHI (assume electron):
      used_SHI_local = 0
   endif

   if (present(hw_phonon) .and. present(Emax_in) ) then
      sigma = Integral_CDF_delta_CS(El_inelast, Mass, Mt, Ekin, CDF, Ip, n_target, identical, used_SHI_local, &
         hw_phonon=hw_phonon, Emax_in=Emax_in)   ! below
   elseif (present(hw_phonon)) then
      sigma = Integral_CDF_delta_CS(El_inelast, Mass, Mt, Ekin, CDF, Ip, n_target, identical, used_SHI_local, hw_phonon=hw_phonon)   ! below
   elseif (present(Emax_in)) then
      sigma = Integral_CDF_delta_CS(El_inelast, Mass, Mt, Ekin, CDF, Ip, n_target, identical, used_SHI_local, Emax_in=Emax_in)   ! below
   else
      sigma = Integral_CDF_delta_CS(El_inelast, Mass, Mt, Ekin, CDF, Ip, n_target, identical, used_SHI_local)   ! below
   endif
   sigma = sigma*(Zeff*Zeff)    ! include effective charge
end function CDF_total_CS_delta


function Integral_CDF_delta_CS(El_inelast, M, mt, E, CDF, Ip, nat, identical, used_SHI, hw_phonon, Emax_in) result(Sigma)
   real(8) Sigma    ! [A^2]
   integer, intent(in) :: El_inelast    ! type of cross section to be used
   real(8), intent(in) :: M, mt ! [kg] indident and target (scatterrer) particles masses
   real(8), intent(in) :: E ! [eV]
   type(Ritchi_CDF), intent(in), target :: CDF  ! CDF coefficients
   real(8), intent(in) :: Ip    ! [eV] ionization potential
   real(8), intent(in) :: nat   ! [1/cm^3] atomic density
   logical, intent(in) :: identical   ! identical particles, true or not
   integer, intent(in) :: used_SHI  ! Z atomic number of SHI, if SHI is used (0 if it is not an SHI)
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8), intent(in), optional :: Emax_in ! [eV] integration limit, used for defining transfered energy
   !---------------------------------------
   real(8) :: Mc2, mtc2, CS, P, Mc22, Emin, Emax, Eeq, dEed, Estart
   real(8) :: a, b  ! coefficients of the linear extrapolation function
   real(8), dimension(:), pointer :: alpha
   real(8), dimension(size(CDF%E0)) :: E0
   integer :: i, Nosc
   
   ! How many delta-functions:
   Nosc = size(CDF%E0(:))
   
   Mc2 = rest_energy(M)   ! module "Relativity"
   mtc2 = rest_energy(mt)   ! module "Relativity"
   ! Parameters of the delta-functions CDF:
!    if (El_inelast == 5) then    ! Single Pole Delta CDF
!       E0 = max(Ip, sqrt(CDF%h_omega_e2))
!       alpha => CDF%h_omega_e2
!    else  ! Delta CDF
      E0 = CDF%E0
   if (used_SHI > 0) then ! it is an SHI, use its renormalized alpha:
      alpha => CDF%alpha_SHI(:,used_SHI)
   else ! it is not an SHI, use electron's renormalized alpha:
      alpha => CDF%alpha
   endif
!    endif
!    E0 = CDF%E0
!    print*, E0, sqrt(CDF%h_omega_e2), (CDF%h_omega_e2), CDF%alpha
   
   CS = 0.0d0   ! to start with
   do i = 1, Nosc
      Emin = W_min(Ip, Mc2, mtc2, E, E0(i)) ! module "CS_integration_limits"

      if (abs(M - mt)/M < 1.0d-6) then  ! identical particles:
         Estart = Ip
      else  ! non-identical
         Estart = minimal_sufficient_E(Ip, Mc2, mtc2) ! module "CS_integration_limits"
      endif
      
      ! Where to switch from the approximate formula to linear extrapolation:
      if (present(hw_phonon)) then
         call find_Wmax_equal_Wmin(Mc2, mtc2, identical, E, Ip, E0(i), Eeq, hw_phonon)   ! module "CS_integration_limits"
      else
         call find_Wmax_equal_Wmin(Mc2, mtc2, identical, E, Ip, E0(i), Eeq)   ! module "CS_integration_limits"
      endif
      
      if (E <= Estart) then ! cannot ionize
         CS = 0.0d0
         P = 0.0d0
      else
         dEed = Eeq/100.0d0    ! shift a little the point at which we switch to linear extrapolation
         
         if (E <= Eeq+dEed) then   ! delta-CDF model does not apply here, replace the CDF with the linear extrapolation cross section
            print*, 'Integral_CDF_delta_CS too low', E, Eeq
            ! Find the coefficients of the linear approximation:
            call Find_linear_a_b(alpha(i), M, Mc2, mtc2, E0(i), Ip, identical, nat, Eeq+dEed, a, b) ! below
            CS = a*E + b    ! linearly approximated cross section
            P = 1.0d0   ! no prefactor needed in this case, it alreaxy was included into the linear extrapolation
         else   ! delta-model works
            if (present(hw_phonon)) then  ! scattering on phonons:
               Emax = W_max(Mc2, mtc2, identical, E, Ip, hw_phonon) ! module "CS_integration_limits"
            else  ! scattering on electrons or atoms:
               Emax = W_max(Mc2, mtc2, identical, E, Ip) ! module "CS_integration_limits"
            endif
            if (present(Emax_in)) then
               if (Emax_in < Emin) then
                  Emax = Emin
               elseif (Emax_in < Emax) then
                  Emax = Emax_in
               endif
            endif
!             CS = CS - (integrated_delta_CDF_CS(alpha(i), Mc2, E0(i), mtc2, Emax, Ip, E) - integrated_delta_CDF_CS(alpha(i), Mc2, E0(i), mtc2, Emin, Ip, E))   ! below
            CS = CS - Total_delta_CDF_CS(alpha(i), Mc2, E0(i), mtc2, Emax, Emin, Ip, E)   ! below
!             ! Prefactor:
!             P = P_prefactor(M, E, nat)   ! module "CDF_delta"
         endif ! (E <= Eeq) 
      endif ! (E <= Ip)
   enddo ! i = 1, Nosc
   
   ! Prefactor:
   P = P_prefactor(M, E, nat)   ! module "CDF_delta"
   
   ! Total cross section:
   Sigma = abs(CS)*P   ! [A^2]

   nullify(alpha)
end function Integral_CDF_delta_CS


function Total_delta_CDF_CS(alpha, Mc2, E0, mtc2, Emax, Emin, Ip, E) result(CS)
   real(8) CS   ! a part of the cross section
   real(8), intent(in) :: alpha, Mc2, E0, mtc2, Ip, E, Emax, Emin
   !------------------------------------------
   real(8) :: Mc22, Wmin, t3, WminE0, eps
   parameter(eps = 1.0d-15)
   
   ! Low integration limit:
   Wmin = W_min(Ip, Mc2, mtc2, E, E0) ! module "CS_integration_limits"
   
   if (Emin > Wmin) then
      Wmin = Emin 
   endif
   
   Mc22 = 2.0d0*Mc2
   if (E <= Ip) then ! cannot ionize or transfer energy
      CS = 0.0d0
   else
      WminE0 = (Wmin-E0)
      if (WminE0 < 0.0d0) then  ! just in case something went wrong
         t3 = 0.0d0
      elseif (WminE0 < eps) then    ! approximate solution, Appendix A [1]
         t3 = dlog(Emax - E0) - ( 2.0d0*dlog(E0) + dlog(Mc2) - dlog(4.0d0*E*mtc2) )
      else  ! regular solution
         !t3 = dlog(abs((Emax - E0)/(Wmin - E0)))
         t3 = dlog(Emax - E0) - dlog(abs(Wmin - E0))
      endif
      !CS = alpha/(g_me_eV * (Mc22 - E0)) * ( (Mc22 - mtc2)*dlog((Mc22 + Emax - E0)/(Mc22 + Wmin  - E0) ) + &
      !         Mc22*(mtc2/E0-1.0d0)*dlog(Emax/Wmin ) -  mtc2*(Mc22/E0-1.0d0)*t3 )
      CS = alpha/(g_me_eV * (Mc22 - E0)) * ( (Mc22 - mtc2)*dlog((Mc22 + Emax - E0)/(Mc22 + Wmin  - E0) ) + &
               Mc22*(mtc2/E0-1.0d0)* ( dlog(Emax) - dlog(Wmin) ) -  mtc2*(Mc22/E0-1.0d0)*t3 )

!      print*, 'Total_delta_CDF_CS', E, Emax, E0, Emax - E0, abs(Wmin - E0)
   endif
end function Total_delta_CDF_CS



function integrated_delta_CDF_CS(alpha, Mc2, E0, mtc2, W, Ip, E) result(CS)
   real(8) CS   ! a part of the cross section
   real(8), intent(in) :: alpha, Mc2, E0, mtc2, W, Ip, E
   !------------------------------------------
   real(8) :: Mc22, Wmin
   ! Low integration limit:
   Wmin = W_min(Ip, Mc2, mtc2, E, E0) ! module "CS_integration_limits"

   if ((W < Wmin) .or. (E <= Ip)) then ! cannot ionize or transfer energy
      CS = 0.0d0
   else
      CS =  integral_CS(alpha, Mc2, mtc2, E0, W)    ! below
   endif
end function integrated_delta_CDF_CS



function energy_loss_delta(El_inelast, E, M, ZSHI, Zeff, Ip, nat, Mt, CDF, identical, used_SHI, hw_phonon) result(Se)
   real(8) :: Se    ! [eV/A] energy loss
   integer, intent(in) :: El_inelast    ! type of cross section to be used
   real(8), intent(in) :: E, M, ZSHI, Zeff, Ip, nat, Mt
   type(Ritchi_CDF), intent(in), target :: CDF	! CDF coefficients
   logical, intent(in) :: identical   ! identical particles, true or not
   integer, intent(in) :: used_SHI ! Z atomic number of SHI (0 if it is not an SHI)
   real(8), intent(in), optional :: hw_phonon   ! [eV]
   !---------------------------------------
   real(8) :: M_units, Mc2, mtc2, S_cur, P, Mc22, Emin, Emax, Eeq, dEed, Estart
   real(8) :: a, b  ! coefficients of the linear extrapolation function
   real(8) :: v, beta, beta2, W_max_nd, b_max
   real(8), dimension(:), pointer :: alpha
   real(8), dimension(size(CDF%E0)) :: E0
   integer :: i, Nosc
   
!    M = M_in*g_amu     ! [kg]
   if (M > g_me*1.1d0) then ! not an electron, assume ion:
      M_units = M/g_amu    ! number of nuclei
   else ! an electron
      M_units = 1.0d0
   endif
   
   ! How many delta-functions:
   Nosc = size(CDF%E0(:))
   
   v = velosity_from_kinetic_energy(E, M, afs=.false.) ! [m/s] module "Relativity"
   beta = beta_factor(v)    ! module "Relativity"
   beta2 = beta*beta
   Mc2 = rest_energy(M)   ! module "Relativity"
   mtc2 = rest_energy(mt)   ! module "Relativity"
   
   if (M > g_me*1.1d0) then ! not an electron, assume ion:
      b_max = 1.2e-5*(M_units)**(1.0d0/3.0d0) ! estimate of the nuclear size [A]
      !W_max_nd = W_from_impact_parameter(Mc2, mtc2, E, dble(ZSHI), 1.0d0, b_max)  ! module "CS_general_tools"
      W_max_nd = 10.0e26
   else
      W_max_nd = 1e25   ! do not use for electrons, only for ions
   endif   
   
   ! Parameters of the delta-functions CDF:
!    if (El_inelast == 5) then    ! Single Pole Delta CDF
! !       E0 = Ip
!       E0 = max(Ip, sqrt(CDF%h_omega_e2))
!       alpha => CDF%h_omega_e2
!    else  ! Delta CDF
      E0 = CDF%E0
      !alpha => CDF%alpha
   if (used_SHI > 0) then ! it is an SHI, use its renormalized alpha:
      alpha => CDF%alpha_SHI(:,used_SHI)
   else ! it is not an SHI, use electron's renormalized alpha:
      alpha => CDF%alpha
   endif
!    endif

   S_cur = 0.0d0   ! to start with
   do i = 1, Nosc
      Emin = W_min(Ip, Mc2, mtc2, E, E0(i)) ! module "CS_integration_limits"

      if (abs(M - mt)/M < 1.0d-6) then  ! identical particles:
         Estart = Ip
      else  ! non-identical
         Estart = minimal_sufficient_E(Ip, Mc2, mtc2) ! module "CS_integration_limits"
      endif
      
      ! Where to switch from the approximate formula to linear extrapolation:
      call  find_Wmax_equal_Wmin(Mc2, mtc2, identical, E, Ip, E0(i), Eeq)   ! module "CS_integration_limits"
      
      if (E <= Estart) then ! cannot ionize
         S_cur = 0.0d0
         P = 0.0d0
      else
         dEed = Eeq/100.0d0    ! shift a little the point at which we switch to linear extrapolation
         
         if (present(hw_phonon)) then  ! scattering on phonons:
            Emax = W_max(Mc2, mtc2, identical, E, Ip, hw_phonon) ! module "CS_integration_limits"   
         else  ! scattering on electrons or atoms:
            Emax = W_max(Mc2, mtc2, identical, E, Ip) ! module "CS_integration_limits"
         endif
         ! Finite nucleus size:
!          if (Emax > W_max_nd) Emax = W_max_nd ! Londhard-Sorsen correction to concider...

         if (E <= Eeq+dEed) then   ! delta-CDF model does not apply here, replace the CDF with the linear extrapolation cross section
            ! Find the coefficients of the linear approximation:
            call Find_linear_a_b(alpha(i), M, Mc2, mtc2, E0(i), Ip, identical, nat, Eeq+dEed, a, b) ! below
            S_cur = 0.5d0*a*(Emax - Estart)    ! linearly approximated cross section
            P = 1.0d0   ! no prefactor needed in this case, it alreaxy was included into the linear extrapolation
         else   ! delta-model works
            !S_cur = S_cur - (integral_CS_x_W(alpha(i), Mc2, E0(i), Estart, mtc2, Emax) - integral_CS_x_W(alpha(i), Mc2, E0(i), Estart, mtc2, Emin)) ! below
            S_cur = S_cur - int_energy_loss(E, alpha(i), M, Mc2, mtc2, E0(i), Emax, Emin)  ! below
            !S_cur = S_cur + (log(1.0d0-beta2) + beta2)*Mc2/g_me_eV*alpha(i)
!             write(*,'(a,e,i,e,e,e,e)') 'alpha:', E, i, S_cur, alpha(i), Mc2, mtc2
            
!             ! Prefactor:
!             P = P_prefactor(M, E, nat)   ! below
         endif ! (E <= Eeq) 
      endif ! (E <= Ip)
   enddo ! i = 1, Nosc
   
   ! Prefactor:
   P = P_prefactor(M, E, nat)   ! below
   ! Energy loss including effective charge:
!    print*, 'energy_loss_delta:', E, S_cur, beta2
   Se = abs(S_cur)*(Zeff*Zeff)*P*nat*1.0d-24   ! [eV/A]
   !- (log(1.0d0-beta2) + beta2)*Mc2/g_me_eV
   
   nullify(alpha)   
end function energy_loss_delta

 
 


! Find the coefficients of the linear approximation:
 subroutine Find_linear_a_b(alpha, M, Mc2, mtc2, E0, Ip, identical, nat, Eeq, a, b)
   real(8), intent(in) :: alpha, M, Mc2, mtc2, E0, Ip
   logical, intent(in) :: identical
   real(8), intent(in) :: nat   ! [1/cm^3] atomic density
   real(8), intent(in) :: Eeq   ! [eV] point where W_min = W_max - there we switch from delta-model to the linear extrapolation
   real(8), intent(out) :: a, b  ! coefficients of the linear extrapolation function
   real(8) :: Wmin_lim, Wmax_lim, CS, P, IpMm
   
   Wmin_lim = W_min(Ip, Mc2, mtc2, Eeq, E0) ! module "CS_integration_limits"
   Wmax_lim = W_max(Mc2, mtc2, identical, Eeq, Ip) ! module "CS_integration_limits"
   P = P_prefactor(M, Eeq, nat)   ! below
   CS = -P * (integral_CS(alpha, Mc2, mtc2, E0, Wmax_lim) - integral_CS(alpha, Mc2, mtc2, E0, Wmin_lim))  ! below
   
   IpMm = minimal_sufficient_E(Ip, Mc2, mtc2)    ! module "CS_integration_limits"
   
   a = CS/(Eeq-IpMm)
   b = -CS*Ip/(Eeq-IpMm)
end subroutine Find_linear_a_b


function integral_CS(alpha, Mc2, mtc2, E0, W) result(CS) ! Eq.(16) [1]
   real(8) CS   ! part of the cross section of scattering (not normalized, prefactor missing)
   real(8), intent(in) :: alpha, Mc2, E0, mtc2, W
   real(8) :: mtc22, Mc22
   Mc22 = 2.0d0*Mc2
   mtc22 = 2.0d0*mtc2
!   CS = alpha/(g_me_eV * (Mc22 - E0)) * ( (Mc22 - mtc2)*log(Mc22 + W - E0) + Mc22*mtc2/E0*( log(W) - log(abs(W - E0))) )
!   CS = alpha/(g_me_eV * (Mc22 - E0)) * ( (Mc22 - mtc2)*log(Mc22 + W - E0) + Mc22*(mtc2/E0-1.0d0)*log(W) -  mtc2*(Mc22/E0-1.0d0)*log(abs(W - E0)) )
   ! Corrected:
   CS = alpha*mtc2/(g_me_eV*(mtc22-E0))*( log(Mc22 + W - E0) + 2.0d0*(mtc2/E0-1.0d0)*log(W) - (mtc22/E0-1.0d0)*log(abs(W - E0)) )
end function integral_CS



pure function int_energy_loss(E, alpha, M, Mc2, mtc2, E0, Wmax, Wmin) result (WCS) ! Eq.(18) [1]
   real(8) :: WCS
   real(8), intent(in) :: E, alpha, M, Mc2, mtc2, Wmin, Wmax, E0
   real(8) :: Mc22, part1, part2
   Mc22 = 2.0d0*Mc2
   !WCS =-alpha*( (mtc2*log(abs((Wmax-E0)/(Wmin-E0))) + (Mc22-mtc2)*log(abs((-Wmax+E0-Mc22)/(-Wmin+E0-Mc22))))/g_me_eV)
   ! Corrected:
   if (abs(Wmin - E0) < 1.0d-14) then
      part1 = log(abs(Wmax-E0))
   else
      part1 = log(abs(Wmax-E0)) - log(abs(Wmin-E0))
   endif
   part2 = log(abs(-Wmax+E0-Mc22)) - log(abs(-Wmin+E0-Mc22))
   WCS =-alpha/g_me_eV*mtc2*( part1 + part2 )
end function int_energy_loss



pure function integral_CS_x_W(alpha, Mc2, E0, Ip, mtc2, W) result (WCS)
   real(8) :: WCS
   real(8), intent(in) :: alpha, Mc2, E0, mtc2, Ip, W
   real(8) :: Mc22, Wmin, IpMm
   ! Low integration limit:
   Wmin = W_min(Ip, Mc2, mtc2, W, E0) ! module "CS_integration_limits"
   IpMm = minimal_sufficient_E(Ip, Mc2, mtc2)    ! module "CS_integration_limits"
   if ((W < Wmin) .or. (W <= IpMm)) then ! cannot ionize or transfer energy
      WCS = 0.0d0
   else
      Mc22 = 2.0d0*Mc2
      WCS = -alpha/g_me_eV*( mtc2*log(abs(W-E0)) + (-mtc2+Mc22)*log(W-E0+Mc22) )
   endif
end function integral_CS_x_W




pure function P_prefactor(M, E, nat) result(P)
   real(8) :: P
   real(8), intent(in) :: M, E  ! mass [kg] and energy [eV] of the projectile
   real(8), intent(in) :: nat   ! [1/cm^3] atomic density
   real(8) :: Mc2, beta, v
   ! beta-factor: v/c
   v = velosity_from_kinetic_energy(E, M, afs=.false.) ! [m/s] module "Relativity"
   beta = beta_factor(v)    ! module "Relativity"
   P = 1.0d24/(g_Pi*g_a0*nat*g_me_eV*(beta*beta)) ! converting into [A^2/eV]
end function P_prefactor




pure subroutine define_alpha(Ai, Gammai, E0i, x_min, alpha, mass_target) ! Eq.(12) [1]
   real(8), intent(in) :: Ai, Gammai, E0i, x_min   ! [eV] the coefficients of Ritchi's CDF, and x_min=Ip or Egap, the starting point
   real(8), intent(out) :: alpha ! coefficient in delta-CDF to be fitted
   real(8), intent(in), optional :: mass_target ! [kg] mass of the target particle (in case it is not an electron)

   alpha = g_Pi*Ai*0.5d0     ! Eq.(12) [1]
   
   if (present(mass_target)) then   ! account for the target particle mass in the sum rule:
      alpha = g_me/mass_target * alpha
   endif
   
end subroutine define_alpha



pure subroutine define_alpha_OLD(Ai, Gammai, E0i, x_min, alpha, mass_target) ! Eq.(12) [1]
   real(8), intent(in) :: Ai, Gammai, E0i, x_min   ! [eV] the coefficients of Ritchi's CDF, and x_min=Ip or Egap, the starting point
   real(8), intent(out) :: alpha ! coefficient in delta-CDF to be fitted
   real(8), intent(in), optional :: mass_target ! [kg] mass of the target particle (in case it is not an electron)
   real(8) :: low_lim, high_lim
   low_lim =  Int_Ritchi_x(Ai,E0i,Gammai,x_min)     ! module "CDF_Ritchi"
   high_lim = Int_Ritchi_x(Ai,E0i,Gammai,1.0d30)   ! module "CDF_Ritchi"
   alpha = high_lim - low_lim   ! normalized to reproduce K-sum rule
   
   if (present(mass_target)) then   ! account for the target particle mass in the sum rule:
      alpha = g_me/mass_target * alpha
   endif
end subroutine define_alpha_OLD
 

 
 
 END MODULE CDF_delta
