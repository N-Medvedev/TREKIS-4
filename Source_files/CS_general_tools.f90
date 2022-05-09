! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains general subroutines useful for calculating cross sections and related quantities
! References:
! [1] Medvedev et al., Phys. Rev. B 82, 125425 (2010)

module CS_general_tools

use Universal_constants
use Objects
use Little_subroutines, only : interpolate_data_single


implicit none

 contains
 
 

pure function W_from_impact_parameter(Mc2, mtc2, E, Z1, Z2, b) result(Etr) ! transferred energy from impact parameter
   real(8) Etr  ! transfered energy
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of projectile and scattering center
   real(8), intent(in) :: E         ! [eV] projectile energy
   real(8), intent(in) :: Z1, Z2    ! [electron charge] Charges of the projectile and scattering center
   real(8), intent(in) :: b         ! [A] impact parameter
   real(8) :: Mr, Mr1, costheta
   Mr = Mc2/mtc2    ! mass ratio
   Mr1 = (Mr+1.0d0)
   
   costheta = Mr1/sqrt(E*b/(g_Ry*g_a0*Z1*Z2)**2 + Mr1**2)
   
   ! Relativistic generalization of Eq.(10) [1]:
   !Etr = 2.0d0*E*(E/mtc2 + 2.0d0*Mr)*Mr1/( (2.0d0*E + Mr1)*(Mr1 + b*E/(g_a0*Z1*Z2*g_Ry)) )
   !Etr = 2.0d0*mtc2*(E*E - Mc2*Mc2)*costheta/((mtc2 + E)**2 - (E-Mc2)**2*costheta)
   Etr = W_from_cos_theta(Mc2, mtc2, E, costheta)    ! below
end function W_from_impact_parameter


pure function W_from_cos_theta(Mc2, mtc2, E, cos_theta) result(W)
   real(8) W    ! [eV]
   real(8), intent(in) :: Mc2, mtc2, E, cos_theta
   real(8) :: cos2, sin2
   cos2 = cos_theta**2
   sin2 = 1.0d0-cos2
   !W = (E+Mc2)*sin2 + mtc2 + cos_theta*sqrt(abs(mtc2**2 - Mc2**2*sin2))
   !W = W * E*(E+2.0d0*Mc2)/( (E+Mc2+mtc2)**2 - E*(E+2*Mc2)*cos2 )
   W = ((E+Mc2)*sin2 + mtc2*(1.0d0-cos_theta))*E*(E+2.0d0*Mc2)/((E+mtc2)**2-E*(E+2.0d0*Mc2)*cos2)
end function W_from_cos_theta



pure function Impact_parameter_from_W(Mc2, mtc2, E, Z1, Z2, W) result(b) ! transferred energy from impact parameter
   real(8) b   ! [A]
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of projectile and scattering center
   real(8), intent(in) :: E         ! [eV] projectile energy
   real(8), intent(in) :: Z1, Z2    ! [electron charge] Charges of the projectile and scattering center
   real(8), intent(in) :: W         ! transfered energy
   real(8) :: Mr, Mr1, A
   Mr = Mc2/mtc2    ! mass ratio
   Mr1 = (Mr+1.0d0)**2
   A = 2.0d0*E*(E/mtc2 + 2.0d0*Mr)*Mr1/(2.0d0*E + Mr1)
   ! Inverse of relativistic generalization of Eq.(10) [1]:
   b = g_a0*Z1*Z2*g_Ry/E*(A/W - Mr1)
end function Impact_parameter_from_W
 
 
 
function find_valence_hole_mass(numpar, DOS, E_DOS) result(Mh)
   real(8) Mh   ! [kg] holes mass
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Density_of_states), intent(in) :: DOS    ! Valence and/or conduction band DOS of the material
   real(8), intent(in) :: E_DOS ! energy of the valence hole counted from top of VB [eV]
   real(8) :: eps, mass_temp
   ! Find holes mass:
   eps = 1.0d-6 ! margin within which mass equals to zero
   if (numpar%H_m_eff > eps) then ! a constant coefficient * me
      Mh = numpar%H_m_eff * g_me
   elseif (abs(numpar%H_m_eff) < eps) then ! equals to user-provided mass (electron, ion, etc.)
      Mh = g_me    ! free electron mass
   else  ! effective mass from DOS
      call interpolate_data_single(DOS%E, DOS%Eff_m, (DOS%E_VB_top - E_DOS), mass_temp) ! module "Little_subroutines"
      Mh = mass_temp * g_me
   endif
end function find_valence_hole_mass
 

pure subroutine get_ranges_from_Se(E, Se, E_start, ranges)
   real(8), dimension(:), intent(in) :: E   ! [eV] energy grid
   real(8), dimension(:), intent(in) :: Se  ! [eV/A] stopping power
   real(8), intent(in) :: E_start   ! [eV] lower limit of integration of Se for definition of range
   real(8), dimension(:), allocatable, intent(inout) :: ranges  ! [A]
   integer :: Nsiz, i
   real(8) :: eps
   eps = 1.0d-10    ! margin or error around zero
   Nsiz = size(E)   ! size of the grid
   if (.not.allocateD(ranges)) allocate(ranges(Nsiz))   ! make sure ranges array is allocated
   ranges = 0.0d0   ! to start with
   do i = 1, Nsiz
      if ((E(i) >= E_start) .and. (Se(i) > eps)) then ! non-zero stopping power
         if (i <= 1) then
            ranges(i) = 0.0d0 + 0.5d0/Se(i) * (E(i+1) - E(i))
         elseif (Se(i-1) < eps) then
            ranges(i) = 0.0d0 + 0.5d0/Se(i) * (E(i+1) - E(i))
         else
            ranges(i) = ranges(i-1) + 0.5d0*(1.0d0/Se(i) + 1.0d0/Se(i-1)) * (E(i) - E(i-1))
         endif
      else
         ranges(i) = 0.0d0  ! just to exclude those ranges, for which there is no stopping
      endif
   enddo
end subroutine get_ranges_from_Se
 
 

pure function get_total_CS_from_partial(sigma, percentage, masking) result(total_CS)
   real(8) total_CS ! total cross section from the partial ones
   real(8), dimension(:), intent(in) :: sigma   ! partial cross sections
   real(8), dimension(:), intent(in) :: percentage  ! weights of the corresponding cross sections
   logical, dimension(:), intent(in), optional :: masking     ! mask array which cross sections to include here
   real(8) :: norm_perc
   if (present(masking)) then
      norm_perc = SUM(percentage(:), MASK = masking)
      total_CS = SUM(sigma(:)*percentage(:), MASK = masking)/norm_perc
   else
      norm_perc = SUM(percentage(:))
      total_CS = SUM(sigma(:)*percentage(:))/norm_perc
   endif
end function get_total_CS_from_partial


function total_CS_from_chennels(Ekin, E_array1, CS_array1, E_array2, CS_array2, E_array3, CS_array3, E_array4, CS_array4, E_array5, CS_array5) result(total_CS)
   real(8) total_CS ! total cross section from the partial ones
   real(8), intent(in) :: Ekin    ! energy of the particle
   real(8), dimension(:), intent(in) :: E_array1 ! energy array
   real(8), dimension(:), intent(in) :: CS_array1   ! partial cross section
   real(8), dimension(:), intent(in), optional :: E_array2, E_array3, E_array4, E_array5 ! energy arrays
   real(8), dimension(:), intent(in), optional :: CS_array2, CS_array3, CS_array4, CS_array5  ! possible other partial arrays
   !---------------------------------------------
   real(8) :: CS
   total_CS = 0.0d0 ! to start with
!    call interpolate_data_single(E_array1, CS_array1, Ekin, CS, x0=0.0d0, y0=0.0d0) ! module "Little_subroutines"
   call interpolate_data_single(E_array1, CS_array1, Ekin, CS) ! module "Little_subroutines"
   total_CS = total_CS + CS ! [A^2]
   
   if (present(E_array2) .and. present(CS_array2)) then
!       call interpolate_data_single(E_array2, CS_array2, Ekin, CS, x0=0.0d0, y0=0.0d0) ! module "Little_subroutines"
      call interpolate_data_single(E_array2, CS_array2, Ekin, CS) ! module "Little_subroutines"
      total_CS = total_CS + CS ! [A^2]   
   endif
   
   if (present(E_array3) .and. present(CS_array3)) then
!       call interpolate_data_single(E_array3, CS_array3, Ekin, CS, x0=0.0d0, y0=0.0d0) ! module "Little_subroutines"
      call interpolate_data_single(E_array3, CS_array3, Ekin, CS) ! module "Little_subroutines"
      total_CS = total_CS + CS ! [A^2]   
   endif
   
   if (present(E_array4) .and. present(CS_array4)) then
!       call interpolate_data_single(E_array4, CS_array4, Ekin, CS, x0=0.0d0, y0=0.0d0) ! module "Little_subroutines"
      call interpolate_data_single(E_array4, CS_array4, Ekin, CS) ! module "Little_subroutines"
      total_CS = total_CS + CS ! [A^2]   
   endif
   
   if (present(E_array5) .and. present(CS_array5)) then
!       call interpolate_data_single(E_array5, CS_array5, Ekin, CS, x0=0.0d0, y0=0.0d0) ! module "Little_subroutines"
      call interpolate_data_single(E_array5, CS_array5, Ekin, CS) ! module "Little_subroutines"
      total_CS = total_CS + CS ! [A^2]   
   endif
end function total_CS_from_chennels



subroutine find_type_of_scattering(i_type, Ekin, E_array1, CS_array1, E_array2, CS_array2, E_array3, CS_array3, E_array4, CS_array4, E_array5, CS_array5)
   integer, intent(out) :: i_type   !  index of the type of scattering: which cross section is realized
   real(8), intent(in) :: Ekin    ! energy of the particle
   real(8), dimension(:), intent(in) :: E_array1 ! energy array
   real(8), dimension(:), intent(in) :: CS_array1   ! partial cross section
   real(8), dimension(:), intent(in), optional :: E_array2, E_array3, E_array4, E_array5 ! energy arrays
   real(8), dimension(:), intent(in), optional :: CS_array2, CS_array3, CS_array4, CS_array5  ! possible other partial arrays
   !---------------------------------------------
   real(8) :: total_CS, CS_cur, CS(5), RN
   integer :: i
   logical :: more_than_one, array_present_2, array_present_3, array_present_4, array_present_5
   array_present_2 = (present(E_array2) .and. present(CS_array2))
   array_present_3 = (present(E_array3) .and. present(CS_array3))
   array_present_4 = (present(E_array4) .and. present(CS_array4))
   array_present_5 = (present(E_array5) .and. present(CS_array5))

   more_than_one = (array_present_2 .or. array_present_3 .or. array_present_4 .or. array_present_5)
   
   if (.not.more_than_one) then ! there is only one CS, nothing to chose from
      i_type = 1
   else ! chose among the cross sections   
      ! Find all cross sections involved:
      CS(:) = 0.0d0 ! to start with
      total_CS = 0.0d0 ! to start with
      call interpolate_data_single(E_array1, CS_array1, Ekin, CS_cur) ! module "Little_subroutines"
      total_CS = total_CS + CS_cur
      CS(1) = total_CS
   
      if (array_present_2) then
         call interpolate_data_single(E_array2, CS_array2, Ekin, CS_cur) ! module "Little_subroutines"
         total_CS = total_CS + CS_cur
         CS(2) = total_CS
      endif
   
      if (array_present_3) then
         call interpolate_data_single(E_array3, CS_array3, Ekin, CS_cur) ! module "Little_subroutines"
         total_CS = total_CS + CS_cur
         CS(3) = total_CS 
      endif
   
      if (array_present_4) then
         call interpolate_data_single(E_array4, CS_array4, Ekin, CS_cur) ! module "Little_subroutines"
         total_CS = total_CS + CS_cur
         CS(4) = total_CS
      endif
   
      if (array_present_5) then
         call interpolate_data_single(E_array5, CS_array5, Ekin, CS_cur) ! module "Little_subroutines"
         total_CS = total_CS + CS_cur
         CS(5) = total_CS
      endif
      
      ! Normalize to get probability:
      if (total_CS > 1.0d-15) then  ! only if thereare non-zero CSs
         CS(:) = CS(:)/total_CS
         ! Sample a random number:
         call random_number(RN)    ! [0,1] intrinsic FORTRAN
         ! Get the cross index of the process from the sampled CSs:
         i = 1
         do while (RN > CS(i))
            i = i + 1
         enddo
         i_type = i
      else  ! if there are no non-zero CS, printout 1 as default
         i_type = 1
      endif
      
      if (total_CS < 1.0d-10) then
         print*, 'find_type_of_scattering', i_type , CS, 'Total:', total_CS
      endif
      
   endif
end subroutine find_type_of_scattering




pure function MFP_from_sigma(sigma, nat) result(lambda)
   real(8) lambda   ! [1/A]
   real(8), intent(in) :: sigma ! [A^2]
   real(8), intent(in) :: nat   ! [cm^-3]
   real(8) :: na
   if (sigma > 1.0d-24) then  ! finite MFP
      na = nat * 1.0d-24    ! [A^-3] converted from [cm^-3]
      lambda = 1.0d0 / (sigma * na)    ! [A] mean free path 
   else
      lambda = 1.0d30   ! infinity
   endif
end function MFP_from_sigma


pure function Time_from_MFP(MFP, V) result(T)
   real(8) T   ! [fs]
   real(8), intent(in) :: MFP ! [A] mean free path
   real(8), intent(in) :: V   ! [A/fs] velosity
   real(8) eps
   eps = 1.0d-6 ! minimal speed [A/fs]
   if ( (MFP > 1.0d20) .or. (V < eps) ) then  ! infinite MFP
      T = 1.0d25
   else
      T = MFP/V
   endif
end function Time_from_MFP


pure function temperature_factor(W, T) result (F)
   real(8) :: F
   real(8), intent(in) :: W ! [eV] transferred energy
   real(8), intent(in) :: T  ! [eV] target temperature
   real(8) :: eps
   eps = 1.0d-4
   if ((T < eps) .or. (W < 1.0d-10)) then    ! zero temperature or no energy
      F = 1.0d0
   else ! non-zero tempreature
      F = 1.0d0/(1.0d0 - exp(-W/T))
   endif
end function temperature_factor



pure function plasmon_energy(Nel, At_Dens, Egap) result(Epl)
   real(8) :: Epl   ! [eV] plasmon energy
   real(8), intent(in) :: Nel           ! valence electron density (per atom)
   real(8), intent(in) :: At_Dens   ! atomic density [1/cm^2]
   real(8), intent(in) :: Egap        ! band gap [eV]
   Epl = sqrt(Nel*At_Dens*1d6*g_h*g_h/(g_me*g_e0) + Egap*Egap)    ! Maximal energy of plasmons [eV]
end function plasmon_energy


pure function remap_energy_step(E, n, G, E0, useold) result(dE)
   real(8) :: dE
   real(8), intent(in) :: E ! [eV] integration point in energy
   integer, intent(in) :: n ! effective number of integration points
   real(8), dimension(:), intent(in), optional :: E0, G ! [eV] E0 and Gamma parameters of CDF
   logical, intent(in), optional :: useold  ! use the old (TREKIS-3) mapping instead of the new one
   logical :: old_mapping
   integer :: i, Nsiz
   real(8) :: E_min, E_cur
   
   ! Check if the user wants to use old mapping or new one:
   old_mapping = .false.
   if (present(useold)) then
      if (useold) then
         old_mapping = .true.
      else
         old_mapping = .false.
      endif
   endif
   
   if (old_mapping) then
      ! the integration step is chosen by this mapping as to provide small steps at low energies
      ! (where there are peaks in CDF), and larger steps at higher energies along smooth decreesing tails:
      dE = (1.0d0/(E+1.0d0) + E)/dble(n)
   else     ! use new mapping
      ! the integration step is chosen small around CDF peaks, but larger away from any of them:
      Nsiz = size(E0)  ! number of CDF oscillators
      E_min = 1d20     ! to start with
      do i = 1, Nsiz
         E_cur = max( (E - E0(i)), G(i) )
         if (E_min > E_cur) E_min = E_cur
!          print*, E, E - E0(i)
      enddo
      dE = E_min/dble(n) 
   endif
end function remap_energy_step


pure function remap_dq(hq, n, W, E0, G) result(dq)
   real(8) dq   ! [eV]
   real(8), intent(in) :: hq    ! [eV] integration point in momentum space (hq)
   real(8), intent(in) :: W     ! [eV] integration point in energy space (hw)
   integer, intent(in) :: n     ! effective number of integration points
   real(8), dimension(:), intent(in) :: E0, G   ! CDF parameters
   integer :: Nsiz, i
   real(8), dimension(size(E0)) :: q_center
   real(8) :: q_min, q_cur
   
   ! Position of the center of the Drude oscillator:
   q_center(:) = W - E0(:)
   
   ! Set the integration step in accordance with the Q:
   Nsiz = size(E0)  ! number of CDF oscillators
   q_min = hq       ! to start with
   do i = 1, Nsiz
      q_cur = max( (hq - q_center(i)), G(i) )
      if (q_min > q_cur) q_min = q_cur
   enddo
   dq = q_min/dble(n)
end function remap_dq



end module CS_general_tools
