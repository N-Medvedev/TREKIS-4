! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019-2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains Monte Carlo routines to simulate positron events
module MC_positron
use Universal_constants
use Objects
use Geometries !,only: m_tollerance_eps
use Little_subroutines, only: interpolate_data_single
use Relativity, only: velosity_from_kinetic_energy, beta_factor, rest_energy, kinetic_energy_from_velosity
use MC_general_tools, only: flight_time_to_boundary, out_of_simulation_box, particle_cross_box_boundary, cos_theta_from_W, &
                            cos_recoil_from_W, check_shell, check_element, total_electrons_in_target, choose_random_shell, &
                            sample_phi, sample_theta, check_phi, transfered_E_from_theta, &
                            deflect_velosity, get_electron_flight_time, get_hole_flight_time, get_photon_flight_time, &
                            get_positron_flight_time, remove_particle_from_MC_array, extend_MC_array, find_the_target, &
                            electron_transmission_probability, reflection_from_surface, define_normal_to_surface, &
                            add_energy_into_MD_array
use CS_general_tools, only: total_CS_from_chennels, MFP_from_sigma, Time_from_MFP, find_type_of_scattering, find_valence_hole_mass
use CS_positrons_inelastic, only: get_inelastic_energy_transfer_pos
use CS_positrons_elastic, only: get_pos_elastic_CS
use CS_electrons_Bremsstrahlung, only: sample_Bremsstrahlung_electron_angles, sample_Bremsstrahlung_theta
use CS_electrons_elastic, only: Mott_sample_mu
use CS_electrons_inelastic, only: CDF_total_CS_nonrel
use CS_positrons_Bremsstrahlung, only: get_energy_transfer_Bremsstrahlung_pos
use CS_positrons_annihilation, only: theta_annih, theta_plus_annih, get_Heitler_CS, ksi_min, get_Heitler_intCS, Eph_from_ksi_annih
use CS_integration_limits, only: W_max, find_Wmax_equal_Wmin
use Dealing_with_DOS, only: select_energy_dos
use SHI_charge_state, only: Equilibrium_charge_SHI

implicit none

 contains

subroutine MC_positron_event(used_target, numpar, N_p, Prtcl, NOP, MC, MD_supce, E_p_at, E_e, E_h)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(inout) :: N_p   ! number of positrons
   type(Positron), dimension(:), allocatable, intent(inout) :: Prtcl  ! the positron to perform some event
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout) :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_p_at, E_e, E_h ! data to pass to MD later
   !----------------------------------
   integer :: i_type   ! type of event: -1=box crossing; 0=boundary; 1=inelastic; 2=elastic; 3=Bremsstrahlung
   
   ! 0) Advance the last-step parameters:
   Prtcl(NOP)%R0(:) = Prtcl(NOP)%R(:)
   Prtcl(NOP)%V0(:) = Prtcl(NOP)%V(:)
   
   ! 1) Move the particle to the end-point:
   Prtcl(NOP)%R(:) = Prtcl(NOP)%R0(:) + (Prtcl(NOP)%V(:)) * (Prtcl(NOP)%ti - Prtcl(NOP)%t0) !* g_ms2Afs
   Prtcl(NOP)%t0 = Prtcl(NOP)%ti
   
   ! 2) Choose a positron event:
   call find_type_of_positron_event(used_target, numpar, Prtcl(NOP), i_type)  ! below
   
   ! 3) Perform the event according to the chosen type:
   select case(i_type)
   case (0) ! target boundary crossing
      call event_positron_target_boundary(used_target, numpar, Prtcl(NOP))   ! below
      
   case (1) ! inelastic
      call event_positron_inelastic(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)   ! below
      
   case (2) ! elastic
      call event_positron_elastic(used_target, numpar, MC, NOP, MD_supce, E_p_at)   ! below
      
   case (3) ! Bremsstrahlung
      call event_positron_Bremsstrahlung(used_target, numpar, MC, NOP)   ! below
      
   case(4)  ! annihilation
      call event_positron_annihilation(N_p, used_target, numpar, MC, NOP, MD_supce, E_e, E_h)   ! below
      
   case default ! simulation box boundary crossing
      call particle_cross_box_boundary(used_target, numpar, N_p, NOP, Prtcl)    ! module "MC_general_tools"
      
   endselect
end subroutine MC_positron_event



!ССССССССССССССССССССССССССССССССССССССССССС
! Positron crossing target boundary:
subroutine event_positron_target_boundary(used_target, numpar, Prtcl)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Positron), intent(inout) :: Prtcl        ! positron as an object
   !------------------------------------------------------
   real(8), dimension(3)  :: R_shift, norm_to_surf, V_perp, V
   real(8) :: Vnorm, Ekin_norm, T, RN, Vabs, Z
   integer :: target_ind
   logical :: transmit
   type(Emission_barrier), pointer :: Em_Barr
   !------------------------------------------------------
   
   ! Find whether an positron crosses the boundary or reflects back:
   ! 1) Kinetic energy towards crossing the boundary:
   call define_normal_to_surface(used_target,  Prtcl, norm_to_surf, 'event_positron_target_boundary')    ! module "MC_general_tools"
   
   ! 2) Find the velosity projection onto the normal to the surface:
   Vnorm = scalar_projection(Prtcl%V, norm_to_surf)     ! module "Geometries"
   Vnorm = Vnorm * g_Afs2ms ! [fs/A] -> [m/s]
   
   ! 3) Get the normal component of the kinetic energy:
   Ekin_norm = kinetic_energy_from_velosity(Vnorm, g_me)    ! module "Relativity"
   
   ! 4) Get the probability of boundary crossing:
   Em_Barr => used_target%Material(Prtcl%in_target)%Surface_barrier
   T = electron_transmission_probability(Em_Barr%Work_func, Em_Barr%gamma, Em_Barr%E1, Ekin_norm)     ! module "MC_general_tools"

   ! 5) Sample transmission vs. vs. reflection:
   transmit = .true.   ! to start with, assume transmission
   RN = 0.0d0
   if (T < 0.99999999d0) then   ! reflection is plausible:
      call random_number(RN)
      if (RN > T) transmit = .false.    ! it will be reflection
   endif

   ! 6) Model reflection or transmission:
   if (transmit) then   ! transmission into another target:
      ! Find into which target this positron enters:
      ! Find which target's boundary the particle is crossing:
      Vabs = SQRT( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
      R_shift = m_tollerance_eps * Prtcl%V(:)/Vabs     ! to place particle inside of the material 
      ! Update particle's material index according to the new material it enters:
      call find_the_target(used_target, Prtcl, R_shift) ! module "MC_general_tools"

      ! 6.a) Shift positron just across the border:
      Prtcl%R(:) = Prtcl%R(:) + m_tollerance_eps * Prtcl%V(:)/Vabs

      ! 6.b) Get the next flight inside the new target (transmitted) or old target (reflected):
      call get_positron_flight_time(used_target, numpar, Prtcl)  ! module "MC_general_tools"
   else ! reflection back
      ! Change velosity according to reflection from the surface with given normal:
      call reflection_from_surface(Prtcl%V, norm_to_surf)   ! module "MC_general_tools"

      ! 6.c) Shift positron just across the border:
      Vabs = SQRT( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
      Prtcl%R(:) = Prtcl%R(:) + m_tollerance_eps * Prtcl%V(:)/Vabs

      ! 6.d) Get the next flight inside the new target (transmitted) or old target (reflected):
      call get_positron_flight_time(used_target, numpar, Prtcl, .true.)  ! module "MC_general_tools"
   endif

   nullify(Em_Barr)
end subroutine event_positron_target_boundary



!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
subroutine event_positron_annihilation(N_p, used_target, numpar, MC, NOP, MD_supce, E_e, E_h)   ! below
   integer, intent(inout) :: N_p   ! number of positrons
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h   ! data to pass to MD later
   !------------------------------------------------------
   type(Positron), pointer :: Prtcl  ! the positron to perform some event
   type(Target_atoms), pointer :: matter
   real(8) :: theta, phi, theta_2, phi_2, polar_angle_h, azimuth_angle_h
   real(8) :: Eph, Eph2, Mh, V_h_abs, V0(3), E_DOS, dE
   real(8) :: Vtot, Vtot0, a0(3), a(3), V_ph(3), V_h(3)
   integer :: KOA, NSH
   
   ! 0) Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Positrons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! 1) Sample first photon emission energy:
   call  find_annihilation_photon_energy(Prtcl%Ekin, Eph) ! below
   
   ! 2) Find the angles of the lower-energy photon emission:
   call  set_photon_annih_angles(Prtcl%Ekin, Eph, theta, phi) ! below
   
    ! 3) Get the cosines of the positron velosity:
   Vtot = sqrt( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
   a0(:) = Prtcl%V(:)/Vtot
   Vtot0 = Vtot ! save it for later setting photon parameters
   
   ! 4) Create first emitted photon:
   MC%N_ph = MC%N_ph + 1
   ! 4a) Get the cosines of the emitting photon velosity:
!    call deflect_velosity(a0(1), a0(2), a0(3), phi, theta, a(1), a(2), a(3))  ! module "MC_general_tools"
   call deflect_velosity(a0(1), a0(2), a0(3), theta, phi, a(1), a(2), a(3))  ! module "MC_general_tools"
!    V_ph(:) = g_cvel * a(:)   ! derection of photon motion [m/s]
   V_ph(:) =g_c_Afs * a(:)   ! derection of photon motion [A/fs]
   ! 4b) Save the photon data into the particle properties:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_ph > size(MC%MC_Photons)) call extend_MC_array(MC%MC_Photons)    ! module "MC_general_tools"
   call make_new_particle(MC%MC_Photons(MC%N_ph), Ekin=Eph, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_ph)    ! module "Objects"
   ! 4c) Define the next photon scattering event:
   call get_photon_flight_time(used_target, numpar, MC%MC_Photons(MC%N_ph))  ! module "MC_general_tools"
   ! Save the data for the last step:
   MC%MC_Photons(MC%N_ph)%R0(:) = MC%MC_Photons(MC%N_ph)%R(:)
   MC%MC_Photons(MC%N_ph)%V0(:) = MC%MC_Photons(MC%N_ph)%V(:)
   
   ! 5) Find the parameters of the hole that is left after electron disappearence:
   ! 5a) Find in which shell this hole appears:
   call choose_random_shell(matter, KOA, NSH)    ! module "MC_general_tools"
   
   ! 6) Energy of the second photon:
   dE = Prtcl%Ekin + 2.0d0*g_me_eV - Eph  ! energy to be given to the photon (subtract the Ip later)
   if (matter%Elements(KOA)%valent(NSH)) then   ! valence hole
      ! Define the energy within DOS where the hole is sitting:
      call select_energy_DOS(matter%DOS%E, matter%DOS%DOS, matter%Integral_DOS_fe, &
                                           matter%DOS%Egap, dE, matter%DOS%alpha_CB, E_DOS)   ! module "Dealing_with_DOS"
      Eph2 = dE - (E_DOS - matter%DOS%Egap)
   else ! core hole
      Eph2 = dE - matter%Elements(KOA)%Ip(NSH)
   endif

   ! 6a) Set the angles of the second photon:
   theta_2 = theta_plus_annih(Prtcl%Ekin, Eph)  ! module "CS_positrons_annihilation"
   phi_2 = g_Pi - phi
   ! Make sure phi_2 is within [0:2Pi]
   call check_phi(phi_2)    ! module "MC_general_tools"
   ! 6b) Create first emitted photon:
   MC%N_ph = MC%N_ph + 1
   ! 6c) Get the cosines of the emitting photon velosity:
!    call deflect_velosity(a0(1), a0(2), a0(3), phi_2, theta_2, a(1), a(2), a(3))  ! module "MC_general_tools"
   call deflect_velosity(a0(1), a0(2), a0(3), theta_2, phi_2, a(1), a(2), a(3))  ! module "MC_general_tools"
!    V_ph(:) = g_cvel * a(:)   ! derection of photon motion [m/s]
   V_ph(:) = g_c_Afs * a(:)   ! derection of photon motion [A/fs]
   ! 6d) Save the photon data into the particle properties:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_ph > size(MC%MC_Photons)) call extend_MC_array(MC%MC_Photons)   ! module "MC_general_tools"
   call make_new_particle(MC%MC_Photons(MC%N_ph), Ekin=Eph2, t0=Prtcl%ti, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_ph)    ! module "Objects"
   ! 6e) Define the next photon scattering event:
   call get_photon_flight_time(used_target, numpar, MC%MC_Photons(MC%N_ph))  ! module "MC_general_tools"
   ! Save the data for the last step:
   MC%MC_Photons(MC%N_ph)%R0(:) = MC%MC_Photons(MC%N_ph)%R(:)
   MC%MC_Photons(MC%N_ph)%V0(:) = MC%MC_Photons(MC%N_ph)%V(:)
   
   ! 7) Create a new hole after electron removal due to annihilation:
   MC%N_h = MC%N_h + 1
   ! for valent hole, sample its velosity:
   if (matter%Elements(KOA)%valent(NSH)) then
      ! Find hole's mass:
      ! Knowing the energy, define the mass accordingly:
      Mh =  find_valence_hole_mass(numpar, matter%DOS, abs(E_DOS)) ! module "MC_general_tools"
      ! Set the velosity of the new hole:
      V_h_abs = velosity_from_kinetic_energy(abs(E_DOS), Mh)    ! [A/fs] module "Relativity"
      polar_angle_h = sample_phi() ! module "MC_general_tools"
      azimuth_angle_h = sample_theta()  ! module "MC_general_tools"
      V0(1) = 1.0d0
      V0(2:3) = 0.0d0
      ! Get the holes velosity in the sampled direction:
!       call  deflect_velosity(V0(1), V0(2), V0(3), polar_angle_h, azimuth_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      call  deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle_h, polar_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      V_h(:) = V_h_abs * V_h(:)   ! set also the absolute value
   else
      V_h(:) = 0.0d0   ! core holes don't fly
   endif
   ! Save parameters of the new particles into array:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_h > size(MC%MC_Holes)) call extend_MC_array(MC%MC_Holes)    ! module "MC_general_tools"
   call make_new_particle(MC%MC_Holes(MC%N_h), Ekin=E_DOS, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_h, &
                                      KOA = KOA, Sh = NSH, valent = matter%Elements(KOA)%valent(NSH) )    ! module "Objects"

   ! 8) Get new holes time of the next event:
   call get_hole_flight_time(used_target, numpar, MC%MC_Holes(MC%N_h), MD_supce, E_h)  ! module "MC_general_tools"
   
   ! 9) Positron disappears:
   call remove_particle_from_MC_array(N_p, NOP, MC%MC_Positrons)   ! module "MC_general_tools"
   
   nullify(Prtcl, matter)
end subroutine event_positron_annihilation



subroutine find_annihilation_photon_energy(Epos, Eph)
   real(8), intent(in) :: Epos  ! [eV] positron energy
   real(8), intent(out) :: Eph   ! [eV] emitted photon energy
   real(8) :: RN, CS_sampled, CS_tot, CS_cur, CS_0
   real(8) :: ksi_low, ksi_high, ksi_cur
   real(8) :: eps
   
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-3
   
   ! Get total annihilation cross section:
   CS_tot = get_Heitler_CS(Epos)  ! module "CS_positrons_annihilation"
   
   ! Sample the cross section:
   call random_number(RN)
   CS_sampled = RN*CS_tot
   
   ! Find the minimal and maximal normalized energy:
   ksi_low = ksi_min(Epos)       ! module "CS_positrons_annihilation"
   ksi_high = 0.5d0 ! before Eq.(3.186) [1]
   
   ! Define the low integartion boundary:
   CS_0 = get_Heitler_intCS(Epos, ksi_low)    ! module "CS_positrons_annihilation"
   
   ! Start finding CS using bisection method:
   ksi_cur = (ksi_low + ksi_high)*0.5d0
   ! Total cross section for the given energy:
   CS_cur = get_Heitler_intCS(Epos, ksi_cur) - CS_0   ! module "CS_positrons_annihilation"
   ! Search by bisection method:
   do while (ABS(CS_cur - CS_sampled)/CS_sampled > eps)
      if (CS_cur > CS_sampled) then
         ksi_high = ksi_cur
      else
         ksi_low = ksi_cur
      endif
      ksi_cur = (ksi_low + ksi_high)/2.0d0
      if (abs(ksi_low - ksi_high) < eps) exit  ! precise enough
      CS_cur = get_Heitler_intCS(Epos, ksi_cur) - CS_0   ! module "CS_positrons_annihilation"
   enddo
   
   ! Convert normalized photon energy back into eV:
   Eph = Eph_from_ksi_annih(Epos, ksi_cur)  ! module "CS_positrons_annihilation"
end subroutine find_annihilation_photon_energy



subroutine set_photon_annih_angles(Epos, Eph, theta, phi)
   real(8), intent(in) :: Epos, Eph ! [eV] positron kinetic energy, and photon emission energy
   real(8), intent(out) :: theta, phi   ! angles of first photon emission
   theta =  theta_annih(Epos, Eph)  ! module "CS_positrons_annihilation"
   phi = sample_phi()   ! module "MC_general_tools"
end subroutine set_photon_annih_angles



!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
! Positron Bremsstrahlung emission subroutines:
subroutine event_positron_Bremsstrahlung(used_target, numpar, MC, NOP)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   !------------------------------------------------------
   type(Positron), pointer :: Prtcl  ! the positron to perform some event
   type(Target_atoms), pointer :: matter
   type(Atom_kind), pointer :: Element
   real(8) :: RN, phi, theta, CS_tot, CS_sampled, CS_sum, CS_elem, N_elem, elem_contrib, CS_cur, Vtot0, Vtot, a0(3), a(3), dE, V_ph(3)
   real(8) :: phi_ph, theta_ph, beta, A_coef, B_coef
   integer :: j, KOA, KOA1
   logical :: found_shl
   
   ! 0) Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Positrons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! 1) Select the element the positron scatters on:
   ! Get the total cross section:
   call interpolate_data_single(matter%Pos_Brems_total%E(:), matter%Pos_Brems_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
   ! Get the random number, to sample partial cross section:
   call random_number(RN)   ! intrinsic FORTRAN subroutine
   ! Realized CS:
   CS_sampled = RN*CS_tot
   ! Find which of the partial cross sections correspond to the sampled one:
   CS_sum = 0.0d0   ! to start with
   N_elem = dble(SUM(matter%Elements(:)%percentage))   ! number of elements in this compound material
   ! Scan through all the atoms' shells to find the one:
   found_shl = .false.  ! to start with
   KOA = matter%N_Elements  ! by default, it is the last element
   if (KOA > 1) then    ! find element
      ! (no need to check, because true by exclusion if all the other elements are not true):
      SRCH:do j =1, matter%N_Elements-1     ! for each element, expect for the last one
         Element => matter%Elements(j)	! all information about this element
         elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         ! For each shell, find its partial cross section:
         call interpolate_data_single(Element%Pos_Brems%E,  Element%Pos_Brems%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
         CS_sum = CS_sum + CS_cur*elem_contrib
         CS_elem = CS_cur  ! save for later use
         ! Check if this shell is the sampled one according to the integral CS:
         call check_element(CS_sum, CS_sampled, j, KOA1, found_shl)    ! module "MC_general_tools"
!          if (found_shl) exit SRCH  ! no need to continue the search if we already found the KOA
         if (found_shl) then
            KOA = KOA1
            exit SRCH  ! no need to continue the search if we already found the KOA
         endif
      enddo SRCH
   else ! only one element, no neet to search
      Element => matter%Elements(KOA)          ! all information about this element
      ! For each shell, find its partial cross section:
      call interpolate_data_single(Element%Pos_Brems%E,  Element%Pos_Brems%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
      CS_elem = CS_cur  ! save for later use
   endif
   
   if ( (KOA > size(matter%Elements(:))) .or. (KOA < 1) ) then
      print*, "Error 1 in event_positron_Bremsstrahlung:", CS_sampled, RN*CS_tot, CS_tot
   endif
   if (isnan(CS_elem) .or. (CS_elem <= 0.0d0)) then
      print*, "Error 2 in event_positron_Bremsstrahlung:", Prtcl%Ekin, CS_elem
   endif
   
   ! Now we know which element the positron scatters on:
   Element => matter%Elements(KOA)  ! all information about this element
   
   ! 2) Get the energy of emitted photon:
   dE = get_energy_transfer_Bremsstrahlung_pos(Prtcl%Ekin, Element, numpar, KOA, CS_elem)    ! module "CS_positrons_Bremsstrahlung"
   
   ! 3) Get the cosines of the positron velosity:
   Vtot = sqrt( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
   a0(:) = Prtcl%V(:)/Vtot
   Vtot0 = Vtot ! save it for later setting photon parameters
   
   ! 4) Get the positron deflection angles (assume the same as for electron):
   call sample_Bremsstrahlung_electron_angles(phi,theta)    ! module "CS_electrons_Bremsstrahlung"
   
   ! 5a) Get the deflection angle (angles of scattering of positron in the laboratory system):
!    call  deflect_velosity(a0(1), a0(2), a0(3), phi, theta, a(1), a(2), a(3))  ! module "MC_general_tools"
   call  deflect_velosity(a0(1), a0(2), a0(3), theta, phi, a(1), a(2), a(3))  ! module "MC_general_tools"
   ! 5b) Change positron energy accordingly:
   Prtcl%Ekin = Prtcl%Ekin - dE ! [eV]
   ! 5c) Change the velosity accordingly:
   Vtot = velosity_from_kinetic_energy(Prtcl%Ekin, g_me)  ! [A/fs] module "Relativity"
   ! Update velosity:
   Prtcl%V(:) = Vtot * a(:)  ! components of the positron velosity
   
   ! 6) Now update the positron parameters:
   call get_positron_flight_time(used_target, numpar, Prtcl)  ! module "MC_general_tools"
   
   ! 7) Create emitted photon:
   MC%N_ph = MC%N_ph + 1
   ! 7a) Get the photon emission angles:
   phi_ph = sample_phi()   ! module "MC_general_tools"
   beta = beta_factor(Vtot0, afs_in = .true.)    ! module "Relativity"
   ! coefficients we don't know yet [1], so set the default ones:
   A_coef = 1.0d0
   B_coef = 0.0d0
   ! Assume emission angle the same as for electron:
   theta_ph = sample_Bremsstrahlung_theta(beta, A_coef, B_coef)   ! module "CS_electrons_Bremsstrahlung"
   ! 7b) Get the cosines of the emitting positron velosity:
!    call deflect_velosity(a0(1), a0(2), a0(3), phi_ph, theta_ph, a(1), a(2), a(3))  ! module "MC_general_tools"
   call deflect_velosity(a0(1), a0(2), a0(3), theta_ph, phi_ph, a(1), a(2), a(3))  ! module "MC_general_tools"
!    V_ph(:) = g_cvel * a(:)   ! derection of photon motion [m/s]
   V_ph(:) = g_c_Afs * a(:)   ! derection of photon motion [A/fs]
   ! 7c) Save the photon data into the particle properties:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_ph > size(MC%MC_Photons)) call extend_MC_array(MC%MC_Photons)    ! module "MC_general_tools"
   call make_new_particle(MC%MC_Photons(MC%N_ph), Ekin=dE, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_ph)    ! module "Objects"
   ! 7d) Define the next photon scattering event:
   call get_photon_flight_time(used_target, numpar, MC%MC_Photons(MC%N_ph))  ! module "MC_general_tools"
   ! Save the data for the last step:
   MC%MC_Photons(MC%N_ph)%R0(:) = MC%MC_Photons(MC%N_ph)%R(:)
   MC%MC_Photons(MC%N_ph)%V0(:) = MC%MC_Photons(MC%N_ph)%V(:)
   
   nullify(Prtcl, matter, Element)
end subroutine event_positron_Bremsstrahlung


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Positron inelastic scattering subroutines:
subroutine event_positron_inelastic(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h   ! data to pass to MD later
   !------------------------------------------------------
   integer :: KOA, NSH
   real(8) :: Ekin, Ekin_new, dE, theta, phi, theta_r, phi_r, eps
   real(8) :: V0(3), V(3), V_tot, V_e_abs, V_e(3), E_DOS, a0(3), a(3)
   real(8) :: polar_angle_h, azimuth_angle_h, Mh, mass_temp, V_h_abs, V_h(3)
   type(Positron), pointer :: Prtcl  ! the positron to perform some event
   type(Target_atoms), pointer :: matter
!    type(Atom_kind), pointer :: Element
   logical :: valent
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Positrons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! 1) Which shell of which element is being ionized:
   call select_shell_ionization_pos(used_target, numpar, MC, NOP, KOA, NSH)  ! below
   
   ! 2) Sample transfered energy in the impact ionization event:
   ! Sample the transfered energy:
    if (KOA == 0) then   ! valence band:
      valent = .true.
      dE = get_inelastic_energy_transfer_pos(Prtcl%Ekin, matter, matter%DOS, numpar, KOA, NSH, matter%DOS%Egap)    ! module "CS_positrons_inelastic"
    else ! core shell:
      valent = .false.
      dE = get_inelastic_energy_transfer_pos(Prtcl%Ekin, matter, matter%DOS, numpar, KOA, NSH, matter%Elements(KOA)%Ip(NSH))    ! module "CS_positrons_inelastic"
   endif
   ! Just in case of numerical inaccuracies, make sure transfered energy is not above the electron energy:
   if (dE > Prtcl%Ekin) dE = Prtcl%Ekin
   
   ! 2a) Scattering angles according to the transferred energy:
   call sample_angles_inelastic_pos(Prtcl%Ekin, dE, theta, phi) ! below
   ! And the angles of emission of a new electron:
   call set_angles_of_new_electron_pos(Prtcl%Ekin, dE, phi, theta_r, phi_r)  ! below
   
   ! 3) Energy an electron is emitted with:
   if (KOA == 0) then   ! valence band:
      ! Sample from where within the valence band it is ionized according to DOS:
      call select_energy_DOS(matter%DOS%E, matter%DOS%DOS, matter%Integral_DOS_fe, &
                                           matter%DOS%Egap, dE, matter%DOS%alpha_CB, E_DOS)   ! module "Dealing_with_DOS"
      Ekin = dE + (E_DOS - matter%DOS%Egap)
   else ! core shell:
      Ekin = dE - matter%Elements(KOA)%Ip(NSH)
      E_DOS = 0.0d0
   endif
   ! Make sure the transfered energy is sufficient:
   if (Ekin < 0.0d0) then   ! insufficient energy is replaced with the minimal allowed:
      Ekin = 0.0d0
      if (KOA == 0) then   ! valence band:
         dE = -E_DOS + matter%DOS%Egap
      else
         dE = matter%Elements(KOA)%Ip(NSH)
      endif
   endif
   
   ! 4) Change the incident positron direction of motion:
   ! Get the cosines of the positron velosity:
   V_tot = sqrt(SUM(Prtcl%V(:)*Prtcl%V(:)))
   V0(:) = Prtcl%V(:)/V_tot
   ! Get the deflection angle (angles of emission of photoelectron in the laboratory system):
   call deflect_velosity(V0(1), V0(2), V0(3), theta, phi, V(1), V(2), V(3))  ! module "MC_general_tools"
   ! Get the NEW absolute electron velosity, to set its components according to the cosines:
   Prtcl%Ekin = Prtcl%Ekin - dE ! update positron energy
   V_e_abs = velosity_from_kinetic_energy(Prtcl%Ekin, g_me)     ! [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the positron velosity
   ! Update velosity:
   Prtcl%V(:) = V_e(:)
   
   ! Get the positron next sampled free path:
   call get_positron_flight_time(used_target, numpar, Prtcl)  ! module "MC_general_tools"
   
   ! 4a) Get the direction of the new electron:
   call deflect_velosity(V0(1), V0(2), V0(3), theta_r, phi_r, V(1), V(2), V(3))  ! module "MC_general_tools"
   V_e_abs = velosity_from_kinetic_energy(Ekin, g_me)     ! [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the electron velosity

   ! 5) Make a new electron:
   MC%N_e = MC%N_e + 1
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_e > size(MC%MC_Electrons)) call extend_MC_array(MC%MC_Electrons)   ! module "MC_general_tools"
   call make_new_particle(MC%MC_Electrons(MC%N_e), Ekin=Ekin, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_e)    ! module "Objects"
   
   ! 6) Get new electrons time of the next event:
   call get_electron_flight_time(used_target, numpar, MC%MC_Electrons(MC%N_e), MD_supce, E_e)  ! module "MC_general_tools"
   
   ! 7) Create a new hole after electron removal:
    MC%N_h = MC%N_h + 1
    ! for valent hole, sample its velosity:
    if (valent) then
       ! Find hole's mass:
       Mh =  find_valence_hole_mass(numpar, used_target%Material(Prtcl%in_target)%DOS, abs(E_DOS)) ! module "MC_general_tools"
       ! Set the velosity of the new hole:
       V_h_abs = velosity_from_kinetic_energy(abs(E_DOS), Mh)    ! [A/fs] module "Relativity"
       polar_angle_h = sample_phi() ! module "MC_general_tools"
       azimuth_angle_h = sample_theta()  ! module "MC_general_tools"
       V0(1) = 1.0d0
       V0(2:3) = 0.0d0
       ! Ger the holes velosity in the sampled direction:   
!        call  deflect_velosity(V0(1), V0(2), V0(3), polar_angle_h, azimuth_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
       call  deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle_h, polar_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
       V_h(:) = V_h_abs * V_h(:)   ! set also the absolute value
    else
       E_DOS = 0.0d0
       V_h(:) = 0.0d0   ! core holes don't fly
    endif
    ! Save parameters of the new particles into array:
    ! in case we have more particles than spaces in the array, extend the array:
    if (MC%N_h > size(MC%MC_Holes)) call extend_MC_array(MC%MC_Holes)   ! module "MC_general_tools"
    call make_new_particle(MC%MC_Holes(MC%N_h), Ekin=abs(E_DOS), t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_h, &
                                      KOA = KOA, Sh = NSH, valent = valent )    ! module "Objects"

   ! 8) Get new holes time of the next event:
   call get_hole_flight_time(used_target, numpar, MC%MC_Holes(MC%N_h), MD_supce, E_h)  ! module "MC_general_tools"
   
9999 nullify(Prtcl, matter)
end subroutine event_positron_inelastic


subroutine select_shell_ionization_pos(used_target, numpar, MC, NOP, KOA, NSH)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   integer, intent(out) :: KOA, NSH ! kind of atom, and the shell number, which absorbs the photon
   !------------------------------------------------------
   type(Target_atoms), pointer :: matter
   type(Positron), pointer :: Prtcl  ! the photon to perform some event
   type(Atom_kind), pointer :: Element
   real(8) :: RN, CS_tot, CS_sampled, CS_cur, CS_sum,  elem_contrib, N_elem
   integer :: i_mat, j, k
   logical :: found_shl, valence_done
   
   ! Starting default values (valence band):
   KOA = 0
   NSH = 0
   valence_done = .false.   ! to start, we did not calculate valence band yet
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Positrons(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! Get the total cross section:
   call interpolate_data_single(matter%Pos_inelastic_total%E(:), matter%Pos_inelastic_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
   ! Get the random number, to sample partial cross section:
   call random_number(RN)   ! intrinsic FORTRAN subroutine
   ! Realized CS:
   CS_sampled = RN*CS_tot

   ! Find which of the partial cross sections correspond to the sampled one:
   CS_sum = 0.0d0   ! to start with
   N_elem = dble(SUM(matter%Elements(:)%percentage))   ! number of elements in this compound material
   ! Scan through all the atoms' shells to find the one:
   SRCH:do j =1, matter%N_Elements	! for each element
      Element => matter%Elements(j)	! all information about this element
      elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
      ! Check the atomic CSs:
      do k = 1, Element%N_shl    ! for all shells in this elements
         
         VAL:if ((numpar%El_inelast /= 2) .and.  (Element%valent(k)) .and. (allocated(matter%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
            if (.not.valence_done) then ! add VB (but only once):
               call interpolate_data_single(matter%Pos_inelastic_valent%E,  matter%Pos_inelastic_valent%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
               CS_sum = CS_sum + CS_cur  ! band contribution
               valence_done = .true.     ! we added VB once (not do double-count it)
            endif
         else VAL
            call interpolate_data_single(Element%Pos_inelastic%E,  Element%Pos_inelastic%Per_shell(k,:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
            CS_sum = CS_sum + CS_cur*elem_contrib   ! elements contribution
         endif VAL
         
!          CS_sum = CS_sum + CS_cur*elem_contrib
         if (.not.Element%valent(k)) then ! core orbitals
            call check_shell(CS_sum, CS_sampled, j, k, KOA, NSH, found_shl)    ! module "MC_general_tools"
         else  ! valence orbitals
            if (allocated(matter%CDF_valence%A)) then    ! valence band
               call check_shell(CS_sum, CS_sampled, 0, 0, KOA, NSH, found_shl)    ! module "MC_general_tools"
            else    ! valence shells
               call check_shell(CS_sum, CS_sampled, j, k, KOA, NSH, found_shl)    ! module "MC_general_tools"
            endif
         endif ! (Element%valent(k)) 
         if (found_shl) exit SRCH  ! no need to continue the search if we already found the KOA nad NSH
      enddo ! k
   enddo SRCH !  j =1, N_elements

   ! Clean up the memory:
   nullify(matter, Prtcl, Element)
end subroutine select_shell_ionization_pos


subroutine sample_angles_inelastic_pos(Ekin, dE, theta, phi)
   real(8), intent(in) :: Ekin, dE  ! [eV] incident and transfered energy
   real(8), intent(out) :: theta, phi   ! scattering angles
   real(8) :: mu
   ! Assume uniform scattering probability into 2Pi for phi:
   phi = sample_phi()    ! module "MC_general_tools"
   ! Assume absolutely elastic scattering for theta (no recoil):
   mu = cos_theta_from_W(Ekin, dE, g_me, g_me)    ! module "MC_general_tools"
   theta = acos(mu)
end subroutine sample_angles_inelastic_pos


pure subroutine set_angles_of_new_electron_pos(Ekin, dE, phi, theta_r, phi_r)  
   real(8), intent(in) :: Ekin, dE  ! [eV] incident and transfered energy
   real(8), intent(in) :: phi   ! angle of deflection of the original particle
   real(8), intent(out) :: theta_r, phi_r   ! emission angles of a new electron
   real(8) :: mu
   ! Assume uniform scattering probability into 2Pi for phi:
   phi_r = phi + g_Pi
   call check_phi(phi_r) ! module "MC_general_tools"
   ! Assume absolutely elastic scattering for theta:
   mu = cos_recoil_from_W(Ekin, dE, g_me, g_me)    ! module "MC_general_tools"
   theta_r = acos(mu)
end subroutine set_angles_of_new_electron_pos


!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! Elastic scattering energy transfer
subroutine event_positron_elastic(used_target, numpar, MC, NOP, MD_supce, E_p_at)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_p_at ! data to pass to MD later
   !------------------------------------------------------
   type(Positron), pointer :: Prtcl  ! the electron to perform some event
   type(Target_atoms), pointer :: matter
   type(Atom_kind), pointer :: Element
   real(8) :: RN, phi, theta, CS_tot, CS_sampled, CS_sum, N_elem, elem_contrib, CS_cur, Vtot, a0(3), a(3), dE, Ekin
   integer :: j, KOA, KOA1
   logical :: found_at
   
   ! 0) Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Positrons(NOP)
   matter => used_target%Material(Prtcl%in_target)

   ! 1) Set the deflection angles:
   ! Sample azimuthal angle randomly within 2Pi:
   phi = sample_phi()   ! module "MC_general_tools"
   
   ! 2,3) Find out transfered energy and deflecting angle:
   call elastic_energy_transfer_p(numpar, Prtcl, matter, KOA, dE, theta)   ! below
   ! Just in case, check that the transfered energy is not too high:
   if (dE > Prtcl%Ekin) dE = Prtcl%Ekin ! [eV]
   
   ! 4) Get the cosines of the electron velosity:
   Vtot = sqrt( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
   a0(:) = Prtcl%V(:)/Vtot
   
   ! 5a) Get the deflection angle (angles of scattering of electron in the laboratory system):
!    call deflect_velosity(a0(1), a0(2), a0(3), phi, theta, a(1), a(2), a(3))  ! module "MC_general_tools"
   call deflect_velosity(a0(1), a0(2), a0(3), theta, phi, a(1), a(2), a(3))  ! module "MC_general_tools"
   
   Ekin = Prtcl%Ekin    ! save for testing
   Prtcl%Ekin = Prtcl%Ekin - dE ! [eV]
   
   ! 5c) Change the velosity accordingly:
   Vtot = velosity_from_kinetic_energy(Prtcl%Ekin, g_me)  ! [A/fs] module "Relativity"
   ! Update velosity:
   Prtcl%V(:) = Vtot * a(:)  ! components of the electron velosity
   
   ! 6) Now update the positron parameters:
   call get_positron_flight_time(used_target, numpar, Prtcl)  ! module "MC_general_tools"
   
   ! 7) Save the atomic in this collision parameters:
   MC%N_at_nrg = MC%N_at_nrg + 1
   if (MC%N_at_nrg > size(MC%MC_Atoms_events)) call extend_MC_array(MC%MC_Atoms_events)    ! module "MC_general_tools"
   ! Save the parameters of this collision (ONLY ENERGY TRANSFER IS SAVED, MOMENTUM NOT DONE YET!)
   call make_new_particle(MC%MC_Atoms_events(MC%N_at_nrg), Ekin=dE, t0=Prtcl%t0,  KOA = KOA, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=a)    ! module "Objects"
   
   ! Save the energy to pass to MD module, if needed:
   if (numpar%DO_MD) then   ! if user requested MD at all
      call add_energy_into_MD_array(numpar, dE, Prtcl%R, MD_supce, E_p_at)   ! module "MC_general_tools"
   endif

   nullify(Prtcl, matter, Element)
end subroutine event_positron_elastic



subroutine elastic_energy_transfer_p(numpar, Prtcl, matter, KOA, dE, theta)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Positron), intent(in) :: Prtcl  ! the electron to perform some event
   type(Target_atoms), intent(in), target :: matter ! material parameters
   integer, intent(out) :: KOA    ! kind of atom
   real(8), intent(out) :: dE, theta    ! [eV] transfered energy, and scattering angle
   !--------------------------------
   real(8) :: CS_tot, RN, CS_sampled, CS_sum, N_elem, elem_contrib, CS_cur, Eeq, max_E0, Zeff
   real(8) :: E_left, E_right, E_cur, eps, mu, hw_phonon, mtc2, Mc2, Se1
   logical :: found_at
   integer :: j, KOA1, Pos_elastic
   type(Atom_kind), pointer :: Element
   
   select case (numpar%Pos_elastic)     ! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   case (1,5)  ! CDF
      eps = 1.0d-3
      E_left = 0.0d0    ! [eV] minimal transferred energy
!       E_right = matter%E_debye    ! [eV] maximal transferred energy
      hw_phonon = maxval( matter%CDF_phonon%E0(:) + matter%CDF_phonon%Gamma(:) )
      mtc2 = rest_energy(matter%Mean_Mass)   ! target rest energy; module "Relativity"
      Mc2 = rest_energy(g_me)   ! incident positron rest energy; module "Relativity"
      E_right = W_max(Mc2, mtc2, .false., Prtcl%Ekin, 1.0d-8, hw_phonon) ! module "CS_integration_limits"
      
      Pos_elastic = numpar%Pos_elastic
      if (Pos_elastic == 1) then
         max_E0 = maxval(matter%CDF_phonon%E0(:))
         call find_Wmax_equal_Wmin(g_me, matter%Mean_Mass, .false., Prtcl%Ekin, 1.0d-6, max_E0, Eeq, hw_phonon)   ! module "CS_integration_limits"
         ! Check if delta-functional CDF works here,  and apply for electrons above chosen threshold:
         if (Prtcl%Ekin < numpar%CDF_Eeq_elast*Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
             Pos_elastic = 4
         endif
      endif ! (Pos_elastic == 1)
      
      ! Sample the cross section:
      call random_number(RN)
      Element => matter%Elements(1) ! unused here, so just to pass something into the subroutine
      call get_pos_elastic_CS(Prtcl%Ekin, matter, Element, numpar, CS_tot)   ! module "CS_electrons_elastic"
      ! Sample the cross section:
      CS_sampled = RN*CS_tot
   
      if ( (RN < 1.0d-10) .or. (CS_sampled < 1.0d-12) ) then   ! no need to sample, it's just the lower limit
         E_cur = E_left
      else ! sample it
         if (Pos_elastic /= 4) then    ! for analytical CSs, use bisection search algorithm:
            ! Start finding CS:
            E_cur = (E_left + E_right)*0.5d0
            call get_pos_elastic_CS(Prtcl%Ekin, matter, Element, numpar, CS_cur, E_max=E_cur)   ! module "CS_electrons_elastic"
            ! Search by bisection method:
!             do while (ABS(CS_cur - CS_sampled) > eps*CS_sampled)
            do while (abs(E_left - E_right) >= eps*E_left)
               if (CS_cur > CS_sampled) then
                  E_right = E_cur
               else
                  E_left = E_cur
               endif
               E_cur = (E_left + E_right)/2.0d0
               if (abs(E_left - E_right) < eps) exit  ! precise enough
               call get_pos_elastic_CS(Prtcl%Ekin, matter, Element, numpar, CS_cur, E_max=E_cur)   ! module "CS_electrons_elastic"
            enddo
!             print*, 'elastic_energy_transfer_p 1', Pos_elastic, CS_sampled, CS_cur, E_cur
         else   ! for numerically integrable CS, use direct search
             ! Target mean atomic number:
             Zeff = 1.0d0 + Equilibrium_charge_SHI(Prtcl%Ekin, g_me, matter%Mean_Z, (matter%Mean_Z-1.0d0), 0, 1.0d0) ! module "SHI_charge_state"
             call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Prtcl%Ekin, g_me, Zeff, 1.0d-10, matter%T_eV, matter%CDF_phonon, matter%Mean_Mass, &
                        matter%DOS%k, matter%DOS%Eff_m, .false., 1.0d0, matter%At_Dens, matter%DOS, numpar%CDF_model, &
                        hw_phonon = hw_phonon, Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
!              print*, 'elastic_energy_transfer_p 2', Pos_elastic, CS_sampled, CS_cur, E_cur
         endif  ! (El_elast /= 4)
      endif ! ( (RN < 1.0d-10) .or. (CS_sampled < 1.0d-12) )
   
      ! Output: sampled transferred energy:
      dE = E_cur
   
      ! Get the corresponding scattering angle:
      mu = cos_theta_from_W(Prtcl%Ekin, dE, g_me, matter%Mean_Mass) ! module "MC_general_tools"
!       print*, 'mu', mu, Prtcl%Ekin, dE, E_left, CS_sampled, CS_tot
      theta = acos(mu)

   case (2)  ! Mott
   
      ! 1) Select the element the electron scatters on:
      ! Get the total cross section:
      call interpolate_data_single(matter%Pos_elastic_total%E(:), matter%Pos_elastic_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
      ! Get the random number, to sample partial cross section:
      call random_number(RN)   ! intrinsic FORTRAN subroutine
      ! Realized CS:
      CS_sampled = RN*CS_tot
      ! Find which of the partial cross sections correspond to the sampled one:
      CS_sum = 0.0d0   ! to start with
      N_elem = dble(SUM(matter%Elements(:)%percentage))   ! number of elements in this compound material
      ! Scan through all the atoms to find the one:
      found_at = .false.   ! to start with
      KOA = matter%N_Elements  ! by default, it is the last element
      ! (no need to check, because true by exclusion if all the other elements are not true):
      SRCH:do j =1, matter%N_Elements-1     ! for each element, expect for the last one
         Element => matter%Elements(j)          ! all information about this element
         elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         ! For each shell, find its partial cross section:
         call interpolate_data_single(Element%Pos_elastic%E,  Element%Pos_elastic%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
         CS_sum = CS_sum + CS_cur*elem_contrib
         ! Check if this shell is the sampled one according to the integral CS:
         call check_element(CS_sum, CS_sampled, j, KOA1, found_at)    ! module "MC_general_tools"
!          if (found_at) exit SRCH  ! no need to continue the search if we already found the KOA
         if (found_at) then
            KOA = KOA1
            exit SRCH  ! no need to continue the search if we already found the KOA
         endif
      enddo SRCH
   
      if ( (KOA > size(matter%Elements(:))) .or. (KOA < 1) ) then
         print*, 'Pos Elast', CS_sampled, RN*CS_tot, CS_tot
      endif
   
      ! Now we know which element the electron scatters on:
      Element => matter%Elements(KOA)  ! all information about this element
   
      ! 2) Sample polar angle according to the differential cross section
      call  get_positron_elastic_polar_angle(Prtcl%Ekin, Element, theta) ! below
         
      ! 3) Change electron energy accordingly:
      dE = transfered_E_from_theta(Prtcl%Ekin, theta, g_me, Element%M) ! module "MC_general_tools"

   case (3) ! DSF
      ! NOT READY YET
   case default ! no scattering
      dE = 0.0d0
      theta = 0.0d0
   end select

   nullify(Element)
end subroutine elastic_energy_transfer_p


subroutine get_positron_elastic_polar_angle(Ekin, Element, theta)
   real(8), intent(in) :: Ekin  ! [eV] positron energy
   type(Atom_kind), intent(in) :: Element   ! parameters of the element the photon scatters on
   real(8), intent(out) :: theta   ! polar angle of photoelectron emission
   real(8) :: RN
   ! Sample the angle:
   call random_number(RN)
   ! The scattering angle is:
   theta = Mott_sample_mu(Ekin, Element%Zat, RN)  ! module "CS_electrons_elastic"
end subroutine get_positron_elastic_polar_angle


!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

subroutine find_type_of_positron_event(used_target, numpar, Prtcl, i_type)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Positron), intent(in) :: Prtcl  ! the Positron to perform some event
   integer, intent(out) :: i_type   ! type of event: -1=box crossing; 0=boundary; 1=inelastic; 2=elastic; 3=Bremsstrahlung; 4=annihilation
   !----------------------------------------------
   real(8) :: eps, MFP, V
   logical :: out
   type(Target_atoms), pointer :: matter
   
   eps = 1.0d-12

!    print*, 'Pos', Prtcl%Ekin, Prtcl%ti, Prtcl%t_sc, Prtcl%R(:)
   
   ! 1) Check if it is a boundary crossing:
   if (abs(Prtcl%ti - Prtcl%t_sc) > eps) then   ! it is a boundary crossing
      ! 2) Check whether it is simulation box crossing or a target crossing:
      out = out_of_simulation_box(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), numpar%box_start_x, numpar%box_end_x, &
               numpar%box_start_y, numpar%box_end_y, numpar%box_start_z, numpar%box_end_z) ! module "MC_general_tools"
      if (out) then ! it is crossing the simulation box boundary
         i_type = -1
      else  ! it is crossing the target boundary
         i_type = 0
      endif
   else ! it is a scattering-type event
      ! Properties of the target material, inside of which the particle is:
      matter => used_target%Material(Prtcl%in_target)
     
!        print*, 'Pos2', Prtcl%Ekin, Prtcl%ti, Prtcl%t_sc, Prtcl%R(:)
     
      ! 3) Find out what kind of scattering event it is:
      call find_type_of_scattering(i_type, Prtcl%Ekin,  matter%Pos_inelastic_total%E, matter%Pos_inelastic_total%Total, &
                    E_array2=matter%Pos_elastic_total%E, CS_array2=matter%Pos_elastic_total%Total, &
                    E_array3=matter%Pos_Brems_total%E, CS_array3=matter%Pos_Brems_total%Total, &
                    E_array4=matter%Pos_annihil_total%E, CS_array4=matter%Pos_annihil_total%Total) ! module "CS_general_tools"
   endif
   nullify(matter)
end subroutine find_type_of_positron_event


   
   
end module MC_positron
