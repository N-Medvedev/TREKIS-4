! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019-2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains Monte Carlo routines to simulate hole events
module MC_hole
use Universal_constants
use Objects
use Geometries
use Little_subroutines, only: interpolate_data_single, print_time_step
use Relativity, only: velosity_from_kinetic_energy, rest_energy
use MC_general_tools, only: flight_time_to_boundary, out_of_simulation_box, select_process_by_time, sample_phi, sample_theta, &
                            cos_theta_from_W, check_phi, cos_recoil_from_W, transfered_E_from_theta, &
                            particle_cross_box_boundary, deflect_velosity, get_electron_flight_time, get_hole_flight_time, &
                            get_photon_flight_time, check_element, remove_particle_from_MC_array, extend_MC_array, &
                            find_the_target, define_normal_to_surface, reflection_from_surface, add_energy_into_MD_array
use CS_integration_limits, only: W_max, find_Wmax_equal_Wmin
use CS_general_tools, only: total_CS_from_chennels, MFP_from_sigma, Time_from_MFP, find_type_of_scattering, find_valence_hole_mass
use CS_electrons_inelastic, only: get_inelastic_energy_transfer, CDF_total_CS_nonrel
use CS_electrons_elastic, only: get_el_elastic_CS, Mott_sample_mu
use CS_holes_elastic, only: get_h_elastic_CS
use CS_holes_inelastic, only: get_inelastic_energy_transfer_h
use Dealing_with_DOS, only: select_energy_DOS
use SHI_charge_state, only: Equilibrium_charge_SHI

implicit none

 contains

subroutine MC_hole_event(used_target, numpar, N_h, Prtcl, NOP, MC, MD_supce, E_h_at, E_e, E_h)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(inout) :: N_h   ! number of holes
   type(Hole), dimension(:), allocatable, intent(inout) :: Prtcl  ! the hole to perform some event
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout) :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_h_at, E_e, E_h ! data to pass to MD later
   !----------------------------------
   integer :: i_type   ! type of event: VALENCE: -1=box crossing; 0=boundary; 1=inelastic; 2=elastic; CORE: -10=Auger; -11=radiative decay
   integer :: INFO  ! to check if procidure executed correctly
   real(8) :: t0, R0(3)    ! save for checking purposes

   R0 = Prtcl(NOP)%R0   ! save for testing
   ! 0) Advance the last-step parameters:
   Prtcl(NOP)%R0(:) = Prtcl(NOP)%R(:)
   Prtcl(NOP)%V0(:) = Prtcl(NOP)%V(:)
   
!    print*, 'MC_hole_event', Prtcl(NOP)%t0, Prtcl(NOP)%ti
   
   ! 1) Move the particle to the end-point:
   if (Prtcl(NOP)%valent) then  ! it can move
      Prtcl(NOP)%R(:) = Prtcl(NOP)%R0(:) + (Prtcl(NOP)%V(:)) * (Prtcl(NOP)%ti - Prtcl(NOP)%t0) !* g_ms2Afs
   else     ! core holes do not move
      Prtcl(NOP)%R(:) = Prtcl(NOP)%R0(:)
   endif
   t0 = Prtcl(NOP)%t0   ! save for testing
!    print*, 'MC_hole_event', Prtcl(NOP)%t0, Prtcl(NOP)%ti, sqrt(Prtcl(NOP)%R(1)*Prtcl(NOP)%R(1) + Prtcl(NOP)%R(2)*Prtcl(NOP)%R(2))
   Prtcl(NOP)%t0 = Prtcl(NOP)%ti
   
   ! 2) Choose a hole event:
   call find_type_of_hole_event(used_target, numpar, Prtcl(NOP), NOP, i_type, N_h)  ! below
   
!    print*, 'MC_hole_event', NOP, i_type
   
!    if ((Prtcl(NOP)%R(3) < 0.0d0 - m_tollerance_eps) .or. (Prtcl(NOP)%R(3) > 10.0d0+m_tollerance_eps)) then
!       print*,  'HOLE OUT OF THE BOX!'
!       print*, NOP, Prtcl(NOP)%R(3), Prtcl(NOP)%V(3)
!     endif
   
   ! 3) Perform the event according to the chosen type:
   select case(i_type)
   case (0) ! target boundary crossing
!        call print_time_step('Hole Before boundary:', dble(NOP), msec=.true.)   ! module "Little_subroutines"
!       call particle_cross_box_boundary(used_target, numpar, N_h, NOP, Prtcl, type_of_periodicity=2) ! module "MC_general_tools"
!       print*, 'event_hole_target_boundary'
      call event_hole_target_boundary(used_target, numpar, Prtcl(NOP), NOP, INFO, MD_supce, E_h)   ! below
      if (INFO /= 0) then   ! somthing went wrong in the boundary scattering:
         print*, 't0=', t0
         print*, 'R0=', R0
         print*, '****************************************'
      endif
      
!        call print_time_step('Hole After boundary:', dble(NOP), msec=.true.)   ! module "Little_subroutines"
!       pause 'MC_hole_event : event_hole_target_boundary'
   ! Here are the events for valence holes:
   case (1) ! inelastic
      call event_hole_inelastic(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)   ! below
   case (2) ! elastic
      call event_hole_elastic(used_target, numpar, MC, NOP, MD_supce, E_h, E_h_at)   ! below
   
   ! And here are the events for core holes:
   case (-10) ! core-hole Auger decay
      call Auger_decay(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)    ! below
   case (-11) ! core-hole radiative decay
      call radiative_decay(used_target, numpar, MC, NOP, MD_supce, E_h)    ! below
   
   case default ! simulation box boundary crossing
!         print*, 'MC_hole_event 4', Prtcl(NOP)%R(3), Prtcl(NOP)%R0(3), Prtcl(NOP)%V(3) * (Prtcl(NOP)%ti - Prtcl(NOP)%t0), &
!                     Prtcl(NOP)%R0(3) + (Prtcl(NOP)%V(3)) * (Prtcl(NOP)%ti - Prtcl(NOP)%t0)    
!         print*, 'MC_hole_event 4', NOP, Prtcl(NOP)%R(3), Prtcl(NOP)%V(3)
!         if ((Prtcl(NOP)%R(3) < 0.0d0) .or. (Prtcl(NOP)%R(3) > 10.0d0)) print*, 'OUT OF THE BOX!'
      call particle_cross_box_boundary(used_target, numpar, N_h, NOP, Prtcl)    ! module "MC_general_tools"
!         print*, 'MC_hole_event 6', NOP, Prtcl(NOP)%R(3), Prtcl(NOP)%V(3)
!         if ((Prtcl(NOP)%R(3) < 0.0d0) .or. (Prtcl(NOP)%R(3) > 10.0d0)) then
!             print*, 'Error in MC_hole_event: HOLE OUT OF THE BOX'
!             print*, NOP, Prtcl(NOP)%R(:), Prtcl(NOP)%V(:)
!         endif
   endselect   
end subroutine MC_hole_event


!ССССССССССССССССССССССССССССССССССССССССССС
! Hole crossing target boundary:
subroutine event_hole_target_boundary(used_target, numpar, Prtcl, NOP, INFO, MD_supce, E_h)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Hole), intent(inout) :: Prtcl        ! hole as an object
   integer, intent(in) :: NOP   ! number of hole
   integer, intent(inout) :: INFO  ! info about errors to pass to main subroutine
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_h ! data to pass to MD later
   !------------------------------------------------------
   real(8), dimension(3)  :: norm_to_surf
   real(8) :: Vabs

   ! Currently, holes can only reflect from the surface, not cross it
   ! Find normal to the surface to reflecto from:
   call define_normal_to_surface(used_target,  Prtcl, norm_to_surf, 'event_hole_target_boundary', INFO)    ! module "MC_general_tools"
   if (INFO /=0) then   !some error occured
      print*, 'Hole: define_normal_to_surface failed for Prtcl:', NOP
      print*, 'R=', Prtcl%R
      print*, 'V=', Prtcl%V
      print*, 'E=', Prtcl%Ekin
      print*, 'ti=', Prtcl%ti
      print*, 't_sc=', Prtcl%t_sc
      print*, 'Inside target #', Prtcl%in_target
      print*, 'n to surface:', norm_to_surf
   endif
   
   ! Change velosity according to reflection from the surface with given normal:
   call reflection_from_surface(Prtcl%V, norm_to_surf)   ! module "MC_general_tools"
   
   ! Shift hole just across the border:
   Vabs = SQRT( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
   if (Vabs >= m_tollerance_eps) then
      !Prtcl%R(:) = Prtcl%R(:) + m_tollerance_eps * Prtcl%V(:)/Vabs
      Prtcl%R(:) = Prtcl%R(:) + m_tollerance_eps * Prtcl%V(:)/abs(Prtcl%V(:))
   else 
      print*, 'Hole has zero velosity'
      Prtcl%R(:) = Prtcl%R(:) + m_tollerance_eps
   endif
   
   ! Find next scattering event on the boundary (but keep old (in-)elastic scattering time):
   call get_hole_flight_time(used_target, numpar, Prtcl, MD_supce, E_h, no_scatternig=.true.)  ! module "MC_general_tools"
end subroutine event_hole_target_boundary


!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Auger decay subroutines:
subroutine Auger_decay(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h   ! data to pass to MD later
   !------------------------------------------------------
   type(Hole), pointer :: Prtcl  ! the hole to perform some event
   type(Target_atoms), pointer :: matter
   type(Atom_kind), pointer :: Element
   integer :: KOA, NSH, i_1, i_2, sh_selected_1, sh_selected_2, Nsiz
   integer :: sh_selected_1_in, sh_selected_2_in, KOA_in
   real(8) :: Ee, phi, theta, V0(3), V_e(3), Vtot, Mh, V_h_abs, polar_angle_h, azimuth_angle_h, V_h(3)
   real(8) :: E_DOS_1, E_DOS_2, E_fin_1, E_fin_2, RN, f_tot, f_sampled, f_sum
   
   ! 0) Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Holes(NOP)    ! all particles parameters
   matter => used_target%Material(Prtcl%in_target)  ! all materials parameters
   KOA = Prtcl%KOA    ! Kind of atom the hole belongs to
   NSH = Prtcl%Sh       ! Shell the hole is in
   Element => matter%Elements(KOA)  ! all atomic data for this element
   Nsiz = size(Element%Ip) ! Number of shells
   
   ! 1) Find the partner shells for Auger decay:
   ! To sample the shells:
   call random_number(RN)   ! intrinsic FORTRAN subroutine
   f_tot = SUM(Element%f_auger(NSH,:,:))    ! to normalize the total probability
   ! Sampled probability to find the shell hole pops up into:
   f_sampled = f_tot * RN
   ! Find the sampled shell:
   f_sum = 0.0d0
   SHLS:do i_1 = 1, Nsiz
      do i_2 = 1, Nsiz
         f_sum = f_sum + Element%f_auger(NSH,i_1,i_2)   ! probabilities
         if (f_sum >= f_sampled) then
            sh_selected_1 = i_1    ! save the first shell number
            sh_selected_2 = i_2    ! save the second shell number
            exit SHLS   ! we found the shell, exit the cycle
         endif
      enddo
   enddo SHLS

   ! In case EADL didn't have data for shell-resolved probabilities,
   ! assume the nearest shells participate in Auger, as the most probable ones:
   if (f_sum < 1.0e-10) then
      sh_selected_1 = Nsiz ! save the first shell number (use as default)
      sh_selected_2 = Nsiz ! save the second shell number (use as default)
      if (NSH < Nsiz-1) then ! no need to search, it's VB:
         SHL2:do i_1 = 1, Nsiz
            if (Element%valent(i_1)) then ! VB
               sh_selected_1 = Nsiz ! save the first shell number
               sh_selected_2 = Nsiz ! save the second shell number
               exit SHL2   ! found shells, no need to continue
            else
               sh_selected_1 = i_1    ! save the first shell number
               if (matter%Elements(KOA)%Ip(sh_selected_1) <= matter%Elements(KOA)%Ip(NSH)) then
                  do i_2 = i_1, Nsiz
                     sh_selected_1 = i_2    ! save the second shell number
                     if (matter%Elements(KOA)%Ip(sh_selected_2) <= &
                        (matter%Elements(KOA)%Ip(NSH) - matter%Elements(KOA)%Ip(sh_selected_1)) ) then
                        exit SHL2   ! found shells, no need to continue
                     endif
                  enddo
               else
                  sh_selected_1 = Nsiz ! use default (valence band)
               endif
            endif
         enddo SHL2
      endif
   endif

   if ((sh_selected_1 <= NSH) .or. (sh_selected_2 <= NSH)) then
      print*, 'Error #1 in (Auger_decay):'
      print*, 'Auger-electrons are from shells deeper than the original hole:'
      print*, NSH, sh_selected_1, sh_selected_2
      print*, f_tot, RN
      print*, f_sum, f_sampled
      print*, '-------------'
      !print*, Element%f_auger(NSH,i_1,i_2)
   endif

   ! In case the first one is a valence hole, sample from where within the valence band it is ionized according to DOS:
   if (Element%valent(sh_selected_1)) then
      !call select_energy_DOS(matter%DOS%E, matter%DOS%DOS, matter%Integral_DOS_fe, &
      !          matter%DOS%Egap, 1.0d6, matter%DOS%alpha_CB, E_DOS_1)   ! module "Dealing_with_DOS"
      call select_energy_DOS(matter%DOS%E, matter%DOS%DOS, matter%Integral_DOS_fe, &
                matter%DOS%Egap, (Element%Ip(NSH)-matter%DOS%Egap), matter%DOS%alpha_CB, E_DOS_1)   ! module "Dealing_with_DOS"
      E_fin_1 = -(E_DOS_1 - matter%DOS%Egap)    ! [eV] energy level where the hole ends up
   else
      E_DOS_1 = 0.0d0 ! no width of a core level
      E_fin_1 = Element%Ip(sh_selected_1)  ! [eV] energy level where the hole ends up
   endif
   ! In case the second one is a valence hole, sample from where within the valence band it is ionized according to DOS:
   if (Element%valent(sh_selected_2)) then
      call select_energy_DOS(matter%DOS%E, matter%DOS%DOS, matter%Integral_DOS_fe, &
                matter%DOS%Egap, (Element%Ip(NSH)-E_fin_1), matter%DOS%alpha_CB, E_DOS_2)   ! module "Dealing_with_DOS"
                !matter%DOS%Egap, 1.0d6, matter%DOS%alpha_CB, E_DOS_2)   ! module "Dealing_with_DOS"
      E_fin_2 = -(E_DOS_2 - matter%DOS%Egap)    ! [eV] energy level where the hole ends up
   else
      E_DOS_2 = 0.0d0 ! no width of a core level
      E_fin_2 = Element%Ip(sh_selected_2)  ! [eV] energy level where the hole ends up
   endif
   
   ! 2) Make a new electron:
   MC%N_e = MC%N_e + 1
   ! 2a) Set the energy of the emitted electron:
   Ee = Element%Ip(NSH) - E_fin_1 - E_fin_2  ! [eV] electron energy
   ! 2b) Set electron angles:
   phi = sample_phi() ! module "MC_general_tools"
   theta = sample_theta()   ! module "MC_general_tools"
   ! 2c) Set electron velosity:
   V0(1) = 1.0d0
   V0(2:3) = 0.0d0
   ! Get the holes velosity in the sampled direction:
   call  deflect_velosity(V0(1), V0(2), V0(3), theta, phi, V_e(1), V_e(2), V_e(3))  ! module "MC_general_tools"
   Vtot = velosity_from_kinetic_energy(Ee, g_me)    ! [A/fs] module "Relativity"
   V_e(:) = Vtot * V_e(:)   ! derection of electron motion
   ! 2d) Save the electron data into the particle properties:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_e > size(MC%MC_Electrons)) call extend_MC_array(MC%MC_Electrons)    ! module "Objects"
   call make_new_particle(MC%MC_Electrons(MC%N_e), Ekin=Ee, t0=Prtcl%t0, &
             generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_e)    ! module "Objects"
   ! 2e) Get new electrons time of the next event:
   call get_electron_flight_time(used_target, numpar, MC%MC_Electrons(MC%N_e), MD_supce, E_e)  ! module "MC_general_tools"
   
   if ((MC%MC_Electrons(MC%N_e)%ti < Prtcl%ti) .or. (MC%MC_Electrons(MC%N_e)%Ekin < 0.0d0)) then
      print*, 'Error #2 in (Auger_decay):'
      print*, 'Electron got negative energy:', MC%MC_Electrons(MC%N_e)%Ekin
      print*, 'N=', MC%N_e, MC%N_h
      print*, MC%MC_Electrons(MC%N_e)%R(:)
      print*, MC%MC_Electrons(MC%N_e)%R0(:)
      print*, MC%MC_Electrons(MC%N_e)%V(:)
      print*, MC%MC_Electrons(MC%N_e)%V0(:)
      print*, MC%MC_Electrons(MC%N_e)%ti - MC%MC_Electrons(MC%N_e)%t0, MC%MC_Electrons(MC%N_e)%ti, MC%MC_Electrons(MC%N_e)%t0
      print*, 'Ip:', Prtcl%KOA, Prtcl%Sh
      print*, 'sh:', sh_selected_1, sh_selected_2
      print*, Element%valent(sh_selected_1), Element%valent(sh_selected_2)
      print*, Element%Ip(NSH), E_fin_1, E_fin_2
      print*, E_DOS_1, E_DOS_2, matter%DOS%Egap
      print*, '-------------'
      !pause 'ERROR #2 PAUSE'
   endif


   ! 3) Make a new hole:
   MC%N_h = MC%N_h + 1
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_h > size(MC%MC_Holes)) then
      nullify(Prtcl, matter)    ! nullify pointers
      call extend_MC_array(MC%MC_Holes)   ! module "MC_general_tools"
      Prtcl => MC%MC_Holes(NOP)    ! all particles parameters
      matter => used_target%Material(Prtcl%in_target)  ! all materials parameters
   endif
   ! 3a) for valent hole, sample its velosity:
   if (Element%valent(sh_selected_2)) then  ! valence hole
      sh_selected_2_in = 0  ! marker of the valence band
      KOA_in = 0
      MC%MC_Holes(MC%N_h)%valent = .true.
      MC%MC_Holes(MC%N_h)%KOA = 0
      MC%MC_Holes(MC%N_h)%Sh = 0
      ! Find hole's mass:
      Mh =  find_valence_hole_mass(numpar, used_target%Material(Prtcl%in_target)%DOS, abs(E_DOS_2)) ! module "MC_general_tools"
      ! Set the velosity of the new hole:
      V_h_abs = velosity_from_kinetic_energy(abs(E_DOS_2), Mh)    ! [A/fs] module "Relativity"
      polar_angle_h = sample_phi() ! module "MC_general_tools"
      azimuth_angle_h = sample_theta()  ! module "MC_general_tools"
      V0(3) = 1.0d0
      V0(1:2) = 0.0d0
      ! Get the holes velosity in the sampled direction:   
!       call  deflect_velosity(V0(1), V0(2), V0(3), polar_angle_h, azimuth_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      call  deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle_h, polar_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      V_h(:) = V_h_abs * V_h(:)   ! set also the absolute value
   else ! core hole
      sh_selected_2_in = sh_selected_2  ! marker of the deep shell
      KOA_in = KOA
      E_DOS_2 = 0.0d0
      V_h(:) = 0.0d0   ! core holes don't fly
   endif
   ! Save parameters of the new particles into array:
   call make_new_particle(MC%MC_Holes(MC%N_h), Ekin=abs(E_DOS_2), t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_h, &
                                      !KOA = KOA, Sh = sh_selected_2, valent = Element%valent(sh_selected_2) )    ! module "Objects"
                                      KOA = KOA_in, Sh = sh_selected_2_in, valent = Element%valent(sh_selected_2) )    ! module "Objects"
   ! 3b) Get new holes time of the next event:
   call get_hole_flight_time(used_target, numpar, MC%MC_Holes(MC%N_h), MD_supce, E_h)  ! module "MC_general_tools"
   
   
   ! 4) Update parameters of the old hole:
   ! 4a) Get the old hole next sampled event:
   Prtcl%Sh = sh_selected_1      ! shell the hole jumped into
   if (Element%valent(Prtcl%Sh)) then   ! valence hole
      Prtcl%valent = .true.
      Prtcl%KOA = 0
      Prtcl%Sh = 0
      ! Find hole's mass:
      Mh =  find_valence_hole_mass(numpar, used_target%Material(Prtcl%in_target)%DOS, abs(E_DOS_1)) ! module "MC_general_tools"
      ! Set the velosity of the new hole:
      V_h_abs = velosity_from_kinetic_energy(abs(E_DOS_1), Mh)    ! [A/fs] module "Relativity"
      polar_angle_h = sample_phi() ! module "MC_general_tools"
      azimuth_angle_h = sample_theta()  ! module "MC_general_tools"
      V0(3) = 1.0d0
      V0(1:2) = 0.0d0
      ! Get the holes velosity in the sampled direction:
!       call  deflect_velosity(V0(1), V0(2), V0(3), polar_angle_h, azimuth_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      call  deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle_h, polar_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      V_h(:) = V_h_abs * V_h(:)   ! set also the absolute value
      Prtcl%V = V_h ! to set new velosity
      Prtcl%V0 = Prtcl%V ! to save
      Prtcl%Ekin = abs(E_DOS_1) ! [eV]
   else ! core hole
      V_h(:) = 0.0d0   ! core holes don't fly
   endif
   ! 4b) New time of the next event:
   call get_hole_flight_time(used_target, numpar, Prtcl, MD_supce, E_h)  ! module "MC_general_tools"
   
   nullify(Prtcl, matter, Element)
end subroutine Auger_decay


!RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
! Radiative decay subroutines:
subroutine radiative_decay(used_target, numpar, MC, NOP, MD_supce, E_h)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_h   ! data to pass to MD later
   !------------------------------------------------------
   type(Hole), pointer :: Prtcl  ! the hole to perform some event
   type(Target_atoms), pointer :: matter
   type(Atom_kind), pointer :: Element
   integer :: KOA, NSH, i, Nsiz, sh_selected
   real(8) :: f_tot, RN, f_sampled, f_sum, E_fin
   real(8) :: Eph, phi, theta, V0(3), V_ph(3), Mh, E_DOS, V_h_abs, polar_angle_h, azimuth_angle_h, V_h(3)
   
   ! 0) Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Holes(NOP)    ! all particles parameters
   matter => used_target%Material(Prtcl%in_target)  ! all materials parameters
   KOA = Prtcl%KOA    ! Kind of atom the hole belongs to
   NSH = Prtcl%Sh       ! Shell the hole is in
   Element => matter%Elements(KOA)  ! all atomic data for this element
   
!    print*, 'radiative_decay0', Prtcl%t0, Prtcl%ti
   
   
   ! 1) Find the shell that participates in the radiative decay:
   ! 1a) Get the random number, to sample partial cross section:
   call random_number(RN)   ! intrinsic FORTRAN subroutine
   
   ! 1b) Total probability to be normalized to 1:
   f_tot = SUM(Element%f_rad(NSH,:))
   ! Sampled probability to find the shell hole pops up into:
   f_sampled = f_tot * RN
   ! 1c) Find the corresponding shell:
   Nsiz = size(Element%Ip) ! Number of shells
   f_sum = 0.0d0    ! to start counting
   do i = NSH, Nsiz
      f_sum = f_sum + Element%f_rad(NSH,i)
      if (f_sum >= f_sampled) then
         sh_selected = i    ! save the shell number
         exit   ! we found the shell, exit the cycle
      endif
   enddo

   ! In case it is a valence hole, sample from where within the valence band it is ionized according to DOS:
   if (Element%valent(sh_selected)) then
      call select_energy_DOS(matter%DOS%E, matter%DOS%DOS, matter%Integral_DOS_fe, &
                                           matter%DOS%Egap, 1.0d6, matter%DOS%alpha_CB, E_DOS)   ! module "Dealing_with_DOS"
      E_fin = -(E_DOS - matter%DOS%Egap)    ! [eV] energy level where the hole ends up
   else
      E_DOS = 0.0d0 ! no width of a core level
      E_fin = Element%Ip(sh_selected)  ! [eV] energy level where the hole ends up
   endif
   
   ! 2) Create emitted photon:
   MC%N_ph = MC%N_ph + 1
   ! 2a) Define emitted photon parameters:
   Eph = Element%Ip(NSH) - E_fin  ! [eV] photon energy
   ! 2b) Set photon angles:
   phi = sample_phi() ! module "MC_general_tools"
   theta = sample_theta()   ! module "MC_general_tools"
   ! 2c) Set photon velosity:
   V0(1) = 1.0d0
   V0(2:3) = 0.0d0
   ! Get the holes velosity in the sampled direction:
   call  deflect_velosity(V0(1), V0(2), V0(3), theta, phi, V_ph(1), V_ph(2), V_ph(3))  ! module "MC_general_tools"
!    V_ph(:) =g_cvel * V_ph(:)   ! derection of photon motion [m/s]
   V_ph(:) =g_c_Afs * V_ph(:)   ! derection of photon motion [m/s]
   ! 2d) Save the photon data into the particle properties:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_ph > size(MC%MC_Photons)) call extend_MC_array(MC%MC_Photons)     ! module "MC_general_tools"
   call make_new_particle(MC%MC_Photons(MC%N_ph), Ekin=Eph, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_ph)    ! module "Objects"
   ! 2e) Define the next photon scattering event:
   call get_photon_flight_time(used_target, numpar, MC%MC_Photons(MC%N_ph))  ! module "MC_general_tools"
   
!    print*, 'radiative_decay1', MC%MC_Photons(MC%N_ph)%t0, MC%MC_Photons(MC%N_ph)%ti
   
   ! 3) Update parameters of the hole:
   Prtcl%Sh = sh_selected      ! shell the hole jumped into
   ! 3a) Get the old hole next sampled event:
   if (Element%valent(Prtcl%Sh)) then   ! valence hole
      Prtcl%valent = .true.
      Prtcl%KOA = 0
      Prtcl%Sh = 0
      ! Find hole's mass:
      Mh =  find_valence_hole_mass(numpar, used_target%Material(Prtcl%in_target)%DOS, abs(E_DOS)) ! module "MC_general_tools"
      ! Set the velosity of the new hole:
      V_h_abs = velosity_from_kinetic_energy(abs(E_DOS), Mh)    ! [A/fs] module "Relativity"
      polar_angle_h = sample_phi() ! module "MC_general_tools"
      azimuth_angle_h = sample_theta()  ! module "MC_general_tools"
      V0(3) = 1.0d0
      V0(1:2) = 0.0d0
      ! Get the holes velosity in the sampled direction:   
!       call  deflect_velosity(V0(1), V0(2), V0(3), polar_angle_h, azimuth_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      call deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle_h, polar_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      V_h(:) = V_h_abs * V_h(:)   ! set also the absolute value
      Prtcl%V = V_h ! to set new velosity
      Prtcl%V0 = Prtcl%V ! to save
      Prtcl%Ekin = abs(E_DOS) ! to save
   endif
   
   ! New time of the next event:
   call get_hole_flight_time(used_target, numpar, Prtcl, MD_supce, E_h)  ! module "MC_general_tools"
   
!    print*, 'radiative_decay2', Prtcl%t0, Prtcl%ti

   nullify(Prtcl, matter, Element)
end subroutine radiative_decay



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Hole inelastic scattering subroutines:
subroutine event_hole_inelastic(used_target, numpar, MC, NOP, MD_supce, E_e, E_h)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h   ! data to pass to MD later
   !------------------------------------------------------
   real(8) :: Ekin, Ekin_new, dE, theta, phi, theta_r, phi_r, eps
   real(8) :: V0(3), V(3), V_tot, V_e_abs, V_e(3), E_DOS, a0(3), a(3)
   real(8) :: polar_angle_h, azimuth_angle_h, Mh, mass_temp, V_h_abs, V_h(3)
   type(Hole), pointer :: Prtcl  ! the electron to perform some event
   type(Target_atoms), pointer :: matter
!    type(Atom_kind), pointer :: Element
   logical :: valent, possible
   real(8) :: h_m
   
   possible = .true.    ! by default, consider scattering possible
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Holes(NOP)
   matter => used_target%Material(Prtcl%in_target)
   ! Get valence hole mass (negative -Ekin is because it's counted from the top of VB):
   h_m =  find_valence_hole_mass(numpar, matter%DOS, Prtcl%Ekin)    ! module "MC_general_tools"
   ! Check just in case of some problems with discret DOS:
   if (h_m <= 0.0d0) then   ! use free electron mass
      h_m = g_me
   endif
   
   ! 1) A hole can only ionize valence band:
   valent = .true.
   ! Sample transfered energy in the impact ionization event:
   dE = get_inelastic_energy_transfer_h(Prtcl%Ekin, matter, matter%DOS, numpar, 0, 0, matter%DOS%Egap)    ! module "CS_holes_inelastic"

!     print*, 'event_hole_inelastic', Prtcl%Ekin, dE
   
   ! 2a) Scattering angles according to the transferred energy:
   call sample_angles_inelastic_hole(Prtcl%Ekin, dE, h_m, theta, phi, possible) ! below
   ! If such an energy transfer is not possible, skip all the subroutine:
   if (.not.possible) goto 99999
   
   ! And the angles of emission of a new electron:
   call hole_set_angles_of_new_electron(Prtcl%Ekin, dE, h_m, phi, theta_r, phi_r)  ! below
   
   ! 3) Energy an electron is emitted with:
   ! Sample from where within the valence band it is ionized according to DOS:
   call select_energy_DOS(matter%DOS%E, matter%DOS%DOS, matter%Integral_DOS_fe, &
                                        matter%DOS%Egap, dE, matter%DOS%alpha_CB, E_DOS)   ! module "Dealing_with_DOS"
   Ekin = dE + (E_DOS - matter%DOS%Egap)

   ! Make sure the transfered energy is sufficient:
   if (Ekin < 0.0d0) then   ! insufficient energy is replaced with the minimal allowed:
      Ekin = 0.0d0
      dE = -E_DOS + matter%DOS%Egap
   endif
   
   ! 4) Change the incident hole direction of motion:
   ! Get the cosines of the hole velosity:
   V_tot = sqrt(SUM(Prtcl%V(:)*Prtcl%V(:)))
   V0(:) = Prtcl%V(:)/V_tot
   ! Get the deflection angle:
   call deflect_velosity(V0(1), V0(2), V0(3), theta, phi, V(1), V(2), V(3))  ! module "MC_general_tools"
   ! Get the NEW absolute hole velosity, to set its components according to the cosines:
   Prtcl%Ekin = Prtcl%Ekin - dE ! update old hole energy
   V_e_abs = velosity_from_kinetic_energy(Prtcl%Ekin, h_m)     ! [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the hole velosity
   ! Update velosity:
   Prtcl%V(:) = V_e(:)
   
!    print*, 'event_hole_inelastic-2:', Prtcl%V(:)
  
   ! 4a) Get the direction of the new electron:
   call deflect_velosity(V0(1), V0(2), V0(3), theta_r, phi_r, V(1), V(2), V(3))  ! module "MC_general_tools"
   V_e_abs = velosity_from_kinetic_energy(Ekin, g_me)     ! [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the electron velosity

   ! 5) Make a new electron:
   MC%N_e = MC%N_e + 1
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_e > size(MC%MC_Electrons)) call extend_MC_array(MC%MC_Electrons)     ! module "MC_general_tools"
   call make_new_particle(MC%MC_Electrons(MC%N_e), Ekin=Ekin, t0=Prtcl%t0, &
             generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_e)    ! module "Objects"
   
   ! 6) Get new electrons time of the next event:
   call get_electron_flight_time(used_target, numpar, MC%MC_Electrons(MC%N_e), MD_supce, E_e)  ! module "MC_general_tools"
   
   if ((MC%MC_Electrons(MC%N_e)%ti < Prtcl%ti) .or. (MC%MC_Electrons(MC%N_e)%Ekin < 0.0d0)) then
      print*, 'Error in (event_hole_inelastic):'
      print*, 'Electron got negative energy:', MC%MC_Electrons(MC%N_e)%Ekin
      print*, 'N=', MC%N_e, MC%N_e
      print*, MC%MC_Electrons(MC%N_e)%R(:)
      print*, MC%MC_Electrons(MC%N_e)%R0(:)
      print*, MC%MC_Electrons(MC%N_e)%V(:)
      print*, MC%MC_Electrons(MC%N_e)%V0(:)
      print*, MC%MC_Electrons(MC%N_e)%ti - MC%MC_Electrons(MC%N_e)%t0, MC%MC_Electrons(MC%N_e)%ti, MC%MC_Electrons(MC%N_e)%t0
   endif


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
      ! Get the holes velosity in the sampled direction:
      call  deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle_h, polar_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      V_h(:) = V_h_abs * V_h(:)   ! set also the absolute value
   else
      E_DOS = 0.0d0
      V_h(:) = 0.0d0   ! core holes don't fly
   endif
   ! Save parameters of the new particles into array:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_h > size(MC%MC_Holes)) then
      nullify(Prtcl, matter)    ! nullify pointers
      call extend_MC_array(MC%MC_Holes)    ! module "MC_general_tools"
      ! Pointers for easier access to scattering particle and target properties:
      Prtcl => MC%MC_Holes(NOP)
      matter => used_target%Material(Prtcl%in_target)
   endif
   call make_new_particle(MC%MC_Holes(MC%N_h), Ekin=abs(E_DOS), t0=Prtcl%t0, &
                          generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_h, &
                          KOA = 0, Sh = 0, valent = valent )    ! module "Objects"

   ! 8) Get new hole's time of the next event:
   call get_hole_flight_time(used_target, numpar, MC%MC_Holes(MC%N_h), MD_supce, E_h)  ! module "MC_general_tools"
   
   ! 9) Get the old hole's time of the next event:
99999  call get_hole_flight_time(used_target, numpar, Prtcl, MD_supce, E_h)  ! module "MC_general_tools"

   nullify(Prtcl, matter)
end subroutine event_hole_inelastic



subroutine sample_angles_inelastic_hole(Ekin, dE, h_m, theta, phi, possible)
   real(8), intent(in) :: Ekin, dE  ! [eV] incident and transfered energy
   real(8), intent(in) :: h_m   ! [kg] holes mass
   real(8), intent(out) :: theta, phi   ! scattering angles
   logical, intent(inout) :: possible   ! is such energy transfer is even possible?
   real(8) :: mu
   ! Assume uniform scattering probability into 2Pi for phi:
   phi = sample_phi()    ! module "MC_general_tools"
   ! Assume absolutely elastic scattering for theta (no recoil):
   mu = cos_theta_from_W(Ekin, dE, h_m, g_me)    ! module "MC_general_tools"
   if (abs(mu) <= 1.0d0) then   ! it is possible
      theta = acos(mu)
      possible = .true.
   else ! such energy transfer is impossible
      theta = 0.0d0
      possible = .false.
!       print*, 'sample_angles_inelastic_hole', mu, Ekin, dE
   endif
end subroutine sample_angles_inelastic_hole


subroutine hole_set_angles_of_new_electron(Ekin, dE, h_m, phi, theta_r, phi_r)  
   real(8), intent(in) :: Ekin, dE  ! [eV] incident and transfered energy
   real(8), intent(in) :: h_m   ! [kg] hole mass
   real(8), intent(in) :: phi   ! angle of deflection of the original particle
   real(8), intent(out) :: theta_r, phi_r   ! emission angles of a new electron
   real(8) :: mu
   ! Assume uniform scattering probability into 2Pi for phi:
   phi_r = phi + g_Pi
   call check_phi(phi_r) ! module "MC_general_tools"
   ! Assume absolutely elastic scattering for theta:
   mu = cos_recoil_from_W(Ekin, dE, h_m, g_me)    ! module "MC_general_tools"
   if (abs(mu) <= 1.0d0) then   ! no problem here, calculate the angle
      theta_r = acos(mu)
   else ! some problem here, sample random angle
      theta_r = sample_theta() ! module "MC_general_tools"
   endif
end subroutine hole_set_angles_of_new_electron



!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! Hole elastic scattering subroutines:

subroutine event_hole_elastic(used_target, numpar, MC, NOP, MD_supce, E_h, E_h_at)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_h, E_h_at   ! data to pass to MD later
   !------------------------------------------------------
   type(Hole), pointer :: Prtcl  ! the electron to perform some event
   type(Target_atoms), pointer :: matter
   type(Atom_kind), pointer :: Element
   real(8) :: RN, phi, theta, CS_tot, CS_sampled, CS_sum, N_elem, elem_contrib, CS_cur, Vtot
   real(8) :: a0(3), a(3), dE, h_m, Ekin
   integer :: j, KOA, KOA1
   logical :: found_at
   
   ! 0) Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_Holes(NOP)
   matter => used_target%Material(Prtcl%in_target)
   ! Get valence hole mass:
   h_m =  find_valence_hole_mass(numpar, matter%DOS, Prtcl%Ekin)    ! module "MC_general_tools"

   ! 1) Set the deflection angles:
   ! Sample azimuthal angle randomly within 2Pi:
   phi = sample_phi()   ! module "MC_general_tools"
   
   ! 2,3) Find out transfered energy and deflecting angle:
   call elastic_energy_transfer_hole(numpar, Prtcl, matter, h_m, KOA, dE, theta)   ! below
   ! Just in case, check that the transfered energy is not too high:
   if (dE > Prtcl%Ekin) dE = Prtcl%Ekin ! [eV]
   
   ! 4) Get the cosines of the electron velosity:
   Vtot = sqrt( SUM( Prtcl%V(:)*Prtcl%V(:) ) )
   a0(:) = Prtcl%V(:)/Vtot
   
   ! 5) Get the deflection angle (angles of scattering of electron in the laboratory system):
!    call deflect_velosity(a0(1), a0(2), a0(3), phi, theta, a(1), a(2), a(3))  ! module "MC_general_tools"
   call deflect_velosity(a0(1), a0(2), a0(3), theta, phi, a(1), a(2), a(3))  ! module "MC_general_tools"
   
   ! 5b) Change hole energy accordingly:  
   Ekin = Prtcl%Ekin    ! save for testing
   Prtcl%Ekin = Prtcl%Ekin - dE ! [eV]
   
   ! 5c) Change the velosity accordingly:
   Vtot = velosity_from_kinetic_energy(Prtcl%Ekin, h_m)  ! [A/fs] module "Relativity"
   Prtcl%V(:) = Vtot * a(:)  ! components of the hole velosity
   
   ! 6) Now update the hole parameters:
   call get_hole_flight_time(used_target, numpar, Prtcl, MD_supce, E_h)  ! module "MC_general_tools"
   
   ! 7) Save the atomic parameters in this collision:
   MC%N_at_nrg = MC%N_at_nrg + 1
   if (MC%N_at_nrg > size(MC%MC_Atoms_events)) call extend_MC_array(MC%MC_Atoms_events)    ! module "MC_general_tools"
   ! Save the parameters of this collision (ONLY ENERGY TRANSFER IS SAVED, MOMENTUM NOT DONE YET!)
   call make_new_particle(MC%MC_Atoms_events(MC%N_at_nrg), Ekin=dE, t0=Prtcl%t0, KOA = KOA, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R)    ! module "Objects"

   ! Save the energy to pass to MD module, if needed:
   if (numpar%DO_MD) then   ! if user requested MD at all
      call add_energy_into_MD_array(numpar, dE, Prtcl%R, MD_supce, E_h_at)   ! module "MC_general_tools"
   endif

   nullify(Prtcl, matter, Element)
end subroutine event_hole_elastic



subroutine elastic_energy_transfer_hole(numpar, Prtcl, matter, h_m, KOA, dE, theta)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Hole), intent(in) :: Prtcl  ! the hole to perform some event
   type(Target_atoms), intent(in), target :: matter ! material parameters
   real(8), intent(in) :: h_m  ! [kg] mass of hole
   integer, intent(out) :: KOA    ! kind of atom
   real(8), intent(out) :: dE, theta    ! [eV] transfered energy, and scattering angle
   !--------------------------------
   real(8) :: CS_tot, RN, CS_sampled, CS_sum, N_elem, elem_contrib, CS_cur, Eeq, max_E0, Zeff
   real(8) :: E_left, E_right, E_cur, eps, mu, hw_phonon, mtc2, Mc2, Se1, CS_tot_test
   logical :: found_at
   integer :: j, KOA1, H_elast 
   type(Atom_kind), pointer :: Element
   
   select case (numpar%H_elast)	! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   case (1,5)  ! CDF
      eps = 1.0d-3
      E_left = 0.0d0! [eV] minimal transferred energy
!       E_right = matter%E_debye    ! [eV] maximal transferred energy
      hw_phonon = maxval( matter%CDF_phonon%E0(:) + matter%CDF_phonon%Gamma(:) )
      mtc2 = rest_energy(matter%Mean_Mass)   ! module "Relativity"
      Mc2 = rest_energy(h_m)   ! incident hole rest energy; module "Relativity"
!       print*, 'elastic_energy_transfer_hole', Mc2, mtc2, Prtcl%Ekin
      E_right = W_max(Mc2, mtc2, .false., Prtcl%Ekin, 1.0d-8, hw_phonon) ! module "CS_integration_limits"
      
      H_elast = numpar%H_elast
      if ((H_elast  == 1) .or. (H_elast  == 5)) then
         max_E0 = maxval(matter%CDF_phonon%E0(:))
         call find_Wmax_equal_Wmin(h_m, matter%Mean_Mass, .false., Prtcl%Ekin, 1.0d-6, max_E0, Eeq, hw_phonon)   ! module "CS_integration_limits"
         ! Check if delta-functional CDF works here,  and apply for electrons above chosen threshold:
         if (Prtcl%Ekin < numpar%CDF_Eeq_elast*Eeq) then   ! switch to nonrelativistic numerically-integrable CDF:
             H_elast  = 4
         endif
      endif ! (H_elast  == 1)
      
      ! Sample the cross section:
      call random_number(RN)
      Element => matter%Elements(1) ! unused here, so just to pass something into the subroutine
!       call get_h_elastic_CS(Prtcl%Ekin, matter, Element, numpar, CS_tot_test)   ! module "CS_holes_elastic"
      ! Use precalculated total cross section:
      call interpolate_data_single(matter%H_elastic_total%E(:), matter%H_elastic_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
!       print*, 'H_elast:', CS_tot_test, CS_tot
      
      ! Sample the cross section:
      CS_sampled = RN*CS_tot
   
      if ( (RN < 1.0d-10) .or. (CS_sampled < 1.0d-12) )then   ! no need to sample, it's just the lower limit
         E_cur = E_left
      else ! sample it 
         if (H_elast /= 4) then    ! for analytical CSs, use bisection search algorithm:
            ! Start finding CS:
            E_cur = (E_left + E_right)*0.5d0
            call get_h_elastic_CS(Prtcl%Ekin, matter, Element, numpar, CS_cur, E_max=E_cur)   ! module "CS_holes_elastic"
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
               call get_h_elastic_CS(Prtcl%Ekin, matter, Element, numpar, CS_cur, E_max=E_cur)   ! module "CS_holes_elastic"
            enddo
!             print*, 'elastic_energy_transfer_hole 1', CS_sampled, CS_cur, E_cur
         else   ! for numerically integrable CS, use direct search
             ! Target mean atomic number:
             Zeff = 1.0d0 + Equilibrium_charge_SHI(Prtcl%Ekin, h_m, matter%Mean_Z, (matter%Mean_Z-1.0d0), 0, 1.0d0) ! module "SHI_charge_state"
             call CDF_total_CS_nonrel(numpar, CS_cur, Se1, Prtcl%Ekin, h_m, Zeff, 1.0d-10, matter%T_eV, matter%CDF_phonon, matter%Mean_Mass, &
                        matter%DOS%k, matter%DOS%Eff_m, .false., 1.0d0, matter%At_Dens, matter%DOS, numpar%CDF_model, &
                        hw_phonon = hw_phonon, Sigma_sampled=CS_sampled, E_sampled=E_cur)  ! module "CS_electrons_inelastic"
!              print*, 'elastic_energy_transfer_hole 2', CS_sampled, CS_cur, E_cur
         endif
      endif
   
      ! Output: sampled transferred energy:
      dE = E_cur
   
      ! Get the corresponding scattering angle:
      mu = cos_theta_from_W(Prtcl%Ekin, dE, h_m, matter%Mean_Mass) ! module "MC_general_tools"
!       print*, 'mu', mu, Prtcl%Ekin, dE, E_left, CS_sampled, CS_tot
      theta = acos(mu)

   case (2)  ! Mott
   
      ! 2) Select the element the hole scatters on:
      ! Get the total cross section:
      call interpolate_data_single(matter%H_elastic_total%E(:), matter%H_elastic_total%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
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
         call interpolate_data_single(Element%H_elastic%E,  Element%H_elastic%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
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
         print*, 'Elast', CS_sampled, RN*CS_tot, CS_tot
      endif
   
      ! Now we know which element the hole scatters on:
      Element => matter%Elements(KOA)  ! all information about this element

      ! 3) Sample polar angle according to the differential cross section
      call get_hole_elastic_polar_angle(Prtcl%Ekin, Element, h_m, theta)    ! below
         
      ! 3) Change hole energy accordingly:
      dE = transfered_E_from_theta(Prtcl%Ekin, theta, h_m, Element%M) ! module "MC_general_tools"
   case (3) ! DSF
      ! NOT READY YET
   case default ! no scattering
      dE = 0.0d0
      theta = 0.0d0
   end select

   nullify(Element)
end subroutine elastic_energy_transfer_hole



subroutine get_hole_elastic_polar_angle(Ekin, Element, mass, theta)
   real(8), intent(in) :: Ekin  ! [eV] photon energy
   type(Atom_kind), intent(in) :: Element   ! parameters of the element the photon scatters on
   real(8), intent(in) :: mass  ! [kg] hole mass
   real(8), intent(out) :: theta   ! polar angle of photoelectron emission
   real(8) :: RN
   ! Sample the angle:
   call random_number(RN)
   ! The scattering angle is:
   theta = Mott_sample_mu(Ekin, Element%Zat, RN, mass)  ! module "CS_electrons_elastic"
end subroutine get_hole_elastic_polar_angle


!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
subroutine find_type_of_hole_event(used_target, numpar, Prtcl, NOP, i_type, N_h)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Hole), intent(in) :: Prtcl  ! the Hole to perform some event
   integer, intent(in) :: NOP, N_h   ! index of particle in the array, total number of holes
   integer, intent(out) :: i_type   ! type of event: VALENCE: -1=box crossing; 0=boundary; 1=inelastic; 2=elastic; CORE: -10=Auger; -11=radiative decay
   !----------------------------------------------
   real(8) :: eps, MFP, V
   integer :: ind
   logical :: out
   type(Target_atoms), pointer :: matter
   
   !eps = 1.0d-12
   eps = 1.0d-10

!     print*, 'Ho', Prtcl%Ekin, Prtcl%ti, Prtcl%t_sc, 'R:', Prtcl%R(:), 'V:', Prtcl%V(:), Prtcl%valent
   
   VAL:if (Prtcl%valent) then   ! transport   
      ! 1) Check if it is a boundary crossing:
      if (abs(Prtcl%ti - Prtcl%t_sc) > eps) then   ! it is a boundary crossing
         ! 2) Check whether it is simulation box crossing or a target crossing:
         out = out_of_simulation_box(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), numpar%box_start_x, numpar%box_end_x, &
               numpar%box_start_y, numpar%box_end_y, numpar%box_start_z, numpar%box_end_z) ! module "MC_general_tools"
         if (out) then ! it is crossing the simulation box boundary
            i_type = -1
         else  ! it is crossing the target boundary
!             print*, 'TRGT cross', Prtcl%R
!             pause 'TRGT pause'
            i_type = 0
         endif
      else ! it is a scattering-type event
         ! Properties of the target material, inside of which the particle is:
         matter => used_target%Material(Prtcl%in_target)
         
!          print*, 'Ho2', Prtcl%Ekin, Prtcl%ti, Prtcl%t_sc, Prtcl%R(:), Prtcl%valent
         
         ! 3) Find out what kind of scattering event it is:
         call find_type_of_scattering(i_type, Prtcl%Ekin,  matter%H_inelastic_total%E, matter%H_inelastic_total%Total, &
                    E_array2=matter%H_elastic_total%E, CS_array2=matter%H_elastic_total%Total) ! module "CS_general_tools"
      endif
   else VAL ! core hole
      ! Sample the process: Auger vs radiative:
      matter => used_target%Material(Prtcl%in_target)
      if (Prtcl%KOA == 0) then  ! inconsistency found
         print*, 'ERROR in find_type_of_hole_event:', Prtcl%valent, Prtcl%KOA, Prtcl%Sh
         print*, 'R=', Prtcl%R
         print*, 'NOP=', NOP, N_h
         print*, 'target:', Prtcl%in_target
         print*, 'times:', Prtcl%ti, Prtcl%t_sc
      endif
      call select_process_by_time(ind, matter%Elements(Prtcl%KOA)%Auger(Prtcl%Sh), time2=matter%Elements(Prtcl%KOA)%Radiat(Prtcl%Sh))   ! module "MC_general_tools"
      select case (ind)
      case (2) ! radiative decay
         i_type = -11
      case default  ! assume Auger
         i_type = -10
      end select
   endif VAL
   
   nullify(matter)
end subroutine find_type_of_hole_event

  
   
end module MC_hole
