! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019-2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains Monte Carlo routines to simulate SHI events
module MC_SHI
use Universal_constants
use Objects
use Geometries !, only: m_tollerance_eps
use Little_subroutines, only: interpolate_data_single
use Relativity, only: velosity_from_kinetic_energy
use MC_general_tools, only: flight_time_to_boundary, out_of_simulation_box, sample_phi, sample_theta, &
                                            cos_recoil_from_W, cos_theta_from_W, check_phi, particle_cross_box_boundary, deflect_velosity, &
                                            get_electron_flight_time, get_hole_flight_time, get_photon_flight_time, get_shi_flight_time, check_shell, &
                                            remove_particle_from_MC_array, extend_MC_array, find_the_target
use CS_general_tools, only: total_CS_from_chennels, MFP_from_sigma, Time_from_MFP, find_type_of_scattering, find_valence_hole_mass
use CS_ion_inelastic, only: get_SHI_inelastic_energy_transfer
use Dealing_with_DOS, only: select_energy_dos

implicit none

 contains

subroutine MC_SHI_event(used_target, numpar, bunch, N_SHI, Prtcl, NOP, MC, MD_supce, E_e, E_h)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   integer, intent(inout) :: N_SHI   ! number of SHIs
   type(SHI), dimension(:), allocatable, intent(inout) :: Prtcl  ! the SHI to perform some event
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout) :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h ! data to pass to MD later
   !--------------------------------
   integer :: i_type   ! type of event: -1=box crossing; 0=boundary; 1=inelastic;
   
   ! 0) Advance the last-step parameters:
   Prtcl(NOP)%R0(:) = Prtcl(NOP)%R(:)
   Prtcl(NOP)%V0(:) = Prtcl(NOP)%V(:)
   
   ! 1) Move the particle to the end-point:
   Prtcl(NOP)%R(:) = Prtcl(NOP)%R0(:) + (Prtcl(NOP)%V(:)) * (Prtcl(NOP)%ti - Prtcl(NOP)%t0) !* g_ms2Afs

   !print*, 'SHI', sqrt(SUM(Prtcl(NOP)%R(:))), Prtcl(NOP)%ti - Prtcl(NOP)%t0, sqrt(SUM(Prtcl(NOP)%V(:)))
!    write(*,'(a,f,f,f,f,f)') "SHI", Prtcl(NOP)%R(:), Prtcl(NOP)%ti, Prtcl(NOP)%ti - Prtcl(NOP)%t0
   
   Prtcl(NOP)%t0 = Prtcl(NOP)%ti

   ! 2) Choose an SHI event:
   call find_type_of_SHI_event(used_target, numpar, Prtcl(NOP), i_type)  ! below
   
   ! 3) Perform the event according to the chosen type:
   select case(i_type)
   case (0) ! target boundary crossing
!       print*, 'event_SHI_target_boundary'
      call event_SHI_target_boundary(used_target, numpar, Prtcl(NOP))   ! below
      
   case (1) ! inelastic
      call event_SHI_inelastic(used_target, numpar, bunch, MC, NOP, MD_supce, E_e, E_h)   ! below
      
   case (2) ! elastic
      ! Not introduced yet...
      
   case default ! simulation box boundary crossing
!     print*, 'MC_SHI_event 4', Prtcl(NOP)%R(3), Prtcl(NOP)%R0(3), Prtcl(NOP)%V(3) * (Prtcl(NOP)%ti - Prtcl(NOP)%t0), &
!                    Prtcl(NOP)%R0(3) + (Prtcl(NOP)%V(3)) * (Prtcl(NOP)%ti - Prtcl(NOP)%t0)          
      call particle_cross_box_boundary(used_target, numpar, N_SHI, NOP, Prtcl)    ! module "MC_general_tools"
!       print*, 'MC_SHI_event 6' 
   endselect
   
end subroutine MC_SHI_event


!ССССССССССССССССССССССССССССССССССССССССССС
! SHI crossing target boundary:
subroutine event_SHI_target_boundary(used_target, numpar, Prtcl)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(SHI), intent(inout) :: Prtcl        ! SHI as an object
   !------------------------------------------------------
   real(8), dimension(3)  :: R_shift
   real(8) :: PrtclV 
   integer :: target_ind
   
   ! Find into which target this SHI enters:
   ! 1) Find which target's boundary the particle is crossing:
   !R_shift = 1.0d-7 * Prtcl%V(:)     ! to place particle inside of the material
   PrtclV = max( SQRT( SUM( Prtcl%V(:)*Prtcl%V(:) ) ), m_tollerance_eps) ! exclude zero
   R_shift = m_tollerance_eps * Prtcl%V(:)/PrtclV     ! to place particle inside of the material
   ! Update particle's material index according to the new material it enters:
   call find_the_target(used_target, Prtcl, R_shift) ! module "MC_general_tools"
   
!    print*, 'SHI target_ind', target_ind
   
   ! Test case (assume vacuum):
!    Prtcl%in_target = 0   ! vacuum
   call get_SHI_flight_time(used_target, numpar, Prtcl)  ! module "MC_general_tools"

end subroutine event_SHI_target_boundary



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SHI inelastic scattering subroutines:
subroutine event_SHI_inelastic(used_target, numpar, bunch, MC, NOP, MD_supce, E_e, E_h)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Radiation_param), dimension(:), intent(in) :: bunch	! incomming radiation
   integer, intent(in) :: NOP   ! index of particle in the array
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e, E_h ! data to pass to MD later
   !------------------------------------------------------
   integer :: KOA, NSH, TSHI, ZSHI_amu, i
   real(8) :: Ekin, Ekin_new, dE, dE_SAVE
   real(8) :: theta, phi, theta_r, phi_r, eps, MSHI_amu, MSHI
   real(8) :: V0(3), V(3), V_tot, V_e_abs, V_e(3), E_DOS, a0(3), a(3)
   real(8) :: polar_angle_h, azimuth_angle_h, Mh, mass_temp, V_h_abs, V_h(3)
   real(8) :: CS_total
   type(SHI), pointer :: Prtcl  ! the electron to perform some event
   type(Target_atoms), pointer :: matter
!    type(Atom_kind), pointer :: Element
   logical :: valent
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_SHIs(NOP)
   matter => used_target%Material(Prtcl%in_target)
   TSHI = Prtcl%KOA ! type of SHI (corresponding to bunch number, where info about the ion is storred)
   ! Set the model parameters:
   MSHI_amu = bunch(TSHI)%Meff  ! Atomic mass of the ion
   MSHI = MSHI_amu*g_amu    ! [kg]
   ZSHI_amu = bunch(TSHI)%Z ! Atomic number of SHI
   
!    print*, 'shi 0', ZSHI_amu, MSHI_amu
   
   ! 1) Which shell of which element is being ionized: (tested, correct)
   call select_shell_ionization_SHI(used_target, numpar, MC, NOP, TSHI, KOA, NSH)  ! below
   
!    if (NSH > 1) then    ! testing
!       print*, 'shi 1', KOA, NSH
!    endif
!    print*, 'shi 1', KOA, NSH
   
   ! 2) Sample transfered energy in the impact ionization event:
   ! Sample the transfered energy: (tested, correct)
   if (KOA == 0) then   ! valence band:
      valent = .true.
!       print*, 'Valence SHI ionization', Prtcl%Ekin
      ! Find the precalculated total cross section:
      call interpolate_data_single(matter%SHI_inelastic_valent(Prtcl%KOA)%E, matter%SHI_inelastic_valent(Prtcl%KOA)%Total, &
                        Prtcl%Ekin, CS_total) ! module "Little_subroutines"
      ! Pass it into energy transfer finding:
      dE = get_SHI_inelastic_energy_transfer(Prtcl%Ekin, ZSHI_amu, MSHI_amu, matter, matter%DOS, numpar, bunch, TSHI, KOA, NSH, &
                        matter%DOS%Egap, CS_total)    ! module "CS_ion_inelastic"
   else ! core shell:
      valent = .false.
      ! Find the precalculated total cross section:
      call interpolate_data_single(matter%Elements(KOA)%SHI_inelastic(Prtcl%KOA)%E, matter%Elements(KOA)%SHI_inelastic(Prtcl%KOA)%Per_shell(NSH,:), &
                        Prtcl%Ekin, CS_total) ! module "Little_subroutines"
      ! Pass it into energy transfer finding:
      dE = get_SHI_inelastic_energy_transfer(Prtcl%Ekin, ZSHI_amu, MSHI_amu, matter, matter%DOS, numpar, bunch, TSHI, KOA, NSH, &
                        matter%Elements(KOA)%Ip(NSH), CS_total)    ! module "CS_ion_inelastic"
   endif
   
!     print*, 'shi 2', valent, dE
!    print*, 'event_SHI_inelastic TEST_dE:', Prtcl%Ekin, dE
!    print*, '------------------------------'

!     dE = 1000.0 ! testing
!     do i = 1, 1000 ! testing
    !dE = dble(i) * 180.0d0 ! testing

   ! 2a) Scattering angles according to the transferred energy:
   call sample_angles_inelastic_SHI(Prtcl%Ekin, dE, MSHI, theta, phi) ! below
   
   ! And the angles of emission of a new electron: (tested, correct)
   call SHI_set_angles_of_new_electron(Prtcl%Ekin, dE, MSHI, phi, theta_r, phi_r)  ! below
   
   ! 3) Energy an electron is emitted with:
   if (KOA == 0) then   ! valence band:
      ! Sample from where within the valence band it is ionized according to DOS:
      call select_energy_DOS(matter%DOS%E, matter%DOS%DOS, matter%Integral_DOS_fe, &
                                           matter%DOS%Egap, dE, matter%DOS%alpha_CB, E_DOS)   ! module "Dealing_with_DOS"
      Ekin = dE + (E_DOS - matter%DOS%Egap) ! E_DOS is counted from the top of VB here!
   else ! core shell:
      E_DOS = 0.0d0
      Ekin = dE - matter%Elements(KOA)%Ip(NSH)
   endif
   
!    if (dE > 150.0d0) then
!       print*, 's1', Ekin, E_DOS
!       print*, 's2', MC%N_e + 1
!       print*, 's3', dE
!    endif
   
   ! Make sure the transfered energy is sufficient:
   if (Ekin < 0.0d0) then   ! insufficient energy is replaced with the minimal allowed:
      print*, "Error in event_SHI_inelastic", Ekin, E_DOS
      Ekin = 0.0d0
      if (KOA == 0) then   ! valence band:
         dE = -E_DOS + matter%DOS%Egap
      else
         dE = matter%Elements(KOA)%Ip(NSH)
      endif
   endif
   
   ! 4) Change the incident ion direction of motion:
   ! Get the cosines of the ion velosity:
   V_tot = sqrt(SUM(Prtcl%V(:)*Prtcl%V(:)))
   V0(:) = Prtcl%V(:)/V_tot
   ! Get the deflection angle:
   call deflect_velosity(V0(1), V0(2), V0(3), theta, phi, V(1), V(2), V(3))  ! module "MC_general_tools"
   ! Get the NEW absolute SHI velosity, to set its components according to the cosines:
   Prtcl%Ekin = Prtcl%Ekin - dE ! update ion energy
   V_e_abs = velosity_from_kinetic_energy(Prtcl%Ekin, MSHI)     ! [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the ion velosity [A/fs]
   
!    write(*,'(a,f,f,f,f,f,f,f,f)') Prtcl%V(:), V_e(:), theta, phi
   
   ! Update new velosity:
   Prtcl%V(:) = V_e(:)
   
   ! Get the SHI next sampled free path:
   call get_SHI_flight_time(used_target, numpar, Prtcl)  ! module "MC_general_tools"
   
   ! 4a) Get the direction of the new electron:
   call deflect_velosity(V0(1), V0(2), V0(3), theta_r, phi_r, V(1), V(2), V(3))  ! module "MC_general_tools"
   V_e_abs = velosity_from_kinetic_energy(Ekin, g_me)     ! [A/fs] module "Relativity"
   V_e(:) = V_e_abs * V(:)  ! components of the electron velosity [A/fs]

!    write(*,'(a,f,f,f,e,e,e)') 'shi:', Ekin, theta_r, phi, V_e(:)

   ! 5) Make a new electron:
   MC%N_e = MC%N_e + 1
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_e > size(MC%MC_Electrons)) call extend_MC_array(MC%MC_Electrons)    ! module "MC_general_tools"
   call make_new_particle(MC%MC_Electrons(MC%N_e), Ekin=Ekin, t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_e)    ! module "Objects"
   
   ! 6) Get new electrons time of the next event:
   call get_electron_flight_time(used_target, numpar, MC%MC_Electrons(MC%N_e), MD_supce, E_e)  ! module "MC_general_tools"

   if (MC%MC_Electrons(MC%N_e)%Ekin < 0.0d0) then
      print*, 'Error in (event_SHI_inelastic):'
      print*, 'Electron got negative energy:', MC%MC_Electrons(MC%N_e)%Ekin
      print*, 'N=', MC%N_e, MC%N_e
      print*, MC%MC_Electrons(MC%N_e)%R(:)
      print*, MC%MC_Electrons(MC%N_e)%R0(:)
      print*, MC%MC_Electrons(MC%N_e)%V(:)
      print*, MC%MC_Electrons(MC%N_e)%V0(:)
      print*, MC%MC_Electrons(MC%N_e)%ti - MC%MC_Electrons(MC%N_e)%t0, MC%MC_Electrons(MC%N_e)%ti, MC%MC_Electrons(MC%N_e)%t0
   endif


!    write(*,'(f,f,f,f)') dE, MC%MC_Electrons(MC%N_e)%ti, MC%MC_Electrons(MC%N_e)%ti * V_e_abs
!    enddo ! testing
!    pause 'event_SHI_inelastic'

   
   ! 7) Create a new hole after electron removal:
   MC%N_h = MC%N_h + 1
   ! for valent hole, sample its velosity:
   if (valent) then
      ! Find hole's mass:
      Mh =  find_valence_hole_mass(numpar, used_target%Material(Prtcl%in_target)%DOS, abs(E_DOS)) ! module "MC_general_tools"

      ! Set the velosity of the new hole:
      V_h_abs = velosity_from_kinetic_energy( abs(E_DOS), Mh)    ! [A/fs] module "Relativity"
      polar_angle_h = sample_phi() ! module "MC_general_tools"
      azimuth_angle_h = sample_theta()  ! module "MC_general_tools"
      V0(3) = 1.0d0
      V0(1:2) = 0.0d0
      ! Get the holes velosity in the sampled direction:   
!       call  deflect_velosity(V0(1), V0(2), V0(3), polar_angle_h, azimuth_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      call  deflect_velosity(V0(1), V0(2), V0(3), azimuth_angle_h, polar_angle_h, V_h(1), V_h(2), V_h(3))  ! module "MC_general_tools"
      V_h(:) = V_h_abs * V_h(:)   ! set also the absolute value
   else
      E_DOS = 0.0d0
      V_h(:) = 0.0d0   ! core holes don't fly
   endif
   
   if (.not.valent .and. (KOA == 0)) then  ! inconsistency found
      print*, 'event_SHI_inelastic', valent, KOA, NSH
   endif
   
   ! Save parameters of the new particles into array:
   ! in case we have more particles than spaces in the array, extend the array:
   if (MC%N_h > size(MC%MC_Holes)) call extend_MC_array(MC%MC_Holes)    ! module "MC_general_tools"
   call make_new_particle(MC%MC_Holes(MC%N_h), Ekin=abs(E_DOS), t0=Prtcl%t0, &
                                      generation=Prtcl%generation+1, in_target=Prtcl%in_target, R=Prtcl%R, V=V_h, &
                                      KOA = KOA, Sh = NSH, valent = valent )    ! module "Objects"

   ! 8) Get new holes time of the next event:
   call get_hole_flight_time(used_target, numpar, MC%MC_Holes(MC%N_h), MD_supce, E_h)  ! module "MC_general_tools"
   
!    print*, 'shi 5', MC%N_h
      
   nullify(Prtcl)
end subroutine event_SHI_inelastic



subroutine select_shell_ionization_SHI(used_target, numpar, MC, NOP, ipart, KOA, NSH)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: NOP   ! index of particle in the array
   integer, intent(in) :: ipart     ! type of SHI
   type(MC_arrays), intent(inout), target :: MC      ! elements of MC array for all particles in one iteration
   integer, intent(out) :: KOA, NSH ! kind of atom, and the shell number, which absorbs the photon
   !------------------------------------------------------
   type(Target_atoms), pointer :: matter
   type(SHI), pointer :: Prtcl  ! the SHI to perform some event
   type(Atom_kind), pointer :: Element
   real(8) :: RN, CS_tot, CS_sampled, CS_cur, CS_sum,  elem_contrib, N_elem
   integer :: i_mat, j, k
   logical :: found_shl, valence_done
   
   ! Starting default values (valence band):
   KOA = 0
   NSH = 0
   valence_done = .false.   ! to start, we did not calculate valence band yet
   
   ! Pointers for easier access to scattering particle and target properties:
   Prtcl => MC%MC_SHIs(NOP)
   matter => used_target%Material(Prtcl%in_target)
   
   ! Get the total cross section:
   call interpolate_data_single(matter%SHI_inelastic_total(ipart)%E(:), matter%SHI_inelastic_total(ipart)%Total(:), Prtcl%Ekin, CS_tot) ! module "Little_subroutines"
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
         
         VAL:if ( (Element%valent(k)) .and. (allocated(matter%CDF_valence%A)) ) then    ! Valence band (not for RBEB atomic model!)
            if (.not.valence_done) then ! add VB (but only once):
               call interpolate_data_single(matter%SHI_inelastic_valent(ipart)%E, matter%SHI_inelastic_valent(ipart)%Total(:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
               CS_sum = CS_sum + CS_cur    ! band contribution
               valence_done = .true.     ! we added VB once (not do double-count it)
            endif
         else VAL
            call interpolate_data_single(Element%SHI_inelastic(ipart)%E, Element%SHI_inelastic(ipart)%Per_shell(k,:), Prtcl%Ekin, CS_cur) ! module "Little_subroutines"
            CS_sum = CS_sum + CS_cur*elem_contrib   ! element contribution
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
end subroutine select_shell_ionization_SHI



subroutine sample_angles_inelastic_SHI(Ekin, dE, h_m, theta, phi)
   real(8), intent(in) :: Ekin, dE  ! [eV] incident and transfered energy
   real(8), intent(in) :: h_m   ! [kg] SHI mass
   real(8), intent(out) :: theta, phi   ! scattering angles
   real(8) :: mu
   ! Assume uniform scattering probability into 2Pi for phi:
   phi = sample_phi()    ! module "MC_general_tools"
   ! Assume abolutely elastic scattering for theta:
   mu = cos_theta_from_W(Ekin, dE, h_m, g_me)    ! module "MC_general_tools"
   theta = acos(mu)
end subroutine sample_angles_inelastic_SHI


pure subroutine SHI_set_angles_of_new_electron(Ekin, dE, h_m, phi, theta_r, phi_r)  
   real(8), intent(in) :: Ekin, dE  ! [eV] incident and transfered energy
   real(8), intent(in) :: h_m   ! [kg] SHI mass
   real(8), intent(in) :: phi   ! angle of deflection of the original particle
   real(8), intent(out) :: theta_r, phi_r   ! emission angles of a new electron
   real(8) :: mu
   ! Assume uniform scattering probability into 2Pi for phi:
   phi_r = phi + g_Pi
   call check_phi(phi_r) ! module "MC_general_tools"
   ! Assume abolutely elastic scattering for theta:
   mu = cos_recoil_from_W(Ekin, dE, h_m, g_me)    ! module "MC_general_tools"
   theta_r = acos(mu)
end subroutine SHI_set_angles_of_new_electron



!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

subroutine find_type_of_SHI_event(used_target, numpar, Prtcl, i_type)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(SHI), intent(in) :: Prtcl  ! the SHI to perform some event
   integer, intent(out) :: i_type   ! type of event: -1=box crossing; 0=boundary; 1=inelastic;
   !----------------------------------------------
   real(8) :: eps, MFP, V
   logical :: out
   type(Target_atoms), pointer :: matter
   
   eps = 1.0d-12

!    print*, 'SHI', Prtcl%Ekin, Prtcl%ti, Prtcl%t_sc, Prtcl%R(:), Prtcl%V(:)
   
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
      ! 3) Find out what kind of scattering event it is (currently only inelastic is included...):
      call find_type_of_scattering(i_type, Prtcl%Ekin,  matter%SHI_inelastic_total(Prtcl%KOA)%E, matter%SHI_inelastic_total(Prtcl%KOA)%Total) ! module "CS_general_tools"
   endif

   nullify(matter)
end subroutine find_type_of_SHI_event

   
   
end module MC_SHI
