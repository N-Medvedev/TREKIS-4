! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! and F. Akhmetov
! in 2019-2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains generally used Monte Carlo routines
! References used in the module:
! [1]  F. Salvat, J. M. Fernandez-Varea, E. Acosta, J. Sempau 
!   "PENELOPE-2014 A Code System for Monte Carlo Simulation of Electron and Photon Transport", OECD (2014)
! [2] R. Rymzhanov et al. Phys. Status Solidi B 252, 159-164 (2015) / DOI 10.1002/pssb.201400130
! [3] M. Azzolini et al. J. Phys. Condens. Matter 31, 055901 (2019)


module MC_general_tools
use Universal_constants
use Objects
use Geometries, only: m_tollerance_eps, which_box_wall_crossing, where_crossing_box_boundary, shift_coordinate_system, &
                    rotate_coordinate_system, intersection_with_sphere, where_crossing_cylinder,  &
                    Cartesian_to_cylindrical, Cartesian_to_spherical, Spherical_to_cartesian
use CS_general_tools, only: total_CS_from_chennels, MFP_from_sigma, Time_from_MFP
use Little_subroutines, only: sample_Poisson, interpolate_data_single, Find_in_array_monoton
use Relativity, only: rest_energy, velosity_from_kinetic_energy

implicit none

 ! this is a function that set a spatial grid:
interface extend_MC_array
   module procedure extend_MC_array_Electrons
   module procedure extend_MC_array_Photons
   module procedure extend_MC_array_Holes
   module procedure extend_MC_array_Positrons
   module procedure extend_MC_array_Atoms
   module procedure extend_MC_array_SHIs
end interface extend_MC_array

interface copy_MC_array
   module procedure copy_MC_array_SHI
   module procedure copy_MC_array_atom
   module procedure copy_MC_array_hole
   module procedure copy_MC_array_positron
   module procedure copy_MC_array_electron
end interface copy_MC_array


public :: extend_MC_array, copy_MC_array


 contains


!MDMDMDMDMDMDMDMDMDMDMDMDMDMDMDMMDMDMD
! Subroutines related to coupling between MC and MD modules:
subroutine add_energy_into_MD_array(numpar, dE, R, MD_supce, E_array)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: dE    ! [eV] transferred energy
   real(8), dimension(3), intent(in) :: R   ! [A] coordinates of the event of the energy transfer
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_array  ! [eV] array that will be passed to MD
   !-----------------------------------
   logical :: in_the_box
   integer :: i_grid(3), i
   real(8) :: R_c, Theta_c, Phi_c, L_c


   ! Check if the event is within the MD supercell:
   in_the_box = check_if_inside_box(R, MD_supce%x_start, MD_supce%x_end, MD_supce%y_start, &
                                    MD_supce%y_end, MD_supce%z_start, MD_supce%z_end) ! below
   ! Only if the event is inside the MD supercell, transfer energy there:
   if (in_the_box) then ! energy can be transfered to MD
      i_grid = 1    ! to start with

      select case (MD_supce%coord_dim)
      case default ! zero dimension
         E_array(1,1,1) = E_array(1,1,1) + dE  ! [eV]

      case (1:3) ! 1,2 or 3d
         do i = 1, MD_supce%coord_dim
            ! Find position of event on the grid:
            select case (MD_supce%coord_type)
            case (1)   ! Cartesian
               select case (MD_supce%coord_axis(i))  ! along which axis
               case (1)   ! X
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i)%grid,  R(1), i_grid(i))	! module "Little_subroutines"
               case (2)   ! Y
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i)%grid,  R(2), i_grid(i))	! module "Little_subroutines"
               case (3)   ! Z
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i)%grid,  R(3), i_grid(i))	! module "Little_subroutines"
               endselect
            case (2)   ! Spherical
               ! Get the spherical coords:
               call Cartesian_to_spherical(R(1), R(2), R(3), R_c, Theta_c, Phi_c)  ! module "Geometries"
               ! Find where on the grid it is:
               select case (MD_supce%coord_axis(i))  ! along which axis
               case (1)   ! R
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i)%grid, R_c, i_grid(i))	! module "Little_subroutines"
               case (2)   ! Theta
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i)%grid, Theta_c, i_grid(i))	! module "Little_subroutines"
               case (3)   ! Phi
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i)%grid, Phi_c, i_grid(i))	! module "Little_subroutines"
               endselect
            case (3)   ! Cylindrical
               call Cartesian_to_cylindrical(R(1), R(2), R(3), R_c, L_c, Phi_c)  ! module "Geometries"
               select case (MD_supce%coord_axis(i))  ! along which axis
               case (1)   ! R
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i)%grid, R_c, i_grid(i))	! module "Little_subroutines"
               case (2)   ! Theta
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i)%grid, L_c, i_grid(i))	! module "Little_subroutines"
               case (3)   ! Phi
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i)%grid, Phi_c, i_grid(i))	! module "Little_subroutines"
               endselect
            endselect

         enddo ! i = 1, MD_supce%coord_dim

         ! Add to the right grid point:
         E_array(i_grid(1),i_grid(2),i_grid(3)) = E_array(i_grid(1),i_grid(2),i_grid(3)) + dE  ! [eV]

      endselect ! (MD_supce%coord_dim)
   else ! outside of the MD supercell
      ! Just for testing of periodic boundaries:
      if ( (numpar%periodic_x == 1) .and. (numpar%periodic_y == 1) .and. (numpar%periodic_z == 1) & ! periodic boundaries
                    .and. ALL(numpar%reset_from_MD(:)) ) then    ! and simulation box conforms with the MD supercell
         print*, 'Trouble in add_energy_into_MD_array, particle outside of the box:'
         print*, 'Coords:', R
         print*, 'Box size X:', MD_supce%x_start, MD_supce%x_end
         print*, 'Box size Y:', MD_supce%y_start, MD_supce%y_end
         print*, 'Box size Z:', MD_supce%z_start, MD_supce%z_end
      endif ! periodic boundaries
   endif ! (in_the_box)
end subroutine add_energy_into_MD_array



pure function check_if_inside_box(R, x_start, x_end, y_start, y_end, z_start, z_end) result(yesno)
   logical yesno    ! is the particle within the box or not?
   real(8), dimension(3), intent(in) :: R
   real(8), intent(in) :: x_start, x_end, y_start, y_end, z_start, z_end
   logical :: insideX, insideY, insideZ
   ! Check along X:
   if ((R(1) >= x_start) .and. (R(1) <= x_end)) then
      insideX = .true.
   else
      insideX = .false.
   endif
   ! Check along Y:
   if ((R(2) >= y_start) .and. (R(2) <= y_end)) then
      insideY = .true.
   else
      insideY = .false.
   endif
   ! Check along Z:
   if ((R(3) >= z_start) .and. (R(3) <= z_end)) then
      insideZ = .true.
   else
      insideZ = .false.
   endif

   yesno = (insideX .and. insideY .and. insideZ)
end function check_if_inside_box




!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Subrtines related to scattering angles:

subroutine deflect_velosity(u0, v0, w0, theta, phi, u, v, w)    ! Eq.(1.131), p.37 [1]
   real(8), intent(in) :: u0, v0, w0    ! cosine directions of the old velosity
   real(8), intent(in) :: theta, phi     ! polar (0,Pi) and azimuthal (0,2Pi) angles
   real(8), intent(out) :: u, v, w       ! new cosine directions
   real(8) :: sin_theta, cos_theta, sin_phi, cos_phi, one_w, eps, sin_t_w, temp
   !eps = 1.0d-8 ! margin of acceptance of w being along Z
   eps = m_tollerance_eps   ! module "Geometries"
   sin_theta = sin(theta)
   cos_theta = cos(theta)
   sin_phi = sin(phi)
   cos_phi = cos(phi)
   if ( abs(abs(w0)-1.0d0) < eps ) then   ! motion parallel to Z
      u = w0*sin_theta*cos_phi
      v = w0*sin_theta*sin_phi
      w = w0*cos_theta
   else ! any other direction of motion
      one_w = sqrt(1.0d0 - w0*w0)
      sin_t_w = sin_theta/one_w
      u = u0*cos_theta + sin_t_w*(u0*w0*cos_phi - v0*sin_phi)
      v = v0*cos_theta + sin_t_w*(v0*w0*cos_phi + u0*sin_phi)
      w = w0*cos_theta - one_w*sin_theta*cos_phi
!       if ( abs(abs(w0)-0.25) < eps*1.0d6 ) print*, 'deflect_velosity 1:', w0, w
   endif
   
   temp = sqrt(u*u + v*v + w*w)
   if (abs(temp-1.0d0) > eps) then  ! renormalize it:
      u = u/temp
      v = v/temp
      w = w/temp
   endif
end subroutine deflect_velosity

 
 
 

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
! Subrtines related to shells of elements:


subroutine check_shell(CS_sum, CS_sampled, i_atom, i_shell, KOA, NSH, found_shl)
   real(8), intent(in) :: CS_sum, CS_sampled    ! current integral CS, and the sampled one
   integer, intent(in) :: i_atom, i_shell  ! current number of atom and shell
   integer, intent(out) :: KOA, NSH     ! to save them
   logical, intent(out) :: found_shl    ! found the shell or not
   if (CS_sum >= CS_sampled) then    ! absorption by this shell of this element
      KOA = i_atom
      NSH = i_shell
      found_shl =.true.
   else
      found_shl =.false.
   endif
end subroutine check_shell


subroutine check_element(CS_sum, CS_sampled, i_atom, KOA, found_shl)
   real(8), intent(in) :: CS_sum, CS_sampled    ! current integral CS, and the sampled one
   integer, intent(in) :: i_atom  ! current number of atom and shell
   integer, intent(out) :: KOA     ! to save them
   logical, intent(out) :: found_shl    ! found the shell or not
   real(8) :: eps
   !eps = 1.0d0 - 1.0d-6     ! precision
   eps = 1.0d0 - m_tollerance_eps   ! module "Geometries"
   if (CS_sum >= CS_sampled * eps) then    ! absorption by this shell of this element
      KOA = i_atom
      found_shl =.true.
   else
      found_shl =.false.
   endif
end subroutine check_element


pure function total_electrons_in_target(Material) result(N_e_tot)
   real(8) N_e_tot  ! average number of electrons per molecule of the target
   type(Target_atoms), intent(in):: Material  !material parameters of a given target
   ! Get how many electrons we have per atom in each element:
   N_e_tot = SUM( dble(Material%Elements(:)%Zat) * dble(Material%Elements(:)%percentage) )
end function total_electrons_in_target



subroutine choose_random_shell(matter, KOA, NSH)
   type(Target_atoms), intent(in), target :: matter ! parameters of the target
   integer, intent(out) :: KOA, NSH ! number of atom, and number fo shell sampled randomly
   !--------------------------------------------------
   type(Atom_kind), pointer :: Element
   real(8) :: RN, N_el_tot, Ne_sampled, Ne_cur, elem_contrib
   integer :: i_atom, i_shell
   logical :: found_shl
   ! number of elements in this compound material:
   N_el_tot = total_electrons_in_target(matter) ! module "MC_general_tools"
   
   ! Get the random number, to sample partial cross section:
   call random_number(RN)   ! intrinsic FORTRAN subroutine
   ! Realized electron:
   Ne_sampled = RN*N_el_tot
   
   ! Find which of the partial cross sections correspond to the sampled one:
   Ne_cur = 0.0d0   ! to start with
   ! Scan through all the atoms and shells to find the one:
   ATM:do i_atom = 1, matter%N_Elements  ! for each element
      Element => matter%Elements(i_atom) ! all information about this element
      elem_contrib = dble(Element%percentage)   ! element contribution to this compound
      ! For each shell, its number of electrons:
      SHL:do i_shell = 1, Element%N_shl
         Ne_cur = Ne_cur + matter%Elements(i_atom)%Ne_shell(i_shell) * elem_contrib
         ! Check if this shell is the sampled one according to the sampled number of electron:
         call check_shell(Ne_cur, Ne_sampled, i_atom, i_shell, KOA, NSH, found_shl)    ! above
         if (found_shl) exit ATM    ! no need to continue the search if we already found the KOA
      enddo SHL
   enddo ATM
   
   nullify(Element)
end subroutine choose_random_shell

 
 
!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
! Subrtines related to boundary crossing:


pure function electron_transmission_probability(Work_func, Em_gamma, Em_E1, Ekin) result(T) 
   real(8) T    ! transmission probability / penetrability
   real(8), intent(in) :: Work_func, Em_gamma, Em_E1   ! precalculated barrier parameters
   real(8), intent(in) :: Ekin  ! electron kinetic energy, perpendicular to the surface component [eV]
   real(8) :: arg, hug
   
   if (Ekin < Work_func) then
      T = 0.0d0
   else
      arg = Em_gamma*(Em_E1 - Ekin)
      hug = log(huge(abs(arg)))
   
      if (arg < -hug) then
         T = 1.0d0
      elseif (arg > hug) then
         T = 0.0d0
      else
!          T = 1.0d0/(1.0d0+exp(Em_gamma*(Em_E1 - Ekin)))   ! Eq.(6) in [2]
         T = 1.0d0/(1.0d0+exp(arg))   ! Eq.(6) in [2]
      endif
   endif
end function electron_transmission_probability


pure function electron_transmission_probability_step(Work_func, Ekin) result(T) 
   real(8) T    ! transmission probability / penetrability
   real(8), intent(in) :: Work_func   ! work function of the material/height of the barrier
   real(8), intent(in) :: Ekin  ! electron kinetic energy, perpendicular to the surface component [eV]
   real(8) :: frac, num, denom, denom2
   
   if (Ekin < Work_func) then
      T = 0.0d0
   else
      frac = Work_func / Ekin
      num = 4.0d0*sqrt(1.0d0 - frac)
      denom = 1.0d0 + sqrt(1.0d0 - frac)
      denom2 = denom*denom
      T = num/denom2 ! Eq.(2) in [3]
   endif
end function electron_transmission_probability_step 




pure subroutine set_all_barriers_parameters(used_target)
   type(Matter), intent(inout) :: used_target   ! parameters of the target
   integer :: i, N_mat
   
   ! How many different materials are there:
   N_mat = size(used_target%Material)
   
   ! For all targets, set their barrier parameters:
   do i = 1, N_mat
      call get_barrier_parameters(used_target%Material(i)%Surface_barrier)  ! below
   enddo
end subroutine set_all_barriers_parameters


! Parameters of the barrier for charged particle barrier crossing:
pure subroutine get_barrier_parameters(Surface_barrier)
   type(Emission_barrier), intent(inout) :: Surface_barrier    ! Parameters of the surface barier for a particle emission from the surface of the material
   real(8) :: work_function, bar_height, Em_L, Em_B, Em_ksi, Em_bb, Em_delta
   real(8) :: Em_E1, Em_gam1, Em_gam2, Em_gamma
   real(8) :: eps
   !eps = 1.0d-10
   eps = 1.0d-3*m_tollerance_eps   ! module "Geometries"
   
   work_function = Surface_barrier%Work_func
   bar_height = Surface_barrier%Bar_height
    
   Em_L = Surface_barrier%Surf_bar*1.0d-10     ! Barrier length     [m]
   Em_B = 2.0d0*bar_height - work_function + 2.0d0*sqrt(bar_height*(bar_height - work_function))   ! Eq.(5) [2]
    
   Em_ksi = 0.5d0*sqrt(8.0d0*g_me*Em_L*Em_L*Em_B*g_e/(g_2Pih*g_2Pih)-1.0d0)     ! Eq.(5) [2]
   Em_bb = cosh(g_2Pi*Em_ksi)
   !Em_delta = g_2Pi*Em_L*sqrt(2.0d0*g_me*g_e)/(g_2Pih)
   Em_delta = Em_L*sqrt(2.0d0*g_me*g_e)/g_h     ! under Eq.(7) [2]
    
   ! Coefficients E1 and gamma from Eqs.(7) in [2]:
   Em_E1 = bar_height + 2.0d0*sqrt(bar_height*(bar_height - work_function))*(acosh(Em_bb)/(Em_delta*(sqrt(bar_height) + sqrt(bar_height - work_function)))-1.0d0)
   Em_gam1 = Em_delta*(sqrt(Em_E1) + sqrt(Em_E1 - work_function))
   Em_gam2 = Em_delta*(sqrt(Em_E1) - sqrt(Em_E1 - work_function))
   if ((Em_E1 - work_function) < eps) Em_E1 = work_function + eps
   Em_gamma = (Em_gam1*sinh(Em_gam1) + 2.0d0*Em_gam2*sinh(Em_gam2))/(sqrt(Em_E1*(Em_E1 - work_function))*(Em_bb + cosh(Em_gam1)))
   
   Surface_barrier%E1 = Em_E1
   Surface_barrier%gamma = Em_gamma
end subroutine get_barrier_parameters


! Place all particles back inside the simulation box:
subroutine place_particles_back_into_box(numpar, Prtcls)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   class(Particle), dimension(:), intent(inout) :: Prtcls  ! the particle that may be outside of the simulation box
   integer :: i, siz
   siz = size(Prtcls)
   do i = 1, siz
      call put_back_into_box(numpar, Prtcls(i)) ! below
   enddo
end subroutine place_particles_back_into_box



! Place the particle back inside the simulation box:
subroutine put_back_into_box(numpar, Prtcl, R)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   class(Particle), intent(inout), optional :: Prtcl  ! the particle that may be outside of the simulation box
   real(8), dimension(3), intent(inout), optional :: R    ! only coordinates of the particles
   real(8) :: X, Y, Z, X0, Y0, Z0
   integer :: d, generation
   
   if (present(Prtcl)) then
      generation = Prtcl%generation ! use real generation of the particle
   else
      generation = 1    ! just to go on with calculations
   endif
   
   if (generation > 0) then   ! exclude external incoming particle
      ! If a particle coordinate is outside of the box, place it back inside via following steps:
      if( numpar%periodic_x == 1 ) then    ! periodic along X
         X0 = numpar%box_end_x - numpar%box_start_x    ! shift coordinate system to start of the box
         if (present(Prtcl)) then
            X = Prtcl%R(1) - numpar%box_start_x   ! shift particle with respect to the box border
         elseif (present(R)) then
            X = R(1) - numpar%box_start_x   ! shift particle with respect to the box border
         else
            print*, 'Error in put_back_into_box subroutine: no coordinates provided'
         endif
         d = floor(X/X0)   ! find how many full box sizes are between particle and box border
         if (present(Prtcl)) then
            Prtcl%R(1) = Prtcl%R(1) - d*X0    ! shift the particle back into the box, and return into coordinate system
            Prtcl%R0(1) = Prtcl%R0(1) - d*X0    ! shift the coordinates of the last step too
         elseif (present(R)) then
            R(1) = R(1) - d*X0    ! shift the particle back into the box, and return into coordinate system
         endif
      endif
      ! Analogously for all other coordinates:
      if( numpar%periodic_y == 1 ) then    ! periodic along Y
         Y0 = numpar%box_end_y - numpar%box_start_y
         if (present(Prtcl)) then
            Y = Prtcl%R(2) - numpar%box_start_y
         elseif (present(R)) then
            Y = R(2) - numpar%box_start_y
         endif
         d = floor(Y/Y0)
         if (present(Prtcl)) then
            Prtcl%R(2) = Prtcl%R(2) - d*Y0
            Prtcl%R0(2) = Prtcl%R0(2) - d*Y0
         elseif (present(R)) then
            R(2) = R(2) - d*Y0
         endif
      endif
      if( numpar%periodic_z == 1 ) then    ! periodic along Z
         Z0 = numpar%box_end_z - numpar%box_start_z
         if (present(Prtcl)) then
            Z = Prtcl%R(3) - numpar%box_start_z
         elseif (present(R)) then
            Z = R(3) - numpar%box_start_z
         endif
         d = floor(Z/Z0)
         if (present(Prtcl)) then
            Prtcl%R(3) = Prtcl%R(3) - d*Z0
            Prtcl%R0(3) = Prtcl%R0(3) - d*Z0
         elseif (present(R)) then
            R(3) = R(3) - d*Z0
         endif
      endif
   endif ! (generation > 0)
end subroutine put_back_into_box



subroutine particle_cross_box_boundary(used_target, numpar, N_particles, NOP, Prtcl, type_of_periodicity)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(inout) :: N_particles   ! number of particles of this type
   integer, intent(in) :: NOP   ! index of particle in the array
   class(Particle), dimension(:), intent(inout) :: Prtcl  ! the photon to perform some event
   integer, intent(in), optional :: type_of_periodicity   ! if defined by the user: 0=absorbing; 1=periodic; 2=reflective; 3=white
   !-----------------------------------------------------
   integer :: Box_wall_index ! 1=top (Z_start), 2= bottom(Z_end), 3=left (Y_start), 4=right (Y_end), 5=back (X_start), 6=front (X_end)
   logical :: particle_removed
   
   particle_removed = .false.   ! to start with, the particle is still there
   
   !1) Find which boundary it is crosing:
   Box_wall_index = which_box_wall_crossing(Prtcl(NOP)%R(:), numpar%box_start_x, numpar%box_end_x, &
                               numpar%box_start_y, numpar%box_end_y, numpar%box_start_z, numpar%box_end_z)  ! module "Geometries"
   ! 2) Check if it is an incomming particle, or an excited from the target:
   if (Prtcl(NOP)%generation == 0) then ! external incomming
      ! 3a) Remove external particle when a box boundary is reached:
      call cross_boundary(numpar, N_particles, NOP, Prtcl, 0, 0, particle_removed) ! below
   else ! internal part of the target that got excited:
      ! 3b) Do something, depending on kind of boundaries defined by the user:
      select case (Box_wall_index)  ! 1=top (Z_start), 2= bottom(Z_end), 3=left (Y_start), 4=right (Y_end), 5=back (X_start), 6=front (X_end)
      case (1:2) ! perpendicular to Z
         if (present (type_of_periodicity)) then    ! enforce this condition, no matter what is provided in numpar%periodic_z
!             print*, 'type_of_periodicity=', type_of_periodicity, NOP
            call cross_boundary(numpar, N_particles, NOP, Prtcl, type_of_periodicity, Box_wall_index, particle_removed) ! below
         else   ! follow what is provided in numpar%periodic_z
!             print*, 'type_of_periodicity=', numpar%periodic_z, NOP
            call cross_boundary(numpar, N_particles, NOP, Prtcl, numpar%periodic_z, Box_wall_index, particle_removed) ! below
         endif
      case (3:4) ! perpendicular to Y
         if (present (type_of_periodicity)) then    ! enforce this condition, no matter what is provided in numpar%periodic_y
            call cross_boundary(numpar, N_particles, NOP, Prtcl, type_of_periodicity, Box_wall_index, particle_removed) ! below
         else   ! follow what is provided in numpar%periodic_y
            call cross_boundary(numpar, N_particles, NOP, Prtcl, numpar%periodic_y, Box_wall_index, particle_removed) ! below
         endif
      case (5:6) ! perpendicular to X
         if (present (type_of_periodicity)) then    ! enforce this condition, no matter what is provided in numpar%periodic_x
            call cross_boundary(numpar, N_particles, NOP, Prtcl, type_of_periodicity, Box_wall_index, particle_removed) ! below
         else   ! follow what is provided in numpar%periodic_x
            call cross_boundary(numpar, N_particles, NOP, Prtcl, numpar%periodic_x, Box_wall_index, particle_removed) ! below
         endif
      end select
      
      ! Update the flight time until the next boundary crossing:
      if (.not.particle_removed) then    ! if not absorbing boundary
!         print*, 'PCBB before:', Prtcl(NOP)%Ekin, Prtcl(NOP)%ti, Prtcl(NOP)%t_sc, 'R:', Prtcl(NOP)%R(:), 'V:', Prtcl(NOP)%V(:)
         call update_flight_time_to_boundary(used_target, numpar, Prtcl(NOP))   ! below
!         print*, 'PCBB after:', Prtcl(NOP)%Ekin, Prtcl(NOP)%ti, Prtcl(NOP)%t_sc, 'R:', Prtcl(NOP)%R(:), 'V:', Prtcl(NOP)%V(:)
      endif
   endif
end subroutine particle_cross_box_boundary


subroutine cross_boundary(numpar, N_particles, NOP, Prtcl, type_of_periodicity, Box_wall_index, particle_removed)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(inout) :: N_particles   ! number of particles of this type
   integer, intent(in) :: NOP   ! index of particle in the array
   class(Particle), dimension(:), intent(inout) :: Prtcl  ! the photon to perform some event
   integer, intent(in) :: type_of_periodicity   ! defined by the user: 0=absorbing; 1=periodic; 2=reflective; 3=white
   integer, intent(in) ::  Box_wall_index ! 1=top (Z_start), 2= bottom(Z_end), 3=left (Y_start), 4=right (Y_end), 5=back (X_start), 6=front (X_end)
   logical, intent(inout) :: particle_removed   ! is particle removed or not
   real(8) :: PrtclV, eps
   
   eps = 1.0d-1 * m_tollerance_eps      ! shift size
   
   select case (type_of_periodicity)
   case (1) ! periodic
      call cross_periodic_boundary(numpar, Prtcl(NOP)%R, Box_wall_index) ! below
!       PrtclV = max( SQRT( SUM( Prtcl(NOP)%V(:)*Prtcl(NOP)%V(:) ) ), m_tollerance_eps) ! exclude zero
!       Prtcl(NOP)%R(:) = Prtcl(NOP)%R(:) + eps * Prtcl(NOP)%V(:)/PrtclV     ! to place particle a minimal distance away from the border
   case (2) ! reflective
      call reflect_from_wall(Box_wall_index, Prtcl(NOP)%V)  ! below
      PrtclV = max( SQRT( SUM( Prtcl(NOP)%V(:)*Prtcl(NOP)%V(:) ) ), m_tollerance_eps) ! exclude zero
      Prtcl(NOP)%R(:) = Prtcl(NOP)%R(:) + eps * Prtcl(NOP)%V(:)/PrtclV     ! to place particle a minimal distance away from the border
   case (3) ! white
      call  scatter_from_white_wall(Box_wall_index, Prtcl(NOP)%V)  ! below
      PrtclV = max( SQRT( SUM( Prtcl(NOP)%V(:)*Prtcl(NOP)%V(:) ) ), m_tollerance_eps) ! exclude zero
      Prtcl(NOP)%R(:) = Prtcl(NOP)%R(:) + eps * Prtcl(NOP)%V(:)/PrtclV     ! to place particle a minimal distance away from the border
   case default ! absorbing boundary
      call remove_particle_from_MC_array(N_particles, NOP, Prtcl)  ! below
      particle_removed = .true. ! particle removed
   end select
end subroutine cross_boundary


subroutine cross_periodic_boundary(numpar, R, Box_wall_index)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(3), intent(inout) :: R   ! coordinates
   integer, intent(in) ::  Box_wall_index ! 1=top (Z_start), 2= bottom(Z_end), 3=left (Y_start), 4=right (Y_end), 5=back (X_start), 6=front (X_end)
   ! Move the particle to the other end of the box:
   select case (Box_wall_index)
   case (1) ! move "up"
      R(3) = numpar%box_end_z - m_tollerance_eps
   case (2)  ! move "down"
      R(3) =  numpar%box_start_z + m_tollerance_eps
   case (3)  ! move "right"
      R(2) = numpar%box_end_y - m_tollerance_eps
   case (4)  ! move "left"
      R(2) = numpar%box_start_y + m_tollerance_eps
   case (5)  ! move "forward"
      R(1) = numpar%box_end_x - m_tollerance_eps
   case (6)  ! move "backward"
      R(1) = numpar%box_start_x + m_tollerance_eps
   end select
end subroutine cross_periodic_boundary



subroutine reflect_from_wall(Box_wall_index, V)
   integer, intent(in) :: Box_wall_index    ! index of the plane
   real(8), dimension(3), intent(inout) :: V  ! particle velosity
   select case (Box_wall_index)
   case (5:6) ! perpendicular to X
      V(1) = -V(1)
   case (3:4) ! perpendicular to Y
      V(2) = -V(2)
   case (1:2) ! perpendicular to Z
      V(3) = -V(3)
   endselect
end subroutine reflect_from_wall


subroutine scatter_from_white_wall(Box_wall_index, V)
   integer, intent(in) :: Box_wall_index    ! index of the plane
   real(8), dimension(3), intent(inout) :: V  ! particle velosity
   real(8) :: V_tot, theta, phi
   ! Set random direction of velocity:
   V_tot = sqrt(sum(V(:)*V(:)))
   theta = sample_theta()   ! below
   phi = sample_phi()   ! below
   call  Spherical_to_cartesian(V_tot, phi, theta, V(1), V(2), V(3)) ! module "Geometries"
   ! Now, make sure it scattered off the box wall, not into it:
   select case (Box_wall_index)  ! 1=top (Z_start), 2= bottom(Z_end), 3=left (Y_start), 4=right (Y_end), 5=back (X_start), 6=front (X_end)
   case (1) ! must move "up"
      V(3) = abs(V(3))
   case (2)  ! must move "down"
      V(3) = -abs(V(3))
   case (3)  ! must move "right"
      V(2) = abs(V(2))
   case (4)  ! must move "left"
      V(2) = -abs(V(2))
   case (5)  ! must move "forward"
      V(1) = abs(V(1))
   case (6)  ! must move "backward"
      V(1) = -abs(V(1))
   end select
end subroutine scatter_from_white_wall


function sample_phi() result(phi)
   real(8) phi    ! uniformly sampled angle around X from 0 to 2Pi
   real(8) :: RN
   call random_number(RN)
   phi = RN * g_2Pi
end function sample_phi


pure subroutine check_phi(phi)
   real(8), intent(inout) :: phi    ! ensure phi to be within (0:2Pi)
   real(8) :: sub_phi   ! value to subtract if needed
   sub_phi = FLOOR(phi/g_2Pi)
   phi = phi - sub_phi*g_2Pi
end subroutine check_phi


function sample_theta() result(theta)
   real(8) theta    ! uniformly sampled angle around Z from 0 to Pi
   real(8) :: RN
   call random_number(RN)
   !theta = g_Pi*(-0.5d0 + RN)
   theta = g_Pi*RN
end function sample_theta



function transfered_E_from_theta(Ekin, theta, M_in, mt) result(dE)
   real(8) dE   ! [eV]  according to Eq.(A.25) in 2015-edition of [1]
   real(8), intent(in) :: Ekin  ! [eV] kinetic energy of the incident particle
   real(8), intent(in) :: theta ! scattering angle
   real(8), intent(in) :: M_in, mt  ! incoming particle and target particle masses
   real(8) :: cos_theta, cos_theta2, sin_theta2, mc2, Mct2, Emc, E2mc, W1, W2, EmcMc
   mc2 = rest_energy(M_in)	! mc^2 [eV],  module "Relativity"
   Mct2 = rest_energy(mt)	! Mc^2 [eV],  module "Relativity"
   cos_theta = cos(theta)
   cos_theta2 = cos_theta*cos_theta
   sin_theta2 = 1.0d0 - cos_theta2
   Emc = Ekin + mc2
   E2mc = Ekin + 2.0d0*mc2
   EmcMc = Emc + Mct2
   W1 = Emc*sin_theta2 + Mct2 - cos_theta * sqrt( Mct2*Mct2 - mc2*mc2*sin_theta2 )
   W2 = Ekin*E2mc / ( EmcMc*EmcMc - Ekin*E2mc*cos_theta2 )
   dE = W1*W2
end function transfered_E_from_theta


function cos_theta_from_W(E, W, M_in, mt) result(mu) ! Checked, identical to TREKIS-3 in classical limit
   real(8) mu   ! cos(theta) scattering angle, Eq.(A.23) in 2015-edition of [1]
   real(8), intent(in) :: E, W  ! [eV] incident and transfered energy
   real(8), intent(in) :: M_in, mt  ! [kg] masses of incoming and target particles
   real(8) :: Erest_in, Erest_t, E2mc, EmW, W1, W2, theta
   Erest_in = rest_energy(M_in)	! module "Relativity"
   Erest_t = rest_energy(mt)	! module "Relativity"
   E2mc = E + 2.0d0*Erest_in
   EmW = E - W
   W1 = E*E2mc - W*(E + Erest_in + Erest_t)
   W2 = E*E2mc*EmW*(E2mc - W)
   
   if (W2 > 0.0d0) then
      mu = W1/sqrt(W2)
   else
      mu = 0.0d0
   endif
   ! Slow electrons scattering on phonons may not satisfy the conservation written for atoms:
   if (abs(mu) > 1.0d0) then    ! sample uniformly
      theta = sample_theta()    ! function above
      mu = cos(theta)
   endif
   
end function cos_theta_from_W



function cos_theta_from_W_NORM(E, W, M_in, mt) result(mu) ! Normalize energies (testing precision)
   real(8) mu   ! cos(theta) scattering angle, Eq.(A.23) in 2015-edition of [1]
   real(8), intent(in) :: E, W  ! [eV] incident and transfered energy
   real(8), intent(in) :: M_in, mt  ! [kg] masses of incoming and target particles
   real(8) :: Erest_in, Erest_t, E2mc, EmW, W1, W2, theta, Emc2, Wmc2
   Erest_in = rest_energy(M_in)	! module "Relativity"
   Erest_t = rest_energy(mt)	! module "Relativity"
   
   Emc2 = E/Erest_in
   Wmc2 = W/Erest_in
   
   E2mc = Emc2 + 2.0d0
   EmW = E2mc - Wmc2
   
   W1 = Emc2*E2mc - Wmc2*(Emc2 + Erest_t/Erest_in + 1.0d0)
   W2 = Emc2*E2mc*EmW*(Emc2 - Wmc2)
   
   if (W2 > 0.0d0) then
      mu = W1/sqrt(W2)
   else
      mu = 0.0d0
   endif
   ! Slow electrons scattering on phonons may not satisfy the conservation written for atoms:
   if (abs(mu) > 1.0d0) then    ! sample uniformly
      theta = sample_theta()    ! function above
      mu = cos(theta)
   endif
   
end function cos_theta_from_W_NORM


pure function cos_recoil_from_W(E, W, M_in, mt) result(mu)
   real(8) mu   ! cos(theta) scattering angle, Eq.(A.24) in 2015-edition of [1]
   real(8), intent(in) :: E, W, M_in, mt
   real(8) :: Erest_in, Erest_t, E2mc, W1, W2
   Erest_in = rest_energy(M_in)	! module "Relativity"
   Erest_t = rest_energy(mt)	! module "Relativity"
   E2mc = 2.0d0*Erest_in
   W1 = W*(E + Erest_in + Erest_t)
   W2 = E*(E+E2mc)*W*(W+2.0d0*Erest_t) ! Note a misprint in Eq.(A.24): should be (W+2mtc^2) !
   if (W2 > 0.0d0) then
      mu = W1/sqrt(W2)
   else
      mu = 0.0d0
   endif
end function cos_recoil_from_W



subroutine flight_time_to_boundary(used_target, numpar, Prtcl, T, neutral)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   class(Particle), intent(inout) :: Prtcl      ! undefined particle as an object
   real(8), intent(out) :: T  ! [fs] time to the crossing point
   logical, intent(in), optional :: neutral ! if it's neutral particle, periodic conditions can be simplified
   !-------------------------------------------
   real(8) :: R0(3), R(3), V(3), t_box, t_target
   t_target = 1.0d25
   ! 1) Check how far it is from the border of the simulation box boundary:
   if (present(neutral)) then   ! skip periodic boundary, if user chose to (mainly for photons):
      call where_crossing_box_boundary(Prtcl%R(:), Prtcl%V(:), numpar%box_start_x, numpar%box_end_x, &
            numpar%box_start_y, numpar%box_end_y, numpar%box_start_z, numpar%box_end_z, &
            Prtcl%generation, numpar%periodic_x, numpar%periodic_y, numpar%periodic_z, t_box, neutral) ! module "Geometries"
   else
      call where_crossing_box_boundary(Prtcl%R(:), Prtcl%V(:), numpar%box_start_x, numpar%box_end_x, &
            numpar%box_start_y, numpar%box_end_y, numpar%box_start_z, numpar%box_end_z, &
            Prtcl%generation, numpar%periodic_x, numpar%periodic_y, numpar%periodic_z, t_box) ! module "Geometries"
   endif
   
!    print*, 'x', numpar%box_start_x, numpar%box_end_x
!    print*, 'y', numpar%box_start_y, numpar%box_end_y
!    print*, 'z', numpar%box_start_z, numpar%box_end_z
!    print*, 't_box=', t_box
   
   ! 2) Figure out, within which target the particle is:
   if (Prtcl%in_target > 0) then ! exclude special case: vacuum, out of any target
      ASSOCIATE (ARRAY => used_target%Geom(Prtcl%in_target)%Set)	! that's the syntax to use when passing polimorphic arrays into subroutines
       select type (ARRAY)	! Depending on the type of used_target%Geom(i)%Set
         type is (Rectangle)
            ! Transform into coordinates of the target center:
            call shift_coordinate_system(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), ARRAY%X, ARRAY%Y, ARRAY%Z, R0(1), R0(2), R0(3))   ! module "Geometries"
            call rotate_coordinate_system(R0(1), R0(2), R0(3), R(1), R(2), R(3), ARRAY%angle_x, ARRAY%angle_y, ARRAY%angle_z)   ! module "Geometries"
            ! Transform velosities into the coordinate system of the target center:
            call rotate_coordinate_system(Prtcl%V(1), Prtcl%V(2), Prtcl%V(3), V(1), V(2), V(3), ARRAY%angle_x, ARRAY%angle_y, ARRAY%angle_z)   ! module "Geometries"
            
!              print*, 'R', R(:)
!              print*, 'V', V(:)
!              print*, ARRAY%Zstart, ARRAY%Zend
!              print*, ARRAY%X, ARRAY%Y, ARRAY%Z
!              pause 'flight_time_to_boundary'
            
            ! Check where the boundary is crossed:
            call where_crossing_box_boundary(R(:),V(:), ARRAY%Xstart, ARRAY%Xend, &
                   ARRAY%Ystart, ARRAY%Yend, ARRAY%Zstart, ARRAY%Zend, &
                   Prtcl%generation, 0, 0, 0, t_target) ! module "Geometries"
                   !Prtcl%generation, numpar%periodic_x, numpar%periodic_y, numpar%periodic_z, t_target) ! module "Geometries"
            
!              print*, 't_target=', t_target, R(:), V(:)
         
         type is (Sphere)
            ! Transform into coordinates of the target center:
            call shift_coordinate_system(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), ARRAY%X, ARRAY%Y, ARRAY%Z, R0(1), R0(2), R0(3))   ! module "Geometries"
            call intersection_with_sphere(R0(1), R0(2), R0(3), Prtcl%V(1), Prtcl%V(2), Prtcl%V(3), ARRAY%R, T)   ! module "Geometries"
         type is (Sphere_segment)
            ! Not done yet! DO NOT USE!!
         type is (Cylinder)
            ! Transform into coordinates of the target center:
            call shift_coordinate_system(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), ARRAY%X, ARRAY%Y, ARRAY%Z, R0(1), R0(2), R0(3))   ! module "Geometries"
            call rotate_coordinate_system(R0(1), R0(2), R0(3), R(1), R(2), R(3), ARRAY%angle_x, ARRAY%angle_y, ARRAY%angle_z)   ! module "Geometries"
            ! Transform velosities into the coordinate system of the target center:
            call rotate_coordinate_system(Prtcl%V(1), Prtcl%V(2), Prtcl%V(3), V(1), V(2), V(3), ARRAY%angle_x, ARRAY%angle_y, ARRAY%angle_z)   ! module "Geometries"
            ! Check where the boundary is crossed:
            call where_crossing_cylinder(R, V, ARRAY%L_start, ARRAY%L_end, ARRAY%R, T)    ! module "Geometries"
         type is (Cylinder_segment)
            ! Not done yet! DO NOT USE!!
       endselect
      END ASSOCIATE
   endif ! (Prtcl%in_target > 0)
   T = min(t_box, t_target) ! the closest surface to be crossed
   
!    print*, 'flight_time_to_boundary', t_box, t_target
!    pause 'flight_time_to_boundary END'
   
end subroutine flight_time_to_boundary


! Define normal to arbitrary surface at the point of particle impact:
 subroutine define_normal_to_surface(used_target,  Prtcl, n, message_in, INFO)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   class(Particle), intent(inout) :: Prtcl      ! undefined particle as an object
   real(8), dimension(3), intent(out) :: n  ! normal vector to the surface at the given point
   character(*), intent(in) :: message_in   ! for checking, message to dysplay
   integer, intent(out), optional :: INFO   ! print out error info
   !----------------------------
   real(8), dimension(3) :: R_shift
   real(8) :: eps, Xc, Yc, Zc, Xcr, Ycr, Zcr, n_abs, PrtclV
   integer :: i, N_mat, ind_boundary
   
   if (present(INFO)) INFO = 0  ! no error at the start
   N_mat = size(used_target%Material)    ! how many different materials we have
   !eps = 1.0d-6 ! precision
   eps = m_tollerance_eps   ! module "Geometries"
   
   ! 1) Find which target's boundary the particle is crossing:
   if (Prtcl%in_target == 0) then ! particle coming from vacuum, find which material it enters:
      !R_shift = 1.0d-7 * Prtcl%V(:)     ! to place particle inside of the material
      PrtclV = max(SQRT( SUM( Prtcl%V(:)*Prtcl%V(:) ) ), m_tollerance_eps) ! exclude zero
      R_shift = m_tollerance_eps * Prtcl%V(:)/PrtclV     ! to place particle inside of the material
      
      ! Update particle's material index according to the new material it enters:
      call find_the_target(used_target, Prtcl, R_shift) ! below
   endif
   
   ! 2) Find which surface it is corssing:
   ASSOCIATE (ARRAY => used_target%Geom(Prtcl%in_target)%Set)	! that's the syntax to use when passing polimorphic arrays into subroutines
       select type (ARRAY)	! Depending on the type of used_target%Geom(i)%Set
         type is (Rectangle)
            
            ! First, transfere coordinates into the center-of-the-target system:
            call shift_coordinate_system(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), ARRAY%X, ARRAY%Y, ARRAY%Z, Xc, Yc, Zc)    ! module "Geometries"

!             print*, 'define_normal_to_surface 0:', Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), Xc, Yc, Zc
            
            ! Then, rotate it accordingly to the target orientation:
            call rotate_coordinate_system(Xc, Yc, Zc, Xcr, Ycr, Zcr, ARRAY%angle_x, ARRAY%angle_y, ARRAY%angle_z)    ! module "Geometries"
            
!             print*, 'define_normal_to_surface 1:', Xc, Yc, Zc
            
            ! To start with:
            n = 0.0d0
            
            ! Check which surface it is crossing:
            if ( (abs(Xcr - ARRAY%Xstart) < eps) .or. (abs(Xcr - ARRAY%Xend) < eps ) ) then ! crossing X
               n(2:3) = 0.0d0
               n(1) = 1.0d0
            endif
            if ( (abs(Ycr - ARRAY%Ystart) < eps) .or. (abs(Ycr - ARRAY%Yend) < eps ) ) then ! crossing Y
               n(1) = 0.0d0
               n(2) = 1.0d0
               n(3) = 0.0d0
            endif
            if ( (abs(Zcr - ARRAY%Zstart) < eps) .or. (abs(Zcr - ARRAY%Zend) < eps ) ) then ! crossing Z
               n(1:2) = 0.0d0
               n(3) = 1.0d0
            endif
            
            ! Check consistency:
            if (all( abs(n(:)) < m_tollerance_eps)) then ! try another way: find the nearest boundary
               if (present(INFO)) INFO = 1  ! particle is too far from a rectangle bondary
!                print*, 'ERROR in define_normal_to_surface: ', trim(adjustl(message_in))
!                print*, 'Cannot find the surface the particles is crossing'
!                print*, 'n=', n
!                print*, 'R(rotated)=', Xcr, Ycr, Zcr
!                print*, 'R(original)', Prtcl%R(:)
!                print*, 'V=', Prtcl%V(:)
!                print*, 'Trying a different way of boundary finding...'
!                pause 'define_normal_to_surface'
               call find_closest_boundary(Xcr, Ycr, Zcr, &
                    ARRAY%Xstart, ARRAY%Ystart, ARRAY%Zstart, ARRAY%Xend, ARRAY%Yend, ARRAY%Zend, ind_boundary)  ! below
               select case (ind_boundary)
               case (1,2)   ! closest boundary is along X
                  n(2:3) = 0.0d0
                  n(1) = 1.0d0
               case (3,4)   ! closest boundary is along Y
                  n(1) = 0.0d0
                  n(2) = 1.0d0
                  n(3) = 0.0d0
               case (5,6)   ! closest boundary is along Z
                  n(1:2) = 0.0d0
                  n(3) = 1.0d0
               endselect
            endif
            
!             print*, 'define_normal_to_surface 2', n, Zcr - ARRAY%Zstart, Zcr - ARRAY%Zend
            
            ! Rotate it to the lab.coordinate system:
            call rotate_coordinate_system(n(1), n(2), n(3), n(1), n(2), n(3), ARRAY%angle_x, ARRAY%angle_y, ARRAY%angle_z)  ! module "Geometries"
            
!             print*, 'define_normal_to_surface END:', n
            
         type is (Sphere)
            ! Vector connecting the point of impact with the center of the sphere:
            n(1) = Prtcl%R(1) - ARRAY%X
            n(2) = Prtcl%R(2) - ARRAY%Y
            n(3) = Prtcl%R(3) - ARRAY%Z
            ! Normalize to 1:
            n(:) = n(:) / ARRAY%R
         
         type is (Sphere_segment)
            ! Not done yet! DO NOT USE!!
         
         type is (Cylinder)
            ! First, transfere coordinates into the center-of-the-target system:
            call shift_coordinate_system(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), ARRAY%X, ARRAY%Y, ARRAY%Z, Xc, Yc, Zc)    ! module "Geometries"
            ! Then, rotate it accordingly to the target orientation:
            call rotate_coordinate_system(Xc, Yc, Zc, Xcr, Ycr, Zcr, ARRAY%angle_x, ARRAY%angle_y, ARRAY%angle_z)    ! module "Geometries"
            ! Check if it is impact of the base of the cyllinder:
            if ( (abs(Zcr - ARRAY%L_start) < eps) .or. (abs(Zcr - ARRAY%L_end) < eps ) ) then ! crossing top or bottom surface
               ! Cylinder base is perpendicular to Z in the cylinder coordinate system:
               n(1:2) = 0.0d0
               n(3) = 1.0d0
            else    ! crossing cyllinder surface:
               ! Vector connecting the point of impact with the center of the cylinder:
               n(1) = Prtcl%R(1) - ARRAY%X
               n(2) = Prtcl%R(2) - ARRAY%Y
               n(3) = 0.0d0
               ! Normalize to 1:
               n(:) = n(:) / ARRAY%R
            endif
            ! Rotate it to the lab.coordinate system:
            call rotate_coordinate_system(n(1), n(2), n(3), n(1), n(2), n(3), ARRAY%angle_x, ARRAY%angle_y, ARRAY%angle_z)  ! module "Geometries"
        
         type is (Cylinder_segment)
            ! Not done yet! DO NOT USE!!
       endselect
      END ASSOCIATE
end subroutine define_normal_to_surface


subroutine find_closest_boundary(X, Y, Z, Xstart, Ystart, Zstart, Xend, Yend, Zend, ind_boundary)
   real(8), intent(in) :: X, Y, Z   ! particle coordinates
   real(8), intent(in) :: Xstart, Ystart, Zstart, Xend, Yend, Zend  ! parallelepiped boundaries
   integer, intent(out) :: ind_boundary ! index of closest boundary: 1=Xstart, 3=Ystart, 5=Zstart, 2=Xend, 4=Yend, 6=Zend
   real(8) :: d, temp
   
   ind_boundary = 1 ! to start with
   d = abs(X - Xstart)
   temp = abs(X - Xend)
   if (temp < d) then
      d = temp
      ind_boundary = 2
   endif
   temp = abs(Y - Ystart)
   if (temp < d) then
      d = temp
      ind_boundary = 3
   endif
   temp = abs(Y - Yend)
   if (temp < d) then
      d = temp
      ind_boundary = 4
   endif
   temp = abs(Z - Zstart)
   if (temp < d) then
      d = temp
      ind_boundary = 5
   endif
   temp = abs(Z - Zend)
   if (temp < d) then
      d = temp
      ind_boundary = 6
   endif
end subroutine find_closest_boundary



subroutine reflection_from_surface(V, n)
   real(8), dimension(3), intent(inout) :: V    ! velosity
   real(8), dimension(3), intent(in) :: n    ! normal defining the surface at the point of scattering
   real(8) :: Vabs0, Vabs1, eps
   real(8), dimension(3) :: Vn, V_save
   eps = 1.0d-10
   V_save = V
   
   ! To make sure the normalization is correct:
   Vabs0 = sqrt( SUM( V(:)*V(:) ) )  ! initial absolute velosity
   
   ! Get the multiplication:
   Vn = SUM( V(:)*n(:) )
   ! Change velosity during reflection:
   V(:) = V(:) - 2.0d0*Vn*n(:)  ! assume absolutely elastic reflection
   
   ! Make sure normalization is conserved:
   Vabs1 = sqrt( SUM( V(:)*V(:) ) )    ! new absolute velosity
   
   if (abs(Vabs1-Vabs0) > eps*Vabs0) then
      V(:) = V(:) / Vabs1*Vabs0    ! renormalize to conserve exactly
!     else
!          print*, 'Error in reflection_from_surface:', Vabs1,  V(:)
!          print*, 'reflection_from_surface:', Vabs1,  V(:)
!          print*, 'V0', Vabs0, V_save
!          print*, 'n', n
!         V(:) = V(:) / Vabs1*Vabs0    ! renormalize to conserve exactly
   endif
end subroutine reflection_from_surface



!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! Subrtines related to targets:

subroutine Find_starting_targets(used_target, numpar, bunch, MC)
   type(Matter), intent(inout) :: used_target	! parameters of the target
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(Radiation_param), dimension(:), intent(inout) :: bunch	! incomming radiation
   type(MC_arrays), dimension(:), intent(inout) :: MC	! all MC arrays for particles: photons, electrons and holes; size equals to number of iterations
   !----------------------------------------
   real(8), dimension(3) :: R_shift
   real(8) :: Vabs
   integer :: i, j, N_size
   N_size = size(bunch)        ! initial size of all MC particle arrays
   
   ! Arrays of particles inside of MC iterations:
   do i = 1, max(numpar%NMC,1)
      if (allocated(MC(i)%MC_Photons)) then
         do j = 1, N_size
            if (MC(i)%MC_Photons(j)%active) then   ! only for active particles
               !R_shift(:) = 1.0d-7 * MC(i)%MC_Photons(j)%V(:)     ! to place particle inside of the material
               Vabs = max( SQRT( SUM(MC(i)%MC_Photons(j)%V(:)*MC(i)%MC_Photons(j)%V(:)) ), m_tollerance_eps)
               R_shift(:) = m_tollerance_eps * MC(i)%MC_Photons(j)%V(:)/Vabs     ! to place particle inside of the material
               ! Update particle's material index according to the new material it enters:
               call find_the_target(used_target, MC(i)%MC_Photons(j), R_shift) ! below
!                print*, 'ph', i, j, MC(i)%MC_Photons(j)%in_target
            endif
         enddo
      endif
      
      if (allocated(MC(i)%MC_Electrons)) then
         do j = 1, N_size
            if (MC(i)%MC_Electrons(j)%active) then   ! only for active particles
               !R_shift(:) = 1.0d-7 * MC(i)%MC_Electrons(j)%V(:)     ! to place particle inside of the material
               Vabs = max( SQRT( SUM(MC(i)%MC_Electrons(j)%V(:)*MC(i)%MC_Electrons(j)%V(:)) ), m_tollerance_eps)
               R_shift(:) = m_tollerance_eps * MC(i)%MC_Electrons(j)%V(:)/Vabs     ! to place particle inside of the material
               
               ! Update particle's material index according to the new material it enters:
               call find_the_target(used_target, MC(i)%MC_Electrons(j), R_shift) ! below
            endif
         enddo
      endif
      
      if (allocated(MC(i)%MC_Positrons)) then
         do j = 1, N_size
            if (MC(i)%MC_Positrons(j)%active) then   ! only for active particles
               !R_shift(:) = 1.0d-7 * MC(i)%MC_Positrons(j)%V(:)     ! to place particle inside of the material
               Vabs = max( SQRT( SUM(MC(i)%MC_Positrons(j)%V(:)*MC(i)%MC_Positrons(j)%V(:)) ), m_tollerance_eps)
               R_shift(:) = m_tollerance_eps * MC(i)%MC_Positrons(j)%V(:)/Vabs     ! to place particle inside of the material
               ! Update particle's material index according to the new material it enters:
               call find_the_target(used_target, MC(i)%MC_Positrons(j), R_shift) ! below
            endif
         enddo
      endif
      
      if (allocated(MC(i)%MC_Holes)) then
         do j = 1, N_size
            if (MC(i)%MC_Holes(j)%active) then   ! only for active particles
               if (MC(i)%MC_Holes(j)%valent) then   ! valence holes can move
                  !R_shift(:) = 1.0d-7 * MC(i)%MC_Holes(j)%V(:)     ! to place particle inside of the material
                  Vabs = max( SQRT( SUM(MC(i)%MC_Holes(j)%V(:)*MC(i)%MC_Holes(j)%V(:)) ), m_tollerance_eps)
                  R_shift(:) = m_tollerance_eps * MC(i)%MC_Holes(j)%V(:)/Vabs     ! to place particle inside of the material
                  ! Update particle's material index according to the new material it enters:
                  call find_the_target(used_target, MC(i)%MC_Holes(j), R_shift) ! below
               else ! no shift for immobile deep-shell holes
                  ! Update particle's material index according to the new material it enters:
                  call find_the_target(used_target, MC(i)%MC_Holes(j)) ! below
               endif
            endif
         enddo
      endif
      
      if (allocated(MC(i)%MC_SHIs)) then
         do j = 1, N_size
            if (MC(i)%MC_SHIs(j)%active) then   ! only for active particles
               !R_shift(:) = 1.0d-7 * MC(i)%MC_SHIs(j)%V(:)     ! to place particle inside of the material
               Vabs = max( SQRT( SUM(MC(i)%MC_SHIs(j)%V(:)*MC(i)%MC_SHIs(j)%V(:)) ), m_tollerance_eps)
               R_shift(:) = m_tollerance_eps * MC(i)%MC_SHIs(j)%V(:)/Vabs     ! to place particle inside of the material
               ! Update particle's material index according to the new material it enters:
               call find_the_target(used_target, MC(i)%MC_SHIs(j), R_shift) ! below
            endif
         enddo
      endif
      
      if (allocated(MC(i)%MC_Atoms_events)) then
         do j = 1, N_size
            if (MC(i)%MC_Atoms_events(j)%active) then   ! only for active particles
               !R_shift(:) = 1.0d-7 * MC(i)%MC_Atoms_events(j)%V(:)     ! to place particle inside of the material
               Vabs = max( SQRT( SUM(MC(i)%MC_Atoms_events(j)%V(:)*MC(i)%MC_Atoms_events(j)%V(:)) ), m_tollerance_eps)
               R_shift(:) = m_tollerance_eps * MC(i)%MC_Atoms_events(j)%V(:)/Vabs     ! to place particle inside of the material
               ! Update particle's material index according to the new material it enters:
               call find_the_target(used_target, MC(i)%MC_Atoms_events(j), R_shift) ! below
            endif
         enddo
      endif
   enddo
end subroutine Find_starting_targets



! Identifying the target within which the particle currently is:
pure subroutine find_the_target(used_target, Prtcl, R_shift, tar_ind) ! Find within which given particle is now
   type(Matter), intent(in) :: used_target   ! parameters of the target
   class(Particle), intent(inout) :: Prtcl      ! undefined particle as an object
   real(8), dimension(3), intent(in), optional :: R_shift   ! move particle deeper by this amount
   integer, intent(out), optional :: tar_ind    ! to print out the target index, instead of updating it in the Prtcl properties
   !-------------------------------------------
   real(8), dimension(3) :: R_add
   integer :: siz, i, target_ind
   logical :: found_target
   siz = size(used_target%Material)    ! how many materials we have in this target
   found_target = .false.   ! to start with
   
   ! Start from the vacuum:
   target_ind = 0
   
   if (present(R_shift)) then   ! shift the particle to these coordinates:
      R_add = Prtcl%R + R_shift
   else ! use the actual particles coordinates
      R_add = Prtcl%R
   endif
   
   TRGT:do i = 1, siz    ! for all materials
      ! Once we know the index, lets use the type for the geometry:
      ASSOCIATE (ARRAY => used_target%Geom(i)%Set)	! that's the syntax to use when passing polimorphic arrays into subroutines
      select type (ARRAY)   ! Depending on the type of used_target%Geom(i)%Set
         type is (Rectangle)
         !if (it_is_within_rectangle(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), ARRAY)) then  ! the particle is within this target
         if (it_is_within_rectangle(R_add(1), R_add(2), R_add(3), ARRAY)) then  ! the particle is within this target
            target_ind = i
            found_target = .true. !  the particle is inside this target
         endif
       type is (Sphere)
         !if (it_is_within_sphere(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), ARRAY)) then  ! the particle is within this target
         if (it_is_within_sphere(R_add(1), R_add(2), R_add(3), ARRAY)) then  ! the particle is within this target
            target_ind = i
            found_target = .true. !  the particle is inside this target
         endif
      type is (Sphere_segment)
         !if (it_is_within_sphere_segment(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), ARRAY)) then  ! the particle is within this target
         if (it_is_within_sphere_segment(R_add(1), R_add(2), R_add(3), ARRAY)) then  ! the particle is within this target
            target_ind = i
            found_target = .true. ! the particle is inside this target
         endif
      type is (Cylinder)
         !if (it_is_within_cylinder(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), ARRAY)) then  ! the particle is within this target
         if (it_is_within_cylinder(R_add(1), R_add(2), R_add(3), ARRAY)) then  ! the particle is within this target
            target_ind = i
            found_target = .true. ! the particle is inside this target
         endif
      type is (Cylinder_segment)
         !if (it_is_within_cylinder_segment(Prtcl%R(1), Prtcl%R(2), Prtcl%R(3), ARRAY)) then  ! the particle is within this target
         if (it_is_within_cylinder_segment(R_add(1), R_add(2), R_add(3), ARRAY)) then  ! the particle is within this target
            target_ind = i
            found_target = .true. ! the particle is inside this target
         endif
      endselect
      END ASSOCIATE
      if (found_target) exit TRGT   ! no need to check further if we found the target within which the particle is
   enddo TRGT
   
   ! Where to pass the found index of the target:
   if (present (tar_ind)) then ! printout the target index
      tar_ind = target_ind
   else ! update particle index
      Prtcl%in_target = target_ind
   endif
end subroutine find_the_target


pure function out_of_simulation_box(X, Y, Z, box_start_x, box_end_x, box_start_y, box_end_y, box_start_z, box_end_z) result(out)
   logical :: out
   real(8), intent(in) :: X, Y, Z   ! particle coordinates
   real(8), intent(in) :: box_start_x, box_end_x, box_start_y, box_end_y, box_start_z, box_end_z    ! box sizes
   logical :: within_x, within_y, within_z
   real(8) :: eps   ! precision
   !eps = 1.0d-7 ! how close to the boundary, to consider that it is crossing
   eps = m_tollerance_eps
   within_x = ((X > box_start_x + eps) .and. (X < box_end_x - eps))
   within_y = ((Y > box_start_y + eps) .and. (Y < box_end_y - eps))
   within_z = ((Z > box_start_z + eps) .and. (Z < box_end_z - eps))
   out = .not.(within_x .and. within_y .and. within_z) ! all coordinates must be within this target to be true
end function out_of_simulation_box


pure function it_is_within_rectangle(X, Y, Z, Geom_object) result(is_or_not)
   logical :: is_or_not ! are the particles coordinates within the given target or not?
   real(8), intent(in) :: X, Y, Z   ! particles coordinates [A]
   type(Rectangle), intent(in) :: Geom_object  ! the target geometry array for the case of rectangle
   !---------------------------------------------------
   real(8) :: Xc, Yc, Zc, Xcr, Ycr, Zcr
   logical :: within_x, within_y, within_z
   
   ! First, transfere coordinates into the center-of-the-target system:
   call shift_coordinate_system(X, Y, Z, Geom_object%X, Geom_object%Y, Geom_object%Z, Xc, Yc, Zc)    ! module "Geometries"
   ! Then, rotate it accordingly to the target orientation:
   call rotate_coordinate_system(Xc, Yc, Zc, Xcr, Ycr, Zcr, Geom_object%angle_x, Geom_object%angle_y, Geom_object%angle_z)    ! module "Geometries"
   ! Now, compare if the coordinates are within the borders of the target:
   within_x = ((Xcr >= Geom_object%Xstart) .and. (Xcr <= Geom_object%Xend))
   within_y = ((Ycr >= Geom_object%Ystart) .and. (Ycr <= Geom_object%Yend))
   within_z = ((Zcr >= Geom_object%Zstart) .and. (Zcr <= Geom_object%Zend))
   ! Combine them together:
   is_or_not = (within_x .and. within_y .and. within_z) ! all coordinates must be within this target to be true
end function it_is_within_rectangle


pure function it_is_within_sphere(X, Y, Z, Geom_object) result(is_or_not)
   logical :: is_or_not ! are the particles coordinates within the given target or not?
   real(8), intent(in) :: X, Y, Z   ! particles coordinates [A]
   type(Sphere), intent(in) :: Geom_object  ! the target geometry array for the case of rectangle
   !---------------------------------------------------
   real(8) :: Xc, Yc, Zc, R
   ! First, transfere coordinates into the center-of-the-target system:
   call shift_coordinate_system(X, Y, Z, Geom_object%X, Geom_object%Y, Geom_object%Z, Xc, Yc, Zc)    ! module "Geometries"
   ! Get the radial coordinate:
   R = sqrt(Xc*Xc + Yc*Yc + Zc*Zc)
   ! Now, compare if the particle is within the borders of the sphere:
   is_or_not = (R <= Geom_object%R)
end function it_is_within_sphere


pure function it_is_within_sphere_segment(X, Y, Z, Geom_object) result(is_or_not)
   logical :: is_or_not ! are the particles coordinates within the given target or not?
   real(8), intent(in) :: X, Y, Z   ! particles coordinates [A]
   type(Sphere_segment), intent(in) :: Geom_object  ! the target geometry array for the case of rectangle
   !---------------------------------------------------
   real(8) :: Xc, Yc, Zc, Xcr, Ycr, Zcr, R, theta, phi
   logical :: within_R, within_theta, within_phi
   
   ! First, transfere coordinates into the center-of-the-target system:
   call shift_coordinate_system(X, Y, Z, Geom_object%X, Geom_object%Y, Geom_object%Z, Xc, Yc, Zc)    ! module "Geometries"
   ! Convert coordinates into spherical ones:
   call Cartesian_to_spherical(Xc, Yc, Zc, R, theta, phi)   ! module "Geometries"
   ! Now, compare if the coordinates are within the borders of the target:
   within_R = ((R >= Geom_object%R_start) .and. (R <= Geom_object%R_end))
   within_phi = ((phi >= Geom_object%phi_start) .and. (theta <= Geom_object%phi_end))
   within_theta = ((theta >= Geom_object%theta_start) .and. (theta <= Geom_object%theta_end))
   ! Combine them together:
   is_or_not = (within_R .and. within_phi .and. within_theta) ! all coordinates must be within this target to be true
end function it_is_within_sphere_segment


pure function it_is_within_cylinder(X, Y, Z, Geom_object) result(is_or_not)
   logical :: is_or_not ! are the particles coordinates within the given target or not?
   real(8), intent(in) :: X, Y, Z   ! particles coordinates [A]
   type(Cylinder), intent(in) :: Geom_object  ! the target geometry array for the case of rectangle
   !---------------------------------------------------
   real(8) :: Xc, Yc, Zc, Xcr, Ycr, Zcr, R, L, theta
   logical :: within_R, within_L
   
   ! First, transfere coordinates into the center-of-the-target system:
   call shift_coordinate_system(X, Y, Z, Geom_object%X, Geom_object%Y, Geom_object%Z, Xc, Yc, Zc)    ! module "Geometries"
   ! Then, rotate it accordingly to the target orientation:
   call rotate_coordinate_system(Xc, Yc, Zc, Xcr, Ycr, Zcr, Geom_object%angle_x, Geom_object%angle_y, Geom_object%angle_z)    ! module "Geometries"
   ! Convert it into cylindrical coordinates:
   call Cartesian_to_cylindrical(Xcr, Ycr, Zcr, R, L, theta)     ! module "Geometries"
   ! Now, compare if the coordinates are within the borders of the target:
   within_R = (R <= Geom_object%R)
   within_L = ((L >= Geom_object%L_start) .and. (L <= Geom_object%L_end))
   ! Combine them together:
   is_or_not = (within_R .and. within_L) ! all coordinates must be within this target to be true
end function it_is_within_cylinder


pure function it_is_within_cylinder_segment(X, Y, Z, Geom_object) result(is_or_not)
   logical :: is_or_not ! are the particles coordinates within the given target or not?
   real(8), intent(in) :: X, Y, Z   ! particles coordinates [A]
   type(Cylinder_segment), intent(in) :: Geom_object  ! the target geometry array for the case of rectangle
   !---------------------------------------------------
   real(8) :: Xc, Yc, Zc, Xcr, Ycr, Zcr, R, L, theta
   logical :: within_R, within_L, within_theta
   
   ! First, transfere coordinates into the center-of-the-target system:
   call shift_coordinate_system(X, Y, Z, Geom_object%X, Geom_object%Y, Geom_object%Z, Xc, Yc, Zc)    ! module "Geometries"
   ! Then, rotate it accordingly to the target orientation:
   call rotate_coordinate_system(Xc, Yc, Zc, Xcr, Ycr, Zcr, Geom_object%angle_x, Geom_object%angle_y, Geom_object%angle_z)    ! module "Geometries"
   ! Convert it into cylindrical coordinates:
   call Cartesian_to_cylindrical(Xcr, Ycr, Zcr, R, L, theta)     ! module "Geometries"
   ! Now, compare if the coordinates are within the borders of the target:
   within_R = ((R >= Geom_object%R_start) .and. (R <= Geom_object%R_end))
   within_L = ((L >= Geom_object%L_start) .and. (L <= Geom_object%L_end))
   within_theta = ((theta >= Geom_object%theta_start) .and. (theta <= Geom_object%theta_end))
   ! Combine them together:
   is_or_not = (within_R .and. within_L .and. within_theta) ! all coordinates must be within this target to be true
end function it_is_within_cylinder_segment



!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! Subroutines related to flight times of particles:


subroutine update_flight_time_to_boundary(used_target, numpar, Prtcl) ! Find within which given particle is now
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   class(Particle), intent(inout) :: Prtcl    ! undefined particle as an object
   !-------------------------------------------
   real(8) :: T  ! time and distance to the cross point
   
   ! Check how far it is from the border of the simulation box boundary:
   call flight_time_to_boundary(used_target, numpar, Prtcl, T) ! module "MC_general_tools"
   ! The time of the next event (for boundary crossing):
   Prtcl%ti = Prtcl%t0 + T  ! [fs]

   ! Check which event is shorter: boundary crossing or scattering:
   if (Prtcl%ti >= Prtcl%t_sc) Prtcl%ti = Prtcl%t_sc    ! the next event will be scattering
   
   if ( Prtcl%ti < Prtcl%t0 ) then
      print*, 'Problem in (update_flight_time_to_boundary)', Prtcl%ti, Prtcl%t0, T
   endif
   
end subroutine update_flight_time_to_boundary


subroutine get_photon_flight_time(used_target, numpar, Prtcl) ! Find within which given particle is now
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Photon), intent(inout) :: Prtcl        ! photon as an object
   !-------------------------------------------
   real(8) :: T  ! time and distance to the cross point
   real(8) :: CS_total, MFP_total  ! total cross section and MFP
   real(8) :: Sampled_path
   real(8) :: V ! [A/fs] velosity
   type(Target_atoms), pointer :: matter
   
   ! Check how far it is from the border of the simulation box boundary:
   if (Prtcl%generation == 0) then  ! it is removed even from periodic boundaries:
      call flight_time_to_boundary(used_target, numpar, Prtcl, T) ! module "MC_general_tools"
   else ! use simple periodic boundaries:
      call flight_time_to_boundary(used_target, numpar, Prtcl, T, neutral=.true.) ! module "MC_general_tools"
   endif
   
   ! The time of the next event (for boundary crossing):
   Prtcl%ti = Prtcl%t0 + T  ! [fs]
   ! Now check the time to a scattering event:
   if (Prtcl%in_target > 0) then    ! it is not in vacuum
      ! Properties of the target material, inside of which the particle is:
      matter => used_target%Material(Prtcl%in_target)
      ! a) get the total cross sections from all possible scattering channels:
      CS_total = total_CS_from_chennels(Prtcl%Ekin, matter%Ph_absorption_total%E, matter%Ph_absorption_total%Total, &
                    E_array2=matter%Ph_coherent_total%E, CS_array2=matter%Ph_coherent_total%Total, &
                    E_array3=matter%Ph_Compton_total%E, CS_array3=matter%Ph_Compton_total%Total, &
                    E_array4=matter%Ph_pair_total%E, CS_array4=matter%Ph_pair_total%Total) ! module "CS_general_tools"
      ! b) get the total mean free path:
      MFP_total = MFP_from_sigma(CS_total, matter%At_Dens)  ! module "CS_general_tools"
      ! c) Sample flight path from the mean free path:
      ! In case user provided photon attenuation length, redefine the time of flight accordingly:
      if ( (numpar%Ph_att_eff >= 0.0d0) .and. (Prtcl%generation == 0)) then ! only for external photons (generation = 0):
         Sampled_path = sample_Poisson(numpar%Ph_att_eff)  ! module "Little_subroutines"   
!          print*, 'get_photon_flight_time, Ph_att_eff', numpar%Ph_att_eff, Sampled_path
      else  ! use default MFP from EPDL:
         Sampled_path = sample_Poisson(MFP_total)  ! module "Little_subroutines"   
      endif
      ! d) get the total flight time:
      T = Time_from_MFP(Sampled_path, g_c_Afs) ! [fs] module "CS_general_tools"
      ! e) The time of the next event (for scattering):
      Prtcl%t_sc = Prtcl%t0 + T  ! [fs]
   else ! not inside of any target
      Prtcl%t_sc = 1.0d25    ! no scattering in vacuum
   endif
   
!    print*, 'get_photon_flight_time, t_sc', Prtcl%ti, Prtcl%t_sc

   ! Check which event is shorter: boundary crossing or scattering:
   if (Prtcl%ti >= Prtcl%t_sc) Prtcl%ti = Prtcl%t_sc    ! the next event will be scattering
      
   ! If user provided cut-off energy, check if photon crossed it:
   call check_cut_off(numpar%Ph_Cutoff, Prtcl)   ! below
   
   nullify(matter)
end subroutine get_photon_flight_time



subroutine get_electron_flight_time(used_target, numpar, Prtcl, MD_supce, E_e, no_scatternig) ! Find within which given particle is now
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Electron), intent(inout) :: Prtcl   ! electron as an object
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_e  ! data to pass to MD later
   logical, optional :: no_scatternig       ! recalculate the scattering event, or leave the old time
   !-------------------------------------------
   real(8) :: T  ! time and distance to the cross point
   real(8) :: CS_total, MFP_total, Sampled_path  ! total cross section and MFP
   real(8) :: V, V1 ! [A/fs] velosity
   real(8) :: eps
   type(Target_atoms), pointer :: matter
   logical :: do_scattering
   
   eps = 1.0d-12    ! energy precision
   
   if (abs(Prtcl%Ekin) < eps) then ! electron with zero energy cannot fly
      Prtcl%ti = 1.0d26
      Prtcl%t_sc = 1.0d25
   else ! non-zero-energy electron can fly
      ! Check how far it is from the border of the simulation box boundary:
      call flight_time_to_boundary(used_target, numpar, Prtcl, T) ! module "MC_general_tools"
   
      ! The time of the next event (for boundary crossing):
      Prtcl%ti = Prtcl%t0 + T  ! [fs]

      do_scattering = .true.   ! by default, recalculate the scattering time
      if (present(no_scatternig) .or. ( abs(Prtcl%Ekin) < eps) ) then  ! electron with zero energy cannot scatter
         if (no_scatternig) do_scattering = .false.
      endif
   
      if (do_scattering) then     ! recalculate scattering event time:
         ! Now check the time to a scattering event:
         if (Prtcl%in_target > 0) then    ! it is not in vacuum
            ! Properties of the target material, inside of which the particle is:
            matter => used_target%Material(Prtcl%in_target)
            ! a) get the total cross sections from all possible scattering channels:
            CS_total = total_CS_from_chennels(Prtcl%Ekin, matter%El_inelastic_total%E, matter%El_inelastic_total%Total, &
                    E_array2=matter%El_elastic_total%E, CS_array2=matter%El_elastic_total%Total, &
                    E_array3=matter%El_Brems_total%E, CS_array3=matter%El_Brems_total%Total) ! module "CS_general_tools"
            ! b) get the total mean free path:
            MFP_total = MFP_from_sigma(CS_total, matter%At_Dens)  ! module "CS_general_tools"
            ! Sample free flight:
            Sampled_path = sample_Poisson(MFP_total)  ! module "Little_subroutines"  
            
            ! c) get the total flight time:
            V = sqrt( SUM(Prtcl%V(:)*Prtcl%V(:)) )  ! electron speed [A/fs]
            V1 = velosity_from_kinetic_energy(Prtcl%Ekin, g_me)    ! module "Relativity"
      
            ! Cross check that velosity is consistent:
            if ( abs(V-V1) > 1.0d-6*V1 ) then
               print*, 'Error in get_electron_flight_time:', V, velosity_from_kinetic_energy(Prtcl%Ekin, g_me) 
            endif
      
            !T =  Time_from_MFP(MFP_total, V) ! [fs] module "CS_general_tools"
            T =  Time_from_MFP(Sampled_path, V) ! [fs] module "CS_general_tools"
            ! The time of the next event (for scattering):
            Prtcl%t_sc = Prtcl%t0 + T  ! [fs]
            
!             write(*,'(f,f,f,f)') Prtcl%Ekin, Sampled_path, MFP_total
            
         else ! not inside of any target
            Prtcl%t_sc = 1.0d25    ! no scattering in vacuum
         endif
!          print*, 'get_electron_flight_time', Prtcl%ti, Prtcl%t_sc
      endif
   
      ! Check which event is shorter: boundary crossing or scattering:
      if (Prtcl%ti >= Prtcl%t_sc) Prtcl%ti = Prtcl%t_sc    ! the next event will be scattering
   endif ! (abs(Prtcl%Ekin) < eps))
   
   ! If user provided cut-off energy, check if electron crossed it:
   ! and save the disappearing electron energy to pass to MD:
   if ((Prtcl%Ekin <= numpar%El_Cutoff) .and. (numpar%DO_MD)) then
      call add_energy_into_MD_array(numpar, Prtcl%Ekin, Prtcl%R, MD_supce, E_e)   ! module "MC_general_tools"
   endif
   ! now remove the particle if needed:
   call check_cut_off(numpar%El_Cutoff, Prtcl)   ! below

   
   nullify(matter)
end subroutine get_electron_flight_time



subroutine get_positron_flight_time(used_target, numpar, Prtcl, no_scatternig) ! Find within which given particle is now
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Positron), intent(inout) :: Prtcl        ! positron as an object
   logical, optional :: no_scatternig       ! recalculate the scattering event, or leave the old time
   !-------------------------------------------
   real(8) :: T  ! time and distance to the cross point
   real(8) :: CS_total, MFP_total, Sampled_path  ! total cross section and MFP
   real(8) :: V ! [A/fs] velosity
   real(8) :: eps
   type(Target_atoms), pointer :: matter
   logical :: do_scattering
   
   eps = 1.0d-12    ! energy precision
   
   if (abs(Prtcl%Ekin) < eps) then ! positron with zero energy cannot fly
      Prtcl%ti = 1.0d26
      Prtcl%t_sc = 1.0d25
   else ! non-zero-energy positron can fly
   
      ! Check how far it is from the border of the simulation box boundary:
      call flight_time_to_boundary(used_target, numpar, Prtcl, T) ! module "MC_general_tools"
      ! The time of the next event (for boundary crossing):
      Prtcl%ti = Prtcl%t0 + T  ! [fs]
      
      do_scattering = .true.   ! by default, recalculate the scattering time
      if (present(no_scatternig) .or. ( abs(Prtcl%Ekin) < eps) ) then  ! electron with zero energy cannot scatter
         if (no_scatternig) do_scattering = .false.
      endif
      
      if (do_scattering) then     ! recalculate scattering event time:
         ! Now check the time to a scattering event:
         if (Prtcl%in_target > 0) then    ! it is not in vacuum
            ! Properties of the target material, inside of which the particle is:
            matter => used_target%Material(Prtcl%in_target)
            ! a) get the total cross sections from all possible scattering channels:
            CS_total = total_CS_from_chennels(Prtcl%Ekin, matter%Pos_inelastic_total%E, matter%Pos_inelastic_total%Total, &
                    E_array2=matter%Pos_elastic_total%E, CS_array2=matter%Pos_elastic_total%Total, &
                    E_array3=matter%Pos_Brems_total%E, CS_array3=matter%Pos_Brems_total%Total, &
                    E_array4=matter%Pos_annihil_total%E, CS_array4=matter%Pos_annihil_total%Total) ! module "CS_general_tools"
            ! b) get the total mean free path:
            MFP_total = MFP_from_sigma(CS_total, matter%At_Dens)  ! module "CS_general_tools"
            ! Sample free flight:
            Sampled_path = sample_Poisson(MFP_total)  ! module "Little_subroutines"  
            ! c) get the total flight time:
            V = sqrt( SUM(Prtcl%V(:)*Prtcl%V(:)) )   ! Positron speed is the speed of light [A/fs]
            !T =  Time_from_MFP(MFP_total, V) ! [fs] module "CS_general_tools"
            T =  Time_from_MFP(Sampled_path, V) ! [fs] module "CS_general_tools"
            ! The time of the next event (for scattering):
            Prtcl%t_sc = Prtcl%t0 + T  ! [fs]
         else ! not inside of any target
            Prtcl%t_sc = 1.0d25    ! no scattering in vacuum
         endif 
      endif ! do_scattering 

      ! Check which event is shorter: boundary crossing or scattering:
      if (Prtcl%ti >= Prtcl%t_sc) Prtcl%ti = Prtcl%t_sc    ! the next event will be scattering
   endif !  (abs(Prtcl%Ekin) < eps))
   
   ! If user provided cut-off energy, check if positron crossed it:
   call check_cut_off(numpar%Pos_Cutoff, Prtcl)   ! below
   
   nullify(matter)
end subroutine get_positron_flight_time



subroutine get_SHI_flight_time(used_target, numpar, Prtcl) ! Find within which given particle is now
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(SHI), intent(inout) :: Prtcl        ! photon as an object
   !-------------------------------------------
   real(8) :: T  ! time and distance to the cross point
   real(8) :: CS_total, MFP_total, Sampled_path  ! total cross section and MFP
   real(8) :: V ! [A/fs] velosity
   real(8) :: eps
   type(Target_atoms), pointer :: matter
   
   ! Prtcl%KOA - is undefined! Define it!!!
   eps = 1.0d-12    ! energy precision
   
   if (abs(Prtcl%Ekin) < eps) then ! ion with zero energy cannot fly
      Prtcl%ti = 1.0d26
      Prtcl%t_sc = 1.0d25
   else ! non-zero-energy ion can fly
   
      ! Check how far it is from the border of the simulation box boundary:
      call flight_time_to_boundary(used_target, numpar, Prtcl, T) ! module "MC_general_tools"
      ! The time of the next event (for boundary crossing):
      Prtcl%ti = Prtcl%t0 + T  ! [fs]
      ! Now check the time to a scattering event:
      if (Prtcl%in_target > 0) then    ! it is not in vacuum
         ! Properties of the target material, inside of which the particle is:
         matter => used_target%Material(Prtcl%in_target)
         ! a) get the total cross sections from all possible scattering channels:
         CS_total = total_CS_from_chennels(Prtcl%Ekin, matter%SHI_inelastic_total(Prtcl%KOA)%E, matter%SHI_inelastic_total(Prtcl%KOA)%Total) ! module "CS_general_tools"
         ! b) get the total mean free path:
         MFP_total = MFP_from_sigma(CS_total, matter%At_Dens)  ! module "CS_general_tools"
         ! Sample free flight:
         Sampled_path = sample_Poisson(MFP_total)  ! module "Little_subroutines"
         ! c) get the total flight time:
         V = sqrt( SUM(Prtcl%V(:)*Prtcl%V(:)) )   ! electron speed is the speed of light [A/fs]
         !T =  Time_from_MFP(MFP_total, V) ! [fs] module "CS_general_tools"
         T =  Time_from_MFP(Sampled_path, V) ! [fs] module "CS_general_tools"
         ! The time of the next event (for scattering):
         Prtcl%t_sc = Prtcl%t0 + T  ! [fs]
      else ! not inside of any target
         Prtcl%t_sc = 1.0d25    ! no scattering in vacuum
      endif

      ! Check which event is shorter: boundary crossing or scattering:
      if (Prtcl%ti >= Prtcl%t_sc) Prtcl%ti = Prtcl%t_sc    ! the next event will be scattering
   endif
   
   ! If user provided cut-off energy, check if SHI crossed it:
   call check_cut_off(numpar%SHI_Cutoff, Prtcl)   ! below
  
   nullify(matter)
end subroutine get_SHI_flight_time

   

subroutine get_hole_flight_time(used_target, numpar, Prtcl, MD_supce, E_h, no_scatternig) ! Find within which given particle is now
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(Hole), intent(inout) :: Prtcl        ! hole as an object
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters for connection between MC and MD modules
   real(8), dimension(:,:,:), intent(inout) :: E_h ! data to pass to MD later
   logical, optional :: no_scatternig       ! recalculate the scattering event, or leave the old time
   !-------------------------------------------
   real(8) :: T  ! time and distance to the cross point
   real(8) :: CS_total, MFP_total, Sampled_path  ! total cross section and MFP
   real(8) :: V ! [A/fs] velosity
   real(8) :: eps
   type(Target_atoms), pointer :: matter
   logical :: valent, do_scattering
   
   eps = 1.0d-12    ! energy precision
   
   ! Is it a valence hole, or a core hole?
   if ((Prtcl%KOA == 0) .and. (Prtcl%Sh == 0) ) then    ! it's valent band
      valent = .true.
   else ! check if it is valent or not
      valent = used_target%Material(Prtcl%in_target)%Elements(Prtcl%KOA)%valent(Prtcl%Sh)
   endif
   
   VAL:if (valent) then ! such a hole is mobile:

      if (abs(Prtcl%Ekin) < eps) then ! hole with zero energy cannot fly
         Prtcl%ti = 1.0d26
         Prtcl%t_sc = 1.0d25
      else ! non-zero-energy hole can fly
         ! Check how far it is from the border of the simulation box boundary:
         call flight_time_to_boundary(used_target, numpar, Prtcl, T) ! module "MC_general_tools"
         ! The time of the next event (for boundary crossing):
         Prtcl%ti = Prtcl%t0 + T  ! [fs]
         
         do_scattering = .true.   ! by default, recalculate the scattering time
         if (present(no_scatternig) .or. ( abs(Prtcl%Ekin) < eps) ) then  ! electron with zero energy cannot scatter
            if (no_scatternig) do_scattering = .false.
         endif
         
         if (do_scattering) then     ! recalculate scattering event time:
            ! Now check the time to a scattering event:
            if (Prtcl%in_target > 0) then    ! it is not in vacuum
               ! Properties of the target material, inside of which the particle is:
               matter => used_target%Material(Prtcl%in_target)
               ! a) get the total cross sections from all possible scattering channels:
               CS_total = total_CS_from_chennels(Prtcl%Ekin, matter%H_inelastic_total%E, matter%H_inelastic_total%Total, &
                    E_array2=matter%H_elastic_total%E, CS_array2=matter%H_elastic_total%Total) ! module "CS_general_tools"
               ! b) get the total mean free path:
               MFP_total = MFP_from_sigma(CS_total, matter%At_Dens)  ! module "CS_general_tools"
               ! Sample free flight:
               Sampled_path = sample_Poisson(MFP_total)  ! module "Little_subroutines"
               ! c) get the total flight time:
               V = sqrt( SUM(Prtcl%V(:)*Prtcl%V(:)) )   ! hole speed is the speed of light [A/fs]
               !T =  Time_from_MFP(MFP_total, V) ! [fs] module "CS_general_tools"
               T =  Time_from_MFP(Sampled_path, V) ! [fs] module "CS_general_tools"
               ! The time of the next event (for scattering):
               Prtcl%t_sc = Prtcl%t0 + T  ! [fs]
            else ! not inside of any target
               Prtcl%t_sc = 1.0d25    ! no scattering in vacuum
            endif ! (Prtcl%in_target > 0)
         endif ! (do_scattering)
         ! Check which event is shorter: boundary crossing or scattering:
         if (Prtcl%ti >= Prtcl%t_sc) Prtcl%ti = Prtcl%t_sc    ! the next event will be scattering
      endif
      
      ! If user provided cut-off energy, check if hole crossed it:
      ! and save the disappearing electron energy to pass to MD:
      if ((Prtcl%Ekin <= numpar%H_Cutoff) .and. (numpar%DO_MD)) then
         if (Prtcl%in_target > 0) then
            ! For a hole, add its potential energy (equal to band gap):
            matter => used_target%Material(Prtcl%in_target)
            call add_energy_into_MD_array(numpar, Prtcl%Ekin + matter%DOS%Egap, Prtcl%R, MD_supce, E_h)   ! module "MC_general_tools"
         else
            print*, 'Trouble in get_hole_flight_time, hole is in vacuum:', Prtcl%in_target, Prtcl%Ekin, Prtcl%R
         endif
      endif
      ! now remove the particle if needed:
      call check_cut_off(numpar%H_Cutoff, Prtcl)   ! below

   else VAL ! core hole
      ! Those decay, but do not move:
      call core_hole_decay_time(used_target%Material(Prtcl%in_target)%Elements(Prtcl%KOA), Prtcl%Sh, T)  ! below
      Prtcl%t_sc = Prtcl%t0 + T     ! count from the current time
      Prtcl%ti = Prtcl%t_sc         ! the next event will be decay      
!       write(*,'(a,f,f,f)') 'get_hole_flight_time', T, Prtcl%t0, Prtcl%t_sc
   endif VAL
   
   nullify(matter)
end subroutine get_hole_flight_time



subroutine core_hole_decay_time(Element, SHL,  t_dec)
   type(Atom_kind), intent(in) :: Element   ! which atoms this target is constructed of
   integer, intent(in) :: SHL   ! index of the shell
   real(8), intent(out) :: t_dec    ! [fs] sampled decay time
   real(8) :: t_tot

   ! 1) Get the total mean decay time from the Auger and the radiative decay times:
   t_tot = 1.0d0/Element%Auger(SHL) + 1.0d0/Element%Radiat(SHL) ! [1/fs] inverse total decay time
   t_tot = 1.0d0/t_tot  ! [fs] total decay time
   
   ! 2) Sample the realized decay time:
   t_dec = sample_Poisson(t_tot)  ! module "Little_subroutines"
end subroutine core_hole_decay_time



subroutine select_process_by_time(ind, time1, time2, time3, time4, time5)
   integer, intent(out) :: ind  ! which process is sampled
   real(8), intent(in) :: time1 ! time of the first process; all times must have the same units, e.g. [fs]
   real(8), intent(in), optional :: time2, time3, time4, time5  ! times of other possible processes
   real(8) :: RN, tot_tim, Inv_sum(5)
   integer :: i
   logical :: more_than_one
   more_than_one = present(time2) .or. present(time3) .or. present(time4) .or. present(time5)
   
   if (.not.more_than_one) then ! there is only one time, nothing to chose from
      ind = 1
   else ! there are more processes to select from:
      tot_tim = 0.0d0   ! just to start
      Inv_sum = 0.0d0   ! just to start
      
      ! Matthiensen's rule for summing up characteristic times (via inversed times):
      tot_tim = tot_tim + 1.0d0/time1
      Inv_sum(1) = tot_tim 
      
      if (present(time2)) then
         tot_tim = tot_tim + 1.0d0/time2
         Inv_sum(2) = tot_tim 
      endif
      
      if (present(time3)) then
         tot_tim = tot_tim + 1.0d0/time3
         Inv_sum(3) = tot_tim 
      endif
      
      if (present(time4)) then
         tot_tim = tot_tim + 1.0d0/time4
         Inv_sum(4) = tot_tim 
      endif
      
      if (present(time5)) then
         tot_tim = tot_tim + 1.0d0/time5
         Inv_sum(5) = tot_tim 
      endif

      ! Normalize to get probability:
      Inv_sum(:) = Inv_sum(:)/tot_tim
      
      ! Get the random number, to sample partial cross section:
      call random_number(RN)   ! intrinsic FORTRAN subroutine
      
      ! Get the index of the process from the sampled time:
      i = 1
      do while (RN > Inv_sum(i))
!          if (i == size(Inv_sum)) exit   ! the last available
         i = i + 1
      enddo
      ind = i   ! the sampled process index
      
   endif

end subroutine select_process_by_time



subroutine check_cut_off(Ecutoff, Prtcl)
   real(8), intent(in) :: Ecutoff   ! [eV] cut off energy
   class(Particle), intent(inout) :: Prtcl    ! undefined particle as an object
!    type(Electron), dimension(:), intent(inout) :: Prtcl    ! undefined particle as an object
   if (Prtcl%Ekin <= Ecutoff) then   ! stop the particle, but don't remove from simulation
      Prtcl%ti = 1.0d24
      Prtcl%t_sc = 1.0d25
      Prtcl%V(:) = 0.0d0
      Prtcl%SV(:) = 0.0d0
      Prtcl%S0(:) = 0.0d0
      Prtcl%V0(:) = 0.0d0
      Prtcl%SV0(:) = 0.0d0
   endif
end subroutine check_cut_off



!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! Dealing with MC arrays:


subroutine renew_atomic_arrays(MC)
   type(MC_arrays), dimension(:), intent(inout) :: MC   ! all MC arrays for all particles; size = number of iterations
   integer i, Nsiz
   Nsiz = size(MC)
   do i = 1, Nsiz
      MC(i)%N_at_nrg = 0  ! restart counting the number of active atomic events of energy transfer
      MC(i)%MC_Atoms_events(:)%active = .false.  ! deactivate all "particles"
      MC(i)%MC_Atoms_events(:)%Ekin = 0.0d0      ! deactivate all "particles"
   enddo
end subroutine renew_atomic_arrays



subroutine remove_particle_from_MC_array(N_particles, i_remove, Prtcl)
   integer, intent(inout) :: N_particles    ! that's how many active particles there were
   integer, intent(in) :: i_remove  ! index of particle to be removed
   class(Particle), dimension(:), intent(inout) :: Prtcl	! undefined particle as an object
   !--------------------------------------------
   integer :: i

!    print*, 'Before:', Prtcl(i_remove)%Ekin, Prtcl(i_remove+1)%Ekin, N_particles
!    print*,  Prtcl(i_remove)%R(:)
!    do i = 1, N_particles
!       if (abs(Prtcl(i)%Ekin-50.0d0) < 1.0d-6) then
!          print*, 'Be', i, Prtcl(i)%Ekin, Prtcl(i)%generation
!       endif
!    enddo

   ! Now one particle is removed:
   do i = i_remove, N_particles-1
      Prtcl(i)%active = Prtcl(i+1)%active
      Prtcl(i)%generation = Prtcl(i+1)%generation
      Prtcl(i)%in_target = Prtcl(i+1)%in_target
      Prtcl(i)%Ekin = Prtcl(i+1)%Ekin 
      Prtcl(i)%t0 = Prtcl(i+1)%t0
      Prtcl(i)%ti = Prtcl(i+1)%ti
      Prtcl(i)%t_sc = Prtcl(i+1)%t_sc
      Prtcl(i)%R(:) = Prtcl(i+1)%R(:) 
      Prtcl(i)%S(:) = Prtcl(i+1)%S(:) 
      Prtcl(i)%V(:) = Prtcl(i+1)%V(:) 
      Prtcl(i)%SV(:) = Prtcl(i+1)%SV(:) 
      Prtcl(i)%R0(:) = Prtcl(i+1)%R0(:) 
      Prtcl(i)%S0(:) = Prtcl(i+1)%S0(:) 
      Prtcl(i)%V0(:) = Prtcl(i+1)%V0(:) 
      Prtcl(i)%SV0(:) = Prtcl(i+1)%SV0(:)
      Prtcl(i)%Mass = Prtcl(i+1)%Mass
      select type (Prtcl)	! which particle array that should become
         type is (Photon)
         type is (Electron)
            Prtcl(i)%A(:) = Prtcl(i+1)%A(:) 
            Prtcl(i)%Force(:) = Prtcl(i+1)%Force(:) 
         type is (Positron)
            Prtcl(i)%A(:) = Prtcl(i+1)%A(:) 
            Prtcl(i)%Force(:) = Prtcl(i+1)%Force(:) 
         type is (Hole)
            Prtcl(i)%A(:) = Prtcl(i+1)%A(:) 
            Prtcl(i)%Force(:) = Prtcl(i+1)%Force(:) 
            Prtcl(i)%KOA = Prtcl(i+1)%KOA 
            Prtcl(i)%Sh = Prtcl(i+1)%Sh 
            Prtcl(i)%valent = Prtcl(i+1)%valent
         type is (Atom)
            Prtcl(i)%Z = Prtcl(i+1)%Z 
            Prtcl(i)%Name = Prtcl(i+1)%Name 
            Prtcl(i)%A(:) = Prtcl(i+1)%A(:)
            Prtcl(i)%Force(:) = Prtcl(i+1)%Force(:)
         type is (SHI)
            Prtcl(i)%Z = Prtcl(i+1)%Z
            Prtcl(i)%Name = Prtcl(i+1)%Name
            Prtcl(i)%Zeff = Prtcl(i+1)%Zeff
            Prtcl(i)%Meff = Prtcl(i+1)%Meff
            Prtcl(i)%A(:) = Prtcl(i+1)%A(:)
            Prtcl(i)%Force(:) = Prtcl(i+1)%Force(:) 
      end select
   enddo
   ! Exclude the last particle:
   call set_default_particle(Prtcl(N_particles))    ! above
   ! And remove it from the counter:
   N_particles = N_particles - 1
   
!    print*, 'After:', Prtcl(i_remove)%Ekin, Prtcl(i_remove+1)%Ekin, N_particles
!    do i = 1, N_particles
!       if (abs(Prtcl(i)%Ekin-50.0d0) < 1.0d-6) then
!          print*, 'Af', i, Prtcl(i)%Ekin, Prtcl(i)%generation
!       endif
!    enddo
end subroutine remove_particle_from_MC_array



pure subroutine extend_MC_array_Electrons(Prtcl)
   type(Electron), dimension(:), allocatable, intent(inout) :: Prtcl    ! all electrons as objects 
   type(Electron), dimension(:), allocatable :: Prtcl_temp     ! temporary aray of electrons
   integer :: siz, siz2, j
   siz = size(Prtcl)
   ! Allocate the temporary array to transiently store data:
   allocate(Prtcl_temp(siz))
   Prtcl_temp(:)%active = Prtcl(:)%active
   Prtcl_temp(:)%Ekin = Prtcl(:)%Ekin 
   Prtcl_temp(:)%Mass = Prtcl(:)%Mass
   Prtcl_temp(:)%t0 = Prtcl(:)%t0
   Prtcl_temp(:)%ti = Prtcl(:)%ti
   Prtcl_temp(:)%t_sc = Prtcl(:)%t_sc
   Prtcl_temp(:)%generation = Prtcl(:)%generation
   Prtcl_temp(:)%in_target = Prtcl(:)%in_target
   Prtcl_temp(:)%R(1) = Prtcl(:)%R(1) 
   Prtcl_temp(:)%R(2) = Prtcl(:)%R(2) 
   Prtcl_temp(:)%R(3) = Prtcl(:)%R(3) 
   Prtcl_temp(:)%S(1) = Prtcl(:)%S(1) 
   Prtcl_temp(:)%S(2) = Prtcl(:)%S(2) 
   Prtcl_temp(:)%S(3) = Prtcl(:)%S(3) 
   Prtcl_temp(:)%V(1) = Prtcl(:)%V(1) 
   Prtcl_temp(:)%V(2) = Prtcl(:)%V(2) 
   Prtcl_temp(:)%V(3) = Prtcl(:)%V(3) 
   Prtcl_temp(:)%SV(1) = Prtcl(:)%SV(1) 
   Prtcl_temp(:)%SV(2) = Prtcl(:)%SV(2) 
   Prtcl_temp(:)%SV(3) = Prtcl(:)%SV(3) 
   Prtcl_temp(:)%R0(1) = Prtcl(:)%R0(1) 
   Prtcl_temp(:)%R0(2) = Prtcl(:)%R0(2)
   Prtcl_temp(:)%R0(3) = Prtcl(:)%R0(3) 
   Prtcl_temp(:)%S0(1) = Prtcl(:)%S0(1) 
   Prtcl_temp(:)%S0(2) = Prtcl(:)%S0(2) 
   Prtcl_temp(:)%S0(3) = Prtcl(:)%S0(3) 
   Prtcl_temp(:)%V0(1) = Prtcl(:)%V0(1) 
   Prtcl_temp(:)%V0(2) = Prtcl(:)%V0(2) 
   Prtcl_temp(:)%V0(3) = Prtcl(:)%V0(3) 
   Prtcl_temp(:)%SV0(1) =Prtcl(:)%SV0(1)
   Prtcl_temp(:)%SV0(2) =Prtcl(:)%SV0(2)
   Prtcl_temp(:)%SV0(3) =Prtcl(:)%SV0(3)

   call copy_MC_array(Prtcl_temp, Prtcl) ! below

   siz2 = 2*siz
   deallocate(Prtcl)
   allocate(Prtcl(siz2))
   
   ! Copy the data back into the arrays:
   Prtcl(1:siz)%active = Prtcl_temp(1:siz)%active
   Prtcl(1:siz)%Ekin = Prtcl_temp(1:siz)%Ekin
   Prtcl(1:siz)%Mass = Prtcl_temp(1:siz)%Mass
   Prtcl(1:siz)%t0 = Prtcl_temp(1:siz)%t0
   Prtcl(1:siz)%ti = Prtcl_temp(1:siz)%ti
   Prtcl(1:siz)%t_sc = Prtcl_temp(1:siz)%t_sc
   Prtcl(1:siz)%generation = Prtcl_temp(1:siz)%generation
   Prtcl(1:siz)%in_target = Prtcl_temp(1:siz)%in_target
   Prtcl(1:siz)%R(1) = Prtcl_temp(1:siz)%R(1)
   Prtcl(1:siz)%R(2) = Prtcl_temp(1:siz)%R(2)
   Prtcl(1:siz)%R(3) = Prtcl_temp(1:siz)%R(3)
   Prtcl(1:siz)%S(1) = Prtcl_temp(1:siz)%S(1) 
   Prtcl(1:siz)%S(2) = Prtcl_temp(1:siz)%S(2) 
   Prtcl(1:siz)%S(3) = Prtcl_temp(1:siz)%S(3) 
   Prtcl(1:siz)%V(1) = Prtcl_temp(1:siz)%V(1) 
   Prtcl(1:siz)%V(2) = Prtcl_temp(1:siz)%V(2) 
   Prtcl(1:siz)%V(3) = Prtcl_temp(1:siz)%V(3) 
   Prtcl(1:siz)%SV(1) = Prtcl_temp(1:siz)%SV(1) 
   Prtcl(1:siz)%SV(2) = Prtcl_temp(1:siz)%SV(2) 
   Prtcl(1:siz)%SV(3) = Prtcl_temp(1:siz)%SV(3) 
   Prtcl(1:siz)%R0(1) = Prtcl_temp(1:siz)%R0(1) 
   Prtcl(1:siz)%R0(2) = Prtcl_temp(1:siz)%R0(2) 
   Prtcl(1:siz)%R0(3) = Prtcl_temp(1:siz)%R0(3) 
   Prtcl(1:siz)%S0(1) = Prtcl_temp(1:siz)%S0(1) 
   Prtcl(1:siz)%S0(2) = Prtcl_temp(1:siz)%S0(2) 
   Prtcl(1:siz)%S0(3) = Prtcl_temp(1:siz)%S0(3) 
   Prtcl(1:siz)%V0(1) = Prtcl_temp(1:siz)%V0(1) 
   Prtcl(1:siz)%V0(2) = Prtcl_temp(1:siz)%V0(2) 
   Prtcl(1:siz)%V0(3) = Prtcl_temp(1:siz)%V0(3) 
   Prtcl(1:siz)%SV0(1) =Prtcl_temp(1:siz)%SV0(1)
   Prtcl(1:siz)%SV0(2) =Prtcl_temp(1:siz)%SV0(2)
   Prtcl(1:siz)%SV0(3) =Prtcl_temp(1:siz)%SV0(3)

   call copy_MC_array(Prtcl(1:siz), Prtcl_temp(1:siz)) ! below
 
 ! All beyond are default to start with:
   do j = siz+1, siz2
      call set_default_particle(Prtcl(j))  ! above
   enddo
   ! clean up:
   deallocate(Prtcl_temp)
end subroutine extend_MC_array_Electrons



pure subroutine extend_MC_array_Photons(Prtcl)
   type(Photon), dimension(:), allocatable, intent(inout) :: Prtcl    ! all electrons as objects 
   type(Photon), dimension(:), allocatable :: Prtcl_temp     ! temporary aray of electrons
   integer :: siz, siz2, j
   siz = size(Prtcl)
   ! Allocate the temporary array to transiently store data:
   allocate(Prtcl_temp(siz))
   
   Prtcl_temp(:)%active = Prtcl(:)%active
   Prtcl_temp(:)%generation = Prtcl(:)%generation
   Prtcl_temp(:)%in_target = Prtcl(:)%in_target
   Prtcl_temp(:)%Ekin = Prtcl(:)%Ekin 
   Prtcl_temp(:)%t0 = Prtcl(:)%t0
   Prtcl_temp(:)%ti = Prtcl(:)%ti
   Prtcl_temp(:)%t_sc = Prtcl(:)%t_sc
   Prtcl_temp(:)%R(1) = Prtcl(:)%R(1) 
   Prtcl_temp(:)%R(2) = Prtcl(:)%R(2) 
   Prtcl_temp(:)%R(3) = Prtcl(:)%R(3) 
   Prtcl_temp(:)%S(1) = Prtcl(:)%S(1) 
   Prtcl_temp(:)%S(2) = Prtcl(:)%S(2) 
   Prtcl_temp(:)%S(3) = Prtcl(:)%S(3) 
   Prtcl_temp(:)%V(1) = Prtcl(:)%V(1) 
   Prtcl_temp(:)%V(2) = Prtcl(:)%V(2) 
   Prtcl_temp(:)%V(3) = Prtcl(:)%V(3) 
   Prtcl_temp(:)%SV(1) = Prtcl(:)%SV(1) 
   Prtcl_temp(:)%SV(2) = Prtcl(:)%SV(2) 
   Prtcl_temp(:)%SV(3) = Prtcl(:)%SV(3) 
   Prtcl_temp(:)%R0(1) = Prtcl(:)%R0(1) 
   Prtcl_temp(:)%R0(2) = Prtcl(:)%R0(2)
   Prtcl_temp(:)%R0(3) = Prtcl(:)%R0(3) 
   Prtcl_temp(:)%S0(1) = Prtcl(:)%S0(1) 
   Prtcl_temp(:)%S0(2) = Prtcl(:)%S0(2) 
   Prtcl_temp(:)%S0(3) = Prtcl(:)%S0(3) 
   Prtcl_temp(:)%V0(1) = Prtcl(:)%V0(1) 
   Prtcl_temp(:)%V0(2) = Prtcl(:)%V0(2) 
   Prtcl_temp(:)%V0(3) = Prtcl(:)%V0(3) 
   Prtcl_temp(:)%SV0(1) =Prtcl(:)%SV0(1)
   Prtcl_temp(:)%SV0(2) =Prtcl(:)%SV0(2)
   Prtcl_temp(:)%SV0(3) =Prtcl(:)%SV0(3)
   Prtcl_temp(:)%Mass = Prtcl(:)%Mass

!    call copy_MC_array(Prtcl_temp, Prtcl) ! below

   siz2 = 2*siz
   deallocate(Prtcl)
   allocate(Prtcl(siz2))
   
   ! Copy the data back into the arrays:
   Prtcl(1:siz)%active = Prtcl_temp(1:siz)%active
   Prtcl(1:siz)%generation = Prtcl_temp(1:siz)%generation
   Prtcl(1:siz)%in_target = Prtcl_temp(1:siz)%in_target
   Prtcl(1:siz)%Ekin = Prtcl_temp(1:siz)%Ekin 
   Prtcl(1:siz)%t0 = Prtcl_temp(1:siz)%t0
   Prtcl(1:siz)%ti = Prtcl_temp(1:siz)%ti
   Prtcl(1:siz)%t_sc = Prtcl_temp(1:siz)%t_sc
   Prtcl(1:siz)%R(1) = Prtcl_temp(1:siz)%R(1)
   Prtcl(1:siz)%R(2) = Prtcl_temp(1:siz)%R(2)
   Prtcl(1:siz)%R(3) = Prtcl_temp(1:siz)%R(3)
   Prtcl(1:siz)%S(1) = Prtcl_temp(1:siz)%S(1) 
   Prtcl(1:siz)%S(2) = Prtcl_temp(1:siz)%S(2) 
   Prtcl(1:siz)%S(3) = Prtcl_temp(1:siz)%S(3) 
   Prtcl(1:siz)%V(1) = Prtcl_temp(1:siz)%V(1) 
   Prtcl(1:siz)%V(2) = Prtcl_temp(1:siz)%V(2) 
   Prtcl(1:siz)%V(3) = Prtcl_temp(1:siz)%V(3) 
   Prtcl(1:siz)%SV(1) = Prtcl_temp(1:siz)%SV(1) 
   Prtcl(1:siz)%SV(2) = Prtcl_temp(1:siz)%SV(2) 
   Prtcl(1:siz)%SV(3) = Prtcl_temp(1:siz)%SV(3) 
   Prtcl(1:siz)%R0(1) = Prtcl_temp(1:siz)%R0(1) 
   Prtcl(1:siz)%R0(2) = Prtcl_temp(1:siz)%R0(2) 
   Prtcl(1:siz)%R0(3) = Prtcl_temp(1:siz)%R0(3) 
   Prtcl(1:siz)%S0(1) = Prtcl_temp(1:siz)%S0(1) 
   Prtcl(1:siz)%S0(2) = Prtcl_temp(1:siz)%S0(2) 
   Prtcl(1:siz)%S0(3) = Prtcl_temp(1:siz)%S0(3) 
   Prtcl(1:siz)%V0(1) = Prtcl_temp(1:siz)%V0(1) 
   Prtcl(1:siz)%V0(2) = Prtcl_temp(1:siz)%V0(2) 
   Prtcl(1:siz)%V0(3) = Prtcl_temp(1:siz)%V0(3) 
   Prtcl(1:siz)%SV0(1) =Prtcl_temp(1:siz)%SV0(1)
   Prtcl(1:siz)%SV0(2) =Prtcl_temp(1:siz)%SV0(2)
   Prtcl(1:siz)%SV0(3) =Prtcl_temp(1:siz)%SV0(3)
   Prtcl(1:siz)%Mass = Prtcl_temp(1:siz)%Mass

!    call copy_MC_array(Prtcl(1:siz), Prtcl_temp(1:siz)) ! below
 
 ! All beyond are default to start with:
   do j = siz+1, siz2
      call set_default_particle(Prtcl(j))  ! above
   enddo
   ! clean up:
   deallocate(Prtcl_temp)
end subroutine extend_MC_array_Photons



pure subroutine extend_MC_array_Holes(Prtcl)
   type(Hole), dimension(:), allocatable, intent(inout) :: Prtcl    ! all electrons as objects 
   type(Hole), dimension(:), allocatable :: Prtcl_temp     ! temporary aray of electrons
   integer :: siz, siz2, j
   siz = size(Prtcl)
   ! Allocate the temporary array to transiently store data:
   allocate(Prtcl_temp(siz))
   
   Prtcl_temp(:)%active = Prtcl(:)%active
   Prtcl_temp(:)%valent = Prtcl(:)%valent
   Prtcl_temp(:)%generation = Prtcl(:)%generation
   Prtcl_temp(:)%in_target = Prtcl(:)%in_target
   Prtcl_temp(:)%Ekin = Prtcl(:)%Ekin 
   Prtcl_temp(:)%t0 = Prtcl(:)%t0
   Prtcl_temp(:)%ti = Prtcl(:)%ti
   Prtcl_temp(:)%t_sc = Prtcl(:)%t_sc
   Prtcl_temp(:)%R(1) = Prtcl(:)%R(1) 
   Prtcl_temp(:)%R(2) = Prtcl(:)%R(2) 
   Prtcl_temp(:)%R(3) = Prtcl(:)%R(3) 
   Prtcl_temp(:)%S(1) = Prtcl(:)%S(1) 
   Prtcl_temp(:)%S(2) = Prtcl(:)%S(2) 
   Prtcl_temp(:)%S(3) = Prtcl(:)%S(3) 
   Prtcl_temp(:)%V(1) = Prtcl(:)%V(1) 
   Prtcl_temp(:)%V(2) = Prtcl(:)%V(2) 
   Prtcl_temp(:)%V(3) = Prtcl(:)%V(3) 
   Prtcl_temp(:)%SV(1) = Prtcl(:)%SV(1) 
   Prtcl_temp(:)%SV(2) = Prtcl(:)%SV(2) 
   Prtcl_temp(:)%SV(3) = Prtcl(:)%SV(3) 
   Prtcl_temp(:)%R0(1) = Prtcl(:)%R0(1) 
   Prtcl_temp(:)%R0(2) = Prtcl(:)%R0(2)
   Prtcl_temp(:)%R0(3) = Prtcl(:)%R0(3) 
   Prtcl_temp(:)%S0(1) = Prtcl(:)%S0(1) 
   Prtcl_temp(:)%S0(2) = Prtcl(:)%S0(2) 
   Prtcl_temp(:)%S0(3) = Prtcl(:)%S0(3) 
   Prtcl_temp(:)%V0(1) = Prtcl(:)%V0(1) 
   Prtcl_temp(:)%V0(2) = Prtcl(:)%V0(2) 
   Prtcl_temp(:)%V0(3) = Prtcl(:)%V0(3) 
   Prtcl_temp(:)%SV0(1) =Prtcl(:)%SV0(1)
   Prtcl_temp(:)%SV0(2) =Prtcl(:)%SV0(2)
   Prtcl_temp(:)%SV0(3) =Prtcl(:)%SV0(3)
   Prtcl_temp(:)%Mass = Prtcl(:)%Mass

   call copy_MC_array(Prtcl_temp, Prtcl) ! below

   siz2 = 2*siz
   deallocate(Prtcl)
   allocate(Prtcl(siz2))
   
   ! Copy the data back into the arrays:
   Prtcl(1:siz)%active = Prtcl_temp(1:siz)%active
   Prtcl(1:siz)%valent = Prtcl_temp(1:siz)%valent
   Prtcl(1:siz)%generation = Prtcl_temp(1:siz)%generation
   Prtcl(1:siz)%in_target = Prtcl_temp(1:siz)%in_target
   Prtcl(1:siz)%Ekin = Prtcl_temp(1:siz)%Ekin 
   Prtcl(1:siz)%t0 = Prtcl_temp(1:siz)%t0
   Prtcl(1:siz)%ti = Prtcl_temp(1:siz)%ti
   Prtcl(1:siz)%t_sc = Prtcl_temp(1:siz)%t_sc
   Prtcl(1:siz)%R(1) = Prtcl_temp(1:siz)%R(1)
   Prtcl(1:siz)%R(2) = Prtcl_temp(1:siz)%R(2)
   Prtcl(1:siz)%R(3) = Prtcl_temp(1:siz)%R(3)
   Prtcl(1:siz)%S(1) = Prtcl_temp(1:siz)%S(1) 
   Prtcl(1:siz)%S(2) = Prtcl_temp(1:siz)%S(2) 
   Prtcl(1:siz)%S(3) = Prtcl_temp(1:siz)%S(3) 
   Prtcl(1:siz)%V(1) = Prtcl_temp(1:siz)%V(1) 
   Prtcl(1:siz)%V(2) = Prtcl_temp(1:siz)%V(2) 
   Prtcl(1:siz)%V(3) = Prtcl_temp(1:siz)%V(3) 
   Prtcl(1:siz)%SV(1) = Prtcl_temp(1:siz)%SV(1) 
   Prtcl(1:siz)%SV(2) = Prtcl_temp(1:siz)%SV(2) 
   Prtcl(1:siz)%SV(3) = Prtcl_temp(1:siz)%SV(3) 
   Prtcl(1:siz)%R0(1) = Prtcl_temp(1:siz)%R0(1) 
   Prtcl(1:siz)%R0(2) = Prtcl_temp(1:siz)%R0(2) 
   Prtcl(1:siz)%R0(3) = Prtcl_temp(1:siz)%R0(3) 
   Prtcl(1:siz)%S0(1) = Prtcl_temp(1:siz)%S0(1) 
   Prtcl(1:siz)%S0(2) = Prtcl_temp(1:siz)%S0(2) 
   Prtcl(1:siz)%S0(3) = Prtcl_temp(1:siz)%S0(3) 
   Prtcl(1:siz)%V0(1) = Prtcl_temp(1:siz)%V0(1) 
   Prtcl(1:siz)%V0(2) = Prtcl_temp(1:siz)%V0(2) 
   Prtcl(1:siz)%V0(3) = Prtcl_temp(1:siz)%V0(3) 
   Prtcl(1:siz)%SV0(1) =Prtcl_temp(1:siz)%SV0(1)
   Prtcl(1:siz)%SV0(2) =Prtcl_temp(1:siz)%SV0(2)
   Prtcl(1:siz)%SV0(3) =Prtcl_temp(1:siz)%SV0(3)
   Prtcl(1:siz)%Mass = Prtcl_temp(1:siz)%Mass

   call copy_MC_array(Prtcl(1:siz), Prtcl_temp(1:siz)) ! below
 
 ! All beyond are default to start with:
   do j = siz+1, siz2
      call set_default_particle(Prtcl(j))  ! above
   enddo
   ! clean up:
   deallocate(Prtcl_temp)
end subroutine extend_MC_array_Holes



pure subroutine extend_MC_array_Positrons(Prtcl)
   type(Positron), dimension(:), allocatable, intent(inout) :: Prtcl    ! all electrons as objects 
   type(Positron), dimension(:), allocatable :: Prtcl_temp     ! temporary aray of electrons
   integer :: siz, siz2, j
   siz = size(Prtcl)
   ! Allocate the temporary array to transiently store data:
   allocate(Prtcl_temp(siz))
   
   Prtcl_temp(:)%active = Prtcl(:)%active
   Prtcl_temp(:)%generation = Prtcl(:)%generation
   Prtcl_temp(:)%in_target = Prtcl(:)%in_target
   Prtcl_temp(:)%Ekin = Prtcl(:)%Ekin 
   Prtcl_temp(:)%t0 = Prtcl(:)%t0
   Prtcl_temp(:)%ti = Prtcl(:)%ti
   Prtcl_temp(:)%t_sc = Prtcl(:)%t_sc
   Prtcl_temp(:)%R(1) = Prtcl(:)%R(1) 
   Prtcl_temp(:)%R(2) = Prtcl(:)%R(2) 
   Prtcl_temp(:)%R(3) = Prtcl(:)%R(3) 
   Prtcl_temp(:)%S(1) = Prtcl(:)%S(1) 
   Prtcl_temp(:)%S(2) = Prtcl(:)%S(2) 
   Prtcl_temp(:)%S(3) = Prtcl(:)%S(3) 
   Prtcl_temp(:)%V(1) = Prtcl(:)%V(1) 
   Prtcl_temp(:)%V(2) = Prtcl(:)%V(2) 
   Prtcl_temp(:)%V(3) = Prtcl(:)%V(3) 
   Prtcl_temp(:)%SV(1) = Prtcl(:)%SV(1) 
   Prtcl_temp(:)%SV(2) = Prtcl(:)%SV(2) 
   Prtcl_temp(:)%SV(3) = Prtcl(:)%SV(3) 
   Prtcl_temp(:)%R0(1) = Prtcl(:)%R0(1) 
   Prtcl_temp(:)%R0(2) = Prtcl(:)%R0(2)
   Prtcl_temp(:)%R0(3) = Prtcl(:)%R0(3) 
   Prtcl_temp(:)%S0(1) = Prtcl(:)%S0(1) 
   Prtcl_temp(:)%S0(2) = Prtcl(:)%S0(2) 
   Prtcl_temp(:)%S0(3) = Prtcl(:)%S0(3) 
   Prtcl_temp(:)%V0(1) = Prtcl(:)%V0(1) 
   Prtcl_temp(:)%V0(2) = Prtcl(:)%V0(2) 
   Prtcl_temp(:)%V0(3) = Prtcl(:)%V0(3) 
   Prtcl_temp(:)%SV0(1) =Prtcl(:)%SV0(1)
   Prtcl_temp(:)%SV0(2) =Prtcl(:)%SV0(2)
   Prtcl_temp(:)%SV0(3) =Prtcl(:)%SV0(3)
   Prtcl_temp(:)%Mass = Prtcl(:)%Mass

   call copy_MC_array(Prtcl_temp, Prtcl) ! below

   siz2 = 2*siz
   deallocate(Prtcl)
   allocate(Prtcl(siz2))
   
   ! Copy the data back into the arrays:
   Prtcl(1:siz)%active = Prtcl_temp(1:siz)%active
   Prtcl(1:siz)%generation = Prtcl_temp(1:siz)%generation
   Prtcl(1:siz)%in_target = Prtcl_temp(1:siz)%in_target
   Prtcl(1:siz)%Ekin = Prtcl_temp(1:siz)%Ekin 
   Prtcl(1:siz)%t0 = Prtcl_temp(1:siz)%t0
   Prtcl(1:siz)%ti = Prtcl_temp(1:siz)%ti
   Prtcl(1:siz)%t_sc = Prtcl_temp(1:siz)%t_sc
   Prtcl(1:siz)%R(1) = Prtcl_temp(1:siz)%R(1)
   Prtcl(1:siz)%R(2) = Prtcl_temp(1:siz)%R(2)
   Prtcl(1:siz)%R(3) = Prtcl_temp(1:siz)%R(3)
   Prtcl(1:siz)%S(1) = Prtcl_temp(1:siz)%S(1) 
   Prtcl(1:siz)%S(2) = Prtcl_temp(1:siz)%S(2) 
   Prtcl(1:siz)%S(3) = Prtcl_temp(1:siz)%S(3) 
   Prtcl(1:siz)%V(1) = Prtcl_temp(1:siz)%V(1) 
   Prtcl(1:siz)%V(2) = Prtcl_temp(1:siz)%V(2) 
   Prtcl(1:siz)%V(3) = Prtcl_temp(1:siz)%V(3) 
   Prtcl(1:siz)%SV(1) = Prtcl_temp(1:siz)%SV(1) 
   Prtcl(1:siz)%SV(2) = Prtcl_temp(1:siz)%SV(2) 
   Prtcl(1:siz)%SV(3) = Prtcl_temp(1:siz)%SV(3) 
   Prtcl(1:siz)%R0(1) = Prtcl_temp(1:siz)%R0(1) 
   Prtcl(1:siz)%R0(2) = Prtcl_temp(1:siz)%R0(2) 
   Prtcl(1:siz)%R0(3) = Prtcl_temp(1:siz)%R0(3) 
   Prtcl(1:siz)%S0(1) = Prtcl_temp(1:siz)%S0(1) 
   Prtcl(1:siz)%S0(2) = Prtcl_temp(1:siz)%S0(2) 
   Prtcl(1:siz)%S0(3) = Prtcl_temp(1:siz)%S0(3) 
   Prtcl(1:siz)%V0(1) = Prtcl_temp(1:siz)%V0(1) 
   Prtcl(1:siz)%V0(2) = Prtcl_temp(1:siz)%V0(2) 
   Prtcl(1:siz)%V0(3) = Prtcl_temp(1:siz)%V0(3) 
   Prtcl(1:siz)%SV0(1) =Prtcl_temp(1:siz)%SV0(1)
   Prtcl(1:siz)%SV0(2) =Prtcl_temp(1:siz)%SV0(2)
   Prtcl(1:siz)%SV0(3) =Prtcl_temp(1:siz)%SV0(3)
   Prtcl(1:siz)%Mass = Prtcl_temp(1:siz)%Mass

   call copy_MC_array(Prtcl(1:siz), Prtcl_temp(1:siz)) ! below
 
 ! All beyond are default to start with:
   do j = siz+1, siz2
      call set_default_particle(Prtcl(j))  ! above
   enddo
   ! clean up:
   deallocate(Prtcl_temp)
end subroutine extend_MC_array_Positrons



pure subroutine extend_MC_array_Atoms(Prtcl)
   type(Atom), dimension(:), allocatable, intent(inout) :: Prtcl    ! all electrons as objects 
   type(Atom), dimension(:), allocatable :: Prtcl_temp     ! temporary aray of electrons
   integer :: siz, siz2, j
   siz = size(Prtcl)
   ! Allocate the temporary array to transiently store data:
   allocate(Prtcl_temp(siz))
   
   Prtcl_temp(:)%active = Prtcl(:)%active
   Prtcl_temp(:)%generation = Prtcl(:)%generation
   Prtcl_temp(:)%in_target = Prtcl(:)%in_target
   Prtcl_temp(:)%Ekin = Prtcl(:)%Ekin 
   Prtcl_temp(:)%t0 = Prtcl(:)%t0
   Prtcl_temp(:)%ti = Prtcl(:)%ti
   Prtcl_temp(:)%t_sc = Prtcl(:)%t_sc
   Prtcl_temp(:)%R(1) = Prtcl(:)%R(1) 
   Prtcl_temp(:)%R(2) = Prtcl(:)%R(2) 
   Prtcl_temp(:)%R(3) = Prtcl(:)%R(3) 
   Prtcl_temp(:)%S(1) = Prtcl(:)%S(1) 
   Prtcl_temp(:)%S(2) = Prtcl(:)%S(2) 
   Prtcl_temp(:)%S(3) = Prtcl(:)%S(3) 
   Prtcl_temp(:)%V(1) = Prtcl(:)%V(1) 
   Prtcl_temp(:)%V(2) = Prtcl(:)%V(2) 
   Prtcl_temp(:)%V(3) = Prtcl(:)%V(3) 
   Prtcl_temp(:)%SV(1) = Prtcl(:)%SV(1) 
   Prtcl_temp(:)%SV(2) = Prtcl(:)%SV(2) 
   Prtcl_temp(:)%SV(3) = Prtcl(:)%SV(3) 
   Prtcl_temp(:)%R0(1) = Prtcl(:)%R0(1) 
   Prtcl_temp(:)%R0(2) = Prtcl(:)%R0(2)
   Prtcl_temp(:)%R0(3) = Prtcl(:)%R0(3) 
   Prtcl_temp(:)%S0(1) = Prtcl(:)%S0(1) 
   Prtcl_temp(:)%S0(2) = Prtcl(:)%S0(2) 
   Prtcl_temp(:)%S0(3) = Prtcl(:)%S0(3) 
   Prtcl_temp(:)%V0(1) = Prtcl(:)%V0(1) 
   Prtcl_temp(:)%V0(2) = Prtcl(:)%V0(2) 
   Prtcl_temp(:)%V0(3) = Prtcl(:)%V0(3) 
   Prtcl_temp(:)%SV0(1) =Prtcl(:)%SV0(1)
   Prtcl_temp(:)%SV0(2) =Prtcl(:)%SV0(2)
   Prtcl_temp(:)%SV0(3) =Prtcl(:)%SV0(3)
   Prtcl_temp(:)%Mass = Prtcl(:)%Mass

   call copy_MC_array(Prtcl_temp, Prtcl) ! below

   siz2 = 2*siz
   deallocate(Prtcl)
   allocate(Prtcl(siz2))
   
   ! Copy the data back into the arrays:
   Prtcl(1:siz)%active = Prtcl_temp(1:siz)%active
   Prtcl(1:siz)%generation = Prtcl_temp(1:siz)%generation
   Prtcl(1:siz)%in_target = Prtcl_temp(1:siz)%in_target
   Prtcl(1:siz)%Ekin = Prtcl_temp(1:siz)%Ekin 
   Prtcl(1:siz)%t0 = Prtcl_temp(1:siz)%t0
   Prtcl(1:siz)%ti = Prtcl_temp(1:siz)%ti
   Prtcl(1:siz)%t_sc = Prtcl_temp(1:siz)%t_sc
   Prtcl(1:siz)%R(1) = Prtcl_temp(1:siz)%R(1)
   Prtcl(1:siz)%R(2) = Prtcl_temp(1:siz)%R(2)
   Prtcl(1:siz)%R(3) = Prtcl_temp(1:siz)%R(3)
   Prtcl(1:siz)%S(1) = Prtcl_temp(1:siz)%S(1) 
   Prtcl(1:siz)%S(2) = Prtcl_temp(1:siz)%S(2) 
   Prtcl(1:siz)%S(3) = Prtcl_temp(1:siz)%S(3) 
   Prtcl(1:siz)%V(1) = Prtcl_temp(1:siz)%V(1) 
   Prtcl(1:siz)%V(2) = Prtcl_temp(1:siz)%V(2) 
   Prtcl(1:siz)%V(3) = Prtcl_temp(1:siz)%V(3) 
   Prtcl(1:siz)%SV(1) = Prtcl_temp(1:siz)%SV(1) 
   Prtcl(1:siz)%SV(2) = Prtcl_temp(1:siz)%SV(2) 
   Prtcl(1:siz)%SV(3) = Prtcl_temp(1:siz)%SV(3) 
   Prtcl(1:siz)%R0(1) = Prtcl_temp(1:siz)%R0(1) 
   Prtcl(1:siz)%R0(2) = Prtcl_temp(1:siz)%R0(2) 
   Prtcl(1:siz)%R0(3) = Prtcl_temp(1:siz)%R0(3) 
   Prtcl(1:siz)%S0(1) = Prtcl_temp(1:siz)%S0(1) 
   Prtcl(1:siz)%S0(2) = Prtcl_temp(1:siz)%S0(2) 
   Prtcl(1:siz)%S0(3) = Prtcl_temp(1:siz)%S0(3) 
   Prtcl(1:siz)%V0(1) = Prtcl_temp(1:siz)%V0(1) 
   Prtcl(1:siz)%V0(2) = Prtcl_temp(1:siz)%V0(2) 
   Prtcl(1:siz)%V0(3) = Prtcl_temp(1:siz)%V0(3) 
   Prtcl(1:siz)%SV0(1) =Prtcl_temp(1:siz)%SV0(1)
   Prtcl(1:siz)%SV0(2) =Prtcl_temp(1:siz)%SV0(2)
   Prtcl(1:siz)%SV0(3) =Prtcl_temp(1:siz)%SV0(3)
   Prtcl(1:siz)%Mass = Prtcl_temp(1:siz)%Mass

   call copy_MC_array(Prtcl(1:siz), Prtcl_temp(1:siz)) ! below
 
 ! All beyond are default to start with:
   do j = siz+1, siz2
      call set_default_particle(Prtcl(j))  ! above
   enddo
   ! clean up:
   deallocate(Prtcl_temp)
end subroutine extend_MC_array_Atoms



pure subroutine extend_MC_array_SHIs(Prtcl)
   type(SHI), dimension(:), allocatable, intent(inout) :: Prtcl    ! all electrons as objects 
   type(SHI), dimension(:), allocatable :: Prtcl_temp     ! temporary aray of electrons
   integer :: siz, siz2, j
   siz = size(Prtcl)
   ! Allocate the temporary array to transiently store data:
   allocate(Prtcl_temp(siz))
   
   Prtcl_temp(:)%active = Prtcl(:)%active
   Prtcl_temp(:)%generation = Prtcl(:)%generation
   Prtcl_temp(:)%in_target = Prtcl(:)%in_target
   Prtcl_temp(:)%Ekin = Prtcl(:)%Ekin 
   Prtcl_temp(:)%t0 = Prtcl(:)%t0
   Prtcl_temp(:)%ti = Prtcl(:)%ti 
   Prtcl_temp(:)%t_sc = Prtcl(:)%t_sc
   Prtcl_temp(:)%R(1) = Prtcl(:)%R(1) 
   Prtcl_temp(:)%R(2) = Prtcl(:)%R(2) 
   Prtcl_temp(:)%R(3) = Prtcl(:)%R(3) 
   Prtcl_temp(:)%S(1) = Prtcl(:)%S(1) 
   Prtcl_temp(:)%S(2) = Prtcl(:)%S(2) 
   Prtcl_temp(:)%S(3) = Prtcl(:)%S(3) 
   Prtcl_temp(:)%V(1) = Prtcl(:)%V(1) 
   Prtcl_temp(:)%V(2) = Prtcl(:)%V(2) 
   Prtcl_temp(:)%V(3) = Prtcl(:)%V(3) 
   Prtcl_temp(:)%SV(1) = Prtcl(:)%SV(1) 
   Prtcl_temp(:)%SV(2) = Prtcl(:)%SV(2) 
   Prtcl_temp(:)%SV(3) = Prtcl(:)%SV(3) 
   Prtcl_temp(:)%R0(1) = Prtcl(:)%R0(1) 
   Prtcl_temp(:)%R0(2) = Prtcl(:)%R0(2)
   Prtcl_temp(:)%R0(3) = Prtcl(:)%R0(3) 
   Prtcl_temp(:)%S0(1) = Prtcl(:)%S0(1) 
   Prtcl_temp(:)%S0(2) = Prtcl(:)%S0(2) 
   Prtcl_temp(:)%S0(3) = Prtcl(:)%S0(3) 
   Prtcl_temp(:)%V0(1) = Prtcl(:)%V0(1) 
   Prtcl_temp(:)%V0(2) = Prtcl(:)%V0(2) 
   Prtcl_temp(:)%V0(3) = Prtcl(:)%V0(3) 
   Prtcl_temp(:)%SV0(1) =Prtcl(:)%SV0(1)
   Prtcl_temp(:)%SV0(2) =Prtcl(:)%SV0(2)
   Prtcl_temp(:)%SV0(3) =Prtcl(:)%SV0(3)
   Prtcl_temp(:)%Mass = Prtcl(:)%Mass

   call copy_MC_array(Prtcl_temp, Prtcl) ! below

   siz2 = 2*siz
   deallocate(Prtcl)
   allocate(Prtcl(siz2))
   
   ! Copy the data back into the arrays:
   Prtcl(1:siz)%active = Prtcl_temp(1:siz)%active
   Prtcl(1:siz)%generation = Prtcl_temp(1:siz)%generation
   Prtcl(1:siz)%in_target = Prtcl_temp(1:siz)%in_target
   Prtcl(1:siz)%Ekin = Prtcl_temp(1:siz)%Ekin 
   Prtcl(1:siz)%t0 = Prtcl_temp(1:siz)%t0
   Prtcl(1:siz)%ti = Prtcl_temp(1:siz)%ti
   Prtcl(1:siz)%t_sc = Prtcl_temp(1:siz)%t_sc
   Prtcl(1:siz)%R(1) = Prtcl_temp(1:siz)%R(1)
   Prtcl(1:siz)%R(2) = Prtcl_temp(1:siz)%R(2)
   Prtcl(1:siz)%R(3) = Prtcl_temp(1:siz)%R(3)
   Prtcl(1:siz)%S(1) = Prtcl_temp(1:siz)%S(1) 
   Prtcl(1:siz)%S(2) = Prtcl_temp(1:siz)%S(2) 
   Prtcl(1:siz)%S(3) = Prtcl_temp(1:siz)%S(3) 
   Prtcl(1:siz)%V(1) = Prtcl_temp(1:siz)%V(1) 
   Prtcl(1:siz)%V(2) = Prtcl_temp(1:siz)%V(2) 
   Prtcl(1:siz)%V(3) = Prtcl_temp(1:siz)%V(3) 
   Prtcl(1:siz)%SV(1) = Prtcl_temp(1:siz)%SV(1) 
   Prtcl(1:siz)%SV(2) = Prtcl_temp(1:siz)%SV(2) 
   Prtcl(1:siz)%SV(3) = Prtcl_temp(1:siz)%SV(3) 
   Prtcl(1:siz)%R0(1) = Prtcl_temp(1:siz)%R0(1) 
   Prtcl(1:siz)%R0(2) = Prtcl_temp(1:siz)%R0(2) 
   Prtcl(1:siz)%R0(3) = Prtcl_temp(1:siz)%R0(3) 
   Prtcl(1:siz)%S0(1) = Prtcl_temp(1:siz)%S0(1) 
   Prtcl(1:siz)%S0(2) = Prtcl_temp(1:siz)%S0(2) 
   Prtcl(1:siz)%S0(3) = Prtcl_temp(1:siz)%S0(3) 
   Prtcl(1:siz)%V0(1) = Prtcl_temp(1:siz)%V0(1) 
   Prtcl(1:siz)%V0(2) = Prtcl_temp(1:siz)%V0(2) 
   Prtcl(1:siz)%V0(3) = Prtcl_temp(1:siz)%V0(3) 
   Prtcl(1:siz)%SV0(1) =Prtcl_temp(1:siz)%SV0(1)
   Prtcl(1:siz)%SV0(2) =Prtcl_temp(1:siz)%SV0(2)
   Prtcl(1:siz)%SV0(3) =Prtcl_temp(1:siz)%SV0(3)
   Prtcl(1:siz)%Mass = Prtcl_temp(1:siz)%Mass

   call copy_MC_array(Prtcl(1:siz), Prtcl_temp(1:siz)) ! below
 
 ! All beyond are default to start with:
   do j = siz+1, siz2
      call set_default_particle(Prtcl(j))  ! above
   enddo
   ! clean up:
   deallocate(Prtcl_temp)
end subroutine extend_MC_array_SHIs


pure subroutine copy_MC_array_SHI(array1, array2)
   type(SHI), dimension(:), intent(inout) :: array1
   type(SHI), dimension(:), intent(inout) :: array2
   array1(:)%Z = array2(:)%Z
   array1(:)%Name = array2(:)%Name
   array1(:)%Zeff = array2(:)%Zeff
   array1(:)%Meff = array2(:)%Meff
   array1(:)%A(1) = array2(:)%A(1) 
   array1(:)%A(2) = array2(:)%A(2) 
   array1(:)%A(3) = array2(:)%A(3) 
   array1(:)%Force(1) = array2(:)%Force(1) 
   array1(:)%Force(2) = array2(:)%Force(2) 
   array1(:)%Force(3) = array2(:)%Force(3) 
end subroutine copy_MC_array_SHI

pure subroutine copy_MC_array_atom(array1, array2)
   type(Atom), dimension(:), intent(inout) :: array1
   type(Atom), dimension(:), intent(inout) :: array2
   array1(:)%Z = array2(:)%Z 
   array1(:)%Name = array2(:)%Name 
   array1(:)%A(1) = array2(:)%A(1) 
   array1(:)%A(2) = array2(:)%A(2) 
   array1(:)%A(3) = array2(:)%A(3) 
   array1(:)%Force(1) = array2(:)%Force(1) 
   array1(:)%Force(2) = array2(:)%Force(2) 
   array1(:)%Force(3) = array2(:)%Force(3) 
end subroutine copy_MC_array_atom
         
pure subroutine copy_MC_array_hole(array1, array2)
   type(Hole), dimension(:), intent(inout) :: array1
   type(Hole), dimension(:), intent(inout) :: array2
   array1(:)%A(1) = array2(:)%A(1) 
   array1(:)%A(2) = array2(:)%A(2) 
   array1(:)%A(3) = array2(:)%A(3) 
   array1(:)%Force(1) = array2(:)%Force(1) 
   array1(:)%Force(2) = array2(:)%Force(2) 
   array1(:)%Force(3) = array2(:)%Force(3) 
   array1(:)%KOA = array2(:)%KOA
   array1(:)%Sh = array2(:)%Sh
!    array1(:)%valent = array2(:)%valent
end subroutine copy_MC_array_hole

pure subroutine copy_MC_array_positron(array1, array2)
   type(Positron), dimension(:), intent(inout) :: array1
   type(Positron), dimension(:), intent(inout) :: array2
   array1(:)%A(1) = array2(:)%A(1) 
   array1(:)%A(2) = array2(:)%A(2) 
   array1(:)%A(3) = array2(:)%A(3) 
   array1(:)%Force(1) = array2(:)%Force(1) 
   array1(:)%Force(2) = array2(:)%Force(2) 
   array1(:)%Force(3) = array2(:)%Force(3) 
end subroutine copy_MC_array_positron

pure subroutine copy_MC_array_electron(array1, array2)
   type(Electron), dimension(:), intent(inout) :: array1
   type(Electron), dimension(:), intent(inout) :: array2
   array1(:)%A(1) = array2(:)%A(1) 
   array1(:)%A(2) = array2(:)%A(2) 
   array1(:)%A(3) = array2(:)%A(3) 
   array1(:)%Force(1) = array2(:)%Force(1) 
   array1(:)%Force(2) = array2(:)%Force(2) 
   array1(:)%Force(3) = array2(:)%Force(3) 
end subroutine copy_MC_array_electron





!--------------------------------------------------------------------------------
! PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

pure subroutine extend_MC_array_poly(Prtcl)
   class(Particle), dimension(:), allocatable, intent(inout) :: Prtcl   ! array of undefined particles
   class(Particle), dimension(:), allocatable :: Prtcl_temp     ! temporary aray of undefined particles
   integer :: siz, siz2, j
   siz = size(Prtcl)
   ! Allocate the temporary array to transiently store data:
   select type (Prtcl)   ! which particle array that should become
   type is (Photon)
      allocate(Photon::Prtcl_temp(siz))
   type is (Electron)
      allocate(Electron::Prtcl_temp(siz))
   type is (Positron)
      allocate(Positron::Prtcl_temp(siz))
   type is (Hole)
      allocate(Hole::Prtcl_temp(siz))
   type is (Atom)
      allocate(Atom::Prtcl_temp(siz))
   type is (SHI)
      allocate(SHI::Prtcl_temp(siz))
   end select
   
   Prtcl_temp(:)%active = Prtcl(:)%active
   Prtcl_temp(:)%generation = Prtcl(:)%generation
   Prtcl_temp(:)%in_target = Prtcl(:)%in_target
   Prtcl_temp(:)%Ekin = Prtcl(:)%Ekin 
   Prtcl_temp(:)%t0 = Prtcl(:)%t0
   Prtcl_temp(:)%ti = Prtcl(:)%ti
   Prtcl_temp(:)%t_sc = Prtcl(:)%t_sc
   Prtcl_temp(:)%R(1) = Prtcl(:)%R(1) 
   Prtcl_temp(:)%R(2) = Prtcl(:)%R(2) 
   Prtcl_temp(:)%R(3) = Prtcl(:)%R(3) 
   Prtcl_temp(:)%S(1) = Prtcl(:)%S(1) 
   Prtcl_temp(:)%S(2) = Prtcl(:)%S(2) 
   Prtcl_temp(:)%S(3) = Prtcl(:)%S(3) 
   Prtcl_temp(:)%V(1) = Prtcl(:)%V(1) 
   Prtcl_temp(:)%V(2) = Prtcl(:)%V(2) 
   Prtcl_temp(:)%V(3) = Prtcl(:)%V(3) 
   Prtcl_temp(:)%SV(1) = Prtcl(:)%SV(1) 
   Prtcl_temp(:)%SV(2) = Prtcl(:)%SV(2) 
   Prtcl_temp(:)%SV(3) = Prtcl(:)%SV(3) 
   Prtcl_temp(:)%R0(1) = Prtcl(:)%R0(1) 
   Prtcl_temp(:)%R0(2) = Prtcl(:)%R0(2)
   Prtcl_temp(:)%R0(3) = Prtcl(:)%R0(3) 
   Prtcl_temp(:)%S0(1) = Prtcl(:)%S0(1) 
   Prtcl_temp(:)%S0(2) = Prtcl(:)%S0(2) 
   Prtcl_temp(:)%S0(3) = Prtcl(:)%S0(3) 
   Prtcl_temp(:)%V0(1) = Prtcl(:)%V0(1) 
   Prtcl_temp(:)%V0(2) = Prtcl(:)%V0(2) 
   Prtcl_temp(:)%V0(3) = Prtcl(:)%V0(3) 
   Prtcl_temp(:)%SV0(1) =Prtcl(:)%SV0(1)
   Prtcl_temp(:)%SV0(2) =Prtcl(:)%SV0(2)
   Prtcl_temp(:)%SV0(3) =Prtcl(:)%SV0(3)
   Prtcl_temp(:)%Mass = Prtcl(:)%Mass
   select type (Prtcl_temp)	! which particle array that should become
      type is (Electron)
         call copy_MC_electron_array(Prtcl_temp, Prtcl) ! below
      type is (Positron)
         call  copy_MC_positron_array(Prtcl_temp, Prtcl)    ! below
      type is (Hole)
         call copy_MC_hole_array(Prtcl_temp, Prtcl)    ! below
      type is (Atom)
         call copy_MC_atom_array(Prtcl_temp, Prtcl)    ! below
      type is (SHI)
         call copy_MC_SHI_array(Prtcl_temp, Prtcl)    ! below
   end select
  
   siz2 = 2*siz
   select type (Prtcl_temp)   ! which particle array that should become
   type is (Photon)
      ! Extend the array:
      deallocate(Prtcl)
      allocate(Photon::Prtcl(siz2))
   type is (Electron)
      allocate(Electron::Prtcl(siz2))
   type is (Positron)
      allocate(Positron::Prtcl(siz2))
   type is (Hole)
      allocate(Hole::Prtcl(siz2))
   type is (Atom)
      allocate(Atom::Prtcl(siz2))
   type is (SHI)
      allocate(SHI::Prtcl(siz2))
   end select
   ! Copy the data back into the arrays:
   Prtcl(1:siz)%active = Prtcl_temp(1:siz)%active
   Prtcl(1:siz)%generation = Prtcl_temp(1:siz)%generation
   Prtcl(1:siz)%in_target = Prtcl_temp(1:siz)%in_target
   Prtcl(1:siz)%Ekin = Prtcl_temp(1:siz)%Ekin 
   Prtcl(1:siz)%t0 = Prtcl_temp(1:siz)%t0
   Prtcl(1:siz)%ti = Prtcl_temp(1:siz)%ti
   Prtcl(1:siz)%t_sc = Prtcl_temp(1:siz)%t_sc
   Prtcl(1:siz)%R(1) = Prtcl_temp(1:siz)%R(1)
   Prtcl(1:siz)%R(2) = Prtcl_temp(1:siz)%R(2)
   Prtcl(1:siz)%R(3) = Prtcl_temp(1:siz)%R(3)
   Prtcl(1:siz)%S(1) = Prtcl_temp(1:siz)%S(1) 
   Prtcl(1:siz)%S(2) = Prtcl_temp(1:siz)%S(2) 
   Prtcl(1:siz)%S(3) = Prtcl_temp(1:siz)%S(3) 
   Prtcl(1:siz)%V(1) = Prtcl_temp(1:siz)%V(1) 
   Prtcl(1:siz)%V(2) = Prtcl_temp(1:siz)%V(2) 
   Prtcl(1:siz)%V(3) = Prtcl_temp(1:siz)%V(3) 
   Prtcl(1:siz)%SV(1) = Prtcl_temp(1:siz)%SV(1) 
   Prtcl(1:siz)%SV(2) = Prtcl_temp(1:siz)%SV(2) 
   Prtcl(1:siz)%SV(3) = Prtcl_temp(1:siz)%SV(3) 
   Prtcl(1:siz)%R0(1) = Prtcl_temp(1:siz)%R0(1) 
   Prtcl(1:siz)%R0(2) = Prtcl_temp(1:siz)%R0(2) 
   Prtcl(1:siz)%R0(3) = Prtcl_temp(1:siz)%R0(3) 
   Prtcl(1:siz)%S0(1) = Prtcl_temp(1:siz)%S0(1) 
   Prtcl(1:siz)%S0(2) = Prtcl_temp(1:siz)%S0(2) 
   Prtcl(1:siz)%S0(3) = Prtcl_temp(1:siz)%S0(3) 
   Prtcl(1:siz)%V0(1) = Prtcl_temp(1:siz)%V0(1) 
   Prtcl(1:siz)%V0(2) = Prtcl_temp(1:siz)%V0(2) 
   Prtcl(1:siz)%V0(3) = Prtcl_temp(1:siz)%V0(3) 
   Prtcl(1:siz)%SV0(1) =Prtcl_temp(1:siz)%SV0(1)
   Prtcl(1:siz)%SV0(2) =Prtcl_temp(1:siz)%SV0(2)
   Prtcl(1:siz)%SV0(3) =Prtcl_temp(1:siz)%SV0(3)
   Prtcl(1:siz)%Mass = Prtcl_temp(1:siz)%Mass
   select type (Prtcl)  ! which particle array that should become
      type is (Electron)
         call copy_MC_electron_array(Prtcl(1:siz), Prtcl_temp(1:siz)) ! below
      type is (Positron)
         call  copy_MC_positron_array(Prtcl(1:siz), Prtcl_temp(1:siz))    ! below
      type is (Hole)
         call copy_MC_hole_array(Prtcl(1:siz), Prtcl_temp(1:siz))    ! below
      type is (Atom)
         call copy_MC_atom_array(Prtcl(1:siz), Prtcl_temp(1:siz))    ! below
      type is (SHI)
         call copy_MC_SHI_array(Prtcl(1:siz), Prtcl_temp(1:siz))    ! below
   end select
   ! All beyond are default to start with:
   do j = siz+1, siz2
      call set_default_particle(Prtcl(j))  ! above
   enddo
   ! clean up:
   deallocate(Prtcl_temp)
end subroutine extend_MC_array_poly


pure subroutine copy_MC_SHI_array(array1, array2)
   type(SHI), dimension(:), intent(inout) :: array1
   class(Particle), dimension(:), intent(inout) :: array2
   select type (array2)   ! which particle array that should become
      type is (SHI)
         array1(:)%Z = array2(:)%Z
         array1(:)%Name = array2(:)%Name
         array1(:)%Zeff = array2(:)%Zeff
         array1(:)%Meff = array2(:)%Meff
         array1(:)%A(1) = array2(:)%A(1) 
         array1(:)%A(2) = array2(:)%A(2) 
         array1(:)%A(3) = array2(:)%A(3) 
         array1(:)%Force(1) = array2(:)%Force(1) 
         array1(:)%Force(2) = array2(:)%Force(2) 
         array1(:)%Force(3) = array2(:)%Force(3) 
   endselect
end subroutine copy_MC_SHI_array

pure subroutine copy_MC_atom_array(array1, array2)
   type(Atom), dimension(:), intent(inout) :: array1
   class(Particle), dimension(:), intent(inout) :: array2
   select type (array2)   ! which particle array that should become
      type is (Atom)
         array1(:)%Z = array2(:)%Z 
         array1(:)%Name = array2(:)%Name 
         array1(:)%A(1) = array2(:)%A(1) 
         array1(:)%A(2) = array2(:)%A(2) 
         array1(:)%A(3) = array2(:)%A(3) 
         array1(:)%Force(1) = array2(:)%Force(1) 
         array1(:)%Force(2) = array2(:)%Force(2) 
         array1(:)%Force(3) = array2(:)%Force(3) 
   endselect
end subroutine copy_MC_atom_array
         
pure subroutine copy_MC_hole_array(array1, array2)
   type(Hole), dimension(:), intent(inout) :: array1
   class(Particle), dimension(:), intent(inout) :: array2
   select type (array2)   ! which particle array that should become
      type is (Hole)
         array1(:)%A(1) = array2(:)%A(1) 
         array1(:)%A(2) = array2(:)%A(2) 
         array1(:)%A(3) = array2(:)%A(3) 
         array1(:)%Force(1) = array2(:)%Force(1) 
         array1(:)%Force(2) = array2(:)%Force(2) 
         array1(:)%Force(3) = array2(:)%Force(3) 
         array1(:)%KOA = array2(:)%KOA
         array1(:)%Sh = array2(:)%Sh 
         array1(:)%valent = array2(:)%valent
   endselect
end subroutine copy_MC_hole_array

pure subroutine copy_MC_positron_array(array1, array2)
   type(Positron), dimension(:), intent(inout) :: array1
   class(Particle), dimension(:), intent(inout) :: array2
   select type (array2)   ! which particle array that should become
      type is (Positron)
         array1(:)%A(1) = array2(:)%A(1) 
         array1(:)%A(2) = array2(:)%A(2) 
         array1(:)%A(3) = array2(:)%A(3) 
         array1(:)%Force(1) = array2(:)%Force(1) 
         array1(:)%Force(2) = array2(:)%Force(2) 
         array1(:)%Force(3) = array2(:)%Force(3) 
   endselect
end subroutine copy_MC_positron_array

pure subroutine copy_MC_electron_array(array1, array2)
   type(Electron), dimension(:), intent(inout) :: array1
   class(Particle), dimension(:), intent(inout) :: array2
   select type (array2)   ! which particle array that should become
      type is (Electron)
         array1(:)%A(1) = array2(:)%A(1) 
         array1(:)%A(2) = array2(:)%A(2) 
         array1(:)%A(3) = array2(:)%A(3) 
         array1(:)%Force(1) = array2(:)%Force(1) 
         array1(:)%Force(2) = array2(:)%Force(2) 
         array1(:)%Force(3) = array2(:)%Force(3) 
   endselect
end subroutine copy_MC_electron_array

   
end module MC_general_tools
