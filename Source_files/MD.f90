! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains general subroutines for molecular dynamics simulations:
! [3] D. Wolf, et al., J. Chem. Phys. 110, 8254 (1999); https://doi.org/10.1063/1.478738
! [4] Hansen et al., J. Phys. Chem. B 116, 5738 (2012); https://doi.org/10.1021/JP300750G

module MD
use Universal_constants
use Objects
use Dealing_with_files, only: close_file
use Geometries, only: Cartesian_to_cylindrical, Cartesian_to_spherical, fit_parabola_to_3points
use Little_subroutines, only: Find_in_array_monoton, print_time_step
use MD_general_tools, only: Verlet_step, get_nearest_neighbors_list, Quenching, find_which_potential, &
                            Fermi_cut_off, d_Fermi_cut_off, Check_MD_periodic_boundaries, rescale_supercell, &
                            get_temperature_from_energy, get_pressure, shortest_distance, Thermostat, Damp_pressure, &
                            rescale_velosities, get_total_energy, transfer_atomic_data
use MD_data_analysis, only: get_total_energies
use MD_Pot_Simple, only: power_potential, d_power_potential, LJ_potential, d_LJ_potential, exp_potential, d_exp_potential
use MD_Pot_Buck, only: Buck_potential, d_Buck_potential, Matsui_potential, d_Matsui_potential
use MD_Pot_Coulomb, only: Bare_Coulomb_pot, d_Bare_Coulomb_pot, Ewalds_Coulomb, Coulomb_Wolf_pot, d_Coulomb_Wolf_pot, &
                            Coulomb_Wolf_self_term
use MD_Pot_SW, only: SW_potential, SW_U_ij, SW_dU_ij, d_SW_potential
use MD_Pot_ZBL, only: ZBL_pot, d_ZBL_pot

implicit none

real(8) :: m_acceleration_factor
character(50) :: m_cohesive_file

parameter (m_cohesive_file = 'OUTPUT_MD_cohesive_energy.txt')   ! file name for cohesive energy
parameter (m_acceleration_factor = g_e * 1.0d-10)   ! factor to convert acceleration from [eV/(A*kg)]->[A/fs^2]

 contains


subroutine MD_step(MD_atoms, MD_supce, MD_pots, numpar, t_cur, dt)
   type(Atom), dimension(:), intent(inout), allocatable :: MD_atoms ! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   type(MD_potential), dimension(:,:), intent(in), allocatable :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   real(8), intent(in) :: t_cur             ! [fs] current time
   real(8), intent(in) :: dt                ! [fs] timestep
   !-------------------------
   real(8) :: dt2, dt_half  ! dt^2, dt^2/2, for speeding up calculations
   integer :: Na,  i

   ! Check if user requested MD at all:
   if (.not.numpar%DO_MD) goto 9998 ! skip entire subroutine
   
   ! Parameters to be used in MD propagator (e.g. Verlet):
   dt_half = dt*0.5d0   ! dt/2
!    dt2 = dt*dt_half     ! dt^2/2
   Na = size(MD_atoms)  ! how many atoms we have

   !------------------------------
   ! I. Make an MD step for all atoms:
   ! The half-step for velocities, and full step for coordinates:
   select case(numpar%MD_integrator)   ! which integrator to use
   case default ! Verlet
      !$omp parallel private (i)
      !$omp do
      do i = 1, Na
         call Verlet_step(MD_atoms(i)%R, MD_atoms(i)%R0, MD_atoms(i)%V, MD_atoms(i)%V0, MD_atoms(i)%A, &
                    dt, dt_half, .true.) ! module "MD_general_tools"
      enddo
      !$omp enddo
      !$omp end parallel
   case (2) ! something else 
      ! not ready
   endselect
   !------------------------------

   ! 1) Place atoms back inside the supercell with periodic boundaries, if needed:
   call Check_MD_periodic_boundaries(MD_atoms, MD_supce)    ! module "MD_general_tools"
!    call print_time_step('periodic_boundaries done:', t_cur, msec=.true.)   ! module "Little_subroutines"

   ! 1.a) Update the list of nearest neighbors after atoms have moved:
   call get_nearest_neighbors_list(MD_atoms, MD_supce, MD_pots, numpar) ! module "MD_general_tools"
!    call print_time_step('nearest_neighbors done:', t_cur, msec=.true.)   ! module "Little_subroutines"

   ! 2) Now, update MD potential and forces acting on all atoms:
   call get_MD_pot_n_force(numpar, MD_atoms, MD_supce, MD_pots)   ! below
!    call print_time_step('pot_n_force done:', t_cur, msec=.true.)   ! module "Little_subroutines"

   ! 3) Get energy transfered from electrons (from MC module):
   call get_energy_from_electrons(MD_atoms, MD_supce)   ! below

   ! 4) Do quenching, if user requested:
   call Quenching(numpar, MD_atoms, t_cur) ! module "MD_general_tools"

   ! 5) Use thermostat, if user requested:
   call Thermostat(numpar, MD_atoms, MD_supce, dt) ! module "MD_general_tools"

   ! 5) Update supercell parameters:
   ! Get average kinetic and potential energies of all atoms:
   call get_total_energies(MD_atoms, MD_supce%Ekin_tot, MD_supce%Epot_tot)    ! module "MD_general_tools"
   ! Get average temeprature:
   call get_temperature_from_energy(dble(Na), MD_supce%Ekin_tot*dble(Na), MD_supce%T_kin, .true., .true.)   ! module "MD_general_tools"
   ! Get average pressure tensor:
   call get_pressure(MD_atoms, MD_supce, MD_supce%Pressure, .true.)   ! module "MD_general_tools"

   !------------------------------
   ! II. Update the first half of velosities for the next Verlet step:
   select case(numpar%MD_integrator)   ! which integrator to use
   case default ! Verlet
      !$omp parallel private (i)
      !$omp do
      do i = 1, Na
         call Verlet_step(MD_atoms(i)%R, MD_atoms(i)%R0, MD_atoms(i)%V, MD_atoms(i)%V0, MD_atoms(i)%A, &
                    dt, dt_half, .false.) ! module "MD_general_tools"
      enddo
      !$omp enddo
      !$omp end parallel
   case (2) ! something else 
      ! not ready
   endselect
   !------------------------------

9998 continue
end subroutine MD_step


subroutine get_energy_from_electrons(MD_atoms, MD_supce)
   type(Atom), dimension(:), intent(inout) :: MD_atoms ! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   !---------------------------------
   real(8) :: E_transferred, Ekin, Epot, R_c, Theta_c, Phi_c, L_c
   integer :: Nat, i_at, i_dim, i, j, k
   integer, dimension(3, size(MD_atoms)) :: At_indices
   logical, dimension(size(MD_atoms)) :: At_mask


   if (MD_supce%coord_dim > 0) then ! at least 1d grid:
      Nat = size(MD_atoms)  ! number of atoms in the supercell
      !$omp workshare
      At_indices = 1    ! to start with
      !$omp end workshare

      ! Distribute all the atos into corresponding cells on the grid of energy transfer:
      !$omp parallel private (i_at, i_dim, R_c, Theta_c, Phi_c, L_c, i, j, k, At_mask, Ekin, Epot, E_transferred)
      !$omp do schedule(dynamic)
      do i_at = 1, Nat ! for all atoms
         do i_dim = 1, MD_supce%coord_dim   ! along all dimensions
            ! Find position of event on the grid:
            select case (MD_supce%coord_type)
            case (1)   ! Cartesian
               select case (MD_supce%coord_axis(i_dim))  ! along which axis
               case (1)   ! X
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i_dim)%grid, MD_atoms(i_at)%R(1), &
                                                At_indices(i_dim, i_at)) ! module "Little_subroutines"
               case (2)   ! Y
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i_dim)%grid, MD_atoms(i_at)%R(2), &
                                                At_indices(i_dim, i_at)) ! module "Little_subroutines"
               case (3)   ! Z
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i_dim)%grid, MD_atoms(i_at)%R(3), &
                                                At_indices(i_dim, i_at)) ! module "Little_subroutines"
               endselect
            case (2)   ! Spherical
               ! Get the spherical coords:
               call Cartesian_to_spherical(MD_atoms(i_at)%R(1), MD_atoms(i_at)%R(2), MD_atoms(i_at)%R(3), &
                                            R_c, Theta_c, Phi_c)  ! module "Geometries"
               ! Find where on the grid it is:
               select case (MD_supce%coord_axis(i_dim))  ! along which axis
               case (1)   ! R
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i_dim)%grid, R_c, At_indices(i_dim, i_at)) ! module "Little_subroutines"
               case (2)   ! Theta
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i_dim)%grid, Theta_c, At_indices(i_dim, i_at)) ! module "Little_subroutines"
               case (3)   ! Phi
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i_dim)%grid, Phi_c, At_indices(i_dim, i_at)) ! module "Little_subroutines"
               endselect
            case (3)   ! Cylindrical
               call Cartesian_to_cylindrical(MD_atoms(i_at)%R(1), MD_atoms(i_at)%R(2), MD_atoms(i_at)%R(3), &
                                                R_c, L_c, Phi_c) ! module "Geometries"
               select case (MD_supce%coord_axis(i_dim))  ! along which axis
               case (1)   ! R
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i_dim)%grid, R_c, At_indices(i_dim, i_at)) ! module "Little_subroutines"
               case (2)   ! Theta
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i_dim)%grid, L_c, At_indices(i_dim, i_at)) ! module "Little_subroutines"
               case (3)   ! Phi
                  call Find_in_array_monoton(MD_supce%MC2MD_dimen(i_dim)%grid, Phi_c, At_indices(i_dim, i_at)) ! module "Little_subroutines"
               endselect
            endselect
         enddo ! i_dim
      enddo ! i_at
      !$omp end do
      !$omp barrier

      ! Now rescale atomc velocities in each subcell:
      !$omp do schedule(dynamic)
      do i = 1, size(MD_supce%MC2MD_dimen(1)%grid)
         do j = 1, size(MD_supce%MC2MD_dimen(2)%grid)
            do k = 1, size(MD_supce%MC2MD_dimen(3)%grid) ! for all axes
               ! Construct the mask for atoms to sort them into subcells:
               At_mask(:) = (At_indices(1,:) == i) .and. (At_indices(2,:) == j) .and. (At_indices(3,:) == k)

               ! First, get the total kinetic energy of atoms inside this subcell:
               call get_total_energy(MD_atoms, Ekin, Epot, At_mask)  ! module "MD_general_tools"

               ! Get the total energy transfered in this subcell:
               E_transferred = MD_supce%E_e_at_from_MC(i,j,k) + MD_supce%E_h_at_from_MC(i,j,k) + &
                      MD_supce%E_p_at_from_MC(i,j,k) + MD_supce%E_e_from_MC(i,j,k) + &
                      MD_supce%E_h_from_MC(i,j,k)

               ! Then, rescale the velocities to add energy to this:
               call rescale_velosities(MD_atoms, Ekin, E_transferred, At_mask) ! module "MD_general_tools"
            enddo ! k
         enddo ! j
      enddo ! i
      !$omp end do
      !$omp end parallel

   else   ! zero-dimensional modeling: all energy is distributed among all atoms homogeneously
      ! Total energy transferred from MC is:
      E_transferred = MD_supce%E_e_at_from_MC(1,1,1) + MD_supce%E_h_at_from_MC(1,1,1) + &
                      MD_supce%E_p_at_from_MC(1,1,1) + MD_supce%E_e_from_MC(1,1,1) + &
                      MD_supce%E_h_from_MC(1,1,1)
      ! Find total kinetic energy of atoms (need it for velocity scaling):
      call get_total_energy(MD_atoms, Ekin, Epot) ! module "MD_general_tools"
!      print*, 'get_energy_from_electrons', Ekin, E_transferred

      ! Rescale atomic velosities to deliver the energy to atoms:
      call rescale_velosities(MD_atoms, Ekin, E_transferred) ! module "MD_general_tools"
   endif
end subroutine get_energy_from_electrons




subroutine prepare_MD_run(MD_atoms, MD_supce, MD_pots, numpar)
   type(Atom), dimension(:), intent(inout), allocatable :: MD_atoms ! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   type(MD_potential), dimension(:,:), intent(in), allocatable :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   !-------------------------
   real(8) :: Nat, dt_half
   integer :: i

   if (numpar%DO_MD) then   ! only if user requested MD

      Nat = dble(size(MD_atoms))
   
      ! Place atoms inside the supercell with periodic boundaries, if needed:
      call Check_MD_periodic_boundaries(MD_atoms, MD_supce)    ! module "MD_general_tools"

      ! Create the list of nearest neighbors after atoms have moved:
      call get_nearest_neighbors_list(MD_atoms, MD_supce, MD_pots, numpar) ! module "MD_general_tools"

      ! Get MD potential and forces acting on all atoms:
      call get_MD_pot_n_force(numpar, MD_atoms, MD_supce, MD_pots)   ! below
   
      ! Get the half of velosities for Verlet algorithm start:
      select case(numpar%MD_integrator)   ! which integrator to use
      case default ! Verlet
         dt_half = numpar%dt_MD*0.5d0   ! dt/2
         !$omp parallel private (i)
         !$omp do
         do i = 1, int(Nat)
            call Verlet_step(MD_atoms(i)%R, MD_atoms(i)%R0, MD_atoms(i)%V, MD_atoms(i)%V0, MD_atoms(i)%A, &
                    numpar%dt_MD, dt_half, .false.) ! module "MD_general_tools"
         enddo
         !$omp enddo
         !$omp end parallel
      case (2) ! something else
         ! not ready
      endselect

      ! Save initial state of the supercell for analysis:
      call transfer_atomic_data(MD_atoms, MD_supce%MD_atoms_0)   ! module "MD_general_tools"

      !------------------------------
      ! Get starting supercell parameters:
      ! Get average kinetic and potential energies of all atoms:
      call get_total_energies(MD_atoms, MD_supce%Ekin_tot, MD_supce%Epot_tot)    ! module "MD_general_tools"
      ! Get average temeprature:
      call get_temperature_from_energy(Nat, MD_supce%Ekin_tot*Nat, MD_supce%T_kin, .true., .true.)   ! module "MD_general_tools"
      ! Get average pressure tensor:
      call get_pressure(MD_atoms, MD_supce, MD_supce%Pressure, .true.)   ! module "MD_general_tools"
   endif ! (numpar%DO_MD)
end subroutine prepare_MD_run



subroutine calculate_cohesive_energy(numpar, MD_atoms, MD_supce, MD_pots) ! cohesive energy calculation
   type(Num_par), intent(inout) :: numpar   ! all numerical parameters
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   !-----------------------------
   integer :: i, FN, Nsiz
   real(8) :: rescal, f_start, f_end, f_step, MD_Ekin, MD_Epot
   real(8), dimension(3,3) :: Press_tens
   character(250) :: file_cohesive
   type(Atom), dimension(size(MD_atoms)) :: MD_atoms_SAVE   ! all atoms in MD as objects
   type(MD_supcell) :: MD_supce_SAVE  ! MD supercell parameters
   
   if (numpar%DO_MD) then
      print*, 'Calculating cohesive energy of MD supercell...'
   else
      print*, 'Cohesive energy cannot be calcualted without MD included'
      goto 9999 ! terminate
   endif

   ! Prepare the file to write out cohesive energy:
   file_cohesive = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_cohesive_file))
   open(newunit = FN, FILE = trim(adjustl(file_cohesive)))
   write(FN, '(a)') '#Rescaling Nearest_neighbors Energy    Pressure_xx Pressure_yy Pressure_zz'
   write(FN, '(a)') '#  A   eV/atom GPa GPa GPa'
   
   ! Save unperturbed supercell:
   MD_supce_SAVE = MD_supce
   MD_atoms_SAVE = MD_atoms
   
   f_start = 0.75d0  ! minimal rescaling factor for the size of supercell
   f_end = 8.0d0    ! maximal rescaling factor for the size of supercell
   f_step = 1.0d0/200.0d0
   Nsiz = CEILING( (f_end - f_start)/f_step )  ! how many steps
   do i = 1, Nsiz
      ! Get the energy for the supercell rescaled by this factor:
      rescal = (f_start + dble(i)*f_step)**(1.0d0/3.0d0)    ! rescaling factor
      ! Rescale supercell size (and correspondingly coordinates of all atoms):
      call rescale_supercell(MD_atoms_SAVE, MD_atoms, MD_supce_SAVE, MD_supce, rescal)  ! module "MD_general_tools"
!       print*, 'FIX', MD_atoms_SAVE(1)%R(:)
!       print*, 'MOV', MD_atoms(1)%R(:)
      ! Update the list of nearest neighbors after rescaling:
      call get_nearest_neighbors_list(MD_atoms, MD_supce, MD_pots, numpar) ! module "MD_general_tools"
      ! Get the potential energy for rescaled supercell:
      call get_MD_pot_n_force(numpar, MD_atoms, MD_supce, MD_pots)    ! below
      ! Get total energy:
      call get_total_energies(MD_atoms, MD_Ekin, MD_Epot)    ! module "MD_data_analysis"
      ! Get pressures:
      call get_pressure(MD_atoms, MD_supce, Press_tens, .true.)     ! module "MD_general_tools"
      ! Printout the cohesive energy:
      write(FN, '(f,f,es,es,es,es)') rescal, minval(numpar%Neighbors_Rij(1,:)), MD_Epot, Press_tens(1,1), Press_tens(2,2), Press_tens(3,3)
      write(*, '(f,f,es,es,es,es)') rescal, minval(numpar%Neighbors_Rij(1,:)), MD_Epot, Press_tens(1,1), Press_tens(2,2), Press_tens(3,3)
   enddo
   
9999   print*, 'It is done, TREKIS terminates'
   call close_file('save',FN=FN)    ! module "Dealing_with_files"
end subroutine calculate_cohesive_energy


subroutine get_MD_pot_n_force(numpar, MD_atoms, MD_supce, MD_pots)
   type(Num_par), intent(inout), target :: numpar   ! all numerical parameters
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   !-----------------------------
   integer :: i, j, Nat, i_pot
   logical :: long_range_present

   Nat = size(MD_atoms) ! all atoms
   MD_supce%P_pot(:,:) = 0.0d0    ! to start with

!    call print_time_step('get_MD_pot_n_force 0:', 0.0, msec=.true.)   ! module "Little_subroutines"

   !------------------------------------
   ! I. All short-rage potentials:

   ! Chose the force calculation scheme:
   select case (numpar%MD_force_ind)
   case default   ! analytical integration
      call analytical_MD_pot_n_force(numpar, MD_atoms, MD_supce, MD_pots)   ! below
   case (1) ! if we use a numerical derivative of the potential to construct a force:
      call numerical_MD_pot_n_force(numpar, MD_atoms, MD_supce, MD_pots) ! below
   end select

   !------------------------------------
   ! II. All long-range potentials (requiring explicit image cells):
   ! (use Ewalds summation for Coulomb interactions):
   ! Find if there is any long-range potential:
   long_range_present = .false. ! to start with
   LRS:do i = 1, size(MD_pots,1)    ! check all pair of atoms
      do j = 1, size(MD_pots,2)
         do i_pot = 1, size(MD_pots(i,j)%Set) ! and all potentials for them
            ASSOCIATE (MDPot => MD_pots(i,j)%Set(i_pot)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(MDPot)
               type is (Coulomb_Ewald)   ! Coulomb as long range potential
                  long_range_present = .true. ! long-range potential found
                  exit LRS    ! exit potential search
               end select
            END ASSOCIATE
         enddo ! i_pot
      enddo ! j
   enddo LRS ! i
   ! If we have long range (Coulomb) potential, use special treatement for it:
   if (long_range_present) then
      call Ewalds_Coulomb(numpar, MD_atoms, MD_supce, MD_pots)  ! module "MD_Pot_Coulomb"
   endif
!    call print_time_step('get_MD_pot_n_force 2:', 2.0, msec=.true.)   ! module "Little_subroutines"

   !------------------------------------
   ! III. Also check if there are additional forces from pressure dampener:
   call Damp_pressure(numpar, MD_atoms, MD_supce) ! module "MD_general_tools"
!    call print_time_step('get_MD_pot_n_force 3:', 3.0, msec=.true.)   ! module "Little_subroutines"

   !------------------------------------
   ! IV. Once total forces are known, get accelerations:
   !$omp parallel private (i)
   !$omp do
   do i = 1, Nat
      MD_atoms(i)%A(:) = MD_atoms(i)%Force(:)/MD_atoms(i)%Mass * m_acceleration_factor  ! [A/fs^2]
   enddo ! i
   !$omp enddo
   !$omp end parallel
   
!    call print_time_step('get_MD_pot_n_force 4:', 4.0, msec=.true.)   ! module "Little_subroutines"
!    print*, 'Total F:',  MD_atoms(1)%Force(1)
!     pause 'get_MD_pot_n_force'
end subroutine get_MD_pot_n_force



subroutine numerical_MD_pot_n_force(numpar, MD_atoms, MD_supce, MD_pots)
   type(Num_par), intent(in), target :: numpar   ! all numerical parameters
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   !-------------------------------
   real(8), dimension(size(MD_atoms),size(MD_atoms)) :: U_ij, U_ij_SAVE  ! repeated contributions for three-body potentials (e.g. SW)
   real(8), dimension(size(MD_atoms),size(MD_atoms)) :: dU_ij ! derivative of three-body potentials (e.g. SW)
   real(8), dimension(size(MD_atoms),3) :: Ux, Uy, Uz
   real(8) :: r, x, y, z, r_disp, r_sign, A, B, C, Force(3), R_ij(3), f_cut, Pot_part, Pot_3bdy
   integer :: Nat, i, j, i_disp, i_axis, KOP1, KOP2, N_pot, i_pot, k
   integer, pointer :: atom_2
   logical :: ion_done

   r_disp = 1.0d-3   ! [A] default displacement for numerical evaluation of the derivative of the potential

   Nat = size(MD_atoms) ! all atoms
   Ux = 0.0d0  ! to start with
   Uy = 0.0d0  ! to start with
   Uz = 0.0d0  ! to start with


   ! If there are three-body_potentials, set and save repeating functions for them:
   call get_repeated_three_body_part(numpar, MD_atoms, MD_pots, MD_supce, U_ij)   ! below
   U_ij_SAVE = U_ij  ! save to reuse

   ! For each pair of interacting atoms:
   !$omp parallel private (i, ion_done, j, atom_2, r, x, y, z, R_ij, KOP1, KOP2, &
   !$                      N_pot, i_pot, f_cut, Pot_part, Pot_3bdy)
   !$omp do
   do i = 1, Nat
      ion_done = .false.    ! to start with
      MD_atoms(i)%U = 0.0d0         ! reset the potential to start with
      ! Find nearest neighbors it is interacting with:
      do j = 1, numpar%Neighbors_Num(i)    ! that's how many atoms the atom "i" interacts with
         atom_2 => numpar%Neighbors_list(i, j)     ! that's index of the atom, which "i" is interacting with
         ! Find te shortest distance between atoms i and atom_2:
         call shortest_distance(i, atom_2, MD_atoms, MD_supce, r, x1=R_ij(1), y1=R_ij(2), z1=R_ij(3)) ! module "MD_general_tools"
         if (r < 1.0d-6) print*, 'Error in numerical_MD_pot_n_force #1: distance too small:', r, R_ij(:)

         ! Find indices of the potential(s) is used for this pair of atoms:
         call find_which_potential(MD_atoms, i, atom_2, MD_pots, KOP1, KOP2)    ! module "MD_general_tools"

         ! How many potentials are used for this pair of atoms:
         Pot_part = 0.0d0  ! to start with
         N_pot = size(MD_pots(KOP1,KOP2)%Set)
         ! Find parameters of all potentials between the two atoms:
         do i_pot = 1, N_pot
            ! Define cut-off function and its derivative, to get total forces:
            f_cut = Fermi_cut_off(r, MD_pots(KOP1,KOP2)%Set(i_pot)%Par%d_cut, MD_pots(KOP1,KOP2)%Set(i_pot)%Par%dd) ! module "MD_general_tools"

            ! Calculate the potential and its derivative by r:
            Pot_3bdy = 0.0d0  ! to start with
            call get_single_potential(MD_pots(KOP1,KOP2)%Set(i_pot)%Par, r, ion_done, Pot_part, &
                     numpar, i, atom_2, R_ij, U_ij, dU_ij, Pot_3bdy, MD_atoms, MD_supce) ! below

            ! Combine the total potential from potential and cut-off function:
            MD_atoms(i)%U = MD_atoms(i)%U + Pot_part*f_cut + Pot_3bdy   ! [eV] add to total potential of this atom

         enddo ! i_pot
      enddo ! j = numpar%Neighbors_Num(i)
      ! Central point for extrapolation is at shift 0:
      Ux(i,2) = MD_atoms(i)%U
      Uy(i,2) = MD_atoms(i)%U
      Uz(i,2) = MD_atoms(i)%U
   enddo ! i
   !$omp enddo
   !$omp end parallel
   nullify(atom_2)


   ! Calcualte potential for displacements in all directions:
   do i_axis = 1, 3   ! x, y, z
      do i_disp = 1, 3   ! -1, 0, +1
         ! Set the direction of displacement:
         select case (i_disp)
         case (1)
            r_sign = -1.0d0
         case (3)
            r_sign = 1.0d0
         case default
            r_sign = 0.0d0
         end select

         if (i_disp /= 2) then   ! there is a displacement

            !$omp parallel private (i, ion_done, U_ij, j, atom_2, r, x, y, z, R_ij, KOP1, KOP2, &
            !$                      N_pot, i_pot, f_cut, Pot_part, Pot_3bdy)
            !$omp do
            do i = 1, Nat
               ion_done = .false.    ! to start with

               U_ij = U_ij_SAVE ! reuse
               call get_repeated_three_body_part(numpar, MD_atoms, MD_pots, MD_supce, U_ij, atom_i=i, disp=r_sign*r_disp, ind_disp=i_axis)   ! below

               ! Find nearest neighbors it is interacting with:
               do j = 1, numpar%Neighbors_Num(i)    ! that's how many atoms the atom "i" interacts with
                  atom_2 => numpar%Neighbors_list(i, j)     ! that's index of the atom, which "i" is interacting with
                  ! Find te shortest distance between atoms i and atom_2:
                  call shortest_distance(i, atom_2, MD_atoms, MD_supce, r, x1=R_ij(1), y1=R_ij(2), z1=R_ij(3), &
                                          disp = r_sign*r_disp, ind_disp = i_axis) ! module "MD_general_tools"
                  if (r < 1.0d-6) print*, 'Error in numerical_MD_pot_n_force #2: distance too small:', r, R_ij(:)

!                   ! Testing:
!                   if ((i == 1) .and. (atom_2 == 5)) print*, R_ij(:)

                  ! Find indices of the potential(s) is used for this pair of atoms:
                  call find_which_potential(MD_atoms, i, atom_2, MD_pots, KOP1, KOP2)    ! module "MD_general_tools"

                  ! How many potentials are used for this pair of atoms:
                  Pot_part = 0.0d0  ! to start with
                  N_pot = size(MD_pots(KOP1,KOP2)%Set)
                  ! Find parameters of all potentials between the two atoms:
                  do i_pot = 1, N_pot
                     ! Define cut-off function and its derivative, to get total forces:
                     f_cut = Fermi_cut_off(r, MD_pots(KOP1,KOP2)%Set(i_pot)%Par%d_cut, MD_pots(KOP1,KOP2)%Set(i_pot)%Par%dd) ! module "MD_general_tools"

                     ! Calculate the potential and its derivative by r:
                     Pot_3bdy = 0.0d0  ! to start with
                     call get_single_potential(MD_pots(KOP1,KOP2)%Set(i_pot)%Par, r, ion_done, Pot_part, &
                              numpar, i, atom_2, R_ij, U_ij, dU_ij, Pot_3bdy, MD_atoms, MD_supce, &
                              disp = r_sign*r_disp, ind_disp = i_axis) ! below

                     ! Combine the total potential from potential and cut-off function:
                     !MD_atoms(i)%U = MD_atoms(i)%U + Pot_part*f_cut + Pot_3bdy   ! [eV] add to total potential of this atom
                     select case (i_axis)
                     case (1) ! X
                        Ux(i,i_disp) = Ux(i,i_disp) + Pot_part*f_cut + Pot_3bdy   ! [eV] add to total potential of this atom
                     case (2) ! Y
                        Uy(i,i_disp) = Uy(i,i_disp) + Pot_part*f_cut + Pot_3bdy   ! [eV] add to total potential of this atom
                     case (3) ! Z
                        Uz(i,i_disp) = Uz(i,i_disp) + Pot_part*f_cut + Pot_3bdy   ! [eV] add to total potential of this atom
                     end select

                  enddo ! i_pot
               enddo ! j = numpar%Neighbors_Num(i)
            enddo ! i = 1, Nat
            !$omp enddo
            !$omp end parallel

         endif ! (i_disp /= 0)
      enddo ! i_disp
   enddo ! i_axis


   ! Calculate numerical derivative of the potential:
   !$omp parallel private (i, A, B, C, Force, k)
   !$omp do
   do i = 1, Nat

      ! Along X:
      call fit_parabola_to_3points(-r_disp, Ux(i,1), 0.0d0, Ux(i,2), r_disp, Ux(i,3), A, B, C)  ! module "Geometries"
      ! Derivative of parabola is linear function:
      Force(1) = -B   ! [eV/A]

      ! Along Y:
      call fit_parabola_to_3points(-r_disp, Uy(i,1), 0.0d0, Uy(i,2), r_disp, Uy(i,3), A, B, C)  ! module "Geometries"
      Force(2) = -B   ! [eV/A]

      ! Along Z:
      call fit_parabola_to_3points(-r_disp, Uz(i,1), 0.0d0, Uz(i,2), r_disp, Uz(i,3), A, B, C)  ! module "Geometries"
      Force(3) = -B   ! [eV/A]

      ! Save the total force on atom i:
      MD_atoms(i)%Force(:) = Force(:)

      ! Also save potential part of pressure in the supercell:
      do k = 1, 3
         MD_supce%P_pot(1,k) = MD_supce%P_pot(1,k) + Force(k)*MD_atoms(i)%R(1)    ! [eV]
         MD_supce%P_pot(2,k) = MD_supce%P_pot(2,k) + Force(k)*MD_atoms(i)%R(2)    ! [eV]
         MD_supce%P_pot(3,k) = MD_supce%P_pot(3,k) + Force(k)*MD_atoms(i)%R(3)    ! [eV]
      enddo

   enddo
   !$omp enddo
   !$omp end parallel

end subroutine numerical_MD_pot_n_force



subroutine analytical_MD_pot_n_force(numpar, MD_atoms, MD_supce, MD_pots)
   type(Num_par), intent(in), target :: numpar   ! all numerical parameters
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   !-------------------
   integer :: i, j, k, Nat, i_pot, N_pot, KOP1, KOP2
   real(8) :: Pot_part, f_cut, d_Pot_part, d_f_cut, tot_derive, Pot_3bdy, R_ij(3)
   real(8) :: cos_a, cos_b, cos_c, Force(3), Force_3bdy(3), Force_3bdy_tmp(3)
   logical :: ion_done
   real(8), pointer :: r, x, y, z
   integer, pointer :: atom_2
   real(8), dimension(size(MD_atoms),size(MD_atoms)) :: U_ij ! repeated contributions for three-body potentials (e.g. SW)
   real(8), dimension(size(MD_atoms),size(MD_atoms)) :: dU_ij ! derivative of three-body potentials (e.g. SW)


   Nat = size(MD_atoms) ! all atoms
   ! If there are three-body_potentials, set and save repeating functions for them:
   call get_repeated_three_body_potential(numpar, MD_atoms, MD_pots, U_ij, dU_ij) ! below

   ! For each pair of interacting atoms:
   !$omp parallel private (i, ion_done, j, atom_2, r, x, y, z, cos_a, cos_b, cos_c, R_ij, KOP1, KOP2, Force, Force_3bdy, Force_3bdy_tmp, &
   !$                      N_pot, i_pot, f_cut, d_f_cut, Pot_part, Pot_3bdy, d_Pot_part, tot_derive, k)
   !$omp do
   do i = 1, Nat
      ion_done = .false.    ! to start with
      MD_atoms(i)%U = 0.0d0         ! reset the potential to start with
      MD_atoms(i)%Force(:) = 0.0d0  ! reset the forces to start with
      MD_atoms(i)%A(:) = 0.0d0  ! reset the accelerations to start with
      ! Find nearest neighbors it is interacting with:
      do j = 1, numpar%Neighbors_Num(i)    ! that's how many atoms the atom "i" interacts with
         atom_2 => numpar%Neighbors_list(i, j)     ! that's index of the atom, which "i" is interacting with
         r => numpar%Neighbors_Rij(i, j)      ! [A] shortest distance between the atoms "i" and "atom_2" (which has number j in the list)
         x => numpar%Neighbors_Xij(i, j, 1)   ! [A] shortest distance between the atoms "i" and "atom_2" along X
         y => numpar%Neighbors_Xij(i, j, 2)   ! [A] shortest distance between the atoms "i" and "atom_2" along Y
         z => numpar%Neighbors_Xij(i, j, 3)   ! [A] shortest distance between the atoms "i" and "atom_2" along Z
         if (r < 1.0d-6) print*, 'Error in analytical_MD_pot_n_force: distance too small:', r, x, y, z
         ! Direction cosines between the two atoms (to calculate forces):
         cos_a = x/r
         cos_b = y/r
         cos_c = z/r
         R_ij = (/ x, y, z /)

         ! Find indices of the potential(s) is used for this pair of atoms:
         call find_which_potential(MD_atoms, i, atom_2, MD_pots, KOP1, KOP2)    ! module "MD_general_tools"

         ! How many potentials are used for this pair of atoms:
         Pot_part = 0.0d0  ! to start with
         Force = 0.0d0  ! to start summing up total forces
         Force_3bdy = 0.0d0   ! three-body contribution to the force
         N_pot = size(MD_pots(KOP1,KOP2)%Set)
         ! Find parameters of all potentials between the two atoms:
         do i_pot = 1, N_pot
            ! Define cut-off function and its derivative, to get total forces:
            f_cut = Fermi_cut_off(r, MD_pots(KOP1,KOP2)%Set(i_pot)%Par%d_cut, MD_pots(KOP1,KOP2)%Set(i_pot)%Par%dd) ! module "MD_general_tools"
            d_f_cut = d_Fermi_cut_off(r, MD_pots(KOP1,KOP2)%Set(i_pot)%Par%d_cut, MD_pots(KOP1,KOP2)%Set(i_pot)%Par%dd) ! module "MD_general_tools"

            ! Calculate the potential and its derivative by r:
            Force_3bdy_tmp = 0.0d0  ! to start with
            Pot_3bdy = 0.0d0  ! to start with
            call get_single_pot_and_force(MD_pots(KOP1,KOP2)%Set(i_pot)%Par, r, ion_done, Pot_part, d_Pot_part, &
                     numpar, i, atom_2, R_ij, MD_atoms, U_ij, dU_ij, Pot_3bdy, Force_3bdy_tmp, MD_pots, MD_supce)    ! below

            ! Combine the total potential from potential and cut-off function:
            MD_atoms(i)%U = MD_atoms(i)%U + Pot_part*f_cut + Pot_3bdy   ! [eV] add to total potential of this atom
            tot_derive = d_Pot_part * f_cut + Pot_part * d_f_cut    ! total derivative of the potential * cut-off function
            ! Get projection of forces (F=-dU(r_ij)/dr_ij):
            Force(1) = Force(1) - tot_derive * cos_a    ! [eV/A] add to total force along X
            Force(2) = Force(2) - tot_derive * cos_b    ! [eV/A] add to total force along Y
            Force(3) = Force(3) - tot_derive * cos_c    ! [eV/A] add to total force along Z
            ! Negative sign before the gradient:
            Force_3bdy(:) = Force_3bdy(:) - Force_3bdy_tmp(:)

!             write(*,'(a,i4,i4,f,f,f,f)'), 'a', i, atom_2, MD_atoms(i)%U, Force_3bdy(:)
!             write(*,'(a,i4,i4,f,f,f,f)'), 'b', i, atom_2, Pot_3bdy, Force_3bdy_tmp(:)

         enddo ! i_pot
         ! Add to the total force on atom i:
         MD_atoms(i)%Force(:) = MD_atoms(i)%Force(:) + Force(:) + Force_3bdy(:)

!          write(*,'(a,i4,i4,f,f,f,f)'), 'a', i, atom_2, MD_atoms(i)%U, MD_atoms(i)%Force(:)
!          write(*,'(a,i4,i4,f,f,f,f,f,f)'), 'b', i, atom_2, Force(:), Force_3bdy(:)

         ! Also save potential part of pressure in the supercell:
         do k = 1, 3
            MD_supce%P_pot(1,k) = MD_supce%P_pot(1,k) + (Force(k) + Force_3bdy(k))*x  ! [eV]
            MD_supce%P_pot(2,k) = MD_supce%P_pot(2,k) + (Force(k) + Force_3bdy(k))*y  ! [eV]
            MD_supce%P_pot(3,k) = MD_supce%P_pot(3,k) + (Force(k) + Force_3bdy(k))*z  ! [eV]
         enddo

      enddo ! j = numpar%Neighbors_Num(i)
   enddo ! i
   !$omp enddo
   !$omp end parallel
   nullify(atom_2, r, x, y, z)
!    call print_time_step('get_MD_pot_n_force 1:', 1.0, msec=.true.)   ! module "Little_subroutines"
end subroutine analytical_MD_pot_n_force




subroutine get_single_potential(MDPar, r, ion_done, Pot_part, &
            numpar, i, j, R_ij, U_ij, dU_ij, Pot_3bdy, MD_atoms, MD_supce, disp, ind_disp) ! short-range potentials
   class(MD_pot), intent(in) :: MDPar   ! MD potential parameters
   real(8), intent(in) :: r ! [A] interatomic distance
   logical, intent(inout) :: ion_done   ! index of ions, to calculate terms that are single-ion specific
   real(8), intent(out) :: Pot_part, Pot_3bdy ! potential and derivative of potential (without cut off function)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: i, j   ! indices of the interacting atoms
   real(8), dimension(3), intent(in) :: R_ij  ! [A] projections of the distance between atoms i and j
   real(8), dimension(:,:), intent(in) :: U_ij, dU_ij ! repeated contributions for three-body potentials (e.g. SW)
   type(Atom), dimension(:), intent(inout) :: MD_atoms	! all atoms in MD as objects
   type(MD_supcell), intent(inout) :: MD_supce  ! MD supercell parameters
!    class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   real(8), intent(in), optional :: disp  ! [A] displacement of atom i
   integer, intent(in), optional :: ind_disp  ! along which axis
   !------------------
   real(8) :: alpha
   select type(MDPar)
   type is (LJ) ! Lennard-Jones
      Pot_part = LJ_potential(MDPar%C12, MDPar%C6, r)   ! module "MD_Pot_Simple"

   type is (Buck) ! Buckingham potential
      Pot_part = Buck_potential(MDPar%A, MDPar%B, MDPar%C, r)   ! module "MD_Pot_Buck"

   type is (Matsui) ! Matsui potential
      Pot_part = Matsui_potential(MDPar%A, MDPar%B, MDPar%C, r)   ! module "MD_Pot_Buck"

   type is (SW) ! Stillinger-Weber potential
      if (present(disp) .and. present(ind_disp)) then
         call SW_potential(MDPar, r, i, j, R_ij, numpar, U_ij, dU_ij, MD_atoms, MD_supce, Pot_part, Pot_3bdy, disp, ind_disp)   ! module "MD_Pot_SW"
      else ! no displacement, regular coordinates:
         call SW_potential(MDPar, r, i, j, R_ij, numpar, U_ij, dU_ij, MD_atoms, MD_supce, Pot_part, Pot_3bdy)   ! module "MD_Pot_SW"
      endif

   type is (Coulomb)    ! Truncated Coulomb according to Wolf et al. [3]
      alpha = 3.0d0/(4.0d0*MDPar%d_cut) ! it's chosen according to optimal value from [4]
      Pot_part = Coulomb_Wolf_pot(MDPar%Z1, MDPar%Z2, r, MDPar%d_cut, alpha)   ! module "MD_Pot_Coulomb"
      ! Subtract self-term from the neutralizing part:
!      if (.not. ion_done) then
!          Pot_part = Pot_part - Coulomg_Wolf_self_term_simple(MDPar%Z1, MDPar%d_cut, alpha)     ! module "MD_Pot_Coulomb"
!          ion_done = .true. ! not to count the same ion again, mark it as done
!       endif

   type is (Coulomb_screened) ! Coulomb potential
      Pot_part = Bare_Coulomb_pot(MDPar%Z1, MDPar%Z2, r)   ! module "MD_Pot_Coulomb"

   type is (Power)   ! Power-law potential
      Pot_part = power_potential(MDPar%C, MDPar%n, r)   ! module "MD_Pot_Simple"

   type is (Exp_law)    ! Exponential potential
      Pot_part = exp_potential(MDPar%C, MDPar%k, r)   ! module "MD_Pot_Simple"

   type is (Soft_Coulomb)   ! Soft Coulomb
    ! not ready

   type is (Morse) ! Morse potential
    ! not ready

   type is (ZBL) ! ZBL
      ! not ready
      !Pot_part = ZBL_pot(MDPar%Z1, MDPar%Z2, r) ! module "ZBL"
      !d_Pot_part = d_ZBL_pot(MDPar%Z1, MDPar%Z2, r) ! module "ZBL"

   class default ! no potential
      Pot_part = 0.0d0
   endselect
end subroutine get_single_potential



subroutine get_single_pot_and_force(MDPar, r, ion_done, Pot_part, d_Pot_part, &
            numpar, i, j, R_ij, MD_atoms, U_ij, dU_ij, Pot_3bdy, Force_3bdy, MD_pots, MD_supce) ! short-range potentials
   class(MD_pot), intent(in) :: MDPar   ! MD potential parameters
   real(8), intent(in) :: r ! [A] interatomic distance
   logical, intent(inout) :: ion_done   ! index of ions, to calculate terms that are single-ion specific
   real(8), intent(out) :: Pot_part, d_Pot_part, Pot_3bdy ! potential and derivative of potential (without cut off function)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(in) :: i, j   ! indices of the interacting atoms
   real(8), dimension(3), intent(in) :: R_ij  ! [A] projections of the distance between atoms i and j
   type(Atom), dimension(:), intent(in) :: MD_atoms
   real(8), dimension(:,:), intent(in) :: U_ij, dU_ij ! repeated contributions for three-body potentials (e.g. SW)
   real(8), dimension(3), intent(out) :: Force_3bdy  ! three-body contribution into the force, if present
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   real(8) :: alpha
   select type(MDPar)
   type is (LJ) ! Lennard-Jones
      Pot_part = LJ_potential(MDPar%C12, MDPar%C6, r)   ! module "MD_Pot_Simple"
      d_Pot_part = d_LJ_potential(MDPar%C12, MDPar%C6, r)   ! module "MD_Pot_Simple"

   type is (Buck) ! Buckingham potential
      Pot_part = Buck_potential(MDPar%A, MDPar%B, MDPar%C, r)   ! module "MD_Pot_Buck"
      d_Pot_part = d_Buck_potential(MDPar%A, MDPar%B, MDPar%C, r)   ! module "MD_Pot_Buck"

   type is (Matsui) ! Matsui potential
      Pot_part = Matsui_potential(MDPar%A, MDPar%B, MDPar%C, r)   ! module "MD_Pot_Buck"
      d_Pot_part = d_Matsui_potential(MDPar%A, MDPar%B, MDPar%C, r)   ! module "MD_Pot_Buck"

   type is (SW) ! Stillinger-Weber potential
      !call SW_potential(MDPar, r, i, j, R_ij, numpar, U_ij, dU_ij, Pot_part, Pot_3bdy)   ! module "MD_Pot_SW"

      call d_SW_potential(MDPar, r, i, j, R_ij, numpar, U_ij, dU_ij, Pot_part, Pot_3bdy, &
                           d_Pot_part, Force_3bdy, MD_atoms, MD_pots, MD_supce)   ! module "MD_Pot_SW"

   type is (Coulomb)    ! Truncated Coulomb according to Wolf et al. [3]
      alpha = 3.0d0/(4.0d0*MDPar%d_cut) ! it's chosen according to optimal value from [4]
      Pot_part = Coulomb_Wolf_pot(MDPar%Z1, MDPar%Z2, r, MDPar%d_cut, alpha)   ! module "MD_Pot_Coulomb"
      d_Pot_part = d_Coulomb_Wolf_pot(MDPar%Z1, MDPar%Z2, r, MDPar%d_cut, alpha)   ! module "MD_Pot_Coulomb"
      ! Subtract self-term from the neutralizing part:
!      if (.not. ion_done) then
!          Pot_part = Pot_part - Coulomg_Wolf_self_term_simple(MDPar%Z1, MDPar%d_cut, alpha)     ! module "MD_Pot_Coulomb"
!          ion_done = .true. ! not to count the same ion again, mark it as done
!       endif

   type is (Coulomb_screened) ! Coulomb potential
      Pot_part = Bare_Coulomb_pot(MDPar%Z1, MDPar%Z2, r)   ! module "MD_Pot_Coulomb"
      d_Pot_part = d_Bare_Coulomb_pot(MDPar%Z1, MDPar%Z2, r)   ! module "MD_Pot_Coulomb"

   type is (Power)   ! Power-law potential
      Pot_part = power_potential(MDPar%C, MDPar%n, r)   ! module "MD_Pot_Simple"
      d_Pot_part = d_power_potential(MDPar%C, MDPar%n, r)   ! module "MD_Pot_Simple"

   type is (Exp_law)    ! Exponential potential
      Pot_part = exp_potential(MDPar%C, MDPar%k, r)   ! module "MD_Pot_Simple"
      d_Pot_part = d_exp_potential(MDPar%C, MDPar%k, r)   ! module "MD_Pot_Simple"

   type is (Soft_Coulomb)   ! Soft Coulomb
    ! not ready

   type is (Morse) ! Morse potential
    ! not ready

   type is (ZBL) ! ZBL
      ! not ready
      !Pot_part = ZBL_pot(MDPar%Z1, MDPar%Z2, r) ! module "ZBL"
      !d_Pot_part = d_ZBL_pot(MDPar%Z1, MDPar%Z2, r) ! module "ZBL"

   class default ! no potential
      Pot_part = 0.0d0
      d_Pot_part = 0.0d0
   endselect
end subroutine get_single_pot_and_force



subroutine get_repeated_three_body_part(numpar, MD_atoms, MD_pots, MD_supce, U_ij, atom_i, disp, ind_disp)
   type(Num_par), intent(in), target :: numpar   ! all numerical parameters
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(Atom), dimension(:), intent(in) :: MD_atoms
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   real(8), dimension(:,:), intent(inout) :: U_ij ! repeated contributions for three-body potentials (e.g. SW)
   integer, intent(in), optional :: atom_i   ! only for this atom
   real(8), intent(in), optional :: disp     ! [A] displacement
   integer, intent(in), optional :: ind_disp   ! along which axis is this displacement
   ! For the moment, only one three-body potential is allowed;
   ! If one wants to add more three-body potentials within the same simulation,
   ! U_ij must be made a 3-d array, including i_pot index for each three-body potential: U_ij(i,j,i_pot)
   !--------------------------------
   logical :: Three_body_pot_present
   integer :: i, j, i_pot, KOP1, KOP2, Nat, k
   real(8) :: r1, U, dU, f_cut, d_f_cut
   real(8) :: r
   integer, pointer :: atom_2

   Three_body_pot_present = .false.
   ! Find if there is any three-body potential:
   CHP:do i = 1, size(MD_pots,1)    ! check all pair of types of atoms
      do j = 1, size(MD_pots,2)
         do i_pot = 1, size(MD_pots(i,j)%Set) ! and all potentials for them
            ASSOCIATE (MDPot => MD_pots(i,j)%Set(i_pot)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(MDPot)
               type is (SW)   ! SW as three-body potential
                  Three_body_pot_present = .true.
                  exit CHP ! found it, no need to continue
               end select
            END ASSOCIATE
         enddo ! i_poy
      enddo ! j
   enddo CHP ! i

   ! If there is, construct the repeating part to be reused:
   if (Three_body_pot_present) then
      Nat = size(MD_atoms) ! all atoms

      if (present(atom_i)) then  ! only for this atom:
         ! Find nearest neighbors it is interacting with:
         !$omp parallel private (j, atom_2, r, i_pot, KOP1, KOP2, f_cut, U)
         !$omp do
         do j = 1, numpar%Neighbors_Num(atom_i)    ! that's how many atoms the atom "i" interacts with
            atom_2 => numpar%Neighbors_list(atom_i, j)     ! that's index of the atom, which "i" is interacting with

            if (present(disp) .and. present(ind_disp)) then ! displaced atom:
               call shortest_distance(atom_i, atom_2, MD_atoms, MD_supce, r, disp = disp, ind_disp = ind_disp) ! module "MD_general_tools"
            else  ! atom in its default position:
               call shortest_distance(atom_i, atom_2, MD_atoms, MD_supce, r) ! module "MD_general_tools"
            endif

            ! Find indices of the potential(s) is used for this pair of atoms:
            call find_which_potential(MD_atoms, atom_i, atom_2, MD_pots, KOP1, KOP2)    ! module "MD_general_tools"

            do i_pot = 1, size(MD_pots(KOP1, KOP2)%Set) ! and all potentials for them
               ! Get the function:
               ASSOCIATE (MDPot => MD_pots(KOP1,KOP2)%Set(i_pot)%Par)
                  select type(MDPot)
                  type is (SW)   ! SW is three-body potential:
                     ! Cut off functions:
                     f_cut = Fermi_cut_off(r, MDPot%d_cut, MDPot%dd) ! module "MD_general_tools"
                     ! Radial functions:
                     U = SW_U_ij(MDPot, r) ! module "MD_Pot_SW"
                     ! Functions including cut off:
                     U_ij(atom_i, atom_2) = U * f_cut

                  end select
               ENDASSOCIATE

            enddo ! i_pot
         enddo ! j
         !$omp enddo
         !$omp end parallel

         ! Symmetrize the interaction:
         do j = 1, numpar%Neighbors_Num(atom_i)    ! that's how many atoms the atom "i" interacts with
            atom_2 => numpar%Neighbors_list(atom_i, j)     ! that's index of the atom, which "i" is interacting with
            U_ij(atom_2, atom_i) = U_ij(atom_i, atom_2)
         enddo
      else  ! for all atoms:
         U_ij = 0.0d0      ! reset

         ! Set it:
         !$omp parallel private (i, j, atom_2, r, i_pot, KOP1, KOP2, f_cut, d_f_cut, U, dU)
         !$omp do
         do i = 1, Nat ! for all atoms
            ! Find nearest neighbors it is interacting with:
            do j = 1, numpar%Neighbors_Num(i)    ! that's how many atoms the atom "i" interacts with
               atom_2 => numpar%Neighbors_list(i, j)     ! that's index of the atom, which "i" is interacting with
!                call shortest_distance(i, atom_2, MD_atoms, MD_supce, r) ! module "MD_general_tools"
               if (present(disp) .and. present(ind_disp)) then ! displaced atom:
                  call shortest_distance(i, atom_2, MD_atoms, MD_supce, r, disp = disp, ind_disp = ind_disp) ! module "MD_general_tools"
               else  ! atom in its default position:
                  call shortest_distance(i, atom_2, MD_atoms, MD_supce, r) ! module "MD_general_tools"
               endif

               ! Find indices of the potential(s) is used for this pair of atoms:
               call find_which_potential(MD_atoms, i, atom_2, MD_pots, KOP1, KOP2)    ! module "MD_general_tools"

               do i_pot = 1, size(MD_pots(KOP1, KOP2)%Set) ! and all potentials for them
                  ! Get the function:
                  ASSOCIATE (MDPot => MD_pots(KOP1,KOP2)%Set(i_pot)%Par)
                     select type(MDPot)
                     type is (SW)   ! SW is three-body potential:
                        ! Cut off functions:
                        f_cut = Fermi_cut_off(r, MDPot%d_cut, MDPot%dd) ! module "MD_general_tools"
                        ! Radial functions:
                        U = SW_U_ij(MDPot, r) ! module "MD_Pot_SW"
                        ! Functions including cut off:
                        U_ij(i, atom_2) = U * f_cut

                     end select
                  ENDASSOCIATE

               enddo ! i_pot
            enddo ! j
         enddo ! i
         !$omp enddo
         !$omp end parallel
      endif ! (present(atom_i))
   else
      if (present(atom_i)) then
         U_ij(atom_i,:) = 0.0d0
      else
         U_ij(:,:) = 0.0d0
      endif
   endif

   nullify(atom_2)
end subroutine get_repeated_three_body_part



subroutine get_repeated_three_body_potential(numpar, MD_atoms, MD_pots, U_ij, dU_ij)
   type(Num_par), intent(in), target :: numpar   ! all numerical parameters
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(Atom), dimension(:), intent(in) :: MD_atoms
   real(8), dimension(:,:), intent(inout) :: U_ij ! repeated contributions for three-body potentials (e.g. SW)
   real(8), dimension(:,:), intent(inout) :: dU_ij ! repeated contributions for derivative of three-body potentials
   ! For the moment, only one three-body potential is allowed;
   ! If one wants to add more three-body potentials within the same simulation,
   ! U_ij must be made a 3-d array, including i_pot index for each three-body potential: U_ij(i,j,i_pot)
   !--------------------------------
   logical :: Three_body_pot_present
   integer :: i, j, i_pot, KOP1, KOP2, Nat, k
   real(8) :: r1, U, dU, f_cut, d_f_cut
   real(8), pointer :: r
   integer, pointer :: atom_2

   Three_body_pot_present = .false.
   ! Find if there is any three-body potential:
   CHP:do i = 1, size(MD_pots,1)    ! check all pair of types of atoms
      do j = 1, size(MD_pots,2)
         do i_pot = 1, size(MD_pots(i,j)%Set) ! and all potentials for them
            ASSOCIATE (MDPot => MD_pots(i,j)%Set(i_pot)%Par) ! this is the sintax we have to use to check the class of defined types
               select type(MDPot)
               type is (SW)   ! SW as three-body potential
                  Three_body_pot_present = .true.
                  exit CHP ! found it, no need to continue
               end select
            END ASSOCIATE
         enddo ! i_poy
      enddo ! j
   enddo CHP ! i

   ! If there is, construct the repeating part to be reused:
   if (Three_body_pot_present) then
      Nat = size(MD_atoms) ! all atoms
      U_ij = 0.0d0      ! reset
      dU_ij = 0.0d0     ! reset

      ! Set it:
      !$omp parallel private (i, j, atom_2, r, i_pot, KOP1, KOP2, f_cut, d_f_cut, U, dU)
      !$omp do
      do i = 1, Nat ! for all atoms
         ! Find nearest neighbors it is interacting with:
         do j = 1, numpar%Neighbors_Num(i)    ! that's how many atoms the atom "i" interacts with
            atom_2 => numpar%Neighbors_list(i, j)     ! that's index of the atom, which "i" is interacting with
            r => numpar%Neighbors_Rij(i, j)      ! [A] shortest distance between the atoms "i" and "atom_2" (which has number j in the list)

            ! Find indices of the potential(s) is used for this pair of atoms:
            call find_which_potential(MD_atoms, i, atom_2, MD_pots, KOP1, KOP2)    ! module "MD_general_tools"

            do i_pot = 1, size(MD_pots(KOP1, KOP2)%Set) ! and all potentials for them
               ! Get the function:
               ASSOCIATE (MDPot => MD_pots(KOP1,KOP2)%Set(i_pot)%Par)
                  select type(MDPot)
                  type is (SW)   ! SW is three-body potential:
                     ! Cut off functions:
                     f_cut = Fermi_cut_off(r, MDPot%d_cut, MDPot%dd) ! module "MD_general_tools"
                     d_f_cut = d_Fermi_cut_off(r, MDPot%d_cut, MDPot%dd) ! module "MD_general_tools"
                     ! Radial functions:
                     U = SW_U_ij(MDPot, r) ! module "MD_Pot_SW"
                     dU = SW_dU_ij(MDPot, r) ! module "MD_Pot_SW"
                     ! Functions including cut off:
                     U_ij(i, atom_2) = U * f_cut
                     dU_ij(i, atom_2) = dU * f_cut + U * d_f_cut

!                      r1 = 1.0d0
!                      do k = 1, 100
!                         r1 = r1 + 4.0d0/100.0d0
!                         print*, k, r1, SW_U_ij(MDPot, r1), SW_dU_ij(MDPot, r1)
!                      enddo
!                      pause 'get_repeated_three_body_potential'

                  end select
               ENDASSOCIATE

            enddo ! i_pot
         enddo ! j
      enddo ! i
      !$omp enddo
      !$omp end parallel
   endif

   nullify(r, atom_2)
end subroutine get_repeated_three_body_potential



end module MD
