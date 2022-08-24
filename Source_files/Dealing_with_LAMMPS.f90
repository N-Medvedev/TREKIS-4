! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2022
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with files in LAMMPS format:
! [1] https://docs.lammps.org/2001/data_format.html
! [2] https://docs.lammps.org/units.html

MODULE Dealing_with_LAMMPS

use Universal_constants
use Objects
use Dealing_with_files
use MD_general_tools, only: find_which_potential

implicit none

 contains


subroutine Write_LAMMPS_input_file(FN, used_target, numpar, MD_supce, MD_atoms, MD_pots, textline)
   integer, intent(in) :: FN    ! file number to write into
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(Num_par), intent(in) :: numpar    ! all numerical parameters
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   type(Atom), dimension(:), intent(inout), allocatable :: MD_atoms ! all atoms in MD as objects
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   character(*), intent(in) :: textline ! first line of LAMMPS input file
   !---------------------------
   integer :: Noft, NofE, N_types, i, j, N_pot, i_pot
   character(20) :: text_var, atomic_type
   logical :: chargedpot, pot_keyword

   ! Write LAMMPS input file:
   ! 1) 1st comment line, text:
   write(FN, '(a)') trim(adjustl(textline))
   ! 2d comment line, use to save info about types and styles:
   ! 1.a) style of atoms (with or without charge)
   chargedpot = .false. ! to start with
   POTS:do i = 1, size(MD_pots,1)
      N_pot = size(MD_pots(i,1)%Set)
      ! Find parameters of all potentials between the two atoms:
      do i_pot = 1, N_pot
         call get_pot_type_keyword(MD_pots(i,1)%Set(i_pot)%Par, pot_keyword)    ! below
         if ( pot_keyword ) then  ! function below
            chargedpot = .true.
            exit POTS   ! we found charged potential, no need to check more
         endif
      enddo
   enddo POTS
   if (chargedpot) then
      text_var = 'charge'
   else
      text_var = 'atomic'
   endif
   atomic_type = text_var   ! save for later
   write(FN, '(a)', advance='no') '# units_style '//trim(adjustl(numpar%LAMMPS_UNITS))// &
                                  ', atom_style '//trim(adjustl(text_var))//', atom types: '
   ! 1.b) Printout types of atoms for reference:
   Noft = size(used_target%Material)    ! number of different target materials
   NofE = 0 ! to start with
   do i = 1, Noft   ! for all targets
      do j = 1, size(used_target%Material(i)%Elements)
         NofE = NofE + 1
         write(text_var, '(i4)') NofE   ! number of atom types
         write(FN, '(a)', advance='no') trim(adjustl(text_var))//'='//trim(adjustl(used_target%Material(i)%Elements(j)%Name))//' '
      enddo
   enddo
   write(FN, '(a)') ! skipline


   ! 2) Number of atoms:
   write(text_var, '(i16)') size(MD_atoms)  ! number of atoms
   write(FN, '(a)') trim(adjustl(text_var))//' atoms'

   ! 3) Number of types atoms (elements):
   ! count diferent types of atoms:
   write(text_var, '(i16)') size(MD_pots,1)  ! number of atom types
   write(FN, '(a)') trim(adjustl(text_var))//' atom types'
   write(FN, '(a)') ! skipline

   ! 4) supercell coordinates [A]:
   write(FN, '(f,f,a)') MD_supce%x_start, MD_supce%x_end, ' xlo xhi'
   write(FN, '(f,f,a)') MD_supce%y_start, MD_supce%y_end, ' ylo yhi'
   write(FN, '(f,f,a)') MD_supce%z_start, MD_supce%z_end, ' zlo zhi'
   write(FN, '(a)') ! skipline

   ! 5) masses:
   write(FN, '(a)') 'Masses'
   write(FN, '(a)') ! skipline
   Noft = size(used_target%Material)    ! number of different target materials
   NofE = 0 ! to start with
   do i = 1, Noft   ! for all targets
      do j = 1, size(used_target%Material(i)%Elements)
         NofE = NofE + 1
         write(text_var, '(i4)') NofE   ! number of atom types
         write(FN, '(a,f)') trim(adjustl(text_var)), used_target%Material(i)%Elements(j)%Mass
      enddo
   enddo
   write(FN, '(a)') ! skipline

   ! 6) atomic coordinates:
   write(FN, '(a)') 'Atoms # '//trim(adjustl(atomic_type))
   write(FN, '(a)') ! skipline
   call write_atomic_coords_LAMMPS(FN, MD_atoms, MD_pots, chargedpot, numpar%LAMMPS_UNITS)  ! below
   write(FN, '(a)') ! skipline

   ! 7) atomic velocities:
   write(FN, '(a)') 'Velocities'
   write(FN, '(a)') ! skipline
   call write_atomic_velocity_LAMMPS(FN, MD_atoms, numpar%LAMMPS_UNITS) ! below

end subroutine Write_LAMMPS_input_file



subroutine write_atomic_coords_LAMMPS(FN, MD_atoms, MD_pots, chargedpot, LAMMPS_UNITS)
   integer, intent(in) :: FN    ! file number to write into
   type(Atom), dimension(:), intent(inout), allocatable :: MD_atoms ! all atoms in MD as objects
   class(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   logical, intent(in) :: chargedpot    ! atomic charge needs to be written or not
   character(*), intent(in) :: LAMMPS_UNITS  ! which LAMMPS units to use
   !-------------------------
   real(8) :: units_coeff, Z1
   integer :: i, Nat, N_pot, i_pot, KOP1, KOP2
   logical :: pot_keyword

   ! Find conversion coefficients into the units used in LAMMPS and defined by the user:
   select case(LAMMPS_UNITS) ! See description in [2]
   case ('real')
      units_coeff = 1.0d0   ! [A] -> [A]
   case ('metal')
      units_coeff = 1.0d0   ! [A] -> [A]
   case ('lj')
      units_coeff = 1.0d0   ! [A] -> [A]
   case ('si')
      units_coeff = 1.0d-10 ! [A] -> [m]
   case ('cgs')
      units_coeff = 1.0d-8  ! [A] -> [cm]
   case ('electron')
      units_coeff = g_a0    ! [A] -> [Bohr]
   case ('micro')
      units_coeff = 1.0d-4  ! [A] -> [mkm]
   case ('nano')
      units_coeff = 1.0d-1  ! [A] -> [nm]
   case default
      units_coeff = 1.0d0  ! [A]
   endselect


   ! Number of atoms:
   Nat = size(MD_atoms)

   ! Defined format: with or without charge:
   if (chargedpot) then ! with charge
      do i = 1, Nat ! for all atoms
         ! Find type of atom:
         call find_which_potential(MD_atoms, i, 1, MD_pots, KOP1, KOP2)   ! module "MD_general_tools"
         N_pot = size(MD_pots(KOP1,KOP2)%Set)
         ! Find parameters of all potentials between the two atoms:
         do i_pot = 1, N_pot    ! Find potential indices that contain charge
            call get_pot_type_keyword(MD_pots(KOP1,KOP2)%Set(i_pot)%Par, pot_keyword, Z1)    ! below
            if (pot_keyword) exit   ! function below
         enddo
         ! Write LAMMPS coordinates:
         write(FN,'(i,i,es,es,es,es)') i, KOP1, Z1, units_coeff*MD_atoms(i)%R(:)
      enddo
   else ! without charge
      do i = 1, Nat ! for all atoms
         ! Find type of atom:
         call find_which_potential(MD_atoms, i, 1, MD_pots, KOP1, KOP2)   ! module "MD_general_tools"
         ! Write LAMMPS coordinates:
         write(FN,'(i,i,es,es,es)') i, KOP1, units_coeff*MD_atoms(i)%R(:)
      enddo
   endif

end subroutine write_atomic_coords_LAMMPS


subroutine write_atomic_velocity_LAMMPS(FN, MD_atoms, LAMMPS_UNITS)
   integer, intent(in) :: FN    ! file number to write into
   type(Atom), dimension(:), intent(inout), allocatable :: MD_atoms ! all atoms in MD as objects
   character(*), intent(in) :: LAMMPS_UNITS  ! which LAMMPS units to use
   !-------------------------
   real(8) :: units_coeff
   integer :: i, Nat

   ! Find conversion coefficients into the units used in LAMMPS and defined by the user:
   select case(LAMMPS_UNITS)    ! See description in [2]
   case ('real')
      units_coeff = 1.0d0   ! [A/fs] -> [A/fs]
   case ('metal')
      units_coeff = 1.0d3   !1.0d0/1.0d-3   ! [A/fs] -> [A/ps]
   case ('lj')
      units_coeff = 1.0d0   ! [A/fs] -> [A/fs]
   case ('si')
      units_coeff = g_Afs2ms   !1.0d-10/1.0d-15 ! [A/fs] -> [m/s]
   case ('cgs')
      units_coeff = 1.0d7   !1.0d-8/1.0d-15  ! [A/fs] -> [cm/s]
   case ('electron')
      units_coeff = g_a0/1.03275d0    ! [A/fs] -> [Bohr/a.t.u.] [1.03275e-15 seconds] [2]
   case ('micro')
      units_coeff = 1.0d-4/1.0d-9  ! [A/fs] -> [um/us]
   case ('nano')
      units_coeff = 1.0d-1/1.0d-6  ! [A/fs] -> [nm/ns]
   case default
      units_coeff = 1.0d0  ! [A/fs]
   endselect

   ! Number of atoms:
   Nat = size(MD_atoms)
   do i = 1, Nat ! for all atoms
      write(FN,'(i,es,es,es)') i, units_coeff*MD_atoms(i)%V(:)
   enddo
end subroutine write_atomic_velocity_LAMMPS



subroutine get_pot_type_keyword(MDPar, chargedpot, Z1)
   logical, intent(out) :: chargedpot
   class(MD_pot), intent(in) :: MDPar   ! MD potential parameters
   real(8), intent(out), optional :: Z1 ! atomic charge
   select type(MDPar)
   type is (LJ) ! Lennard-Jones
      !potkeyword = 'atomic' ! no charges in this potential
      chargedpot = .false.
      if (present(Z1)) Z1 = 0.0d0
   type is (Buck) ! Buckingham potential
      !potkeyword = 'atomic' ! no charges in this potential
      chargedpot = .false.
      if (present(Z1)) Z1 = 0.0d0
   type is (Matsui) ! Matsui potential
      !potkeyword = 'atomic' ! no charges in this potential
      chargedpot = .false.
      if (present(Z1)) Z1 = 0.0d0
   type is (SW) ! Stillinger-Weber potential
      !potkeyword = 'atomic' ! no charges in this potential
      chargedpot = .false.
      if (present(Z1)) Z1 = 0.0d0
   type is (Coulomb)    ! Truncated Coulomb according to Wolf et al. [3]
      !potkeyword = 'charge' ! charges in this potential
      chargedpot = .true.
      if (present(Z1)) Z1 = MDPar%Z1
   type is (Coulomb_screened) ! Coulomb potential
      !potkeyword = 'charge' ! charges in this potential
      chargedpot = .true.
      if (present(Z1)) Z1 = MDPar%Z1
   type is (Power)   ! Power-law potential
      !potkeyword = 'atomic' ! no charges in this potential
      chargedpot = .false.
      if (present(Z1)) Z1 = 0.0d0
   type is (Exp_law)    ! Exponential potential
      !potkeyword = 'atomic' ! no charges in this potential
      chargedpot = .false.
      if (present(Z1)) Z1 = 0.0d0
   type is (Soft_Coulomb)   ! Soft Coulomb
    ! not ready
   type is (Morse) ! Morse potential
    ! not ready
   type is (ZBL) ! ZBL
      ! not ready
   class default ! no pootential - no charge
      !potkeyword = 'atomic' ! no charges in this potential
      chargedpot = .false.
      if (present(Z1)) Z1 = 0.0d0
   endselect
end subroutine get_pot_type_keyword




END MODULE Dealing_with_LAMMPS
