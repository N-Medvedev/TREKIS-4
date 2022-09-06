! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev 
! in 2021
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains routines for analysis of molecular dynamics data:
module MD_data_analysis
use Universal_constants
use Objects
use Little_subroutines, only: Find_in_array_monoton
use MD_general_tools, only: get_total_energies, get_mean_displacements

implicit none


 contains


subroutine analyze_MD_output_data(numpar, MD_supce, MD_pots, MD_atoms, out_data)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(MD_supcell), intent(in) :: MD_supce  ! MD supercell parameters
   type(MD_potential), dimension(:,:), intent(in) :: MD_pots    ! MD potentials for each kind of atom-atom interactions
   type(Atom), dimension(:), intent(inout), allocatable :: MD_atoms ! all atoms in MD as objects
   type(output_data), intent(inout) :: out_data  ! all output data (distributions etc.)
   !--------------------------
   
   if (numpar%DO_MD) then   ! only if MD module is active it makes sense to analyze the data
      ! Get energies in MD:
      call get_total_energies(MD_atoms, out_data%MD_Ekin, out_data%MD_Epot)    ! module "MD_general_tools"

      ! Get mean atomic displacements:
      call get_mean_displacements(MD_supce, MD_pots, MD_atoms, MD_supce%MD_atoms_0, dble(numpar%n_MSD), &
                  out_data%MD_MSD, out_data%MD_MSDP)    ! module "MD_general_tools"

      ! Get pressure:
      ! NOT READY
    
      ! Get distribution of parameters in space on a defined grid:
      ! NOT READY
   else     ! no MD, no data:
      out_data%MD_Ekin = 0.0d0
      out_data%MD_Epot = 0.0d0
   endif
end subroutine analyze_MD_output_data



 
end module MD_data_analysis
