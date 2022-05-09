! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018-2020
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains all global variables

MODULE Variables
use Objects

 implicit none

type(Error_handling) :: g_Err     ! error log
type(Matter) :: g_target            ! parameters of the target
type(Num_par) :: g_numpar      ! numerical parameters
type(Radiation_param), dimension(:), allocatable :: g_bunch     ! incomming radiation
type(MC_arrays), dimension(:), allocatable :: g_MC  ! all MC arrays for particles: photons, electrons and holes; size equals to number of iterations
type(output_data) :: g_output   ! all output data (distributions etc.)
type(Atom), dimension(:), allocatable :: g_MD_atoms     ! all atoms in MD as objects
type(MD_supcell) :: g_MD_supce  ! MD supercell parameters
type(MD_potential), dimension(:,:), allocatable :: g_MD_pots    ! MD potentials for each kind of atom-atom interactions

real(8) :: g_time         ! running time
real(8) :: g_dt_save    ! variable for tracing when to print output files
real(8) :: g_dt_out      ! setting next printout step
integer :: g_i_steps    ! number of timesteps of the main program

integer :: g_ctim(8), g_c1(8)   ! time stamps
real(8) :: g_as1                ! duration of execution
character(200) :: g_text    ! general text variable to be used in the main program

 contains
 
 
END MODULE Variables
