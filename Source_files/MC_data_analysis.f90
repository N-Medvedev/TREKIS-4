! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev 
! and R. Rymzhanov
! in 2020-2021
! 1111111111111111111111111111111111111111111111111111111111111
! References:
! [1] https://www.cs.cornell.edu/courses/cs6630/2015fa/notes/pdf-transform.pdf
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains routines for analysis of Monte Carlo data:
module MC_data_analysis
use Universal_constants
use Objects
use Little_subroutines, only: Find_in_array_monoton
use MC_general_tools, only: renew_atomic_arrays, place_particles_back_into_box, put_back_into_box
use Geometries, only: get_v_theta, m_tollerance_eps

implicit none


real(8) :: m_one_over_halfPi
! Used for testing purposes:
parameter (m_one_over_halfPi = 2.0d0/g_Pi)    ! normalization factor for spherical theta


 contains


!=======================================
! Output data definition:

subroutine analyze_MC_output_data(used_target, numpar, MC, out_data, tim)
   type(Matter), intent(in) :: used_target  ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(MC_arrays), dimension(:), intent(inout) :: MC	! all MC arrays for all particles; size = number of iterations
   type(output_data), intent(inout) :: out_data  ! all output data (distributions etc.)
   real(8), intent(in) :: tim   ! [fs] current time step
   !-----------------------------------------------
   ! Define temporary arrays to average the MC data into distributions:
   ! Energy disributions (spectra):
   real(8), dimension(:), allocatable :: Spectrum_ph, Spectrum_e, Spectrum_p, Spectrum_SHI, Spectrum_mu  ! energy spectra
   real(8), dimension(:,:), allocatable :: Spectrum_h ! VB spectra separate for each target
   ! Velosity theta disributions:
   real(8), dimension(:), allocatable :: Vel_theta_ph, Vel_theta_e, Vel_theta_p, Vel_theta_h, Vel_theta_SHI, Vel_theta_mu  ! velosity theta distributions
   ! Energy spectra vs space in 1d:
   real(8), dimension(:,:), allocatable :: Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X, Spectra_mu_X  ! energy spectra in space along X
   real(8), dimension(:,:), allocatable :: Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y, Spectra_mu_Y  ! energy spectra in space along Y
   real(8), dimension(:,:), allocatable :: Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z, Spectra_mu_Z  ! energy spectra in space along Z
   real(8), dimension(:,:), allocatable :: Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R, Spectra_mu_R  ! energy spectra in space along R
   ! Theta distribution vs space in 1d:
   real(8), dimension(:,:), allocatable :: Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X, Theta_mu_X  ! theta distribution in space along X
   real(8), dimension(:,:), allocatable :: Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y, Theta_mu_Y  ! theta distribution in space along Y
   real(8), dimension(:,:), allocatable :: Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z, Theta_mu_Z  ! theta distribution in space along Z
   real(8), dimension(:,:), allocatable :: Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R, Theta_mu_R  ! theta distribution in space along R
   ! Spatial distributions in 1d:
   real(8), dimension(:), allocatable :: Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_ph_R, Distr_ph_L, Distr_ph_Theta, Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic ! photon
   real(8), dimension(:), allocatable :: Distr_e_X, Distr_e_Y, Distr_e_Z, Distr_e_R, Distr_e_L, Distr_e_Theta, Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic ! electron
   real(8), dimension(:), allocatable :: Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_p_R, Distr_p_L, Distr_p_Theta, Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic ! positron
   real(8), dimension(:,:), allocatable :: Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_h_R, Distr_h_L, Distr_h_Theta, Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic ! hole
   real(8), dimension(:), allocatable :: Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, Distr_SHI_R, Distr_SHI_L, Distr_SHI_Theta, Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic ! SHI
   real(8), dimension(:), allocatable :: Distr_a_X, Distr_a_Y, Distr_a_Z, Distr_a_R, Distr_a_L, Distr_a_Theta, Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic ! Atom
   real(8), dimension(:), allocatable :: Distr_mu_X, Distr_mu_Y, Distr_mu_Z, Distr_mu_R, Distr_mu_L, Distr_mu_Theta, Distr_mu_Rc, Distr_mu_Thetac, Distr_mu_Phic ! muon
   ! Spatial distributions in 2d:
   real(8), dimension(:,:), allocatable :: Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ, Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta
   real(8), dimension(:,:), allocatable :: Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic  ! photon
   real(8), dimension(:,:), allocatable :: Distr_e_XY, Distr_e_YZ, Distr_e_XZ, Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta
   real(8), dimension(:,:), allocatable :: Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic  ! electron
   real(8), dimension(:,:), allocatable :: Distr_p_XY, Distr_p_YZ, Distr_p_XZ, Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta
   real(8), dimension(:,:), allocatable :: Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic  ! positron
   real(8), dimension(:,:,:), allocatable :: Distr_h_XY, Distr_h_YZ, Distr_h_XZ, Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta
   real(8), dimension(:,:,:), allocatable :: Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic  ! hole
   real(8), dimension(:,:), allocatable :: Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ, Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta
   real(8), dimension(:,:), allocatable :: Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic  ! SHI
   real(8), dimension(:,:), allocatable :: Distr_a_XY, Distr_a_YZ, Distr_a_XZ, Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta
   real(8), dimension(:,:), allocatable :: Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic  ! Atom
   real(8), dimension(:,:), allocatable :: Distr_mu_XY, Distr_mu_YZ, Distr_mu_XZ, Distr_mu_RL, Distr_mu_RTheta, Distr_mu_LTheta, &
                                           Distr_mu_RcThc, Distr_mu_RcPhic, Distr_mu_ThcPhic  ! muon
   ! Spatial distributions in 3d:
   real(8), dimension(:,:,:), allocatable :: Distr_ph_XYZ, Distr_ph_RLTheta, Distr_ph_RcThcPhic ! photon
   real(8), dimension(:,:,:), allocatable :: Distr_e_XYZ, Distr_e_RLTheta, Distr_e_RcThcPhic ! electron
   real(8), dimension(:,:,:), allocatable :: Distr_p_XYZ, Distr_p_RLTheta, Distr_p_RcThcPhic ! positron
   real(8), dimension(:,:,:,:), allocatable :: Distr_h_XYZ, Distr_h_RLTheta, Distr_h_RcThcPhic ! hole
   real(8), dimension(:,:,:), allocatable :: Distr_SHI_XYZ, Distr_SHI_RLTheta, Distr_SHI_RcThcPhic ! SHI
   real(8), dimension(:,:,:), allocatable :: Distr_a_XYZ, Distr_a_RLTheta, Distr_a_RcThcPhic ! Atom
   real(8), dimension(:,:,:), allocatable :: Distr_mu_XYZ, Distr_mu_RLTheta, Distr_mu_RcThcPhic ! muon
   ! Spatial energy distributions in 1d:
   real(8), dimension(:), allocatable :: E_Distr_ph_X, E_Distr_ph_Y, E_Distr_ph_Z, E_Distr_ph_R, E_Distr_ph_L, E_Distr_ph_Theta
   real(8), dimension(:), allocatable :: E_Distr_ph_Rc, E_Distr_ph_Thetac, E_Distr_ph_Phic ! photon
   real(8), dimension(:), allocatable :: E_Distr_e_X, E_Distr_e_Y, E_Distr_e_Z, E_Distr_e_R, E_Distr_e_L, E_Distr_e_Theta
   real(8), dimension(:), allocatable :: E_Distr_e_Rc, E_Distr_e_Thetac, E_Distr_e_Phic ! electron
   real(8), dimension(:), allocatable :: E_Distr_p_X, E_Distr_p_Y, E_Distr_p_Z, E_Distr_p_R, E_Distr_p_L, E_Distr_p_Theta
   real(8), dimension(:), allocatable :: E_Distr_p_Rc, E_Distr_p_Thetac, E_Distr_p_Phic ! positron
   real(8), dimension(:,:), allocatable :: E_Distr_h_X, E_Distr_h_Y, E_Distr_h_Z, E_Distr_h_R, E_Distr_h_L, E_Distr_h_Theta
   real(8), dimension(:,:), allocatable :: E_Distr_h_Rc, E_Distr_h_Thetac, E_Distr_h_Phic ! hole
   real(8), dimension(:), allocatable :: E_Distr_SHI_X, E_Distr_SHI_Y, E_Distr_SHI_Z, E_Distr_SHI_R, E_Distr_SHI_L
   real(8), dimension(:), allocatable :: E_Distr_SHI_Theta, E_Distr_SHI_Rc, E_Distr_SHI_Thetac, E_Distr_SHI_Phic ! SHI
   real(8), dimension(:), allocatable :: E_Distr_a_X, E_Distr_a_Y, E_Distr_a_Z, E_Distr_a_R, E_Distr_a_L
   real(8), dimension(:), allocatable :: E_Distr_a_Theta, E_Distr_a_Rc, E_Distr_a_Thetac, E_Distr_a_Phic ! Atom
   real(8), dimension(:), allocatable :: E_Distr_mu_X, E_Distr_mu_Y, E_Distr_mu_Z, E_Distr_mu_R, E_Distr_mu_L, E_Distr_mu_Theta, &
                                         E_Distr_mu_Rc, E_Distr_mu_Thetac, E_Distr_mu_Phic ! muon
   ! Spatial energy distributions in 2d:
   real(8), dimension(:,:), allocatable :: E_Distr_ph_XY, E_Distr_ph_YZ, E_Distr_ph_XZ, E_Distr_ph_RL, E_Distr_ph_RTheta, E_Distr_ph_LTheta
   real(8), dimension(:,:), allocatable :: E_Distr_ph_RcThc, E_Distr_ph_RcPhic, E_Distr_ph_ThcPhic  ! photon
   real(8), dimension(:,:), allocatable :: E_Distr_e_XY, E_Distr_e_YZ, E_Distr_e_XZ, E_Distr_e_RL, E_Distr_e_RTheta, E_Distr_e_LTheta
   real(8), dimension(:,:), allocatable :: E_Distr_e_RcThc, E_Distr_e_RcPhic, E_Distr_e_ThcPhic  ! electron
   real(8), dimension(:,:), allocatable :: E_Distr_p_XY, E_Distr_p_YZ, E_Distr_p_XZ, E_Distr_p_RL, E_Distr_p_RTheta, E_Distr_p_LTheta
   real(8), dimension(:,:), allocatable :: E_Distr_p_RcThc, E_Distr_p_RcPhic, E_Distr_p_ThcPhic  ! positron
   real(8), dimension(:,:,:), allocatable :: E_Distr_h_XY, E_Distr_h_YZ, E_Distr_h_XZ, E_Distr_h_RL, E_Distr_h_RTheta, E_Distr_h_LTheta
   real(8), dimension(:,:,:), allocatable :: E_Distr_h_RcThc, E_Distr_h_RcPhic, E_Distr_h_ThcPhic  ! hole
   real(8), dimension(:,:), allocatable :: E_Distr_SHI_XY, E_Distr_SHI_YZ, E_Distr_SHI_XZ, E_Distr_SHI_RL, E_Distr_SHI_RTheta, E_Distr_SHI_LTheta
   real(8), dimension(:,:), allocatable :: E_Distr_SHI_RcThc, E_Distr_SHI_RcPhic, E_Distr_SHI_ThcPhic  ! SHI
   real(8), dimension(:,:), allocatable :: E_Distr_a_XY, E_Distr_a_YZ, E_Distr_a_XZ, E_Distr_a_RL, E_Distr_a_RTheta, E_Distr_a_LTheta
   real(8), dimension(:,:), allocatable :: E_Distr_a_RcThc, E_Distr_a_RcPhic, E_Distr_a_ThcPhic  ! Atom
   real(8), dimension(:,:), allocatable :: E_Distr_mu_XY, E_Distr_mu_YZ, E_Distr_mu_XZ, E_Distr_mu_RL, E_Distr_mu_RTheta, E_Distr_mu_LTheta, &
                                           E_Distr_mu_RcThc, E_Distr_mu_RcPhic, E_Distr_mu_ThcPhic  ! muon
   ! Spatial energy distributions in 3d:
   real(8), dimension(:,:,:), allocatable :: E_Distr_ph_XYZ, E_Distr_ph_RLTheta, E_Distr_ph_RcThcPhic ! photon
   real(8), dimension(:,:,:), allocatable :: E_Distr_e_XYZ, E_Distr_e_RLTheta, E_Distr_e_RcThcPhic ! electron
   real(8), dimension(:,:,:), allocatable :: E_Distr_p_XYZ, E_Distr_p_RLTheta, E_Distr_p_RcThcPhic ! positron
   real(8), dimension(:,:,:,:), allocatable :: E_Distr_h_XYZ, E_Distr_h_RLTheta, E_Distr_h_RcThcPhic ! hole
   real(8), dimension(:,:,:), allocatable :: E_Distr_SHI_XYZ, E_Distr_SHI_RLTheta, E_Distr_SHI_RcThcPhic ! SHI
   real(8), dimension(:,:,:), allocatable :: E_Distr_a_XYZ, E_Distr_a_RLTheta, E_Distr_a_RcThcPhic ! Atom
   real(8), dimension(:,:,:), allocatable :: E_Distr_mu_XYZ, E_Distr_mu_RLTheta, E_Distr_mu_RcThcPhic ! muon
   ! And regular variables:
   integer :: iter, Nsiz_vel, Nsiz(10), Nsiz1(10), Nsiz2(10), Nsiz3(10), Nspec_siz0(10), Nspec_siz1(10), Nspec_siz2(10), Nspec_siz3(10), &
                  Ntheta_siz0(10), Ntheta_siz1(10), Ntheta_siz2(10), Ntheta_siz3(10)
   integer, dimension(:), allocatable :: Nsiz_VB
   real(8) :: rNMC
   
!    print*, 'analyze_MC_output_data 0'

   ! Nullify the data to start over collecting the distributions:
   call reset_output_arrays(out_data)   ! below

! print*, 'analyze_MC_output_data 1'

   !ttttttttttttttttttttttttttttttttttttttttttttttt
   ! Total numbers:
   call get_total_values(used_target, numpar, MC, out_data) ! below

! print*, 'analyze_MC_output_data 2'
   
   !ddddddddddddddddddddddddddd
   ! Distributions:
   
   ! Set sizes of the arrays:
   call set_sizes_temp_output(numpar, Nsiz, Nsiz1, Nsiz2, Nsiz3, Nsiz_vel, Nspec_siz0, Nspec_siz1, Nspec_siz2, Nspec_siz3, &
                              Ntheta_siz0, Ntheta_siz1, Ntheta_siz2, Ntheta_siz3) ! below
   ! Special grid for valence band holes according to DOS:
   call get_size_VB_output(used_target, Nsiz_VB)    ! below

!    print*, 'analyze_MC_output_data 3'

   ! Allocate the arrays when needed:
   if (.not.allocated(Spectrum_ph)) then    ! all spectra
      call allocate_spectra_arrays(Nsiz, Nsiz_VB, Nspec_siz0, Nspec_siz1, Nspec_siz2, Nspec_siz3, &
            Spectrum_ph, Spectrum_e, Spectrum_p, Spectrum_h, Spectrum_SHI, Spectrum_mu, &
            Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X, Spectra_mu_X, &
            Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y, Spectra_mu_Y, &
            Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z, Spectra_mu_Z, &
            Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R, Spectra_mu_R  ) ! below
   endif
   
!    print*, 'analyze_MC_output_data 4'

   if (.not.allocated(Vel_theta_ph)) then   ! all velosity distributions
      call allocate_vel_theta_arrays(Nsiz_vel, Ntheta_siz0, Ntheta_siz1, Ntheta_siz2, Ntheta_siz3, &
            Vel_theta_ph, Vel_theta_e, Vel_theta_p, Vel_theta_h, Vel_theta_SHI, Vel_theta_mu, &
            Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X, Theta_mu_X, &
            Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y, Theta_mu_Y, &
            Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z, Theta_mu_Z, &
            Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R, Theta_mu_R  ) ! below
   endif
   
!    print*, 'analyze_MC_output_data 5'

   ! Spatial distributions:
   if (.not.allocated(Distr_ph_X)) then
      call allocate_output_arrays(Nsiz, Nsiz1, Nsiz2, Nsiz3, numpar%N_sh_tot, &
       Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_ph_R, Distr_ph_L, Distr_ph_Theta, Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic, &
       Distr_e_X, Distr_e_Y, Distr_e_Z, Distr_e_R, Distr_e_L, Distr_e_Theta, Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic, &
       Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_p_R, Distr_p_L, Distr_p_Theta, Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic, &
       Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_h_R, Distr_h_L, Distr_h_Theta, Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic, &
       Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, Distr_SHI_R, Distr_SHI_L, Distr_SHI_Theta, Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic, &
       Distr_a_X, Distr_a_Y, Distr_a_Z, Distr_a_R, Distr_a_L, Distr_a_Theta, Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic, &
       Distr_mu_X, Distr_mu_Y, Distr_mu_Z, Distr_mu_R, Distr_mu_L, Distr_mu_Theta, Distr_mu_Rc, Distr_mu_Thetac, Distr_mu_Phic, &
       Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ, Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta, &
       Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic, &
       Distr_e_XY, Distr_e_YZ, Distr_e_XZ, Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta, &
       Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic, &
       Distr_p_XY, Distr_p_YZ, Distr_p_XZ, Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta, &
       Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic, &
       Distr_h_XY, Distr_h_YZ, Distr_h_XZ, Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta, &
       Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic, &
       Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ, Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta, &
       Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic, &
       Distr_a_XY, Distr_a_YZ, Distr_a_XZ, Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta, &
       Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic, &
       Distr_mu_XY, Distr_mu_YZ, Distr_mu_XZ, Distr_mu_RL, Distr_mu_RTheta, Distr_mu_LTheta, &
       Distr_mu_RcThc, Distr_mu_RcPhic, Distr_mu_ThcPhic, &
       Distr_ph_XYZ, Distr_ph_RLTheta, Distr_ph_RcThcPhic, &
       Distr_e_XYZ, Distr_e_RLTheta, Distr_e_RcThcPhic, &
       Distr_p_XYZ, Distr_p_RLTheta, Distr_p_RcThcPhic, &
       Distr_h_XYZ, Distr_h_RLTheta, Distr_h_RcThcPhic, &
       Distr_SHI_XYZ, Distr_SHI_RLTheta, Distr_SHI_RcThcPhic, &
       Distr_a_XYZ, Distr_a_RLTheta, Distr_a_RcThcPhic, &
       Distr_mu_XYZ, Distr_mu_RLTheta, Distr_mu_RcThcPhic)    ! below
   endif
   
!    print*, 'analyze_MC_output_data 6'

   ! And for the energy arrays, use the same sizes:
   if (.not.allocated(E_Distr_ph_X)) then
      call allocate_output_arrays(Nsiz, Nsiz1, Nsiz2, Nsiz3, numpar%N_sh_tot, &
       E_Distr_ph_X, E_Distr_ph_Y, E_Distr_ph_Z, E_Distr_ph_R, E_Distr_ph_L, E_Distr_ph_Theta, E_Distr_ph_Rc, E_Distr_ph_Thetac, E_Distr_ph_Phic, &
       E_Distr_e_X, E_Distr_e_Y, E_Distr_e_Z, E_Distr_e_R, E_Distr_e_L, E_Distr_e_Theta, E_Distr_e_Rc, E_Distr_e_Thetac, E_Distr_e_Phic, &
       E_Distr_p_X, E_Distr_p_Y, E_Distr_p_Z, E_Distr_p_R, E_Distr_p_L, E_Distr_p_Theta, E_Distr_p_Rc, E_Distr_p_Thetac, E_Distr_p_Phic, &
       E_Distr_h_X, E_Distr_h_Y, E_Distr_h_Z, E_Distr_h_R, E_Distr_h_L, E_Distr_h_Theta, E_Distr_h_Rc, E_Distr_h_Thetac, E_Distr_h_Phic, &
       E_Distr_SHI_X, E_Distr_SHI_Y, E_Distr_SHI_Z, E_Distr_SHI_R, E_Distr_SHI_L, E_Distr_SHI_Theta, E_Distr_SHI_Rc, E_Distr_SHI_Thetac, E_Distr_SHI_Phic, &
       E_Distr_a_X, E_Distr_a_Y, E_Distr_a_Z, E_Distr_a_R, E_Distr_a_L, E_Distr_a_Theta, E_Distr_a_Rc, E_Distr_a_Thetac, E_Distr_a_Phic, &
       E_Distr_mu_X, E_Distr_mu_Y, E_Distr_mu_Z, E_Distr_mu_R, E_Distr_mu_L, E_Distr_mu_Theta, E_Distr_mu_Rc, E_Distr_mu_Thetac, E_Distr_mu_Phic, &
       E_Distr_ph_XY, E_Distr_ph_YZ, E_Distr_ph_XZ, E_Distr_ph_RL, E_Distr_ph_RTheta, E_Distr_ph_LTheta, &
       E_Distr_ph_RcThc, E_Distr_ph_RcPhic, E_Distr_ph_ThcPhic, &
       E_Distr_e_XY, E_Distr_e_YZ, E_Distr_e_XZ, E_Distr_e_RL, E_Distr_e_RTheta, E_Distr_e_LTheta, &
       E_Distr_e_RcThc, E_Distr_e_RcPhic, E_Distr_e_ThcPhic, &
       E_Distr_p_XY, E_Distr_p_YZ, E_Distr_p_XZ, E_Distr_p_RL, E_Distr_p_RTheta, E_Distr_p_LTheta, &
       E_Distr_p_RcThc, E_Distr_p_RcPhic, E_Distr_p_ThcPhic, &
       E_Distr_h_XY, E_Distr_h_YZ, E_Distr_h_XZ, E_Distr_h_RL, E_Distr_h_RTheta, E_Distr_h_LTheta, &
       E_Distr_h_RcThc, E_Distr_h_RcPhic, E_Distr_h_ThcPhic, &
       E_Distr_SHI_XY, E_Distr_SHI_YZ, E_Distr_SHI_XZ, E_Distr_SHI_RL, E_Distr_SHI_RTheta, E_Distr_SHI_LTheta, &
       E_Distr_SHI_RcThc, E_Distr_SHI_RcPhic, E_Distr_SHI_ThcPhic, &
       E_Distr_a_XY, E_Distr_a_YZ, E_Distr_a_XZ, E_Distr_a_RL, E_Distr_a_RTheta, E_Distr_a_LTheta, &
       E_Distr_a_RcThc, E_Distr_a_RcPhic, E_Distr_a_ThcPhic, &
       E_Distr_mu_XY, E_Distr_mu_YZ, E_Distr_mu_XZ, E_Distr_mu_RL, E_Distr_mu_RTheta, E_Distr_mu_LTheta, &
       E_Distr_mu_RcThc, E_Distr_mu_RcPhic, E_Distr_mu_ThcPhic, &
       E_Distr_ph_XYZ, E_Distr_ph_RLTheta, E_Distr_ph_RcThcPhic, &
       E_Distr_e_XYZ, E_Distr_e_RLTheta, E_Distr_e_RcThcPhic, &
       E_Distr_p_XYZ, E_Distr_p_RLTheta, E_Distr_p_RcThcPhic, &
       E_Distr_h_XYZ, E_Distr_h_RLTheta, E_Distr_h_RcThcPhic, &
       E_Distr_SHI_XYZ, E_Distr_SHI_RLTheta, E_Distr_SHI_RcThcPhic, &
       E_Distr_a_XYZ, E_Distr_a_RLTheta, E_Distr_a_RcThcPhic, &
       E_Distr_mu_XYZ, E_Distr_mu_RLTheta, E_Distr_mu_RcThcPhic  )    ! below
   endif
   
!    print*, 'analyze_MC_output_data 7'

   ! Now, having all arrays needed, go on and fill them with data:
   call sort_data(used_target, MC, numpar, tim, Spectrum_ph, Spectrum_e, Spectrum_p, Spectrum_h, Spectrum_SHI, Spectrum_mu, &
    Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X, Spectra_mu_X, &
    Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y, Spectra_mu_Y, &
    Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z, Spectra_mu_Z, &
    Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R, Spectra_mu_R, &
    Vel_theta_ph, Vel_theta_e, Vel_theta_p, Vel_theta_h, Vel_theta_SHI, Vel_theta_mu, &
    Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X, Theta_mu_X, &
    Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y, Theta_mu_Y, &
    Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z, Theta_mu_Z, &
    Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R, Theta_mu_R, &
    Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_ph_R, Distr_ph_L, Distr_ph_Theta, Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic, &
    Distr_e_X, Distr_e_Y, Distr_e_Z, Distr_e_R, Distr_e_L, Distr_e_Theta, Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic, &
    Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_p_R, Distr_p_L, Distr_p_Theta, Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic, &
    Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_h_R, Distr_h_L, Distr_h_Theta, Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic, &
    Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, Distr_SHI_R, Distr_SHI_L, Distr_SHI_Theta, Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic, &
    Distr_a_X, Distr_a_Y, Distr_a_Z, Distr_a_R, Distr_a_L, Distr_a_Theta, Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic, &
    Distr_mu_X, Distr_mu_Y, Distr_mu_Z, Distr_mu_R, Distr_mu_L, Distr_mu_Theta, Distr_mu_Rc, Distr_mu_Thetac, Distr_mu_Phic, &
    Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ, Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta, &
    Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic, &
    Distr_e_XY, Distr_e_YZ, Distr_e_XZ, Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta, &
    Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic, &
    Distr_p_XY, Distr_p_YZ, Distr_p_XZ, Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta, &
    Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic, &
    Distr_h_XY, Distr_h_YZ, Distr_h_XZ, Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta, &
    Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic, &
    Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ, Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta, &
    Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic, &
    Distr_a_XY, Distr_a_YZ, Distr_a_XZ, Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta, &
    Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic, &
    Distr_mu_XY, Distr_mu_YZ, Distr_mu_XZ, Distr_mu_RL, Distr_mu_RTheta, Distr_mu_LTheta, &
    Distr_mu_RcThc, Distr_mu_RcPhic, Distr_mu_ThcPhic, &
    Distr_ph_XYZ, Distr_ph_RLTheta, Distr_ph_RcThcPhic, &
    Distr_e_XYZ, Distr_e_RLTheta, Distr_e_RcThcPhic, &
    Distr_p_XYZ, Distr_p_RLTheta, Distr_p_RcThcPhic, &
    Distr_h_XYZ, Distr_h_RLTheta, Distr_h_RcThcPhic, &
    Distr_SHI_XYZ, Distr_SHI_RLTheta, Distr_SHI_RcThcPhic,&
    Distr_a_XYZ, Distr_a_RLTheta, Distr_a_RcThcPhic,&
    Distr_mu_XYZ, Distr_mu_RLTheta, Distr_mu_RcThcPhic, &
    E_Distr_ph_X, E_Distr_ph_Y, E_Distr_ph_Z, E_Distr_ph_R, E_Distr_ph_L, &
    E_Distr_ph_Theta, E_Distr_ph_Rc, E_Distr_ph_Thetac, E_Distr_ph_Phic, &
    E_Distr_e_X, E_Distr_e_Y, E_Distr_e_Z, E_Distr_e_R, E_Distr_e_L, E_Distr_e_Theta, E_Distr_e_Rc, E_Distr_e_Thetac, E_Distr_e_Phic, &
    E_Distr_p_X, E_Distr_p_Y, E_Distr_p_Z, E_Distr_p_R, E_Distr_p_L, E_Distr_p_Theta, E_Distr_p_Rc, E_Distr_p_Thetac, E_Distr_p_Phic, &
    E_Distr_h_X, E_Distr_h_Y, E_Distr_h_Z, E_Distr_h_R, E_Distr_h_L, E_Distr_h_Theta, E_Distr_h_Rc, E_Distr_h_Thetac, E_Distr_h_Phic, &
    E_Distr_SHI_X, E_Distr_SHI_Y, E_Distr_SHI_Z, E_Distr_SHI_R, E_Distr_SHI_L, &
    E_Distr_SHI_Theta, E_Distr_SHI_Rc, E_Distr_SHI_Thetac, E_Distr_SHI_Phic, &
    E_Distr_a_X, E_Distr_a_Y, E_Distr_a_Z, E_Distr_a_R, E_Distr_a_L, E_Distr_a_Theta, E_Distr_a_Rc, E_Distr_a_Thetac, E_Distr_a_Phic, &
    E_Distr_mu_X, E_Distr_mu_Y, E_Distr_mu_Z, E_Distr_mu_R, E_Distr_mu_L, E_Distr_mu_Theta, E_Distr_mu_Rc, E_Distr_mu_Thetac, E_Distr_mu_Phic, &
    E_Distr_ph_XY, E_Distr_ph_YZ, E_Distr_ph_XZ, E_Distr_ph_RL, E_Distr_ph_RTheta, E_Distr_ph_LTheta, &
    E_Distr_ph_RcThc, E_Distr_ph_RcPhic, E_Distr_ph_ThcPhic, &
    E_Distr_e_XY, E_Distr_e_YZ, E_Distr_e_XZ, E_Distr_e_RL, E_Distr_e_RTheta, E_Distr_e_LTheta, &
    E_Distr_e_RcThc, E_Distr_e_RcPhic, E_Distr_e_ThcPhic, &
    E_Distr_p_XY, E_Distr_p_YZ, E_Distr_p_XZ, E_Distr_p_RL, E_Distr_p_RTheta, E_Distr_p_LTheta, &
    E_Distr_p_RcThc, E_Distr_p_RcPhic, E_Distr_p_ThcPhic, &
    E_Distr_h_XY, E_Distr_h_YZ, E_Distr_h_XZ, E_Distr_h_RL, E_Distr_h_RTheta, E_Distr_h_LTheta, &
    E_Distr_h_RcThc, E_Distr_h_RcPhic, E_Distr_h_ThcPhic, &
    E_Distr_SHI_XY, E_Distr_SHI_YZ, E_Distr_SHI_XZ, E_Distr_SHI_RL, E_Distr_SHI_RTheta, E_Distr_SHI_LTheta, &
    E_Distr_SHI_RcThc, E_Distr_SHI_RcPhic, E_Distr_SHI_ThcPhic, &
    E_Distr_a_XY, E_Distr_a_YZ, E_Distr_a_XZ, E_Distr_a_RL, E_Distr_a_RTheta, E_Distr_a_LTheta, &
    E_Distr_a_RcThc, E_Distr_a_RcPhic, E_Distr_a_ThcPhic, &
    E_Distr_mu_XY, E_Distr_mu_YZ, E_Distr_mu_XZ, E_Distr_mu_RL, E_Distr_mu_RTheta, E_Distr_mu_LTheta, &
    E_Distr_mu_RcThc, E_Distr_mu_RcPhic, E_Distr_mu_ThcPhic, &
    E_Distr_ph_XYZ, E_Distr_ph_RLTheta, E_Distr_ph_RcThcPhic, &
    E_Distr_e_XYZ, E_Distr_e_RLTheta, E_Distr_e_RcThcPhic, &
    E_Distr_p_XYZ, E_Distr_p_RLTheta, E_Distr_p_RcThcPhic, &
    E_Distr_h_XYZ, E_Distr_h_RLTheta, E_Distr_h_RcThcPhic, &
    E_Distr_SHI_XYZ, E_Distr_SHI_RLTheta, E_Distr_SHI_RcThcPhic, &
    E_Distr_a_XYZ, E_Distr_a_RLTheta, E_Distr_a_RcThcPhic, &
    E_Distr_mu_XYZ, E_Distr_mu_RLTheta, E_Distr_mu_RcThcPhic  )    ! below
   
!    print*, 'analyze_MC_output_data 8'

   ! Get the averaged data, noramlized per number of MC iterations:
    if (numpar%NMC > 0) then
       rNMC = 1.0d0/dble(numpar%NMC) ! convert into double-real type, and take an inverse
    else
       rNMC = 0.0d0
    endif
    ! Spectra:
    out_data%Spectrum_ph = Spectrum_ph * rNMC
    out_data%Spectrum_e = Spectrum_e * rNMC
    out_data%Spectrum_p = Spectrum_p * rNMC
    out_data%Spectrum_h = Spectrum_h * rNMC
    out_data%Spectrum_SHI = Spectrum_SHI * rNMC
    ! Velosity theta distribution:
    out_data%Vel_theta_ph = Vel_theta_ph * rNMC
    out_data%Vel_theta_e = Vel_theta_e * rNMC
    out_data%Vel_theta_p = Vel_theta_p * rNMC
    out_data%Vel_theta_h = Vel_theta_h * rNMC
    out_data%Vel_theta_SHI = Vel_theta_SHI * rNMC
    ! Spectra vs space in 1d:
    out_data%Spectra_ph_X = Spectra_ph_X * rNMC
    out_data%Spectra_e_X = Spectra_e_X * rNMC
    out_data%Spectra_p_X = Spectra_p_X * rNMC
    out_data%Spectra_h_X = Spectra_h_X * rNMC
    out_data%Spectra_SHI_X = Spectra_SHI_X * rNMC
    out_data%Spectra_ph_Y = Spectra_ph_Y * rNMC
    out_data%Spectra_e_Y = Spectra_e_Y * rNMC
    out_data%Spectra_p_Y = Spectra_p_Y * rNMC
    out_data%Spectra_h_Y = Spectra_h_Y * rNMC
    out_data%Spectra_SHI_Y = Spectra_SHI_Y * rNMC
    out_data%Spectra_ph_Z = Spectra_ph_Z * rNMC
    out_data%Spectra_e_Z = Spectra_e_Z * rNMC
    out_data%Spectra_p_Z = Spectra_p_Z * rNMC
    out_data%Spectra_h_Z = Spectra_h_Z * rNMC
    out_data%Spectra_SHI_Z = Spectra_SHI_Z * rNMC
    out_data%Spectra_ph_R = Spectra_ph_R * rNMC
    out_data%Spectra_e_R = Spectra_e_R * rNMC
    out_data%Spectra_p_R = Spectra_p_R * rNMC
    out_data%Spectra_h_R = Spectra_h_R * rNMC
    out_data%Spectra_SHI_R = Spectra_SHI_R * rNMC
    ! Theta distribution vs space in 1d:
    out_data%Theta_ph_X = Theta_ph_X * rNMC
    out_data%Theta_e_X = Theta_e_X * rNMC
    out_data%Theta_p_X = Theta_p_X * rNMC
    out_data%Theta_h_X = Theta_h_X * rNMC
    out_data%Theta_SHI_X = Theta_SHI_X * rNMC
    out_data%Theta_ph_Y = Theta_ph_Y * rNMC
    out_data%Theta_e_Y = Theta_e_Y * rNMC
    out_data%Theta_p_Y = Theta_p_Y * rNMC
    out_data%Theta_h_Y = Theta_h_Y * rNMC
    out_data%Theta_SHI_Y = Theta_SHI_Y * rNMC
    out_data%Theta_ph_Z = Theta_ph_Z * rNMC
    out_data%Theta_e_Z = Theta_e_Z * rNMC
    out_data%Theta_p_Z = Theta_p_Z * rNMC
    out_data%Theta_h_Z = Theta_h_Z * rNMC
    out_data%Theta_SHI_Z = Theta_SHI_Z * rNMC
    out_data%Theta_ph_R = Theta_ph_R * rNMC
    out_data%Theta_e_R = Theta_e_R * rNMC
    out_data%Theta_p_R = Theta_p_R * rNMC
    out_data%Theta_h_R = Theta_h_R * rNMC
    out_data%Theta_SHI_R = Theta_SHI_R * rNMC
    ! Spatial distributions:
    out_data%Distr_ph_X = Distr_ph_X * rNMC
    out_data%Distr_ph_Y = Distr_ph_Y * rNMC
    out_data%Distr_ph_Z = Distr_ph_Z * rNMC
    out_data%Distr_ph_R = Distr_ph_R * rNMC
    out_data%Distr_ph_L = Distr_ph_L * rNMC
    out_data%Distr_ph_Theta = Distr_ph_Theta * rNMC
    out_data%Distr_ph_Rc = Distr_ph_Rc * rNMC
    out_data%Distr_ph_Thetac = Distr_ph_Thetac * rNMC
    out_data%Distr_ph_Phic = Distr_ph_Phic * rNMC
    out_data%Distr_e_X = Distr_e_X * rNMC
    out_data%Distr_e_Y = Distr_e_Y * rNMC
    out_data%Distr_e_Z = Distr_e_Z * rNMC
    out_data%Distr_e_R = Distr_e_R * rNMC
    out_data%Distr_e_L = Distr_e_L * rNMC
    out_data%Distr_e_Theta = Distr_e_Theta * rNMC
    out_data%Distr_e_Rc = Distr_e_Rc * rNMC
    out_data%Distr_e_Thetac = Distr_e_Thetac * rNMC
    out_data%Distr_e_Phic = Distr_e_Phic * rNMC
    out_data%Distr_p_X = Distr_p_X * rNMC
    out_data%Distr_p_Y = Distr_p_Y * rNMC
    out_data%Distr_p_Z = Distr_p_Z * rNMC
    out_data%Distr_p_R = Distr_p_R * rNMC
    out_data%Distr_p_L = Distr_p_L * rNMC
    out_data%Distr_p_Theta = Distr_p_Theta * rNMC
    out_data%Distr_p_Rc = Distr_p_Rc * rNMC
    out_data%Distr_p_Thetac = Distr_p_Thetac * rNMC
    out_data%Distr_p_Phic = Distr_p_Phic * rNMC
    out_data%Distr_h_X = Distr_h_X * rNMC
    out_data%Distr_h_Y = Distr_h_Y * rNMC
    out_data%Distr_h_Z = Distr_h_Z * rNMC
    out_data%Distr_h_R = Distr_h_R * rNMC
    out_data%Distr_h_L = Distr_h_L * rNMC
    out_data%Distr_h_Theta = Distr_h_Theta * rNMC
    out_data%Distr_h_Rc = Distr_h_Rc * rNMC
    out_data%Distr_h_Thetac = Distr_h_Thetac * rNMC
    out_data%Distr_h_Phic = Distr_h_Phic * rNMC
    out_data%Distr_SHI_X = Distr_SHI_X * rNMC
    out_data%Distr_SHI_Y = Distr_SHI_Y * rNMC
    out_data%Distr_SHI_Z = Distr_SHI_Z * rNMC
    out_data%Distr_SHI_R = Distr_SHI_R * rNMC
    out_data%Distr_SHI_L = Distr_SHI_L * rNMC
    out_data%Distr_SHI_Theta = Distr_SHI_Theta * rNMC
    out_data%Distr_SHI_Rc = Distr_SHI_Rc * rNMC
    out_data%Distr_SHI_Thetac = Distr_SHI_Thetac * rNMC
    out_data%Distr_SHI_Phic = Distr_SHI_Phic * rNMC
    ! For atoms, we have to add, because the data are nullified at each time step:
    out_data%Distr_a_X = out_data%Distr_a_X + Distr_a_X * rNMC
    out_data%Distr_a_Y = out_data%Distr_a_Y + Distr_a_Y * rNMC
    out_data%Distr_a_Z = out_data%Distr_a_Z + Distr_a_Z * rNMC
    out_data%Distr_a_R = out_data%Distr_a_R + Distr_a_R * rNMC
    out_data%Distr_a_L = out_data%Distr_a_L + Distr_a_L * rNMC
    out_data%Distr_a_Theta = out_data%Distr_a_Theta + Distr_a_Theta * rNMC
    out_data%Distr_a_Rc = out_data%Distr_a_Rc + Distr_a_Rc * rNMC
    out_data%Distr_a_Thetac = out_data%Distr_a_Thetac + Distr_a_Thetac * rNMC
    out_data%Distr_a_Phic = out_data%Distr_a_Phic+ Distr_a_Phic * rNMC
    out_data%Distr_ph_XY = Distr_ph_XY * rNMC
    out_data%Distr_ph_YZ = Distr_ph_YZ * rNMC
    out_data%Distr_ph_XZ = Distr_ph_XZ * rNMC
    out_data%Distr_ph_RL = Distr_ph_RL * rNMC
    out_data%Distr_ph_RTheta = Distr_ph_RTheta * rNMC
    out_data%Distr_ph_LTheta = Distr_ph_LTheta * rNMC
    out_data%Distr_ph_RcThc = Distr_ph_RcThc * rNMC
    out_data%Distr_ph_RcPhic = Distr_ph_RcPhic * rNMC
    out_data%Distr_ph_ThcPhic = Distr_ph_ThcPhic * rNMC
    out_data%Distr_e_XY = Distr_e_XY * rNMC
    out_data%Distr_e_YZ = Distr_e_YZ * rNMC
    out_data%Distr_e_XZ = Distr_e_XZ * rNMC
    out_data%Distr_e_RL = Distr_e_RL * rNMC
    out_data%Distr_e_RTheta = Distr_e_RTheta * rNMC
    out_data%Distr_e_LTheta = Distr_e_LTheta * rNMC
    out_data%Distr_e_RcThc = Distr_e_RcThc * rNMC
    out_data%Distr_e_RcPhic = Distr_e_RcPhic * rNMC
    out_data%Distr_e_ThcPhic = Distr_e_ThcPhic * rNMC
    out_data%Distr_p_XY = Distr_p_XY * rNMC
    out_data%Distr_p_YZ = Distr_p_YZ * rNMC
    out_data%Distr_p_XZ = Distr_p_XZ * rNMC
    out_data%Distr_p_RL = Distr_p_RL * rNMC
    out_data%Distr_p_RTheta = Distr_p_RTheta * rNMC
    out_data%Distr_p_LTheta = Distr_p_LTheta * rNMC
    out_data%Distr_p_RcThc = Distr_p_RcThc * rNMC
    out_data%Distr_p_RcPhic = Distr_p_RcPhic * rNMC
    out_data%Distr_p_ThcPhic = Distr_p_ThcPhic * rNMC
    out_data%Distr_h_XY = Distr_h_XY * rNMC
    out_data%Distr_h_YZ = Distr_h_YZ * rNMC
    out_data%Distr_h_XZ = Distr_h_XZ * rNMC
    out_data%Distr_h_RL = Distr_h_RL * rNMC
    out_data%Distr_h_RTheta = Distr_h_RTheta * rNMC
    out_data%Distr_h_LTheta = Distr_h_LTheta * rNMC
    out_data%Distr_h_RcThc = Distr_h_RcThc * rNMC
    out_data%Distr_h_RcPhic = Distr_h_RcPhic * rNMC
    out_data%Distr_h_ThcPhic = Distr_h_ThcPhic * rNMC
    out_data%Distr_SHI_XY = Distr_SHI_XY * rNMC
    out_data%Distr_SHI_YZ = Distr_SHI_YZ * rNMC
    out_data%Distr_SHI_XZ = Distr_SHI_XZ * rNMC
    out_data%Distr_SHI_RL = Distr_SHI_RL * rNMC
    out_data%Distr_SHI_RTheta = Distr_SHI_RTheta * rNMC
    out_data%Distr_SHI_LTheta = Distr_SHI_LTheta * rNMC
    out_data%Distr_SHI_RcThc = Distr_SHI_RcThc * rNMC
    out_data%Distr_SHI_RcPhic = Distr_SHI_RcPhic * rNMC
    out_data%Distr_SHI_ThcPhic = Distr_SHI_ThcPhic * rNMC
    ! For atoms, we have to add, because the data a nullified at each time step:
    out_data%Distr_a_XY = out_data%Distr_a_XY+ Distr_a_XY * rNMC
    out_data%Distr_a_YZ = out_data%Distr_a_YZ + Distr_a_YZ * rNMC
    out_data%Distr_a_XZ = out_data%Distr_a_XZ + Distr_a_XZ * rNMC
    out_data%Distr_a_RL = out_data%Distr_a_RL + Distr_a_RL * rNMC
    out_data%Distr_a_RTheta = out_data%Distr_a_RTheta + Distr_a_RTheta * rNMC
    out_data%Distr_a_LTheta = out_data%Distr_a_LTheta + Distr_a_LTheta * rNMC
    out_data%Distr_a_RcThc = out_data%Distr_a_RcThc + Distr_a_RcThc * rNMC
    out_data%Distr_a_RcPhic = out_data%Distr_a_RcPhic + Distr_a_RcPhic * rNMC
    out_data%Distr_a_ThcPhic = out_data%Distr_a_ThcPhic + Distr_a_ThcPhic * rNMC
    out_data%Distr_ph_XYZ = Distr_ph_XYZ * rNMC
    out_data%Distr_ph_RLTheta = Distr_ph_RLTheta * rNMC
    out_data%Distr_ph_RcThcPhic = Distr_ph_RcThcPhic * rNMC
    out_data%Distr_e_XYZ = Distr_e_XYZ * rNMC
    out_data%Distr_e_RLTheta = Distr_e_RLTheta * rNMC
    out_data%Distr_e_RcThcPhic = Distr_e_RcThcPhic * rNMC
    out_data%Distr_p_XYZ = Distr_p_XYZ * rNMC
    out_data%Distr_p_RLTheta = Distr_p_RLTheta * rNMC
    out_data%Distr_p_RcThcPhic = Distr_p_RcThcPhic * rNMC
    out_data%Distr_h_XYZ = Distr_h_XYZ * rNMC
    out_data%Distr_h_RLTheta = Distr_h_RLTheta * rNMC
    out_data%Distr_h_RcThcPhic = Distr_h_RcThcPhic * rNMC
    out_data%Distr_SHI_XYZ = Distr_SHI_XYZ * rNMC
    out_data%Distr_SHI_RLTheta = Distr_SHI_RLTheta * rNMC
    out_data%Distr_SHI_RcThcPhic = Distr_SHI_RcThcPhic * rNMC
    ! For atoms, we have to add, because the data a nullified at each time step:
    out_data%Distr_a_XYZ = out_data%Distr_a_XYZ + Distr_a_XYZ * rNMC
    out_data%Distr_a_RLTheta = out_data%Distr_a_RLTheta + Distr_a_RLTheta * rNMC
    out_data%Distr_a_RcThcPhic = out_data%Distr_a_RcThcPhic + Distr_a_RcThcPhic * rNMC
    ! Spatial energy distributions:
    out_data%E_Distr_ph_X = E_Distr_ph_X * rNMC
    out_data%E_Distr_ph_Y = E_Distr_ph_Y * rNMC
    out_data%E_Distr_ph_Z = E_Distr_ph_Z * rNMC
    out_data%E_Distr_ph_R = E_Distr_ph_R * rNMC
    out_data%E_Distr_ph_L = E_Distr_ph_L * rNMC
    out_data%E_Distr_ph_Theta = E_Distr_ph_Theta * rNMC
    out_data%E_Distr_ph_Rc = E_Distr_ph_Rc * rNMC
    out_data%E_Distr_ph_Thetac = E_Distr_ph_Thetac * rNMC
    out_data%E_Distr_ph_Phic = E_Distr_ph_Phic * rNMC
    out_data%E_Distr_e_X = E_Distr_e_X * rNMC
    out_data%E_Distr_e_Y = E_Distr_e_Y * rNMC
    out_data%E_Distr_e_Z = E_Distr_e_Z * rNMC
    out_data%E_Distr_e_R = E_Distr_e_R * rNMC
    out_data%E_Distr_e_L = E_Distr_e_L * rNMC
    out_data%E_Distr_e_Theta = E_Distr_e_Theta * rNMC
    out_data%E_Distr_e_Rc = E_Distr_e_Rc * rNMC
    out_data%E_Distr_e_Thetac = E_Distr_e_Thetac * rNMC
    out_data%E_Distr_e_Phic = E_Distr_e_Phic * rNMC
    out_data%E_Distr_p_X = E_Distr_p_X * rNMC
    out_data%E_Distr_p_Y = E_Distr_p_Y * rNMC
    out_data%E_Distr_p_Z = E_Distr_p_Z * rNMC
    out_data%E_Distr_p_R = E_Distr_p_R * rNMC
    out_data%E_Distr_p_L = E_Distr_p_L * rNMC
    out_data%E_Distr_p_Theta = E_Distr_p_Theta * rNMC
    out_data%E_Distr_p_Rc = E_Distr_p_Rc * rNMC
    out_data%E_Distr_p_Thetac = E_Distr_p_Thetac * rNMC
    out_data%E_Distr_p_Phic = E_Distr_p_Phic * rNMC
    out_data%E_Distr_h_X = E_Distr_h_X * rNMC
    out_data%E_Distr_h_Y = E_Distr_h_Y * rNMC
    out_data%E_Distr_h_Z = E_Distr_h_Z * rNMC
    out_data%E_Distr_h_R = E_Distr_h_R * rNMC
    out_data%E_Distr_h_L = E_Distr_h_L * rNMC
    out_data%E_Distr_h_Theta = E_Distr_h_Theta * rNMC
    out_data%E_Distr_h_Rc = E_Distr_h_Rc * rNMC
    out_data%E_Distr_h_Thetac = E_Distr_h_Thetac * rNMC
    out_data%E_Distr_h_Phic = E_Distr_h_Phic * rNMC
    out_data%E_Distr_SHI_X = E_Distr_SHI_X * rNMC
    out_data%E_Distr_SHI_Y = E_Distr_SHI_Y * rNMC
    out_data%E_Distr_SHI_Z = E_Distr_SHI_Z * rNMC
    out_data%E_Distr_SHI_R = E_Distr_SHI_R * rNMC
    out_data%E_Distr_SHI_L = E_Distr_SHI_L * rNMC
    out_data%E_Distr_SHI_Theta = E_Distr_SHI_Theta * rNMC
    out_data%E_Distr_SHI_Rc = E_Distr_SHI_Rc * rNMC
    out_data%E_Distr_SHI_Thetac = E_Distr_SHI_Thetac * rNMC
    out_data%E_Distr_SHI_Phic = E_Distr_SHI_Phic * rNMC
    ! For atoms, we have to add, because the data are nullified at each time step:
    out_data%E_Distr_a_X = out_data%E_Distr_a_X + E_Distr_a_X * rNMC
    out_data%E_Distr_a_Y = out_data%E_Distr_a_Y + E_Distr_a_Y * rNMC
    out_data%E_Distr_a_Z = out_data%E_Distr_a_Z + E_Distr_a_Z * rNMC
    out_data%E_Distr_a_R = out_data%E_Distr_a_R + E_Distr_a_R * rNMC
    out_data%E_Distr_a_L = out_data%E_Distr_a_L + E_Distr_a_L * rNMC
    out_data%E_Distr_a_Theta = out_data%E_Distr_a_Theta + E_Distr_a_Theta * rNMC
    out_data%E_Distr_a_Rc = out_data%E_Distr_a_Rc + E_Distr_a_Rc * rNMC
    out_data%E_Distr_a_Thetac = out_data%E_Distr_a_Thetac + E_Distr_a_Thetac * rNMC
    out_data%E_Distr_a_Phic = out_data%E_Distr_a_Phic + E_Distr_a_Phic * rNMC
    out_data%E_Distr_ph_XY = E_Distr_ph_XY * rNMC
    out_data%E_Distr_ph_YZ = E_Distr_ph_YZ * rNMC
    out_data%E_Distr_ph_XZ = E_Distr_ph_XZ * rNMC
    out_data%E_Distr_ph_RL = E_Distr_ph_RL * rNMC
    out_data%E_Distr_ph_RTheta = E_Distr_ph_RTheta * rNMC
    out_data%E_Distr_ph_LTheta = E_Distr_ph_LTheta * rNMC
    out_data%E_Distr_ph_RcThc = E_Distr_ph_RcThc * rNMC
    out_data%E_Distr_ph_RcPhic = E_Distr_ph_RcPhic * rNMC
    out_data%E_Distr_ph_ThcPhic = E_Distr_ph_ThcPhic * rNMC
    out_data%E_Distr_e_XY = E_Distr_e_XY * rNMC
    out_data%E_Distr_e_YZ = E_Distr_e_YZ * rNMC
    out_data%E_Distr_e_XZ = E_Distr_e_XZ * rNMC
    out_data%E_Distr_e_RL = E_Distr_e_RL * rNMC
    out_data%E_Distr_e_RTheta = E_Distr_e_RTheta * rNMC
    out_data%E_Distr_e_LTheta = E_Distr_e_LTheta * rNMC
    out_data%E_Distr_e_RcThc = E_Distr_e_RcThc * rNMC
    out_data%E_Distr_e_RcPhic = E_Distr_e_RcPhic * rNMC
    out_data%E_Distr_e_ThcPhic = E_Distr_e_ThcPhic * rNMC
    out_data%E_Distr_p_XY = E_Distr_p_XY * rNMC
    out_data%E_Distr_p_YZ = E_Distr_p_YZ * rNMC
    out_data%E_Distr_p_XZ = E_Distr_p_XZ * rNMC
    out_data%E_Distr_p_RL = E_Distr_p_RL * rNMC
    out_data%E_Distr_p_RTheta = E_Distr_p_RTheta * rNMC
    out_data%E_Distr_p_LTheta = E_Distr_p_LTheta * rNMC
    out_data%E_Distr_p_RcThc = E_Distr_p_RcThc * rNMC
    out_data%E_Distr_p_RcPhic = E_Distr_p_RcPhic * rNMC
    out_data%E_Distr_p_ThcPhic = E_Distr_p_ThcPhic * rNMC
    out_data%E_Distr_h_XY = E_Distr_h_XY * rNMC
    out_data%E_Distr_h_YZ = E_Distr_h_YZ * rNMC
    out_data%E_Distr_h_XZ = E_Distr_h_XZ * rNMC
    out_data%E_Distr_h_RL = E_Distr_h_RL * rNMC
    out_data%E_Distr_h_RTheta = E_Distr_h_RTheta * rNMC
    out_data%E_Distr_h_LTheta = E_Distr_h_LTheta * rNMC
    out_data%E_Distr_h_RcThc = E_Distr_h_RcThc * rNMC
    out_data%E_Distr_h_RcPhic = E_Distr_h_RcPhic * rNMC
    out_data%E_Distr_h_ThcPhic = E_Distr_h_ThcPhic * rNMC
    out_data%E_Distr_SHI_XY = E_Distr_SHI_XY * rNMC
    out_data%E_Distr_SHI_YZ = E_Distr_SHI_YZ * rNMC
    out_data%E_Distr_SHI_XZ = E_Distr_SHI_XZ * rNMC
    out_data%E_Distr_SHI_RL = E_Distr_SHI_RL * rNMC
    out_data%E_Distr_SHI_RTheta = E_Distr_SHI_RTheta * rNMC
    out_data%E_Distr_SHI_LTheta = E_Distr_SHI_LTheta * rNMC
    out_data%E_Distr_SHI_RcThc = E_Distr_SHI_RcThc * rNMC
    out_data%E_Distr_SHI_RcPhic = E_Distr_SHI_RcPhic * rNMC
    out_data%E_Distr_SHI_ThcPhic = E_Distr_SHI_ThcPhic * rNMC
    ! For atoms, we have to add, because the data are nullified at each time step:
    out_data%E_Distr_a_XY = out_data%E_Distr_a_XY+ E_Distr_a_XY * rNMC
    out_data%E_Distr_a_YZ = out_data%E_Distr_a_YZ + E_Distr_a_YZ * rNMC
    out_data%E_Distr_a_XZ = out_data%E_Distr_a_XZ + E_Distr_a_XZ * rNMC
    out_data%E_Distr_a_RL = out_data%E_Distr_a_RL + E_Distr_a_RL * rNMC
    out_data%E_Distr_a_RTheta = out_data%E_Distr_a_RTheta + E_Distr_a_RTheta * rNMC
    out_data%E_Distr_a_LTheta = out_data%E_Distr_a_LTheta + E_Distr_a_LTheta * rNMC
    out_data%E_Distr_a_RcThc = out_data%E_Distr_a_RcThc + E_Distr_a_RcThc * rNMC
    out_data%E_Distr_a_RcPhic = out_data%E_Distr_a_RcPhic + E_Distr_a_RcPhic * rNMC
    out_data%E_Distr_a_ThcPhic = out_data%E_Distr_a_ThcPhic + E_Distr_a_ThcPhic * rNMC
    out_data%E_Distr_ph_XYZ = E_Distr_ph_XYZ * rNMC
    out_data%E_Distr_ph_RLTheta = E_Distr_ph_RLTheta * rNMC
    out_data%E_Distr_ph_RcThcPhic = E_Distr_ph_RcThcPhic * rNMC
    out_data%E_Distr_e_XYZ = E_Distr_e_XYZ * rNMC
    out_data%E_Distr_e_RLTheta = E_Distr_e_RLTheta * rNMC
    out_data%E_Distr_e_RcThcPhic = E_Distr_e_RcThcPhic * rNMC
    out_data%E_Distr_p_XYZ = E_Distr_p_XYZ * rNMC
    out_data%E_Distr_p_RLTheta = E_Distr_p_RLTheta * rNMC
    out_data%E_Distr_p_RcThcPhic = E_Distr_p_RcThcPhic * rNMC
    out_data%E_Distr_h_XYZ = E_Distr_h_XYZ * rNMC
    out_data%E_Distr_h_RLTheta = E_Distr_h_RLTheta * rNMC
    out_data%E_Distr_h_RcThcPhic = E_Distr_h_RcThcPhic * rNMC
    out_data%E_Distr_SHI_XYZ = E_Distr_SHI_XYZ * rNMC
    out_data%E_Distr_SHI_RLTheta = E_Distr_SHI_RLTheta * rNMC
    out_data%E_Distr_SHI_RcThcPhic = E_Distr_SHI_RcThcPhic * rNMC
    ! For atoms, we have to add, because the data are nullified at each time step:
    out_data%E_Distr_a_XYZ = out_data%E_Distr_a_XYZ+ E_Distr_a_XYZ * rNMC
    out_data%E_Distr_a_RLTheta = out_data%E_Distr_a_RLTheta + E_Distr_a_RLTheta * rNMC
    out_data%E_Distr_a_RcThcPhic = out_data%E_Distr_a_RcThcPhic + E_Distr_a_RcThcPhic * rNMC

! print*, 'analyze_MC_output_data 9'

   ! Deactivate all atomic scattering events as objects, to start over in the next timestep:
   call renew_atomic_arrays(MC) ! module "MC_general_tools"


! print*, 'analyze_MC_output_data END'
end subroutine analyze_MC_output_data



subroutine get_total_values(used_target, numpar, MC, out_data)
   type(Matter), intent(in) :: used_target  ! parameters of the target
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   type(MC_arrays), dimension(:), intent(inout) :: MC	! all MC arrays for all particles; size = number of iterations
   type(output_data), intent(inout) :: out_data  ! all output data (distributions etc.)
   real(8) :: Eph, Ee, Eh_kin, Eh_pot, Ep, Eat, rNMC
   real(8) :: Eph_high, Ee_high, Eh_kin_high, Eh_pot_high, Ep_high
   real(8) :: Nph_high, Ne_high, Nh_high, Np_high
   integer :: i, Nsiz, j
   
   Nsiz = size(MC)
   ! just in case:
   if (Nsiz < 1) print*, 'Error in get_total_values', Nsiz

   if (numpar%NMC > 0) then
      rNMC = 1.0d0/dble(numpar%NMC) ! convert into double-real type, and take an inverse
   else
      rNMC = 1.0d0
   endif
   
   ! Sum all numbers of particles:
   out_data%Nph = SUM(MC(:)%N_ph) * rNMC
   out_data%Ne = SUM(MC(:)%N_e) * rNMC
   out_data%Nh = SUM(MC(:)%N_h) * rNMC
   out_data%Np = SUM(MC(:)%N_p) * rNMC
   
   ! Just to start
   Eph = 0.0d0
   Ee = 0.0d0
   Eh_kin = 0.0d0
   Eh_pot = 0.0d0
   Ep = 0.0d0
   Eat = 0.0d0
   Eph_high = 0.0d0
   Ee_high = 0.0d0
   Eh_kin_high = 0.0d0
   Eh_pot_high = 0.0d0
   Ep_high = 0.0d0
   Nph_high = 0.0d0
   Ne_high = 0.0d0
   Nh_high = 0.0d0
   Np_high = 0.0d0
   
   ! Sum their kinetic energies:
   !$omp parallel private (i, j) shared(Nsiz, Eph, Ee, Eh_kin, Eh_pot, Ep, Eat)
   !$omp do schedule(dynamic) reduction( + : Eph, Ee, Eh_kin, Eh_pot, Ep, Eat, &
   !$    Eph_high, Ee_high, Eh_kin_high, Eh_pot_high, Ep_high, &
   !$    Nph_high, Ne_high, Nh_high, Np_high)
   do i = 1, Nsiz   ! over MC iterations
      Eph = Eph + SUM(MC(i)%MC_Photons(:)%Ekin, MASK = MC(i)%MC_Photons(:)%active)    ! photons
      Ee = Ee + SUM(MC(i)%MC_Electrons(:)%Ekin, MASK = MC(i)%MC_Electrons(:)%active) ! electrons
      Eh_kin = Eh_kin + SUM(MC(i)%MC_Holes(:)%Ekin, MASK = MC(i)%MC_Holes(:)%active) ! kinetic energy of valence holes
      ! potential energy of valence holes:
      Eh_pot = Eh_pot + SUM( used_target%Material(MC(i)%MC_Holes(1:MC(i)%N_h)%in_target)%DOS%Egap , MASK = MC(i)%MC_Holes(1:MC(i)%N_h)%valent )
      do j = 1, MC(i)%N_h   ! Core holes energies (ionization potentials):
         if (.not.MC(i)%MC_Holes(j)%valent) then
            if (MC(i)%MC_Holes(j)%KOA == 0) then  ! inconsistency found
               print*, 'ERROR in get_total_values:', MC(i)%MC_Holes(j)%valent, MC(i)%MC_Holes(j)%KOA, MC(i)%MC_Holes(j)%Sh
               print*, 'R=', MC(i)%MC_Holes(j)%R
               print*, 'Iteration n Prtcl #', i, j
               print*, 'times:', MC(i)%MC_Holes(j)%ti, MC(i)%MC_Holes(j)%t_sc
            endif
            Eh_pot = Eh_pot + used_target%Material(MC(i)%MC_Holes(j)%in_target)%Elements(MC(i)%MC_Holes(j)%KOA)%Ip(MC(i)%MC_Holes(j)%Sh)
         endif
      enddo
      Ep = Ep + SUM(MC(i)%MC_Positrons(:)%Ekin, MASK = MC(i)%MC_Positrons(:)%active) ! positons
      Eat = Eat + SUM(MC(i)%MC_Atoms_events(:)%Ekin, MASK = MC(i)%MC_Atoms_events(:)%active) ! events of elastic energy transfer to atoms

      ! High-energy particles (above cut-off):
      Nph_high = Nph_high + COUNT(MASK = (MC(i)%MC_Photons(:)%active .and. (MC(i)%MC_Photons(:)%Ekin > numpar%Ph_Cutoff)) )    ! photons
      Ne_high = Ne_high + COUNT(MASK = (MC(i)%MC_Electrons(:)%active .and. (MC(i)%MC_Electrons(:)%Ekin > numpar%El_Cutoff)) )  ! electrons
      Nh_high = Nh_high + COUNT(MASK = (MC(i)%MC_Holes(:)%active .and. (MC(i)%MC_Holes(:)%Ekin > numpar%H_Cutoff)) )    ! holes
      Np_high = Np_high + COUNT(MASK = (MC(i)%MC_Positrons(:)%active .and. (MC(i)%MC_Positrons(:)%Ekin > numpar%Pos_Cutoff)) ) ! positrons
      Eph_high = Eph_high + SUM(MC(i)%MC_Photons(:)%Ekin, MASK = (MC(i)%MC_Photons(:)%active .and. &
                (MC(i)%MC_Photons(:)%Ekin > numpar%Ph_Cutoff)) )    ! photons
      Ee_high = Ee_high + SUM(MC(i)%MC_Electrons(:)%Ekin, MASK = (MC(i)%MC_Electrons(:)%active .and. &
                (MC(i)%MC_Electrons(:)%Ekin > numpar%El_Cutoff)) ) ! electrons
      Eh_kin_high = Eh_kin_high + SUM(MC(i)%MC_Holes(:)%Ekin, MASK = (MC(i)%MC_Holes(:)%active .and. &
                (MC(i)%MC_Holes(:)%Ekin > numpar%H_Cutoff)) ) ! kinetic energy of valence holes
      Eh_pot_high = Eh_pot_high + SUM( used_target%Material(MC(i)%MC_Holes(1:MC(i)%N_h)%in_target)%DOS%Egap , &
                MASK = (MC(i)%MC_Holes(1:MC(i)%N_h)%valent  .and. (MC(i)%MC_Holes(:)%Ekin > numpar%H_Cutoff)) )
      do j = 1, MC(i)%N_h   ! Core holes energies (ionization potentials):
         if (.not.MC(i)%MC_Holes(j)%valent) then
            Eh_pot_high = Eh_pot_high + &
            used_target%Material(MC(i)%MC_Holes(j)%in_target)%Elements(MC(i)%MC_Holes(j)%KOA)%Ip(MC(i)%MC_Holes(j)%Sh)
         endif
      enddo
      Ep_high = Ep_high + SUM(MC(i)%MC_Positrons(:)%Ekin, MASK = (MC(i)%MC_Positrons(:)%active .and. &
                (MC(i)%MC_Positrons(:)%Ekin > numpar%Pos_Cutoff)) ) ! positons

!       print*, 'Nsiz', i, Nsiz, rNMC
!       if (Ee < 1.0d-6) then
!          print*, 'get_total_values', Ee, COUNT(MC(i)%MC_Electrons(:)%active)
!       endif
   enddo
   !$omp enddo 
   !$omp end parallel

   out_data%Eph = Eph * rNMC
   out_data%Ee = Ee * rNMC
   out_data%Eh_kin = Eh_kin * rNMC
   out_data%Eh_pot = Eh_pot * rNMC
   out_data%Ep = Ep * rNMC
   out_data%Eat = out_data%Eat + Eat * rNMC ! here we add the data, since atomic data are nullified after each time-step
   ! Energies of high-energy particles (above cut off):
   out_data%Eph_high = Eph_high * rNMC
   out_data%Ee_high = Ee_high * rNMC
   out_data%Eh_kin_high = Eh_kin_high * rNMC
   out_data%Eh_pot_high = Eh_pot_high * rNMC
   out_data%Ep_high = Ep_high * rNMC
   out_data%Nph_high = Nph_high * rNMC
   out_data%Ne_high = Ne_high * rNMC
   out_data%Nh_high = Nh_high * rNMC
   out_data%Np_high = Np_high * rNMC

!    print*, 'Es:', Eph, Ee, Eh_kin, Eh_pot, Ep, Nsiz, SUM(MC(1)%MC_Electrons(:)%Ekin, MASK = MC(1)%MC_Electrons(:)%active) ! electrons
end subroutine get_total_values


subroutine sort_data(used_target, MC, numpar, tim, Spectrum_ph, Spectrum_e, Spectrum_p, Spectrum_h, Spectrum_SHI, Spectrum_mu, &
    Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X, Spectra_mu_X, &
    Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y, Spectra_mu_Y, &
    Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z, Spectra_mu_Z, &
    Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R, Spectra_mu_R, &
    Vel_theta_ph, Vel_theta_e, Vel_theta_p, Vel_theta_h, Vel_theta_SHI, Vel_theta_mu, &
    Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X, Theta_mu_X, &
    Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y, Theta_mu_Y, &
    Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z, Theta_mu_Z, &
    Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R, Theta_mu_R, &
    Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_ph_R, Distr_ph_L, Distr_ph_Theta, Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic, &
    Distr_e_X, Distr_e_Y, Distr_e_Z, Distr_e_R, Distr_e_L, Distr_e_Theta, Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic, &
    Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_p_R, Distr_p_L, Distr_p_Theta, Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic, &
    Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_h_R, Distr_h_L, Distr_h_Theta, Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic, &
    Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, Distr_SHI_R, Distr_SHI_L, Distr_SHI_Theta, Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic, &
    Distr_a_X, Distr_a_Y, Distr_a_Z, Distr_a_R, Distr_a_L, Distr_a_Theta, Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic, &
    Distr_mu_X, Distr_mu_Y, Distr_mu_Z, Distr_mu_R, Distr_mu_L, Distr_mu_Theta, Distr_mu_Rc, Distr_mu_Thetac, Distr_mu_Phic, &
    Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ, Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta, &
    Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic, &
    Distr_e_XY, Distr_e_YZ, Distr_e_XZ, Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta, &
    Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic, &
    Distr_p_XY, Distr_p_YZ, Distr_p_XZ, Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta, &
    Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic, &
    Distr_h_XY, Distr_h_YZ, Distr_h_XZ, Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta, &
    Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic, &
    Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ, Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta, &
    Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic, &
    Distr_a_XY, Distr_a_YZ, Distr_a_XZ, Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta, &
    Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic, &
    Distr_mu_XY, Distr_mu_YZ, Distr_mu_XZ, Distr_mu_RL, Distr_mu_RTheta, Distr_mu_LTheta, &
    Distr_mu_RcThc, Distr_mu_RcPhic, Distr_mu_ThcPhic, &
    Distr_ph_XYZ, Distr_ph_RLTheta, Distr_ph_RcThcPhic, &
    Distr_e_XYZ, Distr_e_RLTheta, Distr_e_RcThcPhic, &
    Distr_p_XYZ, Distr_p_RLTheta, Distr_p_RcThcPhic, &
    Distr_h_XYZ, Distr_h_RLTheta, Distr_h_RcThcPhic, &
    Distr_SHI_XYZ, Distr_SHI_RLTheta, Distr_SHI_RcThcPhic, & 
    Distr_a_XYZ, Distr_a_RLTheta, Distr_a_RcThcPhic, & 
    Distr_mu_XYZ, Distr_mu_RLTheta, Distr_mu_RcThcPhic, &
    E_Distr_ph_X, E_Distr_ph_Y, E_Distr_ph_Z, E_Distr_ph_R, E_Distr_ph_L, &
    E_Distr_ph_Theta, E_Distr_ph_Rc, E_Distr_ph_Thetac, E_Distr_ph_Phic, &
    E_Distr_e_X, E_Distr_e_Y, E_Distr_e_Z, E_Distr_e_R, E_Distr_e_L, E_Distr_e_Theta, E_Distr_e_Rc, E_Distr_e_Thetac, E_Distr_e_Phic, &
    E_Distr_p_X, E_Distr_p_Y, E_Distr_p_Z, E_Distr_p_R, E_Distr_p_L, E_Distr_p_Theta, E_Distr_p_Rc, E_Distr_p_Thetac, E_Distr_p_Phic, &
    E_Distr_h_X, E_Distr_h_Y, E_Distr_h_Z, E_Distr_h_R, E_Distr_h_L, E_Distr_h_Theta, E_Distr_h_Rc, E_Distr_h_Thetac, E_Distr_h_Phic, &
    E_Distr_SHI_X, E_Distr_SHI_Y, E_Distr_SHI_Z, E_Distr_SHI_R, E_Distr_SHI_L, &
    E_Distr_SHI_Theta, E_Distr_SHI_Rc, E_Distr_SHI_Thetac, E_Distr_SHI_Phic, &
    E_Distr_a_X, E_Distr_a_Y, E_Distr_a_Z, E_Distr_a_R, E_Distr_a_L, E_Distr_a_Theta, E_Distr_a_Rc, E_Distr_a_Thetac, E_Distr_a_Phic, &
    E_Distr_mu_X, E_Distr_mu_Y, E_Distr_mu_Z, E_Distr_mu_R, E_Distr_mu_L, E_Distr_mu_Theta, E_Distr_mu_Rc, E_Distr_mu_Thetac, E_Distr_mu_Phic, &
    E_Distr_ph_XY, E_Distr_ph_YZ, E_Distr_ph_XZ, E_Distr_ph_RL, E_Distr_ph_RTheta, E_Distr_ph_LTheta, &
    E_Distr_ph_RcThc, E_Distr_ph_RcPhic, E_Distr_ph_ThcPhic, &
    E_Distr_e_XY, E_Distr_e_YZ, E_Distr_e_XZ, E_Distr_e_RL, E_Distr_e_RTheta, E_Distr_e_LTheta, &
    E_Distr_e_RcThc, E_Distr_e_RcPhic, E_Distr_e_ThcPhic, &
    E_Distr_p_XY, E_Distr_p_YZ, E_Distr_p_XZ, E_Distr_p_RL, E_Distr_p_RTheta, E_Distr_p_LTheta, &
    E_Distr_p_RcThc, E_Distr_p_RcPhic, E_Distr_p_ThcPhic, &
    E_Distr_h_XY, E_Distr_h_YZ, E_Distr_h_XZ, E_Distr_h_RL, E_Distr_h_RTheta, E_Distr_h_LTheta, &
    E_Distr_h_RcThc, E_Distr_h_RcPhic, E_Distr_h_ThcPhic, &
    E_Distr_SHI_XY, E_Distr_SHI_YZ, E_Distr_SHI_XZ, E_Distr_SHI_RL, E_Distr_SHI_RTheta, E_Distr_SHI_LTheta, &
    E_Distr_SHI_RcThc, E_Distr_SHI_RcPhic, E_Distr_SHI_ThcPhic, &
    E_Distr_a_XY, E_Distr_a_YZ, E_Distr_a_XZ, E_Distr_a_RL, E_Distr_a_RTheta, E_Distr_a_LTheta, &
    E_Distr_a_RcThc, E_Distr_a_RcPhic, E_Distr_a_ThcPhic, &
    E_Distr_mu_XY, E_Distr_mu_YZ, E_Distr_mu_XZ, E_Distr_mu_RL, E_Distr_mu_RTheta, E_Distr_mu_LTheta, &
    E_Distr_mu_RcThc, E_Distr_mu_RcPhic, E_Distr_mu_ThcPhic, &
    E_Distr_ph_XYZ, E_Distr_ph_RLTheta, E_Distr_ph_RcThcPhic, &
    E_Distr_e_XYZ, E_Distr_e_RLTheta, E_Distr_e_RcThcPhic, &
    E_Distr_p_XYZ, E_Distr_p_RLTheta, E_Distr_p_RcThcPhic, &
    E_Distr_h_XYZ, E_Distr_h_RLTheta, E_Distr_h_RcThcPhic, &
    E_Distr_SHI_XYZ, E_Distr_SHI_RLTheta, E_Distr_SHI_RcThcPhic, &
    E_Distr_a_XYZ, E_Distr_a_RLTheta, E_Distr_a_RcThcPhic, &
    E_Distr_mu_XYZ, E_Distr_mu_RLTheta, E_Distr_mu_RcThcPhic )
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(MC_arrays), dimension(:), intent(in) :: MC	! all MC arrays for all particles; size = number of iterations
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   ! Energy disributions (spectra):
   real(8), dimension(:), intent(inout) :: Spectrum_ph, Spectrum_e, Spectrum_p, Spectrum_SHI, Spectrum_mu  ! energy spectra
   real(8), dimension(:,:), intent(inout) :: Spectrum_h ! VB holes spectra for each material
   real(8), dimension(:), intent(inout) :: Vel_theta_ph, Vel_theta_e, Vel_theta_p, Vel_theta_h, Vel_theta_SHI, Vel_theta_mu   ! velosity theta distribution
   ! Spectra in 1d space:
   real(8), dimension(:,:), intent(inout) :: Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X, Spectra_mu_X  ! energy spectra in space along X
   real(8), dimension(:,:), intent(inout) :: Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y, Spectra_mu_Y  ! energy spectra in space along Y
   real(8), dimension(:,:), intent(inout) :: Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z, Spectra_mu_Z  ! energy spectra in space along Z
   real(8), dimension(:,:), intent(inout) :: Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R, Spectra_mu_R  ! energy spectra in space along R
   ! Theta distribution in 1d space:
   real(8), dimension(:,:), intent(inout) :: Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X, Theta_mu_X ! theta distr in space along X
   real(8), dimension(:,:), intent(inout) :: Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y, Theta_mu_Y ! theta distr in space along Y
   real(8), dimension(:,:), intent(inout) :: Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z, Theta_mu_Z ! theta distr in space along Z
   real(8), dimension(:,:), intent(inout) :: Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R, Theta_mu_R ! theta distr in space along R

   ! Spatial distributions in 1d:
   real(8), dimension(:), intent(inout) :: Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_ph_R, Distr_ph_L, Distr_ph_Theta, &
                                                        Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic ! photon
   real(8), dimension(:), intent(inout) :: Distr_e_X, Distr_e_Y, Distr_e_Z, Distr_e_R, Distr_e_L, Distr_e_Theta, Distr_e_Rc, &
                                                        Distr_e_Thetac, Distr_e_Phic ! electron
   real(8), dimension(:), intent(inout) :: Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_p_R, Distr_p_L, Distr_p_Theta, Distr_p_Rc, &
                                                        Distr_p_Thetac, Distr_p_Phic ! positron
   real(8), dimension(:,:), intent(inout) :: Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_h_R, Distr_h_L, Distr_h_Theta, &
                                                        Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic ! hole
   real(8), dimension(:), intent(inout) :: Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, Distr_SHI_R, Distr_SHI_L, Distr_SHI_Theta, &
                                                        Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic ! SHI
   real(8), dimension(:), intent(inout) :: Distr_a_X, Distr_a_Y, Distr_a_Z, Distr_a_R, Distr_a_L, Distr_a_Theta, &
                                                        Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic ! Atom
   real(8), dimension(:), intent(inout) :: Distr_mu_X, Distr_mu_Y, Distr_mu_Z, Distr_mu_R, Distr_mu_L, Distr_mu_Theta, Distr_mu_Rc, &
                                                        Distr_mu_Thetac, Distr_mu_Phic ! muon
   ! Spatial distributions in 2d:
   real(8), dimension(:,:), intent(inout) :: Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ, Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta
   real(8), dimension(:,:), intent(inout) :: Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic  ! photon
   real(8), dimension(:,:), intent(inout) :: Distr_e_XY, Distr_e_YZ, Distr_e_XZ, Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta
   real(8), dimension(:,:), intent(inout) :: Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic  ! electron
   real(8), dimension(:,:), intent(inout) :: Distr_p_XY, Distr_p_YZ, Distr_p_XZ, Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta
   real(8), dimension(:,:), intent(inout) :: Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic  ! positron
   real(8), dimension(:,:,:), intent(inout) :: Distr_h_XY, Distr_h_YZ, Distr_h_XZ, Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta
   real(8), dimension(:,:,:), intent(inout) :: Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic  ! hole
   real(8), dimension(:,:), intent(inout) :: Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ, Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta
   real(8), dimension(:,:), intent(inout) :: Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic  ! SHI
   real(8), dimension(:,:), intent(inout) :: Distr_a_XY, Distr_a_YZ, Distr_a_XZ, Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta
   real(8), dimension(:,:), intent(inout) :: Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic  ! Atom
   real(8), dimension(:,:), intent(inout) :: Distr_mu_XY, Distr_mu_YZ, Distr_mu_XZ, Distr_mu_RL, Distr_mu_RTheta, Distr_mu_LTheta, &
                                             Distr_mu_RcThc, Distr_mu_RcPhic, Distr_mu_ThcPhic  ! muon
   ! Spatial distributions in 3d:
   real(8), dimension(:,:,:), intent(inout) :: Distr_ph_XYZ, Distr_ph_RLTheta, Distr_ph_RcThcPhic ! photon
   real(8), dimension(:,:,:), intent(inout) :: Distr_e_XYZ, Distr_e_RLTheta, Distr_e_RcThcPhic ! electron
   real(8), dimension(:,:,:), intent(inout) :: Distr_p_XYZ, Distr_p_RLTheta, Distr_p_RcThcPhic ! positron
   real(8), dimension(:,:,:,:), intent(inout) :: Distr_h_XYZ, Distr_h_RLTheta, Distr_h_RcThcPhic ! hole
   real(8), dimension(:,:,:), intent(inout) :: Distr_SHI_XYZ, Distr_SHI_RLTheta, Distr_SHI_RcThcPhic ! SHI
   real(8), dimension(:,:,:), intent(inout) :: Distr_a_XYZ, Distr_a_RLTheta, Distr_a_RcThcPhic ! Atom
   real(8), dimension(:,:,:), intent(inout) :: Distr_mu_XYZ, Distr_mu_RLTheta, Distr_mu_RcThcPhic ! muon

   ! Spatial energy distributions in 1d:
   real(8), dimension(:), intent(inout) :: E_Distr_ph_X, E_Distr_ph_Y, E_Distr_ph_Z, E_Distr_ph_R, E_Distr_ph_L, E_Distr_ph_Theta
   real(8), dimension(:), intent(inout):: E_Distr_ph_Rc, E_Distr_ph_Thetac, E_Distr_ph_Phic ! photon
   real(8), dimension(:), intent(inout) :: E_Distr_e_X, E_Distr_e_Y, E_Distr_e_Z, E_Distr_e_R, E_Distr_e_L, E_Distr_e_Theta
   real(8), dimension(:), intent(inout) :: E_Distr_e_Rc, E_Distr_e_Thetac, E_Distr_e_Phic ! electron
   real(8), dimension(:), intent(inout) :: E_Distr_p_X, E_Distr_p_Y, E_Distr_p_Z, E_Distr_p_R, E_Distr_p_L, E_Distr_p_Theta
   real(8), dimension(:), intent(inout) :: E_Distr_p_Rc, E_Distr_p_Thetac, E_Distr_p_Phic ! positron
   real(8), dimension(:,:), intent(inout) :: E_Distr_h_X, E_Distr_h_Y, E_Distr_h_Z, E_Distr_h_R, E_Distr_h_L, E_Distr_h_Theta
   real(8), dimension(:,:), intent(inout) :: E_Distr_h_Rc, E_Distr_h_Thetac, E_Distr_h_Phic ! hole
   real(8), dimension(:), intent(inout) :: E_Distr_SHI_X, E_Distr_SHI_Y, E_Distr_SHI_Z, E_Distr_SHI_R, E_Distr_SHI_L
   real(8), dimension(:), intent(inout) :: E_Distr_SHI_Theta, E_Distr_SHI_Rc, E_Distr_SHI_Thetac, E_Distr_SHI_Phic ! SHI
   real(8), dimension(:), intent(inout) :: E_Distr_a_X, E_Distr_a_Y, E_Distr_a_Z, E_Distr_a_R, E_Distr_a_L
   real(8), dimension(:), intent(inout) :: E_Distr_a_Theta, E_Distr_a_Rc, E_Distr_a_Thetac, E_Distr_a_Phic ! Atom
   real(8), dimension(:), intent(inout) :: E_Distr_mu_X, E_Distr_mu_Y, E_Distr_mu_Z, E_Distr_mu_R, E_Distr_mu_L, E_Distr_mu_Theta, &
                                           E_Distr_mu_Rc, E_Distr_mu_Thetac, E_Distr_mu_Phic ! muon
   ! Spatial energy distributions in 2d:
   real(8), dimension(:,:), intent(inout) :: E_Distr_ph_XY, E_Distr_ph_YZ, E_Distr_ph_XZ, E_Distr_ph_RL, E_Distr_ph_RTheta, E_Distr_ph_LTheta
   real(8), dimension(:,:), intent(inout) :: E_Distr_ph_RcThc, E_Distr_ph_RcPhic, E_Distr_ph_ThcPhic  ! photon
   real(8), dimension(:,:), intent(inout) :: E_Distr_e_XY, E_Distr_e_YZ, E_Distr_e_XZ, E_Distr_e_RL, E_Distr_e_RTheta, E_Distr_e_LTheta
   real(8), dimension(:,:), intent(inout) :: E_Distr_e_RcThc, E_Distr_e_RcPhic, E_Distr_e_ThcPhic  ! electron
   real(8), dimension(:,:), intent(inout) :: E_Distr_p_XY, E_Distr_p_YZ, E_Distr_p_XZ, E_Distr_p_RL, E_Distr_p_RTheta, E_Distr_p_LTheta
   real(8), dimension(:,:), intent(inout) :: E_Distr_p_RcThc, E_Distr_p_RcPhic, E_Distr_p_ThcPhic  ! positron
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_h_XY, E_Distr_h_YZ, E_Distr_h_XZ, E_Distr_h_RL, E_Distr_h_RTheta, E_Distr_h_LTheta
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_h_RcThc, E_Distr_h_RcPhic, E_Distr_h_ThcPhic  ! hole
   real(8), dimension(:,:), intent(inout) :: E_Distr_SHI_XY, E_Distr_SHI_YZ, E_Distr_SHI_XZ, E_Distr_SHI_RL, E_Distr_SHI_RTheta, E_Distr_SHI_LTheta
   real(8), dimension(:,:), intent(inout) :: E_Distr_SHI_RcThc, E_Distr_SHI_RcPhic, E_Distr_SHI_ThcPhic  ! SHI
   real(8), dimension(:,:), intent(inout) :: E_Distr_a_XY, E_Distr_a_YZ, E_Distr_a_XZ, E_Distr_a_RL, E_Distr_a_RTheta, E_Distr_a_LTheta
   real(8), dimension(:,:), intent(inout) :: E_Distr_a_RcThc, E_Distr_a_RcPhic, E_Distr_a_ThcPhic  ! Atom
   real(8), dimension(:,:), intent(inout) :: E_Distr_mu_XY, E_Distr_mu_YZ, E_Distr_mu_XZ, E_Distr_mu_RL, E_Distr_mu_RTheta, E_Distr_mu_LTheta, &
                                             E_Distr_mu_RcThc, E_Distr_mu_RcPhic, E_Distr_mu_ThcPhic  ! muon
   ! Spatial energy distributions in 3d:
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_ph_XYZ, E_Distr_ph_RLTheta, E_Distr_ph_RcThcPhic ! photon
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_e_XYZ, E_Distr_e_RLTheta, E_Distr_e_RcThcPhic ! electron
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_p_XYZ, E_Distr_p_RLTheta, E_Distr_p_RcThcPhic ! positron
   real(8), dimension(:,:,:,:), intent(inout) :: E_Distr_h_XYZ, E_Distr_h_RLTheta, E_Distr_h_RcThcPhic ! hole
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_SHI_XYZ, E_Distr_SHI_RLTheta, E_Distr_SHI_RcThcPhic ! SHI
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_a_XYZ, E_Distr_a_RLTheta, E_Distr_a_RcThcPhic ! Atom
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_mu_XYZ, E_Distr_mu_RLTheta, E_Distr_mu_RcThcPhic ! muon
   !-------------------------------------------------------
   integer :: iter


!  print*, 'sort_data 0'

   
   !$omp parallel  &
   !$omp private (iter)

! goto 9999   ! Test

   ! Spectra:
   !$omp do schedule(dynamic) reduction( + : Spectrum_ph, Spectrum_e, Spectrum_p, Spectrum_h, Spectrum_SHI, Spectrum_mu)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
      ! Get photon distributions:
      call get_photon_spectrum(MC(iter), numpar, Spectrum_ph)   ! below
      
      ! Get electron distributions:
      call get_electron_spectrum(MC(iter), numpar, Spectrum_e)   ! below
      
      ! Get positron distributions:
      call get_positron_spectrum(MC(iter), numpar, Spectrum_p)   ! below
      
      ! Get holes distributions:
      call get_hole_spectrum(MC(iter), numpar, Spectrum_h)   ! below
      
      ! Get SHI distributions:
      call get_SHI_spectrum(MC(iter), numpar, Spectrum_SHI)   ! below

      ! Get muon distributions:
      call get_muon_spectrum(MC(iter), numpar, Spectrum_mu)   ! below
   enddo
   !$omp enddo
   !$OMP BARRIER

! goto 9999   ! Test WORKED

   !--------------------------------
   ! Get velosity theta distribution:
   !$omp do schedule(dynamic) reduction( + : Vel_theta_ph, Vel_theta_e, Vel_theta_p, Vel_theta_h, Vel_theta_SHI)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
      ! Get photon velosity theta distributions:
      call get_photon_velotheta(MC(iter), numpar, Vel_theta_ph)   ! below
      
      ! Get electron velosity theta distributions:
      call get_electron_velotheta(MC(iter), numpar, Vel_theta_e)   ! below
      
      ! Get positron velosity theta distributions:
      call get_positron_velotheta(MC(iter), numpar, Vel_theta_p)   ! below
      
      ! Get holes velosity theta distributions:
      call get_hole_velotheta(MC(iter), numpar, Vel_theta_h)   ! below
      
      ! Get SHI velosity theta distributions:
      call get_SHI_velotheta(MC(iter), numpar, Vel_theta_SHI)   ! below
   enddo
   !$omp enddo
   !$OMP BARRIER
   
!    print*, 'sort_data 2', iter
! goto 9999   ! Test FAILED, WORKED 2d time ???

   !--------------------------------
   ! Get spectra vs space along all axes that user requested:
   ! Spectra vs space in 1d:
   !$omp do schedule(dynamic) reduction( + : Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X, &
   !$omp    Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y, &
   !$omp    Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z, &
   !$omp    Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
      call get_spectra_in_space_1d(MC(iter), numpar, tim, Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X, &
                                                    Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y, &
                                                    Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z, &
                                                    Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R) ! below
!        print*, 'sort_data 1.5', iter
   enddo
   !$omp enddo
   !$OMP BARRIER


   !--------------------------------
   ! Get theta distribution vs space along all axes that user requested:
   ! Theta distribution  vs space in 1d:
   !$omp do schedule(dynamic) reduction( + : Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X, &
   !$omp    Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y, &
   !$omp    Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z, &
   !$omp    Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
      call get_theta_in_space_1d(MC(iter), numpar, tim, Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X, &
                                                    Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y, &
                                                    Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z, &
                                                    Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R) ! below
!        print*, 'sort_theta 1.5', iter
   enddo
   !$omp enddo
   !$OMP BARRIER
   
!    print*, 'sort_data 3', iter

   !--------------------------------
   ! Get spatial distributions along all axes that user requested:
   ! Spatial distributions in 1d:
   
   ! Cartesian:
   !$omp do schedule(dynamic) reduction( + : Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_e_X, Distr_e_Y, &
   !$omp  Distr_e_Z, Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, &
   !$omp  Distr_a_X, Distr_a_Y, Distr_a_Z, &
   !$omp  E_Distr_ph_X, E_Distr_ph_Y, E_Distr_ph_Z, E_Distr_e_X, E_Distr_e_Y, &
   !$omp  E_Distr_e_Z, E_Distr_p_X, E_Distr_p_Y, E_Distr_p_Z, E_Distr_h_X, E_Distr_h_Y, E_Distr_h_Z, E_Distr_SHI_X, E_Distr_SHI_Y, E_Distr_SHI_Z, &
   !$omp  E_Distr_a_X, E_Distr_a_Y, E_Distr_a_Z)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
       call get_distribution_1d_Cartesian(used_target, MC(iter), numpar, tim, Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_e_X, Distr_e_Y, &
                  Distr_e_Z, Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, &
                  Distr_a_X, Distr_a_Y, Distr_a_Z, E_Distr_ph_X, E_Distr_ph_Y, E_Distr_ph_Z, E_Distr_e_X, E_Distr_e_Y, &
                  E_Distr_e_Z, E_Distr_p_X, E_Distr_p_Y, E_Distr_p_Z, E_Distr_h_X, E_Distr_h_Y, E_Distr_h_Z, &
                  E_Distr_SHI_X, E_Distr_SHI_Y, E_Distr_SHI_Z, E_Distr_a_X, E_Distr_a_Y, E_Distr_a_Z)    ! below
   enddo
   !$omp end do
   !$OMP BARRIER
   
!    print*, 'sort_data 4', iter

   ! Cylindrical:
   !$omp do schedule(dynamic) reduction( + : Distr_ph_R, Distr_ph_L, Distr_ph_Theta, Distr_e_R, Distr_e_L, Distr_e_Theta, &
   !$omp Distr_p_R, Distr_p_L, Distr_p_Theta, Distr_h_R, Distr_h_L, Distr_h_Theta, Distr_SHI_R, Distr_SHI_L, Distr_SHI_Theta, &
   !$omp Distr_a_R, Distr_a_L, Distr_a_Theta, &
   !$omp E_Distr_ph_R, E_Distr_ph_L, E_Distr_ph_Theta, E_Distr_e_R, E_Distr_e_L, E_Distr_e_Theta, &
   !$omp E_Distr_p_R, E_Distr_p_L, E_Distr_p_Theta, E_Distr_h_R, E_Distr_h_L, E_Distr_h_Theta, &
   !$omp E_Distr_SHI_R, E_Distr_SHI_L, E_Distr_SHI_Theta, E_Distr_a_R, E_Distr_a_L, E_Distr_a_Theta)
   do iter = 1, numpar%NMC  ! analyze data in all iterations   
      call get_distribution_1d_Cylindrical(used_target, MC(iter), numpar, tim, Distr_ph_R, Distr_ph_L, Distr_ph_Theta, Distr_e_R, Distr_e_L, Distr_e_Theta, &
                 Distr_p_R, Distr_p_L, Distr_p_Theta, Distr_h_R, Distr_h_L, Distr_h_Theta, Distr_SHI_R, Distr_SHI_L, Distr_SHI_Theta, &
                 Distr_a_R, Distr_a_L, Distr_a_Theta, &
                 E_Distr_ph_R, E_Distr_ph_L, E_Distr_ph_Theta, E_Distr_e_R, E_Distr_e_L, E_Distr_e_Theta, &
                 E_Distr_p_R, E_Distr_p_L, E_Distr_p_Theta, E_Distr_h_R, E_Distr_h_L, E_Distr_h_Theta, &
                 E_Distr_SHI_R, E_Distr_SHI_L, E_Distr_SHI_Theta, E_Distr_a_R, E_Distr_a_L, E_Distr_a_Theta) ! below
   enddo
   !$omp end do
   !$OMP BARRIER
   
!    print*, 'sort_data 5', iter

   ! Spherical:
   !$omp do schedule(dynamic) reduction( + : Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic, Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic, &
   !$omp Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic,  Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic, Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic, &
   !$omp Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic, &
   !$omp E_Distr_ph_Rc, E_Distr_ph_Thetac, E_Distr_ph_Phic, E_Distr_e_Rc, E_Distr_e_Thetac, E_Distr_e_Phic, &
   !$omp E_Distr_p_Rc, E_Distr_p_Thetac, E_Distr_p_Phic,  E_Distr_h_Rc, E_Distr_h_Thetac, E_Distr_h_Phic, &
   !$omp E_Distr_SHI_Rc, E_Distr_SHI_Thetac, E_Distr_SHI_Phic, E_Distr_a_Rc, E_Distr_a_Thetac, E_Distr_a_Phic)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
      call get_distribution_1d_Spherical(used_target, MC(iter), numpar, tim, Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic, Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic, &
            Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic,  Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic, Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic, &
            Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic, &
            E_Distr_ph_Rc, E_Distr_ph_Thetac, E_Distr_ph_Phic, E_Distr_e_Rc, E_Distr_e_Thetac, E_Distr_e_Phic, &
            E_Distr_p_Rc, E_Distr_p_Thetac, E_Distr_p_Phic,  E_Distr_h_Rc, E_Distr_h_Thetac, E_Distr_h_Phic, &
            E_Distr_SHI_Rc, E_Distr_SHI_Thetac, E_Distr_SHI_Phic, E_Distr_a_Rc, E_Distr_a_Thetac, E_Distr_a_Phic)   ! below
   enddo
   !$omp end do
   !$OMP BARRIER
   
!    print*, 'sort_data 6', iter

   ! Spatial distributions in 2d:
   
   ! Cartesian:
   !$omp do schedule(dynamic) reduction( + : Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ, Distr_e_XY, Distr_e_YZ, Distr_e_XZ, &
   !$omp Distr_p_XY, Distr_p_YZ, Distr_p_XZ,  Distr_h_XY, Distr_h_YZ, Distr_h_XZ, Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ, &
   !$omp Distr_a_XY, Distr_a_YZ, Distr_a_XZ, &
   !$omp E_Distr_ph_XY, E_Distr_ph_YZ, E_Distr_ph_XZ, E_Distr_e_XY, E_Distr_e_YZ, E_Distr_e_XZ, &
   !$omp E_Distr_p_XY, E_Distr_p_YZ, E_Distr_p_XZ, E_Distr_h_XY, E_Distr_h_YZ, E_Distr_h_XZ, &
   !$omp E_Distr_SHI_XY, E_Distr_SHI_YZ, E_Distr_SHI_XZ, E_Distr_a_XY, E_Distr_a_YZ, E_Distr_a_XZ)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
      call get_distribution_2d_Cartesian(used_target, MC(iter), numpar, tim, Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ, Distr_e_XY, Distr_e_YZ, Distr_e_XZ, &
                 Distr_p_XY, Distr_p_YZ, Distr_p_XZ,  Distr_h_XY, Distr_h_YZ, Distr_h_XZ, Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ, &
                 Distr_a_XY, Distr_a_YZ, Distr_a_XZ, &
                 E_Distr_ph_XY, E_Distr_ph_YZ, E_Distr_ph_XZ, E_Distr_e_XY, E_Distr_e_YZ, E_Distr_e_XZ, &
                 E_Distr_p_XY, E_Distr_p_YZ, E_Distr_p_XZ, E_Distr_h_XY, E_Distr_h_YZ, E_Distr_h_XZ, &
                 E_Distr_SHI_XY, E_Distr_SHI_YZ, E_Distr_SHI_XZ, E_Distr_a_XY, E_Distr_a_YZ, E_Distr_a_XZ) ! below
   enddo
   !$omp end do
   !$OMP BARRIER
   
   ! Cylindrical:
   !$omp do schedule(dynamic) reduction( + : Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta,  Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta, &
   !$omp Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta, Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta, Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta, &
   !$omp Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta,&
   !$omp E_Distr_ph_RL, E_Distr_ph_RTheta, E_Distr_ph_LTheta, E_Distr_e_RL, E_Distr_e_RTheta, E_Distr_e_LTheta, &
   !$omp E_Distr_p_RL, E_Distr_p_RTheta, E_Distr_p_LTheta, E_Distr_h_RL, E_Distr_h_RTheta, E_Distr_h_LTheta, &
   !$omp E_Distr_SHI_RL, E_Distr_SHI_RTheta, E_Distr_SHI_LTheta, E_Distr_a_RL, E_Distr_a_RTheta, E_Distr_a_LTheta)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
      call get_distribution_2d_Cylindrical(used_target, MC(iter), numpar, tim, Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta,  Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta, &
                 Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta, Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta, Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta, &
                 Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta, &
                 E_Distr_ph_RL, E_Distr_ph_RTheta, E_Distr_ph_LTheta, E_Distr_e_RL, E_Distr_e_RTheta, E_Distr_e_LTheta, &
                 E_Distr_p_RL, E_Distr_p_RTheta, E_Distr_p_LTheta, E_Distr_h_RL, E_Distr_h_RTheta, E_Distr_h_LTheta, &
                 E_Distr_SHI_RL, E_Distr_SHI_RTheta, E_Distr_SHI_LTheta, E_Distr_a_RL, E_Distr_a_RTheta, E_Distr_a_LTheta)    ! below
   enddo
   !$omp end do
   !$OMP BARRIER
   
   ! Spherical:
   !$omp do schedule(dynamic) reduction( + : Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic, &
   !$omp Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic, Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic, &
   !$omp Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic,  Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic, &
   !$omp Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic, &
   !$omp E_Distr_ph_RcThc, E_Distr_ph_RcPhic, E_Distr_ph_ThcPhic, &
   !$omp E_Distr_e_RcThc, E_Distr_e_RcPhic, E_Distr_e_ThcPhic, E_Distr_p_RcThc, E_Distr_p_RcPhic, E_Distr_p_ThcPhic, &
   !$omp E_Distr_h_RcThc, E_Distr_h_RcPhic, E_Distr_h_ThcPhic, E_Distr_SHI_RcThc, E_Distr_SHI_RcPhic, E_Distr_SHI_ThcPhic,&
   !$omp E_Distr_a_RcThc, E_Distr_a_RcPhic, E_Distr_a_ThcPhic)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
      call get_distribution_2d_Spherical(used_target, MC(iter), numpar, tim, Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic, &
                 Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic, Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic, &
                 Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic,  Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic, &
                 Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic, &
                 E_Distr_ph_RcThc, E_Distr_ph_RcPhic, E_Distr_ph_ThcPhic, &
                 E_Distr_e_RcThc, E_Distr_e_RcPhic, E_Distr_e_ThcPhic, E_Distr_p_RcThc, E_Distr_p_RcPhic, E_Distr_p_ThcPhic, &
                 E_Distr_h_RcThc, E_Distr_h_RcPhic, E_Distr_h_ThcPhic, E_Distr_SHI_RcThc, E_Distr_SHI_RcPhic, E_Distr_SHI_ThcPhic, &
                 E_Distr_a_RcThc, E_Distr_a_RcPhic, E_Distr_a_ThcPhic) ! below
   enddo
   !$omp end do
   !$OMP BARRIER
   
   ! Spatial distributions in 3d:
   
   ! Cartesian:
   !$omp do schedule(dynamic) reduction( + : Distr_ph_XYZ, Distr_e_XYZ, Distr_p_XYZ, Distr_h_XYZ, Distr_SHI_XYZ, Distr_a_XYZ,&
   !$omp E_Distr_ph_XYZ, E_Distr_e_XYZ, E_Distr_p_XYZ, E_Distr_h_XYZ, E_Distr_SHI_XYZ, E_Distr_a_XYZ)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
      call get_distribution_3d_Cartesian(used_target, MC(iter), numpar, tim, &
            Distr_ph_XYZ, Distr_e_XYZ, Distr_p_XYZ, Distr_h_XYZ, Distr_SHI_XYZ, Distr_a_XYZ, &
            E_Distr_ph_XYZ, E_Distr_e_XYZ, E_Distr_p_XYZ, E_Distr_h_XYZ, E_Distr_SHI_XYZ, E_Distr_a_XYZ)  ! Below
   enddo
   !$omp end do
   
   ! Cylindrical:
   !$omp do schedule(dynamic) reduction( + : Distr_ph_RLTheta, Distr_e_RLTheta, Distr_p_RLTheta, Distr_h_RLTheta, Distr_SHI_RLTheta, Distr_a_RLTheta, &
   !$omp E_Distr_ph_RLTheta, E_Distr_e_RLTheta, E_Distr_p_RLTheta, E_Distr_h_RLTheta, E_Distr_SHI_RLTheta, E_Distr_a_RLTheta)
   do iter = 1, numpar%NMC  ! analyze data in all iterations   
      call get_distribution_3d_Cylindrical(used_target, MC(iter), numpar, tim, Distr_ph_RLTheta, Distr_e_RLTheta, &
                 Distr_p_RLTheta, Distr_h_RLTheta, Distr_SHI_RLTheta, Distr_a_RLTheta, &
                 E_Distr_ph_RLTheta, E_Distr_e_RLTheta, E_Distr_p_RLTheta, E_Distr_h_RLTheta, E_Distr_SHI_RLTheta, E_Distr_a_RLTheta)    ! Below
   enddo
   !$omp end do
   !$OMP BARRIER
   
!goto 9999   ! Test FAILED

   ! Spherical:
   !$omp do schedule(dynamic) reduction( + : Distr_ph_RcThcPhic, Distr_e_RcThcPhic, Distr_p_RcThcPhic, Distr_h_RcThcPhic, &
   !$omp Distr_SHI_RcThcPhic, Distr_a_RcThcPhic, &
   !$omp E_Distr_ph_RcThcPhic, E_Distr_e_RcThcPhic, E_Distr_p_RcThcPhic, E_Distr_h_RcThcPhic, E_Distr_SHI_RcThcPhic, E_Distr_a_RcThcPhic)
   do iter = 1, numpar%NMC  ! analyze data in all iterations
      call get_distribution_3d_Spherical(used_target, MC(iter), numpar, tim, Distr_ph_RcThcPhic, Distr_e_RcThcPhic, Distr_p_RcThcPhic, &
                    Distr_h_RcThcPhic, Distr_SHI_RcThcPhic, Distr_a_RcThcPhic, &
                    E_Distr_ph_RcThcPhic, E_Distr_e_RcThcPhic, E_Distr_p_RcThcPhic, E_Distr_h_RcThcPhic, &
                    E_Distr_SHI_RcThcPhic, E_Distr_a_RcThcPhic)  ! Below
   enddo
   !$omp end do
   !$OMP BARRIER
   
! 9999 continue

   !$omp end parallel

!  print*, 'sort_data END'
end subroutine sort_data



!ddddddddddddddddddddddddddddddddddddddddddddddddddddddd
! Distributions in space:
 ! Spatial distributions in 1d:
subroutine get_distribution_1d_Cartesian(used_target, MC, numpar, tim, Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_e_X, Distr_e_Y, &
                 Distr_e_Z, Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, &
                 Distr_a_X, Distr_a_Y, Distr_a_Z, &
                 E_Distr_ph_X, E_Distr_ph_Y, E_Distr_ph_Z, E_Distr_e_X, E_Distr_e_Y, &
                 E_Distr_e_Z, E_Distr_p_X, E_Distr_p_Y, E_Distr_p_Z, E_Distr_h_X, E_Distr_h_Y, E_Distr_h_Z, E_Distr_SHI_X, E_Distr_SHI_Y, E_Distr_SHI_Z, &
                 E_Distr_a_X, E_Distr_a_Y, E_Distr_a_Z)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   ! Density distributions:
   real(8), dimension(:), intent(inout) :: Distr_ph_X, Distr_ph_Y, Distr_ph_Z   ! photon
   real(8), dimension(:), intent(inout) :: Distr_e_X, Distr_e_Y, Distr_e_Z  ! electron
   real(8), dimension(:), intent(inout) :: Distr_p_X, Distr_p_Y, Distr_p_Z  ! positron
   real(8), dimension(:,:), intent(inout) :: Distr_h_X, Distr_h_Y, Distr_h_Z  ! hole
   real(8), dimension(:), intent(inout) :: Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z    ! SHI
   real(8), dimension(:), intent(inout) :: Distr_a_X, Distr_a_Y, Distr_a_Z    ! Atom
   ! Energy distributions:
   real(8), dimension(:), intent(inout) :: E_Distr_ph_X, E_Distr_ph_Y, E_Distr_ph_Z   ! photon
   real(8), dimension(:), intent(inout) :: E_Distr_e_X, E_Distr_e_Y, E_Distr_e_Z  ! electron
   real(8), dimension(:), intent(inout) :: E_Distr_p_X, E_Distr_p_Y, E_Distr_p_Z  ! positron
   real(8), dimension(:,:), intent(inout) :: E_Distr_h_X, E_Distr_h_Y, E_Distr_h_Z  ! hole
   real(8), dimension(:), intent(inout) :: E_Distr_SHI_X, E_Distr_SHI_Y, E_Distr_SHI_Z    ! SHI
   real(8), dimension(:), intent(inout) :: E_Distr_a_X, E_Distr_a_Y, E_Distr_a_Z    ! Atom

   ! Along X:
   if (numpar%grid_par(1)%along_axis) then  ! collect data
      call sort_cartesian(MC%N_ph, MC%MC_Photons(:), numpar, tim, Distr_ph_X, E_Distr_ph_X, 1, neutral=.true.)  ! photons; below
      call sort_cartesian(MC%N_e, MC%MC_Electrons, numpar, tim,Distr_e_X, E_Distr_e_X, 1)  ! electrons; below
      call sort_cartesian(MC%N_p, MC%MC_Positrons, numpar, tim, Distr_p_X, E_Distr_p_X, 1)  ! positrons; below
      call sort_holes_cartesian(used_target, MC%N_h, MC%MC_Holes, numpar, tim, Distr_h_X, E_Distr_h_X, 1)  ! sort holes; below
      call sort_cartesian(MC%N_SHI, MC%MC_SHIs, numpar, tim, Distr_SHI_X, E_Distr_SHI_X, 1)  ! SHI; below
      call sort_cartesian(MC%N_at_nrg, MC%MC_Atoms_events, numpar, tim, Distr_a_X, E_Distr_a_X, 1)  ! Atoms; below
   endif
   
   ! Along Y:
   if (numpar%grid_par(2)%along_axis) then  ! collect data
      call sort_cartesian(MC%N_ph, MC%MC_Photons(:), numpar, tim, Distr_ph_Y, E_Distr_ph_Y, 2, neutral=.true.)  ! photons; below
      call sort_cartesian(MC%N_e, MC%MC_Electrons, numpar, tim,Distr_e_Y, E_Distr_e_Y, 2)  ! electrons; below
      call sort_cartesian(MC%N_p, MC%MC_Positrons, numpar, tim, Distr_p_Y, E_Distr_p_Y, 2)  ! positrons; below
      call sort_holes_cartesian(used_target, MC%N_h, MC%MC_Holes, numpar, tim, Distr_h_Y, E_Distr_h_Y, 2)  ! sort holes; below
      call sort_cartesian(MC%N_SHI, MC%MC_SHIs, numpar, tim, Distr_SHI_Y, E_Distr_SHI_Y, 2)  ! SHI; below
      call sort_cartesian(MC%N_at_nrg, MC%MC_Atoms_events, numpar, tim, Distr_a_Y, E_Distr_a_Y, 2)  ! Atoms; below
   endif
   
   ! Along Z:
   if (numpar%grid_par(3)%along_axis) then  ! collect data
      call sort_cartesian(MC%N_ph, MC%MC_Photons(:), numpar, tim, Distr_ph_Z, E_Distr_ph_Z, 3, neutral=.true.)  ! photons; below
      call sort_cartesian(MC%N_e, MC%MC_Electrons, numpar, tim,Distr_e_Z, E_Distr_e_Z, 3)  ! electrons; below
      call sort_cartesian(MC%N_p, MC%MC_Positrons, numpar, tim, Distr_p_Z, E_Distr_p_Z, 3)  ! positrons; below
      call sort_holes_cartesian(used_target, MC%N_h, MC%MC_Holes, numpar, tim, Distr_h_Z, E_Distr_h_Z, 3)  ! sort holes; below
      call sort_cartesian(MC%N_SHI, MC%MC_SHIs, numpar, tim, Distr_SHI_Z, E_Distr_SHI_Z, 3)  ! SHI; below
      call sort_cartesian(MC%N_at_nrg, MC%MC_Atoms_events, numpar, tim, Distr_a_Z, E_Distr_a_Z, 3)  ! Atoms; below
   endif
end subroutine get_distribution_1d_Cartesian


subroutine get_distribution_1d_Cylindrical(used_target, MC, numpar, tim, Distr_ph_R, Distr_ph_L, Distr_ph_Theta, Distr_e_R, Distr_e_L, Distr_e_Theta, &
                 Distr_p_R, Distr_p_L, Distr_p_Theta, Distr_h_R, Distr_h_L, Distr_h_Theta, Distr_SHI_R, Distr_SHI_L, Distr_SHI_Theta, &
                 Distr_a_R, Distr_a_L, Distr_a_Theta, &
                 E_Distr_ph_R, E_Distr_ph_L, E_Distr_ph_Theta, E_Distr_e_R, E_Distr_e_L, E_Distr_e_Theta, &
                 E_Distr_p_R, E_Distr_p_L, E_Distr_p_Theta, E_Distr_h_R, E_Distr_h_L, E_Distr_h_Theta, &
                 E_Distr_SHI_R, E_Distr_SHI_L, E_Distr_SHI_Theta, E_Distr_a_R, E_Distr_a_L, E_Distr_a_Theta)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   ! Density distributions:
   real(8), dimension(:), intent(inout) :: Distr_ph_R, Distr_ph_L, Distr_ph_Theta   ! photon
   real(8), dimension(:), intent(inout) :: Distr_e_R, Distr_e_L, Distr_e_Theta  ! electron
   real(8), dimension(:), intent(inout) :: Distr_p_R, Distr_p_L, Distr_p_Theta  ! positron
   real(8), dimension(:,:), intent(inout) :: Distr_h_R, Distr_h_L, Distr_h_Theta  ! hole
   real(8), dimension(:), intent(inout) :: Distr_SHI_R, Distr_SHI_L, Distr_SHI_Theta    ! SHI
   real(8), dimension(:), intent(inout) :: Distr_a_R, Distr_a_L, Distr_a_Theta    ! Atom
   ! Energy density distributions:
   real(8), dimension(:), intent(inout) :: E_Distr_ph_R, E_Distr_ph_L, E_Distr_ph_Theta   ! photon
   real(8), dimension(:), intent(inout) :: E_Distr_e_R, E_Distr_e_L, E_Distr_e_Theta  ! electron
   real(8), dimension(:), intent(inout) :: E_Distr_p_R, E_Distr_p_L, E_Distr_p_Theta  ! positron
   real(8), dimension(:,:), intent(inout) :: E_Distr_h_R, E_Distr_h_L, E_Distr_h_Theta  ! hole
   real(8), dimension(:), intent(inout) :: E_Distr_SHI_R, E_Distr_SHI_L, E_Distr_SHI_Theta    ! SHI
   real(8), dimension(:), intent(inout) :: E_Distr_a_R, E_Distr_a_L, E_Distr_a_Theta    ! Atom
   
   ! Along radius R:
   if (numpar%grid_par(8)%along_axis) then  ! collect kinetic energies:
      call sort_radially(MC%N_ph, MC%MC_Photons(:), numpar, tim, Distr_ph_R, E_Distr_ph_R, neutral=.true.)  ! photons; below
      call sort_radially(MC%N_e, MC%MC_Electrons, numpar, tim,Distr_e_R, E_Distr_e_R)  ! electrons; below
      call sort_radially(MC%N_p, MC%MC_Positrons, numpar, tim, Distr_p_R, E_Distr_p_R)  ! positrons; below
!       call sort_radially(MC%N_h, MC%MC_Holes, numpar, tim, Distr_h_R(1,:), E_Distr_h_R(1,:), valent=.true.)  ! valence holes; below
      call sort_holes_radially(used_target, MC%N_h, MC%MC_Holes, numpar, tim, Distr_h_R, E_Distr_h_R)  ! sort holes; below
      call sort_radially(MC%N_SHI, MC%MC_SHIs, numpar, tim, Distr_SHI_R, E_Distr_SHI_R)  ! SHI; below
      call sort_radially(MC%N_at_nrg, MC%MC_Atoms_events, numpar, tim, Distr_a_R, E_Distr_a_R)  ! Atoms; below
   endif

   ! Along depth L:
   if (numpar%grid_par(9)%along_axis) then  ! collect data
      ! NOT READY YET
   endif

   ! Along theta angle Theta:
   if (numpar%grid_par(10)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
end subroutine get_distribution_1d_Cylindrical


subroutine get_distribution_1d_Spherical(used_target, MC, numpar, tim, Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic, Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic, &
                 Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic,  Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic, Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic, &
                 Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic, &
                 E_Distr_ph_Rc, E_Distr_ph_Thetac, E_Distr_ph_Phic, E_Distr_e_Rc, E_Distr_e_Thetac, E_Distr_e_Phic, &
                 E_Distr_p_Rc, E_Distr_p_Thetac, E_Distr_p_Phic,  E_Distr_h_Rc, E_Distr_h_Thetac, E_Distr_h_Phic, &
                 E_Distr_SHI_Rc, E_Distr_SHI_Thetac, E_Distr_SHI_Phic, E_Distr_a_Rc, E_Distr_a_Thetac, E_Distr_a_Phic)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   ! Density distributions:
   real(8), dimension(:), intent(inout) :: Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic ! photon
   real(8), dimension(:), intent(inout) :: Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic ! electron
   real(8), dimension(:), intent(inout) :: Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic ! positron
   real(8), dimension(:,:), intent(inout) :: Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic ! hole
   real(8), dimension(:), intent(inout) :: Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic ! SHI
   real(8), dimension(:), intent(inout) :: Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic ! atom
   ! Eenrgy distributions:
   real(8), dimension(:), intent(inout) :: E_Distr_ph_Rc, E_Distr_ph_Thetac, E_Distr_ph_Phic ! photon
   real(8), dimension(:), intent(inout) :: E_Distr_e_Rc, E_Distr_e_Thetac, E_Distr_e_Phic ! electron
   real(8), dimension(:), intent(inout) :: E_Distr_p_Rc, E_Distr_p_Thetac, E_Distr_p_Phic ! positron
   real(8), dimension(:,:), intent(inout) :: E_Distr_h_Rc, E_Distr_h_Thetac, E_Distr_h_Phic ! hole
   real(8), dimension(:), intent(inout) :: E_Distr_SHI_Rc, E_Distr_SHI_Thetac, E_Distr_SHI_Phic ! SHI
   real(8), dimension(:), intent(inout) :: E_Distr_a_Rc, E_Distr_a_Thetac, E_Distr_a_Phic ! atom
   
   ! Along radius Rc:
   if (numpar%grid_par(14)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
   
   ! Along depth Thetac:
   if (numpar%grid_par(15)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
      
   ! Along theta angle Phic:
   if (numpar%grid_par(16)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
end subroutine get_distribution_1d_Spherical


                 
! Spatial distributions in 2d:
subroutine get_distribution_2d_Cartesian(used_target, MC, numpar, tim, Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ, Distr_e_XY, Distr_e_YZ, Distr_e_XZ, &
                 Distr_p_XY, Distr_p_YZ, Distr_p_XZ,  Distr_h_XY, Distr_h_YZ, Distr_h_XZ, Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ, &
                 Distr_a_XY, Distr_a_YZ, Distr_a_XZ, &
                 E_Distr_ph_XY, E_Distr_ph_YZ, E_Distr_ph_XZ, E_Distr_e_XY, E_Distr_e_YZ, E_Distr_e_XZ, &
                 E_Distr_p_XY, E_Distr_p_YZ, E_Distr_p_XZ, E_Distr_h_XY, E_Distr_h_YZ, E_Distr_h_XZ, &
                 E_Distr_SHI_XY, E_Distr_SHI_YZ, E_Distr_SHI_XZ, E_Distr_a_XY, E_Distr_a_YZ, E_Distr_a_XZ)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   ! Density distributions:
   real(8), dimension(:,:), intent(inout) :: Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ  ! photon
   real(8), dimension(:,:), intent(inout) :: Distr_e_XY, Distr_e_YZ, Distr_e_XZ    ! electron
   real(8), dimension(:,:), intent(inout) :: Distr_p_XY, Distr_p_YZ, Distr_p_XZ ! positron
   real(8), dimension(:,:,:), intent(inout) :: Distr_h_XY, Distr_h_YZ, Distr_h_XZ    ! hole
   real(8), dimension(:,:), intent(inout) :: Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ   ! SHI
   real(8), dimension(:,:), intent(inout) :: Distr_a_XY, Distr_a_YZ, Distr_a_XZ   ! atom
   ! Energy density distributions:
   real(8), dimension(:,:), intent(inout) :: E_Distr_ph_XY, E_Distr_ph_YZ, E_Distr_ph_XZ  ! photon
   real(8), dimension(:,:), intent(inout) :: E_Distr_e_XY, E_Distr_e_YZ, E_Distr_e_XZ    ! electron
   real(8), dimension(:,:), intent(inout) :: E_Distr_p_XY, E_Distr_p_YZ, E_Distr_p_XZ ! positron
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_h_XY, E_Distr_h_YZ, E_Distr_h_XZ    ! hole
   real(8), dimension(:,:), intent(inout) :: E_Distr_SHI_XY, E_Distr_SHI_YZ, E_Distr_SHI_XZ   ! SHI
   real(8), dimension(:,:), intent(inout) :: E_Distr_a_XY, E_Distr_a_YZ, E_Distr_a_XZ   ! atom

   ! Along radius XY:
   if (numpar%grid_par(4)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
   
   ! Along depth XZ:
   if (numpar%grid_par(5)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
      
   ! Along theta angle YZ:
   if (numpar%grid_par(6)%along_axis) then  ! collect data
      ! NOT READY YET
   endif

end subroutine get_distribution_2d_Cartesian



subroutine get_distribution_2d_Cylindrical(used_target, MC, numpar, tim, Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta,  Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta, &
                 Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta, Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta, Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta, &
                 Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta, &
                 E_Distr_ph_RL, E_Distr_ph_RTheta, E_Distr_ph_LTheta, E_Distr_e_RL, E_Distr_e_RTheta, E_Distr_e_LTheta, &
                 E_Distr_p_RL, E_Distr_p_RTheta, E_Distr_p_LTheta, E_Distr_h_RL, E_Distr_h_RTheta, E_Distr_h_LTheta, &
                 E_Distr_SHI_RL, E_Distr_SHI_RTheta, E_Distr_SHI_LTheta, E_Distr_a_RL, E_Distr_a_RTheta, E_Distr_a_LTheta)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   ! Density distributions:
   real(8), dimension(:,:), intent(inout) :: Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta   ! photon
   real(8), dimension(:,:), intent(inout) :: Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta ! electron
   real(8), dimension(:,:), intent(inout) :: Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta ! positron
   real(8), dimension(:,:,:), intent(inout) :: Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta ! hole
   real(8), dimension(:,:), intent(inout) :: Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta ! SHI
   real(8), dimension(:,:), intent(inout) :: Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta ! atom
   ! Energy density distributions:
   real(8), dimension(:,:), intent(inout) :: E_Distr_ph_RL, E_Distr_ph_RTheta, E_Distr_ph_LTheta   ! photon
   real(8), dimension(:,:), intent(inout) :: E_Distr_e_RL, E_Distr_e_RTheta, E_Distr_e_LTheta ! electron
   real(8), dimension(:,:), intent(inout) :: E_Distr_p_RL, E_Distr_p_RTheta, E_Distr_p_LTheta ! positron
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_h_RL, E_Distr_h_RTheta, E_Distr_h_LTheta ! hole
   real(8), dimension(:,:), intent(inout) :: E_Distr_SHI_RL, E_Distr_SHI_RTheta, E_Distr_SHI_LTheta ! SHI
   real(8), dimension(:,:), intent(inout) :: E_Distr_a_RL, E_Distr_a_RTheta, E_Distr_a_LTheta ! atom
    
   ! Along radius R L:
   if (numpar%grid_par(11)%along_axis) then  ! collect data
      call sort_RL(MC%N_ph, MC%MC_Photons(:), numpar, tim, Distr_ph_RL, E_Distr_ph_RL, neutral=.true.)  ! photons; below
      call sort_RL(MC%N_e, MC%MC_Electrons, numpar, tim,Distr_e_RL, E_Distr_e_RL)  ! electrons; below
      call sort_RL(MC%N_p, MC%MC_Positrons, numpar, tim, Distr_p_RL, E_Distr_p_RL)  ! positrons; below
      call sort_holes_RL(used_target, MC%N_h, MC%MC_Holes, numpar, tim, Distr_h_RL, E_Distr_h_RL)  ! sort holes; below
      call sort_RL(MC%N_SHI, MC%MC_SHIs, numpar, tim, Distr_SHI_RL, E_Distr_SHI_RL)  ! SHI; below
      call sort_RL(MC%N_at_nrg, MC%MC_Atoms_events, numpar, tim, Distr_a_RL, E_Distr_a_RL)  ! Atoms; below
   endif
   
   ! Along depth R Theta:
   if (numpar%grid_par(12)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
      
   ! Along theta angle L Theta -- not included at all   
end subroutine get_distribution_2d_Cylindrical


subroutine get_distribution_2d_Spherical(used_target, MC, numpar, tim, Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic, &
                 Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic, Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic, &
                 Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic,  Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic, &
                 Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic, &
                 E_Distr_ph_RcThc, E_Distr_ph_RcPhic, E_Distr_ph_ThcPhic, &
                 E_Distr_e_RcThc, E_Distr_e_RcPhic, E_Distr_e_ThcPhic, E_Distr_p_RcThc, E_Distr_p_RcPhic, E_Distr_p_ThcPhic, &
                 E_Distr_h_RcThc, E_Distr_h_RcPhic, E_Distr_h_ThcPhic, E_Distr_SHI_RcThc, E_Distr_SHI_RcPhic, E_Distr_SHI_ThcPhic, &
                 E_Distr_a_RcThc, E_Distr_a_RcPhic, E_Distr_a_ThcPhic )
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   ! Density distributions:
   real(8), dimension(:,:), intent(inout) :: Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic  ! photon
   real(8), dimension(:,:), intent(inout) :: Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic  ! electron
   real(8), dimension(:,:), intent(inout) :: Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic  ! positron
   real(8), dimension(:,:,:), intent(inout) :: Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic  ! hole
   real(8), dimension(:,:), intent(inout) :: Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic  ! SHI
   real(8), dimension(:,:), intent(inout) :: Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic  ! atom
   ! Energy density distributions:
   real(8), dimension(:,:), intent(inout) :: E_Distr_ph_RcThc, E_Distr_ph_RcPhic, E_Distr_ph_ThcPhic  ! photon
   real(8), dimension(:,:), intent(inout) :: E_Distr_e_RcThc, E_Distr_e_RcPhic, E_Distr_e_ThcPhic  ! electron
   real(8), dimension(:,:), intent(inout) :: E_Distr_p_RcThc, E_Distr_p_RcPhic, E_Distr_p_ThcPhic  ! positron
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_h_RcThc, E_Distr_h_RcPhic, E_Distr_h_ThcPhic  ! hole
   real(8), dimension(:,:), intent(inout) :: E_Distr_SHI_RcThc, E_Distr_SHI_RcPhic, E_Distr_SHI_ThcPhic  ! SHI
   real(8), dimension(:,:), intent(inout) :: E_Distr_a_RcThc, E_Distr_a_RcPhic, E_Distr_a_ThcPhic  ! atom
   
   ! Along radius R Phi:
   if (numpar%grid_par(17)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
   
   ! Along depth R Theta:
   if (numpar%grid_par(18)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
      
   ! Along theta angle Phi Theta -- not included at all   
end subroutine get_distribution_2d_Spherical
                 

! Spatial distributions in 3d:
subroutine get_distribution_3d_Cartesian(used_target, MC, numpar, tim, Distr_ph_XYZ, Distr_e_XYZ, Distr_p_XYZ, Distr_h_XYZ, Distr_SHI_XYZ, Distr_a_XYZ, &
                                                              E_Distr_ph_XYZ, E_Distr_e_XYZ, E_Distr_p_XYZ, E_Distr_h_XYZ, E_Distr_SHI_XYZ, E_Distr_a_XYZ)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   ! Density distribution:
   real(8), dimension(:,:,:), intent(inout) :: Distr_ph_XYZ ! photon
   real(8), dimension(:,:,:), intent(inout) :: Distr_e_XYZ  ! electron
   real(8), dimension(:,:,:), intent(inout) :: Distr_p_XYZ  ! positron
   real(8), dimension(:,:,:,:), intent(inout) :: Distr_h_XYZ  ! hole
   real(8), dimension(:,:,:), intent(inout) :: Distr_SHI_XYZ    ! SHI
   real(8), dimension(:,:,:), intent(inout) :: Distr_a_XYZ    ! atom
   ! Energy density distribution:
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_ph_XYZ ! photon
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_e_XYZ  ! electron
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_p_XYZ  ! positron
   real(8), dimension(:,:,:,:), intent(inout) :: E_Distr_h_XYZ  ! hole
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_SHI_XYZ    ! SHI
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_a_XYZ    ! atom
   ! Along radius XYZ:
   if (numpar%grid_par(7)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
end subroutine get_distribution_3d_Cartesian


subroutine get_distribution_3d_Cylindrical(used_target, MC, numpar, tim, Distr_ph_RLTheta, Distr_e_RLTheta, Distr_p_RLTheta, Distr_h_RLTheta, Distr_SHI_RLTheta, &
                    Distr_a_RLTheta, E_Distr_ph_RLTheta, E_Distr_e_RLTheta, E_Distr_p_RLTheta, E_Distr_h_RLTheta, E_Distr_SHI_RLTheta, E_Distr_a_RLTheta)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   ! Density distributions:
   real(8), dimension(:,:,:), intent(inout) :: Distr_ph_RLTheta ! photon
   real(8), dimension(:,:,:), intent(inout) :: Distr_e_RLTheta  ! electron
   real(8), dimension(:,:,:), intent(inout) :: Distr_p_RLTheta  ! positron
   real(8), dimension(:,:,:,:), intent(inout) :: Distr_h_RLTheta  ! hole
   real(8), dimension(:,:,:), intent(inout) :: Distr_SHI_RLTheta    ! SHI
   real(8), dimension(:,:,:), intent(inout) :: Distr_a_RLTheta    ! atom
   ! Energy density distributions:
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_ph_RLTheta ! photon
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_e_RLTheta  ! electron
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_p_RLTheta  ! positron
   real(8), dimension(:,:,:,:), intent(inout) :: E_Distr_h_RLTheta  ! hole
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_SHI_RLTheta    ! SHI
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_a_RLTheta    ! atom
   ! Along radius RLTheta:
   if (numpar%grid_par(13)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
end subroutine get_distribution_3d_Cylindrical


subroutine get_distribution_3d_Spherical(used_target, MC, numpar, tim, Distr_ph_RcThcPhic, Distr_e_RcThcPhic, Distr_p_RcThcPhic, Distr_h_RcThcPhic, &
                        Distr_SHI_RcThcPhic, Distr_a_RcThcPhic, &
                        E_Distr_ph_RcThcPhic, E_Distr_e_RcThcPhic, E_Distr_p_RcThcPhic, E_Distr_h_RcThcPhic, &
                        E_Distr_SHI_RcThcPhic, E_Distr_a_RcThcPhic)
   type(Matter), intent(in) :: used_target   ! parameters of the target
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   ! Density distributions:
   real(8), dimension(:,:,:), intent(inout) :: Distr_ph_RcThcPhic ! photon
   real(8), dimension(:,:,:), intent(inout) :: Distr_e_RcThcPhic ! electron
   real(8), dimension(:,:,:), intent(inout) :: Distr_p_RcThcPhic ! positron
   real(8), dimension(:,:,:,:), intent(inout) :: Distr_h_RcThcPhic ! hole
   real(8), dimension(:,:,:), intent(inout) :: Distr_SHI_RcThcPhic ! SHI
   real(8), dimension(:,:,:), intent(inout) :: Distr_a_RcThcPhic ! atom
   ! Energy density distributions:
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_ph_RcThcPhic ! photon
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_e_RcThcPhic ! electron
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_p_RcThcPhic ! positron
   real(8), dimension(:,:,:,:), intent(inout) :: E_Distr_h_RcThcPhic ! hole
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_SHI_RcThcPhic ! SHI
   real(8), dimension(:,:,:), intent(inout) :: E_Distr_a_RcThcPhic ! atom
   ! Along radius RThetaPhi:
   if (numpar%grid_par(19)%along_axis) then  ! collect data
      ! NOT READY YET
   endif
end subroutine get_distribution_3d_Spherical



!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

subroutine sort_holes_cartesian(used_target, N_prtcl, MC_Prtcl, numpar, tim, Distr_h_X, E_Distr_h_X, coord_ind)
  type(Matter), intent(in), target :: used_target   ! parameters of the target
   integer, intent(in) :: N_prtcl   ! number of particles of a given kind
   type(Hole), dimension(:), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:,:), intent(inout) :: Distr_h_X, E_Distr_h_X     ! density and energy density distributions for holes in all shells
   integer, intent(in) :: coord_ind ! along which axis: 1=X, 2=Y, 3=Z
   !------------------------------
   real(8) :: S, R_1, R_2
   integer :: i, N_arr
   logical :: anything_to_do

   ! Check if there is any active particle:
   anything_to_do = any(MC_Prtcl(:)%active)

   ACTPAR:if (anything_to_do) then ! there are particles to distribute
      ! To find volume, we define the area according to the simulation box size:
      select case (coord_ind)
      case (1)  ! along X
         R_1 = abs(numpar%box_end_z - numpar%box_start_z)   ! [A]
         R_2 = abs(numpar%box_end_y - numpar%box_start_y)   ! [A]
      case (2)  ! along Y
         R_1 = abs(numpar%box_end_z - numpar%box_start_z)   ! [A]
         R_2 = abs(numpar%box_end_x - numpar%box_start_x)   ! [A]
      case (3)  ! along Z
         R_1 = abs(numpar%box_end_x - numpar%box_start_x)   ! [A]
         R_2 = abs(numpar%box_end_y - numpar%box_start_y)   ! [A]
      end select
      S = R_1 * R_2 ! area [A^2]
      
      ! Go through all electrons and distribute them into arrays:
      do i = 1, N_prtcl
         ! Include only active particles:
         if (MC_Prtcl(i)%active) then
            ! Find to which array to add this hole, according to its element and shell:
            N_arr = hole_number_in_array(used_target, MC_Prtcl(i)%in_target, MC_Prtcl(i)%KOA, MC_Prtcl(i)%Sh)     ! below
            ! Add into the corresponding arrays:
            call add_cartesian_hole(used_target, MC_Prtcl(i), numpar, tim, Distr_h_X(N_arr,:), E_Distr_h_X(N_arr,:), S, coord_ind) ! below
         endif ! (MC_Prtcl(i)%active)
      enddo ! i = 1, N_prtcl
   endif ACTPAR
end subroutine sort_holes_cartesian



subroutine add_cartesian_hole(used_target, MC_Prtcl, numpar, tim, Distr_R, E_Distr_R, S, coord_ind)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Hole), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:), intent(inout) :: Distr_R, E_Distr_R     ! density and energy density distributions
   real(8), intent(in) :: S ! [A^2] area
   integer, intent(in) :: coord_ind ! along which axis: 1=X, 2=Y, 3=Z
   !-----------------------------
   integer :: i_arr, i_ax, N_siz
   real(8) :: R, dR, dV, Rcur(3)
   type(Target_atoms), pointer :: matter
   
   matter => used_target%Material(MC_Prtcl%in_target)
   
   ! Find where the particle is at the time instant "tim":
   if (MC_Prtcl%valent) then    ! only valent holes can move, not core holes
      !Rcur(:) = MC_Prtcl%R0(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
      Rcur(:) = MC_Prtcl%R(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
   else ! core holes are immobile
!       Rcur(:) = MC_Prtcl%R0(:)  ! [A] 
      Rcur(:) = MC_Prtcl%R(:)  ! [A] 
   endif
   
  ! Find where it is on the grid:
   select case (coord_ind)
   case (1)  ! along X
!       call Find_in_array_monoton(numpar%grids(1)%spatial_grid1(:), Rcur(1), i_arr)  ! module "Little_subroutines"
!       ! Define the thickness of the cartesian layer of the grid:
!       if (i_arr == size(numpar%grids(1)%spatial_grid1)) then   ! the last point on the grid
!          dR = numpar%grids(1)%spatial_grid1(i_arr) - numpar%grids(1)%spatial_grid1(i_arr-1)
!       else
!          dR = numpar%grids(1)%spatial_grid1(i_arr+1) - numpar%grids(1)%spatial_grid1(i_arr)
!       endif
      i_ax = 1
      N_siz = size(numpar%grids(i_ax)%spatial_grid1)
      if (Rcur(1) < numpar%grids(i_ax)%spatial_grid1(1)) then   ! R below lower limit
         i_arr = 1
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr)
      else if (Rcur(1) >= numpar%grids(i_ax)%spatial_grid1(N_siz)) then  ! above the max R grid point
         i_arr = N_siz
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      else ! inside the grid
         call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), Rcur(1), i_arr)  ! module "Little_subroutines"
         i_arr = i_arr + 1 ! assign particle to the end of the interval
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      endif
   case (2)  ! along Y
!       call Find_in_array_monoton(numpar%grids(2)%spatial_grid1(:), Rcur(2), i_arr)  ! module "Little_subroutines"
!       ! Define the thickness of the cartesian layer of the grid:
!       if (i_arr == size(numpar%grids(2)%spatial_grid1)) then   ! the last point on the grid
!          dR = numpar%grids(2)%spatial_grid1(i_arr) - numpar%grids(2)%spatial_grid1(i_arr-1)
!       else
!          dR = numpar%grids(2)%spatial_grid1(i_arr+1) - numpar%grids(2)%spatial_grid1(i_arr)
!       endif
      i_ax = 2
      N_siz = size(numpar%grids(i_ax)%spatial_grid1)
      if (Rcur(2) < numpar%grids(i_ax)%spatial_grid1(1)) then   ! R below lower limit
         i_arr = 1
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr)
      else if (Rcur(2) >= numpar%grids(i_ax)%spatial_grid1(N_siz)) then  ! above the max R grid point
         i_arr = N_siz
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      else ! inside the grid
         call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), Rcur(2), i_arr)  ! module "Little_subroutines"
         i_arr = i_arr + 1 ! assign particle to the end of the interval
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      endif
   case (3)  ! along Z
!       call Find_in_array_monoton(numpar%grids(3)%spatial_grid1(:), Rcur(3), i_arr)  ! module "Little_subroutines"
!       ! Define the thickness of the cartesian layer of the grid:
!       if (i_arr == size(numpar%grids(3)%spatial_grid1)) then   ! the last point on the grid
!          dR = numpar%grids(3)%spatial_grid1(i_arr) - numpar%grids(3)%spatial_grid1(i_arr-1)
!       else
!          dR = numpar%grids(3)%spatial_grid1(i_arr+1) - numpar%grids(3)%spatial_grid1(i_arr)
!       endif
      i_ax = 3
      N_siz = size(numpar%grids(i_ax)%spatial_grid1)
      if (Rcur(3) < numpar%grids(i_ax)%spatial_grid1(1)) then   ! R below lower limit
         i_arr = 1
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr)
      else if (Rcur(3) > numpar%grids(i_ax)%spatial_grid1(N_siz)) then  ! above the max R grid point
         i_arr = N_siz
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      else ! inside the grid
         call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), Rcur(3), i_arr)  ! module "Little_subroutines"
         i_arr = i_arr + 1 ! assign particle to the end of the interval
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      endif

   end select
!    print*, numpar%grids(3)%spatial_grid1(i_arr), numpar%grids(3)%spatial_grid1(i_arr+1), R

   ! Define volume:
   dV = dR*S  ! [A^3]
   ! Add particle into the distribution, normalized per volume to get density:
   Distr_R(i_arr) = Distr_R(i_arr) + 1.0d0/dV    ! add a particle into this array [1/A^3]
   ! And corresponding potential energy density:
   if (MC_Prtcl%valent) then    ! for valence holes, use band gap as potential energy
      E_Distr_R(i_arr) = E_Distr_R(i_arr) + (MC_Prtcl%Ekin + matter%DOS%Egap)/dV    ! add a kinetic and potential energies into this array [eV/A^3]
   else ! core shells
      E_Distr_R(i_arr) = E_Distr_R(i_arr) + matter%Elements(MC_Prtcl%KOA)%Ip(MC_Prtcl%Sh)/dV    ! add a ionization potential [eV/A^3]
   endif
   
   nullify(matter)
end subroutine add_cartesian_hole




subroutine sort_holes_radially(used_target, N_prtcl, MC_Prtcl, numpar, tim, Distr_h_R, E_Distr_h_R)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   integer, intent(in) :: N_prtcl   ! number of particles of a given kind
   type(Hole), dimension(:), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:,:), intent(inout) :: Distr_h_R, E_Distr_h_R     ! density and energy density distributions for holes in all shells
   !------------------------------
   real(8) :: L
   integer :: i, N_arr
   logical :: anything_to_do

   ! Check if there is any active particle:
   anything_to_do = any(MC_Prtcl(:)%active)

   ACTPAR:if (anything_to_do) then ! there are particles to distribute
      ! To find volume, we define the depth according to the simulation box size:
      L =  abs(numpar%box_end_z - numpar%box_start_z)   ! [A]
      ! Go through all electrons and distribute them into arrays:
      
!       print*, 'sort_holes_radially', N_prtcl
      
      do i = 1, N_prtcl
         ! Include only active particles:
         if (MC_Prtcl(i)%active) then
            ! Find to which array to add this hole, according to its element and shell:
            N_arr = hole_number_in_array(used_target, MC_Prtcl(i)%in_target, MC_Prtcl(i)%KOA, MC_Prtcl(i)%Sh)     ! below
!             print*, 'sort_holes_radially', N_arr, MC_Prtcl(i)%KOA, MC_Prtcl(i)%Sh
            ! Add into the corresponding arrays:
            call add_radial_hole(used_target, MC_Prtcl(i), numpar, tim, Distr_h_R(N_arr,:), E_Distr_h_R(N_arr,:), L) ! below
         endif ! (MC_Prtcl(i)%active)
      enddo ! i = 1, N_prtcl
   endif ACTPAR
end subroutine sort_holes_radially

    
subroutine sort_holes_RL(used_target, N_prtcl, MC_Prtcl, numpar, tim, Distr_h_RL, E_Distr_h_RL)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   integer, intent(in) :: N_prtcl   ! number of particles of a given kind
   type(Hole), dimension(:), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:,:,:), intent(inout) :: Distr_h_RL, E_Distr_h_RL     ! density and energy density distributions for holes in all shells
   !------------------------------
   real(8) :: L
   integer :: i, N_arr
   logical :: anything_to_do

   ! Check if there is any active particle:
   anything_to_do = any(MC_Prtcl(:)%active)

   ACTPAR:if (anything_to_do) then ! there are particles to distribute
      ! Go through all electrons and distribute them into arrays:
      do i = 1, N_prtcl
         ! Include only active particles:
         if (MC_Prtcl(i)%active) then
            ! Find to which array to add this hole, according to its element and shell:
            N_arr = hole_number_in_array(used_target, MC_Prtcl(i)%in_target, MC_Prtcl(i)%KOA, MC_Prtcl(i)%Sh)     ! below
            ! Add into the corresponding arrays:
            call add_RL_hole(used_target, MC_Prtcl(i), numpar, tim, Distr_h_RL(N_arr,:,:), E_Distr_h_RL(N_arr,:,:)) ! below
         endif ! (MC_Prtcl(i)%active)
      enddo ! i = 1, N_prtcl
   endif ACTPAR
end subroutine sort_holes_RL

pure function hole_number_in_array(used_target, Nt, KOA, SHL) result(i_arr)
   integer i_arr    ! number of the array to be used for sorting this holes data
   type(Matter), intent(in) :: used_target   ! parameters of the target
   integer, intent(in) :: Nt, KOA, SHL  ! number of target, kind of atom, and shell, the given hole is in
   integer :: i, j, k
   i_arr = 1    ! valence band to start counting from
   if ((KOA > 0) .and. (SHL > 0)) then  ! it is indeed a core shell
      do i = 1, Nt  ! how many targets
         if (i < Nt) then   ! not in this target
            do j = 1, size(used_target%Material(i)%Elements(:)) ! add number of all shells in this element
               i_arr = i_arr + used_target%Material(i)%Elements(j)%N_core_shl
            enddo
         else   ! in this target
            do j = 1, KOA
               if (j < KOA) then    ! not in this element
                  i_arr = i_arr + used_target%Material(i)%Elements(j)%N_core_shl    ! add all shells of this element
               else ! in this element
                  i_arr = i_arr + SHL
               endif
            enddo ! j = 1, KOA
         endif ! (i < Nt) 
      enddo ! i = 1, Nt  
   endif  !  ((KOA > 0) .and. (SHL > 0)) 
end function hole_number_in_array



subroutine sort_cartesian(N_prtcl, MC_Prtcl, numpar, tim, Distr_X, E_Distr_X, coord_ind, valent, neutral)
   integer, intent(in) :: N_prtcl   ! number of particles of a given kind
   class(Particle), dimension(:), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:), intent(inout) :: Distr_X, E_Distr_X     ! density and energy density distributions
   integer, intent(in) :: coord_ind ! 1=X, 2=Y, 3=Z, along which coordinate
   logical, intent(in), optional :: valent
   logical, intent(in), optional :: neutral ! is it a neutral particle (photon)?
   !------------------------------
   real(8) :: S, R_1, R_2
   integer :: i
   logical :: anything_to_do

   ! Check if there is any active particle:
   anything_to_do = any(MC_Prtcl(:)%active)

   ACTPAR:if (anything_to_do) then ! there are particles to distribute
      ! To find volume, we define the area according to the simulation box size:
      select case (coord_ind)
      case (1)  ! along X
         R_1 = abs(numpar%box_end_z - numpar%box_start_z)   ! [A]
         R_2 = abs(numpar%box_end_y - numpar%box_start_y)   ! [A]
      case (2)  ! along Y
         R_1 = abs(numpar%box_end_z - numpar%box_start_z)   ! [A]
         R_2 = abs(numpar%box_end_x - numpar%box_start_x)   ! [A]
      case (3)  ! along Z
         R_1 = abs(numpar%box_end_x - numpar%box_start_x)   ! [A]
         R_2 = abs(numpar%box_end_y - numpar%box_start_y)   ! [A]
      end select
      S = R_1 * R_2 ! area [A^2]
      
      ! Go through all electrons and distribute them into arrays:
      do i = 1, N_prtcl
         ! Include only active particles:
         if (MC_Prtcl(i)%active) then
            if (present(valent)) then
               select type(MC_Prtcl)
               type is (Hole)    ! it must be hole, only then attribute "valent" exists
                  if (MC_Prtcl(i)%valent) then
                     call add_cartesian_particle(MC_Prtcl(i), numpar, tim, Distr_X, E_Distr_X, S, coord_ind) ! below
                  endif
               endselect
            else
               if (present(neutral)) then
                  if (neutral) then ! photon etc.
                     call add_cartesian_particle(MC_Prtcl(i), numpar, tim, Distr_X, E_Distr_X, S, coord_ind, neutral) ! below
                  else  ! charged particle
                     call add_cartesian_particle(MC_Prtcl(i), numpar, tim, Distr_X, E_Distr_X, S, coord_ind) ! below
                  endif
               else
                  call add_cartesian_particle(MC_Prtcl(i), numpar, tim, Distr_X, E_Distr_X, S, coord_ind) ! below
               endif
            endif ! (present(valent))
         endif ! (MC_Prtcl(i)%active) 
      enddo ! i = 1, N_prtcl
   endif ACTPAR
end subroutine sort_cartesian



subroutine add_cartesian_particle(MC_Prtcl, numpar, tim, Distr_R, E_Distr_R, S, coord_ind, neutral)
   class(Particle), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:), intent(inout) :: Distr_R, E_Distr_R     ! density and energy density distributions
   real(8), intent(in) :: S ! [A^2] area
   integer, intent(in) :: coord_ind ! 1=X, 2=Y, 3=Z, along which coordinate
   logical, intent(in), optional :: neutral ! is it a neutral particle (photon)?
   !-----------------------------
   integer :: i_arr, N_siz, i_ax
   real(8) :: dR, dV, Rcur(3)

   ! Find where the particle is at the time instant "tim":
   !Rcur(:) = MC_Prtcl%R0(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
   Rcur(:) = MC_Prtcl%R(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
   if (present(neutral)) then   ! make sure periodic boundary conditions are observed:
      if (neutral) call put_back_into_box(numpar, R = Rcur) ! module "MC_general_tools"
   endif
   
!    if (Rcur(3) >= 10.0d0) then
!       print*, 'An electron crossed the simulation box border:'
!       print*, 'add_cartesian_particle', Rcur(3), MC_Prtcl%active, MC_Prtcl%R(3)
!       print*, 'add_cartesian_parti 2:', MC_Prtcl%V(3), tim, MC_Prtcl%t0
!    endif
   
   ! Find where it is on the grid:
   select case (coord_ind)
   case (1)  ! along X
!       call Find_in_array_monoton(numpar%grids(1)%spatial_grid1(:), Rcur(1), i_arr)  ! module "Little_subroutines"
!       ! Define the thickness of the cartesian layer of the grid:
!       if (i_arr == size(numpar%grids(1)%spatial_grid1)) then   ! the last point on the grid
!          dR = numpar%grids(1)%spatial_grid1(i_arr) - numpar%grids(1)%spatial_grid1(i_arr-1)
!       else
!          dR = numpar%grids(1)%spatial_grid1(i_arr+1) - numpar%grids(1)%spatial_grid1(i_arr)
!       endif
      i_ax = 1
      N_siz = size(numpar%grids(i_ax)%spatial_grid1)
      if (Rcur(1) < numpar%grids(i_ax)%spatial_grid1(1)) then   ! R below lower limit
         i_arr = 1
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr)
      else if (Rcur(1) >= numpar%grids(i_ax)%spatial_grid1(N_siz)) then  ! above the max R grid point
         i_arr = N_siz
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      else ! inside the grid
         call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), Rcur(1), i_arr)  ! module "Little_subroutines"
         i_arr = i_arr + 1 ! assign particle to the end of the interval
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      endif

   case (2)  ! along Y
!       call Find_in_array_monoton(numpar%grids(2)%spatial_grid1(:), Rcur(2), i_arr)  ! module "Little_subroutines"
!       ! Define the thickness of the cartesian layer of the grid:
!       if (i_arr == size(numpar%grids(2)%spatial_grid1)) then   ! the last point on the grid
!          dR = numpar%grids(2)%spatial_grid1(i_arr) - numpar%grids(2)%spatial_grid1(i_arr-1)
!       else
!          dR = numpar%grids(2)%spatial_grid1(i_arr+1) - numpar%grids(2)%spatial_grid1(i_arr)
!       endif
      i_ax = 2
      N_siz = size(numpar%grids(i_ax)%spatial_grid1)
      if (Rcur(2) < numpar%grids(i_ax)%spatial_grid1(1)) then   ! R below lower limit
         i_arr = 1
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr)
      else if (Rcur(2) >= numpar%grids(i_ax)%spatial_grid1(N_siz)) then  ! above the max R grid point
         i_arr = N_siz
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      else ! inside the grid
         call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), Rcur(2), i_arr)  ! module "Little_subroutines"
         i_arr = i_arr + 1 ! assign particle to the end of the interval
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      endif

   case (3)  ! along Z
   
      if (isnan(Rcur(3))) then 
         print*, 'Error in add_cartesian_particle #3', Rcur(3), MC_Prtcl%R(:), (MC_Prtcl%V(:)), (tim - MC_Prtcl%t0)
         print*, 'R0', MC_Prtcl%R0(:)
         print*, 'E', MC_Prtcl%Ekin
      endif
      
!       call Find_in_array_monoton(numpar%grids(3)%spatial_grid1(:), Rcur(3), i_arr)  ! module "Little_subroutines"
!       ! Define the thickness of the cartesian layer of the grid:
!       if (i_arr == size(numpar%grids(3)%spatial_grid1)) then   ! the last point on the grid
!          dR = numpar%grids(3)%spatial_grid1(i_arr) - numpar%grids(3)%spatial_grid1(i_arr-1)
!       else
!          dR = numpar%grids(3)%spatial_grid1(i_arr+1) - numpar%grids(3)%spatial_grid1(i_arr)
!       endif
      i_ax = 3
      N_siz = size(numpar%grids(i_ax)%spatial_grid1)
      if (Rcur(3) < numpar%grids(i_ax)%spatial_grid1(1)) then   ! R below lower limit
         i_arr = 1
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr)
      else if (Rcur(3) >= numpar%grids(i_ax)%spatial_grid1(N_siz)) then  ! above the max R grid point
         i_arr = N_siz
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      else ! inside the grid
         call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), Rcur(3), i_arr)  ! module "Little_subroutines"
         i_arr = i_arr + 1 ! assign particle to the end of the interval
         ! Define the thickness of the cylindrical layer of the grid:
         dR = numpar%grids(i_ax)%spatial_grid1(i_arr) - numpar%grids(i_ax)%spatial_grid1(i_arr-1)
      endif

   end select
!    print*, numpar%grids(3)%spatial_grid1(i_arr), numpar%grids(3)%spatial_grid1(i_arr+1), R
   
   ! Define volume:
   dV = S*dR
   
   ! Add particle into the distribution, normalized per volume to get density:
   Distr_R(i_arr) = Distr_R(i_arr) + 1.0d0/dV    ! add a particle into this array [1/A^3]
   
   ! And corresponding energy density:
   E_Distr_R(i_arr) = E_Distr_R(i_arr) + MC_Prtcl%Ekin/dV    ! add a particle energy into this array [eV/A^3]
end subroutine add_cartesian_particle




subroutine sort_radially(N_prtcl, MC_Prtcl, numpar, tim, Distr_R, E_Distr_R, valent, neutral)
   integer, intent(in) :: N_prtcl   ! number of particles of a given kind
   class(Particle), dimension(:), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:), intent(inout) :: Distr_R, E_Distr_R     ! density and energy density distributions
   logical, intent(in), optional :: valent  ! is it a valent hole?
   logical, intent(in), optional :: neutral ! is it a neutral particle (photon)?
   !------------------------------
   real(8) :: L
   integer :: i
   logical :: anything_to_do

   ! Check if there is any active particle:
   anything_to_do = any(MC_Prtcl(:)%active)

   ACTPAR:if (anything_to_do) then ! there are particles to distribute
      ! To find volume, we define the depth according to the simulation box size:
      L =  abs(numpar%box_end_z - numpar%box_start_z)   ! [A]
      ! Go through all electrons and distribute them into arrays:
      do i = 1, N_prtcl
         ! Include only active particles:
         if (MC_Prtcl(i)%active) then
            if (present(valent)) then
               select type(MC_Prtcl)
               type is (Hole)    ! it must be hole, only then attribute "valent" exists
                  if (MC_Prtcl(i)%valent) then
                     call add_radial_particle(MC_Prtcl(i), numpar, tim, Distr_R, E_Distr_R, L) ! below
                  endif
               endselect
            else
               if (present(neutral)) then
                  if (neutral) then ! photons
                     call add_radial_particle(MC_Prtcl(i), numpar, tim, Distr_R, E_Distr_R, L, neutral=.true.) ! below
                  else
                     call add_radial_particle(MC_Prtcl(i), numpar, tim, Distr_R, E_Distr_R, L) ! below
                  endif
               else
                  call add_radial_particle(MC_Prtcl(i), numpar, tim, Distr_R, E_Distr_R, L) ! below
               endif
            endif ! (present(valent))
         endif ! (MC_Prtcl(i)%active) 
      enddo ! i = 1, N_prtcl
   endif ACTPAR
end subroutine sort_radially

    
    
subroutine sort_RL(N_prtcl, MC_Prtcl, numpar, tim, Distr_RL, E_Distr_RL, valent, neutral)
   integer, intent(in) :: N_prtcl   ! number of particles of a given kind
   class(Particle), dimension(:), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:,:), intent(inout) :: Distr_RL, E_Distr_RL     ! density and energy density distributions
   logical, intent(in), optional :: valent  ! is it a valent hole?
   logical, intent(in), optional :: neutral ! is it a neutral particle (photon)?
   !------------------------------
   real(8) :: L
   integer :: i
   logical :: anything_to_do

   ! Check if there is any active particle:
   anything_to_do = any(MC_Prtcl(:)%active)

   ACTPAR:if (anything_to_do) then ! there are particles to distribute
      ! Go through all particles and distribute them into arrays:
      do i = 1, N_prtcl
         ! Include only active particles:
         if (MC_Prtcl(i)%active) then
            if (present(valent)) then
               select type(MC_Prtcl)
               type is (Hole)    ! it must be hole, only then attribute "valent" exists
                  if (MC_Prtcl(i)%valent) then
                     call add_RL_particle(MC_Prtcl(i), numpar, tim, Distr_RL, E_Distr_RL) ! below
                  endif
               endselect
            else
               if (present(neutral)) then   ! photon
                  if (neutral) then ! photon
                     call add_RL_particle(MC_Prtcl(i), numpar, tim, Distr_RL, E_Distr_RL, neutral=.true.) ! below
                  else
                     call add_RL_particle(MC_Prtcl(i), numpar, tim, Distr_RL, E_Distr_RL) ! below
                  endif
               else ! charged particle
                  call add_RL_particle(MC_Prtcl(i), numpar, tim, Distr_RL, E_Distr_RL) ! below
               endif
            endif ! (present(valent))
         endif ! (MC_Prtcl(i)%active) 
      enddo ! i = 1, N_prtcl
   endif ACTPAR
end subroutine sort_RL



subroutine add_radial_hole(used_target, MC_Prtcl, numpar, tim, Distr_R, E_Distr_R, L)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Hole), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:), intent(inout) :: Distr_R, E_Distr_R     ! density and energy density distributions
   real(8), intent(in) :: L ! [A] thickness
   !-----------------------------
   integer :: i_arr, i_ax, N_siz
   real(8) :: R, dR2, dV, Rcur(3)
   type(Target_atoms), pointer :: matter
   
   matter => used_target%Material(MC_Prtcl%in_target)
   i_ax = 8
   N_siz = size(numpar%grids(i_ax)%spatial_grid1)
   
   ! Find where the particle is at the time instant "tim":
   if (MC_Prtcl%valent) then    ! only valent holes can move, not core holes
!       Rcur(:) = MC_Prtcl%R0(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
      Rcur(:) = MC_Prtcl%R(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
   else ! core holes are immobile
!       Rcur(:) = MC_Prtcl%R0(:)  ! [A] 
      Rcur(:) = MC_Prtcl%R(:)  ! [A] 
   endif
   ! Cylindrical radius:
   R = SQRT( Rcur(1)*Rcur(1) + Rcur(2)*Rcur(2) )  ! [A]
   ! Find where it is on the grid:
!    call Find_in_array_monoton(numpar%grids(8)%spatial_grid1(:), R, i_arr)  ! module "Little_subroutines"
! !    print*, numpar%grids(8)%spatial_grid1(i_arr), numpar%grids(8)%spatial_grid1(i_arr+1), R
!    ! Define the thickness of the cylindrical layer of the grid:
!    if (i_arr == size(numpar%grids(8)%spatial_grid1)) then   ! the last point on the grid
!       dR2 = numpar%grids(8)%spatial_grid1(i_arr)*numpar%grids(8)%spatial_grid1(i_arr) - &
!             numpar%grids(8)%spatial_grid1(i_arr-1)*numpar%grids(8)%spatial_grid1(i_arr-1)
!    else if (i_arr == 1) then    ! the first point, assuming it starts from zero
!       dR2 = numpar%grids(8)%spatial_grid1(1)*numpar%grids(8)%spatial_grid1(1)
!    else
!       dR2 = numpar%grids(8)%spatial_grid1(i_arr+1)*numpar%grids(8)%spatial_grid1(i_arr+1) - &
!          numpar%grids(8)%spatial_grid1(i_arr)*numpar%grids(8)%spatial_grid1(i_arr)
!    endif
   if (R < numpar%grids(i_ax)%spatial_grid1(1)) then   ! R below lower limit
      i_arr = 1
      ! Define the thickness of the cylindrical layer of the grid:
      dR2 = numpar%grids(i_ax)%spatial_grid1(1) * numpar%grids(i_ax)%spatial_grid1(1)
   else if (R >= numpar%grids(i_ax)%spatial_grid1(N_siz)) then  ! above the max R grid point
      i_arr = N_siz
      ! Define the thickness of the cylindrical layer of the grid:
      dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr) - &
            numpar%grids(i_ax)%spatial_grid1(i_arr-1)*numpar%grids(i_ax)%spatial_grid1(i_arr-1)
   else ! inside the grid
      call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), R, i_arr)  ! module "Little_subroutines"
      i_arr = i_arr + 1 ! assign particle to the end of the interval
      ! Define the thickness of the cylindrical layer of the grid:
      dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr) - &
            numpar%grids(i_ax)%spatial_grid1(i_arr-1)*numpar%grids(i_ax)%spatial_grid1(i_arr-1)
   endif

   ! Define volume:
   dV = g_Pi*dR2*L  ! [A^3]
   ! Add particle into the distribution, normalized per volume to get density:
   Distr_R(i_arr) = Distr_R(i_arr) + 1.0d0/dV    ! add a particle into this array [1/A^3]
   ! And corresponding potential energy density:
   if (MC_Prtcl%valent) then    ! for valence holes, use band gap as potential energy
      E_Distr_R(i_arr) = E_Distr_R(i_arr) + (MC_Prtcl%Ekin + matter%DOS%Egap)/dV    ! add a kinetic and potential energies into this array [eV/A^3]
   else ! core shells
      E_Distr_R(i_arr) = E_Distr_R(i_arr) + matter%Elements(MC_Prtcl%KOA)%Ip(MC_Prtcl%Sh)/dV    ! add a ionization potential [eV/A^3]
   endif
   
   nullify(matter)
end subroutine add_radial_hole

    
subroutine add_RL_hole(used_target, MC_Prtcl, numpar, tim, Distr_RL, E_Distr_RL)
   type(Matter), intent(in), target :: used_target   ! parameters of the target
   type(Hole), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:,:), intent(inout) :: Distr_RL, E_Distr_RL     ! density and energy density distributions
   !-----------------------------
   integer :: i_arr, j_arr, i_ax, N_siz, N_sizj
   real(8) :: R, dR2, dV, Rcur(3), dL, L
   type(Target_atoms), pointer :: matter
   
   matter => used_target%Material(MC_Prtcl%in_target)
   i_ax = 11
   N_siz = size(numpar%grids(i_ax)%spatial_grid1)
   N_sizj = size(numpar%grids(i_ax)%spatial_grid2)
   
   ! Find where the particle is at the time instant "tim":
   if (MC_Prtcl%valent) then    ! only valent holes can move, not core holes
      Rcur(:) = MC_Prtcl%R(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
   else ! core holes are immobile
      Rcur(:) = MC_Prtcl%R(:)  ! [A] 
   endif
   ! Cylindrical radius:
   R = SQRT( Rcur(1)*Rcur(1) + Rcur(2)*Rcur(2) )  ! [A]
   ! Depth:
   L = Rcur(3)  ![A]
!    ! Find where it is on the grid:
!    call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), R, i_arr)  ! module "Little_subroutines"
!    call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid2(:), L, j_arr)  ! module "Little_subroutines"
! 
!    ! Define the thickness of the cylindrical layer of the grid:
!    if (i_arr == size(numpar%grids(i_ax)%spatial_grid1)) then   ! the last point on the grid
!       dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr) - &
!             numpar%grids(i_ax)%spatial_grid1(i_arr-1)*numpar%grids(i_ax)%spatial_grid1(i_arr-1)
!    else if (i_arr == 1) then    ! the first point, assuming it starts from zero
!       dR2 = numpar%grids(i_ax)%spatial_grid1(1)*numpar%grids(i_ax)%spatial_grid1(1)
!    else
!       dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr+1)*numpar%grids(i_ax)%spatial_grid1(i_arr+1) - &
!          numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr)
!    endif
!      !Define the thickness of the depth layer
!    if (j_arr == size(numpar%grids(i_ax)%spatial_grid2)) then    
!         dL = numpar%grids(i_ax)%spatial_grid2(j_arr) - numpar%grids(i_ax)%spatial_grid2(j_arr-1)
!    else
!         dL = numpar%grids(i_ax)%spatial_grid2(j_arr+1) - numpar%grids(i_ax)%spatial_grid2(j_arr)
!    endif
   
   ! Find where it is on the grid:
   if (R < numpar%grids(i_ax)%spatial_grid1(1)) then   ! R below lower limit
      i_arr = 1
      ! Define the thickness of the cylindrical layer of the grid:
      dR2 = numpar%grids(i_ax)%spatial_grid1(1) * numpar%grids(i_ax)%spatial_grid1(1)
   else if (R >= numpar%grids(i_ax)%spatial_grid1(N_siz)) then  ! above the max R grid point
      i_arr = N_siz
      ! Define the thickness of the cylindrical layer of the grid:
      dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr) - &
            numpar%grids(i_ax)%spatial_grid1(i_arr-1)*numpar%grids(i_ax)%spatial_grid1(i_arr-1)
   else ! inside the grid
      call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), R, i_arr)  ! module "Little_subroutines"
      i_arr = i_arr + 1 ! assign particle to the end of the interval
      ! Define the thickness of the cylindrical layer of the grid:
      dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr) - &
            numpar%grids(i_ax)%spatial_grid1(i_arr-1)*numpar%grids(i_ax)%spatial_grid1(i_arr-1)
   endif
   
   !Define the thickness of the depth layer:
   if (L < numpar%grids(i_ax)%spatial_grid2(1)) then   ! L below lower limit
      j_arr = 1
      dL = numpar%grids(i_ax)%spatial_grid2(j_arr+1) - numpar%grids(i_ax)%spatial_grid2(j_arr)
   else if (L >= numpar%grids(i_ax)%spatial_grid2(N_sizj)) then  ! above the max L grid point
      j_arr = N_sizj
      dL = numpar%grids(i_ax)%spatial_grid2(j_arr) - numpar%grids(i_ax)%spatial_grid2(j_arr-1)
   else ! inside the grid
      call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid2(:), L, j_arr)  ! module "Little_subroutines"
      j_arr = j_arr + 1 ! assign particle to the end of the interval
      dL = numpar%grids(i_ax)%spatial_grid2(j_arr) - numpar%grids(i_ax)%spatial_grid2(j_arr-1)
   endif
   
   ! Define volume:
   dV = g_Pi*dR2*dL  ! [A^3]
   ! Add particle into the distribution, normalized per volume to get density:
   Distr_RL(i_arr, j_arr) = Distr_RL(i_arr, j_arr) + 1.0d0/dV    ! add a particle into this array [1/A^3]
   ! And corresponding potential energy density:
   if (MC_Prtcl%valent) then    ! for valence holes, use band gap as potential energy
      E_Distr_RL(i_arr, j_arr) = E_Distr_RL(i_arr, j_arr) + (MC_Prtcl%Ekin + matter%DOS%Egap)/dV    ! add a kinetic and potential energies into this array [eV/A^3]
   else ! core shells
      E_Distr_RL(i_arr, j_arr) = E_Distr_RL(i_arr, j_arr) + matter%Elements(MC_Prtcl%KOA)%Ip(MC_Prtcl%Sh)/dV    ! add a ionization potential [eV/A^3]
   endif
   
   nullify(matter)
end subroutine add_RL_hole


subroutine add_radial_particle(MC_Prtcl, numpar, tim, Distr_R, E_Distr_R, L, neutral)
   class(Particle), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:), intent(inout) :: Distr_R, E_Distr_R     ! density and energy density distributions
   real(8), intent(in) :: L ! [A] thickness
   logical, intent(in), optional :: neutral ! is it a neutral particle (photon)?
   !-----------------------------
   integer :: i_arr, i_ax, N_siz
   real(8) :: R, dR2, dV, Rcur(3)
   character (15) :: text_prtcl

   i_ax = 8
   N_siz = size(numpar%grids(8)%spatial_grid1)
   
   ! Find where the particle is at the time instant "tim":
!    Rcur(:) = MC_Prtcl%R0(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
   Rcur(:) = MC_Prtcl%R(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
   
   if (present(neutral)) then   ! For photons, make sure periodic boundaries are observed:
      if (neutral) call put_back_into_box(numpar, R = Rcur) ! module "MC_general_tools"
   endif
   
   ! Cylindrical radius:
   R = SQRT( Rcur(1)*Rcur(1) + Rcur(2)*Rcur(2) )  ! [A]
   ! Find where it is on the grid:
!    call Find_in_array_monoton(numpar%grids(8)%spatial_grid1(:), R, i_arr)  ! module "Little_subroutines"
   ! Find where it is on the grid:
   if (R < numpar%grids(i_ax)%spatial_grid1(1)) then   ! R below lower limit
      i_arr = 1
      ! Define the thickness of the cylindrical layer of the grid:
      dR2 = numpar%grids(i_ax)%spatial_grid1(1) * numpar%grids(i_ax)%spatial_grid1(1)
   else if (R >= numpar%grids(i_ax)%spatial_grid1(N_siz)) then  ! above the max R grid point
      i_arr = N_siz
      ! Define the thickness of the cylindrical layer of the grid:
      dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr) - &
            numpar%grids(i_ax)%spatial_grid1(i_arr-1)*numpar%grids(i_ax)%spatial_grid1(i_arr-1)
   else ! inside the grid
      call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), R, i_arr)  ! module "Little_subroutines"
      i_arr = i_arr + 1 ! assign particle to the end of the interval
      ! Define the thickness of the cylindrical layer of the grid:
      dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr) - &
            numpar%grids(i_ax)%spatial_grid1(i_arr-1)*numpar%grids(i_ax)%spatial_grid1(i_arr-1)
   endif
   
!    print*, numpar%grids(8)%spatial_grid1(i_arr), numpar%grids(8)%spatial_grid1(i_arr+1), R
!    ! Define the thickness of the cylindrical layer of the grid:
!    if (i_arr == size(numpar%grids(8)%spatial_grid1)) then   ! the last point on the grid
!       dR2 = numpar%grids(8)%spatial_grid1(i_arr)*numpar%grids(8)%spatial_grid1(i_arr) - &
!             numpar%grids(8)%spatial_grid1(i_arr-1)*numpar%grids(8)%spatial_grid1(i_arr-1)
!    else if (i_arr == 1) then    ! the first point, assuming it starts from zero
!       dR2 = numpar%grids(8)%spatial_grid1(1)*numpar%grids(8)%spatial_grid1(1)
!    else
!       dR2 = numpar%grids(8)%spatial_grid1(i_arr+1)*numpar%grids(8)%spatial_grid1(i_arr+1) - &
!          numpar%grids(8)%spatial_grid1(i_arr)*numpar%grids(8)%spatial_grid1(i_arr)
!    endif

   ! Define volume:
   dV = g_Pi*dR2*L
   ! Add particle into the distribution, normalized per volume to get density:
   Distr_R(i_arr) = Distr_R(i_arr) + 1.0d0/dV    ! add a particle into this array [1/A^3]
   ! And corresponding energy density:
   E_Distr_R(i_arr) = E_Distr_R(i_arr) + MC_Prtcl%Ekin/dV    ! add a particle energy into this array [eV/A^3]

   if ((MC_Prtcl%Ekin < 0.0d0) .or. (dV < 0.0d0)) then
      print*, 'Problem in (add_radial_particle):'
      select type(MC_Prtcl)
      type is (Photon)
         text_prtcl = 'Photon'
      type is (Electron)
         text_prtcl = 'Electron'
      type is (Hole)    ! it must be hole, only then attribute "valent" exists
         text_prtcl = 'Hole'
      type is (SHI)
         text_prtcl = 'SHI'
      end select
      print*, trim(adjustl(text_prtcl))//' with negative energy found: ', MC_Prtcl%Ekin
      print*, MC_Prtcl%active
      print*, MC_Prtcl%R(:)
      print*, MC_Prtcl%R0(:)
      print*, MC_Prtcl%V(:)
      print*, MC_Prtcl%V0(:)
      print*, MC_Prtcl%ti - MC_Prtcl%t0, MC_Prtcl%ti, MC_Prtcl%t0
   endif

end subroutine add_radial_particle


subroutine add_RL_particle(MC_Prtcl, numpar, tim, Distr_RL, E_Distr_RL, neutral)
   class(Particle), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:,:), intent(inout) :: Distr_RL, E_Distr_RL     ! density and energy density distributions
   logical, intent(in), optional :: neutral ! is it a neutral particle (photon)?
   !-----------------------------
   integer :: i_arr, j_arr, i_ax, N_siz, N_sizj
   real(8) :: R, dR2, dV, Rcur(3), dL, L
   
   i_ax = 11
   N_siz = size(numpar%grids(i_ax)%spatial_grid1)
   N_sizj = size(numpar%grids(i_ax)%spatial_grid2)

   ! Find where the particle is at the time instant "tim":
   Rcur(:) = MC_Prtcl%R(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
   if (present(neutral)) then ! for photons, make sure periodic boundaries are observed:
      if (neutral) call put_back_into_box(numpar, R = Rcur) ! module "MC_general_tools"
   endif
   
   ! Cylindrical radius:
   R = SQRT( Rcur(1)*Rcur(1) + Rcur(2)*Rcur(2) )  ! [A]
   ! Depth L:
   L = Rcur(3) ! [A]
   if(L .LT. 0) then
       print*, 'L = ', L
       pause 'Error in subroutine add_RL_particle, module MC_data_analysis'
   endif    
   ! Find where it is on the grid:
!    call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), R, i_arr)  ! module "Little_subroutines"
   if (R < numpar%grids(i_ax)%spatial_grid1(1)) then   ! R below lower limit
      i_arr = 1
      dR2 = numpar%grids(i_ax)%spatial_grid1(1) * numpar%grids(i_ax)%spatial_grid1(1)
   else if (R >= numpar%grids(i_ax)%spatial_grid1(N_siz)) then  ! above the max R grid point
      i_arr = N_siz
      dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr) - &
            numpar%grids(i_ax)%spatial_grid1(i_arr-1)*numpar%grids(i_ax)%spatial_grid1(i_arr-1)
   else ! inside the grid
      call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid1(:), R, i_arr)  ! module "Little_subroutines"
      i_arr = i_arr + 1 ! assign particle to the end of the interval
      dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr) - &
            numpar%grids(i_ax)%spatial_grid1(i_arr-1)*numpar%grids(i_ax)%spatial_grid1(i_arr-1)
   endif
   
   
   ! Include only particles inside borders of L grid:
   !if ((L < numpar%grids(i_ax)%spatial_grid2(1)) .AND. &
   !    (L > numpar%grids(i_ax)%spatial_grid2(size(numpar%grids(i_ax)%spatial_grid2)))) then
       ! Find where it is on the grid:
        !call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid2(:), L, j_arr)  ! module "Little_subroutines"
   ! Find where it is on the grid:
   if (L < numpar%grids(i_ax)%spatial_grid2(1)) then   ! L below lower limit
      j_arr = 1
      dL = numpar%grids(i_ax)%spatial_grid2(j_arr+1) - numpar%grids(i_ax)%spatial_grid2(j_arr)
   else if (L >= numpar%grids(i_ax)%spatial_grid2(N_sizj)) then  ! above the max L grid point
      j_arr = N_sizj
      dL = numpar%grids(i_ax)%spatial_grid2(j_arr) - numpar%grids(i_ax)%spatial_grid2(j_arr-1)
   else ! inside the grid
      call Find_in_array_monoton(numpar%grids(i_ax)%spatial_grid2(:), L, j_arr)  ! module "Little_subroutines"
      j_arr = j_arr + 1 ! assign particle to the end of the interval
      dL = numpar%grids(i_ax)%spatial_grid2(j_arr) - numpar%grids(i_ax)%spatial_grid2(j_arr-1)
   endif
   !endif
        
!    ! Define the thickness of the cylindrical layer of the grid:
!    if (i_arr == N_siz) then   ! the last point on the grid
!       dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr) - &
!             numpar%grids(i_ax)%spatial_grid1(i_arr-1)*numpar%grids(i_ax)%spatial_grid1(i_arr-1)
!    else if (i_arr == 1) then    ! the first point, assuming it starts from zero
!       dR2 = numpar%grids(i_ax)%spatial_grid1(1)*numpar%grids(i_ax)%spatial_grid1(1)
!    else
!       dR2 = numpar%grids(i_ax)%spatial_grid1(i_arr+1)*numpar%grids(i_ax)%spatial_grid1(i_arr+1) - &
!          numpar%grids(i_ax)%spatial_grid1(i_arr)*numpar%grids(i_ax)%spatial_grid1(i_arr)
!    endif
!    !Define the thickness of the depth layer
!    if (j_arr == size(numpar%grids(i_ax)%spatial_grid2)) then        
!         dL = numpar%grids(i_ax)%spatial_grid2(j_arr) - numpar%grids(i_ax)%spatial_grid2(j_arr-1)
!    else
!         dL = numpar%grids(i_ax)%spatial_grid2(j_arr+1) - numpar%grids(i_ax)%spatial_grid2(j_arr)
!    endif

   ! Define volume:
   dV = g_Pi*dR2*dL
   ! Add particle into the distribution, normalized per volume to get density:
   Distr_RL(i_arr, j_arr) = Distr_RL(i_arr, j_arr) + 1.0d0/dV    ! add a particle into this array [1/A^3]
   ! And corresponding energy density:
   E_Distr_RL(i_arr, j_arr) = E_Distr_RL(i_arr, j_arr) + MC_Prtcl%Ekin/dV    ! add a particle energy into this array [eV/A^3]
end subroutine add_RL_particle



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! Velosity distributions:

subroutine get_photon_velotheta(MC, numpar, Vel_theta_ph)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:), intent(inout) :: Vel_theta_ph     ! photon spectrum
   !------------------------------
   real(8) :: d_theta, theta, one_over_N, over_sin_theta
   integer :: i, Nsiz, i_arr
   logical :: anything_to_do
   
   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%vel_theta_grid_par%along_axis) then ! velosity theta data
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Photons(:)%active)

      ACTPAR:if (anything_to_do) then ! there are particles to distribute
!          Nsiz = size(MC%MC_Photons(:))
         
         if (MC%N_ph > 0) one_over_N = 1.0d0 / dble(MC%N_ph) ! to normalize to the number of particles
         
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_ph
            ! Include only active particles:
            if (MC%MC_Photons(i)%active) then
!                print*, MC%N_e, i, theta

               ! Find velosity theta [deg]:
               theta = get_v_theta(MC%MC_Photons(i)%V)   ! module "Geometris"
               ! The distribution by theta needs to be converted into pherical coordiante
               ! to convert from distribution by theta in Cartesian to
               ! theta in Spherical, one need additionally to divide by sin(theta), see page 8 [1]:
               if (abs(theta) > m_tollerance_eps) then
                  over_sin_theta = 1.0d0/sin(theta)
               else
                  over_sin_theta = 1.0d0
               endif
               over_sin_theta = over_sin_theta * m_one_over_halfPi  ! normalization accounting for sin(x)
               
               ! Find where to put in on the given theta grid:
               if (theta < numpar%vel_theta_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  d_theta = numpar%vel_theta_grid(1)
               else if (theta >= numpar%vel_theta_grid(size(numpar%vel_theta_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%vel_theta_grid)
                  d_theta = numpar%vel_theta_grid(i_arr) - numpar%vel_theta_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%vel_theta_grid, theta, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  d_theta = numpar%vel_theta_grid(i_arr) - numpar%vel_theta_grid(i_arr-1)
               endif
                ! add a photon into this array, per energy interval to make distribution:
               Vel_theta_ph(i_arr) = Vel_theta_ph(i_arr) + one_over_N/d_theta * over_sin_theta
!                write(*,'(a,i4,i4,f,f,f)') '(get_photon_velotheta)', i, i_arr, numpar%vel_theta_grid(i_arr), theta, Spectrum_e(i_arr)
            endif
         enddo
         
         ! Normalize per number of particles:
!          Vel_theta_ph = Vel_theta_ph / dble(MC%N_ph)
      endif ACTPAR
   endif SPEC
!    pause 'get_photon_velotheta'
end subroutine get_photon_velotheta


subroutine get_electron_velotheta(MC, numpar, Vel_theta_e)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:), intent(inout) :: Vel_theta_e     ! electron spectrum
   !------------------------------
   real(8) :: d_theta, theta, one_over_N, over_sin_theta
   integer :: i, Nsiz, i_arr
   logical :: anything_to_do
   
   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%vel_theta_grid_par%along_axis) then ! velosity theta data
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Electrons(:)%active)
      ACTPAR:if (anything_to_do) then ! there are particles to distribute
!          Nsiz = size(MC%MC_Electrons(:))

         if (MC%N_e > 0) one_over_N = 1.0d0 / dble(MC%N_e)  ! to normalize to the number of particles
         
         if (MC%N_e - count(MC%MC_Electrons(:)%active) /= 0) then
            print*, 'Problem in get_electron_velotheta:'
            print*, 'number of active electrons is unequal to total number of electrons'
            print*, MC%N_e, count(MC%MC_Electrons(:)%active)
         endif
         
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_e
            ! Include only active particles:
            if (MC%MC_Electrons(i)%active) then
!                print*, MC%N_e, i, theta
               ! Find velosity theta [deg]:
               theta = get_v_theta(MC%MC_Electrons(i)%V)   ! module "Geometris"
               
               ! The distribution by theta needs to be converted into spherical coordiante
               ! to convert from distribution by theta in Cartesian to
               ! theta in Spherical, one need additionally to divide by sin(theta), see page 8 [1]:
               if (abs(theta) > m_tollerance_eps) then
                  over_sin_theta = 1.0d0/sin(theta * g_deg2rad)
               else
                  over_sin_theta = 1.0d0
               endif
               over_sin_theta = over_sin_theta * m_one_over_halfPi  ! normalization accounting for sin(x)
               
!                write(*,'(i4,f,f,f,f)') i, MC%MC_Electrons(i)%V(:), theta
               
               ! Find where to put in on the given theta grid:
               if (theta < numpar%vel_theta_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  d_theta = numpar%vel_theta_grid(1)
               else if (theta >= numpar%vel_theta_grid(size(numpar%vel_theta_grid))) then  ! above the max theta grid point
                  i_arr = size(numpar%vel_theta_grid)
                  d_theta = numpar%vel_theta_grid(i_arr) - numpar%vel_theta_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%vel_theta_grid, theta, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  d_theta = numpar%vel_theta_grid(i_arr) - numpar%vel_theta_grid(i_arr-1)
               endif
               ! add a particle into this array, per theta interval to make distribution
               Vel_theta_e(i_arr) = Vel_theta_e(i_arr) + one_over_N/d_theta * over_sin_theta
!                write(*,'(a,i4,i4,f,f,f,f)') '(get_electron_velotheta)', i, i_arr, numpar%vel_theta_grid(i_arr), theta, Vel_theta_e(i_arr), d_theta
!                pause 'get_electron_velotheta'
            endif
         enddo
         
!          ! Normalize per number of particles:
!         if (MC%N_e > 0) then
!             Vel_theta_e = Vel_theta_e / dble(MC%N_e)
!         else
!           Vel_theta_e = 0.0d0 ! no particle => no distribution
!         endif
      endif ACTPAR
   endif SPEC
!    pause 'get_electron_velotheta'
end subroutine get_electron_velotheta



subroutine get_positron_velotheta(MC, numpar, Vel_theta_p)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:), intent(inout) :: Vel_theta_p     ! positron spectrum
   !------------------------------
   real(8) :: d_theta, theta, one_over_N, over_sin_theta
   integer :: i, Nsiz, i_arr
   logical :: anything_to_do
   
   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%vel_theta_grid_par%along_axis) then ! velosity theta data
      
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Positrons(:)%active)

      ACTPAR:if (anything_to_do) then ! there are particles to distribute
!          Nsiz = size(MC%MC_Positrons(:))
         
         if (MC%N_p > 0) one_over_N = 1.0d0 / dble(MC%N_p)  ! to normalize to the number of particles
         
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_p
            ! Include only active particles:
            if (MC%MC_Positrons(i)%active) then
!                print*, MC%N_e, i, MC%MC_Positrons(i)%Ekin
               ! Find velosity theta [deg]:
               theta = get_v_theta(MC%MC_Positrons(i)%V)   ! module "Geometris"
               ! The distribution by theta needs to be converted into pherical coordiante
               ! to convert from distribution by theta in Cartesian to
               ! theta in Spherical, one need additionally to divide by sin(theta), see page 8 [1]:
               if (abs(theta) > m_tollerance_eps) then
                  over_sin_theta = 1.0d0/sin(theta)
               else
                  over_sin_theta = 1.0d0
               endif
               over_sin_theta = over_sin_theta * m_one_over_halfPi  ! normalization accounting for sin(x)
               
               ! Find where to put in on the given energy grid:
               if (theta < numpar%vel_theta_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  d_theta = numpar%vel_theta_grid(1)
               else if (theta >= numpar%vel_theta_grid(size(numpar%vel_theta_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%vel_theta_grid)
                  d_theta = numpar%vel_theta_grid(i_arr) - numpar%vel_theta_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%vel_theta_grid, theta, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  d_theta = numpar%vel_theta_grid(i_arr) - numpar%vel_theta_grid(i_arr-1)
               endif
                ! add a photon into this array, per energy interval to make distribution:
               Vel_theta_p(i_arr) = Vel_theta_p(i_arr) + one_over_N/d_theta * over_sin_theta
!              write(*,'(a,i4,i4,f,f,f)') '(get_positron_velotheta)', i, i_arr, numpar%vel_theta_grid(i_arr), theta, Spectrum_e(i_arr)
            endif
         enddo
         
         ! Normalize per number of particles:
!          Vel_theta_p = Vel_theta_p / dble(MC%N_p)
      endif ACTPAR
   endif SPEC
!    pause 'get_positron_velotheta'
end subroutine get_positron_velotheta



subroutine get_hole_velotheta(MC, numpar, Vel_theta_h)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:), intent(inout) :: Vel_theta_h     ! positron spectrum
   !------------------------------
   real(8) :: d_theta, theta, one_over_N, over_sin_theta
   integer :: i, Nsiz, i_arr
   logical :: anything_to_do
   
   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%vel_theta_grid_par%along_axis) then ! velosity theta data
      
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Holes(:)%active)

      ACTPAR:if (anything_to_do) then ! there are particles to distribute
!          Nsiz = size(MC%MC_Holes(:))
         
         if (MC%N_h > 0) one_over_N = 1.0d0 / dble(MC%N_h)  ! to normalize to the number of particles
         
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_h
            ! Include only active particles:
            if (MC%MC_Holes(i)%active) then
!                print*, MC%N_e, i, MC%MC_Holes(i)%Ekin
               
               ! Find velosity theta [deg]:
               if (MC%MC_Holes(i)%valent) then   
                  theta = get_v_theta(MC%MC_Holes(i)%V)   ! module "Geometris"
               else ! core holes have no velosity
                  theta = 0.0d0
               endif
               
               ! The distribution by theta needs to be converted into pherical coordiante
               ! to convert from distribution by theta in Cartesian to
               ! theta in Spherical, one need additionally to divide by sin(theta), see page 8 [1]:
               if (abs(theta) > m_tollerance_eps) then
                  over_sin_theta = 1.0d0/sin(theta)
               else
                  over_sin_theta = 1.0d0
               endif
               over_sin_theta = over_sin_theta * m_one_over_halfPi  ! normalization accounting for sin(x)
               
               ! Find where to put in on the given energy grid:
               if (theta < numpar%vel_theta_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  d_theta = numpar%vel_theta_grid(1)
               else if (theta >= numpar%vel_theta_grid(size(numpar%vel_theta_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%vel_theta_grid)
                  d_theta = numpar%vel_theta_grid(i_arr) - numpar%vel_theta_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%vel_theta_grid, theta, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  d_theta = numpar%vel_theta_grid(i_arr) - numpar%vel_theta_grid(i_arr-1)
               endif
                ! add a photon into this array, per energy interval to make distribution:
               Vel_theta_h(i_arr) = Vel_theta_h(i_arr) + one_over_N/d_theta * over_sin_theta
!              write(*,'(a,i4,i4,f,f,f)') '(get_hole_velotheta)', i, i_arr, numpar%vel_theta_grid(i_arr), theta, Spectrum_e(i_arr)
            endif
         enddo
         
         ! Normalize per number of particles:
!          Vel_theta_h = Vel_theta_h / dble(MC%N_h)
      endif ACTPAR
   endif SPEC
!    pause 'get_hole_velotheta'
end subroutine get_hole_velotheta



subroutine get_SHI_velotheta(MC, numpar, Vel_theta_SHI)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:), intent(inout) :: Vel_theta_SHI ! SHI spectrum
   !------------------------------
   real(8) :: d_theta, theta, one_over_N, over_sin_theta
   integer :: i, Nsiz, i_arr
   logical :: anything_to_do
   
   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%NRG_grid_par%along_axis) then ! energy data
      
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_SHIs(:)%active)
   
      ACTPAR:if (anything_to_do) then ! there are particles to distribute
         
         if (MC%N_SHI > 0) one_over_N = 1.0d0 / dble(MC%N_SHI)  ! to normalize to the number of particles
         
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_SHI
            ! Include only active particles:
            if (MC%MC_SHIs(i)%active) then

               ! Find velosity theta:
               theta = get_v_theta(MC%MC_SHIs(i)%V)   ! module "Geometris"
               ! The distribution by theta needs to be converted into pherical coordiante
               ! to convert from distribution by theta in Cartesian to
               ! theta in Spherical, one need additionally to divide by sin(theta), see page 8 [1]:
               if (abs(theta) > m_tollerance_eps) then
                  over_sin_theta = 1.0d0/sin(theta)
               else
                  over_sin_theta = 1.0d0
               endif
               over_sin_theta = over_sin_theta * m_one_over_halfPi  ! normalization accounting for sin(x)
               
               ! Find where to put in on the given energy grid:
               if (theta < numpar%NRG_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  d_theta = numpar%NRG_grid(1)
               else if (theta >= numpar%NRG_grid(size(numpar%NRG_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%NRG_grid)
                  d_theta = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%NRG_grid, theta, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  d_theta = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               endif
               ! add a hole into this array, per energy interval to make distribution:
               Vel_theta_SHI(i_arr) = Vel_theta_SHI(i_arr) + one_over_N/d_theta * over_sin_theta
!                write(*,'(a,i4,i4,f,f,f)') '(get_SHI_velotheta)', i, i_arr, numpar%NRG_grid(i_arr), theta, Spectrum_e(i_arr)
            endif
         enddo
         
         ! Normalize per number of particles:
!          Vel_theta_SHI = Vel_theta_SHI / dble(MC%N_SHI)
!       else ACTPAR
!          Vel_theta_SHI(1) = 0.0d0
      endif ACTPAR
   endif SPEC
!    pause 'get_SHI_velotheta'
end subroutine get_SHI_velotheta


!sssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Spectra:

subroutine get_photon_spectrum(MC, numpar, Spectrum_ph)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:), intent(inout) :: Spectrum_ph     ! electron spectrum
   !------------------------------
   real(8) :: dE, one_over_N
   integer :: i, Nsiz, i_arr
   logical :: anything_to_do
   
   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%NRG_grid_par%along_axis) then ! energy data
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Photons(:)%active)
   
      ACTPAR:if (anything_to_do) then ! there are particles to distribute
!          Nsiz = size(MC%MC_Photons(:))
         
         if (MC%N_ph > 0) one_over_N = 1.0d0 / dble(MC%N_ph)  ! to normalize to the number of particles
         
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_ph
            ! Include only active particles:
            if (MC%MC_Photons(i)%active) then
!                print*, MC%N_e, i, MC%MC_Photons(i)%Ekin
               ! Find where to put in on the given energy grid:
               if (MC%MC_Photons(i)%Ekin < numpar%NRG_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  dE = numpar%NRG_grid(1)
               else if (MC%MC_Photons(i)%Ekin >= numpar%NRG_grid(size(numpar%NRG_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%NRG_grid)
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%NRG_grid, MC%MC_Photons(i)%Ekin, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               endif
               Spectrum_ph(i_arr) = Spectrum_ph(i_arr) + one_over_N/dE    ! add a photon into this array, per energy interval to make distribution
!                write(*,'(a,i4,i4,f,f,f)') '(get_electron_spectrum)', i, i_arr, numpar%NRG_grid(i_arr), MC%MC_Photons(i)%Ekin, Spectrum_e(i_arr)
            endif
         enddo
         
         ! Normalize to the number of particles:
!          Spectrum_ph = Spectrum_ph / dble(MC%N_ph)
      endif ACTPAR
   endif SPEC
!    pause 'get_photon_spectrum'
end subroutine get_photon_spectrum


subroutine get_electron_spectrum(MC, numpar, Spectrum_e)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:), intent(inout) :: Spectrum_e     ! electron spectrum
   !------------------------------
   real(8) :: dE, one_over_N
   integer :: i, Nsiz, i_arr
   logical :: anything_to_do
   
   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%NRG_grid_par%along_axis) then ! energy data
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Electrons(:)%active)
   
      ACTPAR:if (anything_to_do) then ! there are particles to distribute
!          Nsiz = size(MC%MC_Electrons(:))
         
         if (MC%N_e > 0) one_over_N = 1.0d0 / dble(MC%N_e)  ! to normalize to the number of particles
         
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_e
            ! Include only active particles:
            if (MC%MC_Electrons(i)%active) then
!                print*, MC%N_e, i, MC%MC_Electrons(i)%Ekin
!                print*, 'Electrons_s:', allocated(numpar%NRG_grid)
               
               ! Find where to put in on the given energy grid:
               if (MC%MC_Electrons(i)%Ekin < numpar%NRG_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  dE = numpar%NRG_grid(1)
               else if (MC%MC_Electrons(i)%Ekin >= numpar%NRG_grid(size(numpar%NRG_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%NRG_grid)
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%NRG_grid, MC%MC_Electrons(i)%Ekin, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               endif

               Spectrum_e(i_arr) = Spectrum_e(i_arr) + one_over_N/dE    ! add an electron into this array, per energy interval to make distribution
!                write(*,'(a,i4,i4,f,f,f)') '(get_electron_spectrum)', i, i_arr, numpar%NRG_grid(i_arr), MC%MC_Electrons(i)%Ekin, Spectrum_e(i_arr)
            endif
         enddo
         ! Normalize to the number of particles:
!          Spectrum_e = Spectrum_e / dble(MC%N_e)
      endif ACTPAR
   endif SPEC
!    pause 'get_electron_spectrum'
end subroutine get_electron_spectrum




subroutine get_hole_spectrum(MC, numpar, Spectrum_h)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:,:), intent(inout) :: Spectrum_h     ! electron spectrum
   !------------------------------
   real(8) :: dE, one_over_N
   integer :: i, Nsiz, i_targ
   integer :: i_arr
   logical :: anything_to_do

   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%NRG_grid_par%along_axis) then ! energy data
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Holes(:)%active)
   
      ACTPAR:if (anything_to_do) then ! there are particles to distribute
!          Nsiz = size(MC%MC_Holes(:))
         
         if (MC%N_h > 0) one_over_N = 1.0d0 / dble(MC%N_h)  ! to normalize to the number of particles
         
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_h
            ! in which target material this hole is:
            i_targ = MC%MC_Holes(i)%in_target
            ! Include only active particles:
            if (MC%MC_Holes(i)%active) then
!                print*, MC%N_e, i, MC%MC_Holes(i)%Ekin
!                print*, 'Holes_s:', allocated(numpar%NRG_grid_VB)
               
               ! Find where to put in on the given energy grid:
               if (MC%MC_Holes(i)%Ekin < numpar%NRG_grid_VB(1)) then   ! energies below lower limit
                  i_arr = 1
                  dE =numpar%NRG_grid_VB(i_arr+1) - numpar%NRG_grid_VB(i_arr)
               else if (MC%MC_Holes(i)%Ekin >= numpar%NRG_grid_VB(size(numpar%NRG_grid_VB))) then  ! above the max energy grid point
                  i_arr = size(numpar%NRG_grid_VB)
                  dE = numpar%NRG_grid_VB(i_arr) - numpar%NRG_grid_VB(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%NRG_grid_VB, MC%MC_Holes(i)%Ekin, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  dE = numpar%NRG_grid_VB(i_arr) - numpar%NRG_grid_VB(i_arr-1)
               endif
               ! add a hole into this array, per energy interval to make distribution:
               Spectrum_h(i_targ,i_arr) = Spectrum_h(i_targ,i_arr) + one_over_N/dE
               
!              write(*,'(a,i4,i4,f,f,f)') '(get_hole_spectrum)', i, i_arr, numpar%NRG_grid_VB(i_arr), MC%MC_Holes(i)%Ekin, Spectrum_e(i_arr)
            endif
         enddo
         
         ! Normalize to the number of particles:
!          Spectrum_h = Spectrum_h / dble(MC%N_h)
      endif ACTPAR
   endif SPEC

end subroutine get_hole_spectrum



subroutine get_positron_spectrum(MC, numpar, Spectrum_p)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:), intent(inout) :: Spectrum_p     ! electron spectrum
   !------------------------------
   real(8) :: dE, one_over_N
   integer :: i, Nsiz, i_arr
   logical :: anything_to_do
   
   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%NRG_grid_par%along_axis) then ! energy data
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Positrons(:)%active)
   
      ACTPAR:if (anything_to_do) then ! there are particles to distribute
!          Nsiz = size(MC%MC_Positrons(:))
         
         if (MC%N_p > 0) one_over_N = 1.0d0 / dble(MC%N_p)  ! to normalize to the number of particles
         
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_p
            ! Include only active particles:
            if (MC%MC_Positrons(i)%active) then
!                print*, MC%N_e, i, MC%MC_Positrons(i)%Ekin
               ! Find where to put in on the given energy grid:
               if (MC%MC_Positrons(i)%Ekin < numpar%NRG_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  dE = numpar%NRG_grid(1)
               else if (MC%MC_Positrons(i)%Ekin >= numpar%NRG_grid(size(numpar%NRG_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%NRG_grid)
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%NRG_grid, MC%MC_Positrons(i)%Ekin, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               endif
               Spectrum_p(i_arr) = Spectrum_p(i_arr) + one_over_N/dE    ! add a positron into this array, per energy interval to make distribution
!                write(*,'(a,i4,i4,f,f,f)') '(get_positron_spectrum)', i, i_arr, numpar%NRG_grid(i_arr), MC%MC_Electrons(i)%Ekin, Spectrum_e(i_arr)
            endif
         enddo
         
         ! Normalize to the number of particles:
!          Spectrum_p = Spectrum_p / dble(MC%N_p)
      endif ACTPAR
   endif SPEC
!    pause 'get_positron_spectrum'
end subroutine get_positron_spectrum



subroutine get_SHI_spectrum(MC, numpar, Spectrum_SHI)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:), intent(inout) :: Spectrum_SHI ! SHI spectrum
   !------------------------------
   real(8) :: dE, one_over_N
   integer :: i, Nsiz, i_arr
   logical :: anything_to_do
   
   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%NRG_grid_par%along_axis) then ! energy data
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_SHIs(:)%active)
   
      ACTPAR:if (anything_to_do) then ! there are particles to distribute
         
         if (MC%N_SHI > 0) one_over_N = 1.0d0 / dble(MC%N_SHI)  ! to normalize to the number of particles
         
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_SHI
            ! Include only active particles:
            if (MC%MC_SHIs(i)%active) then
!                print*, MC%N_e, i, MC%MC_SHIs(i)%Ekin
!                print*, 'SHI_s:', allocated(numpar%NRG_grid)
               
               ! Find where to put in on the given energy grid:
               if (MC%MC_SHIs(i)%Ekin < numpar%NRG_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  dE = numpar%NRG_grid(1)
               else if (MC%MC_SHIs(i)%Ekin >= numpar%NRG_grid(size(numpar%NRG_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%NRG_grid)
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%NRG_grid, MC%MC_SHIs(i)%Ekin, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               endif
               Spectrum_SHI(i_arr) = Spectrum_SHI(i_arr) + one_over_N/dE    ! add SHI into this array, per energy interval to make distribution
!                write(*,'(a,i4,i4,f,f,f)') '(get_SHI_spectrum)', i, i_arr, numpar%NRG_grid(i_arr), MC%MC_Electrons(i)%Ekin, Spectrum_e(i_arr)
            endif
         enddo
         
         ! Normalize to the number of particles:
!          Spectrum_SHI = Spectrum_SHI / dble(MC%N_SHI)
      endif ACTPAR
   endif SPEC
!    pause 'get_SHI_spectrum'
end subroutine get_SHI_spectrum



subroutine get_muon_spectrum(MC, numpar, Spectrum_mu)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:), intent(inout) :: Spectrum_mu     ! muon spectrum
   !------------------------------
   real(8) :: dE, one_over_N
   integer :: i, Nsiz, i_arr
   logical :: anything_to_do

   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%NRG_grid_par%along_axis) then ! energy data
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Muons(:)%active)

      ACTPAR:if (anything_to_do) then ! there are particles to distribute

         if (MC%N_mu > 0) one_over_N = 1.0d0 / dble(MC%N_mu)  ! to normalize to the number of particles

         ! Go through all muons and distribute them into arrays:
         do i = 1, MC%N_mu
            ! Include only active particles:
            if (MC%MC_Muons(i)%active) then
               ! Find where to put in on the given energy grid:
               if (MC%MC_Muons(i)%Ekin < numpar%NRG_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  dE = numpar%NRG_grid(1)
               else if (MC%MC_Muons(i)%Ekin >= numpar%NRG_grid(size(numpar%NRG_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%NRG_grid)
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%NRG_grid, MC%MC_Muons(i)%Ekin, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               endif

               Spectrum_mu(i_arr) = Spectrum_mu(i_arr) + one_over_N/dE    ! add an electron into this array, per energy interval to make distribution
            endif
         enddo
      endif ACTPAR
   endif SPEC
end subroutine get_muon_spectrum



!sssssssssssssssssssssssssssssssssssssssssssssssss
! Spectra along 1d axis:
subroutine get_spectra_in_space_1d(MC, numpar, tim, Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X, &
                                                    Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y, &
                                                    Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z, &
                                                    Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:,:), intent(inout) :: Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X  ! energy spectra in space along X
   real(8), dimension(:,:), intent(inout) :: Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y  ! energy spectra in space along Y
   real(8), dimension(:,:), intent(inout) :: Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z  ! energy spectra in space along Z
   real(8), dimension(:,:), intent(inout) :: Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R  ! energy spectra in space along R

   ! Along X:
   if (numpar%Spectr_grid_par(1)%along_axis) then  ! collect data
      call get_electron_spectra_1d(MC, numpar, tim, Spectra_e_X) ! below
   endif
   
   ! Along Y:
   if (numpar%Spectr_grid_par(2)%along_axis) then  ! collect data
      call get_electron_spectra_1d(MC, numpar, tim, Spectra_e_Y) ! below
   endif
   
!    print*, 'get_spectra_in_space_1d :', size(Spectra_e_Z,1), size(Spectra_e_Z,2)
   
   ! Along Z:
   if (numpar%Spectr_grid_par(3)%along_axis) then  ! collect data
      call get_electron_spectra_1d(MC, numpar, tim, Spectra_e_Z) ! below 
   endif
   
!    print*, size(Spectra_e_Z,1), size(Spectra_e_Z,1)
!    print*, 'get_spectra_in_space_1d END'
   
end subroutine get_spectra_in_space_1d




subroutine get_electron_spectra_1d(MC, numpar, tim, Spectrum_e)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:,:), intent(inout) :: Spectrum_e     ! electron spectrum
   !------------------------------
   real(8) :: dE
   integer :: i, Nsiz, i_arr, i_space
   logical :: anything_to_do
   
   ! Construct spectrum only if user requested it:
   SPEC:if (numpar%NRG_grid_par%along_axis) then ! energy data are also set
      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Electrons(:)%active)
   
      ACTPAR:if (anything_to_do) then ! there are particles to distribute
         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_e
            ! Include only active particles:
            if (MC%MC_Electrons(i)%active) then
               ! Find where to put in on the given energy grid:
               if (MC%MC_Electrons(i)%Ekin < numpar%NRG_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  dE = numpar%NRG_grid(1)
               else if (MC%MC_Electrons(i)%Ekin >= numpar%NRG_grid(size(numpar%NRG_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%NRG_grid)
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%NRG_grid, MC%MC_Electrons(i)%Ekin, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  !dE = numpar%NRG_grid(i_arr+1) - numpar%NRG_grid(i_arr)
                  dE = numpar%NRG_grid(i_arr) - numpar%NRG_grid(i_arr-1)
               endif
               
               ! Find where the particle is at the time instant "tim":
               ! Along X:
               if (numpar%Spectr_grid_par(1)%along_axis) then  ! collect data
                  i_space = add_cartesian_particle_for_spectra_1d(MC%MC_Electrons(i), numpar, tim, 1)    ! below
                  Spectrum_e(i_arr, i_space) = Spectrum_e(i_arr, i_space) + 1.0d0/dE    ! add an electron into this array, per energy interval to make distribution
               endif
               ! Along Y:
               if (numpar%Spectr_grid_par(2)%along_axis) then  ! collect data
                  i_space = add_cartesian_particle_for_spectra_1d(MC%MC_Electrons(i), numpar, tim, 2)    ! below
                  Spectrum_e(i_arr, i_space) = Spectrum_e(i_arr, i_space) + 1.0d0/dE    ! add an electron into this array, per energy interval to make distribution
               endif
               ! Along Z:
               if (numpar%Spectr_grid_par(3)%along_axis) then  ! collect data
                  i_space = add_cartesian_particle_for_spectra_1d(MC%MC_Electrons(i), numpar, tim, 3)    ! below
!                   print*, 'get_electron_spectra_1d', size(Spectrum_e,1), size(Spectrum_e,2), dE, i_arr, i_space
                  Spectrum_e(i_arr, i_space) = Spectrum_e(i_arr, i_space) + 1.0d0/dE    ! add an electron into this array, per energy interval to make distribution
               endif
               
!                if (i_arr > 1) then  ! Testing
!                   print*, 'spectra_1d', MC%MC_Electrons(i)%Ekin, numpar%NRG_grid(i_arr), numpar%NRG_grid(i_arr-1)
!                   if (i_space > 1) then
!                      print*, MC%MC_Electrons(i)%R(3) + (MC%MC_Electrons(i)%V(3)) * (tim - MC%MC_Electrons(i)%t0) , numpar%Spectr_grid(3)%spatial_grid1(i_space), numpar%Spectr_grid(3)%spatial_grid1(i_space-1)
!                   else
!                      print*, MC%MC_Electrons(i)%R(3) + (MC%MC_Electrons(i)%V(3)) * (tim - MC%MC_Electrons(i)%t0), numpar%Spectr_grid(3)%spatial_grid1(i_space), i_space
!                   endif
!                endif
               
!                write(*,'(a,i4,i4,f,f,f)') '(get_electron_spectra_1d)', i, i_arr, numpar%NRG_grid(i_arr), MC%MC_Electrons(i)%Ekin, Spectrum_e(i_arr)
            endif
         enddo
      endif ACTPAR
   endif SPEC
!    pause 'get_electron_spectra_1d'
end subroutine get_electron_spectra_1d



function add_cartesian_particle_for_spectra_1d(MC_Prtcl, numpar, tim, coord_ind, neutral) result(i_out)
   integer i_out    ! index for the space array
   class(Particle), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   integer, intent(in) :: coord_ind ! 1=X, 2=Y, 3=Z, along which coordinate
   logical, intent(in), optional :: neutral ! is it a neutral particle (photon)?
   !-----------------------------
   integer :: i_arr, Nsiz
   real(8) :: dR, dV, Rcur(3)

!    print*, 'add_cartesian_particle_for_spectra_1d'
   
   ! Find where the particle is at the time instant "tim":
   Rcur(:) = MC_Prtcl%R(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
   if (present(neutral)) then ! for photons, make sure periodic boundaries are observed:
      if (neutral) call put_back_into_box(numpar, R = Rcur) ! module "MC_general_tools"
   endif
   
   ! Find where it is on the grid:
   select case (coord_ind)
   case (1)  ! along X
      if (isnan(Rcur(1))) then
         print*, 'Error in add_cartesian_particle_for_spectra_1d #1', Rcur(1), MC_Prtcl%R(:)
         print*, 'V',  (MC_Prtcl%V(:)), (tim - MC_Prtcl%t0)
         print*, 'R0', MC_Prtcl%R0(:)
         print*, 'E', MC_Prtcl%Ekin
      endif
      Nsiz = size(numpar%Spectr_grid(1)%spatial_grid1)
      if (Rcur(1) < numpar%Spectr_grid(1)%spatial_grid1(1)) then
         i_arr = 1
      elseif (Rcur(1) >= numpar%Spectr_grid(1)%spatial_grid1(Nsiz)) then
         i_arr = Nsiz
      else
         call Find_in_array_monoton(numpar%Spectr_grid(1)%spatial_grid1(:), Rcur(1), i_arr)  ! module "Little_subroutines"
!          i_arr = i_arr + 1
      endif
   case (2)  ! along Y
      if (isnan(Rcur(2))) then
         print*, 'Error in add_cartesian_particle_for_spectra_1d #2', Rcur(2), MC_Prtcl%R(:)
         print*, 'V',  (MC_Prtcl%V(:)), (tim - MC_Prtcl%t0)
         print*, 'R0', MC_Prtcl%R0(:)
         print*, 'E', MC_Prtcl%Ekin
      endif
!       call Find_in_array_monoton(numpar%Spectr_grid(2)%spatial_grid1(:), Rcur(2), i_arr)  ! module "Little_subroutines"
      Nsiz = size(numpar%Spectr_grid(2)%spatial_grid1)
      if (Rcur(2) < numpar%Spectr_grid(2)%spatial_grid1(1)) then
         i_arr = 1
      elseif (Rcur(2) >= numpar%Spectr_grid(2)%spatial_grid1(Nsiz)) then
         i_arr = Nsiz
      else
         call Find_in_array_monoton(numpar%Spectr_grid(2)%spatial_grid1(:), Rcur(2), i_arr)  ! module "Little_subroutines"
!          i_arr = i_arr + 1
      endif
   case (3)  ! along Z
      if (isnan(Rcur(3))) then
         print*, 'Error in add_cartesian_particle_for_spectra_1d #3', Rcur(3), MC_Prtcl%R(:)
         print*, 'V',  (MC_Prtcl%V(:)), (tim - MC_Prtcl%t0)
         print*, 'R0', MC_Prtcl%R0(:)
         print*, 'E', MC_Prtcl%Ekin
      endif
      !call Find_in_array_monoton(numpar%Spectr_grid(3)%spatial_grid1(:), Rcur(3), i_arr)  ! module "Little_subroutines"
      Nsiz = size(numpar%Spectr_grid(3)%spatial_grid1)
      if (Rcur(3) < numpar%Spectr_grid(3)%spatial_grid1(1)) then
         i_arr = 1
      elseif (Rcur(3) >= numpar%Spectr_grid(3)%spatial_grid1(Nsiz)) then
         i_arr = Nsiz
      else
         call Find_in_array_monoton(numpar%Spectr_grid(3)%spatial_grid1(:), Rcur(3), i_arr)  ! module "Little_subroutines"
!          i_arr = i_arr + 1
      endif
   end select
!    print*, numpar%Spectr_grid(3)%spatial_grid1(i_arr), numpar%Spectr_grid(3)%spatial_grid1(i_arr+1), R
   i_out = i_arr
!    print*, 'add_cartesian_particle_for_spectra_1d END'
end function add_cartesian_particle_for_spectra_1d




!ttttttttttttttttttttttttttttttttttttttttttttttttt
! Theta distirbution spatially resolved along 1d axis:
subroutine get_theta_in_space_1d(MC, numpar, tim, Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X, &
                                                    Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y, &
                                                    Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z, &
                                                    Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:,:), intent(inout) :: Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X  ! Theta distribution in space along X
   real(8), dimension(:,:), intent(inout) :: Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y  ! Theta distribution in space along Y
   real(8), dimension(:,:), intent(inout) :: Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z  ! Theta distribution in space along Z
   real(8), dimension(:,:), intent(inout) :: Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R  ! Theta distribution in space along R

   ! Along X:
   if (numpar%Theta_grid_par(1)%along_axis) then  ! collect data
      call get_electron_theta_1d(MC, numpar, tim, Theta_e_X) ! below
   endif

   ! Along Y:
   if (numpar%Theta_grid_par(2)%along_axis) then  ! collect data
      call get_electron_theta_1d(MC, numpar, tim, Theta_e_Y) ! below
   endif

!    print*, 'get_theta_in_space_1d :', size(Spectra_e_Z,1), size(Spectra_e_Z,2)

   ! Along Z:
   if (numpar%Theta_grid_par(3)%along_axis) then  ! collect data
      call get_electron_theta_1d(MC, numpar, tim, Theta_e_Z) ! below
   endif

   ! [ONLY 1d DISTRIBUTION FOR ELECTRONS IS READY SO FAR]

!    print*, size(Spectra_e_Z,1), size(Spectra_e_Z,1)
!    print*, 'get_theta_in_space_1d END'

end subroutine get_theta_in_space_1d


subroutine get_electron_theta_1d(MC, numpar, tim, Theta_e)
   type(MC_arrays), intent(in) :: MC      ! elements of MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   real(8), dimension(:,:), intent(inout) :: Theta_e     ! electron theta distribution
   !------------------------------
   real(8) :: dTheta, theta, over_sin_theta, one_over_N
   integer :: i, Nsiz, i_arr, i_space
   logical :: anything_to_do

   ! Construct theta only if user requested it:
   SPEC:if (numpar%vel_theta_grid_par%along_axis) then ! velosity theta data

      ! Check if there is any active particle:
      anything_to_do = any(MC%MC_Electrons(:)%active)

      ACTPAR:if (anything_to_do) then ! there are particles to distribute

         if (MC%N_e > 0) one_over_N = 1.0d0 / dble(MC%N_e)  ! to normalize to the number of particles

         ! Go through all electrons and distribute them into arrays:
         do i = 1, MC%N_e
            ! Include only active particles:
            if (MC%MC_Electrons(i)%active) then
               ! Find where to put in on the given energy grid:
                              ! Find velosity theta [deg]:
               theta = get_v_theta(MC%MC_Electrons(i)%V)   ! module "Geometris"

               ! The distribution by theta needs to be converted into spherical coordiante
               ! to convert from distribution by theta in Cartesian to
               ! theta in Spherical, one need additionally to divide by sin(theta), see page 8 [1]:
               if (abs(theta) > m_tollerance_eps) then
                  over_sin_theta = 1.0d0/sin(theta * g_deg2rad)
               else
                  over_sin_theta = 1.0d0
               endif
               over_sin_theta = over_sin_theta * m_one_over_halfPi  ! normalization accounting for sin(x)

               ! Find where to put in on the given theta grid:
               if (theta < numpar%vel_theta_grid(1)) then   ! energies below lower limit
                  i_arr = 1
                  dTheta = numpar%vel_theta_grid(1)
               else if (theta >= numpar%vel_theta_grid(size(numpar%vel_theta_grid))) then  ! above the max energy grid point
                  i_arr = size(numpar%vel_theta_grid)
                  dTheta = numpar%vel_theta_grid(i_arr) - numpar%vel_theta_grid(i_arr-1)
               else ! inside the grid
                  call Find_in_array_monoton(numpar%vel_theta_grid, theta, i_arr)  ! module "Little_subroutines"
                  i_arr = i_arr + 1 ! assign particle to the end of the interval
                  dTheta = numpar%vel_theta_grid(i_arr) - numpar%vel_theta_grid(i_arr-1)
               endif

               ! Find where the particle is at the time instant "tim":
               ! Along X:
               if (numpar%Spectr_grid_par(1)%along_axis) then  ! collect data
                  i_space = add_cartesian_particle_for_theta_1d(MC%MC_Electrons(i), numpar, tim, 1)    ! below
                  Theta_e(i_arr, i_space) = Theta_e(i_arr, i_space) + one_over_N/dTheta * over_sin_theta ! add an electron into this array, per energy interval to make distribution
               endif
               ! Along Y:
               if (numpar%Spectr_grid_par(2)%along_axis) then  ! collect data
                  i_space = add_cartesian_particle_for_theta_1d(MC%MC_Electrons(i), numpar, tim, 2)    ! below
                  Theta_e(i_arr, i_space) = Theta_e(i_arr, i_space) + one_over_N/dTheta * over_sin_theta    ! add an electron into this array, per energy interval to make distribution
               endif
               ! Along Z:
               if (numpar%Spectr_grid_par(3)%along_axis) then  ! collect data
                  i_space = add_cartesian_particle_for_theta_1d(MC%MC_Electrons(i), numpar, tim, 3)    ! below
!                   print*, 'get_electron_theta_1d', size(Theta_e,1), size(Theta_e,2), dTheta, i_arr, i_space
                  Theta_e(i_arr, i_space) = Theta_e(i_arr, i_space) + one_over_N/dTheta * over_sin_theta    ! add an electron into this array, per energy interval to make distribution
               endif

!                if (i_arr > 1) then  ! Testing
!                   print*, 'spectra_1d', MC%MC_Electrons(i)%Ekin, numpar%NRG_grid(i_arr), numpar%NRG_grid(i_arr-1)
!                   if (i_space > 1) then
!                      print*, MC%MC_Electrons(i)%R(3) + (MC%MC_Electrons(i)%V(3)) * (tim - MC%MC_Electrons(i)%t0) , numpar%Spectr_grid(3)%spatial_grid1(i_space), numpar%Spectr_grid(3)%spatial_grid1(i_space-1)
!                   else
!                      print*, MC%MC_Electrons(i)%R(3) + (MC%MC_Electrons(i)%V(3)) * (tim - MC%MC_Electrons(i)%t0), numpar%Spectr_grid(3)%spatial_grid1(i_space), i_space
!                   endif
!                endif

!                write(*,'(a,i4,i4,f,f,f)') '(get_electron_theta_1d)', i, i_arr, numpar%NRG_grid(i_arr), MC%MC_Electrons(i)%Ekin, Theta_e(i_arr)
            endif
         enddo
      endif ACTPAR
   endif SPEC
!    pause 'get_electron_theta_1d'
end subroutine get_electron_theta_1d



function add_cartesian_particle_for_theta_1d(MC_Prtcl, numpar, tim, coord_ind, neutral) result(i_out)
   integer i_out    ! index for the space array
   class(Particle), intent(in) :: MC_Prtcl      ! MC array for all particles in one iteration
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(in) :: tim   ! [fs] current time step
   integer, intent(in) :: coord_ind ! 1=X, 2=Y, 3=Z, along which coordinate
   logical, intent(in), optional :: neutral ! is it a neutral particle (photon)?
   !-----------------------------
   integer :: i_arr, Nsiz
   real(8) :: dR, dV, Rcur(3)

!    print*, 'add_cartesian_particle_for_theta_1d'

   ! Find where the particle is at the time instant "tim":
   Rcur(:) = MC_Prtcl%R(:) + (MC_Prtcl%V(:)) * (tim - MC_Prtcl%t0)   ! [A]
   if (present(neutral)) then ! for photons, make sure periodic boundaries are observed:
      if (neutral) call put_back_into_box(numpar, R = Rcur) ! module "MC_general_tools"
   endif

   ! Find where it is on the grid:
   select case (coord_ind)
   case (1)  ! along X
      if (isnan(Rcur(1))) then
         print*, 'Error in add_cartesian_particle_for_theta_1d #1', Rcur(1), MC_Prtcl%R(:)
         print*, 'V',  (MC_Prtcl%V(:)), (tim - MC_Prtcl%t0)
         print*, 'R0', MC_Prtcl%R0(:)
         print*, 'E', MC_Prtcl%Ekin
      endif
      Nsiz = size(numpar%Theta_grid(1)%spatial_grid1)
      if (Rcur(1) < numpar%Theta_grid(1)%spatial_grid1(1)) then
         i_arr = 1
      elseif (Rcur(1) >= numpar%Theta_grid(1)%spatial_grid1(Nsiz)) then
         i_arr = Nsiz
      else
         call Find_in_array_monoton(numpar%Theta_grid(1)%spatial_grid1(:), Rcur(1), i_arr)  ! module "Little_subroutines"
!          i_arr = i_arr + 1
      endif
   case (2)  ! along Y
      if (isnan(Rcur(2))) then
         print*, 'Error in add_cartesian_particle_for_theta_1d #2', Rcur(2), MC_Prtcl%R(:)
         print*, 'V',  (MC_Prtcl%V(:)), (tim - MC_Prtcl%t0)
         print*, 'R0', MC_Prtcl%R0(:)
         print*, 'E', MC_Prtcl%Ekin
      endif
!       call Find_in_array_monoton(numpar%Theta_grid(2)%spatial_grid1(:), Rcur(2), i_arr)  ! module "Little_subroutines"
      Nsiz = size(numpar%Theta_grid(2)%spatial_grid1)
      if (Rcur(2) < numpar%Theta_grid(2)%spatial_grid1(1)) then
         i_arr = 1
      elseif (Rcur(2) >= numpar%Theta_grid(2)%spatial_grid1(Nsiz)) then
         i_arr = Nsiz
      else
         call Find_in_array_monoton(numpar%Theta_grid(2)%spatial_grid1(:), Rcur(2), i_arr)  ! module "Little_subroutines"
!          i_arr = i_arr + 1
      endif
   case (3)  ! along Z
      if (isnan(Rcur(3))) then
         print*, 'Error in add_cartesian_particle_for_theta_1d #3', Rcur(3), MC_Prtcl%R(:)
         print*, 'V',  (MC_Prtcl%V(:)), (tim - MC_Prtcl%t0)
         print*, 'R0', MC_Prtcl%R0(:)
         print*, 'E', MC_Prtcl%Ekin
      endif
      !call Find_in_array_monoton(numpar%Theta_grid(3)%spatial_grid1(:), Rcur(3), i_arr)  ! module "Little_subroutines"
      Nsiz = size(numpar%Theta_grid(3)%spatial_grid1)
      if (Rcur(3) < numpar%Theta_grid(3)%spatial_grid1(1)) then
         i_arr = 1
      elseif (Rcur(3) >= numpar%Theta_grid(3)%spatial_grid1(Nsiz)) then
         i_arr = Nsiz
      else
         call Find_in_array_monoton(numpar%Theta_grid(3)%spatial_grid1(:), Rcur(3), i_arr)  ! module "Little_subroutines"
!          i_arr = i_arr + 1
      endif
   end select
!    print*, numpar%Theta_grid(3)%spatial_grid1(i_arr), numpar%Theta_grid(3)%spatial_grid1(i_arr+1), R
   i_out = i_arr
!    print*, 'add_cartesian_particle_for_spectra_1d END'
end function add_cartesian_particle_for_theta_1d




!eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
! Allocating routines:

subroutine allocate_spectra_arrays(Nsiz, Nsiz_VB, Nspec_siz0, Nspec_siz1, Nspec_siz2, Nspec_siz3, &
    Spectrum_ph, Spectrum_e, Spectrum_p, Spectrum_h, Spectrum_SHI, Spectrum_mu, &
    Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X, Spectra_mu_X, &
    Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y, Spectra_mu_Y, &
    Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z, Spectra_mu_Z, &
    Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R, Spectra_mu_R  )
   integer, dimension(:), intent(in) :: Nsiz, Nsiz_VB, Nspec_siz0, Nspec_siz1, Nspec_siz2, Nspec_siz3    ! sizes
   ! Energy disributions (spectra):
   real(8), dimension(:), allocatable, intent(inout) :: Spectrum_ph, Spectrum_e, Spectrum_p, Spectrum_SHI, Spectrum_mu  ! energy spectra
   real(8), dimension(:,:), allocatable, intent(inout) :: Spectrum_h  ! VB holes in each target
   ! Spectra vs space 1d:
   real(8), dimension(:,:), allocatable, intent(inout) :: Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X, Spectra_mu_X  ! energy spectra in space along X
   real(8), dimension(:,:), allocatable, intent(inout) :: Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y, Spectra_mu_Y  ! energy spectra in space along Y
   real(8), dimension(:,:), allocatable, intent(inout) :: Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z, Spectra_mu_Z  ! energy spectra in space along Z
   real(8), dimension(:,:), allocatable, intent(inout) :: Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R, Spectra_mu_R ! energy spectra in space along R
   
   ! Energy spectra:
   allocate(Spectrum_ph(Nsiz(1)), source = 0.0d0)
   allocate(Spectrum_e(Nsiz(1)), source = 0.0d0)
   allocate(Spectrum_p(Nsiz(1)), source = 0.0d0)
   allocate(Spectrum_h(size(Nsiz_VB),maxval(Nsiz_VB)), source = 0.0d0)
   allocate(Spectrum_SHI(Nsiz(1)), source = 0.0d0)
   allocate(Spectrum_mu(Nsiz(1)), source = 0.0d0)
   ! Spectra vs space 1d:
   allocate(Spectra_ph_X(Nsiz(1), Nspec_siz0(2)), source = 0.0d0)   ! X
   allocate(Spectra_e_X(Nsiz(1), Nspec_siz0(2)), source = 0.0d0)
   allocate(Spectra_p_X(Nsiz(1), Nspec_siz0(2)), source = 0.0d0)
   allocate(Spectra_h_X(Nsiz(1), Nspec_siz0(2)), source = 0.0d0)
   allocate(Spectra_SHI_X(Nsiz(1), Nspec_siz0(2)), source = 0.0d0)
   allocate(Spectra_mu_X(Nsiz(1), Nspec_siz0(2)), source = 0.0d0)
   allocate(Spectra_ph_Y(Nsiz(1), Nspec_siz0(3)), source = 0.0d0)   ! Y
   allocate(Spectra_e_Y(Nsiz(1), Nspec_siz0(3)), source = 0.0d0)
   allocate(Spectra_p_Y(Nsiz(1), Nspec_siz0(3)), source = 0.0d0)
   allocate(Spectra_h_Y(Nsiz(1), Nspec_siz0(3)), source = 0.0d0)
   allocate(Spectra_SHI_Y(Nsiz(1), Nspec_siz0(3)), source = 0.0d0) 
   allocate(Spectra_mu_Y(Nsiz(1), Nspec_siz0(3)), source = 0.0d0)
   allocate(Spectra_ph_Z(Nsiz(1), Nspec_siz0(4)), source = 0.0d0)  ! Z
   allocate(Spectra_e_Z(Nsiz(1), Nspec_siz0(4)), source = 0.0d0)
   allocate(Spectra_p_Z(Nsiz(1), Nspec_siz0(4)), source = 0.0d0)
   allocate(Spectra_h_Z(Nsiz(1), Nspec_siz0(4)), source = 0.0d0)
   allocate(Spectra_SHI_Z(Nsiz(1), Nspec_siz0(4)), source = 0.0d0)
   allocate(Spectra_mu_Z(Nsiz(1), Nspec_siz0(4)), source = 0.0d0)
   allocate(Spectra_ph_R(Nsiz(1), Nspec_siz0(8)), source = 0.0d0)   ! R
   allocate(Spectra_e_R(Nsiz(1), Nspec_siz0(8)), source = 0.0d0)
   allocate(Spectra_p_R(Nsiz(1), Nspec_siz0(8)), source = 0.0d0)
   allocate(Spectra_h_R(Nsiz(1), Nspec_siz0(8)), source = 0.0d0)
   allocate(Spectra_SHI_R(Nsiz(1), Nspec_siz0(8)), source = 0.0d0)
   allocate(Spectra_mu_R(Nsiz(1), Nspec_siz0(8)), source = 0.0d0)
end subroutine allocate_spectra_arrays




subroutine allocate_vel_theta_arrays(Nsiz_vel, Ntheta_siz0, Ntheta_siz1, Ntheta_siz2, Ntheta_siz3, &
            Vel_theta_ph, Vel_theta_e, Vel_theta_p, Vel_theta_h, Vel_theta_SHI, Vel_theta_mu, &
            Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X, Theta_mu_X, &
            Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y, Theta_mu_Y, &
            Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z, Theta_mu_Z, &
            Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R, Theta_mu_R  )
   integer, intent(in) :: Nsiz_vel    ! size
   integer, dimension(:), intent(in) :: Ntheta_siz0, Ntheta_siz1, Ntheta_siz2, Ntheta_siz3
   real(8), dimension(:), allocatable, intent(inout) :: Vel_theta_ph, Vel_theta_e, Vel_theta_p, Vel_theta_h, Vel_theta_SHI, Vel_theta_mu
   ! Theta distribution vs space in 1d:
   real(8), dimension(:,:), allocatable, intent(inout) :: Theta_ph_X, Theta_e_X, Theta_p_X, Theta_h_X, Theta_SHI_X, Theta_mu_X  ! theta distribution in space along X
   real(8), dimension(:,:), allocatable, intent(inout) :: Theta_ph_Y, Theta_e_Y, Theta_p_Y, Theta_h_Y, Theta_SHI_Y, Theta_mu_Y  ! theta distribution in space along Y
   real(8), dimension(:,:), allocatable, intent(inout) :: Theta_ph_Z, Theta_e_Z, Theta_p_Z, Theta_h_Z, Theta_SHI_Z, Theta_mu_Z  ! theta distribution in space along Z
   real(8), dimension(:,:), allocatable, intent(inout) :: Theta_ph_R, Theta_e_R, Theta_p_R, Theta_h_R, Theta_SHI_R, Theta_mu_R  ! theta distribution in space along R
   !------------------------------
   ! Velosity theta distributions:
   allocate(Vel_theta_ph(Nsiz_vel), source = 0.0d0)
   allocate(Vel_theta_e(Nsiz_vel), source = 0.0d0)
   allocate(Vel_theta_p(Nsiz_vel), source = 0.0d0)
   allocate(Vel_theta_h(Nsiz_vel), source = 0.0d0)
   allocate(Vel_theta_SHI(Nsiz_vel), source = 0.0d0)
   allocate(Vel_theta_mu(Nsiz_vel), source = 0.0d0)
   ! Velosity theta distributions vs space 1d:
   allocate(Theta_ph_X(Nsiz_vel, Ntheta_siz0(2)), source = 0.0d0)   ! X
   allocate(Theta_e_X(Nsiz_vel, Ntheta_siz0(2)), source = 0.0d0)
   allocate(Theta_p_X(Nsiz_vel, Ntheta_siz0(2)), source = 0.0d0)
   allocate(Theta_h_X(Nsiz_vel, Ntheta_siz0(2)), source = 0.0d0)
   allocate(Theta_SHI_X(Nsiz_vel, Ntheta_siz0(2)), source = 0.0d0)
   allocate(Theta_mu_X(Nsiz_vel, Ntheta_siz0(2)), source = 0.0d0)
   allocate(Theta_ph_Y(Nsiz_vel, Ntheta_siz0(3)), source = 0.0d0)   ! Y
   allocate(Theta_e_Y(Nsiz_vel, Ntheta_siz0(3)), source = 0.0d0)
   allocate(Theta_p_Y(Nsiz_vel, Ntheta_siz0(3)), source = 0.0d0)
   allocate(Theta_h_Y(Nsiz_vel, Ntheta_siz0(3)), source = 0.0d0)
   allocate(Theta_SHI_Y(Nsiz_vel, Ntheta_siz0(3)), source = 0.0d0)
   allocate(Theta_mu_Y(Nsiz_vel, Ntheta_siz0(3)), source = 0.0d0)
   allocate(Theta_ph_Z(Nsiz_vel, Ntheta_siz0(4)), source = 0.0d0)  ! Z
   allocate(Theta_e_Z(Nsiz_vel, Ntheta_siz0(4)), source = 0.0d0)
   allocate(Theta_p_Z(Nsiz_vel, Ntheta_siz0(4)), source = 0.0d0)
   allocate(Theta_h_Z(Nsiz_vel, Ntheta_siz0(4)), source = 0.0d0)
   allocate(Theta_SHI_Z(Nsiz_vel, Ntheta_siz0(4)), source = 0.0d0)
   allocate(Theta_mu_Z(Nsiz_vel, Ntheta_siz0(4)), source = 0.0d0)
   allocate(Theta_ph_R(Nsiz_vel, Ntheta_siz0(8)), source = 0.0d0)   ! R
   allocate(Theta_e_R(Nsiz_vel, Ntheta_siz0(8)), source = 0.0d0)
   allocate(Theta_p_R(Nsiz_vel, Ntheta_siz0(8)), source = 0.0d0)
   allocate(Theta_h_R(Nsiz_vel, Ntheta_siz0(8)), source = 0.0d0)
   allocate(Theta_SHI_R(Nsiz_vel, Ntheta_siz0(8)), source = 0.0d0)
   allocate(Theta_mu_R(Nsiz_vel, Ntheta_siz0(8)), source = 0.0d0)
end subroutine allocate_vel_theta_arrays



subroutine allocate_output_arrays(Nsiz, Nsiz1, Nsiz2, Nsiz3, Nh_siz, &
    Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_ph_R, Distr_ph_L, Distr_ph_Theta, Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic, &
    Distr_e_X, Distr_e_Y, Distr_e_Z, Distr_e_R, Distr_e_L, Distr_e_Theta, Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic, &
    Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_p_R, Distr_p_L, Distr_p_Theta, Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic, &
    Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_h_R, Distr_h_L, Distr_h_Theta, Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic, &
    Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, Distr_SHI_R, Distr_SHI_L, Distr_SHI_Theta, Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic, &
    Distr_a_X, Distr_a_Y, Distr_a_Z, Distr_a_R, Distr_a_L, Distr_a_Theta, Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic, &
    Distr_mu_X, Distr_mu_Y, Distr_mu_Z, Distr_mu_R, Distr_mu_L, Distr_mu_Theta, Distr_mu_Rc, Distr_mu_Thetac, Distr_mu_Phic, &
    Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ, Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta, &
    Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic, &
    Distr_e_XY, Distr_e_YZ, Distr_e_XZ, Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta, &
    Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic, &
    Distr_p_XY, Distr_p_YZ, Distr_p_XZ, Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta, &
    Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic, &
    Distr_h_XY, Distr_h_YZ, Distr_h_XZ, Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta, &
    Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic, &
    Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ, Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta, &
    Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic, &
    Distr_a_XY, Distr_a_YZ, Distr_a_XZ, Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta, &
    Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic, &
    Distr_mu_XY, Distr_mu_YZ, Distr_mu_XZ, Distr_mu_RL, Distr_mu_RTheta, Distr_mu_LTheta, &
    Distr_mu_RcThc, Distr_mu_RcPhic, Distr_mu_ThcPhic, &
    Distr_ph_XYZ, Distr_ph_RLTheta, Distr_ph_RcThcPhic, &
    Distr_e_XYZ, Distr_e_RLTheta, Distr_e_RcThcPhic, &
    Distr_p_XYZ, Distr_p_RLTheta, Distr_p_RcThcPhic, &
    Distr_h_XYZ, Distr_h_RLTheta, Distr_h_RcThcPhic, &
    Distr_SHI_XYZ, Distr_SHI_RLTheta, Distr_SHI_RcThcPhic, &
    Distr_a_XYZ, Distr_a_RLTheta, Distr_a_RcThcPhic, &
    Distr_mu_XYZ, Distr_mu_RLTheta, Distr_mu_RcThcPhic )
   integer, dimension(:), intent(in) :: Nsiz, Nsiz1, Nsiz2, Nsiz3   ! sizes
   integer, intent(in) :: Nh_siz
   ! Spatial distributions in 1d:
   real(8), dimension(:), allocatable, intent(inout) :: Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_ph_R, Distr_ph_L, Distr_ph_Theta, &
                                                            Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic ! photon
   real(8), dimension(:), allocatable, intent(inout) :: Distr_e_X, Distr_e_Y, Distr_e_Z, Distr_e_R, Distr_e_L, Distr_e_Theta, &
                                                            Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic ! electron
   real(8), dimension(:), allocatable, intent(inout) :: Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_p_R, Distr_p_L, Distr_p_Theta, &
                                                            Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic ! positron
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_h_R, Distr_h_L, Distr_h_Theta, &
                                                            Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic ! hole
   real(8), dimension(:), allocatable, intent(inout) :: Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, Distr_SHI_R, Distr_SHI_L, &
                                                            Distr_SHI_Theta, Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic ! SHI
   real(8), dimension(:), allocatable, intent(inout) :: Distr_a_X, Distr_a_Y, Distr_a_Z, Distr_a_R, Distr_a_L, &
                                                            Distr_a_Theta, Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic ! atom
   real(8), dimension(:), allocatable, intent(inout) :: Distr_mu_X, Distr_mu_Y, Distr_mu_Z, Distr_mu_R, Distr_mu_L, Distr_mu_Theta, &
                                                            Distr_mu_Rc, Distr_mu_Thetac, Distr_mu_Phic ! muon
   ! Spatial distributions in 2d:
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_ph_XY, Distr_ph_YZ, Distr_ph_XZ, Distr_ph_RL, Distr_ph_RTheta, Distr_ph_LTheta
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_ph_RcThc, Distr_ph_RcPhic, Distr_ph_ThcPhic  ! photon
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_e_XY, Distr_e_YZ, Distr_e_XZ, Distr_e_RL, Distr_e_RTheta, Distr_e_LTheta
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_e_RcThc, Distr_e_RcPhic, Distr_e_ThcPhic  ! electron
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_p_XY, Distr_p_YZ, Distr_p_XZ, Distr_p_RL, Distr_p_RTheta, Distr_p_LTheta
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_p_RcThc, Distr_p_RcPhic, Distr_p_ThcPhic  ! positron
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Distr_h_XY, Distr_h_YZ, Distr_h_XZ, Distr_h_RL, Distr_h_RTheta, Distr_h_LTheta
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Distr_h_RcThc, Distr_h_RcPhic, Distr_h_ThcPhic  ! hole
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_SHI_XY, Distr_SHI_YZ, Distr_SHI_XZ, Distr_SHI_RL, Distr_SHI_RTheta, Distr_SHI_LTheta
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_SHI_RcThc, Distr_SHI_RcPhic, Distr_SHI_ThcPhic  ! SHI
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_a_XY, Distr_a_YZ, Distr_a_XZ, Distr_a_RL, Distr_a_RTheta, Distr_a_LTheta
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_a_RcThc, Distr_a_RcPhic, Distr_a_ThcPhic  ! atom
   real(8), dimension(:,:), allocatable, intent(inout) :: Distr_mu_XY, Distr_mu_YZ, Distr_mu_XZ, Distr_mu_RL, Distr_mu_RTheta, Distr_mu_LTheta, &
                                                          Distr_mu_RcThc, Distr_mu_RcPhic, Distr_mu_ThcPhic  ! muon
   ! Spatial distributions in 3d:
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Distr_ph_XYZ, Distr_ph_RLTheta, Distr_ph_RcThcPhic ! photon
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Distr_e_XYZ, Distr_e_RLTheta, Distr_e_RcThcPhic ! electron
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Distr_p_XYZ, Distr_p_RLTheta, Distr_p_RcThcPhic ! positron
   real(8), dimension(:,:,:,:), allocatable, intent(inout) :: Distr_h_XYZ, Distr_h_RLTheta, Distr_h_RcThcPhic ! hole
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Distr_SHI_XYZ, Distr_SHI_RLTheta, Distr_SHI_RcThcPhic ! SHI
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Distr_a_XYZ, Distr_a_RLTheta, Distr_a_RcThcPhic ! atom
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Distr_mu_XYZ, Distr_mu_RLTheta, Distr_mu_RcThcPhic ! muon

   ! 1d data:
   allocate(Distr_ph_X(Nsiz(2)), source = 0.0d0)
   allocate(Distr_e_X(Nsiz(2)), source = 0.0d0)
   allocate(Distr_p_X(Nsiz(2)), source = 0.0d0)
   allocate(Distr_h_X(Nh_siz,Nsiz(2)), source = 0.0d0)
   allocate(Distr_SHI_X(Nsiz(2)), source = 0.0d0)
   allocate(Distr_a_X(Nsiz(2)), source = 0.0d0)
   allocate(Distr_mu_X(Nsiz(2)), source = 0.0d0)
   
   allocate(Distr_ph_Y(Nsiz(3)), source = 0.0d0)
   allocate(Distr_e_Y(Nsiz(3)), source = 0.0d0)
   allocate(Distr_p_Y(Nsiz(3)), source = 0.0d0)
   allocate(Distr_h_Y(Nh_siz,Nsiz(3)), source = 0.0d0)
   allocate(Distr_SHI_Y(Nsiz(3)), source = 0.0d0)
   allocate(Distr_a_Y(Nsiz(3)), source = 0.0d0)
   allocate(Distr_mu_Y(Nsiz(3)), source = 0.0d0)
   
   allocate(Distr_ph_Z(Nsiz(4)), source = 0.0d0)
   allocate(Distr_e_Z(Nsiz(4)), source = 0.0d0)
   allocate(Distr_p_Z(Nsiz(4)), source = 0.0d0)
   allocate(Distr_h_Z(Nh_siz,Nsiz(4)), source = 0.0d0)
   allocate(Distr_SHI_Z(Nsiz(4)), source = 0.0d0)
   allocate(Distr_a_Z(Nsiz(4)), source = 0.0d0)
   allocate(Distr_mu_Z(Nsiz(4)), source = 0.0d0)

   allocate(Distr_ph_R(Nsiz(5)), source = 0.0d0)
   allocate(Distr_e_R(Nsiz(5)), source = 0.0d0)
   allocate(Distr_p_R(Nsiz(5)), source = 0.0d0)
   allocate(Distr_h_R(Nh_siz,Nsiz(5)), source = 0.0d0)
   allocate(Distr_SHI_R(Nsiz(5)), source = 0.0d0)
   allocate(Distr_a_R(Nsiz(5)), source = 0.0d0)
   allocate(Distr_mu_R(Nsiz(5)), source = 0.0d0)
   
   allocate(Distr_ph_L(Nsiz(6)), source = 0.0d0)
   allocate(Distr_e_L(Nsiz(6)), source = 0.0d0)
   allocate(Distr_p_L(Nsiz(6)), source = 0.0d0)
   allocate(Distr_h_L(Nh_siz,Nsiz(6)), source = 0.0d0)
   allocate(Distr_SHI_L(Nsiz(6)), source = 0.0d0)
   allocate(Distr_a_L(Nsiz(6)), source = 0.0d0)
   allocate(Distr_mu_L(Nsiz(6)), source = 0.0d0)
   
   allocate(Distr_ph_Theta(Nsiz(7)), source = 0.0d0)
   allocate(Distr_e_Theta(Nsiz(7)), source = 0.0d0)
   allocate(Distr_p_Theta(Nsiz(7)), source = 0.0d0)
   allocate(Distr_h_Theta(Nh_siz,Nsiz(7)), source = 0.0d0)
   allocate(Distr_SHI_Theta(Nsiz(7)), source = 0.0d0)
   allocate(Distr_a_Theta(Nsiz(7)), source = 0.0d0)
   allocate(Distr_mu_Theta(Nsiz(7)), source = 0.0d0)
   
   allocate(Distr_ph_Rc(Nsiz(8)), source = 0.0d0)
   allocate(Distr_e_Rc(Nsiz(8)), source = 0.0d0)
   allocate(Distr_p_Rc(Nsiz(8)), source = 0.0d0)
   allocate(Distr_h_Rc(Nh_siz,Nsiz(8)), source = 0.0d0)
   allocate(Distr_SHI_Rc(Nsiz(8)), source = 0.0d0)
   allocate(Distr_a_Rc(Nsiz(8)), source = 0.0d0)
   allocate(Distr_mu_Rc(Nsiz(8)), source = 0.0d0)
   
   allocate(Distr_ph_Thetac(Nsiz(9)), source = 0.0d0)
   allocate(Distr_e_Thetac(Nsiz(9)), source = 0.0d0)
   allocate(Distr_p_Thetac(Nsiz(9)), source = 0.0d0)
   allocate(Distr_h_Thetac(Nh_siz,Nsiz(9)), source = 0.0d0)
   allocate(Distr_SHI_Thetac(Nsiz(9)), source = 0.0d0)
   allocate(Distr_a_Thetac(Nsiz(9)), source = 0.0d0)
   allocate(Distr_mu_Thetac(Nsiz(9)), source = 0.0d0)

   allocate(Distr_ph_Phic(Nsiz(10)), source = 0.0d0)
   allocate(Distr_e_Phic(Nsiz(10)), source = 0.0d0)
   allocate(Distr_p_Phic(Nsiz(10)), source = 0.0d0)
   allocate(Distr_h_Phic(Nh_siz,Nsiz(10)), source = 0.0d0)
   allocate(Distr_SHI_Phic(Nsiz(10)), source = 0.0d0)
   allocate(Distr_a_Phic(Nsiz(10)), source = 0.0d0)
   allocate(Distr_mu_Phic(Nsiz(10)), source = 0.0d0)
   
   ! 2d data:
   allocate(Distr_ph_XY(Nsiz1(1),Nsiz2(1)), source = 0.0d0)
   allocate(Distr_e_XY(Nsiz1(1),Nsiz2(1)), source = 0.0d0)
   allocate(Distr_p_XY(Nsiz1(1),Nsiz2(1)), source = 0.0d0)
   allocate(Distr_h_XY(Nh_siz,Nsiz1(1),Nsiz2(1)), source = 0.0d0)
   allocate(Distr_SHI_XY(Nsiz1(1),Nsiz2(1)), source = 0.0d0)
   allocate(Distr_a_XY(Nsiz1(1),Nsiz2(1)), source = 0.0d0)
   allocate(Distr_mu_XY(Nsiz1(1),Nsiz2(1)), source = 0.0d0)
   
   allocate(Distr_ph_XZ(Nsiz1(2),Nsiz2(2)), source = 0.0d0)
   allocate(Distr_e_XZ(Nsiz1(2),Nsiz2(2)), source = 0.0d0)
   allocate(Distr_p_XZ(Nsiz1(2),Nsiz2(2)), source = 0.0d0)
   allocate(Distr_h_XZ(Nh_siz,Nsiz1(2),Nsiz2(2)), source = 0.0d0)
   allocate(Distr_SHI_XZ(Nsiz1(2),Nsiz2(2)), source = 0.0d0)
   allocate(Distr_a_XZ(Nsiz1(2),Nsiz2(2)), source = 0.0d0)
   allocate(Distr_mu_XZ(Nsiz1(2),Nsiz2(2)), source = 0.0d0)
   
   allocate(Distr_ph_YZ(Nsiz1(3),Nsiz2(3)), source = 0.0d0)
   allocate(Distr_e_YZ(Nsiz1(3),Nsiz2(3)), source = 0.0d0)
   allocate(Distr_p_YZ(Nsiz1(3),Nsiz2(3)), source = 0.0d0)
   allocate(Distr_h_YZ(Nh_siz,Nsiz1(3),Nsiz2(3)), source = 0.0d0)
   allocate(Distr_SHI_YZ(Nsiz1(3),Nsiz2(3)), source = 0.0d0)
   allocate(Distr_a_YZ(Nsiz1(3),Nsiz2(3)), source = 0.0d0)
   allocate(Distr_mu_YZ(Nsiz1(3),Nsiz2(3)), source = 0.0d0)

   allocate(Distr_ph_RL(Nsiz1(4),Nsiz2(4)), source = 0.0d0)
   allocate(Distr_e_RL(Nsiz1(4),Nsiz2(4)), source = 0.0d0)
   allocate(Distr_p_RL(Nsiz1(4),Nsiz2(4)), source = 0.0d0)
   allocate(Distr_h_RL(Nh_siz,Nsiz1(4),Nsiz2(4)), source = 0.0d0)
   allocate(Distr_SHI_RL(Nsiz1(4),Nsiz2(4)), source = 0.0d0)
   allocate(Distr_a_RL(Nsiz1(4),Nsiz2(4)), source = 0.0d0)
   allocate(Distr_mu_RL(Nsiz1(4),Nsiz2(4)), source = 0.0d0)

   allocate(Distr_ph_RTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_e_RTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_p_RTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_h_RTheta(Nh_siz,Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_SHI_RTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_a_RTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_mu_RTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   
   allocate(Distr_ph_LTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_e_LTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_p_LTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_h_LTheta(Nh_siz,Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_SHI_LTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_a_LTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   allocate(Distr_mu_LTheta(Nsiz1(5),Nsiz2(5)), source = 0.0d0)
   
   allocate(Distr_ph_RcThc(Nsiz1(6),Nsiz2(6)), source = 0.0d0)
   allocate(Distr_e_RcThc(Nsiz1(6),Nsiz2(6)), source = 0.0d0)
   allocate(Distr_p_RcThc(Nsiz1(6),Nsiz2(6)), source = 0.0d0)
   allocate(Distr_h_RcThc(Nh_siz,Nsiz1(6),Nsiz2(6)), source = 0.0d0)
   allocate(Distr_SHI_RcThc(Nsiz1(6),Nsiz2(6)), source = 0.0d0)
   allocate(Distr_a_RcThc(Nsiz1(6),Nsiz2(6)), source = 0.0d0)
   allocate(Distr_mu_RcThc(Nsiz1(6),Nsiz2(6)), source = 0.0d0)
   
   allocate(Distr_ph_RcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_e_RcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_p_RcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_h_RcPhic(Nh_siz,Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_SHI_RcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_a_RcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_mu_RcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   
   allocate(Distr_ph_ThcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_e_ThcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_p_ThcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_h_ThcPhic(Nh_siz,Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_SHI_ThcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_a_ThcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   allocate(Distr_mu_ThcPhic(Nsiz1(7),Nsiz2(7)), source = 0.0d0)
   
   ! 3d data:
   allocate(Distr_ph_XYZ(Nsiz1(8),Nsiz2(8),Nsiz3(8)), source = 0.0d0)
   allocate(Distr_e_XYZ(Nsiz1(8),Nsiz2(8),Nsiz3(8)), source = 0.0d0)
   allocate(Distr_p_XYZ(Nsiz1(8),Nsiz2(8),Nsiz3(8)), source = 0.0d0)
   allocate(Distr_h_XYZ(Nh_siz,Nsiz1(8),Nsiz2(8),Nsiz3(8)), source = 0.0d0)
   allocate(Distr_SHI_XYZ(Nsiz1(8),Nsiz2(8),Nsiz3(8)), source = 0.0d0)
   allocate(Distr_a_XYZ(Nsiz1(8),Nsiz2(8),Nsiz3(8)), source = 0.0d0)
   allocate(Distr_mu_XYZ(Nsiz1(8),Nsiz2(8),Nsiz3(8)), source = 0.0d0)
   
   allocate(Distr_ph_RLTheta(Nsiz1(9),Nsiz2(9),Nsiz3(9)), source = 0.0d0)
   allocate(Distr_e_RLTheta(Nsiz1(9),Nsiz2(9),Nsiz3(9)), source = 0.0d0)
   allocate(Distr_p_RLTheta(Nsiz1(9),Nsiz2(9),Nsiz3(9)), source = 0.0d0)
   allocate(Distr_h_RLTheta(Nh_siz,Nsiz1(9),Nsiz2(9),Nsiz3(9)), source = 0.0d0)
   allocate(Distr_SHI_RLTheta(Nsiz1(9),Nsiz2(9),Nsiz3(9)), source = 0.0d0)
   allocate(Distr_a_RLTheta(Nsiz1(9),Nsiz2(9),Nsiz3(9)), source = 0.0d0)
   allocate(Distr_mu_RLTheta(Nsiz1(9),Nsiz2(9),Nsiz3(9)), source = 0.0d0)

   allocate(Distr_ph_RcThcPhic(Nsiz1(10),Nsiz2(10),Nsiz3(10)), source = 0.0d0)
   allocate(Distr_e_RcThcPhic(Nsiz1(10),Nsiz2(10),Nsiz3(10)), source = 0.0d0)
   allocate(Distr_p_RcThcPhic(Nsiz1(10),Nsiz2(10),Nsiz3(10)), source = 0.0d0)
   allocate(Distr_h_RcThcPhic(Nh_siz,Nsiz1(10),Nsiz2(10),Nsiz3(10)), source = 0.0d0)
   allocate(Distr_SHI_RcThcPhic(Nsiz1(10),Nsiz2(10),Nsiz3(10)), source = 0.0d0)
   allocate(Distr_a_RcThcPhic(Nsiz1(10),Nsiz2(10),Nsiz3(10)), source = 0.0d0)
   allocate(Distr_mu_RcThcPhic(Nsiz1(10),Nsiz2(10),Nsiz3(10)), source = 0.0d0)
end subroutine allocate_output_arrays



subroutine get_size_VB_output(used_target, N_size_VB, NRG_grid_VB)
   type(Matter), intent(in) :: used_target  ! parameters of the target
   integer, dimension(:), allocatable, intent(out) :: N_size_VB ! size of the grid for VB spectrum
   real(8), dimension(:), allocatable, intent(inout), optional :: NRG_grid_VB  ! energy grid for the VB output data (VB spectra)
   integer :: i, Ntargets, i_DOS
   real(8), parameter :: dE = 0.1d0    ! [eV] step for the grid for VB spectrum

   Ntargets = size(used_target%Material)    ! how many different targets
   ! Allocate sizes for spectra:
   if (allocated(N_size_VB)) deallocate(N_size_VB)
   allocate(N_size_VB(Ntargets)) ! for each target

   ! How many grid points to cover all DOSes:
   i_DOS = CEILING((maxval(used_target%Material(:)%DOS%E_VB_top) - minval(used_target%Material(:)%DOS%E_VB_bottom))/dE)
   N_size_VB = i_DOS    ! set them all equal

   ! Allocate the grid for VB spectra:
   if (present(NRG_grid_VB)) then
      if (.not.allocated(NRG_grid_VB)) then
         allocate(NRG_grid_VB(i_DOS))
         ! Set the grid array values:
         NRG_grid_VB(1) = abs(maxval(used_target%Material(:)%DOS%E_VB_top))
         do i = 2, i_DOS
            NRG_grid_VB(i) = NRG_grid_VB(i-1) + dE
         enddo
      endif
   endif

   ! For all targets, define VB grids according to DOSes:
!    do i = 1, Ntargets
!       N_size_VB = size(used_target%Material(i)%DOS%E)
!    enddo
end subroutine get_size_VB_output



subroutine set_sizes_temp_output(numpar, Nsiz, Nsiz1, Nsiz2, Nsiz3, Nsiz_vel, Nspec_siz0, Nspec_siz1, Nspec_siz2, Nspec_siz3, &
                                    Nsiz_vel_siz0, Nsiz_vel_siz1, Nsiz_vel_siz2, Nsiz_vel_siz3)
   type(Num_par), intent(in) :: numpar   ! all numerical parameters
   integer, intent(out) :: Nsiz_vel
   integer, dimension(:), intent(out) :: Nsiz, Nsiz1, Nsiz2, Nsiz3
   integer, dimension(:), intent(out) :: Nspec_siz0, Nspec_siz1, Nspec_siz2, Nspec_siz3
   integer, dimension(:), intent(out) :: Nsiz_vel_siz0, Nsiz_vel_siz1, Nsiz_vel_siz2, Nsiz_vel_siz3
   !------------------------------
   ! Default values:
   Nsiz = 0
   Nsiz1 = 0
   Nsiz2 = 0
   Nsiz3 = 0
   Nsiz_vel = 0
   Nspec_siz0 = 0
   Nspec_siz1 = 0
   Nspec_siz2 = 0
   Nspec_siz3 = 0
   Nsiz_vel_siz0 = 0
   Nsiz_vel_siz1 = 0
   Nsiz_vel_siz2 = 0
   Nsiz_vel_siz3 = 0

   ! Redefine those that are requested by the user:
   ! Velosity theta distribution:
   if (numpar%vel_theta_grid_par%along_axis) then ! velosity theta data
      Nsiz_vel = size(numpar%vel_theta_grid)
   endif
   ! Energy (spectra):
   if (numpar%NRG_grid_par%along_axis) then ! energy data
      Nsiz(1) = size(numpar%NRG_grid)
   endif
   ! Spectra vs space:
   if (numpar%Spectr_grid_par(1)%along_axis) then ! X
      Nspec_siz0(2) = size(numpar%Spectr_grid(1)%spatial_grid1)
   endif
   if (numpar%Spectr_grid_par(2)%along_axis) then ! Y
      Nspec_siz0(3) = size(numpar%Spectr_grid(2)%spatial_grid1)
   endif
   if (numpar%Spectr_grid_par(3)%along_axis) then ! Z
      Nspec_siz0(4) = size(numpar%Spectr_grid(3)%spatial_grid1)
   endif
   if (numpar%Spectr_grid_par(8)%along_axis) then ! R  (cyllinrical)
      Nspec_siz0(5) = size(numpar%Spectr_grid(8)%spatial_grid1)
   endif
   if (numpar%Spectr_grid_par(9)%along_axis) then ! L
      Nspec_siz0(6) = size(numpar%Spectr_grid(9)%spatial_grid1)
   endif
   if (numpar%Spectr_grid_par(10)%along_axis) then ! Theta (cyllinrical)
      Nspec_siz0(7) = size(numpar%Spectr_grid(10)%spatial_grid1)
   endif
   if (numpar%Spectr_grid_par(14)%along_axis) then ! R (Spherical)
      Nspec_siz0(8) = size(numpar%Spectr_grid(14)%spatial_grid1)
   endif
   if (numpar%Spectr_grid_par(15)%along_axis) then ! Theta (Spherical)
      Nspec_siz0(9) = size(numpar%Spectr_grid(15)%spatial_grid1)
   endif
   if (numpar%Spectr_grid_par(16)%along_axis) then ! Phi (Spherical)
      Nspec_siz0(10) = size(numpar%Spectr_grid(16)%spatial_grid1)
   endif
   if (numpar%Spectr_grid_par(4)%along_axis) then ! XY
      Nspec_siz1(1) = size(numpar%Spectr_grid(4)%spatial_grid1)
      Nspec_siz2(1) = size(numpar%Spectr_grid(4)%spatial_grid2)
   endif
   if (numpar%Spectr_grid_par(5)%along_axis) then ! XZ
      Nspec_siz1(2) = size(numpar%Spectr_grid(5)%spatial_grid1)
      Nspec_siz2(2) = size(numpar%Spectr_grid(5)%spatial_grid2)
   endif
   if (numpar%Spectr_grid_par(6)%along_axis) then ! YZ
      Nspec_siz1(3) = size(numpar%Spectr_grid(6)%spatial_grid1)
      Nspec_siz2(3) = size(numpar%Spectr_grid(6)%spatial_grid2)
   endif
   if (numpar%Spectr_grid_par(11)%along_axis) then ! RL
      Nspec_siz1(4) = size(numpar%Spectr_grid(11)%spatial_grid1)
      Nspec_siz2(4) = size(numpar%Spectr_grid(11)%spatial_grid2)
   endif
   if (numpar%Spectr_grid_par(12)%along_axis) then ! RTheta
      Nspec_siz1(5) = size(numpar%Spectr_grid(12)%spatial_grid1)
      Nspec_siz2(5) = size(numpar%Spectr_grid(12)%spatial_grid2)
   endif
   if (numpar%Spectr_grid_par(17)%along_axis) then ! R Theta    (Spherical)
       Nspec_siz1(6) = size(numpar%Spectr_grid(17)%spatial_grid1)
       Nspec_siz2(6) = size(numpar%Spectr_grid(17)%spatial_grid2)
   endif
   if (numpar%Spectr_grid_par(18)%along_axis) then ! R Phi    (Spherical)
      Nspec_siz1(7) = size(numpar%Spectr_grid(18)%spatial_grid1)
      Nspec_siz2(7) = size(numpar%Spectr_grid(18)%spatial_grid2)
   endif
   if (numpar%Spectr_grid_par(7)%along_axis) then ! X Y Z
      Nspec_siz1(8) = size(numpar%Spectr_grid(7)%spatial_grid1)
      Nspec_siz2(8) = size(numpar%Spectr_grid(7)%spatial_grid2)
      Nspec_siz3(8) = size(numpar%Spectr_grid(7)%spatial_grid3)
   endif
   if (numpar%Spectr_grid_par(13)%along_axis) then ! R L Theta
      Nspec_siz1(9) = size(numpar%Spectr_grid(13)%spatial_grid1)
      Nspec_siz2(9) = size(numpar%Spectr_grid(13)%spatial_grid2)
      Nspec_siz3(9) = size(numpar%Spectr_grid(13)%spatial_grid3)
   endif
   if (numpar%Spectr_grid_par(19)%along_axis) then ! R Theta Phi
      Nspec_siz1(10) = size(numpar%Spectr_grid(19)%spatial_grid1)
      Nspec_siz2(10) = size(numpar%Spectr_grid(19)%spatial_grid2)
      Nspec_siz3(10) = size(numpar%Spectr_grid(19)%spatial_grid3)
   endif

   ! Theta distribution vs space:
   if (numpar%Theta_grid_par(1)%along_axis) then ! X
      Nsiz_vel_siz0(2) = size(numpar%Theta_grid(1)%spatial_grid1)
   endif
   if (numpar%Theta_grid_par(2)%along_axis) then ! Y
      Nsiz_vel_siz0(3) = size(numpar%Theta_grid(2)%spatial_grid1)
   endif
   if (numpar%Theta_grid_par(3)%along_axis) then ! Z
      Nsiz_vel_siz0(4) = size(numpar%Theta_grid(3)%spatial_grid1)
   endif
   if (numpar%Theta_grid_par(8)%along_axis) then ! R  (cyllinrical)
      Nsiz_vel_siz0(5) = size(numpar%Theta_grid(8)%spatial_grid1)
   endif
   if (numpar%Theta_grid_par(9)%along_axis) then ! L
      Nsiz_vel_siz0(6) = size(numpar%Theta_grid(9)%spatial_grid1)
   endif
   if (numpar%Theta_grid_par(10)%along_axis) then ! Theta (cyllinrical)
      Nsiz_vel_siz0(7) = size(numpar%Theta_grid(10)%spatial_grid1)
   endif
   if (numpar%Theta_grid_par(14)%along_axis) then ! R (Spherical)
      Nsiz_vel_siz0(8) = size(numpar%Theta_grid(14)%spatial_grid1)
   endif
   if (numpar%Theta_grid_par(15)%along_axis) then ! Theta (Spherical)
      Nsiz_vel_siz0(9) = size(numpar%Theta_grid(15)%spatial_grid1)
   endif
   if (numpar%Theta_grid_par(16)%along_axis) then ! Phi (Spherical)
      Nsiz_vel_siz0(10) = size(numpar%Theta_grid(16)%spatial_grid1)
   endif
   if (numpar%Theta_grid_par(4)%along_axis) then ! XY
      Nsiz_vel_siz1(1) = size(numpar%Theta_grid(4)%spatial_grid1)
      Nsiz_vel_siz2(1) = size(numpar%Theta_grid(4)%spatial_grid2)
   endif
   if (numpar%Theta_grid_par(5)%along_axis) then ! XZ
      Nsiz_vel_siz1(2) = size(numpar%Theta_grid(5)%spatial_grid1)
      Nsiz_vel_siz2(2) = size(numpar%Theta_grid(5)%spatial_grid2)
   endif
   if (numpar%Theta_grid_par(6)%along_axis) then ! YZ
      Nsiz_vel_siz1(3) = size(numpar%Theta_grid(6)%spatial_grid1)
      Nsiz_vel_siz2(3) = size(numpar%Theta_grid(6)%spatial_grid2)
   endif
   if (numpar%Theta_grid_par(11)%along_axis) then ! RL
      Nsiz_vel_siz1(4) = size(numpar%Theta_grid(11)%spatial_grid1)
      Nsiz_vel_siz2(4) = size(numpar%Theta_grid(11)%spatial_grid2)
   endif
   if (numpar%Theta_grid_par(12)%along_axis) then ! RTheta
      Nsiz_vel_siz1(5) = size(numpar%Theta_grid(12)%spatial_grid1)
      Nsiz_vel_siz2(5) = size(numpar%Theta_grid(12)%spatial_grid2)
   endif
   if (numpar%Theta_grid_par(17)%along_axis) then ! R Theta    (Spherical)
       Nsiz_vel_siz1(6) = size(numpar%Theta_grid(17)%spatial_grid1)
       Nsiz_vel_siz2(6) = size(numpar%Theta_grid(17)%spatial_grid2)
   endif
   if (numpar%Theta_grid_par(18)%along_axis) then ! R Phi    (Spherical)
      Nsiz_vel_siz1(7) = size(numpar%Theta_grid(18)%spatial_grid1)
      Nsiz_vel_siz2(7) = size(numpar%Theta_grid(18)%spatial_grid2)
   endif
   if (numpar%Theta_grid_par(7)%along_axis) then ! X Y Z
      Nsiz_vel_siz1(8) = size(numpar%Theta_grid(7)%spatial_grid1)
      Nsiz_vel_siz2(8) = size(numpar%Theta_grid(7)%spatial_grid2)
      Nsiz_vel_siz3(8) = size(numpar%Theta_grid(7)%spatial_grid3)
   endif
   if (numpar%Theta_grid_par(13)%along_axis) then ! R L Theta
      Nsiz_vel_siz1(9) = size(numpar%Theta_grid(13)%spatial_grid1)
      Nsiz_vel_siz2(9) = size(numpar%Theta_grid(13)%spatial_grid2)
      Nsiz_vel_siz3(9) = size(numpar%Theta_grid(13)%spatial_grid3)
   endif
   if (numpar%Theta_grid_par(19)%along_axis) then ! R Theta Phi
      Nsiz_vel_siz1(10) = size(numpar%Theta_grid(19)%spatial_grid1)
      Nsiz_vel_siz2(10) = size(numpar%Theta_grid(19)%spatial_grid2)
      Nsiz_vel_siz3(10) = size(numpar%Theta_grid(19)%spatial_grid3)
   endif
   
   ! Spatial grids:
   if (numpar%grid_par(1)%along_axis) then ! X
      Nsiz(2) = size(numpar%grids(1)%spatial_grid1)
   endif
   if (numpar%grid_par(2)%along_axis) then ! Y
      Nsiz(3) = size(numpar%grids(2)%spatial_grid1)
   endif
   if (numpar%grid_par(3)%along_axis) then ! Z
      Nsiz(4) = size(numpar%grids(3)%spatial_grid1)
   endif
   if (numpar%grid_par(8)%along_axis) then ! R  (cyllinrical)
      Nsiz(5) = size(numpar%grids(8)%spatial_grid1)
   endif
   if (numpar%grid_par(9)%along_axis) then ! L
      Nsiz(6) = size(numpar%grids(9)%spatial_grid1)
   endif
   if (numpar%grid_par(10)%along_axis) then ! Theta (cyllinrical)
      Nsiz(7) = size(numpar%grids(10)%spatial_grid1)
   endif
   if (numpar%grid_par(14)%along_axis) then ! R (Spherical)
      Nsiz(8) = size(numpar%grids(14)%spatial_grid1)
   endif
   if (numpar%grid_par(15)%along_axis) then ! Theta (Spherical)
      Nsiz(9) = size(numpar%grids(15)%spatial_grid1)
   endif
   if (numpar%grid_par(16)%along_axis) then ! Phi (Spherical)
      Nsiz(10) = size(numpar%grids(16)%spatial_grid1)
   endif
   if (numpar%grid_par(4)%along_axis) then ! XY
      Nsiz1(1) = size(numpar%grids(4)%spatial_grid1)
      Nsiz2(1) = size(numpar%grids(4)%spatial_grid2)
   endif
   if (numpar%grid_par(5)%along_axis) then ! XZ
      Nsiz1(2) = size(numpar%grids(5)%spatial_grid1)
      Nsiz2(2) = size(numpar%grids(5)%spatial_grid2)
   endif
   if (numpar%grid_par(6)%along_axis) then ! YZ
      Nsiz1(3) = size(numpar%grids(6)%spatial_grid1)
      Nsiz2(3) = size(numpar%grids(6)%spatial_grid2)
   endif
   if (numpar%grid_par(11)%along_axis) then ! RL
      Nsiz1(4) = size(numpar%grids(11)%spatial_grid1)
      Nsiz2(4) = size(numpar%grids(11)%spatial_grid2)
   endif
   if (numpar%grid_par(12)%along_axis) then ! RTheta
      Nsiz1(5) = size(numpar%grids(12)%spatial_grid1)
      Nsiz2(5) = size(numpar%grids(12)%spatial_grid2)
   endif
   if (numpar%grid_par(17)%along_axis) then ! R Theta    (Spherical)
       Nsiz1(6) = size(numpar%grids(17)%spatial_grid1)
       Nsiz2(6) = size(numpar%grids(17)%spatial_grid2)
   endif
   if (numpar%grid_par(18)%along_axis) then ! R Phi    (Spherical)
      Nsiz1(7) = size(numpar%grids(18)%spatial_grid1)
      Nsiz2(7) = size(numpar%grids(18)%spatial_grid2)
   endif
   if (numpar%grid_par(7)%along_axis) then ! X Y Z
      Nsiz1(8) = size(numpar%grids(7)%spatial_grid1)
      Nsiz2(8) = size(numpar%grids(7)%spatial_grid2)
      Nsiz3(8) = size(numpar%grids(7)%spatial_grid3)
   endif
   if (numpar%grid_par(13)%along_axis) then ! R L Theta
      Nsiz1(9) = size(numpar%grids(13)%spatial_grid1)
      Nsiz2(9) = size(numpar%grids(13)%spatial_grid2)
      Nsiz3(9) = size(numpar%grids(13)%spatial_grid3)
   endif
   if (numpar%grid_par(19)%along_axis) then ! R Theta Phi
      Nsiz1(10) = size(numpar%grids(19)%spatial_grid1)
      Nsiz2(10) = size(numpar%grids(19)%spatial_grid2)
      Nsiz3(10) = size(numpar%grids(19)%spatial_grid3)
   endif
end subroutine set_sizes_temp_output



subroutine allocate_output(used_target, numpar, out_data)
   type(Matter), intent(in) :: used_target  ! parameters of the target
   type(Num_par), intent(inout), target :: numpar   ! all numerical parameters
   type(output_data), intent(inout) :: out_data  ! all output data (distributions etc.)
   integer :: Nsiz_vel, Nsiz(10), Nsiz1(10), Nsiz2(10), Nsiz3(10), Nspec_siz0(10), Nspec_siz1(10), Nspec_siz2(10), Nspec_siz3(10), &
                  Nsiz_vel_siz0(10), Nsiz_vel_siz1(10), Nsiz_vel_siz2(10), Nsiz_vel_siz3(10)
   integer, dimension(:), allocatable :: Nsiz_VB    ! VB holes for each target material
   integer :: i_DOS
   real(8) :: dE

   out_data%Eph = 0.0d0
   out_data%Ee = 0.0d0
   out_data%Eh_kin = 0.0d0
   out_data%Eh_pot = 0.0d0
   out_data%Ep = 0.0d0
   out_data%Eat = 0.0d0
   ! Energies of high-energy particles (above cut off):
   out_data%Eph_high = 0.0d0
   out_data%Ee_high = 0.0d0
   out_data%Eh_kin_high = 0.0d0
   out_data%Eh_pot_high = 0.0d0
   out_data%Ep_high = 0.0d0
   
   ! Set sizes of the arrays:
   call set_sizes_temp_output(numpar, Nsiz, Nsiz1, Nsiz2, Nsiz3, &
                              Nsiz_vel, Nspec_siz0, Nspec_siz1, Nspec_siz2, Nspec_siz3, &
                              Nsiz_vel_siz0, Nsiz_vel_siz1, Nsiz_vel_siz2, Nsiz_vel_siz3) ! below
   ! Special grid for valence band holes according to DOS:
   call get_size_VB_output(used_target, Nsiz_VB, numpar%NRG_grid_VB)    ! below

   ! Allocate the arrays when needed:
   if (.not.allocated(out_data%Spectrum_ph)) then    ! all spectra
      call allocate_spectra_arrays(Nsiz, Nsiz_VB, Nspec_siz0, Nspec_siz1, Nspec_siz2, Nspec_siz3, &
       out_data%Spectrum_ph, out_data%Spectrum_e, out_data%Spectrum_p, out_data%Spectrum_h, out_data%Spectrum_SHI, out_data%Spectrum_mu, &
       out_data%Spectra_ph_X, out_data%Spectra_e_X, out_data%Spectra_p_X, out_data%Spectra_h_X, out_data%Spectra_SHI_X, out_data%Spectra_mu_X, &
       out_data%Spectra_ph_Y, out_data%Spectra_e_Y, out_data%Spectra_p_Y, out_data%Spectra_h_Y, out_data%Spectra_SHI_Y, out_data%Spectra_mu_Y, &
       out_data%Spectra_ph_Z, out_data%Spectra_e_Z, out_data%Spectra_p_Z, out_data%Spectra_h_Z, out_data%Spectra_SHI_Z, out_data%Spectra_mu_Z, &
       out_data%Spectra_ph_R, out_data%Spectra_e_R, out_data%Spectra_p_R, out_data%Spectra_h_R, out_data%Spectra_SHI_R, out_data%Spectra_mu_R &
            ) ! below
   endif

   ! Allocate velosity theta distributions:
   if (.not.allocated(out_data%Vel_theta_ph)) then   ! all velosity distributions
      call allocate_vel_theta_arrays(Nsiz_vel, Nsiz_vel_siz0, Nsiz_vel_siz1, Nsiz_vel_siz2, Nsiz_vel_siz3, &
            out_data%Vel_theta_ph, out_data%Vel_theta_e, out_data%Vel_theta_p, out_data%Vel_theta_h, out_data%Vel_theta_SHI, out_data%Vel_theta_mu, &
            out_data%Theta_ph_X, out_data%Theta_e_X, out_data%Theta_p_X, out_data%Theta_h_X, out_data%Theta_SHI_X, out_data%Theta_mu_X, &
            out_data%Theta_ph_Y, out_data%Theta_e_Y, out_data%Theta_p_Y, out_data%Theta_h_Y, out_data%Theta_SHI_Y, out_data%Theta_mu_Y, &
            out_data%Theta_ph_Z, out_data%Theta_e_Z, out_data%Theta_p_Z, out_data%Theta_h_Z, out_data%Theta_SHI_Z, out_data%Theta_mu_Z, &
            out_data%Theta_ph_R, out_data%Theta_e_R, out_data%Theta_p_R, out_data%Theta_h_R, out_data%Theta_SHI_R, out_data%Theta_mu_R )    ! below
   endif

   ! Spatial distributions:
   if (.not.allocated(out_data%Distr_ph_X)) then
      call allocate_output_arrays(Nsiz, Nsiz1, Nsiz2, Nsiz3, numpar%N_sh_tot, &
       out_data%Distr_ph_X, out_data%Distr_ph_Y, out_data%Distr_ph_Z, out_data%Distr_ph_R, &
       out_data%Distr_ph_L, out_data%Distr_ph_Theta, out_data%Distr_ph_Rc, out_data%Distr_ph_Thetac, out_data%Distr_ph_Phic, &
       out_data%Distr_e_X, out_data%Distr_e_Y, out_data%Distr_e_Z, out_data%Distr_e_R, out_data%Distr_e_L, &
       out_data%Distr_e_Theta, out_data%Distr_e_Rc, out_data%Distr_e_Thetac, out_data%Distr_e_Phic, &
       out_data%Distr_p_X, out_data%Distr_p_Y, out_data%Distr_p_Z, out_data%Distr_p_R, &
       out_data%Distr_p_L, out_data%Distr_p_Theta, out_data%Distr_p_Rc, out_data%Distr_p_Thetac, out_data%Distr_p_Phic, &
       out_data%Distr_h_X, out_data%Distr_h_Y, out_data%Distr_h_Z, out_data%Distr_h_R, &
       out_data%Distr_h_L, out_data%Distr_h_Theta, out_data%Distr_h_Rc, out_data%Distr_h_Thetac, out_data%Distr_h_Phic, &
       out_data%Distr_SHI_X, out_data%Distr_SHI_Y, out_data%Distr_SHI_Z, out_data%Distr_SHI_R, out_data%Distr_SHI_L, &
       out_data%Distr_SHI_Theta, out_data%Distr_SHI_Rc, out_data%Distr_SHI_Thetac, out_data%Distr_SHI_Phic, &
       out_data%Distr_a_X, out_data%Distr_a_Y, out_data%Distr_a_Z, out_data%Distr_a_R, out_data%Distr_a_L, &
       out_data%Distr_a_Theta, out_data%Distr_a_Rc, out_data%Distr_a_Thetac, out_data%Distr_a_Phic, &
       out_data%Distr_mu_X, out_data%Distr_mu_Y, out_data%Distr_mu_Z, out_data%Distr_mu_R, &
       out_data%Distr_mu_L, out_data%Distr_mu_Theta, out_data%Distr_mu_Rc, out_data%Distr_mu_Thetac, out_data%Distr_mu_Phic, &
       out_data%Distr_ph_XY, out_data%Distr_ph_YZ, out_data%Distr_ph_XZ, out_data%Distr_ph_RL, out_data%Distr_ph_RTheta, &
       out_data%Distr_ph_LTheta, out_data%Distr_ph_RcThc, out_data%Distr_ph_RcPhic, out_data%Distr_ph_ThcPhic, &
       out_data%Distr_e_XY, out_data%Distr_e_YZ, out_data%Distr_e_XZ, out_data%Distr_e_RL, out_data%Distr_e_RTheta, &
       out_data%Distr_e_LTheta, out_data%Distr_e_RcThc, out_data%Distr_e_RcPhic, out_data%Distr_e_ThcPhic, &
       out_data%Distr_p_XY, out_data%Distr_p_YZ, out_data%Distr_p_XZ, out_data%Distr_p_RL, out_data%Distr_p_RTheta, &
       out_data%Distr_p_LTheta, out_data%Distr_p_RcThc, out_data%Distr_p_RcPhic, out_data%Distr_p_ThcPhic, &
       out_data%Distr_h_XY, out_data%Distr_h_YZ, out_data%Distr_h_XZ, out_data%Distr_h_RL, out_data%Distr_h_RTheta, &
       out_data%Distr_h_LTheta, out_data%Distr_h_RcThc, out_data%Distr_h_RcPhic, out_data%Distr_h_ThcPhic, &
       out_data%Distr_SHI_XY, out_data%Distr_SHI_YZ, out_data%Distr_SHI_XZ, out_data%Distr_SHI_RL, &
       out_data%Distr_SHI_RTheta, out_data%Distr_SHI_LTheta, out_data%Distr_SHI_RcThc, &
       out_data%Distr_SHI_RcPhic, out_data%Distr_SHI_ThcPhic, &
       out_data%Distr_a_XY, out_data%Distr_a_YZ, out_data%Distr_a_XZ, out_data%Distr_a_RL, &
       out_data%Distr_a_RTheta, out_data%Distr_a_LTheta, out_data%Distr_a_RcThc, &
       out_data%Distr_a_RcPhic, out_data%Distr_a_ThcPhic, &
       out_data%Distr_mu_XY, out_data%Distr_mu_YZ, out_data%Distr_mu_XZ, out_data%Distr_mu_RL, out_data%Distr_mu_RTheta, &
       out_data%Distr_mu_LTheta, out_data%Distr_mu_RcThc, out_data%Distr_mu_RcPhic, out_data%Distr_mu_ThcPhic, &
       out_data%Distr_ph_XYZ, out_data%Distr_ph_RLTheta, out_data%Distr_ph_RcThcPhic, &
       out_data%Distr_e_XYZ, out_data%Distr_e_RLTheta, out_data%Distr_e_RcThcPhic, &
       out_data%Distr_p_XYZ, out_data%Distr_p_RLTheta, out_data%Distr_p_RcThcPhic, &
       out_data%Distr_h_XYZ, out_data%Distr_h_RLTheta, out_data%Distr_h_RcThcPhic, &
       out_data%Distr_SHI_XYZ, out_data%Distr_SHI_RLTheta, out_data%Distr_SHI_RcThcPhic, &
       out_data%Distr_a_XYZ, out_data%Distr_a_RLTheta, out_data%Distr_a_RcThcPhic, &
       out_data%Distr_mu_XYZ, out_data%Distr_mu_RLTheta, out_data%Distr_mu_RcThcPhic )    ! below
   endif
   ! And for the energy arrays:
   if (.not.allocated(out_data%E_Distr_ph_X)) then
      call allocate_output_arrays(Nsiz, Nsiz1, Nsiz2, Nsiz3, numpar%N_sh_tot, &
       out_data%E_Distr_ph_X, out_data%E_Distr_ph_Y, out_data%E_Distr_ph_Z, out_data%E_Distr_ph_R, &
       out_data%E_Distr_ph_L, out_data%E_Distr_ph_Theta, out_data%E_Distr_ph_Rc, out_data%E_Distr_ph_Thetac, out_data%E_Distr_ph_Phic, &
       out_data%E_Distr_e_X, out_data%E_Distr_e_Y, out_data%E_Distr_e_Z, out_data%E_Distr_e_R, &
       out_data%E_Distr_e_L, out_data%E_Distr_e_Theta, out_data%E_Distr_e_Rc, out_data%E_Distr_e_Thetac, out_data%E_Distr_e_Phic, &
       out_data%E_Distr_p_X, out_data%E_Distr_p_Y, out_data%E_Distr_p_Z, out_data%E_Distr_p_R, &
       out_data%E_Distr_p_L, out_data%E_Distr_p_Theta, out_data%E_Distr_p_Rc, out_data%E_Distr_p_Thetac, out_data%E_Distr_p_Phic, &
       out_data%E_Distr_h_X, out_data%E_Distr_h_Y, out_data%E_Distr_h_Z, out_data%E_Distr_h_R, &
       out_data%E_Distr_h_L, out_data%E_Distr_h_Theta, out_data%E_Distr_h_Rc, out_data%E_Distr_h_Thetac, out_data%E_Distr_h_Phic, &
       out_data%E_Distr_SHI_X, out_data%E_Distr_SHI_Y, out_data%E_Distr_SHI_Z, out_data%E_Distr_SHI_R, &
       out_data%E_Distr_SHI_L, out_data%E_Distr_SHI_Theta, out_data%E_Distr_SHI_Rc, out_data%E_Distr_SHI_Thetac, out_data%E_Distr_SHI_Phic, &
       out_data%E_Distr_a_X, out_data%E_Distr_a_Y, out_data%E_Distr_a_Z, out_data%E_Distr_a_R, &
       out_data%E_Distr_a_L, out_data%E_Distr_a_Theta, out_data%E_Distr_a_Rc, out_data%E_Distr_a_Thetac, out_data%E_Distr_a_Phic, &
       out_data%E_Distr_mu_X, out_data%E_Distr_mu_Y, out_data%E_Distr_mu_Z, out_data%E_Distr_mu_R, &
       out_data%E_Distr_mu_L, out_data%E_Distr_mu_Theta, out_data%E_Distr_mu_Rc, out_data%E_Distr_mu_Thetac, out_data%E_Distr_mu_Phic, &
       out_data%E_Distr_ph_XY, out_data%E_Distr_ph_YZ, out_data%E_Distr_ph_XZ, out_data%E_Distr_ph_RL, &
       out_data%E_Distr_ph_RTheta, out_data%E_Distr_ph_LTheta, &
       out_data%E_Distr_ph_RcThc, out_data%E_Distr_ph_RcPhic, out_data%E_Distr_ph_ThcPhic, &
       out_data%E_Distr_e_XY, out_data%E_Distr_e_YZ, out_data%E_Distr_e_XZ, out_data%E_Distr_e_RL, &
       out_data%E_Distr_e_RTheta, out_data%E_Distr_e_LTheta, &
       out_data%E_Distr_e_RcThc, out_data%E_Distr_e_RcPhic, out_data%E_Distr_e_ThcPhic, &
       out_data%E_Distr_p_XY, out_data%E_Distr_p_YZ, out_data%E_Distr_p_XZ, &
       out_data%E_Distr_p_RL, out_data%E_Distr_p_RTheta, out_data%E_Distr_p_LTheta, &
       out_data%E_Distr_p_RcThc, out_data%E_Distr_p_RcPhic, out_data%E_Distr_p_ThcPhic, &
       out_data%E_Distr_h_XY, out_data%E_Distr_h_YZ, out_data%E_Distr_h_XZ, out_data%E_Distr_h_RL, &
       out_data%E_Distr_h_RTheta, out_data%E_Distr_h_LTheta, &
       out_data%E_Distr_h_RcThc, out_data%E_Distr_h_RcPhic, out_data%E_Distr_h_ThcPhic, &
       out_data%E_Distr_SHI_XY, out_data%E_Distr_SHI_YZ, out_data%E_Distr_SHI_XZ, &
       out_data%E_Distr_SHI_RL, out_data%E_Distr_SHI_RTheta, out_data%E_Distr_SHI_LTheta, &
       out_data%E_Distr_SHI_RcThc, out_data%E_Distr_SHI_RcPhic, out_data%E_Distr_SHI_ThcPhic, &
       out_data%E_Distr_a_XY, out_data%E_Distr_a_YZ, out_data%E_Distr_a_XZ, &
       out_data%E_Distr_a_RL, out_data%E_Distr_a_RTheta, out_data%E_Distr_a_LTheta, &
       out_data%E_Distr_a_RcThc, out_data%E_Distr_a_RcPhic, out_data%E_Distr_a_ThcPhic, &
       out_data%E_Distr_mu_XY, out_data%E_Distr_mu_YZ, out_data%E_Distr_mu_XZ, &
       out_data%E_Distr_mu_RL, out_data%E_Distr_mu_RTheta, out_data%E_Distr_mu_LTheta, &
       out_data%E_Distr_mu_RcThc, out_data%E_Distr_mu_RcPhic, out_data%E_Distr_mu_ThcPhic, &
       out_data%E_Distr_ph_XYZ, out_data%E_Distr_ph_RLTheta, out_data%E_Distr_ph_RcThcPhic, &
       out_data%E_Distr_e_XYZ, out_data%E_Distr_e_RLTheta, out_data%E_Distr_e_RcThcPhic, &
       out_data%E_Distr_p_XYZ, out_data%E_Distr_p_RLTheta, out_data%E_Distr_p_RcThcPhic, &
       out_data%E_Distr_h_XYZ, out_data%E_Distr_h_RLTheta, out_data%E_Distr_h_RcThcPhic, &
       out_data%E_Distr_SHI_XYZ, out_data%E_Distr_SHI_RLTheta, out_data%E_Distr_SHI_RcThcPhic, &
       out_data%E_Distr_a_XYZ, out_data%E_Distr_a_RLTheta, out_data%E_Distr_a_RcThcPhic, &
       out_data%E_Distr_mu_XYZ, out_data%E_Distr_mu_RLTheta, out_data%E_Distr_mu_RcThcPhic  )    ! below
   endif
   
end subroutine allocate_output


subroutine reset_output_arrays(out_data)
   type(output_data), intent(inout) :: out_data  ! all output data (distributions etc.)
    ! Total energies:
    out_data%Eph = 0.0d0
    out_data%Ee = 0.0d0
    out_data%Eh_kin = 0.0d0
    out_data%Eh_pot = 0.0d0
    out_data%Ep = 0.0d0
!     out_data%Eat = 0.0d0  ! do not reset, add up instead
    ! Energies of high-energy particles (above cut off):
    out_data%Eph_high = 0.0d0
    out_data%Ee_high = 0.0d0
    out_data%Eh_kin_high = 0.0d0
    out_data%Eh_pot_high = 0.0d0
    out_data%Ep_high = 0.0d0
    ! Spectra:
    out_data%Spectrum_ph = 0.0d0
    out_data%Spectrum_e = 0.0d0
    out_data%Spectrum_p = 0.0d0
    out_data%Spectrum_h = 0.0d0
    out_data%Spectrum_SHI = 0.0d0
    ! Velosity theta distribution:
    out_data%Vel_theta_ph = 0.0d0
    out_data%Vel_theta_e = 0.0d0
    out_data%Vel_theta_p = 0.0d0
    out_data%Vel_theta_h = 0.0d0
    out_data%Vel_theta_SHI = 0.0d0
    ! Spectra vs space 1d:
    out_data%Spectra_ph_X = 0.0d0
    out_data%Spectra_e_X = 0.0d0
    out_data%Spectra_p_X = 0.0d0
    out_data%Spectra_h_X = 0.0d0
    out_data%Spectra_SHI_X = 0.0d0
    out_data%Spectra_ph_Y = 0.0d0
    out_data%Spectra_e_Y = 0.0d0
    out_data%Spectra_p_Y = 0.0d0
    out_data%Spectra_h_Y = 0.0d0
    out_data%Spectra_SHI_Y = 0.0d0
    out_data%Spectra_ph_Z = 0.0d0
    out_data%Spectra_e_Z = 0.0d0
    out_data%Spectra_p_Z = 0.0d0
    out_data%Spectra_h_Z = 0.0d0
    out_data%Spectra_SHI_Z = 0.0d0
    out_data%Spectra_ph_R = 0.0d0
    out_data%Spectra_e_R = 0.0d0
    out_data%Spectra_p_R = 0.0d0
    out_data%Spectra_h_R = 0.0d0
    out_data%Spectra_SHI_R = 0.0d0
    ! Spatial distributions:
    out_data%Distr_ph_X = 0.0d0
    out_data%Distr_ph_Y = 0.0d0
    out_data%Distr_ph_Z = 0.0d0
    out_data%Distr_ph_R = 0.0d0
    out_data%Distr_ph_L = 0.0d0
    out_data%Distr_ph_Theta = 0.0d0
    out_data%Distr_ph_Rc = 0.0d0
    out_data%Distr_ph_Thetac = 0.0d0
    out_data%Distr_ph_Phic = 0.0d0
    out_data%Distr_e_X = 0.0d0
    out_data%Distr_e_Y = 0.0d0
    out_data%Distr_e_Z = 0.0d0
    out_data%Distr_e_R = 0.0d0
    out_data%Distr_e_L = 0.0d0
    out_data%Distr_e_Theta = 0.0d0
    out_data%Distr_e_Rc = 0.0d0
    out_data%Distr_e_Thetac = 0.0d0
    out_data%Distr_e_Phic = 0.0d0
    out_data%Distr_p_X = 0.0d0
    out_data%Distr_p_Y = 0.0d0
    out_data%Distr_p_Z = 0.0d0
    out_data%Distr_p_R = 0.0d0
    out_data%Distr_p_L = 0.0d0
    out_data%Distr_p_Theta = 0.0d0
    out_data%Distr_p_Rc = 0.0d0
    out_data%Distr_p_Thetac = 0.0d0
    out_data%Distr_p_Phic = 0.0d0
    out_data%Distr_h_X = 0.0d0
    out_data%Distr_h_Y = 0.0d0
    out_data%Distr_h_Z = 0.0d0
    out_data%Distr_h_R = 0.0d0
    out_data%Distr_h_L = 0.0d0
    out_data%Distr_h_Theta = 0.0d0
    out_data%Distr_h_Rc = 0.0d0
    out_data%Distr_h_Thetac = 0.0d0
    out_data%Distr_h_Phic = 0.0d0
    out_data%Distr_SHI_X = 0.0d0
    out_data%Distr_SHI_Y = 0.0d0
    out_data%Distr_SHI_Z = 0.0d0
    out_data%Distr_SHI_R = 0.0d0
    out_data%Distr_SHI_L = 0.0d0
    out_data%Distr_SHI_Theta = 0.0d0
    out_data%Distr_SHI_Rc = 0.0d0
    out_data%Distr_SHI_Thetac = 0.0d0
    out_data%Distr_SHI_Phic = 0.0d0
!     out_data%Distr_a_X = 0.0d0   ! do not reset, add up instead
!     out_data%Distr_a_Y = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_Z = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_R = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_L = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_Theta = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_Rc = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_Thetac = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_Phic = 0.0d0  ! do not reset, add up instead
    out_data%Distr_ph_XY = 0.0d0
    out_data%Distr_ph_YZ = 0.0d0
    out_data%Distr_ph_XZ = 0.0d0
    out_data%Distr_ph_RL = 0.0d0
    out_data%Distr_ph_RTheta = 0.0d0
    out_data%Distr_ph_LTheta = 0.0d0
    out_data%Distr_ph_RcThc = 0.0d0
    out_data%Distr_ph_RcPhic = 0.0d0
    out_data%Distr_ph_ThcPhic = 0.0d0
    out_data%Distr_e_XY = 0.0d0
    out_data%Distr_e_YZ = 0.0d0
    out_data%Distr_e_XZ = 0.0d0
    out_data%Distr_e_RL = 0.0d0
    out_data%Distr_e_RTheta = 0.0d0
    out_data%Distr_e_LTheta = 0.0d0
    out_data%Distr_e_RcThc = 0.0d0
    out_data%Distr_e_RcPhic = 0.0d0
    out_data%Distr_e_ThcPhic = 0.0d0
    out_data%Distr_p_XY = 0.0d0
    out_data%Distr_p_YZ = 0.0d0
    out_data%Distr_p_XZ = 0.0d0
    out_data%Distr_p_RL = 0.0d0
    out_data%Distr_p_RTheta = 0.0d0
    out_data%Distr_p_LTheta = 0.0d0
    out_data%Distr_p_RcThc = 0.0d0
    out_data%Distr_p_RcPhic = 0.0d0
    out_data%Distr_p_ThcPhic = 0.0d0
    out_data%Distr_h_XY = 0.0d0
    out_data%Distr_h_YZ = 0.0d0
    out_data%Distr_h_XZ = 0.0d0
    out_data%Distr_h_RL = 0.0d0
    out_data%Distr_h_RTheta = 0.0d0
    out_data%Distr_h_LTheta = 0.0d0
    out_data%Distr_h_RcThc = 0.0d0
    out_data%Distr_h_RcPhic = 0.0d0
    out_data%Distr_h_ThcPhic = 0.0d0
    out_data%Distr_SHI_XY = 0.0d0
    out_data%Distr_SHI_YZ = 0.0d0
    out_data%Distr_SHI_XZ = 0.0d0
    out_data%Distr_SHI_RL = 0.0d0
    out_data%Distr_SHI_RTheta = 0.0d0
    out_data%Distr_SHI_LTheta = 0.0d0
    out_data%Distr_SHI_RcThc = 0.0d0
    out_data%Distr_SHI_RcPhic = 0.0d0
    out_data%Distr_SHI_ThcPhic = 0.0d0
!     out_data%Distr_a_XY = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_YZ = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_XZ = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_RL = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_RTheta = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_LTheta = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_RcThc = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_RcPhic = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_ThcPhic = 0.0d0  ! do not reset, add up instead
    out_data%Distr_ph_XYZ = 0.0d0
    out_data%Distr_ph_RLTheta = 0.0d0
    out_data%Distr_ph_RcThcPhic = 0.0d0
    out_data%Distr_e_XYZ = 0.0d0
    out_data%Distr_e_RLTheta = 0.0d0
    out_data%Distr_e_RcThcPhic = 0.0d0
    out_data%Distr_p_XYZ = 0.0d0
    out_data%Distr_p_RLTheta = 0.0d0
    out_data%Distr_p_RcThcPhic = 0.0d0
    out_data%Distr_h_XYZ = 0.0d0
    out_data%Distr_h_RLTheta = 0.0d0
    out_data%Distr_h_RcThcPhic = 0.0d0
    out_data%Distr_SHI_XYZ = 0.0d0
    out_data%Distr_SHI_RLTheta = 0.0d0
    out_data%Distr_SHI_RcThcPhic = 0.0d0
!     out_data%Distr_a_XYZ = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_RLTheta = 0.0d0  ! do not reset, add up instead
!     out_data%Distr_a_RcThcPhic = 0.0d0  ! do not reset, add up instead
    ! Energy distributions:
    out_data%E_Distr_ph_X = 0.0d0
    out_data%E_Distr_ph_Y = 0.0d0
    out_data%E_Distr_ph_Z = 0.0d0
    out_data%E_Distr_ph_R = 0.0d0
    out_data%E_Distr_ph_L = 0.0d0
    out_data%E_Distr_ph_Theta = 0.0d0
    out_data%E_Distr_ph_Rc = 0.0d0
    out_data%E_Distr_ph_Thetac = 0.0d0
    out_data%E_Distr_ph_Phic = 0.0d0
    out_data%E_Distr_e_X = 0.0d0
    out_data%E_Distr_e_Y = 0.0d0
    out_data%E_Distr_e_Z = 0.0d0
    out_data%E_Distr_e_R = 0.0d0
    out_data%E_Distr_e_L = 0.0d0
    out_data%E_Distr_e_Theta = 0.0d0
    out_data%E_Distr_e_Rc = 0.0d0
    out_data%E_Distr_e_Thetac = 0.0d0
    out_data%E_Distr_e_Phic = 0.0d0
    out_data%E_Distr_p_X = 0.0d0
    out_data%E_Distr_p_Y = 0.0d0
    out_data%E_Distr_p_Z = 0.0d0
    out_data%E_Distr_p_R = 0.0d0
    out_data%E_Distr_p_L = 0.0d0
    out_data%E_Distr_p_Theta = 0.0d0
    out_data%E_Distr_p_Rc = 0.0d0
    out_data%E_Distr_p_Thetac = 0.0d0
    out_data%E_Distr_p_Phic = 0.0d0
    out_data%E_Distr_h_X = 0.0d0
    out_data%E_Distr_h_Y = 0.0d0
    out_data%E_Distr_h_Z = 0.0d0
    out_data%E_Distr_h_R = 0.0d0
    out_data%E_Distr_h_L = 0.0d0
    out_data%E_Distr_h_Theta = 0.0d0
    out_data%E_Distr_h_Rc = 0.0d0
    out_data%E_Distr_h_Thetac = 0.0d0
    out_data%E_Distr_h_Phic = 0.0d0
    out_data%E_Distr_SHI_X = 0.0d0
    out_data%E_Distr_SHI_Y = 0.0d0
    out_data%E_Distr_SHI_Z = 0.0d0
    out_data%E_Distr_SHI_R = 0.0d0
    out_data%E_Distr_SHI_L = 0.0d0
    out_data%E_Distr_SHI_Theta = 0.0d0
    out_data%E_Distr_SHI_Rc = 0.0d0
    out_data%E_Distr_SHI_Thetac = 0.0d0
    out_data%E_Distr_SHI_Phic = 0.0d0
!     out_data%E_Distr_a_X = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_Y = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_Z = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_R = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_L = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_Theta = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_Rc = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_Thetac = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_Phic = 0.0d0  ! do not reset, add up instead
    out_data%E_Distr_ph_XY = 0.0d0
    out_data%E_Distr_ph_YZ = 0.0d0
    out_data%E_Distr_ph_XZ = 0.0d0
    out_data%E_Distr_ph_RL = 0.0d0
    out_data%E_Distr_ph_RTheta = 0.0d0
    out_data%E_Distr_ph_LTheta = 0.0d0
    out_data%E_Distr_ph_RcThc = 0.0d0
    out_data%E_Distr_ph_RcPhic = 0.0d0
    out_data%E_Distr_ph_ThcPhic = 0.0d0
    out_data%E_Distr_e_XY = 0.0d0
    out_data%E_Distr_e_YZ = 0.0d0
    out_data%E_Distr_e_XZ = 0.0d0
    out_data%E_Distr_e_RL = 0.0d0
    out_data%E_Distr_e_RTheta = 0.0d0
    out_data%E_Distr_e_LTheta = 0.0d0
    out_data%E_Distr_e_RcThc = 0.0d0
    out_data%E_Distr_e_RcPhic = 0.0d0
    out_data%E_Distr_e_ThcPhic = 0.0d0
    out_data%E_Distr_p_XY = 0.0d0
    out_data%E_Distr_p_YZ = 0.0d0
    out_data%E_Distr_p_XZ = 0.0d0
    out_data%E_Distr_p_RL = 0.0d0
    out_data%E_Distr_p_RTheta = 0.0d0
    out_data%E_Distr_p_LTheta = 0.0d0
    out_data%E_Distr_p_RcThc = 0.0d0
    out_data%E_Distr_p_RcPhic = 0.0d0
    out_data%E_Distr_p_ThcPhic = 0.0d0
    out_data%E_Distr_h_XY = 0.0d0
    out_data%E_Distr_h_YZ = 0.0d0
    out_data%E_Distr_h_XZ = 0.0d0
    out_data%E_Distr_h_RL = 0.0d0
    out_data%E_Distr_h_RTheta = 0.0d0
    out_data%E_Distr_h_LTheta = 0.0d0
    out_data%E_Distr_h_RcThc = 0.0d0
    out_data%E_Distr_h_RcPhic = 0.0d0
    out_data%E_Distr_h_ThcPhic = 0.0d0
    out_data%E_Distr_SHI_XY = 0.0d0
    out_data%E_Distr_SHI_YZ = 0.0d0
    out_data%E_Distr_SHI_XZ = 0.0d0
    out_data%E_Distr_SHI_RL = 0.0d0
    out_data%E_Distr_SHI_RTheta = 0.0d0
    out_data%E_Distr_SHI_LTheta = 0.0d0
    out_data%E_Distr_SHI_RcThc = 0.0d0
    out_data%E_Distr_SHI_RcPhic = 0.0d0
    out_data%E_Distr_SHI_ThcPhic = 0.0d0
!     out_data%E_Distr_a_XY = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_YZ = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_XZ = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_RL = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_RTheta = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_LTheta = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_RcThc = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_RcPhic = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_ThcPhic = 0.0d0  ! do not reset, add up instead
    out_data%E_Distr_ph_XYZ = 0.0d0
    out_data%E_Distr_ph_RLTheta = 0.0d0
    out_data%E_Distr_ph_RcThcPhic = 0.0d0
    out_data%E_Distr_e_XYZ = 0.0d0
    out_data%E_Distr_e_RLTheta = 0.0d0
    out_data%E_Distr_e_RcThcPhic = 0.0d0
    out_data%E_Distr_p_XYZ = 0.0d0
    out_data%E_Distr_p_RLTheta = 0.0d0
    out_data%E_Distr_p_RcThcPhic = 0.0d0
    out_data%E_Distr_h_XYZ = 0.0d0
    out_data%E_Distr_h_RLTheta = 0.0d0
    out_data%E_Distr_h_RcThcPhic = 0.0d0
    out_data%E_Distr_SHI_XYZ = 0.0d0
    out_data%E_Distr_SHI_RLTheta = 0.0d0
    out_data%E_Distr_SHI_RcThcPhic = 0.0d0
!     out_data%E_Distr_a_XYZ = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_RLTheta = 0.0d0  ! do not reset, add up instead
!     out_data%E_Distr_a_RcThcPhic = 0.0d0  ! do not reset, add up instead
end subroutine reset_output_arrays



end module MC_data_analysis
