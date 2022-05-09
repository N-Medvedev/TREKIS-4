! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev, R. Rymzhanov, F. Akhmetov
! in 2018-2021
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains all defined objects

MODULE Objects

use Universal_constants

 implicit none



!============================================== 
! Different geometries used to define targets:
type Geom_type	! basic geometric type
   real(8) :: X, Y, Z	! coordinates of its center
end type Geom_type

type, EXTENDS (Geom_type) :: Rectangle ! it is defined by its cartesian coordinates:
   real(8) :: Xstart, Xend	! Length: its beginning and its end along X axis
   real(8) :: Ystart, Yend	! Length: its beginning and its end along Y axis
   real(8) :: Zstart, Zend	! Length: its beginning and its end along Z axis
   real(8) :: angle_x	! its rotation around X-axis
   real(8) :: angle_y	! its rotation around Y-axis
   real(8) :: angle_z	! its rotation around Z-axis
end type Rectangle 

type, EXTENDS (Geom_type) :: Sphere ! it is defined by its spherical coordinates:
   real(8) :: R	! its radius
end type Sphere 

type, EXTENDS (Geom_type) :: Sphere_segment ! it is defined by its spherical coordinates:
   real(8) :: R_start, R_end		! starting and ending points by radius (if it's cut as a spherical layer)
   real(8) :: phi_start, phi_end	! starting and ending points by phi angle within [-Pi/2..Pi/2] (if it's cut as excluded cone)
   real(8) :: theta_start, theta_end	! starting and ending points by theta angle within [0..2*Pi) (if it's cut as a segment)
end type Sphere_segment 

type, EXTENDS (Geom_type) :: Cylinder ! it is defined by its cylindrical coordinates:
   real(8) :: R	! radius
   real(8) :: L_start, L_end	! starting and ending points of the cylinder defining its length
   real(8) :: angle_x	! its rotation around X-axis
   real(8) :: angle_y	! its rotation around Y-axis
   real(8) :: angle_z	! its rotation around Z-axis
end type Cylinder

type, EXTENDS (Geom_type) :: Cylinder_segment ! it is defined by its cylindrical coordinates:
   real(8) :: R_start, R_end	! starting and ending radius (if it's cut as cylindrical layer)
   real(8) :: L_start, L_end	! starting and ending points of the cylinder defining its length
   real(8) :: theta_start, theta_end	! starting and ending points by theta angle within [0..2*Pi) (if it's cut as a segment)
   real(8) :: angle_x	! its rotation around X-axis
   real(8) :: angle_y	! its rotation around Y-axis
   real(8) :: angle_z	! its rotation around Z-axis
end type Cylinder_segment


! Encapsulating format:
type Geom_array
   class(Geom_type), allocatable :: Set	! geometries of different components of target
end type Geom_array


!==============================================
! Complex dielectric function as an object:
type :: Ritchi_CDF 
   ! parameters entering Ritchi-CDF:
   real(8), dimension(:), allocatable :: E0		! parameter E0
   real(8), dimension(:), allocatable :: A			! parameter A
   real(8), dimension(:), allocatable :: Gamma	! parameter gamma
   ! parameters entering delta-CDF:
   real(8), dimension(:), allocatable :: alpha  ! normalization coefficient
   ! parameters entering single-pole delta-CDF:
   real(8), dimension(:), allocatable :: h_omega_e2    ! [eV^2] (h*Omega_e)^2 where Omega_e is the electron plasma frequency
end type Ritchi_CDF

!============================================== 
! Precalculated cross section arrays:
type :: Cross_section
   real(8), dimension(:), allocatable :: E              ! [eV] energy grid
   real(8), dimension(:), allocatable :: Total         ! [A^2] total cross section
   real(8), dimension(:,:), allocatable :: Per_shell    ! number of shell,  [A^2] cross section for chosen shell
   real(8), dimension(:,:), allocatable :: MFP        ! [A] mean free path per shell
   real(8), dimension(:,:), allocatable :: Se           ! [eV/A] energy loss per shell
   real(8), dimension(:), allocatable :: Total_MFP      ! [A] total mean free path
   real(8), dimension(:), allocatable :: Total_Se         !  [eV/A] total energy loss 
   real(8), dimension(:), allocatable :: Total_Range    !  [A] total range
end type Cross_section

!============================================== 
! Parameters of atoms and ions:
type Atom_kind
   character(3) Name	! Chemical element name
   integer :: Zat		! atomic number [electron charge]
   real(8) :: Mass	! [a.m.u]
   real(8) :: M		! [kg]
   integer :: N_shl	! number of shells of the element
   integer :: N_core_shl    ! number of core shells in this element
   real(8) :: percentage	! contribution to the compound
   integer :: NVB		! number of valence electrons according to periodic table
   real(8) :: Pair_R, Pair_nu_inf	! coefficients needed for photon pair creation calculations
   real(8) :: form_a(5)	! coefficients needed to construct atomic form factors
   integer, dimension(:), allocatable :: Shl_dsgnr	! EADL shell designator
   character(11), dimension(:), allocatable :: Shell_name	! names of the shells
   logical, dimension(:), allocatable :: valent	! mark whether this shell is valent (true) or core (false)
   real(8), dimension(:), allocatable :: Ip		! [eV] ionization potentials for all shells
   real(8), dimension(:), allocatable :: Ek	! [eV] mean kinetic energy of all shells
   real(8), dimension(:), allocatable :: Ne_shell	! number of electrons in each shell
   real(8), dimension(:), allocatable :: Radiat	! [fs] radiative-decay times for all shells
   real(8), dimension(:,:), allocatable :: f_rad	! probability of radiative decay from a given shell
   real(8), dimension(:), allocatable :: Auger	! [fs] Auger-decay times for all shells
   real(8), dimension(:), allocatable :: Compton	! Compton profiles for each orbital
   real(8), dimension(:,:,:), allocatable :: f_auger	! probability of Auger decay between two given shells
   type(Ritchi_CDF), dimension(:), allocatable :: CDF	! complex dielectric function parameters for each shell
   ! Cross sections of scattering on this kind of atom:
   type(Cross_section) :: Phot_absorption	! cross sections of photon absorption
   type(Cross_section) :: Phot_compton	! cross sections of photon Compton scattering
   type(Cross_section) :: Phot_pair		! cross sections of pair creation by photon scattering on nucleus
   type(Cross_section) :: Phot_triplet		! cross sections of pair creation by photon scattering on electron
   type(Cross_section) :: Phot_coherent	! cross sections of coherent photon (elastic, Raylegh, Thompson) scattering
   type(Cross_section) :: El_inelastic	! cross sections of electron inelastic scattering
   type(Cross_section) :: El_elastic	! cross sections of electron elastic scattering
   type(Cross_section) :: El_brems	! cross sections of electron Bremsstrahlung
   type(Cross_section), dimension(:), allocatable :: SHI_inelastic	! cross sections of SHIs inelastic scattering
   type(Cross_section) :: Pos_inelastic	! cross sections of positron inelastic scattering
   type(Cross_section) :: Pos_elastic	! cross sections of positron elastic scattering
   type(Cross_section) :: Pos_brems	! cross sections of positron Bremsstrahlung
   type(Cross_section) :: Pos_annihil	! cross sections of positron annihilation
   type(Cross_section) :: H_inelastic	! cross sections of electron inelastic scattering
   type(Cross_section) :: H_elastic	! cross sections of electron elastic scattering
end type Atom_kind

   
!============================================== 
! Barrier of electron emission from the surface of the material:
type :: Emission_barrier
   integer :: barr_type ! type of potential barrier at the surface: 0=step, 1=Eckart
   real(8) :: Work_func		! [eV] Work function
   real(8) :: Surf_bar		! [A] surface barrier length
   real(8) :: Bar_height	! [eV] barrier height for electron emission
   real(8) :: gamma         ! a parameter entering electron transmission probability
   real(8) :: E1            ! a parameter entering electron transmission probability
end type Emission_barrier


!============================================== 
! Density of states of the material:
type :: Density_of_states
   real(8) :: Egap      ! Band gap [eV]
   real(8) :: E_f        ! Fermi energy [eV]
   real(8) :: v_f         ! Fermi velosity -- units as used in one-pole approximation (not SI!)
   real(8) :: E_VB_bottom, E_VB_top ! bottom and top of the valence band [eV], must be in accord with gap and Fermi energies
   real(8) :: E_CB_bottom, E_CB_top ! bottom and top of the conduction band [eV], must be in accord with gap and Fermi energies
   integer :: N_VB_top, N_CB_bottom ! indices of the arrays in DOS for the top of VB and bottom of CB
   real(8) :: alpha_CB  ! coefficient for the free-electron DOS that matches the given DOS at its top
   real(8), dimension(:), allocatable :: E          ! [eV] energy grid within the valence band
   real(8), dimension(:), allocatable :: DOS     ! [1/eV] DOS within the valence band
   real(8), dimension(:), allocatable :: k          ! Wave vector [1/m] 
   real(8), dimension(:), allocatable :: Eff_m   ! [me] Corresponding effective mass of particle in free particle dispersion relation
end type Density_of_states


!============================================== 
! Encapsulating format for all properties of a given target:
type Target_atoms
   character(100) :: Name   ! name of the material of this constituent of the target
   character(100) :: DOS_file   ! name of the file with the valence band DOS
   integer :: N_Elements    ! number of elements in this compound
   real(8) :: Dens          ! Density [g/cm^3]
   real(8) :: At_Dens   ! Atomic density [atoms/cm^3]
   real(8) :: Mean_Mass ! Average atomic mass [kg]
   real(8) :: Mean_Z    ! Average atomic number
   real(8) :: T         ! Temperature [K]
   real(8) :: T_eV      ! Temperature [eV]
   real(8) :: v_sound   ! Speed of sound [m/s]
   real(8) :: E_debye   ! [eV] Debye energy (maximum phonon energy within Debye approximation)
   real(8) :: E_eistein ! [eV] Einstein phonon energy
   real(8) :: N_VB_el   ! number of valence electrons in VB
   real(8) :: me_eff    ! effective mass of electrons in CB for calculation of single-pole CDF
   type(Density_of_states) :: DOS    ! Valence and/or conduction band DOS of the material
   real(8), dimension(:), allocatable :: fe ! electron distribution function
   real(8), dimension(:), allocatable :: Integral_DOS_fe    ! Boltzmann-like integral used to select energy within DOS
   type(Emission_barrier) :: Surface_barrier    ! Parameters of the surface barier for a particle emission from the surface of the material
   type(Atom_kind), dimension(:), allocatable :: Elements   ! which atoms this target is constructed of
   type(Ritchi_CDF) :: CDF_valence      ! complex dielectric function parameters for phonons
   type(Ritchi_CDF) :: CDF_phonon       ! complex dielectric function parameters for phonons
   type(Cross_section) :: Ph_absorption_valent  ! cross sections of photon absorption by the valence band
   type(Cross_section) :: Ph_absorption_total     ! total cross sections of photon absorption in this target
   type(Cross_section) :: Ph_coherent_total       ! total cross sections of coherent photon scattering in this target
   type(Cross_section) :: Ph_Compton_total        ! total cross sections of Compton photon scattering in this target
   type(Cross_section) :: Ph_pair_total           ! total cross sections of pair creation by photon scattering in this target
   type(Cross_section) :: El_inelastic_valent     ! cross sections of electron inelastic scattering on the valence band
   type(Cross_section) :: El_inelastic_total      ! total cross sections of electron inelastic scattering in this target
   type(Cross_section) :: El_elastic_total        ! total cross sections of electron elastic scattering in this target
   type(Cross_section) :: El_Brems_total          ! total cross sections of electron Bremsstrahlung in this target
   type(Cross_section), dimension(:), allocatable :: SHI_inelastic_valent   ! cross sections of SHIs inelastic scattering on the valence band
   type(Cross_section), dimension(:), allocatable :: SHI_inelastic_total    ! total cross sections of SHIs inelastic scattering in this target
   type(Cross_section) :: Pos_inelastic_valent    ! cross sections of positron inelastic scattering on the valence band
   type(Cross_section) :: Pos_inelastic_total     ! total cross sections of positron inelastic scattering in this target
   type(Cross_section) :: Pos_elastic_total       ! total cross sections of positron elastic scattering in this target
   type(Cross_section) :: Pos_Brems_total         ! total cross sections of positron Bremsstrahlung in this target
   type(Cross_section) :: Pos_annihil_total       ! total cross sections of positron annihilation in this target
   type(Cross_section) :: H_inelastic_total       ! total cross sections of valence hole inelast. scatter. in target (= valence band scattering)
   type(Cross_section) :: H_elastic_total         ! total cross sections of valence hole elastic scattering in this target
end type Target_atoms


!============================================== 
! Parameters of the target material:
type Matter
   character(100) :: Name   ! Name of the target to be used for output directories and files
   integer :: NOC           ! number of different components of the target (layers, clusters, etc.)
   integer, dimension(:), allocatable :: Mat_types  ! indexes of the target geometry constituents:
   ! Mat_types = 0 is a rectangle
   ! Mat_types = 1 is a sphere
   ! Mat_types = 2 is a spherical segment
   ! Mat_types = 3 is a cylinder
   ! Mat_types = 4 is a cylindrical segment
   ! Arrays of all types of geometries:
   type(Geom_array), dimension(:), allocatable :: Geom           ! which geometry each target has
   type(Target_atoms), dimension(:), allocatable :: Material     ! material parameters of each target that it's constructed of
end type Matter

!============================================== 
! All parameters of incomming radiation:
type Radiation_param
   integer :: NOP			! Number of particles in the incomming pulse / bunch
   integer :: KOP			! Kind of particle: 0=photon, 1=electron, 2=positron, 3=SHI
   real(8) :: R(3)			! [A] coordinates of impact
   real(8) :: R_spread(3)	! [A] spread (or uncertainty) of coordinates of impact
!    integer :: Space_shape	! Spatial shape: 0 = rectangular, 1 = Gaussian, 2 = SASE
   real(8) :: theta, phi		! [degrees] theta and phi angles of impact
   real(8) :: E			! [eV] energy of radiation particles
   real(8) :: E_spread		! [eV] spread of energy of radiation particles
!    integer :: Spectr_shape	! Spectral shape: 0 = rectangular, 1 = Gaussian, 2 = SASE
   real(8) :: t			! [fs] arrival time of the incomming particle / center of the pulse
   real(8) :: FWHM		! [fs] FWHM-duration of the pulse (ignorred for single particles)
!    integer :: Time_shape	! Temporal shape: 0 = rectangular, 1 = Gaussian, 2 = SASE
   integer :: Z			! in case of SHI, atomic number of ion in periodic table
   character(3) :: Name ! in case of an ion, its element name
   real(8) :: Zeff			! in case of SHI, user can provide its fixed charge [electron charge]
   real(8) :: Meff			! in case of SHI, user can provide its mass [a.m.u]
   integer :: KOA, Sh   ! in case of a hole, specify kind of atom and shell
end type Radiation_param


!==============================================
type output_data
   ! Total values:
   real(8) :: Nph, Ne, Nh, Np, Eph, Ee, Eh_kin, Eh_pot, Ep, Eat  ! total numbers and energies of photons, electrons, holes (kinetic and potential), positrons, atoms
   real(8) :: Nph_high, Ne_high, Nh_high, Np_high, Eph_high, Ee_high, Eh_kin_high, Eh_pot_high, Ep_high  ! numbers and energies of photons, electrons, holes (kinetic and potential), positrons with energies above cut-off
   real(8) :: MD_Ekin, MD_Epot  ! total energies in MD supercell: kinetic and potential
   ! Energy disributions (spectra):
   real(8), dimension(:), allocatable :: Spectrum_ph, Spectrum_e, Spectrum_p, Spectrum_SHI  ! energy spectra
   real(8), dimension(:,:), allocatable :: Spectrum_h ! VB spectra for each target material
   ! Velosity theta disributions:
   real(8), dimension(:), allocatable :: Vel_theta_ph, Vel_theta_e, Vel_theta_p, Vel_theta_h, Vel_theta_SHI
   ! Spectra in 1d space:
   real(8), dimension(:,:), allocatable :: Spectra_ph_X, Spectra_e_X, Spectra_p_X, Spectra_h_X, Spectra_SHI_X  ! energy spectra in space along X
   real(8), dimension(:,:), allocatable :: Spectra_ph_Y, Spectra_e_Y, Spectra_p_Y, Spectra_h_Y, Spectra_SHI_Y  ! energy spectra in space along Y
   real(8), dimension(:,:), allocatable :: Spectra_ph_Z, Spectra_e_Z, Spectra_p_Z, Spectra_h_Z, Spectra_SHI_Z  ! energy spectra in space along Z
   real(8), dimension(:,:), allocatable :: Spectra_ph_R, Spectra_e_R, Spectra_p_R, Spectra_h_R, Spectra_SHI_R  ! energy spectra in space along R
   ! Spatial distributions in 1d:
   real(8), dimension(:), allocatable :: Distr_ph_X, Distr_ph_Y, Distr_ph_Z, Distr_ph_R, Distr_ph_L, &
                                        Distr_ph_Theta, Distr_ph_Rc, Distr_ph_Thetac, Distr_ph_Phic ! photon
   real(8), dimension(:), allocatable :: Distr_e_X, Distr_e_Y, Distr_e_Z, Distr_e_R, Distr_e_L, &
                                        Distr_e_Theta, Distr_e_Rc, Distr_e_Thetac, Distr_e_Phic ! electron
   real(8), dimension(:), allocatable :: Distr_p_X, Distr_p_Y, Distr_p_Z, Distr_p_R, Distr_p_L, &
                                        Distr_p_Theta, Distr_p_Rc, Distr_p_Thetac, Distr_p_Phic ! positron
   real(8), dimension(:,:), allocatable :: Distr_h_X, Distr_h_Y, Distr_h_Z, Distr_h_R, Distr_h_L, &
                                        Distr_h_Theta, Distr_h_Rc, Distr_h_Thetac, Distr_h_Phic ! hole
   real(8), dimension(:), allocatable :: Distr_SHI_X, Distr_SHI_Y, Distr_SHI_Z, Distr_SHI_R, Distr_SHI_L, &
                                        Distr_SHI_Theta, Distr_SHI_Rc, Distr_SHI_Thetac, Distr_SHI_Phic ! SHI
   real(8), dimension(:), allocatable :: Distr_a_X, Distr_a_Y, Distr_a_Z, Distr_a_R, Distr_a_L, &
                                        Distr_a_Theta, Distr_a_Rc, Distr_a_Thetac, Distr_a_Phic ! Atom
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
   ! Spatial distributions in 3d:
   real(8), dimension(:,:,:), allocatable :: Distr_ph_XYZ, Distr_ph_RLTheta, Distr_ph_RcThcPhic ! photon
   real(8), dimension(:,:,:), allocatable :: Distr_e_XYZ, Distr_e_RLTheta, Distr_e_RcThcPhic ! electron
   real(8), dimension(:,:,:), allocatable :: Distr_p_XYZ, Distr_p_RLTheta, Distr_p_RcThcPhic ! positron
   real(8), dimension(:,:,:,:), allocatable :: Distr_h_XYZ, Distr_h_RLTheta, Distr_h_RcThcPhic ! hole
   real(8), dimension(:,:,:), allocatable :: Distr_SHI_XYZ, Distr_SHI_RLTheta, Distr_SHI_RcThcPhic ! SHI
   real(8), dimension(:,:,:), allocatable :: Distr_a_XYZ, Distr_a_RLTheta, Distr_a_RcThcPhic ! Atom
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
   ! Spatial energy distributions in 3d:
   real(8), dimension(:,:,:), allocatable :: E_Distr_ph_XYZ, E_Distr_ph_RLTheta, E_Distr_ph_RcThcPhic ! photon
   real(8), dimension(:,:,:), allocatable :: E_Distr_e_XYZ, E_Distr_e_RLTheta, E_Distr_e_RcThcPhic ! electron
   real(8), dimension(:,:,:), allocatable :: E_Distr_p_XYZ, E_Distr_p_RLTheta, E_Distr_p_RcThcPhic ! positron
   real(8), dimension(:,:,:,:), allocatable :: E_Distr_h_XYZ, E_Distr_h_RLTheta, E_Distr_h_RcThcPhic ! hole
   real(8), dimension(:,:,:), allocatable :: E_Distr_SHI_XYZ, E_Distr_SHI_RLTheta, E_Distr_SHI_RcThcPhic ! SHI
   real(8), dimension(:,:,:), allocatable :: E_Distr_a_XYZ, E_Distr_a_RLTheta, E_Distr_a_RcThcPhic ! Atom
end type output_data



!==============================================
! Parameter of the output grids etc.:
type grid_params
   ! index defined in subroutine 'read_output_grid_coord' in module "Read_numerical_parameters":
   logical :: along_axis       ! does the user want to printout the data along this axes combinations?
   logical :: log_scale(3)     ! printout spacial grids along 3 axes in log or linear scale, for each grid index
   real(8) :: gridstep(3)      ! [A] grid step along each used axis in case of linear scales
   real(8) :: gridstart(3), gridend(3)    ! [A] starting and ending points of the grid along each axis, for each index
end type grid_params

type grids_sets
   real(8), dimension(:), allocatable :: spatial_grid1  ! along the first coordinate (X, or R, or Rc)
   real(8), dimension(:), allocatable :: spatial_grid2  ! along the second coordinate (Y, ot L, or Theta_c)
   real(8), dimension(:), allocatable :: spatial_grid3  ! along the third coordinate (Z, ot Theta, or Phi_c)
end type grids_sets


!==============================================
! Paramets used for gnuplpot:
type gnu_par
   logical :: do_gnuplot
   character(10) :: gnu_extension
   character(100) :: gnu_terminal
end type gnu_par



!==============================================
! MD classical potentials:
type MD_pot ! properties of all MD potentials
   character(50) :: Name
   real(8) :: d_cut ! [A] cut off distance
   real(8) :: dd    ! [A] cut off smearing parameter
end type MD_pot

! ZBL:
type, extends (MD_pot) :: ZBL
    ! no additional parameters required
end type ZBL

! Lenard-Jones:
type, extends (MD_pot) :: LJ
   real(8) :: C12	! [eV*A^12] Lenard-Jones C12
   real(8) :: C6	! [eV*A^6] Lenard-Jones C6
end type LJ

! Power law:
type, extends (MD_pot) :: Power ! C*r^n
   real(8) :: C, n
end type Power

! Exponential law:
type, extends (MD_pot) :: Exp_law ! C*exp(-k*r)
   real(8) :: C, k
end type Exp_law

! Buckingham:
type, extends (MD_pot) :: Buck
   real(8) :: A, B, C
end type Buck

! Matsui:
type, extends (MD_pot) :: Matsui
   real(8) :: A, B, C
end type Matsui

! Truncated Coulomb (following Wolf et al.):
type, extends (MD_pot) :: Coulomb
   real(8) :: Z1, Z2
end type Coulomb

! Short-range Coulomb (softly cut, like screened):
type, extends (MD_pot) :: Coulomb_screened
   real(8) :: Z1, Z2
end type Coulomb_screened

! Long-range Coulomb (1/r, requiring special treatement according to Ewalds):
type, extends (MD_pot) :: Coulomb_Ewald
   real(8) :: Z1, Z2
end type Coulomb_Ewald

! Soft Coulomb:
type, extends (MD_pot) :: Soft_Coulomb
   real(8) :: Z1, Z2
   real(8) :: r0    ! [A] k*Z1*Z2*e^2/(r+r0)
end type Soft_Coulomb

! Morse:
type, extends (MD_pot) :: Morse
   ! https://en.wikipedia.org/wiki/Morse_potential
   real(8) :: De    ! [eV] well depth
   real(8) :: re    ! [A] equilibrium bond distance
   real(8) :: a     ! [1/A] width of the potential "alpha"
end type Morse

!---------------------
! Parameters of MD combined within an object:
type MD_pot_set
   class(MD_pot), allocatable :: Par ! parameters of the MD potential of the defined class
end type MD_pot_set

type MD_potential
   ! To account for a possibility of multiple potentials of different types
   type(MD_pot_set), dimension(:), allocatable :: Set   ! set of MD parameters
   character(3) :: El1, El2 ! names of elements which this potential is for
end type MD_potential



!==============================================
! Data type containing all numerical parameters used in the code:
type Num_par
   character(1) :: path_sep	! path separator for the environment
   ! OUTPUT FILE NAMES AND NUMBERS:
   character(200) :: input_path	! path to the folder with all input data
   character(200) :: output_path	! path to the folder with all output data
   character(200) :: FILE_parameters	! name of the file with output parameters
   integer :: FN_parameters		! file number of the file with output parameters
   character(200) :: FILE_communication	! name of the file thru which user can communicate with the program
   integer :: FN_communication	! file thru which user can communicate with the program
   integer :: MOD_TIME ! time when the communication.txt file was last modified
   character(200), dimension(:), allocatable :: FILE_totals     ! name of the file with total numbers
   integer, dimension(:), allocatable :: FN_totals              ! file number with total numbers
   character(200), dimension(:), allocatable :: FILE_totals_cutoff  ! name of the file with total numbers (above cutoff)
   integer, dimension(:), allocatable :: FN_totals_cutoff           ! file number with total numbers (above cutoff)

   ! MD files:
   character(200) :: FILE_MD_totals    ! name of the file with total numbers in MD
   integer :: FN_MD_totals             ! file number with total numbers in MD
   character(200) :: FILE_MD_average   ! name of the file with average values in MD
   integer :: FN_MD_average            ! file number with average values in MD
   character(200) :: FILE_MD_XYZ       ! name of the file with atomic coordinates from MD in XYZ
   integer :: FN_MD_XYZ                ! file number with atomic coordinates from MD in XYZ
   character(200) :: FILE_MD_V_XYZ     ! name of the file with atomic velocities from MD in XYZ
   integer :: FN_MD_V_XYZ              ! file number with atomic velocities from MD in XYZ   
   character(200) :: FILE_MCMD_info    ! name of the file with energy transfered from MC to MD
   integer :: FN_MCMD_info             ! file number with energy transfered from MC to MD
   character(200) :: FILE_MD_LAMMPS       ! name of the file with atomic coordinates from MD in LAMMPS format
   integer :: FN_MD_LAMMPS                ! file number with atomic coordinates from MD in LAMMPS format

   ! Files with spectra:
   character(200) :: FILE_spectrum_ph    ! name of the file with photonic spectrum
   integer :: FN_spectrum_ph     ! file number with photonic spectrum
   character(200) :: FILE_spectrum_e    ! name of the file with electronic spectrum
   integer :: FN_spectrum_e     ! file number with electronic spectrum
   character(200), dimension(:), allocatable :: FILE_spectrum_h    ! name of the file with holes spectrum
   integer, dimension(:), allocatable :: FN_spectrum_h     ! file number with holes spectrum
   character(200) :: FILE_spectrum_p    ! name of the file with positronic spectrum
   integer :: FN_spectrum_p     ! file number with positronic spectrum
   character(200) :: FILE_spectrum_SHI    ! name of the file with SHI spectrum
   integer :: FN_spectrum_SHI     ! file number with SHI spectrum
   ! name of the file with velosity distribution by theta
   character(200) :: FILE_vel_theta_ph, FILE_vel_theta_e, FILE_vel_theta_h, FILE_vel_theta_p, FILE_vel_theta_SHI
   integer :: FN_vel_theta_ph, FN_vel_theta_e, FN_vel_theta_h, FN_vel_theta_p, FN_vel_theta_SHI
   ! Files with spectra vs space:
   integer :: FN_spectrum_ph_X, FN_spectrum_e_X, FN_spectrum_p_X, FN_spectrum_SHI_X, FN_spectrum_h_X
   integer :: FN_spectrum_ph_Y, FN_spectrum_e_Y, FN_spectrum_p_Y, FN_spectrum_SHI_Y, FN_spectrum_h_Y
   integer :: FN_spectrum_ph_Z, FN_spectrum_e_Z, FN_spectrum_p_Z, FN_spectrum_SHI_Z, FN_spectrum_h_Z
   ! File numbers with spatial distributions:
   ! Cartesian:
   integer :: FN_car_1d_X_ph, FN_car_1d_X_e, FN_car_1d_X_p, FN_car_1d_X_SHI, FN_car_1d_X_a   ! densities along X
   integer :: FN_car_1d_Y_ph, FN_car_1d_Y_e, FN_car_1d_Y_p, FN_car_1d_Y_SHI, FN_car_1d_Y_a   ! densities along Y
   integer :: FN_car_1d_Z_ph, FN_car_1d_Z_e, FN_car_1d_Z_p, FN_car_1d_Z_SHI, FN_car_1d_Z_a   ! densities along Z
   integer :: FN_car_1d_X_E_ph, FN_car_1d_X_E_e, FN_car_1d_X_E_p, FN_car_1d_X_E_SHI, FN_car_1d_X_E_a   ! energy densities along X
   integer :: FN_car_1d_Y_E_ph, FN_car_1d_Y_E_e, FN_car_1d_Y_E_p, FN_car_1d_Y_E_SHI, FN_car_1d_Y_E_a   ! energy densities along Y
   integer :: FN_car_1d_Z_E_ph, FN_car_1d_Z_E_e, FN_car_1d_Z_E_p, FN_car_1d_Z_E_SHI, FN_car_1d_Z_E_a   ! energy densities along Z
   integer, dimension(:), allocatable :: FN_car_1d_X_h, FN_car_1d_X_E_h   ! densities and doses of holes in all shells of all elements along X
   integer, dimension(:), allocatable :: FN_car_1d_Y_h, FN_car_1d_Y_E_h   ! densities and doses of holes in all shells of all elements along Y
   integer, dimension(:), allocatable :: FN_car_1d_Z_h, FN_car_1d_Z_E_h   ! densities and doses of holes in all shells of all elements along Z
   ! Cyllindrical:
    !1d   
   integer :: FN_cyl_1d_R_ph, FN_cyl_1d_R_e, FN_cyl_1d_R_p, FN_cyl_1d_R_SHI, FN_cyl_1d_R_a   ! densities
   integer :: FN_cyl_1d_R_E_ph, FN_cyl_1d_R_E_e, FN_cyl_1d_R_E_p, FN_cyl_1d_R_E_SHI, FN_cyl_1d_R_E_a   ! energy densities
   integer, dimension(:), allocatable :: FN_cyl_1d_R_h, FN_cyl_1d_R_E_h   ! densities and energy densities of holes in all shells of all elements
    !2d
   integer :: FN_cyl_2d_RL_ph, FN_cyl_2d_RL_e, FN_cyl_2d_RL_p, FN_cyl_2d_RL_SHI, FN_cyl_2d_RL_a   ! densities
   integer :: FN_cyl_2d_RL_E_ph, FN_cyl_2d_RL_E_e, FN_cyl_2d_RL_E_p, FN_cyl_2d_RL_E_SHI, FN_cyl_2d_RL_E_a   ! energy densities
   integer, dimension(:), allocatable :: FN_cyl_2d_RL_h, FN_cyl_2d_RL_E_h   ! densities and energy densities of holes in all shells of all elements
   ! OUTPUT PRINTOUT:
   type(gnu_par) :: gnupl       ! parameters for gnuplotting
   logical :: printout_DOS      ! user defines to printout analyzed DOS and related parameters or not
   logical :: printout_MFPs     ! user defines to printout analytical particles mean free paths or not
   logical :: printout_ranges   ! user defines to printout analytical particles ranges or not
   logical :: print_MD_R_xyz, print_MD_V_xyz    ! printout all atoms coordinates, velocities
   logical :: print_MD_LAMMPS   ! printout all atoms coordinates, velocities as input file for LAMMPS
   character(20) :: LAMMPS_UNITS       ! Units according to LAMMPS format: https://docs.lammps.org/units.html
   ! Units are only used to write LAMMPS input file (as one of TREKIS-4 outpit files), but TREKIS itself use its own units

   type(grid_params), dimension(19) :: grid_par   ! all the parameters of printout spatial grids
   type(grids_sets), dimension(19) :: grids          ! all spatial grids for the output data
   type(grid_params), dimension(19) :: Spectr_grid_par   ! all the parameters of printout Energy vs spatial grids
   type(grids_sets), dimension(19) :: Spectr_grid  ! Space grid in 1d for spectra calculations
   type(grid_params) :: NRG_grid_par   ! all the parameters of printout energy grid
   real(8), dimension(:), allocatable :: NRG_grid  ! energy grid for the output data (spectra)
   real(8), dimension(:), allocatable :: NRG_grid_VB  ! energy grid for the VB output data (VB spectra)
   type(grid_params) :: vel_theta_grid_par   ! all the parameters of printout velosity theta distribution
   real(8), dimension(:), allocatable :: vel_theta_grid  ! particles velosity distribution by theta: Vz/V
   ! NUMERICS OF THE TARGET PARAMETERS:
   integer :: N_sh_tot  ! total number of core shells in the target material
   ! NUMERICAL PARAMETERS:
   integer :: NOMP	! number of nodes for OpenMP
   integer :: NMC	! number of MC iterations
   real(8) :: dt_MD	! [fs] time step for MD
   real(8), dimension(:), allocatable :: dt_MD_grid    ! user-provided changes in the time grid
   real(8), dimension(:), allocatable :: dt_MD_reset_grid  ! user-provided grid for changing MD time-step
   integer :: i_dt          ! counter for MD timegrid
   real(8) :: dt_printout	! [fs] how often to print out the data into files
   integer :: i_dt_printout ! counter for timegrid to ptrintou data
   real(8), dimension(:), allocatable :: dt_out_grid    ! user-provided time grid
   real(8), dimension(:), allocatable :: dt_reset_grid  ! user-provided grid for changing time-step
   real(8) :: t_start   ! [fs] starting time of the simulation
   real(8) :: t_total	! [fs] when to stop simulation
   logical :: print_each_step   ! printout time each MD timestep or not?
   ! SIMULATION BOX PARAMETERS:
   logical :: reset_from_MD(3)  ! if the simulation box boundaries should be taken MD supercell
   real(8) :: box_start_x, box_end_x    ! [A] coordinates of the left and right ends of the simulation box along X
   real(8) :: box_start_y, box_end_y    ! [A] coordinates of the left and right ends of the simulation box along Y
   real(8) :: box_start_z, box_end_z    ! [A] coordinates of the left and right ends of the simulation box along Z
   integer :: periodic_x, periodic_y, periodic_z    ! Bondary conditions along x,y,z: 0=free, 1=periodic
   ! MODULES AND MODELS:
   logical :: DO_MC		! Activate MC module or not
   logical :: DO_MD		! Activate MD module or not
   logical :: DO_TTM	! Activate TTM module or not
   integer :: MC_vs_MD 	! MC or MD target model: 0=MC, 1=MD
   logical :: recalculate_MFPs  ! Force program to recalculate cross sections and MFPs even if there are files containing precalculated ones
   ! MODELS FOR ELECTRONS:
   integer :: El_forces		! Include forces and fields among electrons: 0=exclude, 1=include as MD, 2=include as mean field
   integer :: El_inelast		! inelastic scattering: 0=excluded, 1=CDF, 2=RBEB
   integer :: El_elastic		! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   integer :: El_Brems		! Bremsstrahlung: 0=excluded, 1= BHW
   integer :: El_Cheren		! Cherenkov radiation: 0=excluded, 1= ...
   real(8) :: El_Cutoff		! [eV] Cut-off energy (electrons with lower energies are excluded from calculation)
   ! MODELS FOR POSITRONS:
   integer :: Pos_inelast		! inelastic scattering: 0=excluded, 1=CDF, 2=RBEB
   integer :: Pos_elastic		! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   integer :: Pos_Brems		! Bremsstrahlung: 0=excluded, 1= BHW
   integer :: Pos_annih		! Annihilation: 0=excluded, 1= Heitler
   real(8) :: Pos_Cutoff		! [eV] Cut-off energy (electrons with lower energies are excluded from calculation)
   ! MODELS FOR PHOTONS:
   integer :: Ph_absorb	! Photoabsorption CSs: 0=excluded, 1=CDF, 2=EPDL97
   integer :: Ph_Compton	! Compton effect: 0=excluded, 1= ...
   integer :: Ph_Thomson	! Thomson / Rayleigh scattering: 0=excluded, 1= ...
   integer :: Ph_Pairs		! Electron-positron pair creation: 0=excluded, 1=included, 2=included with Landau-Pomeranchuk-Migdal
   integer :: Ph_Nucl		! Photonuclear Physics: 0=excluded, 1=included
   real(8) :: Ph_Cutoff		! [eV] Cut-off energy (electrons with lower energies are excluded from calculation)
   real(8) :: Ph_att_eff     ! [A] Effective photon attenuation length (<0 means do not use effective one, use real one from EPDL)
   ! MODELS FOR SHI:
   integer :: SHI_inelast		! inelastic scattering: 0=excluded, 1:3=CDF Delta, 4=CDF Ritchie
   integer :: SHI_ch_st		! Charge state: 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande, 4=fixed Zeff, 5=charge exchange
   integer :: SHI_ch_shape	! SHI charge shape: 0=point-like charge; 1=Brandt-Kitagawa ion
   real(8) :: SHI_Cutoff		! [eV] Cut-off energy (electrons with lower energies are excluded from calculation)
   ! MODEL PARAMETERS FOR CDF:
   integer :: CDF_model	! CDF model: 0 = Drude / Lindhard CDF, 1=Mermin CDF, 2=Full conserving CDF
   integer :: CDF_dispers	! target dispersion relation: 1=free electron, 2=plasmon-pole, 3=Ritchie
   integer :: CDF_m_eff	! effective mass [in me] (0=effective mass from DOS of VB; -1=free-electron)
   integer :: CDF_plasmon	! Include plasmon integration limit (0=no, 1=yes)
   real(8) :: CDF_Eeq_factor    ! coeff. where to switch from Ritchie to Delta CDF: E = CDF_Eeq_factor * Wmin (INELASTIC)
   real(8) :: CDF_Eeq_elast     ! coeff. where to switch from Ritchie to Delta CDF: E = CDF_Eeq_elast * Wmin (ELASTIC)
   integer :: CDF_elast_Zeff    ! use Zeff or Z=1
   integer :: CDF_int_n_inelastE ! effective number of grid points for inelastic cross section integration over energy
   integer :: CDF_int_n_inelastQ ! effective number of grid points for inelastic cross section integration over momentum
   integer :: CDF_int_n_elastE ! effective number of grid points for elastic cross section integration over energy
   integer :: CDF_int_n_elastQ ! effective number of grid points for elastic cross section integration over momentum
   ! MODELS FOR HOLES:
   integer :: H_Auger		! Auger decays: 0=excluded, 1=provided in cdf-file, 2=EADL
   integer :: H_Radiat		! Radiative decays: 0=excluded, 1=EADL
   real(8) :: H_m_eff		! [m_e] effective mass of a VB-hole
   integer :: H_inelast		! inelastic scattering: 0=excluded, 1=CDF, 2=RBEB
   integer :: H_elast		! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   real(8) :: H_Cutoff		! [eV] Cut-off energy (electrons with lower energies are excluded from calculation)
   ! MD NUMERICAL PARAMETERS:
   integer :: MD_integrator     ! index of MD integrator (currently only Verlet)
   logical :: recenter(3)       ! recenter MD supercell along 3 axes or not?
   logical :: print_MC_MD_energy    ! printout energies passed from MC to MD
   logical :: use_thermostat    ! use thermostat for the MD supercell or not?
   integer :: thermostat_inx    ! index of thermostat to use: 0=Berendsen; 1=NOT READY
   logical :: damp_pressure     ! use damping of pressure wave or not?
   logical :: do_quenching      ! use quenching or not?
   logical :: do_cohesive       ! calculate cohesive energy instead of full MD run
   real(8) :: t_quench_start, dt_quench ! [fs] when to start quenching, and how often to use
   real(8) :: t_quench_run      ! [fs] running time to trace when to quench
   integer, dimension(:), allocatable :: Neighbors_Num      ! number of nearest neighbors for each atom
   integer, dimension(:,:), allocatable :: Neighbors_list   ! list of nearest neighbors for each atom
   real(8), dimension(:,:), allocatable :: Neighbors_Rij    ! distance between atoms in the nearest neighbors list [A]
   real(8), dimension(:,:,:), allocatable :: Neighbors_Xij  ! distance between atoms along 3 axes in the nearest neighbors list [A]
   real(8), dimension(:), allocatable :: Ewalds_kfac        ! exp(-k*b)/k^2 for Ewalds summation
end type Num_par


!==============================================
! Types used in MC and MD:
type :: Particle    ! basic class for all particles
   logical :: active    ! switch: is particle participating in calculation or not
   real(8) :: Ekin      ! [eV] kinetic energy
   real(8) :: Mass      ! [kg] (effective) mass of this particle
   real(8) :: t0    ! time of the last event (scattering, decay, etc.)
   real(8) :: ti    ! time of the next event (scattering, decay, etc.)
   real(8) :: t_sc  ! time of the next scattering event (excluding a possibility of a border crossing)
   integer :: generation    ! which generation it belongs: 0=incomming, 1=created by incomming ones, 2=secondary generated... etc.
   integer :: in_target     ! in which target this particle is now
   ! coordinates and velosities:
   real(8), dimension(3) :: R	! [A] coordinates (x,y,z)
   real(8), dimension(3) :: S	! [a.u.] relative coordinates (sx, sy, sz)
   real(8), dimension(3) :: V	! [A/fs] velosities (Vx, Vy, Vz)
   real(8), dimension(3) :: SV	! [a.u.] relative velosities (SVx, SVy, SVz)
   ! on the last time-step:
   real(8), dimension(3) :: R0	! [A] coordinates
   real(8), dimension(3) :: S0	! [a.u.] relative coordinates
   real(8), dimension(3) :: V0	! [A/fs] velosities
   real(8), dimension(3) :: SV0	! [a.u.] relative velosities
end type Particle

type, EXTENDS (Particle) :: Photon ! photon as an object, contains no additional info
end type Photon

type, EXTENDS (Particle) :: Electron ! electron as an object, contains the same into as photon + more:
   real(8) :: U     ! [eV] potential energy
   real(8), dimension(3) :: A       ! [A^2/fs] accelerations
   real(8), dimension(3) :: Force   ! [eV/fs] forces
end type Electron

type, EXTENDS (Electron) :: Positron    ! positron as an object
end type Positron

type, EXTENDS (Electron) :: Hole    ! hole as an object
   integer :: KOA	! kind of atom
   integer :: Sh	! shell in which this hole sits
   logical :: valent    ! is it a valent hole (true) of core hole (false)
end type Hole

type, EXTENDS (Electron) :: Atom    ! atom for MD as an object
   integer :: KOA		! kind of atom in the compound (according to Atom_kind type above)
   integer :: Z		! atomic number in the periodic table
   character(3) :: Name	! abbreviation of the atom according to periodic table
end type Atom

type, EXTENDS (Atom) :: SHI     ! swift heavy ion as an object
   real(8) :: Zeff	! effective charge [electron charge]
   real(8) :: Meff	! user-defined mass [amu]
end type SHI

type, EXTENDS (Atom) :: MacroAtom    ! Macroparticle for MD/Coarse-grained/PIC etc. simulations
   ! Macroparticles could be molecules, clusters, macroparticles in PIC, etc., any coarse-grained lump of atoms
   real(8) :: Radius    ! [A] radius of the macroparticle
end type MacroAtom

!==============================================
! Parameters for exchange energy between MC and MD:
type :: MCMD_grid
   real(8), dimension(:), allocatable :: grid    ! grids along 3 axes
end type MCMD_grid

! Parameters of the MD supercell / simulation box:
type :: MD_supcell
   real(8) :: x_start, x_end    ! [A] coordinates of the start and end of the box along X
   real(8) :: y_start, y_end    ! [A] coordinates of the start and end of the box along Y
   real(8) :: z_start, z_end    ! [A] coordinates of the start and end of the box along Z
   real(8), dimension(3) :: vect        ! [A] orthogonal supercell sizes along X, Y, Z
   real(8) :: V                 ! [A^3] volume of the supercell
   integer, dimension(3) :: boundary    ! index for boundary conditions: 0=free, 1=periodic
   real(8), dimension(3,3) :: SC_vec    ! [A] matrix of supercell vectors (must be consistent with coordinates): UNUSED
   ! Derived from MD and evolving properties of the supercell:
   real(8) :: Ekin_tot, Epot_tot    ! average kinetic and potential energy of all atoms [eV/atom]
   real(8) :: T_kin ! average kinetic temperature of all atoms [K]
   real(8), dimension(3,3) :: Pressure, P_kin, P_pot ! average pressure tensor [GPa], its kinetic and potential contributions
   ! Thermostat parameters:
   real(8) :: thermo_dx_start, thermo_dx_end   ! thickness of layer where thermostat is used along X
   real(8) :: thermo_dy_start, thermo_dy_end   ! thickness of layer where thermostat is used along Y
   real(8) :: thermo_dz_start, thermo_dz_end   ! thickness of layer where thermostat is used along Z
   real(8) :: therm_tau ! [fs] characteristic time of cooling
   real(8) :: Bath_T    ! [K] bath temperature
   ! Pressure wave dampener:
   real(8) :: pressdamp_dx_start, pressdamp_dx_end   ! thickness of layer where pressure waves are damped along X
   real(8) :: pressdamp_dy_start, pressdamp_dy_end   ! thickness of layer where pressure waves are damped along Y
   real(8) :: pressdamp_dz_start, pressdamp_dz_end   ! thickness of layer where pressure waves are damped along Z
   real(8) :: press_tau ! [fs] characteristic time of pressure dampening
   ! Parameters of connection between MC and MD:
   integer :: coord_type  ! type of coordinates system grid: 1=cartesian; 2=spherical; 3=cylindrical
   integer :: coord_dim   ! dimensionality of the grid for transfer of info between MC and MD
   integer :: coord_axis(3)  ! index of the axis: 1=X or R or R; 2=Y or Theta or L; 3=Z or Phi or Phi (for Cart, Spher, or Cyl)
   real(8) :: step(3)  ! grid step [A] or [deg]
   type(MCMD_grid), dimension(3) :: MC2MD_dimen ! grids for info exchange between MC and MD
   real(8), dimension(:,:,:), allocatable :: E_e_at_from_MC  ! energy elastically transferred from electrons to atoms from within MC module
   real(8), dimension(:,:,:), allocatable :: E_h_at_from_MC  ! energy elastically transferred from holes to atoms from within MC module
   real(8), dimension(:,:,:), allocatable :: E_p_at_from_MC  ! energy elastically transferred from positrons to atoms from within MC module
   real(8), dimension(:,:,:), allocatable :: E_e_from_MC   ! energy of electrons below cut-off from within MC module
   real(8), dimension(:,:,:), allocatable :: E_h_from_MC   ! energy of holes below cut-off from within MC module
end type MD_supcell


!==============================================
! Arrays for all types of MC particles collected together:
type :: MC_arrays
   ! Number of particles currently active:
   integer :: N_ph    ! number of active photons
   integer :: N_e     ! number of active electrons
   integer :: N_p     ! number of active positrons
   integer :: N_h     ! number of active holes
   integer :: N_SHI   ! number of active SHIs
   integer :: N_at_nrg  ! number of elastic scattering events transfering energy to atoms
   ! Arrays for all MC particles:
   type(Photon), dimension(:), allocatable :: MC_Photons        ! all photons as objects
   type(Electron), dimension(:), allocatable :: MC_Electrons    ! all electrons as objects
   type(Positron), dimension(:), allocatable :: MC_Positrons    ! all positrons as objects
   type(Hole), dimension(:), allocatable :: MC_Holes    ! all holes as objects
   type(SHI), dimension(:), allocatable :: MC_SHIs      ! all SHIs as objects
   type(Atom), dimension(:), allocatable :: MC_Atoms_events     ! all elastic energy transfer events as objects
end type MC_arrays


!==============================================
! Error handling as an object:
type :: Error_handling
   logical Err		! indicates that an error occured
   integer Err_Num	! assign a number to an error
   integer File_Num		! number of the file with error log
   character(200) File_Name	! name of the file with the error log
   character(200) Err_descript	! describes more details about the error
end type Error_handling
!==============================================


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 contains
! How to write a log about an error:
subroutine Save_error_details(Err_name, Err_num, Err_data)
   type(Error_handling) :: Err_name
   integer, intent(in) :: Err_num
   character(*), intent(in) :: Err_data
   integer FN
   FN = Err_name%File_Num
   Err_name%Err = .true.
   Err_name%Err_Num = Err_num
   Err_name%Err_descript = Err_data
   write(FN, '(a,i2,1x,a)') 'Error #', Err_name%Err_Num, trim(adjustl(Err_name%Err_descript))
   ! Descriptors of errors:
   ! Err_num = 1 means error in opening a file
   ! Err_num = 2 means error in reading a file
   ! Err_num = 3 means error in a database
   ! Err_num = 4 means error in format somewhere with the Periodic Table
   ! Err_num = 5 means error in a subroutine call
   ! Err_num = 6 means error in format of input
end subroutine Save_error_details


pure subroutine allocate_cross_section(CS, Ngrid, Nshl)
   type(Cross_section), intent(inout) :: CS	! cross section variable
   integer, intent(in) :: Ngrid, Nshl		! number of grid points, number of shell in this element
   if (.not.allocated(CS%E)) allocate(CS%E(Ngrid))	! [eV] energy grid
   if (.not.allocated(CS%Total)) allocate(CS%Total(Ngrid))	! [A^2] total cross section
   if (.not.allocated(CS%Per_shell)) allocate(CS%Per_shell(Nshl,Ngrid))	! number of shell,  [A^2] cross section for chosen shell
end subroutine allocate_cross_section



pure subroutine set_default_particle_array(Prtcl, typ, siz)
   class(Particle), dimension(:), allocatable, intent(inout) :: Prtcl	! undefined particle as an object
   character(*), intent(in) :: typ	! type of the particle
   integer, intent(in) :: siz	! size of the array
   if (.not.allocated(Prtcl)) then	! it's the first time, define the array
      select case (typ)	! which particle array that should become
      case ('PHOTON', 'Photon', 'photon', 'ph')
          allocate(Photon::Prtcl(siz)) ! make it an array of photon, size SYZ
      case ('ELECTRON', 'Electron', 'electron', 'e')
         allocate(Electron::Prtcl(siz)) ! make it an array of Electrons, size SYZ
      case ('POSITRON', 'Positron', 'positron', 'p')
         allocate(Positron::Prtcl(siz)) ! make it an array of Positrons, size SYZ
      case ('HOLE', 'Hole', 'hole', 'h')
         allocate(Hole::Prtcl(siz)) ! make it an array of Holes, size SYZ
      case ('ATOM', 'Atom', 'atom', 'a')
         allocate(Atom::Prtcl(siz)) ! make it an array of Atoms, size SYZ
      case ('SHI', 'ION', 'Ion', 'ion', 'shi')
         allocate(SHI::Prtcl(siz)) ! make it an array of SHIs, size SYZ
      end select
   endif
end subroutine set_default_particle_array


pure subroutine make_new_particle(Prtcl, Ekin, Mass, t0, ti, t_sc, generation, in_target, R, S, V, SV, R0, S0, V0, SV0, &
                                    Force, KOA, Sh, valent, Z, Name, Zeff, Meff)
   class(Particle), intent(inout) :: Prtcl	! undefined particle as an object
   real(8), intent(in), optional :: Ekin      ! [eV] kinetic energy
   real(8), intent(in), optional :: Mass     ! [kg] (effective) mass of this particle
   real(8), intent(in), optional :: t0    ! time of the last event (scattering, decay, etc.)
   real(8), intent(in), optional :: ti     ! time of the next event (scattering, decay, etc.)
   real(8), intent(in), optional :: t_sc  ! time of the next scattering event (excluding a possibility of a border crossing)
   integer, intent(in), optional :: generation    ! which generation it belongs: 0=incomming, 1=created by incomming ones, 2=secondary generated... etc.
   integer, intent(in), optional :: in_target      ! in which target this particle is now
   real(8), dimension(3), intent(in), optional :: R		! [A] coordinates (x,y,z)
   real(8), dimension(3), intent(in), optional :: S		! [a.u.] relative coordinates (sx, sy, sz)
   real(8), dimension(3), intent(in), optional :: V		! [A/fs] velosities (Vx, Vy, Vz)
   real(8), dimension(3), intent(in), optional :: SV	! [a.u.] relative velosities (SVx, SVy, SVz)
   real(8), dimension(3), intent(in), optional :: R0	! [A] coordinates
   real(8), dimension(3), intent(in), optional :: S0	! [a.u.] relative coordinates
   real(8), dimension(3), intent(in), optional :: V0	! [A/fs] velosities
   real(8), dimension(3), intent(in), optional :: SV0	! [a.u.] relative velosities
   real(8), dimension(3), intent(in), optional :: Force	! [eV/fs] forces
   integer, intent(in), optional :: KOA		! kind of atom in the compound (according to Atom_kind type above)
   integer, intent(in), optional :: Sh	! shell in which this hole sits
   logical, intent(in), optional :: valent    ! is it a valent hole (true) of core hole (false)
   integer, intent(in), optional :: Z		! atomic number in the periodic table
   character(3), intent(in), optional :: Name	! abbreviation of the atom according to periodic table
   real(8), intent(in), optional :: Zeff	! effective charge [electron charge]
   real(8), intent(in), optional :: Meff	! user-defined mass [amu]
   !--------------------------------------------------
   ! To start with, make a default one:
   call set_default_particle(Prtcl) ! below
   ! Activate it:
   Prtcl%active = .true.
   if (present(Ekin)) Prtcl%Ekin = Ekin
   if (present(Mass)) Prtcl%Mass = Mass
   if (present(t0)) Prtcl%t0 = t0
   if (present(ti)) Prtcl%ti = ti
   if (present(t_sc)) Prtcl%t_sc = t_sc
   if (present(generation)) Prtcl%generation = generation
   if (present(in_target)) Prtcl%in_target = in_target
   if (present(R)) Prtcl%R = R
   if (present(S)) Prtcl%S = S
   if (present(V)) Prtcl%V = V
   if (present(SV)) Prtcl%SV = SV
   if (present(R0)) then
      Prtcl%R0 = R0
   else
      Prtcl%R0(:) = Prtcl%R(:)
   endif
   if (present(S0)) Prtcl%S0 = S0
   if (present(V0)) then
      Prtcl%V0 = V0
   else
      Prtcl%V0(:) = Prtcl%V(:)
   endif
   if (present(SV0)) Prtcl%SV0 = SV0
   
   select type (Prtcl)	! which particle array that should become
      type is (Photon)
      type is (Electron)
         if (present(Force)) Prtcl%Force = Force
      type is (Positron)
         if (present(Force)) Prtcl%Force = Force
      type is (Hole)
         if (present(Force)) Prtcl%Force = Force
         if (present(KOA)) Prtcl%KOA = KOA
         if (present(Sh)) Prtcl%Sh = Sh
         if (present(valent)) Prtcl%valent = valent
      type is (Atom)
         if (present(Force)) Prtcl%Force = Force
         if (present(KOA)) Prtcl%KOA = KOA
         if (present(Z)) Prtcl%Z = Z
         if (present(Name)) Prtcl%Name = Name
      type is (SHI)
         if (present(Force)) Prtcl%Force = Force
         if (present(KOA)) Prtcl%KOA = KOA
         if (present(Z)) Prtcl%Z = Z
         if (present(Name)) Prtcl%Name = Name
         if (present(Zeff)) Prtcl%Zeff = Zeff
         if (present(Meff)) Prtcl%Meff = Meff
   end select
end subroutine make_new_particle


pure subroutine set_default_particle(Prtcl)
   class(Particle), intent(inout) :: Prtcl	! undefined particle as an object
   ! Default values to start with:
   Prtcl%active = .false.   ! by default, a particle is excluded from simulations
   Prtcl%generation = -1    ! has not been generated yet
   Prtcl%in_target = 1      ! by default, let it be in the first target
   Prtcl%Ekin = 0.0d0     ! [eV] kinetic energy
   Prtcl%t0 = 0.0d0         ! [fs] starting time
   Prtcl%ti = 1.0d22        ! [fs] time of next event (scattering, decay, etc.)
   Prtcl%t_sc = 1.0d23        ! [fs] time of next event (scattering, decay, etc.)
   Prtcl%R(:) = 0.0d0       ! [A] coordinates (x,y,z)
   Prtcl%S(:) = 0.0d0       ! [a.u.] relative coordinates (sx, sy, sz)
   Prtcl%V(:) = 0.0d0       ! [A/fs] velosities (Vx, Vy, Vz)
   Prtcl%SV(:) = 0.0d0     ! [a.u.] relative velosities (SVx, SVy, SVz)
   Prtcl%R0(:) = Prtcl%R(:)     ! [A] coordinates
   Prtcl%S0(:) = 0.0d0     ! [a.u.] relative coordinates
   Prtcl%V0(:) = Prtcl%V(:)     ! [A/fs] velosities
   Prtcl%SV0(:)= 0.0d0    ! [a.u.] relative velosities
   
   select type (Prtcl)	! which particle array that should become
      type is (Photon)
         Prtcl%Mass = 0.0d0     ! [kg] photon mass(less)
      type is (Electron)
         Prtcl%Mass = g_me     ! [kg] free electron rest mass
         Prtcl%A(:) = 0.0d0        ! [A^2/fs] accelerations
         Prtcl%Force(:) = 0.0d0  ! [eV/fs] forces
      type is (Positron)
         Prtcl%Mass = g_me     ! [kg] free electron rest mass
         Prtcl%A(:) = 0.0d0        ! [A^2/fs] accelerations
         Prtcl%Force(:) = 0.0d0  ! [eV/fs] forces
      type is (Hole)
         Prtcl%Mass = g_me     ! [kg] free electron rest mass
         Prtcl%A(:) = 0.0d0        ! [A^2/fs] accelerations
         Prtcl%Force(:) = 0.0d0  ! [eV/fs] forces
         Prtcl%KOA = 0            ! kind of atom
         Prtcl%Sh = 0               ! shell in which this hole sits
         Prtcl%valent = .false. ! valent of core hole
      type is (Atom)
         Prtcl%Mass = g_Mp     ! [kg] proton rest mass (by default, before set by the user)
         Prtcl%Z = 1                 ! proton by default
         Prtcl%KOA = 0            ! kind of atom
         Prtcl%Name = 'H'         ! abbreviation of the atom according to periodic table
         Prtcl%A(:) = 0.0d0        ! [A^2/fs] accelerations
         Prtcl%Force(:) = 0.0d0  ! [eV/fs] forces
      type is (SHI)
         Prtcl%Mass = g_Mp     ! [kg] proton rest mass (by default, before set by the user)
         Prtcl%Z = 1                  ! proton by default
         Prtcl%KOA = 0            ! kind of atom
         Prtcl%Name = 'H'          ! abbreviation of the atom according to periodic table
         Prtcl%Zeff = 1              ! effective charge [electron charge]
         Prtcl%Meff = 1             ! user-defined mass [a.m.u.]
         Prtcl%A(:) = 0.0d0        ! [A^2/fs] accelerations
         Prtcl%Force(:) = 0.0d0  ! [eV/fs] forces
   end select
end subroutine set_default_particle


END MODULE Objects
