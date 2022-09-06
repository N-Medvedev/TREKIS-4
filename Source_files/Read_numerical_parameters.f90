! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018-2020
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to read input files:

MODULE Read_numerical_parameters
use Universal_constants
use Objects
use Dealing_with_files, only: read_file, close_file, Count_lines_in_file, Count_columns_in_file
use Little_subroutines, only: create_grid


implicit none

! this is a function that read a grid from a user provided file:
interface read_grid_from_file
   module procedure read_grid_from_file_1d
   module procedure read_grid_from_file_2d
end interface read_grid_from_file

! private :: ! hides items not listed on public statement 
public :: read_grid_from_file



! All paths to input data and databases are collected here within this module:
character(50) :: m_input_folder

! All input folders / directories:
parameter(m_input_folder = 'INPUT_DATA')			! folder with all input data (other folders and files)


 contains

!==============================
! Read the file with all numerical parameters:
subroutine read_num_pars(FN, File_name, numpar, Err)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !----------------------
   integer :: Reason, count_lines, i, N_rp, temp
   character(200) :: Error_descript, temp_ch
   logical :: read_well
   
   count_lines = 0	! to start counting lines in the file

   !==================================================
   ! Skip the first line: it is only an indicator. May be used for comments.
   read(FN,*,IOSTAT=Reason)     ! ::: NUMERICAL PARAMETERS :::
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! number of MC iterations
   read(FN,*,IOSTAT=Reason) numpar%NMC
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! number of threads for parallel calculations with OpenMP (1 if nonparrelelized)
   read(FN,*,IOSTAT=Reason) numpar%NOMP
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! [fs] Time-step for MD
   read(FN,*,IOSTAT=Reason) temp_ch  ! check whether there is time grid provided
   !read(FN,*,IOSTAT=Reason) numpar%dt_MD
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   else
      call set_MD_step_grid(temp_ch, numpar, read_well, Error_descript)  ! below
      if (.not. read_well) then
         call Save_error_details(Err, 2, Error_descript)    ! module "Objects"
         goto 9998
      endif
   endif
   
   
   ! [fs] How often to print out the data into files
   read(FN,*,IOSTAT=Reason) temp_ch  ! check whether there is time grid provided
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      write(Error_descript,'(a)') trim(adjustl(Error_descript))//'. Setting printout timestep equal to simulation step.'
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      numpar%dt_printout = numpar%dt_MD ! printout data every timestep by default
   else
      call set_time_grid(temp_ch, numpar)   ! below
   endif
   
   ! [fs] when to stop simulation
   read(FN,*,IOSTAT=Reason) numpar%t_start
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif

   ! [fs] when to stop simulation
   read(FN,*,IOSTAT=Reason) numpar%t_total
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Activate MC module or not:
   read(FN,*,IOSTAT=Reason) temp
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   if (temp /= 0) then	! perform MC calculations:
      numpar%DO_MC = .true.
   else	! do NOT perform MC calculations:
      numpar%DO_MC = .false.
   endif 
   
   ! Activate MD module or not:
   read(FN,*,IOSTAT=Reason) temp
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   if (temp /= 0) then	! perform MD calculations:
      numpar%DO_MD = .true.
   else	! do NOT perform MD calculations:
      numpar%DO_MD = .false.
   endif 
   
   ! Activate TTM module or not:
   read(FN,*,IOSTAT=Reason) temp
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   if (temp /= 0) then	! perform MD calculations:
      numpar%DO_TTM = .true.
   else	! do NOT perform MD calculations:
      numpar%DO_TTM = .false.
   endif
   
   ! Recalculate cross sections and MFPs (T=true; F=false):
   read(FN,*,IOSTAT=Reason) numpar%recalculate_MFPs
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif


   ! [A] coordinates of the left and right ends of the simulation box along X, and whether reset it via MD:
   read(FN,*,IOSTAT=Reason) numpar%box_start_x, numpar%box_end_x, numpar%reset_from_MD(1)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif

   ! [A] coordinates of the left and right ends of the simulation box along Y, and whether reset it via MD:
   read(FN,*,IOSTAT=Reason) numpar%box_start_y, numpar%box_end_y, numpar%reset_from_MD(2)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! [A] coordinates of the left and right ends of the simulation box along Z, and whether reset it via MD:
   read(FN,*,IOSTAT=Reason) numpar%box_start_z, numpar%box_end_z, numpar%reset_from_MD(3)
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif

   read(FN,*,IOSTAT=Reason) numpar%periodic_x, numpar%periodic_y, numpar%periodic_z    ! Bondary conditions along x,y,z: 0=free, 1=periodic
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   !==================================================
   ! Skip line
   read(FN,*,IOSTAT=Reason) !::: MODELS FOR ELECTRONS :::
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! MC or MD target model: 0=MC, 1=MD
   read(FN,*,IOSTAT=Reason) numpar%MC_vs_MD
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Include forces and fields among electrons: 0=exclude, 1=include
   read(FN,*,IOSTAT=Reason) numpar%El_forces
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! inelastic scattering: 0=excluded, 1=relativ.CDF, 2=RBEB, 3=delta, 4=nonrelativ.CDF (DO NOT USE!), 5=SPdelta
   read(FN,*,IOSTAT=Reason) numpar%El_inelast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! elastic scattering: 0=excluded, 1=CDF, 2=Mott, 3=DSF
   read(FN,*,IOSTAT=Reason) numpar%El_elastic
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Bremsstrahlung: 0=excluded, 1= ...
   read(FN,*,IOSTAT=Reason) numpar%El_Brems
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Cherenkov radiation: 0=excluded, 1= ...
   read(FN,*,IOSTAT=Reason) numpar%El_Cheren
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! [eV] Cut-off energy (electrons with lower energies are excluded from calculation):
   read(FN,*,IOSTAT=Reason) numpar%El_Cutoff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   
   !==================================================
   ! Skip line
   read(FN,*,IOSTAT=Reason) !::: MODELS FOR PHOTONS :::
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Photoabsorption CSs: 0=excluded, 1=CDF, 2=EPDL97:
   read(FN,*,IOSTAT=Reason) numpar%Ph_absorb
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Compton effect: 0=excluded, 1= ...
   read(FN,*,IOSTAT=Reason) numpar%Ph_Compton
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Thomson / Rayleigh scattering: 0=excluded, 1= ...
   read(FN,*,IOSTAT=Reason) numpar%Ph_Thomson
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Electron-positron pair creation: 0=excluded, 1=included, 2=included with 6) Landau-Pomeranchuk-Migdal suppression effect:
   read(FN,*,IOSTAT=Reason) numpar%Ph_Pairs
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Photonuclear Physics: 0=excluded, 1=included:
   read(FN,*,IOSTAT=Reason) numpar%Ph_Nucl
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! [eV] Cut-off energy (photons with lower energies are excluded from calculation):
   read(FN,*,IOSTAT=Reason) numpar%Ph_Cutoff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! [A] Effective photon attenuation length (<0 means do not use effective one, use real one from EPDL):
   read(FN,*,IOSTAT=Reason) numpar%Ph_att_eff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   !==================================================
   ! Skip line:
   read(FN,*,IOSTAT=Reason) !::: MODELS FOR SHI :::
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif

   ! SHI inelastic scattering: 0=excluded, 1:3=delta, 4=nonrelativ.CDF (DO NOT USE!), 5=SPdelta
   read(FN,*,IOSTAT=Reason) numpar%SHI_inelast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Charge state: 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande, 4=fixed Zeff, 5=charge exchange:
   read(FN,*,IOSTAT=Reason) numpar%SHI_ch_st
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! SHI charge shape: 0=point-like charge; 1=Brandt-Kitagawa ion:
   read(FN,*,IOSTAT=Reason) numpar%SHI_ch_shape
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! [eV] Cut-off energy (SHIs with lower energies are excluded from calculation):
   read(FN,*,IOSTAT=Reason) numpar%SHI_Cutoff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   !==================================================
   ! Skip line:
   read(FN,*,IOSTAT=Reason) !::: MODELS FOR POSITRONS :::
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Positron inelastic scattering: 0=excluded, 1=delta
   read(FN,*,IOSTAT=Reason) numpar%Pos_inelast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Positron elastic scattering: 0=excluded, 1=delta
   read(FN,*,IOSTAT=Reason) numpar%Pos_elastic
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Positron Bremsstrahlung scattering: 0=excluded, 1=delta
   read(FN,*,IOSTAT=Reason) numpar%Pos_Brems
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Positron annihilation: 0=excluded, 1= Heitler
   read(FN,*,IOSTAT=Reason) numpar%Pos_annih
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! [eV] Cut-off energy (electrons with lower energies are excluded from calculation):
   read(FN,*,IOSTAT=Reason) numpar%Pos_Cutoff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   !==================================================
   ! Skipe line
   read(FN,*,IOSTAT=Reason) !::: MODEL PARAMETERS FOR CDF :::
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! CDF model: 0 = Drude / Lindhard CDF, 1=Mermin CDF, 2=Full conserving CDF:
   read(FN,*,IOSTAT=Reason) numpar%CDF_model
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! target dispersion relation: 0=free electron, 1=plasmon-pole, 2=Ritchie; effective mass [in me] (0=effective mass from DOS of VB; -1=free-electron):
   read(FN,*,IOSTAT=Reason) numpar%CDF_dispers, numpar%CDF_m_eff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Include plasmon integration limit (0=no, 1=yes):
   read(FN,*,IOSTAT=Reason) numpar%CDF_plasmon
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Coefficient where to switch from Ritchie to Delta CDF: E = k * Wmin (INELASTIC):
   read(FN,*,IOSTAT=Reason) numpar%CDF_Eeq_factor
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Coefficient where to switch from Ritchie to Delta CDF: E = k * Wmin (ELASTIC):
   read(FN,*,IOSTAT=Reason) numpar%CDF_Eeq_elast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
    ! Flag to use for target atoms Zeff (set 0), or Z=1 (set 1):
   read(FN,*,IOSTAT=Reason) numpar%CDF_elast_Zeff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! effective number of grid points for inelastic cross section integration over energy
   read(FN,*,IOSTAT=Reason) numpar%CDF_int_n_inelastE
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! effective number of grid points for inelastic cross section integration over momentum
   read(FN,*,IOSTAT=Reason) numpar%CDF_int_n_inelastQ
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! effective number of grid points for elastic cross section integration over energy
   read(FN,*,IOSTAT=Reason) numpar%CDF_int_n_elastE
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! effective number of grid points for elastic cross section integration over momentum
   read(FN,*,IOSTAT=Reason) numpar%CDF_int_n_elastQ
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   
   !==================================================
   ! Skip line
   read(FN,*,IOSTAT=Reason) !::: MODELS FOR HOLES :::
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Auger decays:  0=excluded, 1=EADL
   read(FN,*,IOSTAT=Reason) numpar%H_Auger
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Radiative decays: 0=excluded, 1=EADL
   read(FN,*,IOSTAT=Reason) numpar%H_Radiat
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! [me] effective valence hole mass in units of electron mass
   read(FN,*,IOSTAT=Reason) numpar%H_m_eff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Valence hole inelastic scattering: 0=excluded, 1=CDF, 2=BEB, 3=Delta
   read(FN,*,IOSTAT=Reason) numpar%H_inelast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Valence hole elastic scattering: 0=excluded, 1=Mott
   read(FN,*,IOSTAT=Reason) numpar%H_elast
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! [eV] Cut-off energy (holes with lower energies are excluded from calculation):
   read(FN,*,IOSTAT=Reason) numpar%H_Cutoff
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   !==================================================
   ! Skip line
   read(FN,*,IOSTAT=Reason)  ! ::: MD MODEL PARAMETERS :::
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Quenching in MD:
   read(FN,*,IOSTAT=Reason) numpar%do_quenching, numpar%t_quench_start, numpar%dt_quench
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   numpar%t_quench_run = 0.0d0  ! starting
   
   !==================================================
   ! Skip line
   read(FN,*,IOSTAT=Reason)  ! ::: OUTPUT DATA :::
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Which format to use for gnuplot figures:
   read(FN,*,IOSTAT=Reason) numpar%gnupl%gnu_extension
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Does user require to printout DOS of materials?
   read(FN,*,IOSTAT=Reason) numpar%printout_DOS
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Does user require to printout MFPs in materials?
   read(FN,*,IOSTAT=Reason) numpar%printout_MFPs
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif

   ! Does user require to printout Ranges in materials?
   read(FN,*,IOSTAT=Reason) numpar%printout_ranges
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      call Save_error_details(Err, 2, Error_descript)	! module "Objects"
      goto 9998
   endif
   
   ! Which spatial grids to use for printing output data:
   call read_output_grid_coord(FN, File_name, numpar, Err, count_lines)  ! below
   
9998 continue
end subroutine read_num_pars



subroutine set_MD_step_grid(File_name, numpar, read_well_out, Error_descript)
   character(*), intent(in) :: File_name    ! file name with input data
   type(Num_par), intent(inout), target :: numpar	! all numerical parameters
   logical, intent(inout) :: read_well_out
   character(*), intent(inout) :: Error_descript
   !-----------------------------------------
   character(200) :: Path, full_file_name
   logical :: file_exist, read_well
   integer :: FN, Nsiz, count_lines, Reason, i

   read_well_out = .false.   ! to start with
   Path = trim(adjustl(m_input_folder))//numpar%path_sep    ! where to find the file with the data
   full_file_name = trim(adjustl(Path))//trim(adjustl(File_name))   ! to read the file from the INPUT_DATA directory
   inquire(file=trim(adjustl(full_file_name)),exist=file_exist) ! check if input file is there
   if (file_exist) then ! try to read it, if there is a grid provided
      open(newunit = FN, FILE = trim(adjustl(full_file_name)), status = 'old', action='read')
      ! Find the grid size from the file:
      call Count_lines_in_file(FN, Nsiz) ! module "Dealing_with_files"
      ! Knowing the size, create the grid-array and read from the file:
      if (allocated(numpar%dt_MD_reset_grid)) deallocate(numpar%dt_MD_reset_grid) ! make sure it's possible to allocate
      allocate(numpar%dt_MD_reset_grid(Nsiz)) ! allocate it
      if (allocated(numpar%dt_MD_grid)) deallocate(numpar%dt_MD_grid) ! make sure it's possible to allocate
      allocate(numpar%dt_MD_grid(Nsiz)) ! allocate it

      ! Read data on the grid from the file:
      count_lines = 0   ! just to start counting lines in the file
      do i = 1, Nsiz    ! read grid line by line from the file
         read(FN,*,IOSTAT=Reason) numpar%dt_MD_reset_grid(i), numpar%dt_MD_grid(i)    ! grid data from the file
         call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
         if (.not. read_well) then ! something wrong with the user-defined grid
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            goto 9993  ! couldn't read the data, exit the cycle
         endif
      enddo
      read_well_out = .true.    ! we read the grid from the file well
      numpar%i_dt = 0   ! to start with
      numpar%dt_MD = numpar%dt_MD_grid(1)   ! to start with

!       print*, 'set_MD_step_grid:'
!       print*, numpar%dt_MD_reset_grid(:)
!       print*, 'dt=', numpar%dt_MD_grid(:)

9993 call close_file('close', FN=FN) ! module "Dealing_with_files"
   else ! If there is no input file, check if the teimstep is provided instead:
      count_lines = 0   ! just to start counting lines in the file
      read(File_name,*,IOSTAT=Reason) numpar%dt_MD    ! grid data from the file
      call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
      if (read_well) read_well_out = .true.    ! we read the grid from the file well
      numpar%i_dt = -1   ! to mark that the reset option is unused
   endif ! file_exist
!    pause 'set_MD_step_grid'
end subroutine set_MD_step_grid



subroutine set_time_grid(File_name, numpar)
   character(*), intent(in) :: File_name    ! file name with input data
   type(Num_par), intent(inout), target :: numpar	! all numerical parameters
   !-----------------------------------------
   character(200) :: Path, full_file_name
   logical :: file_exist, read_well, setgrid
   integer :: FN, Nsiz, Ncol, count_lines, Reason, i
   
   setgrid = .false.   ! to start with 
   Path = trim(adjustl(m_input_folder))//numpar%path_sep    ! where to find the file with the data
   full_file_name = trim(adjustl(Path))//trim(adjustl(File_name))   ! to read the file from the INPUT_DATA directory
   inquire(file=trim(adjustl(full_file_name)),exist=file_exist) ! check if input file is there
   if (file_exist) then ! try to read it, if there is a grid provided
      open(newunit = FN, FILE = trim(adjustl(full_file_name)), status = 'old', action='read')
      ! Find the grid size from the file:
      call Count_lines_in_file(FN, Nsiz) ! module "Dealing_with_files"
      ! Knowing the size, create the grid-array and read from the file:
      if (allocated(numpar%dt_out_grid)) deallocate(numpar%dt_out_grid) ! make sure it's possible to allocate
      allocate(numpar%dt_out_grid(Nsiz)) ! allocate it
      ! Check if there is a column with resetting timestep:
      call Count_columns_in_file(FN, Ncol) ! module "Dealing_with_files"
      if (Ncol > 1) then    ! set the grid for resetting timestep of the simulation
         if (allocated(numpar%dt_reset_grid)) deallocate(numpar%dt_reset_grid) ! make sure it's possible to allocate
         allocate(numpar%dt_reset_grid(Nsiz)) ! allocate it
      endif
      ! Read data on the grid from the file:
      count_lines = 0   ! just to start counting lines in the file
      do i = 1, Nsiz    ! read grid line by line from the file
          if (Ncol > 1) then ! read two Count_columns
             read(FN,*,IOSTAT=Reason) numpar%dt_out_grid(i), numpar%dt_reset_grid(i)    ! grid data from the file
          else  ! read one column
             read(FN,*,IOSTAT=Reason) numpar%dt_out_grid(i)    ! grid data from the file
          endif
          call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
!           print*, i, numpar%dt_out_grid(i), numpar%dt_reset_grid(i)

          if (.not. read_well) then ! something wrong with the user-defined grid
             write(6,'(a,i3)') 'In the file '//trim(adjustl(full_file_name))//' could not read line ', count_lines
             write(6,'(a)') 'Using default grid instead.'
             deallocate(numpar%dt_out_grid)
             if (allocated(numpar%dt_reset_grid)) deallocate(numpar%dt_reset_grid)
             goto 9994  ! couldn't read the data, exit the cycle
          endif
      enddo
      setgrid = .true.    ! we read the grid from the file well

9994 call close_file('close', FN=FN) ! module "Dealing_with_files"
   else ! If there is no input file, check if the timestep is provided instead:
      count_lines = 0   ! just to start counting lines in the file
      read(File_name,*,IOSTAT=Reason) numpar%dt_printout    ! grid data from the file
      call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
      if (read_well) setgrid = .true.    ! we read the grid from the file well
!       print*, 'Set constant dt_out:' , numpar%dt_printout
   endif ! file_exist
   
   ! In case something went wrong while reading the file, set default:
   if (.not.setgrid) then! if the number is not provided in a correct format, use MD timestep
      write(6,'(a,i3)') 'In the file '//trim(adjustl(full_file_name))//' could not read line ', count_lines
      write(6,'(a)') 'Using MD temistep for printout instead.'
      numpar%dt_printout = numpar%dt_MD ! printout data every timestep by default
   endif
end subroutine set_time_grid


subroutine read_output_grid_coord(FN, File_name, numpar, Err, count_lines)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   type(Num_par), intent(inout), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   integer, intent(inout) :: count_lines    ! where are we reading file right now
   !----------------------
   type(grid_params), pointer:: grid_par     ! grid parameters
   logical :: read_well , grid_created     ! was there an error reading file?
   real(8) :: max_grid_size
   integer :: Reason, i, N_rp, temp, i_ax, grid_ind
   character(200) :: Error_descript
   character(30) :: temp_ch, temp_ch2

   ! to start with:
   temp_ch = ''
   read_well = .true.
   numpar%print_each_step = .false.    ! to start with the default value (do not printout each MD step timing)
   numpar%MD_force_ind = 0 ! Default value of the MD force calculator
   
   ! Define default parameters for the grids:
   call set_default_grid_params(numpar) ! below
   
   ! Check what user defined:
   GR:do while (read_well)	! check all possible entries (valence band and phonons)
      read(FN,*,IOSTAT=Reason) temp_ch  ! check whether there are grids for printout
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then
         exit	! nothing more to read from the file
      else  !  (.not. read_well)
         
         select case (trim(adjustl(temp_ch)))

         !============================================
         ! Printout atomic coordinates from MD in XYZ format:
         case ('Displacement', 'DISPLACEMENT', 'displacement', 'MSD', 'msd')
            read(FN,*,IOSTAT=Reason) numpar%n_MSD   ! power of mean displacement to print out (set integer N: <u^N>-<u0^N>)
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"

         !============================================
         ! Printout atomic coordinates from MD in XYZ format:
         case ('Numerical_potential', 'Numeric_potential', 'Numeric_pot', 'numeric_pot', 'numeric_force', 'Numeric_force')
            numpar%MD_force_ind = 1 ! calculate forces as numerical derivative of the potential (SLOW)

         !============================================
         ! Printout atomic coordinates from MD in XYZ format:
         case ('Cohesive', 'COHESIVE', 'cohesive', 'Cohesive_energy', 'COHESIVE_ENERGY', 'cohesive_energy')
            numpar%do_cohesive = .true.    ! calculate cohesive energy (instead of full TREKIS run)
         
         !============================================
         ! Printout atomic coordinates from MD in XYZ format:
         case ('Print_XYZ', 'print_XYZ', 'Print_xyz', 'print_xyz', 'PrintXZY', 'PRINT_XYZ')
            numpar%print_MD_R_xyz = .true.    ! printout MD atomic coordinates in XYZ

         !============================================
         ! Printout atomic velosities from MD in XYZ format:
         case ('Print_V_XYZ', 'print_V_XYZ', 'Print_v_xyz', 'print_v_xyz', 'PrintVXZY', 'PRINT_V_XYZ')
            numpar%print_MD_V_xyz = .true.    ! printout MD atomic velosities in XYZ

         !============================================
         ! Create LAMMPS input file with atmic state at the final instant of time:
         case ('Print_LAMMPS', 'print_LAMMPS', 'Print_Lammps', 'print_lammps', 'PrintLAMMPS', 'PRINT_LAMMPS')
            numpar%print_MD_LAMMPS = .true.    ! create input file for LAMMPS at the final time instant
            read(FN,*,IOSTAT=Reason) numpar%LAMMPS_UNITS   ! Read LAMMPS units to be used for LAMMPS input file
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"

         !============================================
         ! Printout energy passed from MC to MD:
         case ('Print_MC_MD_energy', 'print_MC_MD_energy', 'PRINT_MC_MD_ENERGY', 'print_mc_md_energy', 'Print_MC_MD_Energy')
            numpar%print_MC_MD_energy = .true.    ! printout MC-MD energy transfer

         !============================================
         ! Printout each timestep:
         case ('PrintTheta', 'print_theta', 'print_Theta', 'Print_theta', 'Print_Theta', 'Printtheta', 'printtheta', &
                'printTheta')    ! printout particles theta-distribution
            grid_par => numpar%vel_theta_grid_par   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes; theta=acos(Vz/V)
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default velosity theta grid instead'
               
               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, 1)  ! below
               
               ! Create grid:
               if (grid_par%log_scale(1)) then   ! log-scale grid along theta
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%vel_theta_grid, 1)  ! module "Little_subroutines"
               else    ! linear scale along E
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%vel_theta_grid, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%vel_theta_grid, grid_created)  ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, 1)  ! below

                  ! Create grid:
                  if (grid_par%log_scale(1)) then   ! log-scale grid along theta
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%vel_theta_grid, 1)  ! module "Little_subroutines"
                  else    ! linear scale along E
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%vel_theta_grid, 0)  ! module "Little_subroutines"
                  endif
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
            
         !============================================
         ! Printout each timestep:
         case ('PrintMDSteps', 'printMDsteps', 'PrintMDsteps', 'Print_MD_Steps', 'Print_MD_steps', 'print_MD_steps', &
                'Print_each_step', 'print_each_step')    ! printout each MD step timeing
            numpar%print_each_step = .true.    ! printout each MD step timing

         !============================================
         ! Set energy grid for spectra printout:
         case ('ENERGY', 'Energy', 'NRG', 'energy', 'nrg', 'E', 'e')    ! Energy grid parameters
            grid_par => numpar%NRG_grid_par   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes
            
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default energy grid instead'
               
               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, 1)  ! below
               
               ! Create grid:
               if (grid_par%log_scale(1)) then   ! log-scale grid along E
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%NRG_grid, 1)  ! module "Little_subroutines"
               else    ! linear scale along E
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%NRG_grid, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%NRG_grid, grid_created)  ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, 1)  ! below

                  ! Create grid:
                  if (grid_par%log_scale(1)) then   ! log-scale grid along E
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%NRG_grid, 1)  ! module "Little_subroutines"
                  else    ! linear scale along E
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%NRG_grid, 0)  ! module "Little_subroutines"
                  endif
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
            
         ! Set spatial grids for distributions printout:
         case ('Spectra_x', 'SPECTRA_X', 'SPECTRA_x', 'spectra_x', 'spectra_X', 'Spectra_X')    ! Spectra along X-axis (Cartesian)
            i_ax = 1    ! index for this type of printout
            grid_ind = 1 ! X

            grid_par => numpar%Spectr_grid_par(i_ax)   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes
            
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default Spactra X-grid instead'
               
               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below

               ! Create grid:
               if (grid_par%log_scale(1)) then   ! log-scale grid along X
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
               else    ! linear scale along R
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%Spectr_grid(i_ax)%spatial_grid1, grid_created)    ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below
                  !Last three parameters are:   grid_type, grid_dim, grid_ind
                  !  grid_type     ! grid type: 0=cartesian, 1=cylindrical,  2=spherical
                  !  grid_dim      ! dimension: 1=1d, 2=2d, 3=3d
                  !  grid_ind       ! which axis: 1=x or R or R;   2=y or L or Theta;  3=z or Theta or Phi  (for Cartesian or Cyllindrical or Spherical)
                  
                  ! Create grid:
                  if (grid_par%log_scale(1)) then   ! log-scale grid along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
                  else    ! linear scale along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
                  endif
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
         
         case ('Spectra_y', 'SPECTRA_Y', 'SPECTRA_y', 'spectra_Y', 'spectra_y', 'Spectra_Y')    ! Spectra along Y-axis (Cartesian)
            i_ax = 2    ! index for this type of printout
            grid_ind = 1    ! Y
            
            grid_par => numpar%Spectr_grid_par(i_ax)   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes
            
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default Spactra Y-grid instead'

               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below
               
               ! Create grid:
               if (grid_par%log_scale(1)) then   ! log-scale grid along Y
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
               else    ! linear scale along R
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%Spectr_grid(i_ax)%spatial_grid1, grid_created)    ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below
                  !Last three parameters are:   grid_type, grid_dim, grid_ind
                  
                  ! Create grid:
                  if (grid_par%log_scale(1)) then   ! log-scale grid along R
                     call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
                  else    ! linear scale along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
                  endif
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
            
            
         case ('Spectra_z', 'SPECTRA_Z', 'SPECTRA_z', 'spectra_Z', 'spectra_z', 'Spectra_Z')    ! Spectra along Z-axis (Cartesian)
            i_ax = 3    ! index for this type of printout
            grid_ind = 1    ! Z
            
            grid_par => numpar%Spectr_grid_par(i_ax)   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes
            
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default Spactra Z-grid instead'

               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below
               
               ! Create grid:
               if (grid_par%log_scale(1)) then   ! log-scale grid along Y
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
               else    ! linear scale along R
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%Spectr_grid(i_ax)%spatial_grid1, grid_created)    ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below
                  !Last three parameters are:   grid_type, grid_dim, grid_ind

                  ! Create grid:
                  if (grid_par%log_scale(1)) then   ! log-scale grid
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
                  else    ! linear scale along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
                  endif
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
            
!             print*, numpar%Spectr_grid(3)%along_axis
!             print*, size(numpar%Spectr_grid(i_ax)%spatial_grid1)
!             pause 'READ GRID SPECTRUM'
!             print*, 'GRID Z',  numpar%Spectr_grid(i_ax)%spatial_grid1(:)
!             pause 

         case ('Spectra_r', 'Spectra_R', 'SPECTRA_r', 'S_PECTRA_R', 'spectra_r')    ! Spctra only along Radius (Cylindric)
            i_ax = 8    ! index for this type of printout
            grid_par => numpar%Spectr_grid_par(i_ax)   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes
            
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default Spactra R-grid instead'
               
               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, 1)  ! below

               ! Create grid:
               if (grid_par%log_scale(1)) then   ! log-scale grid along R
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(8)%spatial_grid1, 1)  ! module "Little_subroutines"
               else    ! linear scale along R
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(8)%spatial_grid1, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%Spectr_grid(8)%spatial_grid1, grid_created)    ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, 1)  ! below
                  ! Create grid:
                  if (grid_par%log_scale(1)) then   ! log-scale grid along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(8)%spatial_grid1, 1)  ! module "Little_subroutines"
                  else    ! linear scale along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%Spectr_grid(8)%spatial_grid1, 0)  ! module "Little_subroutines"
                  endif
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
         
         
         !============================================   
         ! Set spatial grids for distributions printout:
         case ('x', 'X')    ! only along X-axis (Cartesian)
            i_ax = 1    ! index for this type of printout
            grid_ind = 1 ! X
            
            grid_par => numpar%grid_par(i_ax)   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes
            
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default X-grid instead'
               
               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below

               ! Create grid:
               if (grid_par%log_scale(1)) then   ! log-scale grid along X
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
               else    ! linear scale along R
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%grids(i_ax)%spatial_grid1, grid_created)    ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below
                  !Last three parameters are:   grid_type, grid_dim, grid_ind
                  !  grid_type     ! grid type: 0=cartesian, 1=cylindrical,  2=spherical
                  !  grid_dim      ! dimension: 1=1d, 2=2d, 3=3d
                  !  grid_ind       ! which axis: 1=x or R or R;   2=y or L or Theta;  3=z or Theta or Phi  (for Cartesian or Cyllindrical or Spherical)
                  
                  ! Create grid:
                  if (grid_par%log_scale(1)) then   ! log-scale grid along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
                  else    ! linear scale along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
                  endif
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
            
         case ('y', 'Y')    ! only along Y-axis (Cartesian)
            i_ax = 2    ! index for this type of printout
            grid_ind = 1    ! Y
            
            grid_par => numpar%grid_par(i_ax)   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes
            
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default Y-grid instead'

               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below
               
               ! Create grid:
               if (grid_par%log_scale(1)) then   ! log-scale grid along Y
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
               else    ! linear scale along R
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%grids(i_ax)%spatial_grid1, grid_created)    ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below
                  !Last three parameters are:   grid_type, grid_dim, grid_ind
                  
                  ! Create grid:
                  if (grid_par%log_scale(1)) then   ! log-scale grid along R
                     call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
                  else    ! linear scale along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
                  endif
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
            
            
         case ('z', 'Z')    ! only along Z-axis (Cartesian)
            i_ax = 3    ! index for this type of printout
            grid_ind = 1    ! Z
            
            grid_par => numpar%grid_par(i_ax)   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes
            
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default Z-grid instead'

               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below
               
               ! Create grid:
               if (grid_par%log_scale(1)) then   ! log-scale grid along Z
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
               else    ! linear scale along R
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%grids(i_ax)%spatial_grid1, grid_created)    ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 0, 1, grid_ind)  ! below
                  !Last three parameters are:   grid_type, grid_dim, grid_ind
                  
                  ! Create grid:
                  if (grid_par%log_scale(1)) then   ! log-scale grid along R
                     call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
                  else    ! linear scale along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
                  endif
!                   print*, 'TEST READING:', grid_ind, allocated(numpar%grids(i_ax)%spatial_grid1)
!                   pause 'create_grid'
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
            
            
         case ('xy', 'XY', 'Xy', 'Yx')  ! only along X and Y-axes (Cartesian)
            i_ax = 4    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('xz', 'XZ', 'Xz', 'Zx')   ! only along X and Z-axes (Cartesian)
            i_ax = 5    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('yz', 'YZ', 'Yz', 'Zy')   ! only along Y and Z-axes (Cartesian)
            i_ax = 6    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('xyz', 'XYZ', 'Xyz', 'XYz', 'yxz', 'YXZ', 'zyx', 'ZYX', 'xzy', 'XZY', 'zxy', 'ZXY')     ! along X, Y and Z-axes (Cartesian)
            i_ax = 7    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('r', 'R')    ! only along Radius (Cylindric)
            i_ax = 8    ! index for this type of printout
            grid_par => numpar%grid_par(i_ax)   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes
            
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default R-grid instead'
               
               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, 1)  ! below

               ! Create grid:
               if (grid_par%log_scale(1)) then   ! log-scale grid along R
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(8)%spatial_grid1, 1)  ! module "Little_subroutines"
               else    ! linear scale along R
                  call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(8)%spatial_grid1, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%grids(8)%spatial_grid1, grid_created)    ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, 1)  ! below
                  ! Create grid:
                  if (grid_par%log_scale(1)) then   ! log-scale grid along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(8)%spatial_grid1, 1)  ! module "Little_subroutines"
                  else    ! linear scale along R
                     call create_grid(grid_par%gridstart(1), grid_par%gridend(1), grid_par%gridstep(1), numpar%grids(8)%spatial_grid1, 0)  ! module "Little_subroutines"
                  endif
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
            
         case ('L', 'l', 'Depth', 'depth')   ! only along L (Cylindric)
            i_ax = 9    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('phi', 'Phi', 'PHI')       ! only along phi (Cylindric)
            i_ax = 10    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('rl', 'RL', 'Rl', 'lr', 'LR', 'Lr')      ! only along Radius and Length (Cylindric)
            i_ax = 11    ! index for this type of printout
            
            grid_par => numpar%grid_par(i_ax)   ! just to access easier
            grid_par%along_axis = .true.     ! the user whats to printout the data along this axes
            
            !Reading grid for the firts dimension - R
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            grid_ind = 1        ! R
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default R-grid instead'

               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, grid_ind)  ! below
               
               ! Create grid:
               if (grid_par%log_scale(grid_ind)) then   ! log-scale grid along Z
                  call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(grid_ind), grid_par%gridstep(grid_ind), numpar%grids(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
               else    ! linear scale along R
                  call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(grid_ind), grid_par%gridstep(grid_ind), numpar%grids(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%grids(i_ax)%spatial_grid1, grid_created)    ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, 1)  ! below
                  !Last three parameters are:   grid_type, grid_dim, grid_ind
                  
                  ! Create grid:
                  if (grid_par%log_scale(grid_ind)) then   ! log-scale grid along R
                     call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(grid_ind), grid_par%gridstep(grid_ind), numpar%grids(i_ax)%spatial_grid1, 1)  ! module "Little_subroutines"
                  else    ! linear scale along R
                     call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(grid_ind), grid_par%gridstep(grid_ind), numpar%grids(i_ax)%spatial_grid1, 0)  ! module "Little_subroutines"
                  endif
!                   print*, 'TEST READING:', grid_ind, allocated(numpar%grids(i_ax)%spatial_grid1)
!                   pause 'create_grid'
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
            
            !Reading grid along second dimension - L
            read(FN,*,IOSTAT=Reason) temp_ch2   ! Read the grid parameter: either file name with the grid, or its type
            call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
            grid_ind = 2    ! L
            if (.not. read_well) then   ! user did not provide correct parameters
               write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
               write(6,'(a)') trim(adjustl(Error_descript))
               write(6,'(a)') ' Using default L-grid instead'

               backspace ( FN ) ! to read the line again from the file into a proper variable
               ! Read the parameters of the grid:
               call read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, grid_ind)  ! below
               
               ! Create grid:
               if (grid_par%log_scale(grid_ind)) then   ! log-scale grid along L
                  call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(grid_ind), grid_par%gridstep(grid_ind), numpar%grids(i_ax)%spatial_grid2, 1)  ! module "Little_subroutines"
               else    ! linear scale along R
                  call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(grid_ind), grid_par%gridstep(grid_ind), numpar%grids(i_ax)%spatial_grid2, 0)  ! module "Little_subroutines"
               endif
            else ! (.not. read_well)   ! if user provided grid parameters
               call read_grid_from_file(temp_ch2, numpar, numpar%grids(i_ax)%spatial_grid2, grid_created)    ! below
               ! Check the grid was not read from the file:
               if (.not.grid_created) then   ! create a default grid
                  backspace ( FN ) ! to read the line again from the file into a proper variable
                  ! Read the parameters of the grid:    Last three parameters are:   grid_type, grid_dim, grid_ind
                  call  read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, 1, 1, grid_ind)  ! below
                  
                  ! Create grid:
                  if (grid_par%log_scale(grid_ind)) then   ! log-scale grid along L
                     call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(grid_ind), grid_par%gridstep(grid_ind), numpar%grids(i_ax)%spatial_grid2, 1)  ! module "Little_subroutines"
                  else    ! linear scale along L
                     call create_grid(grid_par%gridstart(grid_ind), grid_par%gridend(grid_ind), grid_par%gridstep(grid_ind), numpar%grids(i_ax)%spatial_grid2, 0)  ! module "Little_subroutines"
                  endif
               endif    ! (.not.grid_created) 
            endif   ! (.not. read_well)
            
         case ('rtheta', 'RTHETA', 'RTheta', 'Rtheta', 'thetar', 'THETAR', 'THETAr')      ! only along Radius and Theta (Cylindric)
            i_ax = 12    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('rltheta', 'RLTHETA', 'RLTheta', 'Rltheta', 'lrtheta', 'LRTHETA', 'Lrtheta')      ! along Radius, Length and Theta (Cylindric)
            i_ax = 13    ! index for this type of printout
            print*, 'NOT READY YET!'           
            
         case ('rs', 'RS', 'Rs', 'rS')    ! only along Radius (Spherical)
            i_ax = 14    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('thetas', 'Thetas', 'THETAS', 'THETAs')    ! only along theta (Spherical)
            i_ax = 15    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('phis', 'Phis', 'PHIS', 'PHIs') ! only along phi (Spherical)
            i_ax = 16    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('rphis', 'rphiS', 'RPhis', 'RPHIS', 'RPhiS') ! only along R and Phi (Spherical)
            i_ax = 17    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('rthetas', 'rthetaS', 'RThetas', 'RTHETAS', 'RThetaS') ! only along R and theta (Spherical)
            i_ax = 18    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         case ('rthetaphis', 'rthetaphiS', 'RThetaPhis', 'RTHETAPHIS', 'RThetaPhiS') ! along R, theta and phi (Spherical)
            i_ax = 19    ! index for this type of printout
            print*, 'NOT READY YET!'
            
         end select
         
      endif !  (.not. read_well)
   enddo GR

9997 continue
   
!     print*, 'GRID E:', numpar%NRG_grid_par%gridstart(1), numpar%NRG_grid_par%gridend(1), numpar%NRG_grid_par%gridstep(1)
!     do i = 1, size(numpar%NRG_grid(:))
!        print*, i, numpar%NRG_grid(i)
!     enddo
!     print*, 'GRID R:', numpar%grid_par(8)%gridstart(1), numpar%grid_par(8)%gridend(1), numpar%grid_par(8)%gridstep(1)
!     do i = 1, size(numpar%grids(8)%spatial_grid(1,:))
!        print*, i, numpar%grids(8)%spatial_grid(1,i)
!     enddo
!     pause 'read_output_grid_coord'

   nullify(grid_par)
end subroutine read_output_grid_coord


subroutine read_grid_parameters(FN, File_name, count_lines, grid_par, read_well, grid_type, grid_dim, grid_ind)
   integer, intent(in) :: FN	! file number to read from
   character(*), intent(in) :: File_name	! file name with input data
   integer, intent(inout) :: count_lines    ! where are we reading file right now
   type(grid_params), intent(inout):: grid_par     ! grid parameters
   logical, intent(inout) :: read_well  ! did we read it well
   integer, intent(in) :: grid_type     ! grid type: 0=cartesian, 1=cylindrical,  2=spherical
   integer, intent(in) :: grid_dim      ! dimension: 1=1d, 2=2d, 3=3d
   integer, intent(in) :: grid_ind       ! which axis: 1=x or R or R;   2=y or L or Theta;  3=z or Theta or Phi  (for Cartesian or Cyllindrical or Spherical)
   !---------------------------------------------------   
   character(200) :: Error_descript
   integer :: Reason, ind, dims
   
   ! Check if the dimension provided is correct:
   if ((grid_ind > 3) .or. (grid_ind < 1)) then
      ! by default, it will be:
      select case (grid_type)
      case default  ! Cartesian:
         ind = 3   ! along Z
      case (1)  ! Cyllindrical
         ind = 1    ! radial
      case (2)  ! Spherical
         ind = 1    ! radial
      endselect
   else
      ind = grid_ind
   endif
   
   if ((grid_dim > 3) .or. (grid_dim < 1)) then
      ! by default, it will be 1 d:
      dims = 1
   else
      dims= grid_dim
   endif
   
   ! Read the lines in predefined order:
   read(FN,*,IOSTAT=Reason) grid_par%log_scale(ind)  ! printout spatial data in log-scale?
   call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
   if (.not. read_well) then    ! incorrect parameterization provided
      write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
      write(6,'(a)') trim(adjustl(Error_descript))
      write(6,'(a)') ' Using default grid instead'
      goto 9995
   endif
   
   ! Order of variables for:
   ! a) Cartesian: X, Y, Z
   ! b) Cyllindrical: R, L, theta
   ! c) Spherical: R, Theta, Phi
   if (dims == 3) then   ! all dimensions:
      read(FN,*,IOSTAT=Reason) grid_par%gridstart(1),  grid_par%gridend(1),  grid_par%gridstep(1)  ! start, end, step
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then    ! incorrect parameterization provided
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         write(6,'(a)') trim(adjustl(Error_descript))
         write(6,'(a)') ' Using default grid instead'
         goto 9995
      endif
      
      read(FN,*,IOSTAT=Reason) grid_par%gridstart(2),  grid_par%gridend(2),  grid_par%gridstep(2)  ! start, end, step
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then    ! incorrect parameterization provided
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         write(6,'(a)') trim(adjustl(Error_descript))
         write(6,'(a)') ' Using default grid instead'
         goto 9995
      endif
      
      read(FN,*,IOSTAT=Reason) grid_par%gridstart(3),  grid_par%gridend(3),  grid_par%gridstep(3)  ! start, end, step
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then    ! incorrect parameterization provided
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         write(6,'(a)') trim(adjustl(Error_descript))
         write(6,'(a)') ' Using default grid instead'
         goto 9995
      endif
         
   elseif (dims == 2) then    ! 2d grid: xy=1, xz=2, yz=3;   RL=1, RTheta=2, LTheta=3;   RTheta=1, RPhi=2, ThetaPhi=3;
      
      select case (ind)
      case default   !  xy, or RL, or RTheta
         read(FN,*,IOSTAT=Reason) grid_par%gridstart(1),  grid_par%gridend(1),  grid_par%gridstep(1)  ! start, end, step
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then    ! incorrect parameterization provided
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            write(6,'(a)') trim(adjustl(Error_descript))
            write(6,'(a)') ' Using default grid instead'
            goto 9995
         endif
      
         read(FN,*,IOSTAT=Reason) grid_par%gridstart(2),  grid_par%gridend(2),  grid_par%gridstep(2)  ! start, end, step
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then    ! incorrect parameterization provided
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            write(6,'(a)') trim(adjustl(Error_descript))
            write(6,'(a)') ' Using default grid instead'
            goto 9995
         endif
      case (2)   ! xz, or RTheta, or RPhi
         read(FN,*,IOSTAT=Reason) grid_par%gridstart(1),  grid_par%gridend(1),  grid_par%gridstep(1)  ! start, end, step
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then    ! incorrect parameterization provided
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            write(6,'(a)') trim(adjustl(Error_descript))
            write(6,'(a)') ' Using default grid instead'
            goto 9995
         endif
      
         read(FN,*,IOSTAT=Reason) grid_par%gridstart(3),  grid_par%gridend(3),  grid_par%gridstep(3)  ! start, end, step
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then    ! incorrect parameterization provided
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            write(6,'(a)') trim(adjustl(Error_descript))
            write(6,'(a)') ' Using default grid instead'
            goto 9995
         endif
      case (3)   ! yz, or LTheta, or ThetaPhi
         read(FN,*,IOSTAT=Reason) grid_par%gridstart(2),  grid_par%gridend(2),  grid_par%gridstep(2)  ! start, end, step
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then    ! incorrect parameterization provided
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            write(6,'(a)') trim(adjustl(Error_descript))
            write(6,'(a)') ' Using default grid instead'
            goto 9995
         endif
            
         read(FN,*,IOSTAT=Reason) grid_par%gridstart(3),  grid_par%gridend(3),  grid_par%gridstep(3)  ! start, end, step
         call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
         if (.not. read_well) then    ! incorrect parameterization provided
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            write(6,'(a)') trim(adjustl(Error_descript))
            write(6,'(a)') ' Using default grid instead'
            goto 9995
         endif
      end select
      
   else   ! assume 1d
   
      read(FN,*,IOSTAT=Reason) grid_par%gridstart(ind),  grid_par%gridend(ind),  grid_par%gridstep(ind)  ! start, end, step
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
      if (.not. read_well) then    ! incorrect parameterization provided
         write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
         write(6,'(a)') trim(adjustl(Error_descript))
         write(6,'(a)') ' Using default grid instead'
         goto 9995
      endif
   endif  ! (grid_dim == 3)

9995 continue
end subroutine read_grid_parameters                  




subroutine read_grid_from_file_1d(file_name, numpar, E_grid, it_worked)    ! below
   character(*), intent(in) :: file_name    ! file name to check and read the grid from
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), dimension(:), allocatable, intent(inout) :: E_grid  ! grid to create
   logical, intent(out) ::  it_worked   ! did we read it from the file well?
   !-------------------------------------------------
   character(300) :: Path, full_file_name
   integer :: FN, Nsiz, Reason, count_lines, i
   logical :: file_exist, read_well
   
   it_worked = .false.  ! to start with: no grid by default
   
   Path = trim(adjustl(m_input_folder))//numpar%path_sep    ! where to find the file with the data
   full_file_name = trim(adjustl(Path))//trim(adjustl(file_name))   ! to read the file from the INPUT_DATA directory
   inquire(file=trim(adjustl(full_file_name)),exist=file_exist) ! check if input file is there
   if (file_exist) then ! try to read it, if there is a grid provided
      open(newunit = FN, FILE = trim(adjustl(full_file_name)), status = 'old', action='read')
      ! Find the grid size from the file:
      call Count_lines_in_file(FN, Nsiz) ! module "Dealing_with_files"
      ! Knowing the size, create the grid-array and read from the file:
      if (allocated(E_grid)) deallocate(E_grid) ! make sure it's possible to allocate
      allocate(E_grid(Nsiz))
      ! Read data on the grid from the file:
      count_lines = 0   ! just to start counting lines in the file
      do i = 1, Nsiz    ! read grid line by line from the file
          read(FN,*,IOSTAT=Reason) E_grid(i)    ! grid data from the file
          call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
          if (.not. read_well) then ! something wrong with the user-defined grid
             write(6,'(a,i3)') 'In the file '//trim(adjustl(full_file_name))//' could not read line ', count_lines
             write(6,'(a)') 'Using default spatial grid instead.'
             deallocate(E_grid)
             goto 9996  ! couldn't read the data, exit the cycle
          endif
      enddo
      it_worked = .true.    ! we read the grid from the file well

9996 call close_file('close', FN=FN) ! module "Dealing_with_files"
   endif ! file_exist
end subroutine read_grid_from_file_1d

   
   

subroutine read_grid_from_file_2d(file_name, numpar, spatial_grid, ind, it_worked)    ! below
   character(*), intent(in) :: file_name    ! file name to check and read the grid from
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   real(8), dimension(:,:), allocatable, intent(inout) :: spatial_grid  ! grid to create
   integer, intent(in) :: ind   ! which dimensions of the grid
   logical, intent(out) ::  it_worked   ! did we read it from the file well?
   !-------------------------------------------------
   character(200) :: Path, full_file_name
   integer :: FN, Nsiz, Reason, count_lines, i
   logical :: file_exist, read_well
   
   it_worked = .false.  ! to start with: no grid by default
   
   Path = trim(adjustl(m_input_folder))//numpar%path_sep    ! where to find the file with the data
   full_file_name = trim(adjustl(Path))//trim(adjustl(file_name))   ! to read the file from the INPUT_DATA directory
   inquire(file=trim(adjustl(full_file_name)),exist=file_exist) ! check if input file is there
   if (file_exist) then ! try to read it, if there is a grid provided
      open(newunit = FN, FILE = trim(adjustl(full_file_name)), status = 'old', action='read')
      ! Find the grid size from the file:
      call Count_lines_in_file(FN, Nsiz) ! module "Dealing_with_files"
      ! Knowing the size, create the grid-array and read from the file:
      if (allocated(spatial_grid)) deallocate(spatial_grid) ! make sure it's possible to allocate
      allocate(spatial_grid(3,Nsiz))
      ! Read data on the grid from the file:
      count_lines = 0   ! just to start counting lines in the file
      do i = 1, Nsiz    ! read grid line by line from the file
          read(FN,*,IOSTAT=Reason) spatial_grid(ind,i)    ! grid data from the file
          call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
          if (.not. read_well) then ! something wrong with the user-defined grid
             write(6,'(a,i3)') 'In the file '//trim(adjustl(full_file_name))//' could not read line ', count_lines
             write(6,'(a)') 'Using default spatial grid instead.'
             deallocate(spatial_grid)
             goto 9996  ! couldn't read the data, exit the cycle
          endif
      enddo
      it_worked = .true.    ! we read the grid from the file well

9996 call close_file('close', FN=FN) ! module "Dealing_with_files"
   endif ! file_exist
end subroutine read_grid_from_file_2d



subroutine set_default_grid_params(numpar)
   type(Num_par), intent(inout) :: numpar	! all numerical parameters
   
   ! Set default MD parameters:
   numpar%do_cohesive = .false.         ! don't calculate cohesive energy
   numpar%print_MD_R_xyz = .false.      ! don't printout XYZ file with atomic coordinates
   numpar%print_MD_V_xyz = .false.      ! don't printout XYZ file with atomic velocities
   numpar%print_MC_MD_energy = .false.  ! don't printout MC-MD energy transfer
   numpar%print_MD_LAMMPS = .false.     ! don't create LAMMPS input data
   numpar%n_MSD = 1.0d0    ! power of mean displacement to print out (set integer N: <u^N>-<u0^N>)
   
   ! Set default energy grid parameters:
   numpar%NRG_grid_par%along_axis = .false.         ! does the user want to printout spectra?
   numpar%NRG_grid_par%log_scale(1) = .true.       ! printout energy grid in log or linear scale
   numpar%NRG_grid_par%log_scale(2:3) = .false.    ! dummy argument, unused
   numpar%NRG_grid_par%gridstep(1) = 1.0d0         ! [eV] grid step
   numpar%NRG_grid_par%gridstart(1) = 0.1d0        ! [eV] to start energy axis
   numpar%NRG_grid_par%gridend(1) =  1.0d11        ! [eV] ending point
   
   ! Set default parameters for velosity theta distribution:
   numpar%vel_theta_grid_par%along_axis = .false.        ! does the user want to printout spectra?
   numpar%vel_theta_grid_par%log_scale(1) = .false.      ! printout energy grid in log or linear scale
   numpar%vel_theta_grid_par%log_scale(2:3) = .false.    ! dummy argument, unused
   numpar%vel_theta_grid_par%gridstep(1) = 1.0d0         ! [deg] grid step
   numpar%vel_theta_grid_par%gridstart(1) = 0.0d0        ! [deg] to start theta axis
   numpar%vel_theta_grid_par%gridend(1) =  180.0d0       ! [deg] ending point
  
   ! Set defaults spatial grid parameters:
   numpar%grid_par(:)%along_axis = .false.         ! does the user what to printout the data along this axes combinations?
   numpar%grid_par(:)%log_scale(1) = .false.       ! printout spatial grids along 3 axes in log or linear scale, for each grid index
   numpar%grid_par(:)%log_scale(2) = .false.       ! printout spatial grids along 3 axes in log or linear scale, for each grid index
   numpar%grid_par(:)%log_scale(3) = .false.       ! printout spatial grids along 3 axes in log or linear scale, for each grid index
   
   ! Set defaults cartesian grids:
   numpar%grid_par(1:7)%gridstep(1) = 1.0d0            ! [A] grid step along each used axis in case of linear scales
   numpar%grid_par(1:7)%gridstep(2) = 1.0d0            ! [A] grid step along each used axis in case of linear scales
   numpar%grid_par(1:7)%gridstep(3) = 1.0d0            ! [A] grid step along each used axis in case of linear scales
   numpar%grid_par(1:7)%gridstart(1) = numpar%box_start_x      ! [A] to start X-axis for the spatial grid at zero
   numpar%grid_par(1:7)%gridend(1) =  numpar%box_end_x        ! [A] ending points
   numpar%grid_par(1:7)%gridstart(2) = numpar%box_start_y      ! [A] to start Y-axis for the spatial grid at zero
   numpar%grid_par(1:7)%gridend(2) =  numpar%box_end_y        ! [A] ending points
   numpar%grid_par(1:7)%gridstart(3) = numpar%box_start_z      ! [A] to start Z-axis for the spatial grid at zero
   numpar%grid_par(1:7)%gridend(3) =  numpar%box_end_z        ! [A] ending points
   
   ! Set defaults cyllindrical grids:
   numpar%grid_par(8:13)%gridstep(1) = 1.0d0            ! [A] grid step along each used axis in case of linear scales
   numpar%grid_par(8:13)%gridstep(2) = numpar%box_end_z/100.0d0            ! [A] grid step along each used axis in case of linear scales
   numpar%grid_par(8:13)%gridstep(3) = g_2Pi/100.0d0            ! grid step along each used axis in case of linear scales
   numpar%grid_par(8:13)%gridstart(1) = 0.0d0      ! [A] Radius starting
   numpar%grid_par(8:13)%gridend(1) = numpar%box_end_z      ! [A] Radius ending
   numpar%grid_par(8:13)%gridstart(2) = 0.0d0      ! [A] L
   numpar%grid_par(8:13)%gridend(2) = numpar%box_end_z      ! [A] L
   numpar%grid_par(8:13)%gridstart(3) = 0.0d0      ! [A] theta
   numpar%grid_par(8:13)%gridend(3) = g_2Pi      ! [A] theta
   
   ! Set defaults spherical grids: phi=[0,Pi] counted from Z; theta=[0,2*Pi] within (X,Y) plane
   numpar%grid_par(14:19)%gridstep(1) = 1.0d0      ! [A] R
   numpar%grid_par(14:19)%gridstep(2) = g_2Pi/360.0d0      ! theta
   numpar%grid_par(14:19)%gridstep(3) = g_Pi/180.0d0      ! phi
   numpar%grid_par(14:19)%gridstart(1) = 0.0d0      ! [A] Radius starting
   numpar%grid_par(14:19)%gridend(1) = numpar%box_end_z      ! [A] Radius ending
   numpar%grid_par(14:19)%gridstart(2) = 0.0d0      ! [A] theta
   numpar%grid_par(14:19)%gridend(2) = g_2Pi      ! [A] theta
   numpar%grid_par(14:19)%gridstart(3) = 0.0d0      ! [A] phi
   numpar%grid_par(14:19)%gridend(3) = g_Pi      ! [A] phi
   
end subroutine set_default_grid_params


END MODULE Read_numerical_parameters
