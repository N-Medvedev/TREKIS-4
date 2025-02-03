! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with optical data and convert them into CDF:
MODULE CDF_get_from_data
use Universal_constants
use Objects
use CDF_Ritchi
use CDF_delta, only : define_alpha, define_alpha_OLD
use CS_general_tools, only: MFP_from_sigma
use Read_input_data, only: read_CDF_from_file, m_folder_normalized_CDF, m_input_folder, m_optical_data, m_fitted_coeffs, m_folder_materials, m_photon_absorption
use Little_subroutines, only: extend_array_size_by_one, interpolate_data_single
use Dealing_with_files, only: Count_lines_in_file, read_file
! use CS_integration_limits, only: W_min

implicit none
 contains

 
subroutine get_CDF(numpar, used_target, Err)
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   type(Matter), intent(inout), target :: used_target	! parameters of the target
   type(Error_handling), intent(inout) :: Err	! error log
   !---------------------------------------------------------------
   real(8), dimension(:), allocatable :: lambda, CDF_data
   real(8) :: Omega, ksum, fsum, sigma, elem_contrib, sigma_cur, Wmin
   integer :: i, j, k, m, Nat, Nsiz, FN, FN2, N_CDF, Reason, count_lines, N_elem
   character(200) :: folder, folder_with_cdf, file_with_cdf, command, file_with_coefs, Path_valent, File_name
   logical :: file_exist, read_well
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   path_sep => numpar%path_sep
   folder = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_normalized_CDF))
   
   ! For all target materials:
   TRGT:do i = 1, used_target%NOC		! number of different components of the target (layers, clusters, etc.)   
      ! Total number of elements in this target:            
      N_elem = dble(SUM(used_target%Material(i)%Elements(:)%percentage))	! number of elements in this compound material
                  
      ELMNT:do j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent
         Element => used_target%Material(i)%Elements(j)	! all information about this element
         ! Check if folder for this element already exists:
         folder_with_cdf = trim(adjustl(folder))//path_sep//trim(adjustl(Element%Name))
         inquire(DIRECTORY=trim(adjustl(folder_with_cdf)),exist=file_exist)    ! check if input file excists
         if (.not.file_exist) then	! to make sure that such a folder is present (even if empty)
            ! Create a new directory for output files:
            command='mkdir '//trim(adjustl(folder_with_cdf))	! to create a folder use this command
            CALL system(command)	! create the folder
         endif
         !-------------------------------------
         ! Get MFPs and CDFs for each shell:
         SHL:do m = 1, Element%N_shl		! for all shells in this elements
            Nsiz = size(Element%Phot_absorption%E)	! how many grid points for the cross section we have
            if (.not.allocated(lambda)) allocate(lambda(Nsiz))
            if (.not.allocated(CDF_data)) allocate(CDF_data(Nsiz))
            
            ! Check file with optical data:
            file_with_cdf = trim(adjustl(folder_with_cdf))//path_sep//trim(adjustl(m_optical_data))//'_'//trim(adjustl(Element%Shell_name(m)))//'.dat'
            inquire(file=trim(adjustl(file_with_cdf)),exist=file_exist) ! check if this file is there
            if (file_exist) then	! just read from the file, no need to recalculate:
               open(newunit = FN, FILE = trim(adjustl(file_with_cdf)),action='read')
               ! Get the MFP and CDF points:
               count_lines = 0
               do k = 1, Nsiz
                  E => Element%Phot_absorption%E(k)	! photon energy [eV]
                  read(FN,'(es,es,es)',IOSTAT=Reason) E, lambda(k), CDF_data(k)
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not.read_well) then
                     close(FN)	! redo the file
                     goto 8392
                  endif
               enddo
               close(FN)
            else	! no such file => create it
8392        open(newunit = FN, FILE = trim(adjustl(file_with_cdf)))
               ! Get the MFP and CDF points:
               GRD:do k = 1, Nsiz
                  E => Element%Phot_absorption%E(k)	! photon energy [eV]
                  if (Element%Phot_absorption%Per_shell(m,k) > 0.0d0) then
                     lambda(k) = 1.0d0/(Element%Phot_absorption%Per_shell(m,k)) ! MFP is normlaized to atomic density of 1 [atom/A^3]
                     ! To convert for particular material: lambda = lambda/(At_density*1.0d24) : where At_density is in [atom/cm^3] and 1.0d24 [cm^3/A^3]
                  else
                     lambda(k) = 1.0d25	! infinity
                  endif
                  ! Loss function is accordingly normalized to the same density: 
                  CDF_data(k) = g_cvel*1.0d10*g_h/(E*g_e*lambda(k))
                  write(FN,'(es,es,es)') E, lambda(k), CDF_data(k)
               enddo GRD
               close(FN)
            endif ! (file_exist)
            
            ! Allocate CDF arrays:
            if (.not.allocated(Element%CDF)) allocate(Element%CDF(Element%N_shl))
            ! Check if the files with fitted CDF coefficients already exist:
            file_with_coefs = trim(adjustl(folder_with_cdf))//path_sep//trim(adjustl(m_fitted_coeffs))//'_'//trim(adjustl(Element%Shell_name(m)))//'.dat'
            inquire(file=trim(adjustl(file_with_coefs)),exist=file_exist) ! check if this file is there
            if (file_exist) then	! just read from the file, no need to recalculate:
               if (.not.Element%valent(m)) then
                  print*, 'Reading CDF coefficients for '//trim(adjustl(Element%Name))//', shell '//trim(adjustl(Element%Shell_name(m)))
               endif
               open(newunit = FN2, FILE = trim(adjustl(file_with_coefs)))	! file with coefficients
               call Count_lines_in_file(FN2, N_CDF, skip_lines=1)	! modul "Dealing_with_files"
               ! Now as we know how many fitting functions are for this shell, allocate the arrays:
               allocate(Element%CDF(m)%A(N_CDF))
               allocate(Element%CDF(m)%E0(N_CDF))
               allocate(Element%CDF(m)%Gamma(N_CDF))
               ! Read the file with CDF coefficients:
               read(FN2,*) ! Skip the descriptor line: 'A	E0	Gamma'
               !do k = 1, N_CDF	! all coefficients
               !   read(FN2,*) Element%CDF(m)%A(k), Element%CDF(m)%E0(k), Element%CDF(m)%Gamma(k)
               !enddo
               call read_CDF_from_file(FN2, Element%CDF(m), N_CDF, read_well)	! module "Read_input_data"
               if (.not.read_well) then
                  print*, ' Could not read CDF coefficients from file '//trim(adjustl(file_with_coefs))
                  close(FN2)	! redo the file
                  goto 8393
               endif
            else	! Fit the oscillator CDF to those data:
8393        print*, 'Coefficients of CDF are going to be fitted - they need CHECKING before use!'
               print*, 'Fitting '//trim(adjustl(Element%Name))//' shell '//trim(adjustl(Element%Shell_name(m)))
               call fit_oscillators_CDF(Element%Phot_absorption%E, CDF_data, Element%CDF(m), Element%Ne_shell(m))		! see below
               ! And save them in the output files:
!                file_with_coefs = trim(adjustl(folder_with_cdf))//path_sep//trim(adjustl(m_fitted_coeffs))//'_'//trim(adjustl(Element%Shell_name(m)))//'.dat'
               open(newunit = FN2, FILE = trim(adjustl(file_with_coefs)))	! file with coefficients
               write(FN2,'(a)') 'A	E0	Gamma'
               do k = 1, size(Element%CDF(m)%A)	! all coefficients
                  write(FN2,'(f,f,f)') Element%CDF(m)%A(k), Element%CDF(m)%E0(k), Element%CDF(m)%Gamma(k)
               enddo
            endif ! (file_exist)
            close(FN2)
            
            ! Renormalize the CDF coefficients to the proper material density:
            Element%CDF(m)%A = Element%CDF(m)%A/1.0d24*used_target%Material(i)%At_Dens
            
            ! Check sum rules of the obtained coefficients
            !Omega =  w_plasma(1d6*1.0d24)	! module "CDF_Ritchi"
            Omega =  w_plasma(1d6*used_target%Material(i)%At_Dens)	! module "CDF_Ritchi"
            
            ! Minimal energy sufficient for ionization is 2*Ip (testing):
!             call sum_rules(Element%CDF(m)%A,  Element%CDF(m)%E0, Element%CDF(m)%Gamma, 2.0d0*Element%Ip(m), Omega, ksum, fsum)	! module "CDF_Ritchi"
            
            call sum_rules(Element%CDF(m)%A,  Element%CDF(m)%E0, Element%CDF(m)%Gamma, Element%Ip(m), Omega, ksum, fsum)	! module "CDF_Ritchi"
            if (.not.Element%valent(m)) then ! only for core orbitals
               write(*,'(a,f10.3,a,f10.3,a,f12.5)') 'K-sum rule:', ksum, ' Ne=', Element%Ne_shell(m), ' F-sum rule:', fsum
               print*, '------------------------'
            endif
            ! Save the coefficient of the electronic plasma frequency:
            if (.not.allocated(Element%CDF(m)%h_omega_e2)) allocate(Element%CDF(m)%h_omega_e2(size(Element%CDF(m)%A)))
            Element%CDF(m)%h_omega_e2 = (g_h/g_e)*(g_h/g_e) * Omega * Element%Ne_shell(m)    ! [eV^2]
            Element%CDF(m)%h_omega_e2 = g_Pi * 0.5d0* Element%CDF(m)%h_omega_e2 ! with coefficient to convert to "alpha"
            
            ! Define MFPs for photoabsorption for each shell:
            if (.not. allocated(Element%Phot_absorption%MFP)) allocate(Element%Phot_absorption%MFP(size(Element%Phot_absorption%Per_shell,1),size(Element%Phot_absorption%Per_shell,2)))
            ! Take into account the density of atoms of this particular kind:
            elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
            do k = 1, Nsiz  ! for all energy grid points:
               Element%Phot_absorption%MFP(m,k) = MFP_from_sigma(Element%Phot_absorption%Per_shell(m,k), used_target%Material(i)%At_Dens)    ! module "CS_general_tools"
               Element%Phot_absorption%MFP(m,k) = Element%Phot_absorption%MFP(m,k) / elem_contrib
            enddo
         enddo SHL
         ! Clear up the arrays to be reused for the next element:
         if (allocated(lambda)) deallocate(lambda)
         if (allocated(CDF_data)) deallocate(CDF_data)
      enddo ELMNT
      
      ! Path to the files with the material parameters:
      Path_valent = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
      Path_valent = trim(adjustl(Path_valent))//path_sep//trim(adjustl(used_target%Material(i)%Name))

      ! Check if valence band needs single-pole approximation:
      if ( (.not.allocated(used_target%Material(i)%CDF_valence%A)) .or. (numpar%El_inelast == 5)) then
         print*, 'Do not have CDF for valence band in '//trim(adjustl(used_target%Material(i)%Name))
         print*, 'Using single-pole approximation instead:'
         ! Allocate the single pole parameters:
         if (allocated(used_target%Material(i)%CDF_valence%A)) deallocate(used_target%Material(i)%CDF_valence%A)
         if (allocated(used_target%Material(i)%CDF_valence%E0)) deallocate(used_target%Material(i)%CDF_valence%E0)
         if (allocated(used_target%Material(i)%CDF_valence%Gamma)) deallocate(used_target%Material(i)%CDF_valence%Gamma)
         if (allocated(used_target%Material(i)%CDF_valence%h_omega_e2)) deallocate(used_target%Material(i)%CDF_valence%h_omega_e2)
         allocate(used_target%Material(i)%CDF_valence%A(1))
         allocate(used_target%Material(i)%CDF_valence%E0(1))
         allocate(used_target%Material(i)%CDF_valence%Gamma(1))
         allocate(used_target%Material(i)%CDF_valence%h_omega_e2(1))
         ! Set them according to the single-pole approximation:
         ! E0 is set at plasma frequency (with electron effective mass):
         !Omega =  w_plasma( 1d6*used_target%Material(i)%At_Dens/dble(SUM(used_target%Material(i)%Elements(:)%percentage)))	! module "CDF_Ritchi"
         Omega = w_plasma( 1d6*used_target%Material(i)%At_Dens/dble(SUM(used_target%Material(i)%Elements(:)%percentage)), &
                            Mass=used_target%Material(i)%me_eff*g_me )  ! module "CDF_Ritchi"
          
         used_target%Material(i)%CDF_valence%E0(1) = sqrt((g_h/g_e)*(g_h/g_e) * Omega * used_target%Material(i)%N_VB_el)  ! [eV]
         ! Gamma set equal to E0 without effective mass (empirical approximation):
         used_target%Material(i)%CDF_valence%Gamma(1) = used_target%Material(i)%CDF_valence%E0(1)
         ! A is set vie normalization (sum rule):
         used_target%Material(i)%CDF_valence%A(1) = 1.0d0   ! just to get sum rule to renormalize below
         ! Get plasmon omega for free-electrons:
         Omega = w_plasma( 1d6*used_target%Material(i)%At_Dens/dble(SUM(used_target%Material(i)%Elements(:)%percentage)) )  ! module "CDF_Ritchi"
         
         call sum_rules(used_target%Material(i)%CDF_valence%A,  used_target%Material(i)%CDF_valence%E0, used_target%Material(i)%CDF_valence%Gamma, used_target%Material(i)%DOS%Egap, Omega, ksum, fsum)	! module "CDF_Ritchi"
         used_target%Material(i)%CDF_valence%A(1) = used_target%Material(i)%N_VB_el/ksum
         ! Parameterse used for delta-CDF:
         used_target%Material(i)%CDF_valence%h_omega_e2 = (g_h/g_e)*(g_h/g_e) * Omega * used_target%Material(i)%N_VB_el ! [eV^2]
         used_target%Material(i)%CDF_valence%h_omega_e2 = g_Pi*0.5d0*used_target%Material(i)%CDF_valence%h_omega_e2 ! with coefficient to convert to "alpha"
         ! Now we can recalculate sum-rules:
!          call sum_rules(used_target%Material(i)%CDF_valence%A,  used_target%Material(i)%CDF_valence%E0, used_target%Material(i)%CDF_valence%Gamma, used_target%Material(i)%DOS%Egap, Omega, ksum, fsum)	! module "CDF_Ritchi"
!          write(*,'(a,f6.3,a,f6.3,a,f6.3)') 'K-sum rule:', ksum, ' Ne=', used_target%Material(i)%N_VB_el, ' F-sum rule:', fsum
!          print*, '------------------------'
      endif ! if (allocated(used_target%Material(i)%CDF_valence%A))
      
      ! Deal with the valence band
      if (allocated(used_target%Material(i)%CDF_valence%A)) then    ! just to check
         ! Save the total valence band cross section:
         File_name = trim(adjustl(Path_valent))//path_sep//trim(adjustl(m_photon_absorption))//'_valence_CDF_Ritchie.dat'
         if(.not.allocated(used_target%Material(i)%Ph_absorption_valent%Total)) then   ! valence band has not been defined yet, so do that
            Nsiz = size(used_target%Material(i)%Ph_absorption_valent%E)
            allocate(used_target%Material(i)%Ph_absorption_valent%Total(Nsiz))
            allocate(used_target%Material(i)%Ph_absorption_valent%Total_MFP(Nsiz))
            ! Check file with CSs:
            inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
            !if (file_exist) then	! just read from the file, no need to recalculate:
            if ( (file_exist) .and. (.not.numpar%recalculate_MFPs) ) then ! just read from the file, no need to recalculate:
               open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
               ! Get the cross section points:
               count_lines = 0
               do m = 1, Nsiz	! for all energy grid points:
                  E => used_target%Material(i)%Ph_absorption_valent%E(m)	! photon energy [eV]
                  read(FN,'(es,es,es)', IOSTAT=Reason) E, used_target%Material(i)%Ph_absorption_valent%Total(m), used_target%Material(i)%Ph_absorption_valent%Total_MFP(m) !, used_target%Material(i)%Ph_absorption_valent%Total_Se(m)
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"
                  if (.not.read_well) then
                  close(FN)	! redo the file
                  goto 8391
                  endif
               enddo
            else	! no such file => create it
8391           open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
               if (.not.allocated(lambda)) allocate(lambda(Nsiz))
               lambda = 1.0d30  ! to start with
               do m = 1, Nsiz	! for all energy grid points:
                  E => used_target%Material(i)%Ph_absorption_valent%E(m)	! photon energy [eV]
                  if (E >= used_target%Material(i)%DOS%Egap) then   ! absorption is possible
                     sigma = get_photoabsorption_from_CDF(used_target%Material(i)%CDF_valence, used_target%Material(i)%At_Dens, E) ! below
                     used_target%Material(i)%Ph_absorption_valent%Total(m) = sigma   ! [A^2]
                  else  ! absorption is not possible
                     used_target%Material(i)%Ph_absorption_valent%Total(m) = 0.0d0   ! [A^2]
                  endif
                  lambda(m) = MFP_from_sigma(sigma, used_target%Material(i)%At_Dens)    ! module "CS_general_tools"
                  used_target%Material(i)%Ph_absorption_valent%Total_MFP(m) = lambda(m)  ! [1/A]
               enddo ! m = 1, Nsiz
               
               ! Combine valence CDF (at low energy) with the atomic one (at the high energy):
               do m = Nsiz, 1, -1   ! scan all energy grid points, starting from the top
                  E => used_target%Material(i)%Ph_absorption_valent%E(m)	! photon energy [eV]
                  ! Only in the region where we need to use the atomic CS:
                  ! Get sigma from the atomic CSs:
                  sigma = 0.0d0 ! to start with
                  do j =1, used_target%Material(i)%N_Elements	! for each element
                     Element => used_target%Material(i)%Elements(j)	! all information about this element
                     elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
                     !sigma = sigma + SUM(Element%Phot_absorption%Per_shell(:,m), MASK = Element%valent)*elem_contrib	! [A^2]
                     ! Get the cross section for this element on this grid point (which may not coincide with the grid for this element due to element-specific grid points):
                     do k = 1, Element%N_shl		! for all shells in this elements
                        if (Element%valent(k)) then
                           ! For each shell, find its partial cross section:
                           call interpolate_data_single(Element%Phot_absorption%E,  Element%Phot_absorption%Per_shell(k,:), E, sigma_cur) ! module "Little_subroutines"
                           ! Add it into the total valent cross section:
                        else    ! exclude non-valent shells
                           sigma_cur = 0.0d0
                        endif
                        sigma = sigma + sigma_cur*elem_contrib
                     enddo ! k
                  enddo ! j
                  ! At sufficiently low energy (100 eV), a peak in valence CDF should cross the atomic one:
                  if ( (E < 100.0d0) .and. (sigma < used_target%Material(i)%Ph_absorption_valent%Total(m))) then ! replace atomic sigma with valent one:
                     exit ! the rest we leave equal to the valence CS precalculated above
                  endif
                  used_target%Material(i)%Ph_absorption_valent%Total(m) = sigma   ! [A^2]
!                   print*, j, E, Element%valent
                  lambda(m) = MFP_from_sigma(sigma, used_target%Material(i)%At_Dens)    ! module "CS_general_tools"
                  used_target%Material(i)%Ph_absorption_valent%Total_MFP(m) = lambda(m)  ! [1/A]
               enddo ! m
               
               ! Write the combined valence CDF into the file:
               do m = 1, Nsiz	! for all energy grid points:
                  write(FN,'(es,es,es)') used_target%Material(i)%Ph_absorption_valent%E(m), used_target%Material(i)%Ph_absorption_valent%Total(m), lambda(m) !, Se
               enddo ! m = 1, Nsiz
               ! Clean up after usage:
               if (allocated(lambda)) deallocate(lambda)
            endif ! (file_exist)
            close(FN)
         endif !  (.not.allocated(used_target%Material(i)%Ph_absorption_valent%E))

         ! Print out k-sum rule:
         print*, 'Valence band in '//trim(adjustl(used_target%Material(i)%Name))
         Omega =  w_plasma( 1d6*used_target%Material(i)%At_Dens/dble(SUM(used_target%Material(i)%Elements(:)%percentage)))	! module "CDF_Ritchi"
         call sum_rules(used_target%Material(i)%CDF_valence%A,  used_target%Material(i)%CDF_valence%E0, used_target%Material(i)%CDF_valence%Gamma, used_target%Material(i)%DOS%Egap, Omega, ksum, fsum)	! module "CDF_Ritchi"
         write(*,'(a,f10.3,a,f10.3,a,f12.5)') 'K-sum rule:', ksum, ' Ne=', used_target%Material(i)%N_VB_el, ' F-sum rule:', fsum
         print*, '------------------------'
         ! Save the coefficient of the electronic plasma frequency:
         if (.not.allocated(used_target%Material(i)%CDF_valence%h_omega_e2)) then
            allocate(used_target%Material(i)%CDF_valence%h_omega_e2(size(used_target%Material(i)%CDF_valence%A)))
         endif
         used_target%Material(i)%CDF_valence%h_omega_e2 = (g_h/g_e)*(g_h/g_e) * Omega * used_target%Material(i)%N_VB_el ! [eV^2]
         used_target%Material(i)%CDF_valence%h_omega_e2 = g_Pi * 0.5d0* used_target%Material(i)%CDF_valence%h_omega_e2 ! with coefficient to convert to "alpha"
!          !------------------------------------------------
!          ! Print out k-sum rule for phonons:
!          print*, 'Phonons in '//trim(adjustl(used_target%Material(i)%Name))
!          Omega = w_plasma( 1d6*used_target%Material(i)%At_Dens/dble(SUM(used_target%Material(i)%Elements(:)%percentage)), &
!                                          Mass=used_target%Material(i)%Mean_Mass )  ! module "CDF_Ritchi"
!          call sum_rules(used_target%Material(i)%CDF_phonon%A,  used_target%Material(i)%CDF_phonon%E0, &
!                                 used_target%Material(i)%CDF_phonon%Gamma, 1.0d-8, Omega, ksum, fsum)	! module "CDF_Ritchi"
!          print*, 'K-sum rule:', ksum, ' Ne=', dble(SUM(used_target%Material(i)%Elements(:)%percentage))
!          print*, '------------------------'
      else  ! if (allocated(used_target%Material(i)%CDF_valence%A)) then
         print*, 'Do not have CDF for valence band in '//trim(adjustl(used_target%Material(i)%Name))
         !print*, 'Can not calculate sum rules...'
         print*, 'Using single-pole approximation instead:'
         print*, '------------------------'
      endif ! if (allocated(used_target%Material(i)%CDF_valence%A)) then
      
      ! Also deal with phonons if possible:
      if (allocated(used_target%Material(i)%CDF_phonon%A) .and. (numpar%El_elastic == 1)) then
         ! Print out k-sum rule for phonons:
         print*, 'Phonons in '//trim(adjustl(used_target%Material(i)%Name))
         Omega = w_plasma( 1d6*used_target%Material(i)%At_Dens/dble(SUM(used_target%Material(i)%Elements(:)%percentage)), &
                                         Mass=used_target%Material(i)%Mean_Mass )  ! module "CDF_Ritchi"
         call sum_rules(used_target%Material(i)%CDF_phonon%A,  used_target%Material(i)%CDF_phonon%E0, &
                                used_target%Material(i)%CDF_phonon%Gamma, 1.0d-8, Omega, ksum, fsum)	! module "CDF_Ritchi"
         write(*,'(a,f10.3,a,f10.3,a,f12.5)') 'K-sum rule:', ksum, ' Na=', dble(SUM(used_target%Material(i)%Elements(:)%percentage)), ' F-sum rule:', fsum
         print*, '------------------------'
      else
         print*, 'Do not have CDF for phonons in '//trim(adjustl(used_target%Material(i)%Name))
         print*, 'Using single-pole approximation instead'
         !print*, 'Can not calculate sum rules...'
         ! Allocate the single pole parameters:
         if (allocated(used_target%Material(i)%CDF_phonon%A)) deallocate(used_target%Material(i)%CDF_phonon%A)
         if (allocated(used_target%Material(i)%CDF_phonon%E0)) deallocate(used_target%Material(i)%CDF_phonon%E0)
         if (allocated(used_target%Material(i)%CDF_phonon%Gamma)) deallocate(used_target%Material(i)%CDF_phonon%Gamma)
         if (allocated(used_target%Material(i)%CDF_phonon%h_omega_e2)) deallocate(used_target%Material(i)%CDF_phonon%h_omega_e2)
         allocate(used_target%Material(i)%CDF_phonon%A(1))
         allocate(used_target%Material(i)%CDF_phonon%E0(1))
         allocate(used_target%Material(i)%CDF_phonon%Gamma(1))
         allocate(used_target%Material(i)%CDF_phonon%h_omega_e2(1))
         ! Set them according to the single-pole approximation:
         ! E0 is set at plasma frequency:
         !used_target%Material(i)%CDF_phonon%E0(1) = sqrt((g_h/g_e)*(g_h/g_e) * Omega)
         !used_target%Material(i)%CDF_phonon%E0(1) = used_target%Material(i)%Edebye
         ! Set it at [2 x Einstein frequency]:
         used_target%Material(i)%CDF_phonon%E0(1) = 2.0d0 * used_target%Material(i)%E_eistein   ! [eV]

         ! Gamma set equal to 0.5*E0 (empirical approximation):
         used_target%Material(i)%CDF_phonon%Gamma(1) = used_target%Material(i)%CDF_phonon%E0(1) * 0.5d0

         ! A is set vie normalization (sum rule):
         used_target%Material(i)%CDF_phonon%A(1) = 1.0d0   ! just to get sum rule to renormalize below
         Omega = w_plasma( 1d6*used_target%Material(i)%At_Dens/dble(SUM(used_target%Material(i)%Elements(:)%percentage)), &
                                         Mass=used_target%Material(i)%Mean_Mass )  ! module "CDF_Ritchi"
         call sum_rules(used_target%Material(i)%CDF_phonon%A,  used_target%Material(i)%CDF_phonon%E0, &
                                used_target%Material(i)%CDF_phonon%Gamma, 1.0d-8, Omega, ksum, fsum)	! module "CDF_Ritchi"
         used_target%Material(i)%CDF_phonon%A(1) = dble(SUM(used_target%Material(i)%Elements(:)%percentage))/ksum
         ! And renormalize per mass:
         !used_target%Material(i)%CDF_phonon%A(1) = used_target%Material(i)%CDF_phonon%A(1) * sqrt(g_me/used_target%Material(i)%Mean_Mass)
         ! Parameterse used for delta-CDF:
         used_target%Material(i)%CDF_phonon%h_omega_e2 = (g_h/g_e)*(g_h/g_e)*Omega*dble(SUM(used_target%Material(i)%Elements(:)%percentage)) ! [eV^2]
         used_target%Material(i)%CDF_phonon%h_omega_e2 = g_Pi*0.5d0*used_target%Material(i)%CDF_phonon%h_omega_e2 ! with coefficient to convert to "alpha"
         ! Now we can recalculate sum-rules:
         call sum_rules(used_target%Material(i)%CDF_phonon%A,  used_target%Material(i)%CDF_phonon%E0, &
                                used_target%Material(i)%CDF_phonon%Gamma, 1.0d-8, Omega, ksum, fsum)	! module "CDF_Ritchi"
         write(*,'(a,f10.3,a,f10.3,a,f12.5)') 'K-sum rule:', ksum, ' Na=', dble(SUM(used_target%Material(i)%Elements(:)%percentage)), ' F-sum rule:', fsum
         print*, '------------------------'
      endif

      ! And the total CS:
      N_elem = dble(SUM(used_target%Material(i)%Elements(:)%percentage))	! number of elements in this compound material
      allocate(used_target%Material(i)%Ph_absorption_total%Total(size(used_target%Material(i)%Ph_absorption_total%E)))
      allocate(used_target%Material(i)%Ph_absorption_total%Total_MFP(size(used_target%Material(i)%Ph_absorption_total%E)))
      used_target%Material(i)%Ph_absorption_total%Total(:) = 0.0d0
      used_target%Material(i)%Ph_absorption_total%Total_MFP(:) = 0.0d0
      ! Sum up CSs from each element:
      do j =1, used_target%Material(i)%N_Elements	! for each element
         Element => used_target%Material(i)%Elements(j)	! all information about this element
         elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         do m = 1, Nsiz	! for all energy grid points:
            E => used_target%Material(i)%Ph_absorption_total%E(m)	! photon energy [eV]
            ! Core shells:
            !used_target%Material(i)%Ph_absorption_total%Total(m) = used_target%Material(i)%Ph_absorption_total%Total(m) + &
            !        SUM(Element%Phot_absorption%Per_shell(:,m), MASK = .not.Element%valent )*elem_contrib	! [A^2]
            ! Get sigma from the atomic CSs:
            sigma = 0.0d0 ! to start with
            ! Get the cross section for this element on this grid point (which may not coincide with the grid for this element due to element-specific grid points):
            do k = 1, Element%N_shl		! for all shells in this elements
               if (.not.Element%valent(k)) then ! core orbitals
                  ! For each shell, find its partial cross section:
                  call interpolate_data_single(Element%Phot_absorption%E,  Element%Phot_absorption%Per_shell(k,:), E, sigma_cur) ! module "Little_subroutines"
               else    ! exclude valent shells
                  sigma_cur = 0.0d0
               endif
               sigma = sigma + sigma_cur*elem_contrib
            enddo ! k
            used_target%Material(i)%Ph_absorption_total%Total(m) = used_target%Material(i)%Ph_absorption_total%Total(m) + sigma
            ! Add valence band / shells:
            if (allocated(used_target%Material(i)%CDF_valence%A)) then    ! valence band is defined and used in the CS model (RBEB excluded)
               if (j == 1) then ! add valence band only once - when studying the first element of the material
                  used_target%Material(i)%Ph_absorption_total%Total(m) = used_target%Material(i)%Ph_absorption_total%Total(m) + used_target%Material(i)%Ph_absorption_valent%Total(m)   ! [A^2]
               endif
            else ! valent atomic levels
               sigma = 0.0d0 ! to start with
               ! Get the cross section for this element on this grid point (which may not coincide with the grid for this element due to element-specific grid points):
               do k = 1, Element%N_shl		! for all shells in this elements
                  if (Element%valent(k)) then ! valence orbitals
                     ! For each shell, find its partial cross section:
                     call interpolate_data_single(Element%Phot_absorption%E,  Element%Phot_absorption%Per_shell(k,:), E, sigma_cur) ! module "Little_subroutines"
                  else    ! exclude non-valent shells
                     sigma_cur = 0.0d0
                  endif
                  sigma = sigma + sigma_cur*elem_contrib
               enddo ! k
               used_target%Material(i)%Ph_absorption_total%Total(m) = used_target%Material(i)%Ph_absorption_total%Total(m) + sigma
            endif ! (allocated(used_target%Material(i)%CDF_valence%A)) 
         enddo ! m = 1, Nsiz
      enddo !  j =1, N_elements
      ! And total MFP for photoabsorption:
      do m = 1, Nsiz	! for all energy grid points:
            ! Get the total MFP from the total CS:
            used_target%Material(i)%Ph_absorption_total%Total_MFP(m) = MFP_from_sigma(used_target%Material(i)%Ph_absorption_total%Total(m), used_target%Material(i)%At_Dens)    ! module "CS_general_tools"
      enddo
      ! Save it into the file:
      File_name = trim(adjustl(Path_valent))//path_sep//trim(adjustl(m_photon_absorption))//'_total_CDF_Ritchie.dat'
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
      if (.not.file_exist) then	! only create it if file does not exist
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            ! Get the total MFP from the total CS:
            write(FN,'(es,es,es)') used_target%Material(i)%Ph_absorption_total%E(m), used_target%Material(i)%Ph_absorption_total%Total(m), &
                                                  used_target%Material(i)%Ph_absorption_total%Total_MFP(m) !, used_target%Material(i)%Ph_absorption_total%Total_Se(m)
         enddo ! m = 1, Nsiz
         close(FN)
      endif
      
      
      ! If delta-model of CDF is used (which is always used for SHIs), get alpha coefficients:
      ! For all shells:
      do j = 1, used_target%Material(i)%N_Elements	! for all elements from this target consituent
         Element => used_target%Material(i)%Elements(j)	! all information about this element
         ! Get MFPs and CDFs for each shell:
         do m = 1, Element%N_shl		! for all shells in this elements
            if (.not.allocated(Element%CDF(m)%alpha)) then
               allocate(Element%CDF(m)%alpha(size(Element%CDF(m)%A)))
            endif
            !if (.not.allocated(Element%CDF(m)%alpha_SHI)) then
            !   allocate(Element%CDF(m)%alpha_SHI(size(Element%CDF(m)%A)))
            !endif
            if (.not.Element%valent(m)) then ! only for core orbitals
               print*, trim(adjustl(Element%Name)), ' '//trim(adjustl(Element%Shell_name(m))), ' CDF parameters for delta-function:'
            endif
            do k = 1, size(Element%CDF(m)%A)   ! for all oscillators describing this shell's CDF:
!                 call define_alpha(Element%CDF(m)%A(k), Element%CDF(m)%Gamma(k), Element%CDF(m)%E0(k), Element%Ip(m), Element%CDF(m)%alpha(k))  ! module "CDF_delta"
!                print*, k, 'NEW  E0:', Element%CDF(m)%E0(k)
!                print*, k, 'NEW alpha:', Element%CDF(m)%alpha(k)
               
!                call define_alpha_OLD(Element%CDF(m)%A(k), Element%CDF(m)%Gamma(k), Element%CDF(m)%E0(k), Element%Ip(m), Element%CDF(m)%alpha(k))  ! module "CDF_delta"
               ! Account for the fact that integration limit in the cross section starts from Wmin=c*Ip:
               call define_alpha_OLD(Element%CDF(m)%A(k), Element%CDF(m)%Gamma(k), Element%CDF(m)%E0(k), Element%CDF(m)%E0(k), Element%CDF(m)%alpha(k))  ! module "CDF_delta"

               if (.not.Element%valent(m)) then ! only for core orbitals
                  print*, k, 'E0:', Element%CDF(m)%E0(k)
                  print*, k, 'alpha:', Element%CDF(m)%alpha(k)
               endif
               !! Save for SHI:
               !Element%CDF(m)%alpha_SHI(k) = Element%CDF(m)%alpha(k)

            enddo
            if (.not.Element%valent(m)) print*, '------------------------'
         enddo ! m = 1, Element%N_shl
      enddo ! j = 1, used_target%Material(i)%N_Elements
      ! And also valence band, if present:
      if (allocated(used_target%Material(i)%CDF_valence%A)) then
         if (.not.allocated(used_target%Material(i)%CDF_valence%alpha)) then
            allocate(used_target%Material(i)%CDF_valence%alpha(size(used_target%Material(i)%CDF_valence%A)))
         endif
         !if (.not.allocated(used_target%Material(i)%CDF_valence%alpha_SHI)) then
         !   allocate(used_target%Material(i)%CDF_valence%alpha_SHI(size(used_target%Material(i)%CDF_valence%A)))
         !endif
         print*, 'Valence band CDF parameters for delta-function:'
         do k = 1, size(used_target%Material(i)%CDF_valence%A)   ! for all oscillators describing valence CDF
!             call define_alpha(used_target%Material(i)%CDF_valence%A(k), used_target%Material(i)%CDF_valence%Gamma(k), used_target%Material(i)%CDF_valence%E0(k), used_target%Material(i)%DOS%Egap, used_target%Material(i)%CDF_valence%alpha(k))  ! module "CDF_delta"
!             print*, k, 'E0:', used_target%Material(i)%CDF_valence%E0(k)
!             print*, k, 'alpha:', used_target%Material(i)%CDF_valence%alpha(k)
            
               call define_alpha_OLD(used_target%Material(i)%CDF_valence%A(k), used_target%Material(i)%CDF_valence%Gamma(k), &
                      used_target%Material(i)%CDF_valence%E0(k), used_target%Material(i)%DOS%Egap, used_target%Material(i)%CDF_valence%alpha(k))  ! module "CDF_delta"
               ! Normalize it to the molecular density:
!                ! Account for the fact that integration limit in the cross section starts from Wmin=2*Ip:
!                call define_alpha_OLD(used_target%Material(i)%CDF_valence%A(k), used_target%Material(i)%CDF_valence%Gamma(k), &
!                 used_target%Material(i)%CDF_valence%E0(k), 1.5d0*used_target%Material(i)%DOS%Egap, used_target%Material(i)%CDF_valence%alpha(k))  ! module "CDF_delta"
               print*, k, 'E0:', used_target%Material(i)%CDF_valence%E0(k)
               print*, k, 'alpha:', used_target%Material(i)%CDF_valence%alpha(k)
               !! Save for SHI:
               !used_target%Material(i)%CDF_valence%alpha_SHI(k) = used_target%Material(i)%CDF_valence%alpha(k)

         enddo
         print*, '------------------------'
      endif
      ! And also phonons, if present:
      if (allocated(used_target%Material(i)%CDF_phonon%A)) then
         if (.not.allocated(used_target%Material(i)%CDF_phonon%alpha)) then
            allocate(used_target%Material(i)%CDF_phonon%alpha(size(used_target%Material(i)%CDF_phonon%A)))
         endif
         !if (.not.allocated(used_target%Material(i)%CDF_phonon%alpha_SHI)) then
         !   allocate(used_target%Material(i)%CDF_phonon%alpha_SHI(size(used_target%Material(i)%CDF_phonon%A)))
         !endif

         print*, 'Phonons CDF parameters for delta-function:'
         do k = 1, size(used_target%Material(i)%CDF_phonon%A)   ! for all oscillators describing valence CDF
            ! For scattering on atoms, also account for the mass factor in the sum rule:
!             call define_alpha(used_target%Material(i)%CDF_phonon%A(k), used_target%Material(i)%CDF_phonon%Gamma(k), &
!                                       used_target%Material(i)%CDF_phonon%E0(k), used_target%Material(i)%DOS%Egap, &
!                                       used_target%Material(i)%CDF_phonon%alpha(k), mass_target = used_target%Material(i)%Mean_Mass)  ! module "CDF_delta"
!             print*, k, 'E0:', used_target%Material(i)%CDF_phonon%E0(k), used_target%Material(i)%CDF_phonon%A(k)
!             print*, k, 'alpha:', used_target%Material(i)%CDF_phonon%alpha(k) !, used_target%Material(i)%Mean_Mass

            call define_alpha_OLD(used_target%Material(i)%CDF_phonon%A(k), used_target%Material(i)%CDF_phonon%Gamma(k), &
                                      used_target%Material(i)%CDF_phonon%E0(k), 1.0d-10, &
!                                       used_target%Material(i)%CDF_phonon%E0(k), minval(used_target%Material(i)%CDF_phonon%E0(:)), &
                                      used_target%Material(i)%CDF_phonon%alpha(k), mass_target = used_target%Material(i)%Mean_Mass)  ! module "CDF_delta"
            print*, k, 'E0,A:', used_target%Material(i)%CDF_phonon%E0(k), used_target%Material(i)%CDF_phonon%A(k)
            print*, k, 'alpha:', used_target%Material(i)%CDF_phonon%alpha(k) !, used_target%Material(i)%Mean_Mass
            !! Save for SHI:
            !used_target%Material(i)%CDF_phonon%alpha_SHI(k) = used_target%Material(i)%CDF_phonon%alpha(k)
            
         enddo
         ! Save the coefficient of the atomic plasma frequency (phonons):
!          Omega =  w_plasma( 1d6*used_target%Material(i)%At_Dens/dble(SUM(used_target%Material(i)%Elements(:)%percentage)))	! module "CDF_Ritchi"
!          if (.not.allocated(used_target%Material(i)%CDF_phonon%h_omega_e2)) then
!             allocate(used_target%Material(i)%CDF_phonon%h_omega_e2(size(used_target%Material(i)%CDF_phonon%A)))
!          endif
!          used_target%Material(i)%CDF_phonon%h_omega_e2 = (g_h/g_e)*(g_h/g_e) * Omega * used_target%Material(i)%N_VB_el ! [eV^2]
!          used_target%Material(i)%CDF_phonon%h_omega_e2 = g_Pi * 0.5d0* used_target%Material(i)%CDF_phonon%h_omega_e2 ! with coefficient to convert to "alpha"
         print*, '------------------------'
      endif
      
   enddo TRGT

! clean up the memory:
   if (allocated(lambda)) deallocate(lambda)
   if (allocated(CDF_data)) deallocate(CDF_data)
   nullify(path_sep, Element, E)
end subroutine get_CDF



 function get_photoabsorption_from_CDF(CDF, nat, E) result (sigma)
   real(8) sigma    ! [A^2] cross section of photoabsorption
   type(Ritchi_CDF), intent(in) :: CDF	! complex dielectric function parameters for each shell
   real(8), intent(in) :: nat   ! [1/cm^3] atomic density
   real(8), intent(in) :: E ! [eV] photon energy
   real(8) :: Im, coef
   ! Ger CDF function:
   Im = Diel_func(CDF%A, CDF%E0, CDF%Gamma, E, 0.0d0, g_me)     ! module "CDF_Ritchi"
   coef = g_e*1.0d14 ! g_e from [eV] -> [J]; 1/10^2 from  [m/c] -> [cm/s];  10^16 from [cm^2] -> [A^2]
   sigma = Im*E/(nat*g_h*g_cvel)*coef   ! [A^2]
end function get_photoabsorption_from_CDF




subroutine fit_oscillators_CDF(E, CDF_data, CDF, Ne)
   real(8), dimension(:), intent(in) :: E, CDF_data	! [eV] energy, CDF data points
   type(Ritchi_CDF), intent(inout) :: CDF	! complex dielectric function parameters for each shell
   real(8), intent(in) :: Ne	! number of electrons per this shell, to check sum rules
   !----------------------------
   type(Ritchi_CDF) :: CDF_fit	! complex dielectric function parameters
   type(Ritchi_CDF) :: CDF_temp	! complex dielectric function parameters
   type(Ritchi_CDF), dimension(:), allocatable :: der_CDF	! derivatives of CDF by parameters
   type(Ritchi_CDF) :: lambda
   real(8), dimension(size(CDF_data)) :: temp_CDF_data	! CDF data points to work with
   real(8) :: MSD, eps, MSD0, estimator, Omega, ksum, fsum, inv_CDF
   integer :: INFO, i, Nfit, Niter, N_data_points, j_count, j, i_nonzero, i_start, Nmin, N_iter_max, N_interval, N_start
   real(8), dimension(:), allocatable :: f_cdf	! calculated values of CDF fitted function on the grid of E
   integer, dimension(:), allocatable :: peaks_pos	! locations of peaks in CDF optical data
   integer, dimension(:), allocatable :: minima_pos	! locations of minima in CDF optical data
   real(8), dimension(size(E)) :: ri	! residuals
   real(8), dimension(size(E)) :: weights	! wieghts of data points
   logical :: one_more_try
   !----------------------------
   one_more_try = .true.
   eps = 1.0d-8	! precision of mean square deviation (MSD) required
   weights(:) = E(:)**3		! to fit the tails, make their weights larger
!    weights(:) = 1.0d0	! to fit the tails, make their weights larger
   
   ! To use for fitting:
!    log_CDF_data = log(CDF_data)
   
   N_iter_max = INT(1e5)
   N_data_points = size(CDF_data)
9999   allocate(f_cdf(N_data_points))
   temp_CDF_data = CDF_data	! to start working with
   !----------------------------
   ! Use non-linear Gauss-Newton method to fit the data points with given functions. It proceeds in a few steps:
   ! 1) Estimate, how may oscillator functions are required to fit the data. Assume it is equal to number fo peaks in data + 1.
   call find_number_of_peaks(CDF_data, Nfit, peaks_pos)	! function below
   ! Locations of minima on the right from the first peak
   call find_number_of_minima(CDF_data, peaks_pos(1), Nmin, minima_pos)	! function below
   
   ! Allocate arrays and objects used for fitting:
   if (.not.allocated(CDF_fit%A)) allocate(CDF_fit%A(Nfit))
   if (.not.allocated(CDF_fit%Gamma)) allocate(CDF_fit%Gamma(Nfit))
   if (.not.allocated(CDF_fit%E0)) allocate(CDF_fit%E0(Nfit))
   allocate(CDF_temp%A(Nfit))
   allocate(CDF_temp%Gamma(Nfit))
   allocate(CDF_temp%E0(Nfit))
   allocate(der_CDF(N_data_points))
   do i = 1, N_data_points
      allocate(der_CDF(i)%A(Nfit))
      allocate(der_CDF(i)%Gamma(Nfit))
      allocate(der_CDF(i)%E0(Nfit))   
   enddo
   allocate(CDF%A(Nfit))
   allocate(CDF%Gamma(Nfit))
   allocate(CDF%E0(Nfit))
   allocate(lambda%A(Nfit))
   allocate(lambda%Gamma(Nfit))
   allocate(lambda%E0(Nfit))
   !----------------------------
   ! 2) Estimate initial guess for the parameters A, E, gamma, according to some assumptions.
   do i = 1, Nfit
!       if (one_more_try) then
       if (i <= size(peaks_pos)) then
         CDF_fit%E0(i) = E(peaks_pos(i))	! E0 is expected to be around a maxima
         CDF_fit%Gamma(i) = E(peaks_pos(i))	! Set starting gamma as equal to E0 (just to try it out at the first iteration)
         CDF_fit%A(i) =  estimate_A(CDF_data(peaks_pos(i)), E(peaks_pos(i)), CDF_fit%E0(i), CDF_fit%Gamma(i))	! module "CDF_Ritchi"
       else	! special one for the tail:
         CDF_fit%E0(i) = CDF_fit%E0(i-1)*3.0d0	! E0 is expected to be around a maxima
         CDF_fit%Gamma(i) = CDF_fit%E0(i)	! Set starting gamma as equal to E0 (just to try it out at the first iteration)
         CDF_fit%A(i) = -CDF_fit%A(i-1)/2.0d0 !/(CDF_fit%E0(i) - CDF_fit%E0(i-1))
       endif
!       endif
      write(*,'(a,i3,f,f,f)') 'CDF0: ', i, CDF_fit%A(i), CDF_fit%E0(i), CDF_fit%Gamma(i)
   enddo
   
   do i = 1, Nfit
      lambda%A(i) = transfer(minval(ABS(CDF_fit%A(:))),1.0d0)
      if (lambda%A(i) < 1.0d-8) lambda%A(i) = 1.0d-1	! just in case if it's excidentally zero
      lambda%E0(i) = 0.05d0*transfer(minval(ABS(CDF_fit%E0(:))),1.0d0)
      if (lambda%E0(i) < 1.0d-8) lambda%E0(i) = 1.0d-1	! just in case if it's excidentally zero
      lambda%Gamma(i) = 0.1d0*transfer(minval(ABS(CDF_fit%Gamma(:))),1.0d0)
      if (lambda%Gamma(i) < 1.0d-8) lambda%Gamma(i) = 1.0d-1	! just in case if it's excidentally zero
   enddo
   ! In case convergence was not reached from the first time, make steps smaller
   if (.not.one_more_try) then	! if it's not the first time, make steps 10 times finer:
      lambda%A(:) = 0.1d0*lambda%A(:)
      lambda%E0(:) = 0.1d0*lambda%E0(:)
      lambda%Gamma(:) = 0.1d0*lambda%Gamma(:)
   endif
   
   ! Evaluate which part of the data makes sense to consider (significantly large):
   i_start = 1
   do while (CDF_data(i_start) < 1.0d-16)	! find the first meaningful point
      i_start = i_start  + 1
   enddo
   i_nonzero = peaks_pos(1)
   do while (CDF_data(i_nonzero) > CDF_data(peaks_pos(1))*1.0d-12)	! find the last meaningful point
      i_nonzero = i_nonzero + 1
   enddo
!    print*, 'i_start', i_start
!    eps = eps*(CDF_data(peaks_pos(1)))**2
   
   !--------------------------------------
   ! Split the data into intervals according to peaks inside:
   N_interval = 0	! start from the first one
   do N_interval = 1, size(minima_pos)	! for all intervals
    ! Next interval should know that the previous one was already fitted:
    if (N_interval > 1) then
       do i = 1, N_data_points
         f_cdf(i) = Diel_func(CDF_fit%A(N_interval-1:N_interval-1), CDF_fit%E0(N_interval-1:N_interval-1), CDF_fit%Gamma(N_interval-1:N_interval-1), E(i), 0.0d0, g_me)	! module "CDF_Ritchi"
       enddo
       temp_CDF_data(:) = residuals(temp_CDF_data, f_cdf)	! new function with fitted function subtracted
    endif
    Niter = 0
    estimator = 1.0d0	! starting point
    MSD0 = 1.0d10
    eps = eps*(CDF_data(peaks_pos(min(N_interval,size(peaks_pos)))))
    
    do while ( (estimator > eps) .or. (Niter < 1000) ) ! do as long as the fitting quality is reached
      ! Set the values of corresponding cdf function:
      do i = 1, N_data_points
         f_cdf(i) = Diel_func(CDF_fit%A(N_interval:N_interval), CDF_fit%E0(N_interval:N_interval), CDF_fit%Gamma(N_interval:N_interval), E(i), 0.0d0, g_me)	! module "CDF_Ritchi"
      enddo
      
      !----------------------------
      ! 3) Evaluate mean square deviation (MSD) of this function from the data?
      ri(:) = residuals(temp_CDF_data, f_cdf)	! see below
!       ! Fit log of the function, to better fit the tail:
!       ri(:) = residuals(log_CDF_data, f_cdf)	! see below
!       if (N_interval == 1) then
         N_start = i_start
!       else
!          N_start = minima_pos(N_interval-1)
!       endif
      MSD = SUM(ri(N_start:minima_pos(N_interval))*ri(N_start:minima_pos(N_interval))*weights(N_start:minima_pos(N_interval)))	! mean square deviation as residuals squared
      
      estimator = ABS((MSD - MSD0)/MSD0)
      if (estimator <= eps) then	! make at least 1000 iterations
          if (Niter > 1000) then
             print*, 'Exiting by the condition of (estimator <= eps)'
             exit	! close enough, consider it done
          endif
       endif
      MSD0 = MSD	! save for the next iteration to compare with
!       print*, 'N_start', N_start, minima_pos(N_interval), MSD, MSD0
      !----------------------------
      ! 4) Evaluate derivatives by parameters:
      do i = N_start, minima_pos(N_interval)
!          do j = 1, Nfit
            j = min(N_interval,size(CDF_fit%A))
            der_CDF(i)%A(j) = dCDF_dA(CDF_fit%A(j), CDF_fit%E0(j), CDF_fit%Gamma(j), E(i))		! module "CDF_Ritchi"
            der_CDF(i)%E0(j) = dCDF_dE0(CDF_fit%A(j), CDF_fit%E0(j), CDF_fit%Gamma(j), E(i))		! module "CDF_Ritchi"
            der_CDF(i)%Gamma(j) = dCDF_dGamma(CDF_fit%A(j), CDF_fit%E0(j), CDF_fit%Gamma(j), E(i))	! module "CDF_Ritchi"
!          enddo
      enddo
     
      !----------------------------
      ! 5) Evaluate next iteration of parameters by steepest discent method:
      CDF_temp%A = 0.0d0
      CDF_temp%E0 = 0.0d0
      CDF_temp%Gamma = 0.0d0
!       do j = 1, Nfit
!       j = N_interval
         j = min(N_interval,size(CDF_fit%A))
         do i = N_start, minima_pos(N_interval)
             CDF_temp%A(j) = CDF_temp%A(j) + der_CDF(i)%A(j) * ri(i)
             CDF_temp%E0(j) = CDF_temp%E0(j) + der_CDF(i)%E0(j) * ri(i)
             CDF_temp%Gamma(j) = CDF_temp%Gamma(j) + der_CDF(i)%Gamma(j) * ri(i)
         enddo
         CDF_fit%A(j) = CDF_fit%A(j) + lambda%A(j) * CDF_temp%A(j)
         CDF_fit%E0(j) = CDF_fit%E0(j) + lambda%E0(j) * CDF_temp%E0(j)
         CDF_fit%Gamma(j) = CDF_fit%Gamma(j) + lambda%Gamma(j) * CDF_temp%Gamma(j)
         if (CDF_fit%Gamma(j) < 0.0d0) then	! transfere the sign from Gamma to A
            CDF_fit%Gamma(j) = -CDF_fit%Gamma(j)
            CDF_fit%A(j) = -CDF_fit%A(j)
         endif
!       enddo
      
      Niter = Niter + 1
      if (Niter > N_iter_max) then
         print*, 'Exit due to #iterations:', Niter
         exit	! too many iterations, time to stop this
      endif
    enddo	! while (estimator > eps)
   enddo	! N_interval
   !----------------------------
   ! Check sum rules of the obtained coefficients
   Omega =  w_plasma(1d6*1.0d24)	! module "CDF_Ritchi"
   call sum_rules(CDF_fit%A,  CDF_fit%E0, CDF_fit%Gamma, E(peaks_pos(1)), Omega, ksum, fsum)	! module "CDF_Ritchi"
   write(*,'(a,f6.3,a,f6.3,a,f6.3)') 'Sum rules:', ksum, ' Ne=', Ne, ' F-sum rule:', fsum
   print*, 'iterations=', Niter
   print*, 'estimator=', estimator
   print*, '------------------------'
   
   do i = 1, Nfit
       write(*,'(a,i3,f,f,f)') 'COEF: ', i, CDF_fit%A(i), CDF_fit%E0(i), CDF_fit%Gamma(i)
   enddo
   
   ! Negative sum rule indicate that convergence was not achieved, redo the fitting:
   if ((ksum <= 0.0d0) .or. ( ABS(ksum-Ne)/Ne > 0.5d0 )) then	! deviation is more then 50%
      if (one_more_try) then
!          if (allocated(CDF_fit%A)) deallocate(CDF_fit%A)
!          if (allocated(CDF_fit%Gamma)) deallocate(CDF_fit%Gamma)
!          if (allocated(CDF_fit%E0)) deallocate(CDF_fit%E0)
         if (allocated(CDF_temp%A)) deallocate(CDF_temp%A)
         if (allocated(CDF_temp%Gamma)) deallocate(CDF_temp%Gamma)
         if (allocated(CDF_temp%E0)) deallocate(CDF_temp%E0)
         if (allocated(der_CDF)) deallocate(der_CDF)
         if (allocated(CDF%A)) deallocate(CDF%A)
         if (allocated(CDF%Gamma)) deallocate(CDF%Gamma)
         if (allocated(CDF%E0)) deallocate(CDF%E0)
         if (allocated(lambda%A)) deallocate(lambda%A)
         if (allocated(lambda%E0)) deallocate(lambda%E0)
         if (allocated(lambda%Gamma)) deallocate(lambda%Gamma)
         if (allocated(f_cdf)) deallocate(f_cdf)
         if (allocated(peaks_pos)) deallocate(peaks_pos)
         one_more_try = .false.
         N_iter_max = N_iter_max*2	! try longer
         print*, 'LET ME TRY ONE MORE TIME...'
         goto 9999 ! do over
      endif
   endif
!    PAUSE 'CHECK'
   
   CDF = CDF_fit	! save for output
   !----------------------------
   ! Clean up at the end:
   if (allocated(f_cdf)) deallocate(f_cdf)
   if (allocated(peaks_pos)) deallocate(peaks_pos)
   if (allocated(CDF_fit%A)) deallocate(CDF_fit%A)
   if (allocated(CDF_fit%Gamma)) deallocate(CDF_fit%Gamma)
   if (allocated(CDF_fit%E0)) deallocate(CDF_fit%E0)
   if (allocated(CDF_temp%A)) deallocate(CDF_temp%A)
   if (allocated(CDF_temp%Gamma)) deallocate(CDF_temp%Gamma)
   if (allocated(CDF_temp%E0)) deallocate(CDF_temp%E0)
   if (allocated(der_CDF)) deallocate(der_CDF)
   if (allocated(lambda%A)) deallocate(lambda%A)
   if (allocated(lambda%E0)) deallocate(lambda%E0)
   if (allocated(lambda%Gamma)) deallocate(lambda%Gamma)
end subroutine fit_oscillators_CDF






subroutine fit_oscillators_CDF_OLD(E, CDF_data, CDF, Ne)
   real(8), dimension(:), intent(in) :: E, CDF_data	! [eV] energy, CDF data points
   type(Ritchi_CDF), intent(inout) :: CDF	! complex dielectric function parameters for each shell
   real(8), intent(in) :: Ne	! number of electrons per this shell, to check sum rules
   !----------------------------
   type(Ritchi_CDF) :: CDF_fit	! complex dielectric function parameters
   type(Ritchi_CDF) :: CDF_temp	! complex dielectric function parameters
   type(Ritchi_CDF), dimension(:), allocatable :: der_CDF	! derivatives of CDF by parameters
   type(Ritchi_CDF) :: lambda
   real(8), dimension(size(CDF_data)) :: log_CDF_data	! log of CDF data points
   real(8) :: MSD, eps, MSD0, estimator, Omega, ksum, fsum, inv_CDF
   integer :: INFO, i, Nfit, Niter, N_data_points, j_count, j, i_nonzero, i_start, Nmin, N_iter_max
   real(8), dimension(:), allocatable :: f_cdf	! calculated values of CDF fitted function on the grid of E
   integer, dimension(:), allocatable :: peaks_pos	! locations of peaks in CDF optical data
   integer, dimension(:), allocatable :: minima_pos	! locations of minima in CDF optical data
   real(8), dimension(size(E)) :: ri	! residuals
   real(8), dimension(size(E)) :: weights	! wieghts of data points
   logical :: one_more_try
   !----------------------------
   one_more_try = .true.
   eps = 1.0d-8	! precision of mean square deviation (MSD) required
   weights(:) = exp(sqrt(E(:)))	! to fit the tails, make their weights larger
!    weights(:) = 1.0d0	! to fit the tails, make their weights larger
   
   ! To use for fitting:
!    log_CDF_data = log(CDF_data)
   
   N_iter_max = INT(1e6)
   N_data_points = size(CDF_data)
9999   allocate(f_cdf(N_data_points))
   !----------------------------
   ! Use non-linear Gauss-Newton method to fit the data points with given functions. It proceeds in a few steps:
   ! 1) Estimate, how may oscillator functions are required to fit the data. Assume it is equal to number fo peaks in data + 1.
   call find_number_of_peaks(CDF_data, Nfit, peaks_pos)	! function below
   ! Locations of minima on the right from the first peak
   call find_number_of_minima(CDF_data, peaks_pos(1), Nmin, minima_pos)	! function below
   
   ! Allocate arrays and objects used for fitting:
   allocate(CDF_fit%A(Nfit))
   allocate(CDF_fit%Gamma(Nfit))
   allocate(CDF_fit%E0(Nfit))
   allocate(CDF_temp%A(Nfit))
   allocate(CDF_temp%Gamma(Nfit))
   allocate(CDF_temp%E0(Nfit))
   allocate(der_CDF(N_data_points))
   do i = 1, N_data_points
      allocate(der_CDF(i)%A(Nfit))
      allocate(der_CDF(i)%Gamma(Nfit))
      allocate(der_CDF(i)%E0(Nfit))   
   enddo
   allocate(CDF%A(Nfit))
   allocate(CDF%Gamma(Nfit))
   allocate(CDF%E0(Nfit))
   allocate(lambda%A(Nfit))
   allocate(lambda%Gamma(Nfit))
   allocate(lambda%E0(Nfit))
   !----------------------------
   ! 2) Estimate initial guess for the parameters A, E, gamma, according to some assumptions.
   do i = 1, Nfit
      if (i <= size(peaks_pos)) then
         CDF_fit%E0(i) = E(peaks_pos(i))	! E0 is expected to be around a maxima
         CDF_fit%Gamma(i) = E(peaks_pos(i))	! Set starting gamma as equal to E0 (just to try it out at the first iteration)
         CDF_fit%A(i) =  estimate_A(CDF_data(peaks_pos(i)), E(peaks_pos(i)), CDF_fit%E0(i), CDF_fit%Gamma(i))	! module "CDF_Ritchi"
      else	! special one for the tail:
         CDF_fit%E0(i) = CDF_fit%E0(i-1)*3.0d0	! E0 is expected to be around a maxima
         CDF_fit%Gamma(i) = CDF_fit%E0(i)	! Set starting gamma as equal to E0 (just to try it out at the first iteration)
         CDF_fit%A(i) = -CDF_fit%A(i-1)/2.0d0 !/(CDF_fit%E0(i) - CDF_fit%E0(i-1))
      endif
      write(*,'(a,i3,f,f,f)') 'CDF0: ', i, CDF_fit%A(i), CDF_fit%E0(i), CDF_fit%Gamma(i)
   enddo
   
   do i = 1, Nfit
      lambda%A(i) = transfer(minval(ABS(CDF_fit%A(:))),1.0d0)
      if (lambda%A(i) < 1.0d-8) lambda%A(i) = 1.0d-1	! just in case if it's excidentally zero
      lambda%E0(i) = 0.05d0*transfer(minval(ABS(CDF_fit%E0(:))),1.0d0)
      if (lambda%E0(i) < 1.0d-8) lambda%E0(i) = 1.0d-1	! just in case if it's excidentally zero
      lambda%Gamma(i) = 0.1d0*transfer(minval(ABS(CDF_fit%Gamma(:))),1.0d0)
      if (lambda%Gamma(i) < 1.0d-8) lambda%Gamma(i) = 1.0d-1	! just in case if it's excidentally zero
   enddo
   ! In case convergence was not reached from the first time, make steps smaller
   if (.not.one_more_try) then	! if it's not the first time, make steps 10 times finer:
      lambda%A(:) = 0.1d0*lambda%A(:)
      lambda%E0(:) = 0.1d0*lambda%E0(:)
      lambda%Gamma(:) = 0.1d0*lambda%Gamma(:)
   endif
   
   ! Evaluate which part of the data makes sense to consider (significantly large):
   i_start = 1
   do while (CDF_data(i_start) < 1.0d-16)	! find the first meaningful point
      i_start = i_start  + 1
   enddo
   i_nonzero = peaks_pos(1)
   do while (CDF_data(i_nonzero) > CDF_data(peaks_pos(1))*1.0d-12)	! find the last meaningful point
      i_nonzero = i_nonzero + 1
   enddo
   
!    eps = eps*(CDF_data(peaks_pos(1)))**2
   eps = eps*(CDF_data(peaks_pos(1)))
   
   Niter = 0
   estimator = 1.0d0	! starting point
   MSD0 = 1.0d10
   do while ( (estimator > eps) .or. (Niter < 1000) ) ! do as long as the fitting quality is reached
      ! Set the values of corresponding cdf function:
      do i = 1, N_data_points
         f_cdf(i) = Diel_func(CDF_fit%A, CDF_fit%E0, CDF_fit%Gamma, E(i), 0.0d0, g_me)	! module "CDF_Ritchi"
      enddo
      ! Fit log of the function to better reproduce tails:
!       f_cdf = log(f_cdf)
      
      !----------------------------
      ! 3) Evaluate mean square deviation (MSD) of this function from the data?
      ri(:) = residuals(CDF_data, f_cdf)	! see below
!       ! Fit log of the function, to better fit the tail:
!       ri(:) = residuals(log_CDF_data, f_cdf)	! see below
      MSD = SUM(ri(i_start:i_nonzero)*ri(i_start:i_nonzero)*weights(i_start:i_nonzero))	! mean square deviation as residuals squared
      estimator = ABS((MSD - MSD0)/MSD0)
      if (estimator <= eps) then	! make at least 1000 iterations
          if (Niter > 1000) then
             print*, 'Exiting by the condition of (estimator <= eps)'
             exit	! close enough, consider it done
          endif
       endif
      MSD0 = MSD	! save for the next iteration to compare with
      !----------------------------
      ! 4) Evaluate derivatives by parameters:
      do i = 1, N_data_points
         do j = 1, Nfit
            der_CDF(i)%A(j) = dCDF_dA(CDF_fit%A(j), CDF_fit%E0(j), CDF_fit%Gamma(j), E(i))		! module "CDF_Ritchi"
            der_CDF(i)%E0(j) = dCDF_dE0(CDF_fit%A(j), CDF_fit%E0(j), CDF_fit%Gamma(j), E(i))		! module "CDF_Ritchi"
            der_CDF(i)%Gamma(j) = dCDF_dGamma(CDF_fit%A(j), CDF_fit%E0(j), CDF_fit%Gamma(j), E(i))	! module "CDF_Ritchi"
         enddo
!          ! For fitting log of the function:
!          inv_CDF = Diel_func(CDF_fit%A, CDF_fit%E0, CDF_fit%Gamma, E(i), 0.0d0)		! module "CDF_Ritchi"
!          der_CDF(i)%A(:) = der_CDF(i)%A(:)/inv_CDF
!          der_CDF(i)%E0(:) = der_CDF(i)%E0(:)/inv_CDF
!          der_CDF(i)%Gamma(:) = der_CDF(i)%Gamma(:)/inv_CDF
      enddo
     
      !----------------------------
      ! 5) Evaluate next iteration of parameters by steepest discent method:
      CDF_temp%A = 0.0d0
      CDF_temp%E0 = 0.0d0
      CDF_temp%Gamma = 0.0d0
      do j = 1, Nfit
         do i = i_start, i_nonzero
             CDF_temp%A(j) = CDF_temp%A(j) + der_CDF(i)%A(j) * ri(i)
             CDF_temp%E0(j) = CDF_temp%E0(j) + der_CDF(i)%E0(j) * ri(i)
             CDF_temp%Gamma(j) = CDF_temp%Gamma(j) + der_CDF(i)%Gamma(j) * ri(i)
         enddo
         CDF_fit%A(j) = CDF_fit%A(j) + lambda%A(j) * CDF_temp%A(j)
         CDF_fit%E0(j) = CDF_fit%E0(j) + lambda%E0(j) * CDF_temp%E0(j)
         CDF_fit%Gamma(j) = CDF_fit%Gamma(j) + lambda%Gamma(j) * CDF_temp%Gamma(j)
         if (CDF_fit%Gamma(j) < 0.0d0) then	! transfere the sign from Gamma to A
            CDF_fit%Gamma(j) = -CDF_fit%Gamma(j)
            CDF_fit%A(j) = -CDF_fit%A(j)
         endif
      enddo
      
      Niter = Niter + 1
      if (Niter > N_iter_max) then
         print*, 'Exit due to #iterations:', Niter
         exit	! too many iterations, time to stop this
      endif
   enddo	! while (estimator > eps)
   !----------------------------
   ! Check sum rules of the obtained coefficients
   Omega =  w_plasma(1d6*1.0d24)	! module "CDF_Ritchi"
   call sum_rules(CDF_fit%A,  CDF_fit%E0, CDF_fit%Gamma, E(peaks_pos(1)), Omega, ksum, fsum)	! module "CDF_Ritchi"
   write(*,'(a,f6.3,a,f6.3,a,f6.3)') 'Sum rules:', ksum, ' Ne=', Ne, ' F-sum rule:', fsum
   print*, 'iterations=', Niter
   print*, 'estimator=', estimator
   print*, '------------------------'
   
   do i = 1, Nfit
       write(*,'(a,i3,f,f,f)') 'COEF: ', i, CDF_fit%A(i), CDF_fit%E0(i), CDF_fit%Gamma(i)
   enddo
   
   ! Negative sum rule indicate that convergence was not achieved, redo the fitting:
   if ((ksum <= 0.0d0) .or. ( ABS(ksum-Ne)/Ne > 0.5d0 )) then	! deviation is more then 50%
      if (one_more_try) then
         if (allocated(CDF_fit%A)) deallocate(CDF_fit%A)
         if (allocated(CDF_fit%Gamma)) deallocate(CDF_fit%Gamma)
         if (allocated(CDF_fit%E0)) deallocate(CDF_fit%E0)
         if (allocated(CDF_temp%A)) deallocate(CDF_temp%A)
         if (allocated(CDF_temp%Gamma)) deallocate(CDF_temp%Gamma)
         if (allocated(CDF_temp%E0)) deallocate(CDF_temp%E0)
         if (allocated(der_CDF)) deallocate(der_CDF)
         if (allocated(CDF%A)) deallocate(CDF%A)
         if (allocated(CDF%Gamma)) deallocate(CDF%Gamma)
         if (allocated(CDF%E0)) deallocate(CDF%E0)
         if (allocated(lambda%A)) deallocate(lambda%A)
         if (allocated(lambda%E0)) deallocate(lambda%E0)
         if (allocated(lambda%Gamma)) deallocate(lambda%Gamma)
         if (allocated(f_cdf)) deallocate(f_cdf)
         if (allocated(peaks_pos)) deallocate(peaks_pos)
         one_more_try = .false.
         N_iter_max = N_iter_max*2	! try longer
         print*, 'LET ME TRY ONE MORE TIME...'
         goto 9999 ! do over
      endif
   endif
!    PAUSE 'CHECK'
   
   CDF = CDF_fit	! save for output
   !----------------------------
   ! Clean up at the end:
   if (allocated(f_cdf)) deallocate(f_cdf)
   if (allocated(peaks_pos)) deallocate(peaks_pos)
   if (allocated(CDF_fit%A)) deallocate(CDF_fit%A)
   if (allocated(CDF_fit%Gamma)) deallocate(CDF_fit%Gamma)
   if (allocated(CDF_fit%E0)) deallocate(CDF_fit%E0)
   if (allocated(CDF_temp%A)) deallocate(CDF_temp%A)
   if (allocated(CDF_temp%Gamma)) deallocate(CDF_temp%Gamma)
   if (allocated(CDF_temp%E0)) deallocate(CDF_temp%E0)
   if (allocated(der_CDF)) deallocate(der_CDF)
   if (allocated(lambda%A)) deallocate(lambda%A)
   if (allocated(lambda%E0)) deallocate(lambda%E0)
   if (allocated(lambda%Gamma)) deallocate(lambda%Gamma)
end subroutine fit_oscillators_CDF_OLD



pure subroutine find_number_of_peaks(y_data, Npeaks, peaks_pos)
   real(8), dimension(:), intent(in) :: y_data	! array of data points
   integer, intent(out) :: Npeaks	! number of peaks in the data
   integer, dimension(:), allocatable, intent(out) :: peaks_pos	! positions of the peaks
   integer :: i, N, last_pos
   logical :: lp2 ,lp1, rp1 ,rp2
   N = size (y_data)	! size of array
   Npeaks = 0
   last_pos = -5	! to start
   if (N > 5) then	! assume no more then 1 peak for points less then 5
      do i = 3, N-3	! search thru all the data points
         ! peak is defined as local maximum with at least 5 points:
         ! if two points are increasing:
!          lp2 = (y_data(i-2) <= y_data(i-1))	! second left point
         lp2 = .true.
         lp1 = (y_data(i-1) <= y_data(i))	! first left point
         ! and the next two points are decreasing:
         rp1 = (y_data(i) > y_data(i+1))	! first right point
         rp2 = (y_data(i+1) > y_data(i+2))	! second right point
         ! then it's a local maximum:
         if (lp2 .and. lp1 .and. rp1 .and. rp2) then ! it's probably a peak
            if (i > last_pos + 5) then	! consider it as a new peak
               Npeaks = Npeaks + 1
               call extend_array_size_by_one(peaks_pos)	! module "Little_subroutines"
               peaks_pos(Npeaks) = i	! save location of the peak
               last_pos = i
            else	! don't consider it as a new peak, let it be the same peak since they are too close
            endif
         endif
      enddo
   endif
   if (Npeaks <= 0) then
      Npeaks = 1	! at least one peak
      call extend_array_size_by_one(peaks_pos)	! module "Little_subroutines"
      peaks_pos(Npeaks) = transfer(maxloc(y_data),1)	! save location of the peak
   endif
end subroutine find_number_of_peaks



pure subroutine find_number_of_minima(y_data, y_start, Nminima, minima_pos)
   real(8), dimension(:), intent(in) :: y_data	! array of data points
   integer, intent(in) :: y_start	! where to start counting from
   integer, intent(out) :: Nminima	! number of peaks in the data
   integer, dimension(:), allocatable, intent(out) :: minima_pos	! positions of the peaks
   integer :: i, N, last_pos
   logical :: lp2 ,lp1, rp1 ,rp2
   N = size (y_data)	! size of array
   Nminima = 0
   last_pos = -5	! to start
   if (N > 5) then	! assume no more then 1 minimum for points less then 5
      do i = y_start+2, N-3	! search thru all the data points
         ! minimum is defined as local minimum with at least 5 points:
         ! if two points are increasing:
         lp2 = (y_data(i-2) >= y_data(i-1))	! second left point
                  lp2 = .true.
         lp1 = (y_data(i-1) >= y_data(i))	! first left point
         ! and the next two points are decreasing:
         rp1 = (y_data(i) < y_data(i+1))	! first right point
         !rp2 = (y_data(i+1) < y_data(i+2))	! second right point
         rp2 = .true.
         ! then it's a local minimum:
         if (lp2 .and. lp1 .and. rp1 .and. rp2) then ! it's a minimum
            if (i > last_pos + 5) then	! consider it as a new minimum
               Nminima = Nminima + 1
               call extend_array_size_by_one(minima_pos)	! module "Little_subroutines"
               minima_pos(Nminima) = i	! save location of the minimum
               last_pos = i	! save the last peak position
            else	! don't consider it as a new minimum, let it be the same minimum since they are too close
            endif
         endif
      enddo
   endif
   ! Always consider the last point of the array as a minimum:
   Nminima = Nminima + 1	! at least one minimum
   call extend_array_size_by_one(minima_pos)	! module "Little_subroutines"
   minima_pos(Nminima) = size(y_data)	! save location of the minimum
end subroutine find_number_of_minima



pure function mean_square_deviation(y_data, func) result (MSD)
   real(8), dimension(:), intent(in) :: y_data, func
   real(8) :: MSD
   real(8), dimension(size(y_data)) :: ri
   ri = residuals(y_data, func)	! see below
   MSD = SUM(ri*ri)	! mean square deviation as residuals squared
end function mean_square_deviation


pure function residuals(y_data, func) result (r)
   real(8), dimension(:), intent(in) :: y_data, func
   real(8), dimension(size(y_data)) :: r	! residuals
   r = y_data - func
end function residuals




END MODULE CDF_get_from_data
