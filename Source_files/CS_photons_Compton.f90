! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains cross sections for phoon scattering: Compton
! References used in the module:
! [1]  F. Salvat, J. M. Fernandez-Varea, E. Acosta, J. Sempau "PENELOPE – A Code System for Monte Carlo Simulation of Electron and Photon Transport", OECD (2001)


module CS_photons_Compton
use Universal_constants
use Objects
use Read_input_data, only: m_input_folder, m_photon_CS, m_photon_Compton, m_folder_materials
use Dealing_with_files, only: read_file
use CS_general_tools, only: MFP_from_sigma
use Little_subroutines, only: interpolate_data_single

implicit none

real(8) :: m_Compton_factor, m_Cr, m_two_third, m_Thom_factor, m_five_sixteenth

parameter (m_Compton_factor = g_h/(g_me*g_e*g_e/(4.0d0*g_Pi*g_e0)))
parameter (m_Cr = 1.0093d0)	! [1], Page 59, under Eq.(2.92)
parameter (m_two_third = 2.0d0/3.0d0)
parameter (m_Thom_factor = 20.6074d0)	! [1] Constant for Eq.(2.5)
parameter (m_five_sixteenth = 5.0d0/16.0d0)

 contains


!---------------------------------------------------------------------------------
! Compton subroutines:

subroutine get_photon_Compton(Material, numpar, Err)
   type(Target_atoms), dimension(:), intent(inout), target :: Material	!material parameters of each target that it's constructed of
   type(Num_par), intent(in), target :: numpar	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error log
   !--------------------------------
   real(8) :: sigma, sigma_cur, N_elem, elem_contrib
   integer :: i, j, k, N_targets, N_elements, Ngrid, N_shells, FN, m
   integer :: Reason, count_lines, Nsiz
   character(300) :: Folder_with_CS, Path, command, Model_name, File_name, Path_total
   logical :: file_exist, read_well
   real(8), pointer :: E
   character, pointer :: path_sep
   type(Atom_kind), pointer :: Element
   
   write(*, '(a)', advance='no') ' Obtaining photon Compton scattering cross sections...'
   
   path_sep => numpar%path_sep
   Path = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_photon_CS))	! Electron CSs are storred in this folder
   
   N_targets = size(Material)	! that's how many different targets user specified
   
   ! Get electron MFPs:
   TRGT:do i = 1, N_targets	! for each target
      N_elem = dble(SUM(Material(i)%Elements(:)%percentage))	! number of elements in this compound material
      N_elements = size(Material(i)%Elements)	! that's how many different elements are in this target
      LMNT:do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         ! Check if folder for this element already exists:
         Folder_with_CS = trim(adjustl(Path))//path_sep//trim(adjustl(Element%Name))
         inquire(DIRECTORY=trim(adjustl(Folder_with_CS)),exist=file_exist)    ! check if input file excists
         if (.not.file_exist) then	! to make sure that such a folder is present (even if empty)
            ! Create a new directory for output files:
            command='mkdir '//trim(adjustl(Folder_with_CS))	! to create a folder use this command
            CALL system(command)	! create the folder
         endif
         !-------------------------------------
          ! Set part of the file name corresponding to the model used for Compton CS:
         select case (numpar%Ph_Compton)
         case (1)	! PENELOPE
            Model_name = 'PENELOPE'
         case default	! exclude
            Model_name = 'NO'
         end select
         
         ! For Compton CSs use the same grid as for photon absorption:
         Ngrid = size(Element%Phot_absorption%E)
         allocate(Element%Phot_compton%E(Ngrid))
         Element%Phot_compton%E = Element%Phot_absorption%E
         allocate(Element%Phot_compton%Total(Ngrid))
         allocate(Element%Phot_compton%Per_shell(Element%N_shl,Ngrid))
         allocate(Element%Phot_compton%MFP(Element%N_shl,Ngrid))
         N_shells = Element%N_shl
         allocate(Element%Phot_compton%Total_MFP(Ngrid))
         
         ! Calculate total cross sections for all shells of this element:
         count_lines = 0	! to count lines
         do k = 1, N_shells
            ! Check file with CSs:
            File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_photon_Compton))//'_'//trim(adjustl(Element%Shell_name(k)))//'_'//trim(adjustl(Model_name))//'.dat'
            inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
!             if (file_exist) then	! just read from the file, no need to recalculate:
            if ((file_exist) .and. (.not.numpar%recalculate_MFPs)) then    ! just read from the file, no need to recalculate:
               open(newunit = FN, FILE = trim(adjustl(File_name)),action='read')
               ! Get the Compton total cross section:
               do m = 1, Ngrid	! for all energy grid points:
                  E => Element%Phot_compton%E(m)	! photon energy [eV]
                  read(FN,'(es,es,es)', IOSTAT=Reason) E, Element%Phot_compton%Per_shell(k,m), Element%Phot_compton%MFP(k,m)
                  call read_file(Reason, count_lines, read_well)	! module "Dealing_wil_files"
                  if (.not.read_well) then
                     close(FN)	! redo the file
                     goto 8390
                  endif
!                   print*, k, m, Element%Phot_compton%E(m), Element%Phot_compton%Per_shell(k,m)
               enddo
            else	! no such file => create it
8390        open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
               do m = 1, Ngrid	! for all energy grid points:
                  E => Element%Phot_compton%E(m)	! photon energy [eV]
                  call get_Compton_CS(E, Element, numpar, k, sigma)	! see below
!                   write(*,'(a,es,es)') trim(adjustl(Element%Name))//' '//trim(adjustl(Element%Shell_name(k))), E, sigma
                  Element%Phot_compton%Per_shell(k,m) = sigma	! [A^2]
                  Element%Phot_compton%MFP(k,m) = MFP_from_sigma(sigma,  1.0d24)    ! module "CS_general_tools"
                  write(FN,'(es,es,es)') E, sigma, Element%Phot_compton%MFP(k,m)
               enddo ! m = 1, Ngrid
            endif ! (file_exist)
            close(FN)
            ! Normalize MFPs to material density:
            ! Take into account the density of atoms of this particular kind:
            elem_contrib = dble(Element%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
            Element%Phot_compton%MFP(k,:) = Element%Phot_compton%MFP(k,:) * 1.0d24/(Material(i)%At_Dens * elem_contrib)
         enddo ! k = 1, N_shells
         ! And the total CS for this element:
         do m = 1, Ngrid	! for all energy grid points:
            Element%Phot_compton%Total(m) = SUM(Element%Phot_compton%Per_shell(:,m))	! [A^2]
            ! Get also total MFP:
            Element%Phot_compton%Total_MFP(m) = MFP_from_sigma(Element%Phot_compton%Total(m),  Material(i)%At_Dens)    ! module "CS_general_tools"
         enddo ! m = 1, Ngrid
         File_name = trim(adjustl(Folder_with_CS))//path_sep//trim(adjustl(m_photon_Compton))//'_total_'//trim(adjustl(Model_name))//'.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there
         if (.not.file_exist .or. numpar%recalculate_MFPs) then	! only create it if file does not exist
            open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
            do m = 1, Ngrid	! for all energy grid points:
               write(FN,'(es,es,es)') Element%Phot_compton%E(m), Element%Phot_compton%Total(m), Element%Phot_compton%Total_MFP(m)
            enddo ! m = 1, Ngrid
            close(FN)
         endif
      enddo LMNT
      
      
      ! And the total CS:
      Path_total = trim(adjustl(m_input_folder))//path_sep//trim(adjustl(m_folder_materials))
      Path_total = trim(adjustl(Path_total))//path_sep//trim(adjustl(Material(i)%Name))
      File_name = trim(adjustl(Path_total))//path_sep//trim(adjustl(m_photon_Compton))//'_total_'//trim(adjustl(Model_name))//'.dat'
      Nsiz = size(Material(i)%Ph_absorption_total%E)   ! grid for total cross section is different from core-shells grid (but the same as valence grid)
      allocate(Material(i)%Ph_Compton_total%E(Nsiz))
      Material(i)%Ph_Compton_total%E = Material(i)%Ph_absorption_total%E
      allocate(Material(i)%Ph_Compton_total%Total(Nsiz))
      allocate(Material(i)%Ph_Compton_total%Total_MFP(Nsiz))
      Material(i)%Ph_Compton_total%Total(:) = 0.0d0
      Material(i)%Ph_Compton_total%Total_MFP(:) = 0.0d0
      ! Sum up CSs from each element:
      do j =1, N_elements	! for each element
         Element => Material(i)%Elements(j)	! all information about this element
         elem_contrib = dble(Material(i)%Elements(j)%percentage)/N_elem   ! element contribution to this compound (e.g. in SiO2: it is 1/3 of Si, 2/3 of O)
         do m = 1, Nsiz	! for all energy grid points:
            E => Material(i)%Ph_Compton_total%E(m)    ! electron energy [eV]
            ! Core shells:
            sigma = 0.0d0 ! to start with
            ! Get the cross section for this element on this grid point (which may not coincide with the grid for this element due to element-specific grid points):
            do k = 1, Element%N_shl		! for all shells in this elements
               call interpolate_data_single(Element%Phot_compton%E,  Element%Phot_compton%Per_shell(k,:), E, sigma_cur) ! module "Little_subroutines"
               ! Summing over all shells:
               sigma = sigma + sigma_cur*elem_contrib
            enddo ! k
            Material(i)%Ph_Compton_total%Total(m) = Material(i)%Ph_Compton_total%Total(m) + sigma   ! [A^2]
         enddo ! m = 1, Nsiz
      enddo !  j =1, N_elements
      ! Get also MFP:
      do m = 1, Ngrid	! for all energy grid points:
         Material(i)%Ph_Compton_total%Total_MFP(m) = MFP_from_sigma(Material(i)%Ph_Compton_total%Total(m),  Material(i)%At_Dens)    ! module "CS_general_tools"
      enddo
      ! Save the total cross section into the file:
      inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if this file is there      
      if (.not.file_exist .or. numpar%recalculate_MFPs) then	! only create it if file does not exist
         open(newunit = FN, FILE = trim(adjustl(File_name)),action='write')
         do m = 1, Nsiz	! for all energy grid points:
            write(FN,'(es,es,es)') Material(i)%Ph_Compton_total%E(m), Material(i)%Ph_Compton_total%Total(m), Material(i)%Ph_Compton_total%Total_MFP(m)
         enddo ! m = 1, Nsiz
         close(FN)
      endif
      
   enddo TRGT
   
   write(*, '(a)') ' Done.'
!    print*, 'Photon Compton scattering cross sections are obtained.'
   nullify(path_sep, Element, E)
end subroutine get_photon_Compton



! Compton scattering cross section per shell of an element:
pure subroutine get_Compton_CS(Ee, Element, numpar, k, sigma)
   real(8), intent(in) :: Ee	! [eV] photon energy
   type(Atom_kind), intent(in) :: Element	! data for this element
   type(Num_par), intent(in) :: numpar	! all numerical parameters
   integer, intent(in) :: k	! numbers of shell
   real(8), intent(out) :: sigma	! [A^2] cross section
   select case (numpar%Ph_Compton)
   case (1)	! PENELOPE
      sigma = Compton_CS(Ee, Element, k)	! see below
   case default	! exclude
      sigma = 0.0d0
   end select
end subroutine get_Compton_CS



! Compton total cross section according to
! [1] see Chapter 2.3
pure function Compton_CS(Ee, Element, k, mu_max_in) result(sigma)	! see below
   real(8) :: sigma	! [A^2] cross section
   real(8), intent(in) :: Ee	! [eV] photon energy
   type(Atom_kind), intent(in), target :: Element	! data for this element
   integer, intent(in) :: k	! index of shell
   real(8), intent(in), optional :: mu_max_in   ! upper limit of integration
   !------------------------------------
   real(8) :: dSigma, dSigma0, dSigma_mid, dS, eps, sig_step	! d Sigma / s Omega
   real(8) :: mu, mu0, dmu, dmu_half, mu_mid, mu_min, mu_max, dmu_max, dmu_min
   integer :: i, Ngrid
   real(8) :: Ip, mu_Ip
   
   Ip = Element%Ip(k)
   if (Ee > Ip) then	! Photon with such energy can participate in Compton scattering
   
      eps = 1.0d-15	! precision limit
      dS = 0.01d0	! maximal allowed change in dSigma per step
      Ngrid = 100	! grid point for integration over mu

!       mu_min = -1.0d0   ! cos(-Pi)
      ! Account for the ionization potential:
      mu_Ip = Compton_mu_from_Ec(Ee, Ip)    ! below
      mu_min = max(-1.0d0, mu_Ip)
      
      if (present(mu_max_in)) then
         if (mu_max_in <= 1.0d0) then
            mu_max = mu_max_in
         else
            mu_max = 1.0d0  ! cos(0)
         endif
      else
         mu_max = 1.0d0 ! cos(0)
      endif
      mu = mu_min   ! starting from cos(-Pi)
      dmu_max = (1.0d0 - mu_min)/dble(Ngrid)	! maximal allowed integration step
      dmu_min = dmu_max/100.0d0
      dmu = dmu_max	! start with it, and later reduce if needed
      dmu_half = dmu*0.50d0
      ! Start integration of differential cross section to obtain the total one:
      sigma = 0.0d0
      dSigma0 =  Compton_DCS(Ee, Element, k, mu)	! below
      i = 0   
      do while (mu < mu_max)
         i = i + 1	! steps counter
         ! Simpson 3/8 method of integration:
         mu0 = mu
         mu = mu + dmu
         if (mu > mu_max) then	! if by chance we exceeded the limit
            mu = mu_max
            dmu = mu_max - mu0
            dmu_half = 0.5d0*dmu
         endif
         dSigma = Compton_DCS(Ee, Element, k, mu)	! below
         ! Adaptive step: if it's too large, reduce it:
         if ((dSigma > eps) .and. (dSigma0 > eps)) then	! makes sense to take care of integration
            ! Check if it's too large change per step:
!             print*, 'Sigma', dSigma, dSigma0
            do while ((ABS(dSigma - dSigma0) > dS*dSigma0) .and. (dmu > dmu_min))
               dmu = dmu*0.5d0	! reduce timestep
               dmu_half = 0.5d0*dmu	! adjust half-step correspondingly
               mu = mu0 + dmu
               if (mu > mu_max) exit
               dSigma = Compton_DCS(Ee, Element, k, mu)	! below
               if (dmu < dmu_min) exit	! too smal step
            enddo
            ! Check if it's too small change per step (inefficient):
            do while ((ABS(dSigma - dSigma0) < 0.5d0*dS*dSigma0) .and. (dmu < dmu_max))
               dmu = dmu + dmu_min	! increase timestep
               dmu_half = 0.5d0*dmu	! adjust half-step correspondingly
               mu = mu0 + dmu
               if (mu > mu_max) exit
               dSigma = Compton_DCS(Ee, Element, k, mu)	! below
               if (dmu > dmu_max) exit	! too large step
            enddo
         else	! reset the step
            dmu = dmu_max
            dmu_half = dmu*0.5d0
            mu = mu0 + dmu
         endif
         ! Now proceed with the integration:
         if (mu > mu_max) then	! if by chance we exceeded the limit
            mu = mu_max
            dmu = mu_max - mu0
            dmu_half = dmu*0.50d0
         endif
         mu_mid =  mu0 + dmu_half
         dSigma_mid = Compton_DCS(Ee, Element, k, mu_mid)	! below
      
         ! Add up the contributions according to Simpson-3/8 scheme:
         sig_step = dmu/6.0d0*(dSigma0 + 4.0d0*dSigma_mid + dSigma)
         sigma = sigma + sig_step
!          write(*,'(a,f,f,es,es)') 'S', mu, dmu, dSigma, sigma
         ! Save the data for the next point of integration:
         dSigma0 = dSigma
      enddo
   else	! such photon energy is insufficient to Compton-scatter
      sigma = 0.0d0
   endif
   sigma = 2.0d0*g_Pi*sigma	! integrated over phi = (0,2Pi)
end function Compton_CS   


! [1], Eq.(2.52) or Eq.(2.46) in 2015-version of [1]:
pure function Compton_DCS(E, Element, k, mu) result(dSigma)	! see below
   real(8) :: dSigma	! differential cross section: [A^2/rad]
   real(8), intent(in) :: E	! [eV] photon energy
   type(Atom_kind), intent(in), target :: Element	! data for this element
   integer, intent(in) :: k	! numbers of shell
   real(8), intent(in) :: mu	! cos(teta)
   !-----------------------------
   real(8) :: Ommu, Ec, Pmax, EIp, sqrP, temp
   real(8) :: Ni, sgnP, expAi, expins, Jo, EcE
   real(8) :: Ip
   Ip = Element%Ip(k)
   if (E > Ip) then	! it can participate in Compton scattering
      ! Get Compton energy:
      Ommu = 1.0d0 - mu
      Ec = g_me_eV*E / (g_me_eV + E*Ommu)	! [1] Eq.(2.32)
      ! Get Compton maximum transferred momentum:
      EIp = E - Ip
      temp = E*EIp*Ommu
      sqrP = 2.0d0*temp + Ip*Ip
      Pmax = g_e*(temp - g_me_eV*Ip) / (g_cvel*sqrt(sqrP))	! [1] Eq.(2.49) // or (2.43) in (2015) version of [1]
      ! Get an integral of the Compton profile:
      Jo = Element%Compton(k) * m_Compton_factor	! Compton profile for P=0 read from the database
      sgnP = sign(1.0d0,Pmax)	! sign(Pmax)
      expins = 1.0d0 + 2.0d0*sgnP*Jo*Pmax
      expAi = 0.5d0*( 1.0d0 - expins*expins )
      Ni = 0.50d0*(1.0d0 + sgnP*(1.0d0 - exp(expAi) ) )	! [1] Eq.(2.58)
      ! Get differential Compton cross section:
      EcE = Ec/E
      dSigma = 0.5d0*g_r0*g_r0*EcE*EcE*( EcE + 1.0d0/EcE - 1.0d0 + mu*mu) * Ni * Element%Ne_shell(k)	! [1] Eq.(2.52), but per shell (without summation over shells) // (2.46)
   else	! it can not participate in Compton scattering
      dSigma = 0.0d0
   endif
end function Compton_DCS



pure function Compton_Ec(E, mu) result(Ec)
   real(8) Ec	! Compton energy
   real(8), intent(in) :: E	! [eV] photon energy
   real(8), intent(in) :: mu	! cos(teta)
   Ec = g_me_eV*E / (g_me_eV + E*(1.0d0 - mu))	! [1] Eq.(2.32) // (2.26) in 2015
end function Compton_Ec


pure function Compton_mu_from_Ec(E, Ec) result(mu)
   real(8) mu   ! cos(theta)
   real(8), intent(in) :: E ! [eV] photon energy
   real(8), intent(in) :: Ec ! [eV] energy transfer
   if (Ec > 1.0-6) then
      mu = 1.0d0 - g_me_eV/E * (E/Ec - 1.0d0)  ! inverse of Eq.(2.26) [1]
   else
      mu = -1.0d6
   endif
end function Compton_mu_from_Ec


pure function Compton_electron_emission_angle(Eph, EphC, mu) result (mue)
   real(8) mue  ! cos(theta_e)  ! Eq.(2.72) in 2015-edition of [1]
   real(8), intent(in) :: Eph   ! [eV] incoming photon energy
   real(8), intent(in) :: EphC  ! [eV] Compton photon energy
   real(8), intent(in) :: mu    ! cos(theta) scattering angle of the photon
   mue = (Eph - EphC * mu)/sqrt(Eph*Eph + EphC*EphC - 2.0d0*Eph*EphC*mu)
end function Compton_electron_emission_angle


end module CS_photons_Compton
