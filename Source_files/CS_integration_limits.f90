! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2019-2020
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains subroutines for integration limits for cross sections
! References used:
! [1] N.Medvedev, A.E.Volkov, J. PHys. D (2020)    DOI: https://doi.org/10.1088/1361-6463/ab7c09
! [2] F.Salvat, J.M. Fernandez-Varea, J.Sempau, "PENELOPE-2014: A Code System for Monte Carlo Simulation of Electron and Photon Transport" (2014)

MODULE CS_integration_limits
use Relativity, only : rest_energy

implicit none
 contains


!-------------------------------------------------------------
! Choice of the integration limit depending on the model used:

pure subroutine define_integration_limits(E_start, E_end, CS_method, Emin, Emax, E0_min, E0_max)
   real(8), intent(out) :: E_start, E_end
   integer, intent(in) :: CS_method
   real(8), intent(in) :: Emin, Emax
   real(8), intent(in) :: E0_min, E0_max
   !------------------
   select case (abs(CS_method))
   case default   ! Old default, don't redefine
      E_start = Emin
      E_end = Emax

   case (3) ! New default, exclude empty space
      E_start = max(Emin, E0_min)
      E_end = min(Emax, E0_max)
   end select

end subroutine define_integration_limits


!-------------------------------------------------------------
! Relativistic transferred momentum limits:
pure function Q_min(M, mt, E, W) result(Qmin)   ! Eq.(A.31) [2]
   real(8) Qmin ! [eV] energy corresponding to the minimum transfered momentum
   real(8), intent(in) :: M, mt ! [kg] mass of the incident and target particle
   real(8), intent(in) :: E, W  ! [eV] incident energy and transfered energy
   real(8) :: Mc2, mtc2, Emc2, term1, term2, term3
   Mc2 = rest_energy(M)  ! module "Relativity"
   Emc2 = E + 2.0d0*Mc2
   mtc2 = rest_energy(mt)  ! module "Relativity"
   term1 = sqrt(E*Emc2)
   term2 = sqrt( (E-W)*(Emc2-W) )
   term3 = term1 - term2
   Qmin = sqrt( term3*term3 + mtc2*mtc2 ) - mtc2
end function Q_min


pure function Q_max(M, mt, E, W) result(Qmax)    ! Eq.(A.31) [2]
   real(8) Qmax ! [eV] energy corresponding to the minimum transfered momentum
   real(8), intent(in) :: M, mt ! [kg] mass of the incident and target particle
   real(8), intent(in) :: E, W  ! [eV] incident energy and transfered energy
   real(8) :: Mc2, mtc2, Emc2, term1, term2, term3
   Mc2 = rest_energy(M)  ! module "Relativity"
   Emc2 = E + 2.0d0*Mc2
   mtc2 = rest_energy(mt)  ! module "Relativity"
   term1 = sqrt(E*Emc2)
   term2 = sqrt( (E-W)*(Emc2-W) )
   term3 = term1 + term2
   Qmax = sqrt( term3*term3 + mtc2*mtc2 ) - mtc2
end function Q_max



!-------------------------------------------------------------
! Non-relativistic transferred momentum limits:
pure function Q_min_nonrel(E, W, M, mt) result(Qmin)
   real(8) Qmin ! [eV] energy corresponding to the minimum transfered momentum
   real(8), intent(in) :: E, W  ! [eV] incident energy and transfered energy
   real(8), intent(in), optional :: M, mt   ! masses of the incident particle and target particle
   Qmin = sqrt(E) - sqrt(E-W)   ! [sqrt(eV)]
   Qmin = Qmin*Qmin ! [eV]
   if (present(M) .and. present(mt)) Qmin = Qmin*M/mt
end function Q_min_nonrel


pure function Q_max_nonrel(E, W, M, mt, hw_phonon) result(Qmax)
   real(8) Qmax ! [eV] energy corresponding to the minimum transfered momentum
   real(8), intent(in) :: E, W  ! [eV] incident energy and transfered energy
   real(8), intent(in), optional :: M, mt   ! masses of the incident particle and target particle
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   Qmax = sqrt(E) + sqrt(E-W)    ! [sqrt(eV)]
   Qmax = Qmax*Qmax ! [eV]
   if (present(M) .and. present(mt)) Qmax = Qmax*M/mt
   if (present(hw_phonon)) Qmax = max(Qmax,hw_phonon)
end function Q_max_nonrel


!-------------------------------------------------------------
! Limits of transferred energy:
pure function W_min(Ip, Mc2, mtc2, E, E0) result(Wmin)
   real(8) :: Wmin  ! minimal transfered energy [eV]
   real(8), intent(in) :: Ip    ! [eV] ionization potential, or band gap, or the limiting low energy (which may also be zero)
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident and target particle
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in), optional :: E0 ! [eV] position of maximum of delta-model CDF
   !-----------------------
   real(8) :: E0min, mtc2_2, arg
   
!    Wmin = Ip
   mtc2_2 = mtc2*mtc2
   if (E0/E*Mc2/mtc2 > 1.0d-6) then ! exact nonrelativistic expression, Eq.(24) [1]
      arg = 1.0d0 - E0/E*(1.0d0+Mc2/mtc2)
      if (arg > 0.0d0) then ! sqrt is defined:
         E0min = mtc2/(Mc2+mtc2)*(E0 + 2.0d0*Mc2/(Mc2+mtc2)*E*(1.0d0 - sqrt(arg) ) )
      else  ! sqrt is undefined, ionization is impossible
         E0min = E
      endif
   else   ! approximation would do:
      E0min = E0*(1.0d0+0.25d0*E0/E*Mc2/mtc2)    ! Last eq. in Appendix A [1]
   endif
   Wmin = max(Ip,E0min)
end function W_min


pure function W_max(Mc2, mtc2, identical, E, E0, hw_phonon) result(Wmax)    ! Eq.(15) [1]
   real(8) :: Wmax  ! maximal transfered energy [eV]
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident particles, and the scattering center
   logical, intent(in) :: identical   ! identical particles, true or not
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in) :: E0    ! [eV] oscilator energy (ionization potential)
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8) :: nu, A, B, C, D, F, eps, mtc22, E02, E2, DDFC
   eps = 1.0d-6 ! within what margin we consider particles masses equal
   
   mtc22 = 2.0d0*mtc2
   E2 = 2.0d0*E
   E02 = 2.0d0*E0

   nu = E + 2.0d0*Mc2
   A = mtc22 + E + nu - E02
   B = E0*E0 - mtc22*E0 - E2*nu
   C = A*A - 2.0d0*E2*nu
   D = A*B + E2*nu*(E + nu)
   F = B*B - 4.0d0*E*E*nu*nu
   
   DDFC = D*D - F*C
   if (DDFC > 0.0d0) then   ! and only then the SQRT is >0
      Wmax = (-D + sqrt(abs(DDFC)) )/C   

      if (identical) then   ! elecron-electron, etc.
         Wmax = (Wmax + E0)*0.50d0
      endif
   else ! undefined transfer energy, event is impossible
      Wmax = E
   endif
   
   if (present(hw_phonon)) then
      Wmax = max(Wmax, hw_phonon)
   endif
end function W_max



pure function W_max_nonrel(Mc2, mtc2, identical, E, E0, hw_phonon) result(Wmax) ! Eq.(24) [1]
   real(8) :: Wmax  ! maximal nonrelativistic transfered energy [eV]
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident particles, and the scattering center
   logical, intent(in) :: identical   ! identical particles, true or not
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in) :: E0    ! [eV] oscilator energy (ionization potential)
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8) :: eps, Emtc22, EMmt, Mm, Wmax1
!    eps = 1.0d-6 ! within what margin we consider particles masses equal

   Wmax = mtc2/(Mc2+mtc2)*(E0 + 2.0d0*Mc2/(Mc2+mtc2)*E*(1.0d0 + sqrt(abs(1.0d0 - E0/E*(1.0d0+Mc2/mtc2))) ) )
   
   if (identical) then   ! elecron-electron, etc.
      Wmax =  (Wmax + E0)*0.50d0
   endif
   
   if (present(hw_phonon)) Wmax = max(Wmax, hw_phonon)
end function W_max_nonrel


pure function W_max_nonrel_free(Mc2, mtc2, identical, E, E0, hw_phonon) result(Wmax) ! Eq.(24) [1]
   real(8) :: Wmax  ! maximal nonrelativistic transfered energy [eV]
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident particles, and the scattering center
   logical, intent(in) :: identical   ! identical particles, true or not
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in) :: E0    ! [eV] oscilator energy (ionization potential)
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8) :: eps, Emtc22, EMmt, Mm, Wmax1

   Mm = Mc2+mtc2
   Wmax = 4.0d0*mtc2*Mc2/(Mm*Mm) * E    ! free particle limit
   
   if (identical) then   ! electron-electron, etc.
      Wmax =  (Wmax + E0)*0.50d0
   endif
   
   if (present(hw_phonon)) Wmax = max(Wmax, hw_phonon)
end function W_max_nonrel_free



pure function W_max_OLD(Mc2, mtc2, identical, E, Ip, hw_phonon) result(Wmax)
   real(8) :: Wmax  ! maximal transfered energy [eV]
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident particles, and the scattering center
   logical, intent(in) :: identical   ! identical particles, true or not
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in) :: Ip    ! [eV] minimum energy (ionization potential)
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8) :: eps, Emtc22, EMmt, Mm, Wmax1
!    eps = 1.0d-6 ! within what margin we consider particles masses equal
   if (identical) then   ! elecron-electron, etc.
      Wmax = (E+Ip)*0.50d0
   else ! not identical particles
      Emtc22 =  2.0d0*mtc2*E
      Mm = Mc2+mtc2
      Wmax = Emtc22*(E + Mc2*2.0d0) / (Emtc22 + Mm*Mm)    ! Eq.(A.25) [2]
      if (present(hw_phonon)) then  ! it's scattering on atoms/phonons, so may be it's a phonon:
         if (Wmax < hw_phonon) Wmax = hw_phonon
      endif
      ! Classical limit:
!       Wmax1 = 4.0d0*Mc2*mtc2/(Mm*Mm)*E
!       write(*,'(a,e,e,e,e,e)') 'Wmax', E, Wmax1, Wmax, Mc2, mtc2
   endif
end function W_max_OLD


pure function W_min_OLD(Ip, Mc2, mtc2, E, E0) result(Wmin)
   real(8) :: Wmin  ! minimal transfered energy [eV]
   real(8), intent(in) :: Ip    ! [eV] ionization potential, or band gap, or the limiting low energy (which may also be zero)
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident and target particle
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in), optional :: E0 ! [eV] position of maximum of delta-model CDF
   !-----------------------
   real(8) :: E0min
   real(8) :: mtc2_2
   
   Wmin = Ip

   if (present(E0)) then
      if (abs(Mc2 - mtc2)/Mc2 < 1.0d-6) then   ! e.g. electron-electron
         E0min = E0*(1.0d0-0.25d0*E0/E)    ! delta-model CDF
      else ! e.g. ion-electron
           if (E < 1.0d9*E0) then ! exact expression:
              mtc2_2 = mtc2*mtc2
              E0min = (-2.0d0*Mc2*sqrt(E)*(Mc2*sqrt(E)-sqrt(Mc2*E0*mtc2-E0*mtc2_2+mtc2_2*E))/(Mc2-mtc2)+2.0d0*Mc2*E-E0*mtc2)/(Mc2-mtc2)
           else   ! approximation would do:
              E0min = E0*(1.0d0-0.25d0*E0/E*Mc2/mtc2)    ! delta-model CDF
           endif
      endif
      Wmin = max(Wmin,E0min)
   endif
end function W_min_OLD





pure subroutine find_Wmax_equal_Wmin(Mc2, mtc2, identical, E, Ip, E0, Eeq, hw_phonon)   ! Eq.(19) [1]
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident particles, and the scattering center
   logical, intent(in) :: identical   ! identical particles, true or not
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in) :: Ip    ! [eV] minimum energy (ionization potential)
   real(8), intent(in) :: E0    ! [eV]
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8), intent(out) :: Eeq  ! point where Wmax=Wmin
   !----------------------------
   real(8) :: mtc2_2, Mc2_2, Eeq1, Mm, Mm4, adj_coef
  
   if (.not.identical) then ! scattering between two different particles (ion-electron, electron-ion, positron-electron, etc.)
      Eeq =E0*(mtc2+Mc2)/mtc2
      if (present(hw_phonon)) then  ! scattering on phonon:
         Eeq = max (Eeq, 0.25d0*E0*(-2.0d0*hw_phonon+E0)/(-hw_phonon+E0))
      endif
   else ! scattering between identical particles (electron-electron etc.)
      Eeq = E0*(1.0d0+0.75d0*sqrt(2.0d0))
   endif
end subroutine find_Wmax_equal_Wmin



 subroutine find_Wmax_equal_Wmin_OLD(Mc2, mtc2, identical, E, Ip, E0, Eeq, hw_phonon)
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident particles, and the scattering center
   logical, intent(in) :: identical   ! identical particles, true or not
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in) :: Ip    ! [eV] minimum energy (ionization potential)
   real(8), intent(in) :: E0    ! [eV]
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8), intent(out) :: Eeq  ! point where Wmax=Wmin
   !----------------------------
   real(8) :: mtc2_2, Mc2_2, Eeq1, Mm, Mm4, adj_coef
  
   if (.not.identical) then ! scattering between two different particles (ion-electron, electron-ion, positron-electron, etc.)
      if (present(hw_phonon)) then  ! scattering on phonon:
         Eeq =0.25d0*E0*(-2.0d0*hw_phonon+E0)/(-hw_phonon+E0)
      else ! scattering on a particle
         mtc2_2 = mtc2*mtc2
         Mc2_2 = Mc2*Mc2         
         Mm = Mc2+mtc2
         Mm4 = Mm*Mm*Mm*Mm
         adj_coef = 1.28d0  ! shift the crossing point by this much
         
         !Eeq = adj_coef*E0*(Mc2/mtc2)*(Mm4/abs(4.0d0*Mc2*Mc2*Mc2*Mc2 - Mm4))
         Eeq = E0*(Mc2/mtc2)*(Mm4/abs(4.0d0*Mc2*Mc2 - Mm4))
         
!            Eeq1 = 0.125d0*E0*Mm*Mm/(Mc2*mtc2)   ! just for comparison
!            write(*,'(a,e,e,e,f)') 'Eeq', E, Eeq, Eeq1, E0
      endif
   else ! scattering between identical particles (electron-electron etc.)
      Eeq = 5.0d0/4.0d0*E0 - Ip/2.0d0 + 0.25d0*sqrt(17.0d0*E0*E0 - 12.0d0*Ip*E0 + 4.0d0*Ip*Ip)
   endif
end subroutine find_Wmax_equal_Wmin_OLD
 
 
 
! For a given minimal transferred energy, what is ion energy that is sufficient to transfer such an amount:
pure function minimal_sufficient_E(Ip, Mc2, mtc2) result(Emin)
   real(8) :: Emin
   real(8), intent(in) :: Ip    ! [eV] transferred energy
   real(8), intent(in) :: Mc2 ! [eV] ion mass
   real(8), intent(in) :: mtc2  ! [eV] target particle mass (electron)
   real(8) :: Mc2_2, mtc2_2
   
   if (abs(Mc2-mtc2)/Mc2 < 1.0d-6) then ! identical
      Emin = Ip
   else ! not identical
      Mc2_2 = Mc2*Mc2
      mtc2_2 = mtc2*mtc2

      Emin = 0.50d0*(Ip*Mc2 - 2.0d0*mtc2*Mc2 + 2.0d0*Ip*mtc2 + sqrt(-Ip*Ip*Mc2_2 + 4.0d0*mtc2_2*Mc2_2 - &
               6.0d0*mtc2_2*Mc2*Ip + 2.0d0*Ip*Ip*mtc2_2 + 2.0d0*Ip*Mc2*Mc2_2))/(Mc2 - Ip)
   endif
end function minimal_sufficient_E


 
END MODULE CS_integration_limits
