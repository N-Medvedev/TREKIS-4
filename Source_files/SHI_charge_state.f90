! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2018-2020
! 1111111111111111111111111111111111111111111111111111111111111
! Module contains calculations of effective charge state of an ion
module SHI_charge_state
use Universal_constants
use Relativity, only : velosity_from_kinetic_energy

implicit none

 contains


pure function Equilibrium_charge_SHI(Ekin, Mass, ZSHI, Zmean, Kind_Zeff, fixed_Zeff) result (Zeff)  ! Equilibrium charge
   real(8) Zeff	! effective SHI state
   real(8), intent(in) :: Ekin	! [eV] kinetic energy of SHI
   real(8), intent(in) :: Mass  ! [kg] SHI mass
   real(8), intent(in) :: ZSHI	! SHI atomic number
   real(8), intent(in) :: Zmean	! mean atomic number of elements in the target
   integer, intent(in) :: Kind_Zeff	! model for effective charge
   real(8), intent(in) :: fixed_Zeff	! for the case of user-provided fixed charge of SHI
   !--------------------------
   real(8) x, x2, x4, Zt, Zp, vp, vpvo, c1, c2
   !vp = sqrt(2.0d0*Ekin*g_e/(Mass*g_amu))  ! SHI velocity
   vp = velosity_from_kinetic_energy(Ekin, Mass, afs=.false.)     ! module "Relativity"
   
   
   !Zp = dble(ZSHI) ! SHI atomic number
   Zp = (ZSHI) ! SHI atomic number
   select case (Kind_Zeff)   ! 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande;
      case (1)   ! Original Bohr:
         Zeff = Zp*(1.0d0-exp(-(vp/g_v0/(Zp**(0.66666666d0)))))
      case (2)   ! Nikolaev, Dmitriev, Phys. Lett. 28A, 277 (1968):
         c1 = 0.6d0      ! k
         c2 = 0.45d0     ! alpha
         Zeff = Zp*(1.0d0 + (vp/(Zp**c2*g_v0*4.0d0/3.0d0))**(-1.0d0/c1))**(-c1)
      case (3) ! Schiwietz et al, NIMB 225, 4 (2004):
         Zt = Zmean	! mean atomic number of target atoms
         c1 = 1.0d0 - 0.26d0*exp(-Zt/11.0d0 - (Zt-Zp)*(Zt-Zp)/9.0d0)
         vpvo = Zp**(-0.543d0)*vp/g_v0
         c2 = 1.0d0 + 0.03d0*vpvo*dlog(Zt)
         x = c1*(vpvo/c2/1.54d0)**(1.0d0 + 1.83d0/Zp)
         x2 = x*x
         x4 = x2*x2
         Zeff = Zp*(8.29d0*x + x4)/(0.06d0/x + 4.0d0 + 7.4d0*x + x4)
      case (4)
         Zeff = fixed_Zeff
      case default ! Barkas:
         Zeff = Zp*(1.0d0-exp(-(vp*125.0d0/g_cvel/(Zp**(0.66666666d0)))))
   endselect
end function Equilibrium_charge_SHI

 

end module SHI_charge_state
