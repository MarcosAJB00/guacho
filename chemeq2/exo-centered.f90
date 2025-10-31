!=======================================================================
!> @file exoplanet.f90
!> @brief Exoplanet problem module
!> @author M. Schneiter, C. Villarreal  D'Angelo, A. Esquivel
!> @date 4/May/2016

! Copyright (c) 2016 Guacho Co-Op
!
! This file is part of Guacho-3D.
!
! Guacho-3D is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief Exoplanet module
!> @details Problem Module for exoplanet

module exoplanet

  use parameters
  use globals, only: rank
  implicit none

  TYPE body
    real :: x,y,z
    real :: mass,radius,amdot
    real :: rho
    real :: temp
    real :: bfield
    real :: WR,WT,WV,WD
    real :: S0H, S0He1, S0He3
  END TYPE

  type(body) star,planet,barycenter

  real :: torb    !< planet: orbital period
  real :: rorb    !<  orbital radius
  real :: omegap  !< planet: angular velocity
!  real :: rho_p   !< density at rp
  real :: cs2     !< speed of sound at rp
  real :: Kconst  !< constant tripathi eq.19
  real :: rho0, R0         !< Inner region
  real :: rho_out, Rout   !< Outer region
  real :: Rbound           !< Boundary conditions
  real :: dens_c1, dens_c2 !< constants for density computation
  real :: energy_c1        !< constant for energy computation
  real :: alpha, ymol, fHII, fHeII, fHII_s, fHeII_s

  character (len=128) :: fileout
  integer             :: unitout
  real                :: RMdot  ! Radius at which the mass loss
                                ! rate is computed

contains

!=======================================================================

!> @brief Module initialization
!> @details Here the parameters of the Star are initialized, and scaled
!! to code units

subroutine init_exo()

  use constants, only : Rg, Msun, Yr, Ggrav, pi, mjup, au, day,Rjup, rsun,Kb, amh
  use parameters, only: mu, xphys

  implicit none

  !----------------STAR PARAMETERS ------------------
  star%mass      = 0.86*msun!1.119*msun
  star%radius    = 0.76*rsun!1.155*rsun
!  star%amdot     = 0.0 !1.3E-16*msun/yr    ! Stellar Mass Loss rate (g s^-1)
  star%WT        = 1.4e6 !1.35e6            ! Stellar temperature (K)
  !star%WR        = 8.16*star%radius        ! wind launching radius
  star%WV        = 550.0e5                  ! Stellar wind velocity (cm/s)
  star%WD        = 5.5e-21                  ! Stellar wind density  (g cm^-3)
  star%bfield    = 0.0                      ! Stellar magnetic field (g)
  star%S0H       = 1.3e14                   ! soft_euv lyman weinn flux
  star%S0He1     = 2.6e13                   ! hard_euv (1/cm2/s)
  star%S0He3     = 4.1e15                   ! soft_euv **(mid-UV)**
  !----------------PLANET PARAMETERS------------------
  planet%mass    = 0.0787*mjup
  planet%amdot   = 0.0                       ! Planetary Mass Loss rate (g/s)
  planet%rho     = 3.6e-15
  planet%radius  = 0.4466*rjup!1.359*Rjup
  planet%temp    = 1.0e3                    ! Planets temperature
!  planet%WT      = 0.0                       ! Planets temperature
!  planet%WR      = 0.0                       ! Planetary wind radius (cm)
!  planet%WV      = 0.0                       ! Planets wind velocity (cm/s)
!  planet%WD      = 0.0                       ! Planetary wind density
  planet%bfield  = 0.0                       ! Planetary magnetic field (g)

  alpha    = 1.9e-11!         1./3.       ! pop fraction between singlet and triplet state in he
  ymol     = 0.1!0.1           ! % of He from total density. ntH/ntHe = xmol/ymol. with xmol beinng % of H.
  fHII     = 1.7e-3       ! HII fraction in planetary atmosphere
  fHeII    = 1.9e-9       ! HeII fraction in planetary atmosphere
  fHII_s   = 1.0          ! HII  fraction for solar wind
  fHeII_s  = 1.0          ! HeII fraction for solar wind

  !ORBITAL PARAMETERS
  rorb = 0.0532*AU
  torb = 4.887802443*day

  cs2 = gamma*Rg/1.*planet%temp
  Kconst = planet%rho**(1.0 - gamma)*cs2
  Rbound = 0.75*planet%radius

  dens_c2 = (gamma - 1.0)/gamma*planet%mass*Ggrav/Kconst
  if (rank == master) print*, 'dens_c2=', dens_c2
  dens_c1 = planet%rho**(gamma - 1.0) - dens_c2/planet%radius

  R0 = 0.5*planet%radius
  rho0 = (dens_c1 + dens_c2/R0)**(1.0/(gamma - 1.0))
  if (rank == master) print*,'rho0 [g]=', rho0

  dens_c1 = rho0**(gamma - 1.0) - dens_c2/R0

  rho_out = planet%rho/10.0
  Rout = dens_c2/(rho_out**(gamma - 1.0) - dens_c1)
  if (rank == master) print*,'Rout [Rp]=',Rout/planet%radius

  energy_c1 = Kconst/(gamma-1.0)

! scaled parameters
  rorb = rorb/rsc
  torb = torb/tsc
  omegap = 2.0*pi/torb

  star%x = 0.0
  star%y = 0.0
  star%z = -rorb

  planet%x = 0.00
  planet%y = 0.00
  planet%z = 0.00

  barycenter%x = 0.0
  barycenter%y = 0.0
  barycenter%z = -rorb*star%mass/(star%mass + planet%mass)

  !wind launching radius viewed from the star
  star%WR = rorb-xphys/rsc/2
  if (rank == master) print*,'Rsw [Rs]=', star%WR*rsc/star%radius
  star%WV = star%WV/vsc

  unitout = 15
  RMdot   = 2.0 * Rout
  if (rank == master) then
    write(fileout,'(a)')  trim(outputpath)//'BIN/mass_loss_rate.dat'
    open(unit=unitout,file=fileout,status='unknown',access='append')
  end if

end subroutine init_exo

!=======================================================================

!> @brief Inject sources of wind
!> @details Imposes the sources of wond from the star and planet
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserver variables
!> @param real [time] time : current integration timr
  !--------------------------------------------------------------------
subroutine impose_exo(u,time)

  use constants, only : pi, amh, Kb
  use globals, only : coords, dx, dy, dz

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: x,y,z
  real :: xp,yp,zp,rp
  real :: xs,ys,zs,rs
  real :: velx,vely,velz,dens
  !real :: bx,by,bz,phi
  real :: neutral_h,ion_h,neutral_heS,neutral_heM,ion_he, electrons, ntot
  real :: passive, nht, nhet

#ifdef BFIELD
  real :: cpi
#endif
  integer :: i,j,k
  integer :: global_i,global_j,global_k
  logical :: set_cell
  real    :: energy
  real    :: vp, tp, fion, rr

  do i=nxmin,nxmax
    do j=nymin,nymax
      do k=nzmin,nzmax

        set_cell = .false.

        global_i = i + coords(0)*nx
        global_j = j + coords(1)*ny
        global_k = k + coords(2)*nz

        ! Position measured from the centre of the grid (planet)
        x = (float(global_i - nxtot/2) - 0.5)*dx
        y = (float(global_j - nytot/2) - 0.5)*dy
        z = (float(global_k - nztot/2) - 0.5)*dz

        xp = x - Planet%x
        yp = y - Planet%y
        zp = z - Planet%z

        ! Distance from the centre of the planet
        rp = sqrt(xp**2 + yp**2 + zp**2)*rsc     ! cgs

        nHt   = (1. - ymol)/ (1. + 3.*ymol)  !total hydrogen density
        nHet  =        ymol/ (1. + 3.*ymol)    !total he density

        ! IF INSIDE THE PLANET
        if(rp <= R0) then ! IF INSIDE RMASK

          set_cell = .true.

          dens = rho0     !total rho
          velx = 0.0
          vely = 0.0
          velz = 0.0

          neutral_h   = (1. - fHII) * nHt * dens/ rhosc        !neutral density of H [g/cm3]
          ion_h       = fHII * nHt * dens / rhosc

          neutral_heS = (1. - fHeII)* nHet * dens/ rhosc/ (1+alpha)       !neutral density of He singlet [g/cm3]
          neutral_heM = alpha*neutral_heS         ! neutral density of He triplet [g/cm3]
          ion_he      = fHeII * nHet * dens/ rhosc
          electrons   = ion_h + ion_he

          energy      = energy_c1*(dens)**gamma / vsc2 / rhosc
          passive     = 1000*dens/ rhosc
          dens        = dens / rhosc

        elseif(rp <= Rout)then

          set_cell = .true.

          velx = 0.0
          vely = 0.0
          velz = 0.0

          dens        = (dens_c1 + dens_c2/rp)**(1.0/(gamma - 1.0))

          neutral_h   = (1. - fHII) * nHt * dens/ rhosc    !neutral density of H [g/cm3]
          ion_h       = fHII * nHt * dens / rhosc

          neutral_heS = (1. - fHeII) * nHet * dens/ rhosc /(1+alpha) !neutral density of He singlet [g/cm3]
          neutral_heM = alpha*neutral_heS     ! neutral density of He triplet [g/cm3]
          ion_he      = fHeII * nHet * dens/ rhosc
          electrons   = ion_h + ion_he

          energy      = energy_c1*(dens)**gamma /vsc2 / rhosc
          passive     = 1000*dens / rhosc
          dens        = dens / rhosc

        else

          !set_cell = .true.

          !velx =0.0 ! vp*xp/rp
          !vely =0.0 ! vp*yp/rp
          !velz =0.0 ! vp*zp/rp

          !dens        = rho_out/10000.

          !neutral_h   = (1. - fHII_s) * nHt * dens/rhosc
          !ion_h       = fHII_s * nHt * dens/rhosc

          !neutral_heS = (1. - fHeII_s) * nHet * dens /(1+alpha)/rhosc
          !neutral_heM = alpha * neutral_heS
          !ion_he      = fHeII_s * nHet *dens /rhosc
          !electrons   = ion_h + ion_he

          !energy      = energy_c1*rho_out**gamma /vsc2/rhosc
          !passive     = -1000*dens/rhosc
          !dens        = dens/rhosc

          !          energy      = 0.5*dens*vp**2 + cv*(2.*dens-neutral_h)*tp
          
          !Stellar Wind
          set_cell = .true.
         
          xp = x - barycenter%x
          yp = y - barycenter%y
          zp = z - barycenter%z
          !
          xs = x - star%x
          ys = y - star%y
          zs = z - star%z
          rs = sqrt(xs**2+ys**2+zs**2)
          
          velx = star%WV*xs/rs - omegap*zp
          vely = star%WV*ys/rs
          velz = star%WV*zs/rs + omegap*xp   !code units
          
          dens = star%WD/rhosc!*star%WR**2/rs**2
          
          neutral_h = (1. - fHII_s) * nHt * dens
          ion_h     = fHII_s * nHt * dens
          
          neutral_heS = (1. - fHeII_s) * nHet * dens /(1+alpha)
          neutral_heM = alpha * neutral_heS
          ion_he      = fHeII_s * nHet *dens
          electrons   = ion_h + ion_he
          
          ntot = neutral_h + ion_h + neutral_heS + neutral_heM +ion_he
          
          energy = 0.5*(dens)*(velx**2 + vely**2 + velz**2)
          energy = energy + cv*(ntot + electrons)*star%WT*Kb/rhosc/vsc2
          
          passive = -1000*dens

        end if

        if(set_cell)then

          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz
          u(5,i,j,k) = energy

#ifdef BFIELD
          u(6,i,j,k) = bx
          u(7,i,j,k) = by
          u(8,i,j,k) = bz
          u(5,i,j,k) = u(5,i,j,k) + 0.5*(bx**2 + by**2 + bz**2)  !Magnetica
#endif

          if (passives) then
#ifdef PASSIVES
            !density of neutrals
            u(neqdyn+1,i,j,k) = neutral_h    !/rhosc        ! neutral h
            u(neqdyn+2,i,j,k) = ion_h        !/rhosc           ! ion h
            u(neqdyn+3,i,j,k) = neutral_heS  !/rhosc      ! neutral singlet he
            u(neqdyn+4,i,j,k) = neutral_heM  !/rhosc       ! neutral triplet he
            u(neqdyn+5,i,j,k) = ion_he       !/rhosc          ! ion he
            u(neqdyn+6,i,j,k) = electrons    !/rhosc
            !   passive scalar (tag)
            u(neqdyn+7,i,j,k) = passive      !/rhosc
#endif
          endif
        endif

      end do
    end do
  end do

end subroutine impose_exo

subroutine boundary_exo(u,time)

  use constants, only : pi,Kb
  use globals, only : coords, dx, dy, dz

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: x,y,z
  real :: xp,yp,zp,rp
  real :: xs,ys,zs,rs
  real :: velx,vely,velz,dens
  !real :: bx,by,bz,temp,phi
  real :: neutral_h,ion_h,neutral_heS,neutral_heM,ion_he, electrons, ntot
  real :: passive, nHt, nHet,fion

#ifdef BFIELD
  real :: cpi
#endif
  integer :: i,j,k
  integer :: global_i,global_j,global_k
  logical :: set_cell
  real    :: energy

  do i=nxmin,nxmax
    do j=nymin,nymax
      do k=nzmin,nzmax

        set_cell = .false.

        global_i = i + coords(0)*nx
        global_j = j + coords(1)*ny
        global_k = k + coords(2)*nz

        ! Position measured from the centre of the grid (planet)
        x = (float(global_i - nxtot/2)-0.5)*dx
        y = (float(global_j - nytot/2)-0.5)*dy
        z = (float(global_k - nztot/2)-0.5)*dz

        xp = x - Planet%x
        yp = y - Planet%y
        zp = z - Planet%z

        ! Distance from the centre of the planet
        rp = sqrt(xp**2 + yp**2 + zp**2)*rsc

        nHt   = (1.-ymol)/ (1. + 3.*ymol )  !total hydrogen density
        nHet  =      ymol/ (1. + 3.*ymol )  !total he density

        ! IF INSIDE THE PLANET
        if(rp <= R0) then ! IF INSIDE RMASK

          set_cell = .true.

          velx = 0.0
          vely = 0.0
          velz = 0.0
          dens = rho0

          neutral_h   = (1. - fHII) * nHt * dens/rhosc    !neutral density of H [g/cm3]
          ion_h       = fHII * nHt * dens/rhosc

          neutral_heS = (1. - fHeII) * nHet * dens/rhosc /(1+alpha) !neutral density of He singlet [g/cm3]
          neutral_heM = alpha*neutral_heS     ! neutral density of He triplet [g/cm3]
          ion_he      = fHeII * nHet * dens/rhosc
          electrons   = ion_h + ion_he

          energy = energy_c1*(dens)**gamma /rhosc/vsc2
          passive = 1000*dens/rhosc

          dens = dens/rhosc


        elseif(rp <= Rbound) then ! IF INSIDE RMASK

          set_cell = .true.

          velx = 0.0
          vely = 0.0
          velz = 0.0
          dens = (dens_c1 + dens_c2/rp)**(1.0/(gamma - 1.0))

          neutral_h   = (1. - fHII) * nHt * dens/rhosc     !neutral density of H [g/cm3]
          ion_h       = fHII * nHt * dens / rhosc

          neutral_heS = (1. - fHeII) * nHet * dens/rhosc /(1+alpha)  !neutral density of He singlet [g/cm3]
          neutral_heM = alpha*neutral_heS     ! neutral density of He triplet [g/cm3]
          ion_he      = fHeII * nHet * dens/rhosc
          electrons   = ion_h + ion_he

          energy      = energy_c1*(dens)**gamma /rhosc/vsc2
          passive     = 1000*dens/rhosc
          dens        = dens/rhosc

        endif

         if((global_i <= 1).AND.(u(2,i,j,k) > 0))then
           set_cell = .true.

           velx = 0.0
           vely = 0.0!u(3,i,j,k)*vsc
           velz = 0.0!u(4,i,j,k)*vsc

           dens = rho_out*1.0E-5
           energy = energy_c1*rho_out**gamma /rhosc/vsc2

           neutral_h   = (1. - fHII_s) * nHt * dens/rhosc     !neutral density of H [g/cm3]
           ion_h       = fHII_s * nHt * dens /rhosc

           neutral_heS = (1. - fHeII_s) * nHet * dens/rhosc /(1+alpha)  !neutral density of He singlet [g/cm3]
           neutral_heM = alpha*neutral_heS     ! neutral density of He triplet [g/cm3]
           ion_he      = fHeII_s * nHet * dens/rhosc
           electrons   = ion_h + ion_he

           passive     = -1000*dens/rhosc
           dens        = dens/rhosc

         endif

!         if((global_i >= nxtot).AND.(u(2,i,j,k) > 0))then
!           set_cell = .true.
!
!           velx = 0.0
!           vely = 0.0! u(3,i,j,k)*vsc
!           velz = 0.0!u(4,i,j,k)*vsc
!
!           dens = rho_out*1.0E-5
!           energy = energy_c1*rho_out**gamma /rhosc/vsc2
!
!           neutral_h   = (1. - fHII_s) * nHt * dens/rhosc     !neutral density of H [g/cm3]
!           ion_h       = fHII_s * nHt * dens /rhosc
!
!           neutral_heS = (1. - fHeII_s) * nHet * dens/rhosc /(1+alpha)  !neutral density of He singlet [g/cm3]
!           neutral_heM = alpha*neutral_heS     ! neutral density of He triplet [g/cm3]
!           ion_he      = fHeII_s * nHet * dens/rhosc
!           electrons   = ion_h + ion_he
!
!           passive     = -1000*dens/rhosc
!           dens        = dens/rhosc
!
!         endif
!
         if((global_j <= 1).AND.(u(3,i,j,k) > 0))then
           set_cell = .true.

           velx = 0.0!u(2,i,j,k)*vsc
           vely = 0.0
           velz = 0.0!u(4,i,j,k)*vsc

           dens = rho_out*1.0E-5
           energy = energy_c1*rho_out**gamma /rhosc/vsc2

           neutral_h   = (1. - fHII_s) * nHt * dens/rhosc     !neutral density of H [g/cm3]
           ion_h       = fHII_s * nHt * dens /rhosc

           neutral_heS = (1. - fHeII_s) * nHet * dens/rhosc /(1+alpha)  !neutral density of He singlet [g/cm3]
           neutral_heM = alpha*neutral_heS     ! neutral density of He triplet [g/cm3]
           ion_he      = fHeII_s * nHet * dens/rhosc
           electrons   = ion_h + ion_he

           passive     = -1000*dens/rhosc
           dens        = dens/rhosc

         endif

         if((global_j >= nytot).AND.(u(3,i,j,k) < 0))then
           set_cell = .true.

           velx = 0.0!u(2,i,j,k)*vsc
           vely = 0.0
           velz = 0.0!u(4,i,j,k)*vsc

           dens = rho_out*1.0E-5
           energy = energy_c1*rho_out**gamma /rhosc/vsc2

           neutral_h   = (1. - fHII_s) * nHt * dens/rhosc     !neutral density of H [g/cm3]
           ion_h       = fHII_s * nHt * dens /rhosc

           neutral_heS = (1. - fHeII_s) * nHet * dens/rhosc /(1+alpha)  !neutral density of He singlet [g/cm3]
           neutral_heM = alpha*neutral_heS     ! neutral density of He triplet [g/cm3]
           ion_he      = fHeII_s * nHet * dens/rhosc
           electrons   = ion_h + ion_he

           passive     = -1000*dens/rhosc
           dens        = dens/rhosc

         endif

         if((global_k >= nztot).AND.(u(4,i,j,k) < 0))then
           set_cell = .true.

           velx = 0.0!u(2,i,j,k)*vsc
           vely = 0.0!u(3,i,j,k)*vsc
           velz = 0.0

           dens = rho_out*1.0E-5
           energy = energy_c1*rho_out**gamma /rhosc/vsc2

           neutral_h   = (1. - fHII_s) * nHt * dens/rhosc     !neutral density of H [g/cm3]
           ion_h       = fHII_s * nHt * dens /rhosc

           neutral_heS = (1. - fHeII_s) * nHet * dens/rhosc /(1+alpha)  !neutral density of He singlet [g/cm3]
           neutral_heM = alpha*neutral_heS     ! neutral density of He triplet [g/cm3]
           ion_he      = fHeII_s * nHet * dens/rhosc
           electrons   = ion_h + ion_he

           passive     = -1000*dens/rhosc
           dens        = dens/rhosc

         endif
!
!         if((global_k <= 1).AND.(u(4,i,j,k) < 0))then
!           set_cell = .true.
!
!           velx = 0.0!u(2,i,j,k)*vsc
!           vely = 0.0!u(3,i,j,k)*vsc
!           velz = 0.0
!
!           dens = rho_out*1.0E-5
!           energy = energy_c1*rho_out**gamma /rhosc/vsc2
!
!           neutral_h   = (1. - fHII_s) * nHt * dens/rhosc     !neutral density of H [g/cm3]
!           ion_h       = fHII_s * nHt * dens /rhosc
!
!           neutral_heS = (1. - fHeII_s) * nHet * dens/rhosc /(1+alpha)  !neutral density of He singlet [g/cm3]
!           neutral_heM = alpha*neutral_heS     ! neutral density of He triplet [g/cm3]
!           ion_he      = fHeII_s * nHet * dens/rhosc
!           electrons   = ion_h + ion_he
!
!           passive     = -1000*dens/rhosc
!           dens        = dens/rhosc
!
!         endif
       !Stellar Wind at z= 0 and x=nxtot
        xp = x - barycenter%x
        yp = y - barycenter%y
        zp = z - barycenter%z
        
        xs = x - star%x
        ys = y - star%y
        zs = z - star%z
        rs = sqrt(xs**2+ys**2+zs**2)
        
       if(global_k <= 1 .or. global_i >= nxtot) then
          set_cell = .true.
        
          velx = star%WV*xs/rs -omegap*zp
          vely = star%WV*ys/rs
          velz = star%WV*zs/rs + omegap*xp   !code units
          dens = star%WD/rhosc
        
          neutral_h = (1. - fHII_s) * nHt * dens!/rhosc
          ion_h     = fHII_s * nHt * dens !ion density of H [g/cm3]
        
          neutral_heS = (1.- fHeII_s) * nHet * dens /(1+alpha)
          neutral_heM = alpha * neutral_heS
          ion_he      = fHeII_s * nHet *dens !ion density of He [g/cm3]
          electrons   = ion_h + ion_he
        
        
          ntot = neutral_h + ion_h + neutral_heS + neutral_heM +ion_he
          energy = 0.5*dens*(velx**2 + vely**2 + velz**2)
          energy = energy + cv*(ntot+electrons)*star%WT*Kb/vsc2/rhosc
       
          passive = -1000*dens
          energy = energy

        endif


        if(set_cell)then

          ! dens = dens/rhosc
          ! energy = energy/rhosc/vsc2

          u(1,i,j,k) = dens
          u(2,i,j,k) = dens*velx
          u(3,i,j,k) = dens*vely
          u(4,i,j,k) = dens*velz
          u(5,i,j,k) = energy

#ifdef BFIELD
          u(6,i,j,k) = bx
          u(7,i,j,k) = by
          u(8,i,j,k) = bz
          u(5,i,j,k) = u(5,i,j,k) + 0.5*(bx**2 + by**2 + bz**2)  !Magnetica
#endif

          if (passives) then
#ifdef PASSIVES
            u(neqdyn+1,i,j,k) = neutral_h   !/rhosc    ! neutral h
            u(neqdyn+2,i,j,k) = ion_h       !/rhosc    ! ion h
            u(neqdyn+3,i,j,k) = neutral_heS !/rhosc    ! neutral singlet he
            u(neqdyn+4,i,j,k) = neutral_heM !/rhosc    ! neutral triplet he
            u(neqdyn+5,i,j,k) = ion_he      !/rhosc    ! ion he
            u(neqdyn+6,i,j,k) = electrons   !/rhosc
            u(neqdyn+7,i,j,k) = passive     !/rhosc    !passive scalar
#endif
          endif
        endif

      end do
    end do
  end do

end subroutine boundary_exo

!=======================================================================

subroutine get_central_mass(u, Mass)

  use constants, only : pi
  use globals, only : coords, dx, dy, dz, rank
  implicit none
  real, intent (in) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (out) :: Mass
  real :: x, y, z, rads
  integer ::  i, j, k

  Mass = 0.0

  do k=1,nz
    do j=1,ny
      do i=1,nx

        ! Position measured from the centre of the grid
        x=(real(i+coords(0)*nx-nxtot/2)-0.5)*dx
        y=(real(j+coords(1)*ny-nytot/2)-0.5)*dy
        z=(real(k+coords(2)*nz-nztot/2)-0.5)*dz

        ! Distance from the centre of the first source
        rads=sqrt( (x-planet%x)**2 + (y-planet%y)**2 + (z-planet%z)**2 )

        ! If inside the radius defined
        if( rads <= RMdot ) then

          if(rads == 0.) rads=dx*0.10
          !dens = dw(1)*rw(1)**2/rads**2

          Mass = Mass + u(1,i,j,k)*dx*dy*dz * rhosc *rsc**3

        end if
      end do
    end do
  end do
end subroutine get_central_mass

!======================================================================
!  Write mass acretion rate
subroutine  output_mdot(Mini, Mfin, dt_CFL)
  use parameters,  only : tsc
  use globals,     only : time, currentIteration
  use constants,   only : hr
  implicit none
  real, intent(in )  :: Mini, Mfin, dt_CFL
  real :: mass_loss_rate

  mass_loss_rate = (Mfin - Mini)/(dt_CFL*tsc)

  write(unitout,'(i0,4es15.7)') currentIteration, time*tsc, mass_loss_rate, &
                                dt_CFL*tsc, (Mfin - Mini)
  call flush(unitout)


end subroutine output_mdot
!======================================================================

end module exoplanet

!=======================================================================
