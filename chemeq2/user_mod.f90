!=======================================================================
!> @file user_mod.f90
!> @brief User input module
!> @author C. Villarreal, M. Schneiter, A. Esquivel
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
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief User imput module
!> @details  This is an attempt to have all input neede from user in a
!! single file
!!!  This module should load additional modules (i.e. star, jet, sn), to
!!  impose initial and boundary conditions (such as sources)

module user_mod

  ! load auxiliary modules
  use exoplanet

  implicit none

contains

!> @brief Initializes variables in the module, as well as other
!! modules loaded by user.
!! @n It has to be present, even if empty
subroutine init_user_mod()

  implicit none
  !  initialize modules loaded by user
  call init_exo()

end subroutine init_user_mod

!=====================================================================

!> @brief Here the domain is initialized at t=0
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param real [in] time : time in the simulation (code units)

subroutine initial_conditions(u)

  use parameters, only : neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         passives, neqdyn
!  use globals,    only: coords, dx ,dy ,dz

  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)

  call impose_exo(u,0.0)

end subroutine initial_conditions

!=====================================================================

!> @brief User Defined Boundary conditions
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param real [in] time : time in the simulation (code units)
!> @param integer [in] order : order (mum of cells to be filled in case
!> domain boundaries are being set)

subroutine impose_user_bc(u,order)

  use parameters, only:  neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals,    only: time
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: order

  !  In this case the boundary is the same for 1st and second order)
  if (order >= 1) then
    call boundary_exo(u,time)
  end if

end subroutine impose_user_bc

!=======================================================================

!> @brief User Defined source terms
!> This is a generic interrface to add a source term S in the equation
!> of the form:  dU/dt+dF/dx+dG/dy+dH/dz=S
!> @param real [in] pp(neq) : vector of primitive variables
!> @param real [inout] s(neq) : vector with source terms, has to add to
!>  whatever is there, as other modules can add their own sources
!> @param integer [in] i : cell index in the X direction
!> @param integer [in] j : cell index in the Y direction
!> @param integer [in] k : cell index in the Z direction

subroutine get_user_source_terms(pp, s, i, j , k)

  use constants,  only : Ggrav
  use parameters, only : nx, ny, nz, nxtot, nytot, nztot, rsc, vsc2, neqdyn,beta_pressure
  use globals,    only : dx, dy, dz, coords
  use exoplanet
  use radpress

  implicit none
  integer, intent(in) :: i,j,k
  real, intent(in)    :: pp(neq)
  real, intent(inout) :: s(neq)
  integer, parameter  :: nb=3
  real    :: x(nb),y(nb),z(nb), GM(nb), rad2(nb)
  integer :: index
  real    :: xc ,yc, zc
  real    :: GradPhi(nb), OmegaSq
  real    :: rsoft
  real    :: v, fracv, frac_neutro

  rsoft = (dx*0.1)**2

  GM(2) = Ggrav*Star%mass/rsc/vsc2
  GM(1) = Ggrav*Planet%mass/rsc/vsc2

  !   get cell position
  xc = (float(i + coords(0)*nx - nxtot/2) - 0.5)*dx
  yc = (float(j + coords(1)*ny - nytot/2) - 0.5)*dy
  zc = (float(k + coords(2)*nz - nztot/2) - 0.5)*dz

  ! calculate distance from the sources
  ! planet
  x(1) = xc - Planet%x
  y(1) = yc - Planet%y
  z(1) = zc - Planet%z
  rad2(1) = x(1)**2 + y(1)**2 + z(1)**2

  if(rad2(1) < rsoft)then
    rad2(1) = rsoft
  endif

  ! star
  x(2) = xc - Star%x
  y(2) = yc - Star%y
  z(2) = zc - Star%z
  rad2(2) = x(2)**2 + y(2)**2 + z(2)**2

  if(rad2(2) < rsoft)then
    rad2(2) = rsoft
  endif

  ! barycenter
  x(3) = xc - Barycenter%x
  y(3) = 0.0 ! porque queremos la distancia en el plano orbital x,z
  z(3) = zc - Barycenter%z
  rad2(3) = x(3)**2 + y(3)**2 + z(3)**2
  if(rad2(3) < rsoft)then
    rad2(3) = rsoft
  endif

  if ( beta_pressure ) then
      beta(i,j,k) = 0.
      !  do only outside the planet
      if( rad2(1) >= planet%radius**2 ) then

        !Each cell feels a given pressure proportional to the neutrals fraction
        frac_neutro = pp(neqdyn+1)/pp(1)
        !  Radial velocity in km s^-1 at the stellar frame
        v = (((pp(2)+omegap*rorb)*x(2) + pp(3)*y(2) + pp(4)*z(2))/sqrt(rad2(2)))* (vsc)

        fracv = (v-vr(1))/(vr(Nr)-vr(1))*Nr
        index = int(fracv)+1

        if (index < 1) then
          index = 1
        else if ( index > Nr-1 ) then
          index = Nr-1
        end if
        !Linear interpolation for Beta
        Beta(i,j,k) = (Br(index) + (v-vr(index))*(Br(index+1)-Br(index))/(vr(index+1)-vr(index)))*frac_neutro
        !Update scale factor GM
        GM(2)=GM(2)*(1.-Beta(i,j,k))

      end if
    endif


  OmegaSq =  ( GM(2) + GM(1) )/rorb**3

  GradPhi(1) = GM(1)*x(1)/rad2(1)**1.5 + GM(2)*x(2)/rad2(2)**1.5 - OmegaSq*x(3)
  GradPhi(2) = GM(1)*y(1)/rad2(1)**1.5 + GM(2)*y(2)/rad2(2)**1.5 - OmegaSq*y(3)
  GradPhi(3) = GM(1)*z(1)/rad2(1)**1.5 + GM(2)*z(2)/rad2(2)**1.5 - OmegaSq*z(3)

  !  update source terms with gravity
  s(2)= s(2) - pp(1)*GradPhi(1) - 2*pp(1) * sqrt(OmegaSq) * pp(4)
  s(3)= s(3) - pp(1)*GradPhi(2) ! 0
  s(4)= s(4) - pp(1)*GradPhi(3) + 2*pp(1) * sqrt(OmegaSq) * pp(2)

  ! energy
  s(5)= s(5) - pp(1)*(pp(2)*GradPhi(1) + pp(3)*GradPhi(2) + pp(4)*GradPhi(3))

end subroutine get_user_source_terms


!=======================================================================

end module user_mod

!=======================================================================
