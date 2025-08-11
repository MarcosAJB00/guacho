!=======================================================================
!> @file chemistry.f90
!> @brief chemistry  module
!> @author A. Castellanos, P. Rivera A. Rodriguez, A. Raga  and A. Esquivel
!> @date 10/Mar/2016

! Copyright (c) 2016 A. Esquivel et al.
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

!> @brief chemistry module
!> @details module to solve the chemical/ionic network.
module chemistry

  use network
  implicit none
  integer :: failed_convergence
  integer :: failed_conservation

contains

  !=======================================================================
  !> @brief Advances the chemistry network
  !> @details Advances the chemistry network on the entire domain
  !> (except ghost cells), updates primitives and conserved variables
  !> in globals
  subroutine update_chem()

    use parameters, only : neq, neqdyn, nx, ny, nz, tsc, rhosc, rsc,           &
                           nxtot, nytot, nztot, n_spec, n1_chem
    use globals,    only : u, primit, dt_CFL, coords, dx, dy, dz, rank
    use network,    only : n_elem, iHI, iHII, iHeIS, iHeIM, iHeII, iH, iHe
    use hydro_core, only : u2prim
    use difradHe,   only : phHI, phHeIS, phHeIM
    use exoplanet,  only : Planet,Rbound

    implicit none
    real :: dt_seconds, T, y(n_spec), y0(n_elem)
    integer :: i, j, k, l
    real    :: radp, xp, yp, zp, xx,yy,zz

    dt_seconds = dt_CFL*tsc
    failed_convergence = 0
    failed_conservation = 0

    do k=1,nz
      do j=1,ny
        do i=1,nx

          !   get the primitives (and T)
          call u2prim(u(:,i,j,k),primit(:,i,j,k),T)
          y(1:n_spec) = primit(n1_chem: n1_chem+n_spec-1,i,j,k)
          y0(iH     ) = y(iHI)   + y(iHII)
          y0(iHe    ) = y(iHeIS) + y(iHeIM) + y(iHeII)

          !  update the passive primitives (should not work in single precision)

          ! Position measured from the centre of the grid (planet)
          xx = ( float(i+coords(0)*nx-nxtot/2) - 0.5 )*dx
          yy = ( float(j+coords(1)*ny-nytot/2) - 0.5 )*dy
          zz = ( float(k+coords(2)*nz-nztot/2) - 0.5 )*dz

          xp = xx - Planet%x
          yp = yy - Planet%y
          zp = zz - Planet%z
          ! Distance from the centre of the planet
          radp=sqrt(xp**2+yp**2+zp**2)*rsc

          ! if outside the planet and only planetary material
          if( radp > Rbound .and. u(neqdyn+7,i,j,k)>0) then !0.5*Planet%radius) then

            call chemstep(y, y0, T, dt_seconds, phHI(i,j,k), phHeIS(i,j,k),     &
                                                phHeIM(i,j,k) )

          !  if (failed_conservation > 0) then
          !    u(neqdyn+7,i,j,k) = 222
          !    primit(neqdyn+7,i,j,k) = u(neqdyn+7,i,j,k)
          !  endif
          !  if (failed_convergence > 0 ) then
          !    u(neqdyn+7,i,j,k) = 111
          !    primit(neqdyn+7,i,j,k) = u(neqdyn+7,i,j,k)
          !  endif
          end if

          !  update the primitives and conserved variables
          do l = 0, n_spec-1
            u     (n1_chem+l, i, j, k ) = max( y(l+1), 0. )
            primit(n1_chem+l, i, j, k ) = u(n1_chem+l,i,j,k)
          end do

           ! "correct" the density
          u(1,i,j,k)      = ( y(iHI) + y(iHII) )                               &
                          + 4*( y(iHeIS) + y(iHeIM) + y(iHeII) )
          primit(1,i,j,k) = u(1,i,j,k)

        end do
      end do
    end do

    !if (failed_convergence > 0) write(*, '(a,i3,a,i0,a)') 'in rank: ',rank,    &
    !', chemistry convergence failed in ',failed_convergence,'cells'

    !if (failed_conservation > 0) write(*, '(a,i3,a,i0,a)') 'in rank: ',rank,    &
    !', density conservation failed in ',failed_conservation,'cells'

  end subroutine update_chem

  !=======================================================================
  !> @brief Advances the chemistry network in one cell
  !> @details Advances the chemistry network on the in one cell
  !> @param real [inout] y(n_spec) : number densities of the species
  !> to be updated by the chemistry
  !> @param real [in] y[n_elem] : total number density of each of the
  !> elements involved in the reactions
  !> @param real [in] T : Temperature [K]
  !> @param real [in] deltt : time interval (from the hydro, in seconds)
  subroutine chemstep(y, y0, T, deltt, phHI, phHeIS, phHeIM)

    use linear_system
    use parameters, only : n_spec
    use network,    only : n_reac, n_elem, get_reaction_rates,                 &
                           derv, get_jacobian, n_nequ, check_no_conservation
    use globals, only:rank

    implicit none
    real (kind=8), intent(inout) :: y(n_spec)
    real (kind=8), intent(in) :: y0(n_elem), T, deltt, phHI, phHeIS, phHeIM
    real (kind=8) :: dtm
    real (kind=8) :: y1(n_spec),yin(n_spec), y0_in(n_elem), yt(n_spec)
    real (kind=8) :: rate(n_reac),dydt(n_spec),jacobian(n_spec,n_spec)
    integer, parameter  :: niter=100     ! number of iterations
    real,    parameter  :: tol = 0.005   ! convergence tolerance
    integer :: n,i
    logical :: failed_cons

    n=0
    dtm=1./deltt
    yin(:)   = y(:)
    y0_in(:) = y0(:)

    call get_reaction_rates(rate, T, phHI, phHeIS, phHeIM)
    failed_cons = .false.

    do while ( n <= niter )

      !  initial guess for Newton-Raphson

      !if ( check_no_conservation(y,y0_in) ) then
      !  failed_cons = .true.
      !  call nr_init(y,y0_in)
      !  print('y0_in'= y0_in)
      !end if

      call derv(y,rate,dydt,y0)
      call get_jacobian(y,jacobian,rate)

      do i=1,n_nequ
        jacobian(i,i)= jacobian(i,i) - dtm
        dydt(i)      = dydt(i) - ( y(i)-yin(i) )*dtm
      end do
      y1(:) = -dydt(:)

      call linsys(jacobian, y1, n_spec)
      y(:) = y(:) + y1(:)

      !  exit the loop if converged
      yt(:)=y1(:)/( y(:) + 1e-10 )
      if(all(abs(yt(:)) <= tol )) then
        y(:)=max( y(:), 1.e-10 )
        exit
      end if

      n=n+1

    end do

    if (failed_cons) then
      failed_conservation = failed_conservation + 1
    endif

    if (n >= niter) then
      failed_convergence = failed_convergence + 1
    end if

    return

  end subroutine chemstep

  !=======================================================================

  end module chemistry
