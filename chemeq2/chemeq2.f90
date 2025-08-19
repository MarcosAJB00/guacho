!=======================================================================
!> @file chemistry.f90
!> @brief chemistry  module
!> @author Marcos Baracchi, C. Villareal and A. Esquivel
!> @date 11/Ago/2025

! Copyright (c) 2024 Guacho Co Op.
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
!> @details module to solve the chemical/ionic network following an
!> adaptation of the method in Mott D.R. & Oran E. 2001
!> (CHEMEQ2: A Solver for the Stiff Ordinary Differential Equations of Chemical
!?> Kinetics)
module chemeq2

  use parameters, only : n => n_spec
  use network

  implicit none
  real    :: epsmin, sqreps, epscl, epsmax, dtmin, tstart
  integer :: itermax
  real    :: ymin(n)

contains

  !=======================================================================
  !> @brief Updates the chemistry network
  !> @details Advances the chemistry network with the chemeq2 algorithm
  subroutine update_chemeq2()

    use parameters, only: neq, nx, ny,nz, tsc, n1_chem
    use globals,    only : u, primit, dt_CFL, coords, dx, dy, dz, rank

    implicit none
    real    :: dt_seconds, T, y(n)
    integer :: i, j, k

    dt_seconds = dt_CFL*tsc

    !  loop over the entire domain
    do k=1,nz
      do j=1,ny
        do i=1,nx

          ! for the exoplanet we can restrict here not to do the chemistry
          ! inside the planetary boundary conditions

          !   get the primitives (and T)
          call u2prim(u(:,i,j,k),primit(:,i,j,k),T)

          !  copy densities to y vector
          y(1:n) = primit(n1_chem: n1_chem + n -1 )

          !  advances the y's
          call chemeq2solve( dt_seconds, y, T )

        end do
      end do
    end do

    !  Aditional conservation check/correction can/should be made here

  end subroutine update_chemeq2

  !=======================================================================
  !> @brief Main driver for the module
  !> @details Advances the rate equations by a time increment dtg
  subroutine chemeq2solve(dtg, y, T)

    implicit none
    real,    intent(in)    :: dtg, T
    real,    intent(inout) :: y(n)
    integer :: i, gcount
    integer :: iter
    real    :: ts, tn, tfd
    real    :: ymn(n)
    real    :: q(n), d(n), rtaus(n), y1(n)
    real    :: alpha, qs(n)
    real    :: ys(n), y0(n), rtau(n)
    real    :: scr1, scr2, scrarray(n)
    real    :: dt, dto
    real    :: rswitch
    real    :: scrtch, ascr, eps
    real    :: rtaui, rtaub, qt, pb, rteps
    !! ym1, ym2 and stab are only used for the stability check on dt
    !! real    :: ym1(n), ym2(n), stab

    tn = 0.0

    !  store and limit to 'ymin' the initial values
    q(:)  = 0.0
    p(:)  = 0.0
    y0(:) = y(:)
    do i = 1, n
      y(i)  = max( y(i), ymin(i) )
    end do

    !  obtain the derivatives of the initial values
    !call gsub(y, q, d, tn + tstart)
    gcount = gcount + 1

    ! Estimate the initial stepsize
    ! Strongly increasing functions (q >> d assumed here) use a stepsize
    ! proportional to the step neede fpr the function to reach equilibrium
    ! where as functions decreasing or in equilibrium use a stepsize
    ! proportional to the characteristic steosize of the function.
    ! Comnvergence of the integration scheme is likely since the smallest
    ! estimate is chosen for the initial stepsize.

    ! initial stepsize
    scrtch = 1.0e-25
    do i = 1, n
      ascr   = abs( q(i) )
      scr2   = sign( 1.0/y(i), 0.1*epsmin*ascr - d(i) )
      scr1   = scr2 * d(i)
      scrtch = max( scr1, -abs( ascr - d(i) ) * scr2, scrtch )
    end do
    dt = min( sqreps/scrtch, dtg )

    ! The starting values are stored
    100 continue
    ts = tn

    do i = 1,  n
      rtau(i)  = dt*d(i)/y(i)
      ys(i)    = y(i)
      qs(i)    = q(i)
      rtaus(i) = rtau(i)
    end do

    ! find the predictor terms

    101 continue

    do i = 1, n

      !  prediction
      rtaui = rtau(i)

      !  note the one of two approximations for alpha is chosen:
      !  1) Pade b for all rtaui (see supporting NRL memo report)
      !     or
      !  2) Pade for rtaui <= rswitch (see supporting NRL memo report
      !     (Mott et al. 2002))

      ! Option 1): Pade b
      alpha = (180.0 + rtaui*(60.0 + rtaui*(11.0 + rtaui) ) ) &
            / (360.0 + rtaui*(60.0 + rtaui*(12.0 + rtaui) ) )

      ! Option 2): Pade a
      !  if ( rtaui <= rswitch) then
      !    alpha = ( 840.0 + rtaui*(140.0 + rtaui*(20.0 + rtaui) ) ) &
      !          / (1680.0 + 40.0 * rtaui**2)
      !  else
      !    alpha = 1.0 - 1.0/rtaui
      !  end if

      scrarray(i) = (q (i) - d(i) ) / (1.0 + alpha*rtaui )

    end do

    iter = 1
    do while(iter <= itermax)

      !  limit decreasing functions to their minimum values
      do i = 1, n
        !!ym2(i) = ym1(i)
        !!ym1(i) = y(i)
        y(i) = max( ys(i)+ dt*scrarray(i), ymin(i) )
      end do

      if ( iter == 1 ) then
        ! The first corrector step advances the time (tentatively) and saves the
        ! initial predictor value as y1 for the timestep check later
        tn = ts + dt
        y1(:) = y(:)
      end if

      ! evaluate derivatives for the corrector
      ! call gsub(y, q, d, tn + tstart)
      gcount = gcount + 1
      eps    = 1.0e-10

      do i = 1, n

        rtaub = 0.5*( rtaus(i) + dt*d(i)/y(i) )

        ! Same options for calculating alpha as in predictor
        ! Option 1): Pade b
        alpha = (180.0 + rtaub*(60.0 + rtaub*(11.0 + rtaub) ) ) &
              / (360.0 + rtaub*(60.0 + rtaub*(12.0 + rtaub) ) )

        ! Option 2): Pade a
        !  if ( rtaub <= rswitch) then
        !    alpha = ( 840.0 + rtaub*(140.0 + rtaub*(20.0 + rtaub) ) ) &
        !          / (1680.0 + 40.0 * rtaub**2)
        !  else
        !    alpha = 1.0 - 1.0/rtaub
        !  end if

        qt = qs(i)*( 1.0 - alpha ) + q(i)*alpha
        pb = rtaub/dt
        scrarray(i) = ( qt - ys(i)*pb ) / ( 1.0 - alpha*rtaub )

      end do

      iter = iter + 1

    end do

    ! Calculate new f, check for convergence, and limit decreasing functions.
    ! The order of the operations in this loop is important
    do i = 1, n
      scr2 = max( ys(i) + dt*scrarray(i), 0.0 )
      scr1 = abs( scr2 - y1(i) )
      y(i) = max( scr2, ymin(i) )
      !!ym2(i) = ym1(i)
      !!ym1(i) = y(i)

      if( 0.25*(ys(i) + y(i)) >  ymin(i) ) then
        scr1 = scr1/y(i)
        eps  = max(0.5*(scr1+min(abs(q(i)-d(i))/(q(i)+d(i)+1e-30),scr1)),eps)
      end if
    end do

    eps = eps*epscl

    !  print out diagnostics if stepsize becomes to small
    if (dt < dtmin + 1.0e-16*tn) then
      print*, 'chemeq error; stepsize too small'

      ! call error diagnostic routine
      ! call chemer()
    end if

    ! Check for convergence

    ! The following section is used for the stability check
    stab = 0.01
    if (itermax >= 3) then
      do i = 1, n
        stab = max( stab, abs(y(i)-ym(i))/( abs(ym1(i)-ym2(i))+1.0e-20*y(i) ) )
      end do
    end if

    if  (eps <= epsmax ) then
    !!if ( (eps <= epsmax) .and. (stab <= 1) ) then
      ! Valid step. Return if dtg has been reached
      if (dtg <= tn*tfd ) return
    else
      ! Invalid step; reset tn to ts
      tn =ts
    end if

    ! Perform stepsize modifications
    ! estimate sqrt(eps) by Newton iteration
    rteps = 0.5*(   eps +    1.0    )
    rteps = 0.5*( rteps + eps/rteps )
    rteps = 0.5*( rteps + eps/rteps )

    dto = dt
    dt  = min( dt*(1.0/rteps + 0.005), tfd*( dtg-tn ) )
    !!dt  = min( dt*(1.0/rteps + 0.05), tfd*( dtg-tn ), dto/(stab+0.001) )

    ! begin new stwp if previous dot converged
    if (eps > epsmax) then
    !!if (eps > epsmax .or. stab > 1) then
      rcount = rcount + 1

      ! After an unsuccessful steo the initial timescales don't change, but dt
      ! doesm requiring rtaus to be scaled by the ratio of the new and old
      ! timesteps
      dto = dt/dto
      rtaus(:) = rtaus(:)*dto

      ! unsuccessful steps return to line 101 so that the initial source terms
      ! do not get recalculated
      go to 101
    end if

    ! Successful step; get the source terms for the next step and continue back
    ! to 100
    !call gsub(y, q, d, tn + tstart)
    gcount = gcount + 1
    go to 100

    return

    end subroutine chemeq2solve

  !=======================================================================
  !> @brief Parameter setup
  !> @details Resets local parameters if their value is > 0, otherwise the
  !> defaults values are returned
  !> @param real epsmn : maximum relative error, default value is 1.0e-2
  !> @param real epsmx : this number provides the basis for deciding whether
  !>                     convergence can be aqchieved without step size
  !>                     reduction. if (eps/epsmin > epsmx) further reduction
  !>                     is applied. Default value: 10.0
  !> @param real dtmn  : smallest stepsize allowed. Default: 1.0e-15
  !> @param real tnot  : initial value of time
  !> @param real ymn(n): floor value for y
  !> @param integer itermx : number of times the corrector is applied
  !>                         Default = 1
  subroutine chemsp(epsmn, epsmx, dtmn, tnot, ymn, itermx)

    implicit none
    real,    intent(in) :: epsmn, epsmx, dtmn, tnot, ymn(n)
    integer, intent(in) :: itermx
    integer :; i

    epsmin = 1.0e02
    if (epsmn > 0.0) then
      epsmin = epsmn
      sqreps = 5.0 * sqrt(epsmin)
    end if

    epscl  = 1.0/epsmin
    epsmax = 10.0
    if (epsmx > 0) epsmax = epsmx

    dtmin = 1.0e-15
    if (dtm > 0) dtmin = dtm
    tstart = tnot

    itermax = 1
    if (itermx > 0) itermax = itermx

    do i=1,n
      ymn(i) = 1.0e-20
      if ( ymn(i) > 0.0 ) ymin(i) = ymn(i)
    end do

  end subroutine chemsp
  !=======================================================================
end module chemeq2
