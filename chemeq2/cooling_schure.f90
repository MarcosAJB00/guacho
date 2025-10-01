!=======================================================================
!> @file cooling_schure.f90
!> @brief Cooling module with cooling curves of Schure et al. 2009
!> @author Alejandro Esquivel
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

!> @brief Cooling module with cooling curves of Schure et al. 2009
!> @details Cooling module cooling curves of Schure et al. 2009, A&A, 508,751
!> @n The location of the tables is assumed to be in
!> src/cool_lib/coolingSKKKV.tab

module cooling_schure

  implicit none
  real (kind=8), allocatable :: cooltab(:,:)

  !> The following selects the column of the ionization fraction used for low T
  !> dmc_f = 1 : f_i = 1e-4
  !>         2 : f_i = 1e-3
  !>         3 : f_i = 1e-2
  !>         4 : f_i = 1e-1
  integer, parameter         :: dmc_f = 1


contains

  !=======================================================================
  !> @brief Initializes the DMC cooling
  !> @details Declares variables and reads table
  subroutine init_cooling_schure()

    implicit none

    allocate(cooltab(2,180))
    call read_table_schure()

  end subroutine init_cooling_schure

  !=======================================================================
  !> @brief Reads the cooling curve table
  !> @details Reads the cooling curve table,
  !! the location is assumed in /src/cool_lib/coolingSKKKV.tab
  subroutine read_table_schure()

    use parameters, only : workdir, master
    use globals, only : rank
#ifdef MPIP
    use mpi
#endif
    implicit none
    integer :: i, err
    real (kind=8) :: data(5)  ! Table has five columns


    if(rank == master) then
      open(unit=10,file= trim(workdir)//'../src/cool_lib/coolingSKKKV.tab',    &
           status='old')
!      open(unit=10,file= trim(workdir)//'./schure_tab2.txt',status='old') !111L
!      read(10,*) data2(:)
!      open(unit=11,file= trim(workdir)//'./schure_tab3.txt',status='old') ! 78
!      read(11,*) data3(:)

      do i=1,180
        read(10,*) data(:)
        cooltab(1,i)= 10.0**data(1)
        cooltab(2,i)= 10.0**data(1+dmc_f)
      end do
      close(unit=10)
    endif
#ifdef MPIP
    call mpi_bcast(cooltab,360,mpi_double_precision,0,mpi_comm_world,err)
#endif

end subroutine read_table_schure

  !=======================================================================
  !> @brief Returns the cooling coefficient interpolating the table
  !> @param real [in] T : Temperature K
  function get_lambda(T)

    implicit none
    real , intent(in) :: T
    integer           :: if1
    real, parameter   :: deltaTemp=0.04 !  spacing of T in tables
    real, parameter   :: logTmin = 1.0 !  log of minimum value of temp
    real (kind=8)     :: get_lambda, T0, T1, C0, C1


    if(T.gt.1e8) then
      get_lambda=0.21e-26*sqrt(T)     !free free cooling
    else
      if1=int(( log10(T)- logTmin) /deltaTemp) +1
      if (if1 < 1) then
        get_lambda = cooltab(2,1)
        return
      end if
      T0=cooltab(1,if1)
      c0=cooltab(2,if1)
      T1=cooltab(1,if1+1)
      c1=cooltab(2,if1+1)
      get_lambda=(c1-c0)*(T-T0)/(T1-T0)+c0
    end if

  end function get_lambda

  !=======================================================================
  subroutine cooling_rates_cen(T,ne,nh,nhII,nhe,nheM,nheII,cool_rate)

    implicit none
    real, intent(in)   :: T, ne, nh, nhe,nheM, nheII, nhII
    real, intent (out) :: cool_rate
    real :: cita_h, cita_heS, cita_heII, cita_heM, lya_cool, eta_h, eta_he, &
            w_heII, si_h, si_heI, si_heII, teta

   !collisional ionization cooling (erg/cm3/s)
    !h
    cita_h   = 1.27e-21*sqrt(T)  * exp(-157809.1/T) *ne*nh
    !he I
    cita_heS = 9.38e-22*sqrt(T)  * exp(-285335.4/T) *ne*nhe
    !he II
    cita_heII= 4.95e-22*sqrt(T)  * exp(-631515.0/T) *ne*nheII
    !he (2^3S)
    !cita_heM = 5.01e-27*T**(-0.1687)/ (1+sqrt(T/1e5)) * exp(-55338.0/T) *ne*ne*nheII
    !different from Cen1992 because we included the triplet density explicitly
    cita_heM = 6.41e-21*sqrt(T)  * exp(-55338.0/T) *ne*nheM

    !lyman alpha cooling (Black 1981)
    lya_cool = 7.5e-19*exp(-118348/T) *ne*nh

    ! recombination cooling
    !hII
    eta_h  = 8.70e-27*sqrt(T)*(T/1e3)**(-0.2)/(1 + (T/1e6)**0.7) *ne*nhII
    !heII
    eta_he = 1.55e-26*T**0.3647 *ne*nheII

    !dielectronic recombination cooling (He)
    w_heII = 1.24e-13*T**(-1.5)*exp(-470000/T)*(1 + 0.3*exp(-94000/T)) *ne*nheII

    !colisional excitation cooling
    !h (all n)
    si_h    = 7.5e-19 * exp(-118348/T) * ne*nh
    !he II (n=2)
    si_heII = 5.54e-17*T**(-0.397) * exp(-473638/T) *ne*nheII
    !he I (n=2,3,4 triplets)
    !si_heI  = 9.10e-27*T**(-0.1687) / (1+sqrt(T/1e5)) * exp(-13179/T) *ne*ne*nheII
    si_heI  = 1.16e-20*sqrt(T) * exp(-13179/T) * ne * nheM

    !Bremsstrahlung cooling of all ions:
    !gff=1.5
!    teta = 1.42e-27 * 1.5 *sqrt(T)* (nhII + nheII ) *ne!+ 4*nheIII) * ne
    teta = 1.42e-27 * 1.5 *sqrt(T)* (nhII + nheII ) *ne!+ 4*nheIII) * ne

    cool_rate = cita_h + cita_heS + cita_heII + cita_heM + eta_h + eta_he &
              + w_heII + si_h + si_heII + si_heI + teta + lya_cool
!    cool_rate = cita_h + eta_h !+ si_h  + teta + lya_cool

   end subroutine cooling_rates_cen

  !=======================================================================
  !> @brief High level wrapper to apply cooling with CHIANTI tables
  !> @details High level wrapper to apply cooling with CHIANTI tables
  !> @n cooling is applied in the entire domain and updates both the
  !! conserved and primitive variables
  subroutine coolingschure()

    use parameters, only : nx, ny, nz, cv, Psc, tsc, dif_rad, mhd, n1_chem,   &
                           dif_rad, dif_radHe, nxtot, nytot, nztot, rsc, neqdyn
    use constants,  only : Kb
    use globals,    only : u, primit, dt_CFL, coords, dx, dy, dz
    use hydro_core, only : u2prim
    use difradHe,   only : phHI, phHeIS, phHeIM
    use network,    only : iHI, iHII, iHeIS, iHeIM, iHeII, ie
    use exoplanet,  only : Planet, Rbound

    implicit none
    real                 :: T , dens
    real, parameter      :: Tmin=1000.
    real (kind=8)        :: Lambda0, emtauC, gain, loss
    integer              :: i, j, k, iiHI, iiHeIS, iiHeIM, iiHII, iiHeII, iie
    real                 :: dt_seconds, ch_factor
    real                 :: ne, nh, nhe, nheM, nheII, nhII
    real                 :: radp, xp, yp, zp, xx, yy, zz

    dt_seconds = dt_CFL*tsc
    iiHI   = n1_chem - 1 + iHI    ! neqdyn + 1
    iiHII  = n1_chem - 1 + iHII   ! neqdyn + 2
    iiHeIS = n1_chem - 1 + iHeIS  ! neqdyn + 3
    iiHeIM = n1_chem - 1 + iHeIM  ! neqdyn + 4
    iiHeII = n1_chem - 1 + iHeII  ! neqdyn + 5
    iie    = n1_chem - 1 + ie     ! neqdyn + 6

    do k=1,nz
      do j=1,ny
        do i=1,nx

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
          if( radp > Rbound .and. u(neqdyn+7,i,j,k)>0) then

              !   get the primitives (and T)
              call u2prim(u(:,i,j,k),primit(:,i,j,k),T)

              if(T < 100) then
                print*, 'Temperature, T=', T, primit(1,i,j,k), primit(5,i,j,k)
              end if

              nh    = u(iiHI   ,i,j,k)
              nhII  = u(iiHII  ,i,j,k)
              nhe   = u(iiHeIS ,i,j,k) + u(iiHeIM ,i,j,k)
              nheM  = u(iiHeIM ,i,j,k)
              nheII = u(iiHeII ,i,j,k)
              ne    = u(iie    ,i,j,k)

              if(T > Tmin) then

                !Lambda0=get_lambda(T)  ![erg cm3 /s]
                call cooling_rates_cen(T,ne,nh,nhII,nhe,nheM,nheII,loss) ! [erg/cm3/s]

                if (dif_radHe) then
                  !  energy per photo ionization from solar_spectrum.py (in erg)
                  gain = phHI(i,j,k)   * u(iiHI   ,i,j,k) * 3.8e-12              & ! [erg/cm3/s]
                       + phHeIS(i,j,k) * u(iiHeIS ,i,j,k) * 1.9e-11              &
                       + phHeIM(i,j,k) * u(iiHeIM ,i,j,k) * 9.6e-13

                else

                  gain=0.0

                end if

                dens=primit(1,i,j,k)

                Lambda0 = loss/dens**2 ! [erg cm3 /s]
                !  e^{-dt/Tau}=e^{-2. L0 dt/(3 n K T)}
                emtauC = exp( -2.0*dt_seconds*dens*Lambda0/(3.0*Kb*T) )
                !  this is the Temperature factor of change
                ch_factor = (gain/(dens**2*Lambda0))*(1.0-emtauC)+emtauC

                !  limit changes to avoid catastrophic cooling
                ch_factor = max(ch_factor,0.5)
                ch_factor = min(ch_factor,2.0)

                if ( ch_factor < 0 ) print*, 'cooling factor:', ch_factor

                !  apply cooling to primitive and conserved variables
                primit(5,i,j,k)=primit(5,i,j,k)*ch_factor

                !  update total energy density
                u(5,i,j,k) = cv*primit(5,i,j,k)                                    &
                             + 0.5*primit(1,i,j,k)*(  primit(2,i,j,k)**2           &
                                                    + primit(3,i,j,k)**2           &
                                                    + primit(4,i,j,k)**2  )
#ifdef BFIELD
              if (mhd) then
                u(5,i,j,k) = u(5,i,j,k) + 0.5*(  primit(6,i,j,k)**2                  &
                                               + primit(7,i,j,k)**2                  &
                                               + primit(8,i,j,k)**2  )
              end if
#endif

              end if
          endif
        end do
      end do
    end do


  end subroutine coolingschure

  !======================================================================

end module cooling_schure
