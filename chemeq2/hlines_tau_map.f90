!> @brief tau_utilities
!> @details Utilities to compute the Lyman-@f and H \alpha @f$ absorption

module tau_utilities

use exoplanet, only : init_exo, barycenter, omegap

#ifdef MPIP
  integer :: nps,N_X,N_Y,N_Z
#endif

contains

!> @brief Initializes data
!> @details Initializes data, MPI and other stuff

subroutine init_tau()

!  Initializes MPI, data arrays, etc
use parameters
use globals, only : u, dx, dy, dz, coords, rank, left, right, &
                    top, bottom, out, in, rank, comm3d

implicit none

#ifdef MPIP
  integer :: err
  integer, dimension(0:ndim-1) :: dims
  logical, dimension(0:ndim-1) :: period
  logical :: perx=.false., pery=.false., perz=.false.
#endif

  !initializes MPI
#ifdef MPIP

  if ((bc_left   == BC_PERIODIC).and.(bc_right == BC_PERIODIC)) perx = .true.
  if ((bc_bottom == BC_PERIODIC).and.(bc_top   == BC_PERIODIC)) pery = .true.
  if ((bc_out    == BC_PERIODIC).and.(bc_in    == BC_PERIODIC)) perz = .true.

  period(0) = perx
  period(1) = pery
  period(2) = perz

  call mpi_init (err)
  call mpi_comm_rank (mpi_comm_world,rank,err)
  call mpi_comm_size (mpi_comm_world,nps,err)

  !if (nps.ne.np) then
    ! print*, 'processor number (',nps,') is not equal to pre-defined number (',np,')'
    ! call mpi_finalize(err)
    ! stop
  !endif

  N_X = 1
  N_Y = 1
  N_Z = 1

  dims(0) = MPI_NBX/N_X
  dims(1) = MPI_NBY/N_Y
  dims(2) = MPI_NBZ/N_Z

  if(rank.eq.master) then
    print*,dims(0),dims(1),dims(2)
  endif

#else
  rank = 0
  coords(:) = 0
#endif

  if(rank.eq.master) then
     print '(a)' ,"*******************************************"
     print '(a)' ,"                        _                 *"
     print '(a)' ,"  __   _   _  __ _  ___| |__   ___    3   *"
     print '(a)' ," / _ `| | | |/ _` |/ __| '_ \ / _ \    D  *"
     print '(a)' ,"| (_| | |_| | (_| | (__| | | | (_) |      *"
     print '(a)' ," \__, |\__,_|\__,_|\___|_| |_|\___/       *"
     print '(a)' ," |___/                                    *"
  endif

#ifdef MPIP
  if(rank.eq.master) then
     print '(a,i3,a)','*    running with mpi in', nps , ' processors    *'
     print '(a)' ,'*******************************************'
     print '(a)', 'Calculating hydrogen lines tau map'
  end if
  call mpi_cart_create(mpi_comm_world,ndim,dims,period,.true.,comm3d,err)
  call mpi_comm_rank(comm3d, rank, err)
  call mpi_cart_coords(comm3d, rank, ndim, coords, err)
  print '(a,i3,a,3i4)', 'processor ', rank                              &
       ,' ready w/coords',coords(0),coords(1),coords(2)
  !call mpi_cart_shift(comm3d, 0, 1, left  , right, err)
  !call mpi_cart_shift(comm3d, 1, 1, bottom, top  , err)
  !call mpi_cart_shift(comm3d, 2, 1, out   , in   , err)
  call mpi_barrier(mpi_comm_world, err)
  !
#else
  print '(a)' ,'*******************************************'
  print '(a)' ,'*     running on a single processor       *'
  print '(a)' ,'*******************************************'
  print '(a)', 'Calculating Lyman and H Alpha Tau'
#endif

  !   grid spacing
  dx = xmax/nxtot
  dy = ymax/nytot
  dz = zmax/nztot

  !   allocate big arrays in memory
  allocate(u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )

  call init_exo()

end subroutine init_tau

!=======================================================================

!> @brief reads data from file
!> @details reads data from file
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param integer [in] itprint : number of output
!> @param string [in] filepath : path where the output is

subroutine read_data(u,itprint,filepath)

  use parameters, only : np, neq, nxmin, nxmax, nymin, nymax, &
                         nzmin, nzmax, MPI_NBY, MPI_NBZ
  use globals, only : rank, comm3d, coords
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: itprint
  character (len=128), intent(in) :: filepath
  integer :: unitin, ip, err
  character (len=128) file1
  character           :: byte_read
  character, parameter  :: lf = char(10)
  integer :: nxp, nyp, nzp, x0p, y0p, z0p,&
             mpi_xp, mpi_yp, mpi_zp,neqp,neqdynp, nghostp, num
  real :: dxp, dyp, dzp, scal(3), cvp

  take_turns: do ip = 0,nps-1
    if (rank == ip) then
#ifdef MPIP
      num = (coords(0)*MPI_NBY + coords(1))*MPI_NBZ + coords(2)
      !print*,rank,coords(0),coords(1),coords(2),num
      write(file1,'(a,i3.3,a,i3.3,a)')  &
            trim(filepath)//'BIN/points',num,'.',itprint,'.bin'
      unitin = rank + 10
#else
      write(file1,'(a,i3.3,a)')         &
            trim(filepath)//'BIN/points',itprint,'.bin'
      unitin=10
#endif
      open(unit=unitin,file=file1,status='unknown',access='stream', &
           convert='LITTLE_ENDIAN')

      !   discard the ascii header
      do while (byte_read /= achar(255) )
        read(unitin) byte_read
        !print*, byte_read
      end do

      !  read bin header, sanity check to do
      read(unitin) byte_read
      read(unitin) byte_read
      read(unitin) nxp, nyp, nzp
      read(unitin) dxp, dyp, dzp
      read(unitin) x0p, y0p, z0p
      read(unitin) mpi_xp, mpi_yp, mpi_zp
      read(unitin) neqp, neqdynp
      read(unitin) nghostp
      read(unitin) scal(1:3)
      read(unitin) cvp
      read(unitin) u(:,:,:,:)
      close(unitin)

      print'(i3,a,a)',rank,' read file:',trim(file1)

#ifdef MPIP
      call mpi_barrier(comm3d,err)
#endif
    end if

  end do take_turns

end subroutine read_data

!=======================================================================
!> @brief gets position of a cell
!> @details Returns the position and spherical radius calculated with
!! respect to  the center of the grid
!> @param integer [in] i : cell index in the x direction
!> @param integer [in] j : cell index in the y direction
!> @param integer [in] k : cell index in the z direction
!> @param real    [in] x : x position in the grid
!> @param real    [in] y : y position in the grid
!> @param real    [in] z : z position in the grid
subroutine getXYZ(i,j,k,x,y,z)
  use globals,    only : dx, dy, dz, coords
  use parameters, only : nx, ny, nz, nxtot, nytot, nztot
  implicit none
  integer, intent(in)  :: i, j, k
  real,    intent(out) :: x, y, z
  !planet centered position
  x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
  y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
  z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz
end subroutine getXYZ

!=======================================================================
!> @brief Rotation around the X axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid
subroutine rotation_x(theta,x,y,z,xn,yn,zn)
  ! rotation around the x axis by an angle theta
  implicit none
  real, intent(in ) :: theta, x, y, z
  reaL, intent(out) :: xn, yn, zn
  xn =   x
  yn =   y*cos(theta) - z*sin(theta)
  zn =   y*sin(theta) + z*cos(theta)
end subroutine rotation_x

!=======================================================================

!> @brief Rotation around the Y axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid
subroutine rotation_y(theta,x,y,z,xn,yn,zn)
  implicit none
  real, intent(in ) :: theta, x, y, z
  real, intent(out) :: xn, yn, zn
  xn =   x*cos(theta) + z*sin(theta)
  yn =   y
  zn = - x*sin(theta) + z*cos(theta)
end subroutine rotation_y

!=======================================================================
!> @brief Rotation around the Z axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid
subroutine rotation_z(theta,x,y,z,xn,yn,zn)
  implicit none
  real, intent(in ) :: theta, x, y, z
  real, intent(out) :: xn, yn, zn
  xn =   x*cos(theta) - y*sin(theta)
  yn =   x*sin(theta) + y*cos(theta)
  zn =   z
end subroutine rotation_z

!=======================================================================
!> @brief Fill target map
!> @details Fills the target map of one MPI block
!> @param integer [in] nxmap : Number of X cells in target
!> @param integer [in] nymap : Number of Y cells in target
!> @param real [in] u(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax) :
!! conserved variables
!> @param real [out] map(nxmap,mymap) : Target map
!> @param real [in] dxT : target pixel width
!> @param real [in] dyT : target pixel height
!> @param real [in] thetax : Rotation around X
!> @param real [in] thetay : Rotation around Y
!> @param real [in] thetaz : Rotation around Z
subroutine fill_map(nxmap,nymap,nvmap,vmin,vmax,u,map_ha,map_lya,dxT,dyT,&
                   theta_x,theta_y,theta_z)!,n2p_matrix)!,n2_save)

  use constants, only : clight
  use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         neq, nx, ny, nz, vsc, rsc,nztot, nxtot,nytot,neqdyn
  use globals, only : dz, coords
  use hydro_core, only : u2prim

  implicit none

  integer, intent(in) :: nxmap,nymap,nvmap
  real, intent(in)  :: vmin, vmax
  real, intent(inout)  :: u(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax)
  real , intent(in) :: dxT, dyT, theta_x, theta_y, theta_z
!  real, intent(in)  :: n2p_matrix(200,200)
!  real, intent(inout)  :: n2_save(nxtot,nytot,nztot )
  real, intent(out) :: map_ha(nxmap,nymap,nvmap)
  real, intent(out) :: map_lya(nxmap,nymap,nvmap)
  integer :: i,j,k, iobs, jobs
  real :: x,y,z,xn,yn,zn, vx, vy, vz,vxn, vyn, vzn
  real :: T, prim(neq), profile(nvmap),n2pp
  real, parameter :: sigmaHA = 0.017014, lambdaHA = 6.56279e-5 !(c/nu0=lambda)
  real, parameter :: sigmaLA = 0.01105 , lambdaLA = 1.215668e-5 !(c/nu0=lambda)

  do k=1,nz
     do j=1,ny
        do i=1,nx

          !  obtain original position
          call getXYZ(i,j,k, x,y,z)

          !convert to stellar reference frame
          x = x - barycenter%x
          y = y - barycenter%y
          z = z - barycenter%z

          !  do the rotation of the coordinates
          call rotation_y(theta_y,x,y,z,xn,yn,zn)
          call rotation_x(theta_x,xn,yn,zn,x,y,z)
          call rotation_z(theta_z,x,y,z,xn,yn,zn)

          ! This is the position on the target (centered)
          ! Integration is along Z
          iobs = xn/dxT + nxmap/2
          jobs = yn/dyT + nymap/2

          !  make sure the result lies in the map bounds
          if( (iobs >= 1    ).and.(jobs >= 1    ).and. &
              (iobs <= nxmap).and.(jobs <= nymap)) then

            !only do to the side facing the observer
            if(zn > 0) then

              !  get the velocity in cm/s
              call u2prim(u(:,i,j,k),prim,T)

              vx = prim(2)
              vy = prim(3)
              vz = prim(4)

              !convert to inercial frame
              vx = vx + omegap*z
              vz = vz - omegap*x   !cm/s

              !  obtain the LOS velocity
              call rotation_y(theta_y,vx,vy,vz,vxn,vyn,vzn)
              call rotation_x(theta_x,vxn,vyn,vzn,vx,vy,vz)
              call rotation_z(theta_z,vx,vy,vz,vxn,vyn,vzn)

              vxn = vxn*vsc
              vyn = vyn*vsc
              vzn = vzn*vsc

              ! call subroutine to interpolate the value in the n2 level
              !call calc_pop((prim(1)-prim(neqdyn+1)),T,n2p_matrix,n2pp)
!              call calc_pop(prim(neqdyn+6),T,n2p_matrix,n2pp)  ! Sep 4th, 2024

              !n2_save(i+coords(0)*nx,j+coords(1)*ny,k+coords(2)*nz)=n2pp*prim(neqdyn+1)

              !  calculate the line profile function
              !call phigauss(T, -vzn,vmin,vmax,nvmap,profile)
              !call phivoigt(T, -vzn,vmin,vmax,nvmap,profile,lambdaHA, 4.4101e7)
               !if (prim(neqdyn+2)<0) then
!                    call phivoigt(T, -vzn,vmin,vmax,nvmap,profile,lambdaHA, 4.4101e7)
!                    map_ha(iobs,jobs,:)= map_ha(iobs,jobs,:) + &
!                                 dz*rsc*n2pp*prim(neqdyn+1)*sigmaHA*lambdaHA*profile(:)

                    call phivoigt(T, -vzn,vmin,vmax,nvmap,profile,lambdaLA, 6.27e8)
                    map_lya(iobs,jobs,:)= map_lya(iobs,jobs,:) + &
!                                 dz*rsc*(1.-n2pp)*prim(neqdyn+1)*sigmaLA*lambdaLA*profile(:)
                                 dz*rsc*prim(neqdyn+1)*sigmaLA*lambdaLA*profile(:)
               !end if
            end if
          end if
      end do
    end do
  end do

end subroutine fill_map

!=======================================================================

!> @brief Writes projection to file
!> @details Writes projection to file
!> @param integer [in] itprint : number of output
!> @param string [in] filepath : path where to write
!> @param integer [in] nxmap : Number of X cells in target
!> @param integer [in] nymap : Number of Y cells in target
!> @param integer [in] nvmap : Number of velocity channels
!> @param real [in] map(nxmap,mymap) : Target map f.read_reals(kind).reshape(nxmap,nymap,nvmap,or

subroutine  write_maps(itprint,i,filepath,nxmap,nymap,nvmap,map_ha,map_lya)

  use parameters, only: nxtot,nytot,nztot

  implicit none
  integer, intent(in) :: nxmap, nymap,nvmap,itprint,i
  character (len=128), intent(in) :: filepath
  real, intent(in)    :: map_ha(nxmap,nymap,nvmap)
  real, intent(in)    :: map_lya(nxmap,nymap,nvmap)
!  real, intent(in)    :: n2_save(nxtot,nytot,nztot)
  character (len=128) file1
  character (len=128) file2
  integer ::  unitout

!  write(file1,'(a,i3.3,a,i3.3,a)')  trim(filepath)//'hlines/HA_tau-',itprint,'_',i,'.bin'
!  unitout=11
!  open(unit=unitout,file=file1,status='unknown',form='unformatted', &
!       convert='LITTLE_ENDIAN')

!  write (unitout) map_ha(:,:,:)
!  close(11)
!  print'(a,a)'," wrote file:",trim(file1)

  write(file2,'(a,i3.3,a,i3.3,a)')  trim(filepath)//'hlines/LA_tau-',itprint,'_',i,'.bin'
  open(unit=12,file=file2,status='unknown',form='unformatted', &
      convert='LITTLE_ENDIAN')

  write (12) map_lya(:,:,:)
  close(12)
  print'(a,a)'," wrote file:",trim(file2)

end subroutine write_maps

!=======================================================================

!> @brief This routine computes a gaussian line profile
!> @details This routine computes a gaussian line profile

subroutine phigauss(T,vzn,vmin,vmax,nvmap,profile)

  use constants, only: amh, pi, kB, clight
  implicit none
  real, intent(in) :: T, vzn, vmin, vmax
  integer, intent(in) :: nvmap
  real, intent(out) :: profile(nvmap)
  integer :: i
  real :: coef, dv, vr

  profile(:)=0.
  dv=(vmax-vmin)/float(nvmap)

  coef=amh/(2.*kB*T)

  do i=1,nvmap
     vr=(float(i)-0.5)*dv+vmin
     profile(i)=sqrt(coef/pi)*exp(-coef*((vr-vzn)**2) )
  end do

end subroutine phigauss

!=======================================================================

!> @brief This routine computes a voigt line profile
!> @details This routine computes a voigt line profile

subroutine phivoigt(T,vzn,vmin,vmax,nvmap,profile, wvl, Ga)
  use constants, only: amh, pi, kB, clight
  implicit none
  real, intent(in) :: T, vzn,vmin, vmax, wvl, Ga
  integer, intent(in) :: nvmap
  real, intent(out) :: profile(nvmap)
  integer :: i
  real :: sigm,gamm,dv,vr,ccs,a0,a,v,coef
  !real, parameter :: lambdaHA=6.56279e-5 !(c/nu0=lambda)
  complex (kind = 8) :: W

  a0 = 0.0529e-6 !Radio de Bohr [cm]
  ccs = pi*(2.0*a0)**2

  profile(:) = 0.
  dv = (vmax-vmin)/float(nvmap)

  sigm = sqrt(2.0*kB*T/amh)
  gamm = Ga*wvl!*rho_n*ccs/(2.0*pi)*sigm
  coef = 1./sqrt(pi)/sigm

  do i=1,nvmap
    vr=(float(i)-0.5)*dv + vmin
    v = (vr - vzn)/sigm
    a = gamm/4.0/pi/sigm
    call humlicek(a,v,W)
    profile(i) = coef*REAL(W)
  end do
end subroutine phivoigt

subroutine humlicek(a,v,W)
  implicit none
  real    (kind = 8), intent(in)  :: a,v
  complex (kind = 8), intent(out) :: W
  complex (kind = 8)              :: z, u
  real    (kind = 8)              :: s

  z = cmplx(a, -v)
  s = abs(v) + a

  if (s >= 15.0) then
    !* --- Approximation in region I --        -------------- *!
    W = (z * 0.5641896) / (0.5 + (z * z))
  else if (s >= 5.5) then
    !* --- Approximation in region II --       -------------- *!
    u = z * z
    W = (z * (1.410474 + u*0.5641896)) / (0.75 + (u*(3.0 + u)))
  else if (a >= 0.195*abs(v) - 0.176) then
    !* --- Approximation in region III --      -------------- *!
    W = (16.4955 + z*(20.20933 + z*(11.96482 + z*(3.778987 + &
         0.5642236*z)))) / &
         (16.4955 + z*(38.82363 + z*(39.27121 + z*(21.69274 + &
         z*(6.699398 + z)))))
  else
    !* --- Approximation in region IV --       -------------- *!
    u = z * z
    W = exp(u) - (z*(36183.31 - u*(3321.99 - u*(1540.787 - &
         u*(219.031 - u*(35.7668 - u*(1.320522 - u*0.56419)))))) / &
         (32066.6 - u*(24322.84 - u*(9022.228 - u*(2186.181 - &
         u*(364.2191 - u*(61.57037 - u*(1.841439 - u))))))))
  endif
end subroutine humlicek

subroutine read_tab(t_array,ne_array, n2p_matrix)
  implicit none

  real*8, allocatable,intent(out) :: t_array(:),ne_array(:),n2p_matrix(:,:)
  integer :: i

! read table generated with chiantipy
  open(1,file='./H2pop_chianti_lya_T8150.bin',status='unknown',access='stream')

  read(1) i
  !print*,i

  allocate(t_array(i),ne_array(i),n2p_matrix(i,i))

  read(1) ne_array
  read(1) t_array
  read(1) n2p_matrix
  close(1)
end subroutine read_tab

subroutine calc_pop(ne,te,n2p_matrix,n2pp)

  implicit none
  real, intent(in) :: ne, te, n2p_matrix(200,200)
  real, intent(out):: n2pp
  integer :: ite, ine
  integer, parameter :: N=200
  real,parameter ::  T_low=1e2, T_high=1e7, N_low=1e1, N_high=1e9
  real :: n0, n1, t0, t1, c0, c1, c2, c3, f, g

! find the location of the table corresponding to the cell value te and ne
! ex for te: T_low*(T_high/T_low)**((i-1)/(N-1)) = Te
  if ((te >= T_high) .or. (te <= T_low) ) then
    print*,'te out of range'
    stop
  endif
  ite = int(1 + (N-1)*log10(te/T_low)/log10(T_high/T_low))
  if(ite >200) print*,'ite is out of range'
!same for ne
  if ((ne >= N_high) .or. (ne <= N_low) ) then
    !print*,'ne out of range',ne,N_high,N_low
    n2pp = 0.0
    return
  endif
  ine = int(1 + (N-1)*log10(ne/N_low)/log10(N_high/N_low))
  if(ine >200) print*,'ine is out of range'
! get the interpolated value of n2 pop from the cell values
!if T is above the ionization value for n2p what happend???

  t0 = T_low*(T_high/T_low)**((ite-1.)/N)
  t1 = T_low*(T_high/T_low)**((ite    )/N)
  n0 = N_low*(N_high/N_low)**((ine-1.)/N)
  n1 = N_low*(N_high/N_low)**((ine    )/N)

  c0 = n2p_matrix(ine  , ite)
  c1 = n2p_matrix(ine+1, ite)
  c2 = n2p_matrix(ine  , ite+1)
  c3 = n2p_matrix(ine+1, ite+1)

  f = log10(te/t0)/log10(t1/t0)
  g = log10(ne/n0)/log10(n1/n0)
  n2pp = c0*(c1/c0)**f*(c2/c0)**g*(c3*c0/(c1*c2))**(f*g)
  ! if (n2pp<1e-30) then
  !   print*,'n2p=',n2pp
  !   print*,'ne=',ne,'te=',te
  !   print*,'ine=',ine,'ite=',ite
  ! endif
end subroutine calc_pop

!=======================================================================

end module tau_utilities

!=======================================================================

!> @brief Computes the Ly-alpha apbsorption
!> @details Computes the Ly-alpha apbsorption
!! @n It rotates the data along each of the coordinates axis
!! by an amount @f$ \theta_x, \theta_y, \theta_z @f$, and the LOS
!! is along the Z axis

program hlines_tau_map

  use constants,  only : pi
  use parameters, only : xmax,master, mpi_real_kind,np, &
                         nxtot,nytot,nztot,outputpath
  use globals, only : u, rank, comm3d,coords
  use tau_utilities

  implicit none

#ifdef MPIP
  include "mpif.h"
#endif

  character (len=128) :: filepath
  integer :: err
  integer :: itprint
  !
  real, parameter :: theta_x = 0.0
  real            :: theta_y, x, thetay
  real, parameter :: theta_z = 0.00
  ! thetay array dimensions and values
  integer, parameter :: nsteps=80
  real, parameter    :: thetay_max=40, thetay_min=-40
  integer            :: step, i
  !   map and its dimensions
  integer, parameter :: nvmap=250
  integer            :: nxmap, nymap
  real               :: dxT, dyT, vmin,vmax
  real, allocatable  :: map_ha(:,:,:), map1(:,:,:)
  real, allocatable  :: map_lya(:,:,:), map2(:,:,:)
!  real               :: n2_save(nxtot,nytot,nztot)
  integer            :: baston,ierr,stat(MPI_STATUS_SIZE)
  real,allocatable   :: n2p_matrix(:,:), te_array(:), ne_array(:)
  integer            :: cx, cy, cz
  integer :: isnap_x,isnap_y,isnap_z

  ! read table of values generated with chianti
!  call read_tab(te_array,ne_array,n2p_matrix)

  nxmap = nxtot !100 !nxtot/2 + 100
  nymap = nytot !100 ! nytot

!  allocate(map_ha(nxmap,nymap,nvmap) ,map1(nxmap,nymap,nvmap))
  allocate(map_lya(nxmap,nymap,nvmap),map2(nxmap,nymap,nvmap))

  ! initializes program
  call init_tau()

  cx = coords(0)
  cy = coords(1)
  cz = coords(2)

  !  minumim and maximum of the velocity channels
  vmin = -300.0e5
  vmax =  300.0e5

  !  Target pixel size, relative to the simulation
  dxT = xmax/float(nxtot)
  dyT = dxT

  ! chose output (fix later to input form screen)
  filepath = trim(outputpath)
  itprint = 100

  loop_over_theta_y : do i=0,nsteps

    theta_y = (thetay_max - thetay_min)*real(i)/real(nsteps) + thetay_min
    theta_y = theta_y*pi/180.0

    !  read ph and u from file
   ! if(rank > 0)then
   !   call MPI_RECV(baston,1,MPI_INT,rank-1,66,MPI_COMM_WORLD,stat,ierr)
   ! end if
   ! call read_data(u,itprint,filepath)
   ! if(rank < np-1)then
   !   call MPI_SEND(baston,1,MPI_INT,rank+1,66,MPI_COMM_WORLD,ierr)
   ! end if

    !  resets map
!    map_ha(:,:,:)  = 0.0
!    map1(:,:,:)    = 0.0
    map_lya(:,:,:) = 0.0
    map2(:,:,:)    = 0.0

    do isnap_x = 0, N_X-1
      do isnap_y = 0, N_Y-1
       do isnap_z = 0, N_Z-1

         coords(0) = cx*N_X + isnap_x
         coords(1) = cy*N_Y + isnap_y
         coords(2) = cz*N_Z + isnap_z

        call read_data(u,itprint,filepath)

        !if (rank == master) then
        !   print'(a)', 'Calculating projection with angles of rotaton'
        !   print'(f6.2,a,f6.2,a,f6.2,a)',theta_x*180./pi,'° around X, ' &
        !                                ,theta_y*180./pi,'° around Y, '&
        !                                ,theta_z*180./pi,'° around Z, '
        !end if

        !  add info to the map
        call fill_map(nxmap,nymap,nvmap,vmin,vmax,u,map_ha,map_lya,dxT,dyT,theta_x,theta_y,theta_z)!,n2p_matrix)!, n2_save )
       end do
     end do
    end do

    !  sum all the partial sums
!     call mpi_reduce(map_ha,map1,nxmap*nymap*nvmap, mpi_real_kind, mpi_sum, master, comm3d, err)
     call mpi_reduce(map_lya,map2,nxmap*nymap*nvmap, mpi_real_kind, mpi_sum, master, comm3d, err)

    !  write result
    if (rank == master) then
      call write_maps(itprint,i,filepath,nxmap,nymap,nvmap,map1,map2)!,n2_save)
    end if

  end do loop_over_theta_y

  if (rank == master) print*, 'my work here is done, have a  nice day'
#ifdef MPIP
  call mpi_finalize(err)
#endif
  !
  stop

end program hlines_tau_map

!=======================================================================
