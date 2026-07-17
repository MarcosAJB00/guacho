!> @brief tau_utilities
!> @details Utilities to compute the He 10830 (3P-2 metastable) absorption
!!          map, leyendo directamente las primitivas (densidad de He*,
!!          temperatura, rapidez radial) generadas a partir de un modelo
!!          1D esferizado, sin pasar por variables conservadas ni EOS.

module tau_utilities

use exoplanet, only : init_exo, barycenter, omegap

contains

!> @brief Initializes data (modo serie, sin MPI)
subroutine init_tau()

  use parameters
  use globals, only : dx, dy, dz

  implicit none

  print '(a)' ,"*******************************************"
  print '(a)' ,"                        _                 *"
  print '(a)' ,"  __   _   _  __ _  ___| |__   ___    3   *"
  print '(a)' ," / _ `| | | |/ _` |/ __| '_ \ / _ \    D  *"
  print '(a)' ,"| (_| | |_| | (_| | (__| | | | (_) |      *"
  print '(a)' ," \__, |\__,_|\__,_|\___|_| |_|\___/       *"
  print '(a)' ," |___/                                    *"
  print '(a)' ,'*******************************************'
  print '(a)' ,'*     running on a single processor       *'
  print '(a)' ,'*******************************************'
  print '(a)', 'Calculating Helium 10830 Tau Map (primitive input)'

  !   grid spacing
  dx = xmax/nxtot
  dy = ymax/nytot
  dz = zmax/nztot

  call init_exo()

end subroutine init_tau

!=======================================================================

!> @brief Lee las tres primitivas desde archivos planos (sin header)
!> @details Cada archivo es un dump binario crudo de nx*ny*nz real(kind=8)es,
!!          en orden Fortran (columna-mayor, i varía más rápido).
!!          densidad: n_He* en cm^-3. temperatura: K. velocidad
!!          radial en cm/s (escalar, se vectoriza en fill_map_prim).

subroutine read_data_prim(dens_he, temp, vmod, filepath)

  use parameters, only : nx, ny, nz

  implicit none
  real(kind=8), intent(out) :: dens_he(nx,ny,nz), temp(nx,ny,nz), vmod(nx,ny,nz)
  character (len=128), intent(in) :: filepath
  integer :: unitin

  open(newunit=unitin, file=trim(filepath)//'density.bin', status='old', &
       access='stream', convert='LITTLE_ENDIAN')
  read(unitin) dens_he(:,:,:)
  close(unitin)
  print '(a)', ' read file: density.bin'

  open(newunit=unitin, file=trim(filepath)//'temperature.bin', status='old', &
       access='stream', convert='LITTLE_ENDIAN')
  read(unitin) temp(:,:,:)
  close(unitin)
  print '(a)', ' read file: temperature.bin'

  open(newunit=unitin, file=trim(filepath)//'velocity.bin', status='old', &
       access='stream', convert='LITTLE_ENDIAN')
  read(unitin) vmod(:,:,:)
  close(unitin)
  print '(a)', ' read file: velocity.bin'

end subroutine read_data_prim

!=======================================================================
!> @brief gets position of a cell (centrada en el planeta = 0,0,0 de la caja)
subroutine getXYZ(i,j,k,x,y,z)
  use globals,    only : dx, dy, dz
  use parameters, only : nxtot, nytot, nztot
  implicit none
  integer, intent(in)  :: i, j, k
  real(kind=8),    intent(out) :: x, y, z
  x=(float(i-nxtot/2)+0.5)*dx
  y=(float(j-nytot/2)+0.5)*dy
  z=(float(k-nztot/2)+0.5)*dz
end subroutine getXYZ

!=======================================================================
subroutine rotation_x(theta,x,y,z,xn,yn,zn)
  implicit none
  real(kind=8), intent(in ) :: theta, x, y, z
  real(kind=8), intent(out) :: xn, yn, zn
  xn =   x
  yn =   y*cos(theta) - z*sin(theta)
  zn =   y*sin(theta) + z*cos(theta)
end subroutine rotation_x

!=======================================================================
subroutine rotation_y(theta,x,y,z,xn,yn,zn)
  implicit none
  real(kind=8), intent(in ) :: theta, x, y, z
  real(kind=8), intent(out) :: xn, yn, zn
  xn =   x*cos(theta) + z*sin(theta)
  yn =   y
  zn = - x*sin(theta) + z*cos(theta)
end subroutine rotation_y

!=======================================================================
subroutine rotation_z(theta,x,y,z,xn,yn,zn)
  implicit none
  real(kind=8), intent(in ) :: theta, x, y, z
  real(kind=8), intent(out) :: xn, yn, zn
  xn =   x*cos(theta) - y*sin(theta)
  yn =   x*sin(theta) + y*cos(theta)
  zn =   z
end subroutine rotation_z

!=======================================================================
!> @brief Fill target map usando primitivas directamente (sin u2prim)
subroutine fill_map_prim(nxmap,nymap,nvmap,vmin,vmax,dens_he,temp,vmod,map_he,&
                          dxT,dyT,theta_x,theta_y,theta_z)

  use constants,  only : clight, pi, emass, echarge
  use parameters, only : nx, ny, nz, vsc, rsc, tsc
  use globals,    only : dz

  implicit none

  integer, intent(in) :: nxmap,nymap,nvmap
  real(kind=8),    intent(in) :: vmin, vmax
  real(kind=8),    intent(in) :: dens_he(nx,ny,nz), temp(nx,ny,nz), vmod(nx,ny,nz)
  real(kind=8),    intent(in) :: dxT, dyT, theta_x, theta_y, theta_z
  real(kind=8),   intent(out) :: map_he(nxmap,nymap,nvmap)
  integer :: i,j,k, iobs, jobs
  real(kind=8)    :: x,y,z,xn,yn,zn, vx, vy, vz, vxn, vyn, vzn
  real(kind=8)    :: xp,yp,zp, rad
  real(kind=8)    :: prof_3pj0(nvmap),prof_3pj1(nvmap), prof_3pj2(nvmap), cte
  real(kind=8), parameter :: Temp_floor = 10.0
  real(kind=8) :: Tcell, dcell

  real(kind=8), parameter :: f3pj0 = 5.9902e-2, lambda3pj0 = 1082.909e-7, Aij= 1.0216e7 !1/s
  real(kind=8), parameter :: f3pj1 = 1.7974e-1, lambda3pj1 = 1083.025e-7 !cm
  real(kind=8), parameter :: f3pj2 = 2.9958e-1, lambda3pj2 = 1083.034e-7
  real(kind=8)            :: sigma3pj0, sigma3pj1, sigma3pj2
  real(kind=8)            :: shift_1, shift_2

  cte = pi*echarge**2/emass/clight
  sigma3pj0 = cte*f3pj0 !0.0015940154597622735
  sigma3pj1 = cte*f3pj1 !0.004782951132477565
  sigma3pj2 = cte*f3pj2 !0.007971940025968781

  ! shift of the lines with respect to the j=0 transition
  shift_1 = clight*(lambda3pj2 - lambda3pj0)/lambda3pj2
  shift_2 = clight*(lambda3pj2 - lambda3pj1)/lambda3pj2

  do k=1,nz
     do j=1,ny
        do i=1,nx

          !  posición centrada en el planeta (origen de la caja)
          call getXYZ(i,j,k, x,y,z)

          !  guardamos esta posición planet-centered: es la que define
          !  la dirección radial del viento esferizado (1D -> 3D)
          xp = x
          yp = y
          zp = z

          !  convertimos a marco estelar para la proyección/rotación
          x = x - barycenter%x
          y = y - barycenter%y
          z = z - barycenter%z

          call rotation_y(theta_y,x,y,z,xn,yn,zn)
          call rotation_x(theta_x,xn,yn,zn,x,y,z)
          call rotation_z(theta_z,x,y,z,xn,yn,zn)

          iobs = xn/dxT + nxmap/2
          jobs = yn/dyT + nymap/2

          if( (iobs >= 1    ).and.(jobs >= 1    ).and. &
              (iobs <= nxmap).and.(jobs <= nymap)) then

            if(zn > 0) then

              !  vector velocidad: viento radial centrado en el planeta,
              !  usando vmod(i,j,k) [cm/s] y el versor (xp,yp,zp)/rad
              rad = sqrt(xp*xp + yp*yp + zp*zp)
              if (rad > 0.0) then
                vx = (vmod(i,j,k)*xp/rad) / vsc
                vy = (vmod(i,j,k)*yp/rad) / vsc
                vz = (vmod(i,j,k)*zp/rad) / vsc
              else
                vx = 0.0 ; vy = 0.0 ; vz = 0.0
              end if

              ! convertir al marco inercial
              vx = vx + omegap*z !tsc agragado por mi
              vz = vz - omegap*x !tsc agragado por mi

              call rotation_y(theta_y,vx,vy,vz,vxn,vyn,vzn)
              call rotation_x(theta_x,vxn,vyn,vzn,vx,vy,vz)
              call rotation_z(theta_z,vx,vy,vz,vxn,vyn,vzn)

              vxn = vxn*vsc
              vyn = vyn*vsc
              vzn = vzn*vsc

              Tcell = max(temp(i,j,k), Temp_floor)
              dcell = max(dens_he(i,j,k), 0.0) !antes 1e-30

              !call phivoigt(Tcell,-vzn,vmin,vmax,nvmap,prof_3pj0,lambda3pj0,Aij,0.0)!shift_1)
              call phivoigt(Tcell,-vzn,vmin,vmax,nvmap,prof_3pj1,lambda3pj1,Aij,0.0)!shift_2) 
              !call phivoigt(Tcell,-vzn,vmin,vmax,nvmap,prof_3pj2,lambda3pj2,Aij,0.0)

              map_he(iobs,jobs,:) = map_he(iobs,jobs,:) + dz*rsc*dcell* &
                                    !( sigma3pj0*lambda3pj0*prof_3pj0(:))! + &
                                    ( sigma3pj1*lambda3pj1*prof_3pj1(:) )!+ &
                                    !( sigma3pj2*lambda3pj2*prof_3pj2(:) )

            end if
          end if

        end do
     end do
  end do

end subroutine fill_map_prim

!=======================================================================
!> @brief Writes projection to file (un archivo por paso de rotación i)
subroutine write_maps(i,filepath,nxmap,nymap,nvmap,map_he)

  implicit none
  integer, intent(in) :: nxmap, nymap, nvmap, i
  real(kind=8), intent(in)    :: map_he(nxmap,nymap,nvmap)
  character (len=128), intent(in) :: filepath
  character (len=128) :: file1
  integer ::  unitout

  write(file1,'(a,i3.3,a)') trim(filepath)//'helines/HE_tau_3pj1_',i,'.bin'
  unitout=11
  open(unit=unitout,file=file1,status='unknown',form='unformatted', &
       convert='LITTLE_ENDIAN')
  write (unitout) map_he(:,:,:)
  close(unitout)
  print '(a,a)'," wrote file:",trim(file1)

end subroutine write_maps

!=======================================================================
subroutine phigauss(T,vzn,vmin,vmax,nvmap,profile)
  use constants, only: amh, pi, kB, clight
  implicit none
  real(kind=8), intent(in) :: T, vzn, vmin, vmax
  integer, intent(in) :: nvmap
  real(kind=8), intent(out) :: profile(nvmap)
  integer :: i
  real(kind=8) :: coef, dv, vr

  profile(:)=0.
  dv=(vmax-vmin)/float(nvmap)
  coef=amh/(2.*kB*T)

  do i=1,nvmap
     vr=(float(i)-0.5)*dv+vmin
     profile(i)=sqrt(coef/pi)*exp(-coef*((vr-vzn)**2) )
  end do

end subroutine phigauss

!=======================================================================
subroutine phivoigt(T,vzn,vmin,vmax,nvmap,profile, wvl, Ga, shift)
  use constants, only: amh, pi, kB, clight
  implicit none
  real(kind=8), intent(in)    :: T, vzn,vmin, vmax, wvl, Ga,shift
  real(kind=8), intent(out)   :: profile(nvmap)
  integer, intent(in) :: nvmap
  integer :: i
  real(kind=8)    :: sigm,gamm,dv,vr,a,v,coef
  complex (kind = 8) :: W

  profile(:) = 0.
  dv = (vmax-vmin)/float(nvmap)
  sigm = sqrt(2.0*kB*T/(4*amh))
  gamm = Ga*wvl
  coef = 1./sqrt(pi)/sigm

  do i=1,nvmap
    vr=(float(i)-0.5)*dv + vmin
    v = (vr - vzn + shift)/sigm
    a = gamm/4.0/pi/sigm
    call humlicek(a,v,W)
    profile(i) = coef*real(W)
  end do
end subroutine phivoigt

subroutine humlicek(a,v,W)
  implicit none
  real(kind=8)    , intent(in)  :: a,v
  complex (kind = 8), intent(out) :: W
  complex (kind = 8)              :: z, u
  real(kind=8)                 :: s

  z = cmplx(a, -v)
  s = abs(v) + a

  if (s >= 15.0) then
    W = (z * 0.5641896) / (0.5 + (z * z))
  else if (s >= 5.5) then
    u = z * z
    W = (z * (1.410474 + u*0.5641896)) / (0.75 + (u*(3.0 + u)))
  else if (a >= 0.195*abs(v) - 0.176) then
    W = (16.4955 + z*(20.20933 + z*(11.96482 + z*(3.778987 + &
         0.5642236*z)))) / &
         (16.4955 + z*(38.82363 + z*(39.27121 + z*(21.69274 + &
         z*(6.699398 + z)))))
  else
    u = z * z
    W = exp(u) - (z*(36183.31 - u*(3321.99 - u*(1540.787 - &
         u*(219.031 - u*(35.7668 - u*(1.320522 - u*0.56419)))))) / &
         (32066.6 - u*(24322.84 - u*(9022.228 - u*(2186.181 - &
         u*(364.2191 - u*(61.57037 - u*(1.841439 - u))))))))
  endif
end subroutine humlicek

!=======================================================================

end module tau_utilities

!=======================================================================

!> @brief Computes the He 10830 absorption a lo largo de una secuencia
!! de fases orbitales (rotación en theta_y), a partir de un único
!! conjunto de primitivas (densidad de He*, T, rapidez radial).

program helines_tau_map

  use constants,  only : pi
  use parameters, only : xmax, nx, ny, nz, nxtot, nytot
  use tau_utilities

  implicit none

  character (len=128) :: filepath

  real(kind=8), parameter    :: theta_x = 0.0
  real(kind=8)               :: theta_y
  real(kind=8), parameter    :: theta_z = 0.00
  integer, parameter :: nsteps=4
  real(kind=8), parameter    :: thetay_max=5.0, thetay_min=-5.0
  integer            :: i

  integer, parameter :: nvmap=250
  integer            :: nxmap, nymap
  real(kind=8)               :: dxT, dyT, vmin, vmax

  real(kind=8), allocatable  :: dens_he(:,:,:), temp(:,:,:), vmod(:,:,:)
  real(kind=8), allocatable  :: map_he(:,:,:)

  nxmap = nxtot
  nymap = nytot

  allocate(dens_he(nx,ny,nz), temp(nx,ny,nz), vmod(nx,ny,nz))
  allocate(map_he(nxmap,nymap,nvmap))

  call init_tau()

  vmin = -25.0e5
  vmax =  25.0e5

  dxT = xmax/float(nxtot)
  dyT = dxT

  filepath = './'

  call read_data_prim(dens_he, temp, vmod, filepath)

  loop_over_theta_y : do i=0,nsteps

    theta_y = (thetay_max - thetay_min)*real(i)/real(nsteps) + thetay_min
    theta_y = theta_y*pi/180.0

    map_he(:,:,:) = 0.0

    call fill_map_prim(nxmap,nymap,nvmap,vmin,vmax,dens_he,temp,vmod,map_he,&
                        dxT,dyT,theta_x,theta_y,theta_z)

    call write_maps(i,filepath,nxmap,nymap,nvmap,map_he)

  end do loop_over_theta_y

  print*, 'my work here is done, have a nice day'

  stop

end program helines_tau_map

!=======================================================================