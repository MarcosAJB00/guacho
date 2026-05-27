program test_1D

  use omp_lib
  use solver
  implicit none

  character(len=50)  :: label
  character(len=200) :: perfil_densidad
  character(len=200) :: perfil_temperatura

  real(kind=8) :: T, Rp, time, dtime
  real(kind=8) :: FluxHI, FluxHeIS, FluxHeIM

  integer, parameter :: ns = 6
  integer :: i, j, ios, n_points, total_lines, itermx

  real(kind=8) :: ymn(ns), ti, tf, y(ns)
  real(kind=8) :: epsmn, epsmx, dtmn, tnot, prt

  real(kind=8), allocatable :: r(:), rho(:), T_profile(:)
  real(kind=8), allocatable :: y0_HI(:), y0_HII(:), y0_HeIS(:)
  real(kind=8), allocatable :: y0_HeIM(:), y0_HeII(:), y0_e(:)
  real(kind=8), allocatable :: yi_all(:,:)

  character(len=100) :: filename

  real(kind=8) :: mp, R_jupiter, r_planet, dr
  real(kind=8) :: T_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc
  real(kind=8) :: tau_global(3), tau_loc(3), y0(3)
  real(kind=8) :: t1, t2
  real(kind=8) :: r_dummy
  real(kind=8) :: nHt, nHet, alpha, ymol, fHII, fHeII, n_tot

  character(len=200) :: line

  ! Arrays para guardar tau y phi precalculados para cada celda
  ! (permiten separar la barrera serial de taus del loop paralelo de quimica)
  real(kind=8), allocatable :: tau_all(:,:)   ! (n_points, 3)
  real(kind=8), allocatable :: phi_all(:,:)   ! (n_points, 3)

  ! Array local por thread para las especies durante la integracion
  real(kind=8), allocatable :: y_local(:,:)   ! (n_points, ns)

! =======================================================================
! LEER INPUTS
! =======================================================================
  open(99,file='./inputs/inputs.dat',status='old',action='read')

  read(99,*) label, Rp
  read(99,*) label, FluxHI
  read(99,*) label, FluxHeIS
  read(99,*) label, FluxHeIM
  read(99,*) label, ymol
  read(99,*) label, fHII
  read(99,*) label, fHeII
  read(99,*) label, alpha
  read(99,*) label, dtime
  read(99,*) label, time
  read(99,*) label, perfil_densidad
  read(99,*) label, perfil_temperatura

  close(99)

  print *, 'Inputs:'
  print *, 'Rp = ', Rp, ' R_jupiter'
  print *, 'FluxHI = ', FluxHI, ' cm^-2 s^-1'
  print *, 'FluxHeIS = ', FluxHeIS, ' cm^-2 s^-1'
  print *, 'FluxHeIM = ', FluxHeIM, ' cm^-2 s^-1'
  print *, 'He/H ratio = ', ymol
  print *, 'H ionization fraction = ', fHII
  print *, 'He ionization fraction = ', fHeII
  print *, 'He metastable fraction = ', alpha
  print *, 'dtime = ', dtime, ' s'
  print *, 'time = ', time, ' s'
  print *, 'Threads disponibles: ', omp_get_max_threads()

  call cpu_time(t1)

! =======================================================================
! CONTAR LINEAS PERFIL DENSIDAD
! =======================================================================
  open(10,file='./inputs/'//trim(perfil_densidad),status='old',action='read')
  print *, 'Reading density profile from: ', trim(perfil_densidad)

  total_lines = 0
  do
    read(10,'(A)',iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) /= '#') total_lines = total_lines + 1
  end do

  close(10)

  n_points = total_lines

  allocate(r(n_points), rho(n_points), T_profile(n_points))
  allocate(y0_HI(n_points), y0_HII(n_points))
  allocate(y0_HeIS(n_points), y0_HeIM(n_points))
  allocate(y0_HeII(n_points), y0_e(n_points))
  allocate(yi_all(n_points, ns))
  allocate(tau_all(n_points, 3))
  allocate(phi_all(n_points, 3))
  allocate(y_local(n_points, ns))

! =======================================================================
! LEER PERFIL DE DENSIDAD
! =======================================================================
  open(10,file='./inputs/'//trim(perfil_densidad),status='old',action='read')

  j = 0
  do
    read(10,'(A)',iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) == '#') cycle
    j = j + 1
    read(line,*) r(j), rho(j)
  end do

  close(10)

! =======================================================================
! LEER PERFIL DE TEMPERATURA
! =======================================================================
  open(11,file='./inputs/'//trim(perfil_temperatura), &
      status='old',action='read')

  j = 0
  do
    read(11,'(A)',iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) == '#') cycle
    j = j + 1
    read(line,*) r_dummy, T_profile(j)
  end do
  close(11)

  mp        = 1.67d-24
  R_jupiter = 7.1492d9
  r_planet  = Rp * R_jupiter
  dr        = r_planet*(r(n_points)-r(1))/real(n_points,8)

  print *, 'dr = ', dr, ' cm'
  print *, 'r_planet = ', r_planet, ' cm'
  call flush(6)

! =======================================================================
! PERFILES INICIALES  (igual que el original)
! =======================================================================
  open(30,file='./output/perfiles_iniciales.txt',status='replace')
  write(30,*) '# j y0_HI y0_HII y0_HeIS y0_HeIM y0_HeII y0_e'

  do j = 1, n_points
    n_tot = rho(j)/mp

    nHt  = (1d0 - ymol)/(1d0 + 3d0*ymol)
    nHet = ymol/(1d0 + 3d0*ymol)

    y0_HI(j)  = (1d0 - fHII) * nHt * n_tot
    y0_HII(j) = fHII * nHt * n_tot

    y0_HeIS(j) = (1d0 - fHeII) * nHet * n_tot / (1d0 + alpha)
    y0_HeIM(j) = alpha * y0_HeIS(j)
    y0_HeII(j) = fHeII * nHet * n_tot

    y0_e(j) = y0_HII(j) + y0_HeII(j)

    write(30,'(I6,1X,6ES16.8)') j, &
        y0_HI(j), y0_HII(j), &
        y0_HeIS(j), y0_HeIM(j), &
        y0_HeII(j), y0_e(j)
  end do

  close(30)

! =======================================================================
! CONDICION INICIAL (t = 0)
!
! FISICA: la radiacion viene de afuera, por eso se recorre j=1..n_points
! con i = n_points-j+1 (de afuera hacia adentro) acumulando tau_global.
! Esto es serial por definicion: tau de celda i depende de todas las
! celdas mas externas. NO se puede paralelizar esta barrida.
! =======================================================================
  t = 0.0d0
  tau_global = (/ 0.d0, 0.d0, 0.d0 /)

  ! --- Barrida serial de taus (afuera -> adentro) ---
  do j = 1, n_points
    i = n_points - j + 1

    y0 = (/ y0_HI(i), y0_HeIS(i), y0_HeIM(i) /)
    call local_taus(y0, dr, tau_global, tau_loc)
    tau_global = tau_loc

    ! Guardar tau y phi de cada celda para usarlos luego
    tau_all(j,:) = tau_loc
    call phis_rate(FluxHI, FluxHeIS, FluxHeIM, tau_loc, &
                   phi_all(j,1), phi_all(j,2), phi_all(j,3))
  end do

  ! --- Escritura condicion inicial (serial, una sola vez) ---
  do j = 1, n_points
    i = n_points - j + 1

    write(filename,'(A,I0,A)') './output/species/cell_', j, '.dat'
    open(unit=100, file=filename, status='replace', action='write')
    write(100,*) '# t r HI HII HeIS HeIM HeII e'
    write(100,'(ES12.4,ES12.4,6ES15.6)') t, r(i), &
        y0_HI(i), y0_HII(i), y0_HeIS(i), &
        y0_HeIM(i), y0_HeII(i), y0_e(i)
    close(100)

    write(filename,'(A,I0,A)') './output/phis/cell_', j, '.dat'
    open(unit=200, file=filename, status='replace', action='write')
    write(200,*) '# t r phiHI phiHeIS phiHeIM tauHI tauHeIS tauHeIM'
    write(200,'(ES12.4,ES12.4,6ES15.6)') t, r(i), &
        phi_all(j,1), phi_all(j,2), phi_all(j,3), &
        tau_all(j,1), tau_all(j,2), tau_all(j,3)
    close(200)
  end do

! =======================================================================
! LOOP PRINCIPAL EN EL TIEMPO
! =======================================================================
  do while (t < time)

    tf = dtime

    ! -------------------------------------------------------------------
    ! PASO 1: Barrida serial de taus (afuera -> adentro)
    ! No se puede evitar: tau(j) = tau(j-1) + dtau(j)
    ! Es rapida comparada con la quimica, no es el cuello de botella.
    ! -------------------------------------------------------------------
    tau_global = (/ 0.d0, 0.d0, 0.d0 /)

    do j = 1, n_points
      i = n_points - j + 1

      y0 = (/ y0_HI(i), y0_HeIS(i), y0_HeIM(i) /)
      call local_taus(y0, dr, tau_global, tau_loc)
      tau_global = tau_loc

      tau_all(j,:) = tau_loc
      call phis_rate(FluxHI, FluxHeIS, FluxHeIM, tau_loc, &
                     phi_all(j,1), phi_all(j,2), phi_all(j,3))
    end do

    ! -------------------------------------------------------------------
    ! PASO 2: Integracion quimica — PARALELIZADO con OpenMP
    !
    ! Cada celda es independiente dado que ya se conocen phi y tau.
    ! Todas las variables del loop son PRIVATE por thread.
    ! Los arrays compartidos (y0_*, phi_all, etc.) no tienen conflictos
    ! porque cada j accede a un indice diferente.
    ! -------------------------------------------------------------------
    !$OMP PARALLEL DO                                              &
    !$OMP   DEFAULT(NONE)                                         &
    !$OMP   SHARED(n_points, tf, T_profile, r,               &
    !$OMP          y0_HI, y0_HII, y0_HeIS, y0_HeIM,              &
    !$OMP          y0_HeII, y0_e, phi_all, y_local)               &
    !$OMP   PRIVATE(j, i, y, T_loc,                              &
    !$OMP           epsmn, epsmx, dtmn, tnot, itermx, ymn, prt)
    do j = 1, n_points
      i = n_points - j + 1

      ! Parametros del solver (privados por thread)
      epsmn  = 1e-5
      epsmx  = 0.0d0
      dtmn   = 0.0d0
      tnot   = 0.0d0
      itermx = 5
      ymn(:) = 1e-20
      prt    = 1.0d0

      call chemsp(epsmn, epsmx, dtmn, tnot, itermx, ns, ymn, prt)

      T_loc = T_profile(i)

      y(1) = y0_HI(i)
      y(2) = y0_HII(i)
      y(3) = y0_HeIS(i)
      y(4) = y0_HeIM(i)
      y(5) = y0_HeII(i)
      y(6) = y0_e(i)

      call chemeq2solve(tf, y, ns, T_loc, &
                        phi_all(j,1), phi_all(j,2), phi_all(j,3))

      ! Guardar resultado en buffer (cada j escribe en su propio i, sin conflicto)
      y_local(j,:) = y(:)

    end do
    !$OMP END PARALLEL DO

    ! -------------------------------------------------------------------
    ! PASO 3: Actualizar condiciones y escribir resultados (serial)
    !
    ! Se hace fuera del paralelo para evitar race conditions en I/O
    ! y para mantener la actualizacion de y0_* ordenada.
    ! -------------------------------------------------------------------
    do j = 1, n_points
      i = n_points - j + 1

      ! Actualizar condiciones iniciales para el proximo paso de tiempo
      y0_HI(i)   = y_local(j,1)
      y0_HII(i)  = y_local(j,2)
      y0_HeIS(i) = y_local(j,3)
      y0_HeIM(i) = y_local(j,4)
      y0_HeII(i) = y_local(j,5)
      y0_e(i)    = y_local(j,6)

      ! Escritura species
      write(filename,'(A,I0,A)') './output/species/cell_', j, '.dat'
      open(unit=100, file=filename, status='old', position='append', action='write')
      write(100,'(ES12.4,ES12.4,6ES15.6)') t+tf, r(i), y_local(j,:)
      close(100)

      ! Escritura phis
      write(filename,'(A,I0,A)') './output/phis/cell_', j, '.dat'
      open(unit=200, file=filename, status='old', position='append', action='write')
      write(200,'(ES12.4,ES12.4,6ES15.6)') t+tf, r(i), &
          phi_all(j,1), phi_all(j,2), phi_all(j,3), &
          tau_all(j,1), tau_all(j,2), tau_all(j,3)
      close(200)

    end do

    t = t + tf
    print *, 't = ', t, ' s of ', time, ' s'
    call flush(6)

  end do

  call cpu_time(t2)
  print *, 'Finished test_1D, Time = ', t2-t1, ' s'

contains

! =======================================================================
  subroutine phis_rate(FluxHI, FluxHeIS, FluxHeIM, taus, phiHI, phiHeIS, phiHeIM)
! =======================================================================
  implicit none
  real(8), intent(in)  :: FluxHI, FluxHeIS, FluxHeIM, taus(3)
  real(8), intent(out) :: phiHI, phiHeIS, phiHeIM
  real(8) :: a0_H, a0_HeIS, a0_HeIM

  a0_H    = 3.87d-18
  a0_HeIS = 4.03d-18
  a0_HeIM = 3.59d-18

  phiHI   = a0_H   * FluxHI   * exp(-taus(1))
  phiHeIS = a0_HeIS * FluxHeIS * exp(-taus(2))
  phiHeIM = a0_HeIM * FluxHeIM * exp(-taus(3))

  end subroutine phis_rate


! =======================================================================
  subroutine local_taus(y0, dr, taus_global, taus_loc)
! =======================================================================
  implicit none
  real(8), intent(in)  :: y0(3)
  real(8), intent(in)  :: dr
  real(8), intent(in)  :: taus_global(3)
  real(8), intent(out) :: taus_loc(3)
  real(8) :: a0_H, a0_HeIS, a0_HeIM

  a0_H    = 3.87d-18
  a0_HeIS = 4.03d-18
  a0_HeIM = 3.59d-18

  taus_loc(1) = a0_H   * y0(1) * dr + taus_global(1)
  taus_loc(2) = a0_HeIS * y0(2) * dr + taus_global(2)
  taus_loc(3) = a0_HeIM * y0(3) * dr + taus_global(3)

  end subroutine local_taus


end program test_1D
