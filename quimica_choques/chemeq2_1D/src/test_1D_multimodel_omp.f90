program test_1D_multimodel

  use omp_lib
  use solver
  implicit none

  character(len=50)  :: label
  character(len=200) :: line

  real(kind=8) :: Rp, time, dtime
  real(kind=8) :: FluxHI, FluxHeIS, FluxHeIM
  real(kind=8) :: mp, t1, t2
  real(kind=8) :: ymol, fHII, fHeII, alpha

  integer, parameter :: ns = 6
  integer :: i, ios, n_points, total_lines, modelo_len

! =======================================================================
! LEER INPUTS
! =======================================================================
  open(99, file='./inputs/inputs.dat', status='old', action='read')

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
  read(99,*) label, modelo_len

  close(99)

  print *, 'Inputs:'
  print *, 'FluxHI   = ', FluxHI,   ' cm^-2 s^-1'
  print *, 'FluxHeIS = ', FluxHeIS, ' cm^-2 s^-1'
  print *, 'FluxHeIM = ', FluxHeIM, ' cm^-2 s^-1'
  print *, 'He/H ratio          = ', ymol
  print *, 'H  ionization frac  = ', fHII
  print *, 'He ionization frac  = ', fHeII
  print *, 'He metastable frac  = ', alpha
  print *, 'dtime = ', dtime, ' s'
  print *, 'time  = ', time,  ' s'
  print *, 'Cantidad de modelos = ', modelo_len
  print *, 'Threads disponibles = ', omp_get_max_threads()

  mp = 1.67d-24

  call cpu_time(t1)

! =======================================================================
! CONTAR LINEAS: una sola vez con el modelo 1
! =======================================================================
  open(10, file='../shock_1D/uniform_output/model_1.dat', &
       status='old', action='read')

  total_lines = 0
  do
    read(10,'(A)', iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) /= '#') total_lines = total_lines + 1
  end do
  close(10)

  n_points = total_lines
  print *, 'n_points = ', n_points

! =======================================================================
! LOOP SOBRE MODELOS — PARALELIZADO
!
! Cada modelo es independiente: perfil propio, arrays propios,
! archivo de salida propio. SCHEDULE(DYNAMIC) es importante si los
! modelos tienen distinto costo computacional (distinto tiempo de
! convergencia de la quimica), para balancear la carga entre threads.
! =======================================================================
  !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
  do i = 1, modelo_len
    call run_model(i, n_points, ns, dtime, time, mp, &
                   FluxHI, FluxHeIS, FluxHeIM,       &
                   ymol, fHII, fHeII, alpha)
  end do
  !$OMP END PARALLEL DO

  call cpu_time(t2)
  print *, ''
  print *, 'Finished. Total time = ', t2-t1, ' s'

contains

! =======================================================================
  subroutine run_model(model_id, n_points, ns, dtime, time, mp, &
                       FluxHI, FluxHeIS, FluxHeIM,              &
                       ymol, fHII, fHeII, alpha)
!
! Corre un modelo completo de forma autocontenida.
! Al ser una subrutina, todas sus variables son automaticamente
! locales (stack) → thread-safe sin necesidad de clausulas PRIVATE.
! =======================================================================
    use omp_lib
    use solver
    implicit none

    integer,      intent(in) :: model_id, n_points, ns
    real(kind=8), intent(in) :: dtime, time, mp
    real(kind=8), intent(in) :: FluxHI, FluxHeIS, FluxHeIM
    real(kind=8), intent(in) :: ymol, fHII, fHeII, alpha

    ! Variables locales — cada thread tiene su propia copia
    integer :: j, ic, itermx, ios
    real(kind=8) :: tf, dr, t_run, T_loc, n_tot
    real(kind=8) :: epsmn, epsmx, dtmn, tnot, prt
    real(kind=8) :: ymn(ns), y(ns), y0(3)
    real(kind=8) :: tau_global(3), tau_loc(3)
    real(kind=8) :: phiHI_loc, phiHeIS_loc, phiHeIM_loc
    real(kind=8) :: nHt, nHet
    character(len=200) :: filename, line

    real(kind=8), allocatable :: r(:), rho(:), T_profile(:)
    real(kind=8), allocatable :: y0_HI(:), y0_HII(:), y0_HeIS(:)
    real(kind=8), allocatable :: y0_HeIM(:), y0_HeII(:), y0_e(:)

    ! Unit de I/O unico por modelo para evitar conflictos entre threads
    integer :: unit_in, unit_out
    unit_in  = 1000 + model_id
    unit_out = 2000 + model_id

    print '(A,I4,A,I4,A,I2)', &
        ' Starting model ', model_id, ' of ', modelo_len, &  ! modelo_len es shared
        '  (thread ', omp_get_thread_num(), ')'

    allocate(r(n_points), rho(n_points), T_profile(n_points))
    allocate(y0_HI(n_points),  y0_HII(n_points))
    allocate(y0_HeIS(n_points), y0_HeIM(n_points))
    allocate(y0_HeII(n_points), y0_e(n_points))

    ! -----------------------------------------------------------------
    ! Leer perfil
    ! -----------------------------------------------------------------
    write(filename,'(A,I0,A)') &
        '../shock_1D/uniform_output/model_', model_id, '.dat'
    open(unit=unit_in, file=trim(filename), status='old', action='read')

    j = 0
    do
      read(unit_in,'(A)', iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#') cycle
      j = j + 1
      read(line,*) r(j), rho(j), T_profile(j)
    end do
    close(unit_in)

    dr = (r(n_points) - r(1)) / real(n_points, 8)

    ! -----------------------------------------------------------------
    ! Condiciones iniciales
    ! -----------------------------------------------------------------
    nHt  = (1d0 - ymol) / (1d0 + 3d0*ymol)
    nHet = ymol          / (1d0 + 3d0*ymol)

    do j = 1, n_points
      n_tot = rho(j) / mp

      y0_HI(j)  = (1d0 - fHII)  * nHt  * n_tot
      y0_HII(j) = fHII           * nHt  * n_tot

      y0_HeIS(j) = (1d0 - fHeII) * nHet * n_tot / (1d0 + alpha)
      y0_HeIM(j) = alpha * y0_HeIS(j)
      y0_HeII(j) = fHeII * nHet  * n_tot

      y0_e(j) = y0_HII(j) + y0_HeII(j)
    end do

    ! -----------------------------------------------------------------
    ! Evolucion temporal
    ! -----------------------------------------------------------------
    t_run = 0.0d0

    do while (t_run < time)

      tf         = dtime
      tau_global = (/ 0.d0, 0.d0, 0.d0 /)

      do j = 1, n_points
        !ic = n_points - j + 1 ! Barrer de afuera hacia adentro
        ic = j                   ! Barrer de adentro hacia afuera

        epsmn  = 1e-5
        epsmx  = 0.0d0
        dtmn   = 0.0d0
        tnot   = 0.0d0
        itermx = 5
        ymn(:) = 1e-20
        prt    = 1.0d0

        call chemsp(epsmn, epsmx, dtmn, tnot, itermx, ns, ymn, prt)

        y0 = (/ y0_HI(ic), y0_HeIS(ic), y0_HeIM(ic) /)
        call local_taus(y0, dr, tau_global, tau_loc)
        tau_global = tau_loc

        call phis_rate(FluxHI, FluxHeIS, FluxHeIM, tau_loc, &
                       phiHI_loc, phiHeIS_loc, phiHeIM_loc)

        T_loc = T_profile(ic)

        y(1) = y0_HI(ic)
        y(2) = y0_HII(ic)
        y(3) = y0_HeIS(ic)
        y(4) = y0_HeIM(ic)
        y(5) = y0_HeII(ic)
        y(6) = y0_e(ic)

        call chemeq2solve(tf, y, ns, T_loc, &
                          phiHI_loc, phiHeIS_loc, phiHeIM_loc)

        y0_HI(ic)   = y(1)
        y0_HII(ic)  = y(2)
        y0_HeIS(ic) = y(3)
        y0_HeIM(ic) = y(4)
        y0_HeII(ic) = y(5)
        y0_e(ic)    = y(6)

      end do

      t_run = t_run + tf

    end do

    ! -----------------------------------------------------------------
    ! Guardar perfil final
    ! Cada thread escribe en su propio archivo: sin race conditions
    ! -----------------------------------------------------------------
    write(filename,'(A,I0,A)') './output/final_model_', model_id, '.dat'
    open(unit=unit_out, file=trim(filename), status='replace', action='write')
    write(unit_out,'(A)') &
      '# r[cm]  HI  HII  HeIS  HeIM  HeII  e  phiHI  phiHeIS  phiHeIM  tauHI  tauHeIS  tauHeIM'

    tau_global = (/ 0.d0, 0.d0, 0.d0 /)

    do j = 1, n_points
      !ic = n_points - j + 1 ! Barrer de afuera hacia adentro
      ic = j                 ! Barrer de adentro hacia afuera

      y0 = (/ y0_HI(ic), y0_HeIS(ic), y0_HeIM(ic) /)
      call local_taus(y0, dr, tau_global, tau_loc)
      tau_global = tau_loc

      call phis_rate(FluxHI, FluxHeIS, FluxHeIM, tau_loc, &
                     phiHI_loc, phiHeIS_loc, phiHeIM_loc)

      write(unit_out,'(13ES15.6)') r(ic), &
          y0_HI(ic), y0_HII(ic), y0_HeIS(ic), &
          y0_HeIM(ic), y0_HeII(ic), y0_e(ic), &
          phiHI_loc, phiHeIS_loc, phiHeIM_loc, &
          tau_loc(1), tau_loc(2), tau_loc(3)
    end do

    close(unit_out)

    deallocate(r, rho, T_profile)
    deallocate(y0_HI, y0_HII, y0_HeIS, y0_HeIM, y0_HeII, y0_e)

    print '(A,I4,A,I2)', &
        ' Done model ', model_id, '  (thread ', omp_get_thread_num(), ')'

  end subroutine run_model


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

    phiHI   = a0_H    * FluxHI   * exp(-taus(1))
    phiHeIS = a0_HeIS * FluxHeIS * exp(-taus(2))
    phiHeIM = a0_HeIM * FluxHeIM * exp(-taus(3))

  end subroutine phis_rate


! =======================================================================
  subroutine local_taus(y0, dr, taus_global, taus_loc)
! =======================================================================
    implicit none
    real(8), intent(in)  :: y0(3), dr, taus_global(3)
    real(8), intent(out) :: taus_loc(3)
    real(8) :: a0_H, a0_HeIS, a0_HeIM

    a0_H    = 3.87d-18
    a0_HeIS = 4.03d-18
    a0_HeIM = 3.59d-18

    taus_loc(1) = min(a0_H    * y0(1) * dr + taus_global(1), 140.d0)
    taus_loc(2) = min(a0_HeIS * y0(2) * dr + taus_global(2), 140.d0)
    taus_loc(3) = min(a0_HeIM * y0(3) * dr + taus_global(3), 140.d0)

  end subroutine local_taus

end program test_1D_multimodel
