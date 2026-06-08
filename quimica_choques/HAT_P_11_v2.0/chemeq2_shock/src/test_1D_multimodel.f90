program test_1D_multimodel

  use solver
  implicit none

  character(len=50)  :: label
  character(len=200) :: filename, line

  real(kind=8) :: Rp, time, dtime
  real(kind=8) :: FluxHI, FluxHeIS, FluxHeIM

  integer, parameter :: ns = 6
  integer :: i, j, ic, ios, n_points, total_lines, itermx, modelo_len

  real(kind=8) :: ymn(ns), tf, y(ns)
  real(kind=8) :: epsmn, epsmx, dtmn, tnot, prt

  real(kind=8), allocatable :: r(:), rho(:), T_profile(:)
  real(kind=8), allocatable :: y0_HI(:), y0_HII(:), y0_HeIS(:)
  real(kind=8), allocatable :: y0_HeIM(:), y0_HeII(:), y0_e(:)

  real(kind=8) :: mp, dr, t_run
  real(kind=8) :: T_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc
  real(kind=8) :: tau_global(3), tau_loc(3), y0(3)
  real(kind=8) :: t1, t2
  real(kind=8) :: nHt, nHet, alpha, ymol, fHII, fHeII, n_tot

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

  mp = 1.67d-24   ! g

  call cpu_time(t1)

! =======================================================================
! CONTAR LINEAS: se hace una sola vez con el modelo 1
! (todos los modelos tienen la misma grilla espacial)
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
! LOOP SOBRE MODELOS
! =======================================================================
  do i = 1, modelo_len

    print *, ''
    print *, 'Running model ', i, ' of ', modelo_len

    ! Alocar arrays (se desalocan al final de cada modelo)
    allocate(r(n_points), rho(n_points), T_profile(n_points))
    allocate(y0_HI(n_points),  y0_HII(n_points))
    allocate(y0_HeIS(n_points), y0_HeIM(n_points))
    allocate(y0_HeII(n_points), y0_e(n_points))

    ! -------------------------------------------------------------------
    ! Leer perfil de densidad y temperatura del modelo i
    ! -------------------------------------------------------------------
    write(filename,'(A,I0,A)') '../shock_1D/uniform_output/model_', i, '.dat'
    open(10, file=trim(filename), status='old', action='read')
    print *, '  Reading: ', trim(filename)

    j = 0
    do
      read(10,'(A)', iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#') cycle
      j = j + 1
      read(line,*) r(j), rho(j), T_profile(j)
    end do
    close(10)

    ! dr calculado después de leer r
    dr = (r(n_points) - r(1)) / real(n_points, 8)
    print *, '  dr = ', dr, ' cm'

    ! -------------------------------------------------------------------
    ! Condiciones iniciales de las especies quimicas
    ! -------------------------------------------------------------------
    do j = 1, n_points
      n_tot = rho(j) / mp

      nHt  = (1d0 - ymol) / (1d0 + 3d0*ymol)
      nHet = ymol          / (1d0 + 3d0*ymol)

      y0_HI(j)  = (1d0 - fHII)  * nHt  * n_tot
      y0_HII(j) = fHII           * nHt  * n_tot

      y0_HeIS(j) = (1d0 - fHeII) * nHet * n_tot / (1d0 + alpha)
      y0_HeIM(j) = alpha * y0_HeIS(j)
      y0_HeII(j) = fHeII * nHet  * n_tot

      y0_e(j) = y0_HII(j) + y0_HeII(j)
    end do

    ! -------------------------------------------------------------------
    ! Evolucion temporal de la quimica
    ! -------------------------------------------------------------------
    t_run      = 0.0d0
    tau_global = (/ 0.d0, 0.d0, 0.d0 /)

    do while (t_run < time)

      tf         = dtime
      tau_global = (/ 0.d0, 0.d0, 0.d0 /)

      ! Barrer de afuera hacia adentro (ic = indice de celda)
      do j = 1, n_points
        ic = n_points - j + 1

        epsmn  = 1e-5
        epsmx  = 0.0d0
        dtmn   = 0.0d0
        tnot   = 0.0d0
        itermx = 5
        ymn(:) = 1e-20
        prt    = 1.0d0

        call chemsp(epsmn, epsmx, dtmn, tnot, itermx, ns, ymn, prt)

        ! Tau acumulado hasta esta celda
        y0 = (/ y0_HI(ic), y0_HeIS(ic), y0_HeIM(ic) /)
        call local_taus(y0, dr, tau_global, tau_loc)
        tau_global = tau_loc

        ! Flujos ionizantes locales
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

        ! Actualizar para el siguiente paso temporal
        y0_HI(ic)   = y(1)
        y0_HII(ic)  = y(2)
        y0_HeIS(ic) = y(3)
        y0_HeIM(ic) = y(4)
        y0_HeII(ic) = y(5)
        y0_e(ic)    = y(6)

      end do

      t_run = t_run + tf

    end do

    ! -------------------------------------------------------------------
    ! Guardar perfil FINAL de especies, phis y taus para este modelo.
    ! Se hace una pasada final de taus con las abundancias evolucionadas.
    ! -------------------------------------------------------------------
    write(filename,'(A,I0,A)') './output/final_model_', i, '.dat'
    open(unit=50, file=trim(filename), status='replace', action='write')
    write(50,'(A)') &
      '# r[cm]  HI  HII  HeIS  HeIM  HeII  e  phiHI  phiHeIS  phiHeIM  tauHI  tauHeIS  tauHeIM'

    tau_global = (/ 0.d0, 0.d0, 0.d0 /)

    do j = 1, n_points
      ic = n_points - j + 1

      y0 = (/ y0_HI(ic), y0_HeIS(ic), y0_HeIM(ic) /)
      call local_taus(y0, dr, tau_global, tau_loc)
      tau_global = tau_loc

      call phis_rate(FluxHI, FluxHeIS, FluxHeIM, tau_loc, &
                     phiHI_loc, phiHeIS_loc, phiHeIM_loc)

      write(50,'(13ES15.6)') r(ic), &
          y0_HI(ic), y0_HII(ic), y0_HeIS(ic), &
          y0_HeIM(ic), y0_HeII(ic), y0_e(ic), &
          phiHI_loc, phiHeIS_loc, phiHeIM_loc, &
          tau_loc(1), tau_loc(2), tau_loc(3)
    end do

    close(50)
    print *, '  Saved: ', trim(filename)

    ! Desalocar para el siguiente modelo
    deallocate(r, rho, T_profile)
    deallocate(y0_HI, y0_HII, y0_HeIS, y0_HeIM, y0_HeII, y0_e)

  end do

  call cpu_time(t2)
  print *, ''
  print *, 'Finished. Total time = ', t2-t1, ' s'

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

    taus_loc(1) = a0_H    * y0(1) * dr + taus_global(1)
    taus_loc(2) = a0_HeIS * y0(2) * dr + taus_global(2)
    taus_loc(3) = a0_HeIM * y0(3) * dr + taus_global(3)

  end subroutine local_taus

end program test_1D_multimodel