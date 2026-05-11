program test_1D
  
  use solver
  implicit none

  character(len=50)  :: label
  character(len=200) :: perfil_densidad

  real(8) :: T, Rp, time, dtime
  real(8) :: FluxHI, FluxHeIS, FluxHeIM
  real(8) :: He_H_ratio
  real(8) :: h_ion_frac, he_ion_frac, he_metastable_frac

  integer, parameter :: ns = 6
  integer :: i, j, ios, n_points, total_lines, itermx

  real(kind=8)    :: ymn(ns), ti, tf, y(ns) !, yi(ns)
  real(kind=8)    :: epsmn, epsmx, dtmn, tnot, prt

  real(kind=8), allocatable :: r(:), rho(:)
  real(kind=8), allocatable :: y0_HI(:), y0_HII(:), y0_HeIS(:)
  real(kind=8), allocatable :: y0_HeIM(:), y0_HeII(:), y0_e(:)
  real(kind=8), allocatable :: yi_all(:,:)

  integer, allocatable :: unit_out(:), unit_phys(:)
  character(len=100) :: filename

  real(kind=8) :: mp, R_jupiter, r_planet, dr
  real(kind=8) :: T_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc
  real(kind=8) :: tau_global(3), tau_loc(3), y0(3)
  real(kind=8) :: t1, t2

  character(len=200) :: line

! LEER INPUTS
  open(99,file='./inputs/inputs.dat',status='old',action='read')

  read(99,*) label, T
  read(99,*) label, Rp
  read(99,*) label, FluxHI
  read(99,*) label, FluxHeIS
  read(99,*) label, FluxHeIM
  read(99,*) label, He_H_ratio
  read(99,*) label, h_ion_frac
  read(99,*) label, he_ion_frac
  read(99,*) label, he_metastable_frac
  read(99,*) label, dtime
  read(99,*) label, time
  read(99,*) label, perfil_densidad

  close(99)

  print *, 'Inputs:'
  print *, 'T = ', T, ' K'
  print *, 'Rp = ', Rp, ' R_jupiter'
  print *, 'FluxHI = ', FluxHI, ' cm^-2 s^-1' 
  print *, 'FluxHeIS = ', FluxHeIS, ' cm^-2 s^-1'
  print *, 'FluxHeIM = ', FluxHeIM, ' cm^-2 s^-1'
  print *, 'He/H ratio = ', He_H_ratio
  print *, 'H ionization fraction = ', h_ion_frac
  print *, 'He ionization fraction = ', he_ion_frac 
  print *, 'He metastable fraction = ', he_metastable_frac  
  print *, 'dtime = ', dtime, ' s'
  print *, 'time = ', time, ' s'

  call cpu_time(t1)

! CONTAR LINEAS PERFIL DENSIDAD
  open(10,file=trim(perfil_densidad),status='old',action='read')

  total_lines = 0
  do
    read(10,'(A)',iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) /= '#') total_lines = total_lines + 1
  end do

  close(10)

  n_points = total_lines

  allocate(r(n_points), rho(n_points))
  allocate(y0_HI(n_points), y0_HII(n_points))
  allocate(y0_HeIS(n_points), y0_HeIM(n_points))
  allocate(y0_HeII(n_points), y0_e(n_points))
  allocate(yi_all(n_points, ns))
  allocate(unit_out(n_points), unit_phys(n_points))

  ! LEER PERFIL DE DENSIDAD
  open(10,file=trim(perfil_densidad),status='old',action='read')

  j = 0
  do
    read(10,'(A)',iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) == '#') cycle

    j = j + 1
    read(line,*) r(j), rho(j)
  end do

  close(10)

  mp        = 1.67d-24
  R_jupiter = 7.1492d9
  r_planet  = 1.38 * R_jupiter
  dr        = r_planet*(r(n_points)-r(1))/real(n_points,8)
  T_loc = T

! PERFILES INICIALES
  open(30,file='./output/perfiles_iniciales.txt',status='replace')

  write(30,*) '# j y0_HI y0_HII y0_HeIS y0_HeIM y0_HeII y0_e'

  do j=1,n_points

    y0_HI(j)   = rho(j)*(1d0-He_H_ratio)*(1d0-h_ion_frac)/mp
    y0_HII(j)  = rho(j)*(1d0-He_H_ratio)*h_ion_frac/mp

    y0_HeIS(j) = rho(j)*He_H_ratio*(1d0-he_ion_frac)*(1d0-he_metastable_frac)/mp
    y0_HeIM(j) = rho(j)*He_H_ratio*(1d0-he_ion_frac)*he_metastable_frac/mp
    y0_HeII(j) = rho(j)*He_H_ratio*he_ion_frac/mp

    y0_e(j)    = y0_HII(j) + y0_HeII(j)

    write(30,'(I6,1X,6ES16.8)') j, y0_HI(j), y0_HII(j), y0_HeIS(j), &
                                y0_HeIM(j), y0_HeII(j), y0_e(j)
  end do

  close(30)

  ! LOOP PRINCIPAL
  t = 0.0

  ! ESCRIBIR CONDICIÓN INICIAL (t = 0)
  tau_global = (/ 0.d0, 0.d0, 0.d0 /)

  do j = 1, n_points

    i = n_points - j + 1

    ! Tau inicial
    y0 = (/ y0_HI(i), y0_HeIS(i), y0_HeIM(i) /)
    call local_taus(y0, dr, tau_global, tau_loc)

    tau_global = tau_loc

    ! Phi inicial
    call phis_rate(FluxHI, FluxHeIS, FluxHeIM, tau_loc, &
                 phiHI_loc, phiHeIS_loc, phiHeIM_loc)

    ! ===== species =====
    write(filename,'(A,I0,A)') './output/species/cell_', j, '.dat'
    open(unit=100, file=filename, status='replace', action='write')

    write(100,*) '# t r HI HII HeIS HeIM HeII e'
    write(100,'(ES12.4, ES12.4,6ES15.6)') t, r(i), &
        y0_HI(i), y0_HII(i), y0_HeIS(i), &
        y0_HeIM(i), y0_HeII(i), y0_e(i)

    close(100)

    ! ===== phis =====
    write(filename,'(A,I0,A)') './output/phis/cell_', j, '.dat'
    open(unit=200, file=filename, status='replace', action='write')

    write(200,*) '# t r phiHI phiHeIS phiHeIM tauHI tauHeIS tauHeIM'
    write(200,'(ES12.4, ES12.4,6ES15.6)') t, r(i), &
        phiHI_loc, phiHeIS_loc, phiHeIM_loc, &
        tau_global(1), tau_global(2), tau_global(3)

    close(200)

  end do

  do while (t < time)

    tf = dtime
    tau_global = (/ 0.d0, 0.d0, 0.d0 /)

    do j = 1, n_points

      ti = 0.0
      epsmn = 1e-5
      epsmx = 0.0
      dtmn  = 0.0
      tnot  = ti
      itermx = 5
      ymn(:) = 1e-20
      prt = 1.0

      call chemsp(epsmn, epsmx, dtmn, tnot, itermx, ns, ymn, prt)

      ! Índice invertido (radiación desde afuera)
      i = n_points - j + 1

      !tau local
      y0 = (/ y0_HI(i), y0_HeIS(i), y0_HeIM(i) /)
      call local_taus(y0, dr, tau_global, tau_loc)
      tau_global = tau_loc

      ! Phi  local
      call phis_rate(FluxHI, FluxHeIS, FluxHeIM, tau_loc, &
                    phiHI_loc, phiHeIS_loc, phiHeIM_loc)

      ! Inicializar estado químico
      y(1) = y0_HI(i)
      y(2) = y0_HII(i)
      y(3) = y0_HeIS(i)
      y(4) = y0_HeIM(i)
      y(5) = y0_HeII(i)
      y(6) = y0_e(i)

      ! Integración química
      call chemeq2solve(tf, y, ns, T_loc, &
                      phiHI_loc, phiHeIS_loc, phiHeIM_loc, j)

      ! ESCRITURA RESULTADOS 
      ! species 
      write(filename,'(A,I0,A)') './output/species/cell_', j, '.dat'
      open(unit=100, file=filename, status='old', position='append', action='write')
      write(100,'(ES12.4, ES12.4,6ES15.6)') t+tf, r(i), y(:)
      close(100)

      ! phis 
      write(filename,'(A,I0,A)') './output/phis/cell_', j, '.dat'
      open(unit=200, file=filename, status='old', position='append', action='write')
      write(200,'(ES12.4, ES12.4,6ES15.6)') t+tf, r(i), &
                phiHI_loc, phiHeIS_loc, phiHeIM_loc, &
                tau_global(1), tau_global(2), tau_global(3)
      close(200)

      ! Actualizar condiciones iniciales para el siguiente paso de tiempo
      y0_HI(i)   = y(1)
      y0_HII(i)  = y(2)
      y0_HeIS(i) = y(3)
      y0_HeIM(i) = y(4)
      y0_HeII(i) = y(5)
      y0_e(i)    = y(6)

    end do

    t = t + tf
    print *, 't = ', t, ' s of ', time, ' s'

  end do

  call cpu_time(t2)
  print *, 'Finished test_1D_1, Time = ', t2-t1, ' s'

contains

subroutine phis_rate(FluxHI, FluxHeIS, FluxHeIM, taus, phiHI, phiHeIS, phiHeIM)
  !esta subrutina calcula el valor de phi necesarios para el chemeq solver en funcion del tau local 
  !calculado a partir de las abundancias iniciales

  implicit none
  real(8), intent(in)  :: FluxHI, FluxHeIS, FluxHeIM, taus(3)
  real(8), intent(out) :: phiHI, phiHeIS, phiHeIM 
  real(8) :: a0_H, a0_HeIS, a0_HeIM
  
  a0_H    = 3.87d-18 !seccion eficaz HI
  a0_HeIS = 4.03d-18
  a0_HeIM = 3.59d-18

  
  phiHI   = a0_H*FluxHI*exp(-taus(1))
  phiHeIS = a0_HeIS*FluxHeIS*exp(-taus(2))
  phiHeIM = a0_HeIM*FluxHeIM*exp(-taus(3))

end subroutine phis_rate


subroutine local_taus(y0, dr, taus_global, taus_loc)
  !Esta rutina calcula el tau de cada dr para las abundancias de ese dr
  !El tau local esta compuesto por la suma de los taus de los dr anteriores (taus global) + dtau local
  
  implicit none
  real(8), intent(in)  :: y0(3)
  real(8), intent(in)  :: dr !tamaño de cada celda, en principio constante
  real(8), intent(in)  :: taus_global(3)
  real(8), intent(out) :: taus_loc(3)
  real(8) :: a0_H, a0_HeIS, a0_HeIM

  a0_H    = 3.87d-18 !seccion eficaz HI
  a0_HeIS = 4.03d-18 !
  a0_HeIM = 3.59d-18 !

  taus_loc(1) = a0_H*y0(1)*dr    + taus_global(1)
  taus_loc(2) = a0_HeIS*y0(2)*dr + taus_global(2)
  taus_loc(3) = a0_HeIM*y0(3)*dr + taus_global(3)

end subroutine local_taus


end program test_1D
