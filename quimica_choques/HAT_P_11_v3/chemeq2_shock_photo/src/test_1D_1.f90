program test_1D
  
  use solver
  implicit none

  character(len=50)  :: label
  !character(len=200) :: modelo
  integer :: modelo_len
  !character(len=200) :: perfil_temperatura 

  real(kind = 8) :: T, Rp, time, dtime
  real(kind = 8) :: FluxHI, FluxHeIS, FluxHeIM

  integer, parameter :: ns = 6
  integer :: i, j, ios, n_points, total_lines, itermx

  real(kind=8)    :: ymn(ns), ti, tf, y(ns) !, yi(ns)
  real(kind=8)    :: epsmn, epsmx, dtmn, tnot, prt

  real(kind=8), allocatable :: r(:), rho(:), T_profile(:)
  real(kind=8), allocatable :: y0_HI(:), y0_HII(:), y0_HeIS(:)
  real(kind=8), allocatable :: y0_HeIM(:), y0_HeII(:), y0_e(:)
  real(kind=8), allocatable :: yi_all(:,:)

  integer, allocatable :: unit_out(:), unit_phys(:)
  character(len=100) :: filename

  real(kind=8) :: mp, R_jupiter, r_planet, dr
  real(kind=8) :: T_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc
  real(kind=8) :: tau_global(3), tau_loc(3), y0(3)
  real(kind=8) :: t1, t2
  real(kind=8) :: r_dummy
  real(kind=8) :: nHt, nHet, alpha, ymol, fHII, fHeII, n_tot

  character(len=200) :: line

! LEER INPUTS
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
  read(99,*) label, modelo_len
  !read(99,*) label, perfil_temperatura

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
  print *, 'Cantidad de modelos = ', modelo_len
  
  call cpu_time(t1)

! CONTAR LINEAS PERFIL DENSIDAD, solo para el primer modelo, el resto es igual
  open(10,file='../shock_1D/uniform_output/model_1.dat',status='old',action='read')
  print *, 'Reading density profile from: ', trim(modelo)

  total_lines = 0
  do
    read(10,'(A)',iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) /= '#') total_lines = total_lines + 1
  end do

  close(10)

  n_points = total_lines

  mp        = 1.67d-24        ! protón mass en g
  !El r ya esta en cm
  dr = (r(n_points)-r(1))/real(n_points,8)
  
  print *, 'dr = ', dr, ' cm'
  print *, 'r_planet = ', r_planet, ' cm'
  call flush(6)

  
  do i = 1, modelo_len

    print *, 'Running model ', i, ' of ', modelo_len

    ! Perfiles finales para cada modelo
    write(filename,'(A,I0,A)') './output/finales', i, '.dat'
    open(unit=100, file=filename, status='replace', action='write')
    write(100,*) '# model x HI HII HeIS HeIM HeII e phiHI phiHeIS phiHeIM tauHI tauHeIS tauHeIM'  
    close(100)

    allocate(r(n_points), rho(n_points), T_profile(n_points))
    allocate(y0_HI(n_points), y0_HII(n_points))
    allocate(y0_HeIS(n_points), y0_HeIM(n_points))
    allocate(y0_HeII(n_points), y0_e(n_points))
    allocate(yi_all(n_points, ns))
    allocate(unit_out(n_points), unit_phys(n_points))

    ! LEER PERFIL DE DENSIDAD Y TEMPERATURA
    write(filename,'(A,I0,A)') '../shock_1D/uniform_output/model_', i, '.dat'
    open(10,file=filename,status='old',action='read')
    print *, 'Reading density and temperature profiles from: ', trim(filename)

    j = 0
    do
      read(10,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#') cycle

      j = j + 1
      read(line,*) r(j), rho(j), T_profile(j)
    end do

    close(10)


    do j=1,n_points
      ! Densidad total numérica
      n_tot = rho(j)/mp   ! o dens(j)/rhosc si rho ya viene escalada

      ! Fracciones H/He consistentes
      nHt  = (1d0 - ymol)/(1d0 + 3d0*ymol)
      nHet = ymol/(1d0 + 3d0*ymol)

      ! Hidrógeno
      y0_HI(j)  = (1d0 - fHII) * nHt * n_tot
      y0_HII(j) = fHII * nHt * n_tot

      ! Helio
      y0_HeIS(j) = (1d0 - fHeII) * nHet * n_tot / (1d0 + alpha)
      y0_HeIM(j) = alpha * y0_HeIS(j)
      y0_HeII(j) = fHeII * nHet * n_tot

      ! Electrones 
      y0_e(j) = y0_HII(j) + y0_HeII(j)

    end do
    
    ! LOOP PRINCIPAL
    t = 0.0

    ! ESCRIBIR CONDICIÓN INICIAL (t = 0)
    tau_global = (/ 0.d0, 0.d0, 0.d0 /)

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

        ! temperatura local
        T_loc = T_profile(i)

        ! Inicializar estado químico
        y(1) = y0_HI(i)
        y(2) = y0_HII(i)
        y(3) = y0_HeIS(i)
        y(4) = y0_HeIM(i)
        y(5) = y0_HeII(i)
        y(6) = y0_e(i)

        ! Integración química
        call chemeq2solve(tf, y, ns, T_loc, &
                        phiHI_loc, phiHeIS_loc, phiHeIM_loc)

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
      !print *, 't = ', t, ' s of ', time, ' s'
      call flush(6)

    end do

    open(unit=100, file=filename, status='old', position='append', action='write')
    write(100,'(ES12.4, ES12.4,6ES15.6)') r(i), y(:), &
                  phiHI_loc, phiHeIS_loc, phiHeIM_loc, &
                  tau_global(1), tau_global(2), tau_global(3)
    close(100)

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
