program shock_grid 
  implicit none 
  real(kind=8), parameter :: chiH = 2.179d-11   ! erg,  energia de ionizacion del H (13.6eV)
  real(kind=8), parameter :: kB   = 1.38d-16    ! erg/K
  real(kind=8), parameter :: mp   = 1.67d-24    ! g
  real(kind=8), parameter :: mu   = 1.3d0       ! masa media (1 para H puro, 1.3 para 10% He)
  real(kind=8), parameter :: T_corte = 1.0d2 ! K, temperatura de corte para terminar el modelo
  integer :: i, j, z, k, n_model
  character(len=100) :: filename
  
  integer, parameter :: n0_points = 5
  integer, parameter :: T0_points = 5
  integer, parameter :: u0_points = 5
  integer, parameter :: y0_points = 1

  real(kind=8) :: n0_grid(n0_points), T0_grid(T0_points)
  real(kind=8) :: u0_grid(u0_points), y0_grid(y0_points)

!# r[Rp] density[g/cm^3] temperature[K] velocity[cm/s]
!2.5580000000e+00 1.2897914462e+08 3.8776549256e+03 3.9891866155e+05
!3.5580000000e+00 4.3128190849e+07 3.1180117755e+03 6.1651556661e+05
!4.5580000000e+00 2.1083157393e+07 2.6204164922e+03 7.6838252852e+05
!5.5580000000e+00 1.2285928831e+07 2.2681740462e+03 8.8669059850e+05
!6.5580000000e+00 7.9338270250e+06 2.0030205610e+03 9.8617052484e+05 


  n0_grid = (/ 1.2897914462d+08, 4.3128190849d+07, 2.1083157393d+07, 1.2285928831d+07, 7.9338270250d+06 /) ! cm^-3
  T0_grid = (/ 3.8776549256d+03, 3.1180117755d+03, 2.6204164922d+03, 2.2681740462d+03, 2.0030205610d+03 /) ! K
  u0_grid = (/ 3.9891866155d+05, 6.1651556661d+05, 7.6838252852d+05, 8.8669059850d+05, 9.8617052484d+05 /) ! cm/s
  y0_grid = (/ 0.95d0  /) ! fraccion de ionizacion inicial
  
  n_model = 0

  open(10, file='./output/model_list.dat', status='replace', action='write')
  write(10,'(A)') '# model_number  T0[K]  n0[cm^-3]  u0[cm/s]  y0'
  !Barrido de la grilla de condiciones iniciales
  do i = 1, T0_points
    !do j = 1, n0_points
      !do k = 1, u0_points
        !do z = 1, y0_points
    !numero de modelo para nombrar los archivos de salida
    n_model = n_model + 1
    !printeo en pantalla el modelo y guardo las condiciones iniciales del mismo
    print *, 'Running model ', n_model, 'of', T0_points!*n0_points*u0_points*y0_points
      
    write(10,*) n_model, T0_grid(i), n0_grid(i), u0_grid(i), y0_grid(1)
          
    !corro el modelo con las condiciones iniciales dadas y guardo los resultados
    call run_model(y0_grid(1), n0_grid(i), u0_grid(i), T0_grid(i), n_model)
 
  !      end do
  !    end do
  !  end do  
  end do

  close(10)

contains

!=======================================================================
subroutine coef_c(T, c) ! Coeficiente de ionizacion colisional del H
  implicit none
  real(kind=8), intent(in)  :: T
  real(kind=8), intent(out) :: c

  c = 5.83d-11 * (T**0.5d0) * exp(-157800.0d0/T)

end subroutine coef_c
!=======================================================================

!=======================================================================
subroutine coef_alpha(T, alpha) ! Coeficiente de recombinacion del H
  implicit none
  real(kind=8), intent(in)  :: T
  real(kind=8), intent(out) :: alpha

  alpha = 3.69d-10 * (T**(-0.79d0))

end subroutine coef_alpha
!=======================================================================

!=======================================================================
subroutine run_model(y0i, n0i, u0i, T0i, model_number)
  ! Aqui iria el codigo que corre el modelo con las condiciones dadas
  ! y devuelve los resultados. Por ahora solo es un placeholder.
  ! inputs: 
  ! - y0i: fraccion de ionizacion de H
  ! - n0i: densidad de particulas (cm^-3)
  ! - u0i: velocidad del fluido (cm/s)
  ! - T0i: temperatura (K)
  ! - model_number: numero del modelo a correr (integer)
  ! outputs: perfil de T y rho.

  implicit none

  real(kind=8), intent(in) :: y0i, n0i, u0i, T0i
  integer, intent(in) :: model_number
  
  real(kind=8) :: y, n, u, T
  real(kind=8) :: rho, P, h, x, dx, G, L
  real(kind=8) :: rho0, u0, P0, n0, ne
  real(kind=8) :: c, alpha, dhdx, dydx

  integer :: contador
  character(len=100) :: filename

  !guardo las condiciones iniciales en variables locales para no modificar los inputs
  y = y0i
  n = n0i
  u = u0i
  T = T0i

  ! Condiciones iniciales calciladas a partir de los inputs
  rho = n0i*mu*mp      ! g/cm^3
  P   = n0i*kB*T0i       ! dyn/cm^2
  h   = (u0i**2)/2.0d0 + (5.0d0/2.0d0)*(P/rho)  ! erg/g

  ! Otras condiciones iniciales y variables necesarias para el modelo
  x        = 0.0d0
  dx       = 1.0d3   ! cm, paso inicial
  G        = 0.0d0   ! calentamiento 

  ! Perdida de energia inicial
  call coef_c(T0i, c)
  L = n0i**2 * y0i*(1.0d0 - y0i)*c*chiH                             &
    + n0i**2 * y0i*(6.1d-19)*(T0i**(-0.63d0))                       &
            *(1.0d0 - exp(-(T0i/1.0d5)**1.63d0))

  ! Guardo las condiciones iniciales
  write(filename,'(A,I0,A)') './output/output_', model_number, '.dat'
  open(20, file=filename, status='replace', action='write')
  write(20,'(A)') '# x[cm]  rho[g/cm^3]   T[K]    y    n[cm^-3]    u[cm/s]    P[dyn/cm^2]    h[erg/g]    L[dyn/cm^3]'
  write(20,*) x, rho, T0i, y, n, u, P, h, L

  contador = 0
  !Si mi temperatura inicial es mayor a la de corte, 
  !entonces el contador arranca en 1 para que se corte el loop la primera vez que T < T_corte,
  !, sino en 0, entonces el loop se corta la segunda vez que T < T_corte
  if (T0i > T_corte) contador = contador + 1
  
  ! Loop principal: termina cuando T < T_corte dos veces consecutivas 
  ! El primero es para que entre al loop, el segundo para que termine
  do while (contador < 2)

    if (T <= T_corte) contador = contador + 1

    call coef_c(T, c)
    call coef_alpha(T, alpha)

    dhdx = (G - L) / (rho*u)
    dydx = (n*y/u) * ((1.0d0 - y)*c - y*alpha)

    ! Achico dx si los cambios relativos superan el 1%
    do while (abs(dhdx*dx/h) > 0.01d0 .or. abs(dydx*dx/y) > 0.01d0)
      dx = 0.5d0*dx
    end do

    ! Avanzar
    x = x + dx
    h = h + dhdx*dx
    y = y + dydx*dx

    ! Aumento dx si los cambios relativos son menores al 0.1%
    if (abs(dhdx*dx/h) < 0.001d0 .and. abs(dydx*dx/y) < 0.001d0) then
      dx = 2.0d0*dx
    end if

    ! Guardar estado del paso anterior 
    rho0 = rho
    u0   = u
    P0   = P
    n0   = n

    ! Actualizar variables termodinamicas 
    rho = ( (5.0d0/2.0d0)*(P0 + rho0*u0**2)                   &
          + sqrt( ((5.0d0/2.0d0)*(P0 + rho0*u0**2))**2        &
                  - 8.0d0*h*rho0**2*u0**2 ) )                  &
          / (2.0d0*h)

    P  = P0 + rho0*u0**2*(1.0d0 - rho0/rho)
    u  = rho0*u0/rho
    n  = rho/(mu*mp)
    ne = y*n
    T  = P/((n + ne)*kB)

    ! Perdida radiativa
    L = n**2 * y*(1.0d0 - y)*c*chiH                           &
      + n**2 * y*(6.1d-19)*(T**(-0.63d0))                     &
              *(1.0d0 - exp(-(T/1.0d5)**1.63d0))

    write(20,*) x, rho, T, y, n, u, P, h, L

  end do
  
  close(20)
  

end subroutine run_model



!=======================================================================

end program shock_grid