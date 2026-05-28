program shock

  implicit none

  real(kind=8) :: chiH, kB, mp, mu, x, dx
  real(kind=8) :: h, y, rho, u, T, P, n, ne
  real(kind=8) :: rho0, u0, P0, n0
  real(kind=8) :: G, L, dhdx, dydx, c, alpha
  integer      :: contador

  ! Constantes fisicas
  chiH = 2.179d-11   ! erg,  energia de ionizacion del H
  kB   = 1.38d-16    ! erg/K
  mp   = 1.67d-24    ! g
  mu   = 1.0d0       ! masa media (1 para H puro, 1.3 para 10% He)

  ! Condiciones iniciales 
  u   = 50.0d0*1.0d5 ! cm/s
  n   = 100.0d0      ! cm^-3
  rho = n*mu*mp      ! g/cm^3
  y   = 1.0d-4       ! fraccion de ionizacion inicial
  T   = 100.0d0      ! K
  P   = n*kB*T       ! dyn/cm^2
  h   = (u**2)/2.0d0 + (5.0d0/2.0d0)*(P/rho)  ! erg/g

  x        = 0.0d0
  dx       = 1.0d9   ! cm, paso inicial
  G        = 0.0d0   ! calentamiento 
  contador = 0

  ! Perdida de energia inicial
  call coef_c(T, c)
  L = n**2 * y*(1.0d0 - y)*c*chiH                             &
    + n**2 * y*(6.1d-19)*(T**(-0.63d0))                       &
            *(1.0d0 - exp(-(T/1.0d5)**1.63d0))

  ! Guardo las condiciones iniciales
  open(20, file='output.dat', status='unknown', action='write')
  write(20,'(A)') '# x[cm]  h[erg/g]  y  rho[g/cm^3]  u[cm/s]  T[K]  P[dyn/cm^2]  n[cm^-3]'
  write(20,*) x, h, y, rho, u, T, P, n

  ! Loop principal: termina cuando T < 1e3 K dos veces consecutivas 
  ! El primero es para que entre al loop, el segundo para que termine
  do while (contador < 2)

    if (T <= 1.0d3) contador = contador + 1

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

    write(20,*) x, h, y, rho, u, T, P, n

  end do

  close(20)

contains

! =======================================================================
  subroutine coef_c(T, c) ! Coeficiente de ionizacion colisional del H
    implicit none
    real(kind=8), intent(in)  :: T
    real(kind=8), intent(out) :: c

    c = 5.83d-11 * (T**0.5d0) * exp(-157800.0d0/T)

  end subroutine coef_c
! =======================================================================

! =======================================================================
  subroutine coef_alpha(T, alpha) ! Coeficiente de recombinacion del H
    implicit none
    real(kind=8), intent(in)  :: T
    real(kind=8), intent(out) :: alpha

    alpha = 3.69d-10 * (T**(-0.79d0))

  end subroutine coef_alpha
! =======================================================================

end program shock