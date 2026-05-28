program shock 

real(kind=8) :: G, L, rho, u, T, y0i, dhdx, dydx
real(kind =8) :: u0i, rho0i, mu, mp, n0i, nHIIi, T0i, chiH, kB, P0i, h0i

chiH = 2.179d-11 ! erg, energia de ionizacion de H

u0i = 50*1d5 ! cm/s
!rho0 = 1e-12 ! g/cm^3
mu = 1.0 ! masa media, 1 para H, 1.3 para 10% de He
mp = 1.67d-24 ! g

n0i = 100 ! cm^-3
rho0i = n0i*mu*mp ! g/cm^3
!nHII0 = 0.076*n0i ! ionizacion inicial
!y0 = nHII0/n0i
y0i = 1.0d-4 ! ionizacion inicial
nHII0i = y0i*n0i ! ionizacion inicial
T0i = 100.0 ! K


kB = 1.38d-16 ! erg/K
P0i = n0i*kB*T0i
h0i = (u0i**2)/2 + (5/2)*(P0i/rho0i) !entalpia inicial

call evolve(h0i, y0i, rho0i, u0i, T0i, P0i, n0i)

contains


subroutine coef_c(T, c)
  real(kind=8), intent(in) :: T
  real(kind=8), intent(out) :: c
  c = 5.83d-11*(T**0.5) * exp(-157800/T)
end subroutine coef_c

subroutine coef_alpha(T, alpha)
  real(kind=8), intent(in) :: T
  real(kind=8), intent(out) :: alpha
  alpha = 3.69d-10*(T**(-0.79d0))
end subroutine coef_alpha

subroutine evolve(h0i, y0i, rho0i, u0i, T0i, P0i, n0i)
  real(kind=8) :: x, dx, h, y, P, n, ne
  real(kind=8), intent(in) :: h0i, y0i, rho0i, u0i, T0i, P0i, n0i
  real(kind=8) :: h0, y0, rho0, u0, T0, P0
  real(kind=8) :: G, L, rho, u, T, ny, dhdx, dydx, L0
  real(kind=8) :: c, alpha
  integer :: contador 
  contador = 0
  x = 0.0
  dx = 1.0d9 ! cm

  !A Fortran no le gusta pisar inputs
  h0 = h0i
  y0 = y0i
  rho0 = rho0i
  u0 = u0i
  T0 = T0i
  P0 = P0i
  n0 = n0i


  !Definimos las primeras variables para el primer paso.

  h=h0
  y=y0
  rho = rho0
  u = u0
  T = T0
  P = P0
  n = n0
  ! Calculamos la pérdida de energía inicial
  call coef_c(T, c)
  L0 = n**2 * y * (1.d0 - y) * c * chiH + &
    (n**2) * y * (6.1d-19) * T**(-0.63d0) * &
    (1.d0 - exp( -(T/1.d5)**1.63d0 ))
  L = L0
  G = 0.0 !Y así para siempre


  open(20, file='output.dat', status='unknown', action='write')
  write(20,*) '# x [cm] h [erg/g] y rho [g/cm^3] u [cm/s] T [K] P [dyn/cm^2] n [cm^-3]'
  write(20,*) x, h, y, rho, u, T, P, n



  !Iniciamos la iteración
  do while (contador < 2)
    if (T <= 1.0d3) then
        contador = contador + 1
    end if
    call coef_c(T, c)
    call coef_alpha(T, alpha)

    dhdx = (G - L)/(rho*u)
    dydx = (n*y/u)*((1-y)*c - y*alpha)
    
   !Con el dx inicial (a ojo) se calcula el tamaño de los cambios en ese diferencial.
   !Si el cambio es muy grande, se reduce el dx a la mitad, y así sucesivamente hasta que los cambios sean menores al 1%.
    do while (abs(dhdx*dx/h) > 0.01 .or. abs(dydx*dx/y) > 0.01)
        dx = 0.5*dx
    end do
    !Ahora que tenemos un dx como la gente, se calculan las variables
    x = x + dx
    h = h + dhdx*dx
    y = y + dydx*dx


    !Si resulta que el dx inicial era tan pequeño que los cambios son menores al 0.1%, se puede aumentar el dx a la mitad,
    ! y así sucesivamente hasta que los cambios sean mayores al 0.1% (pero menores al 1% por el control anterior).
    if (abs(dhdx*dx/h) < 0.001 .and. abs(dydx*dx/y) < 0.001) then
        dx = 2.0*dx
    end if
    
    !Actualizamos las variables

    rho = ( (5.d0/2.d0)*(P0 + rho0*u0**2) + &
           sqrt( ((5.d0/2.d0)*(P0 + rho0*u0**2))**2 - &
           8.d0*h*rho0**2*u0**2 ) ) / (2.d0*h)

    P = P0 + rho0*u0**2 * (1.d0 - rho0/rho)

    u = rho0*u0 / rho

    n = rho/(mu*mp)
    ne = y*n ! los ionizados = e-
    T = P/((n+ne)*kB)

    L = n**2 * y * (1.d0 - y) * c * chiH + &
    (n**2) * y * (6.1d-19) * T**(-0.63d0) * &
    (1.d0 - exp( -(T/1.d5)**1.63d0 ))


    
    !Escribimos
    write(20,*) x, h, y, rho, u, T, P, n

    !Actualizamos las iniciales
    rho0 = rho
    u0 = u
    T0 = T
    P0 = P
    h0 = h
    y0 = y
    L0 = L
    n0 = n
    


  end do
  close(20)
end subroutine evolve

end program shock