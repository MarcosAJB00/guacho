program shock 

real(kind=8) :: G, L, rho, u, T, ny, y0, dhdx, dydx
real(kind =8) :: u0, rho0, mu, mp, n0, nHII, T0, chiH, kB, P0, h0

chiH = 2.179d-11 ! erg, energia de ionizacion de H

u0 = 1e6 ! cm/s
rho0 = 1e-12 ! g/cm^3
mu = 1.0 ! masa media, 1 para H, 1.3 para 10% de He
mp = 1.67d-24 ! g
n0 = mu*mp*rho0
nHII0 = 0.076*n0 ! ionizacion inicial
y0 = nHII0/n0
T0 = 1e4

kB = 1.38d-16 ! erg/K
P0 = n0*kB*T0
h0 = (u0**2)/2 + (5/2)*(P0/rho0) !entalpia inicial

call evolve(h0, y0, rho0, u0, T0, P0)

contains


subroutine c(T)
  real(kind=8) :: T
  c = 5.83d-11*(T**0.5) * exp(-157800/T)
end subroutine c

subroutine alpha(T)
  real(kind=8) :: T
  alpha = 3.69d-10*(T**-0.79d0)
end subroutine alpha

subroutine evolve(h0, y0, rho0, u0, T0, P0)
  real(kind=8) :: x, dx
  real(kind=8) :: h, y,
  real(kind=8) :: G, L, rho, u, T, ny, dhdx, dydx, L0
  
  x = 0.0
  dx = 1 ! cm

  h=h0
  y=y0
  rho = rho0
  u = u0
  T = T0
  P = P0
  n = n0
  ! Calculamos la pérdida de energía inicial
  call c(T)
  L0 = n**2 * y * (1.d0 - y) * c * chiH + &
    (n**2) * y * (6.1d-19) * T**(-0.63d0) * &
    (1.d0 - exp( -(T/1.d5)**1.63d0 ))
  L = L0
  G = 0.0 !Y así para siempre


  open(20, file='output.dat', status='unknown', action='write')
  write(20,*) x, h, y, rho, u, T, P, n
  close(20)


  !Iniciamos la iteración
  while (T<1e4)

    call c(T)
    call alpha(T)

    dhdx = (G - L)/(rho*u)
    dydx = (n*y/u)*((1-y)*c - y*alpha)
    
   !Con el dx inicial (a ojo) se calcula el tamaño de los cambios en ese diferencial.
   !Si el cambio es muy grande, se reduce el dx a la mitad, y así sucesivamente hasta que los cambios sean menores al 1%.
    while (abs(dhdx*dx/h) > 0.01 or abs(dydx*dx/y) > 0.01)
        dx *= 0.5
    end while
    !Ahora que tenemos un dx como la gente, se calculan las variables
    x = x + dx
    h = h + dhdx*dx
    y = y + dydx*dx


    !Si resulta que el dx inicial era tan pequeño que los cambios son menores al 0.1%, se puede aumentar el dx a la mitad,
    ! y así sucesivamente hasta que los cambios sean mayores al 0.1% (pero menores al 1% por el control anterior).
    if (abs(dhdx*dx/h) < 0.001 and abs(dydx*dx/y) < 0.001):
        dx *= 2.0*dx
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
    open(20, file='output.dat', status='unknown', action='write')
    write(20,*) x, h, y, rho, u, T, P, n
    close(20)

    !Actualizamos las iniciales
    rho0 = rho
    u0 = u
    T0 = T
    P0 = P
    h0 = h
    y0 = y
    L0 = L
    n0 = n
    


  end while
  
end subroutine evolve

end program shock