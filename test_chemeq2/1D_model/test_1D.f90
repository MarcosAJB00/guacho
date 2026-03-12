program test_1D
  
  use solver

  implicit none
  integer, parameter :: ns = 6
  real    :: ymn(ns), ti, tf
  !real(8)    :: y(ns), y_f_paper(ns), epsil(ns), yi(ns)
  real    :: y(ns), yi(ns)
  real    :: epsmn, epsmx, dtmn, tnot, prt
  integer :: itermx
  integer :: i,j

  integer, parameter :: n_points = 10

  real(8) :: Temp_list(n_points)
  real(8) :: y0_HI(n_points), y0_HeIS(n_points), y0_HeIM(n_points)
  real(8) :: T_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc, Flux_loc
  real(8) :: tau_global(3), tau_loc(3), y_new(3)
  real(8) :: y0(3), y0_HI_loc, y0_HeIS_loc, y0_HeIM_loc
  real(8) :: dr

  !Perfil de temperatura. Necesario para chemeq
  Temp_list = (/ 2000., 5000., 7000., 6500., 5500., &
               5000., 4000., 3700., 3500., 3400. /)

  !Perfil de abundancias neutras iniciales. Necesario para calcular los tau's iniciales, phi's y para chemeq
  y0_HI = (/ 1., 1., 1., 1., 1., &
               1., 1., 1., 1., 1. /)
  y0_HeIS = (/ 1., 1., 1., 1., 1., &
               1., 1., 1., 1., 1. /)
  y0_HeIM = (/ 1., 1., 1., 1., 1., &
               1., 1., 1., 1., 1. /)

  !aca tambien van los perfiles de los ionizados y e. Necesarios para chemeq
  
  !Flujo de la estrella que llega al planeta. Necesario para los phi's
  Flux_loc = 1e10

  !Tamaño de cada dr. Necesario para calcular los tau's
  dr = 1.0

  !inicializacion del tau global. Necesario para calcular el tau del primer (atras para delante) dr
  tau_global = (/ 0., 0., 0. /)
 
  do j=1,n_points

    ti=0.0
    tf = 1000.0
  
    epsmn = 1e-5
    epsmx = 0.0
    dtmn = 0.0
    tnot = ti
    itermx = 5
    ymn(:) = 1e-20
    prt = 0.0

    call chemsp(epsmn, epsmx, dtmn, tnot, itermx ,ns, ymn, prt) 


    T_loc = Temp_list(n_points - j + 1) !arranca por el ultimo
    y0_HI_loc = y0_HI(n_points - j + 1)
    y0_HeIS_loc = y0_HeIS(n_points - j + 1)
    y0_HeIM_loc = y0_HeIM(n_points - j + 1)

    y0 = (/ y0_HI_loc, y0_HeIS_loc, y0_HeIM_loc /)

    call init_local_taus(y0, dr, tau_global, tau_loc)

    call phis_rate(Flux_loc,tau_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc)


  !Numero de las especies
  !iHI- ---> 1
  !iHII ---> 2
  !iHeIS --> 3
  !iHeIM --> 4
  !iHeII --> 5
  ! ie ----> 6

  !Inicializacion de las especies
  !call nr_init(y, y0)???
    yi(1) = y0_HI_loc!1.498e9
    yi(2) = 0.01!2.55e6
    yi(3) = y0_HeIS_loc!1.669e8
    yi(4) = y0_HeIM_loc!3.17e-3
    yi(5) = 0.01!3.17e-1
    yi(6) = 0.0!yi(2) + yi(5)
  

  !Copy initial values to y
    do i  =1,ns
      y(i) = yi(i)
    end do

  !Solve the system
    call chemeq2solve(tf, y, ns, T_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc, j)

  !solo me interesan las especies neutras para actualizar los taus  
    y_new = (/ y(1), y(3), y(4) /) 
  !Actualiza los valores de tau con los y's de chemeq y lo acumulo con los de los anteriores dr 
    call update_global_taus(y_new, dr, tau_global)

  end do !inicial 

 
contains

subroutine phis_rate(Flux, taus, phi_H, phi_HeIS, phi_HeIM)
  !esta subrutina calcula el valor de phi necesarios para el chemeq solver en funcion del tau local 
  !calculado a partir de las abundancias iniciales

  implicit none
  real(8), intent(in)  :: Flux, taus(3)
  real(8), intent(out) :: phi_H, phi_HeIS, phi_HeIM 
  real(8) :: a0_H, a0_HeIS, a0_HeIM
  
  a0_H    = 3.87d-18 !seccion eficaz HI
  a0_HeIS = 4.03d-18
  a0_HeIM = 3.59d-18

  
  phi_H    = a0_H*Flux*exp(-taus(1))
  phi_HeIS = a0_HeIS*Flux*exp(-taus(2))
  phi_HeIM = a0_HeIM*Flux*exp(-taus(3))

end subroutine phis_rate


subroutine init_local_taus(y0, dr, taus_global, taus_loc)
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

end subroutine init_local_taus


subroutine update_global_taus(y_new, dr, taus_global)
  !Esta subrutina actualiza el tau global. El tau global es la suma de los taus de cada dr calculado a partir de los y's finales
  !optenido por chemeq en cada dr

  implicit none
  real(8), intent(in) :: y_new(3)
  real(8), intent(in)  :: dr !tamaño de cada celda, en principio constante
  real(8), intent(inout) :: taus_global(3)
  real(8) :: a0_H, a0_HeIS, a0_HeIM

  a0_H    = 3.87d-18 !seccion eficaz HI
  a0_HeIS = 4.03d-18 !
  a0_HeIM = 3.59d-18 !

  taus_global(1) = a0_H*y_new(1)*dr    + taus_global(1)
  taus_global(2) = a0_HeIS*y_new(2)*dr + taus_global(2)
  taus_global(3) = a0_HeIM*y_new(3)*dr + taus_global(3)

end subroutine update_global_taus


end program test_1D
