program test_1D
  
  use solver

  implicit none
  integer, parameter :: ns = 6
  !real(8), parameter :: alpha = 0.5!1.9e-11
  real    :: ymn(ns), ti, tf
  !real(8)    :: y(ns), y_f_paper(ns), epsil(ns), yi(ns)
  real    :: y(ns), yi(ns)
  real    :: epsmn, epsmx, dtmn, tnot, prt
  integer :: itermx
  integer :: i,j

  integer :: n_points

  real(8), allocatable :: r(:), Temp_list(:)
  real(8), allocatable :: y0_HI(:), y0_HeIS(:), y0_HeIM(:)
  real(8), allocatable :: y0_HII(:), y0_HeII(:), y0_e(:)

  ! variables temporales para leer columnas
  real(8) :: rho, v, p, heat, cool
  real(8) :: nhi, nhii, nhei, nheii, nheiii, nheiTR

  real(8) :: T_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc, Flux_loc
  real(8) :: tau_global(3), tau_loc(3), y_new(3)
  real(8) :: y0(3)
  real(8) :: dr
  integer :: ios
  character(len=200) :: line
  integer :: total_lines
  integer :: unit_status 

  open(10,file='Hydro_ioniz_adv.txt',status='old')

  total_lines = 0
  do
    read(10,'(A)',iostat=ios) line
    if (ios /= 0) exit
    total_lines = total_lines + 1
  end do

  close(10)

  n_points = total_lines - 100 !descarto los primeros 100 puntos (mucho ruido)

  allocate(r(n_points), Temp_list(n_points))
  allocate(y0_HI(n_points), y0_HeIS(n_points), y0_HeIM(n_points))
  allocate(y0_HII(n_points), y0_HeII(n_points), y0_e(n_points))

  open(10,file='Hydro_ioniz_adv.txt',status='old')
  do j=1,100
    read(10,*,iostat=ios)
    if (ios /= 0) exit
  end do

  do j=1,n_points
   read(10,*,iostat=ios) r(j), rho, v, p, Temp_list(j), heat, cool
     if (ios /= 0) then
      print*, 'EOF reached at j = ', j
      exit
     end if
  end do

  open(11,file='Ion_species_adv.txt',status='old')
  do j=1,100
    read(11,*,iostat=ios)
    if (ios /= 0) exit
  end do

  do j=1,n_points
    read(11,*) r(j), nhi, nhii, nhei, nheii, nheiii, nheiTR

    y0_HI(j)   = nhi
    y0_HeIS(j) = nhei
    y0_HeIM(j) = nheiTR*1e5
    y0_HII(j)  = nhii
    y0_HeII(j) = nheii
    y0_e(j)    = nhii + nheii 

  end do
  close(11)
  
  !Flujo de la estrella que llega al planeta. Necesario para los phi's
  Flux_loc = !520055.63 ! [photons cm^-2 s^-1] ????

  !Tamaño de cada dr. Necesario para calcular los tau's
  dr = (r(n_points) - r(1))/real(n_points)

  !inicializacion del tau global. Necesario para calcular el tau del primer (atras para delante) dr
  tau_global = (/ 0., 0., 0. /)

  unit_status = 20

  open(unit_status, file='status.log', status='replace', action='write')

  write(unit_status,*) '# STATUS FILE FOR test_1D_1'
  write(unit_status,*) '# alpha        = ', alpha
  write(unit_status,*) '# total_lines  = ', total_lines
  write(unit_status,*) '# n_points     = ', n_points
  write(unit_status,*) '# Flux         = ', Flux_loc
  write(unit_status,*) '# dr           = ', dr
  write(unit_status,*) '#'
  write(unit_status,*) '# Columns: r  T  phiHI  phiHeIS  phiHeIM  tau_HI  tau_HeIS  tau_HeIM'
  write(unit_status,*)
 
  do j=1,n_points

    ti=0.0
    tf = 1000.0
  
    epsmn = 1e-5!1e-5
    epsmx = 0.0
    dtmn = 0.0
    tnot = ti
    itermx = 5
    ymn(:) = 1e-20
    prt = 0.0

    call chemsp(epsmn, epsmx, dtmn, tnot, itermx ,ns, ymn, prt) 

    T_loc = Temp_list(n_points - j + 1) !arranca por el ultimo
    y0 = (/ y0_HI(n_points - j + 1), y0_HeIS(n_points - j + 1), y0_HeIM(n_points - j + 1) /)

    call init_local_taus(y0, dr, tau_global, tau_loc)

    tau_global = tau_loc !actualizo el tau global para el siguiente dr

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
    yi(1) = y0_HI(n_points - j + 1)
    yi(2) = y0_HII(n_points - j + 1)
    yi(3) = y0_HeIS(n_points - j + 1)
    yi(4) = y0_HeIM(n_points - j + 1)
    yi(5) = y0_HeII(n_points - j + 1)
    yi(6) = y0_e(n_points - j + 1)!yi(2) + yi(5)
  

  !Copy initial values to y
    do i  =1,ns
      y(i) = yi(i)
    end do

  !Solve the system
    call chemeq2solve(tf, y, ns, T_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc, j)

  !solo me interesan las especies neutras para actualizar los taus  
    !y_new = (/ y(1), y(3), y(4) /) 
  !Actualiza los valores de tau con los y's de chemeq y lo acumulo con los de los anteriores dr 
    !call update_global_taus(y_new, dr, tau_global)

    write(unit_status,'(ES12.4,7ES15.6)') r(n_points-j+1), T_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc, &
                       tau_global(1), tau 
    !flush(unit_status)

  end do !inicial 

  close(unit_status)

 
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
