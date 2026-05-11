program test_1D
  
  use solver
  implicit none

!====================================================
! INPUTS
!====================================================
  character(len=50)  :: label
  character(len=200) :: perfil_densidad

  real(8) :: T, Rp, time
  real(8) :: FluxHI, FluxHeIS, FluxHeIM
  real(8) :: He_H_ratio
  real(8) :: h_ion_frac, he_ion_frac, he_metastable_frac

!====================================================
! GENERALES
!====================================================
  integer, parameter :: ns = 6
  integer :: i, j, ios, n_points, total_lines, itermx, unit_status

  real(kind=8)    :: ymn(ns), ti, tf, y(ns), yi(ns)
  real(kind=8)    :: epsmn, epsmx, dtmn, tnot, prt

  real(kind=8), allocatable :: r(:), rho(:)
  real(kind=8), allocatable :: y0_HI(:), y0_HII(:), y0_HeIS(:)
  real(kind=8), allocatable :: y0_HeIM(:), y0_HeII(:), y0_e(:)

  real(kind=8) :: mp, R_jupiter, r_planet, dr
  real(kind=8) :: T_loc, phiHI_loc, phiHeIS_loc, phiHeIM_loc
  real(kind=8) :: tau_global(3), tau_loc(3), y0(3)
  real(kind=8) :: t1, t2
  !real(kind=8) :: y_new(3)

  character(len=200) :: line

!====================================================
! LEER INPUTS
!====================================================
  open(99,file='inputs.dat',status='old',action='read')

  read(99,*) label, T
  read(99,*) label, Rp
  read(99,*) label, FluxHI
  read(99,*) label, FluxHeIS
  read(99,*) label, FluxHeIM
  read(99,*) label, He_H_ratio
  read(99,*) label, h_ion_frac
  read(99,*) label, he_ion_frac
  read(99,*) label, he_metastable_frac
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

!====================================================
  call cpu_time(t1)

!====================================================
! CONTAR LINEAS PERFIL DENSIDAD
!====================================================
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

!====================================================
! LEER PERFIL
!====================================================
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

!====================================================
! CONSTANTES
!====================================================
  mp        = 1.67d-24
  R_jupiter = 7.1492d9
  r_planet  = 1.38 * R_jupiter
  dr        = r_planet*(r(n_points)-r(1))/real(n_points,8)
  T_loc = T

!====================================================
! PERFILES INICIALES
!====================================================
  open(30,file='perfiles_iniciales.txt',status='replace')

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

!====================================================
! STATUS FILE
!====================================================
  tau_global = (/ 0., 0., 0. /)

  unit_status = 20
  open(unit_status,file='status.log',status='replace')

  write(unit_status,*) '# STATUS FILE FOR test_1D_1'
  write(unit_status,*) '# total_lines = ', total_lines
  write(unit_status,*) '# n_points    = ', n_points
  write(unit_status,*) '# FluxHI      = ', FluxHI
  write(unit_status,*) '# FluxHeIS    = ', FluxHeIS
  write(unit_status,*) '# FluxHeIM    = ', FluxHeIM
  write(unit_status,*) '# dr          = ', dr/R_jupiter, ' R_jupiter'
  write(unit_status,*) '# r_P         = ', r_planet/R_jupiter, ' R_jupiter'
  write(unit_status,*) '#'
  write(unit_status,*) '# Columns: r T phiHI phiHeIS phiHeIM tauHI tauHeIS tauHeIM'
  write(unit_status,*)

!====================================================
! LOOP PRINCIPAL
!====================================================
  do j=1,n_points

    print *, 'Running cell ', j, ' of ', n_points

    ti = 0.0
    tf = real(time)

    epsmn = 1e-5
    epsmx = 0.0
    dtmn  = 0.0
    tnot  = ti
    itermx = 5
    ymn(:) = 1e-20
    prt = 1.0

    call chemsp(epsmn, epsmx, dtmn, tnot, itermx, ns, ymn, prt)

    y0 = (/ y0_HI(n_points-j+1), y0_HeIS(n_points-j+1), y0_HeIM(n_points-j+1) /)

    call init_local_taus(y0, dr, tau_global, tau_loc)
    
    tau_global = tau_loc !actualizo el tau global para el siguiente dr

    call phis_rate(FluxHI, FluxHeIS, FluxHeIM, tau_loc, &
                   phiHI_loc, phiHeIS_loc, phiHeIM_loc)

    yi(1)=y0_HI(n_points-j+1)
    yi(2)=y0_HII(n_points-j+1)
    yi(3)=y0_HeIS(n_points-j+1)
    yi(4)=y0_HeIM(n_points-j+1)
    yi(5)=y0_HeII(n_points-j+1)
    yi(6)=y0_e(n_points-j+1)

    do i=1,ns
        y(i)=yi(i)
    end do

    call chemeq2solve(tf,y,ns,T_loc,phiHI_loc,phiHeIS_loc,phiHeIM_loc,j)
    
    !prueba
    !y_new = (/ y(1), y(3), y(4) /) !solo me interesan estas especies para el tau
    !call update_global_taus(y_new, dr, tau_global)
    !tau_global = tau_global  
  
    !escribo el status de cada dr
    write(unit_status,'(ES12.4,7ES15.6)') r(n_points-j+1), T_loc, &
         phiHI_loc, phiHeIS_loc, phiHeIM_loc, &
         tau_global(1), tau_global(2), tau_global(3)

  end do

  close(unit_status)

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
