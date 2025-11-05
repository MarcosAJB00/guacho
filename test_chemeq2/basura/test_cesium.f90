program cesium
  use solver        ! módulo donde está chemeq2 y chemsp
  implicit none

  integer, parameter :: ns = 7, na = 5, mxcase = 9
  real :: y(ns), yi(ns), yf(ns), ymin_local(ns)
  real :: epsil(ns), eps(mxcase)
  real :: ti, tf, deltat, tscale, cput, epsmin_local, rms
  integer :: i, icase, inlp
  external :: csdfe

  !--------------------------------------------
  ! Inicialización de parámetros
  !--------------------------------------------
  ymin_local = 1.0d-20
  eps = (/ 1.0d-1, 5.0d-2, 1.0d-2, 5.0d-3, 1.0d-3, 5.0d-4, 1.0d-4, 5.0d-5, 1.0d-5 /)
  tscale = 1.0d0 / 1024.0d0
  inlp   = 1
  ti     = 0.0d0
  tf     = 1000.0d0
  deltat = (tf - ti) / inlp

  !--------------------------------------------
  ! Valores iniciales y finales (del paper)
  !--------------------------------------------
  yi = (/ 5.20d2, 6.20d2, 1.0d12, 0.0d0, 3.6d14, 1.4d15, 1.0d2 /)
  yf = (/ 2.59139492061d4, 7.55718460300d4, 1.53194051722d3, 9.99999923516d11, &
         3.59000000051d14, 1.4d15, 4.96578968239d4 /)

  print *, '==============================================================='
  print *, '   CHEMEQ2 Cesium test problem (Section 10.2 example)'
  print *, '==============================================================='

  !--------------------------------------------
  ! Bucle principal de casos (epsmin distintos)
  !--------------------------------------------
  do icase = 1, mxcase
     epsmin_local = eps(icase)
     print *, '---------------------------------------------------------------'
     print *, 'CASE ', icase, '   epsmin =', epsmin_local

     ! Llamada a chemsp usando variable local
     call chemsp(epsmin_local, 0.0e0, 0.0e0, ti, ymin_local, 1)

     ! Reset iniciales
     y = yi
     cput = 0.0d0

     ! Integración (dtime/CPU se omite)
     call chemeq2solve(deltat, y, na)

     ! Calcular densidad electrónica final
     y(7) = y(2) - y(1)

     ! Errores relativos
     do i = 1, ns
        epsil(i) = abs(y(i) - yf(i)) / max(1.0d-20, min(y(i), yf(i)))
     end do

     rms = sqrt(sum(epsil(1:ns-1)**2) / real(ns-1,8))

     ! Mostrar resultados resumidos
     print '(A)', ' Species   Y_init        Y_final       Y_sol         RelErr'
     do i = 1, ns
        print '(I3,3E15.6,E10.3)', i, yi(i), yf(i), y(i), epsil(i)
     end do
     print *, 'RMS Error = ', rms
  end do

end program cesium

