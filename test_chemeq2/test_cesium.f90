program test_cesium
  
  use solver

  implicit none
 ! real    :: epsmin, sqreps, epscl, epsmax, dtmin, tstart
 ! integer :: itermax
  integer, parameter :: n_spec = 7
  real    :: ymin_loc(n_spec)
  real    :: y(n_spec),tnot,dtg, yf(n_spec), y_f_paper(n_spec)

  tnot=0.0
  call chemsp(epsmin, epsmax, dtmin,tnot, ymin_loc, itermax) !agregada por el Ale
  print*, 'epsmin=', epsmin, ' epsmax=', epsmax, ' dtmin=', dtmin, ' itermax=', itermax

  !Numero de las especies
  !o2- ---> 1
  !cs+ ---> 2
  !cs ----> 3
  !cso2 --> 4
  !o2 ----> 5
  !n2 ----> 6
  ! e ----> 7

  !Inicializacion de las especies

  y(1) = 5.2e2
  y(2) = 6.2e2
  y(3) = 1e12
  y(4) = 0.0
  y(5) = 3.6e14
  y(6) = 1.4e15
  y(7) = 1e2

  !Tiempo total de integracion 
  dtg = 1000.0

  call chemeq2solve(dtg, y, n_spec)
  
  !Calculate final electron density from densities of other charges species
  y(7) = y(2) - y(1)
  
  !Final results
  yf(:) = y(:)


  !Values from paper for comparison
  y_f_paper(1) = 2.5913e4
  y_f_paper(2) = 7.5571e4
  y_f_paper(3) = 1.5319e3
  y_f_paper(4) = 1.00e12
  y_f_paper(5) = 3.590e14
  y_f_paper(6) = 1.400e15
  y_f_paper(7) = 4.9665e4

  print *, 'Resultados a t =', dtg, ' s:'
  print '(A8,1PE12.4)', 'O2-  =', yf(1), '  (paper:', y_f_paper(1),')'
  print '(A8,1PE12.4)', 'Cs+  =', yf(2), '  (paper:', y_f_paper(2),')'
  print '(A8,1PE12.4)', 'Cs   =', yf(3), '  (paper:', y_f_paper(3),')'
  print '(A8,1PE12.4)', 'CsO2 =', yf(4), '  (paper:', y_f_paper(4),')'
  print '(A8,1PE12.4)', 'O2   =', yf(5), '  (paper:', y_f_paper(5),')'
  print '(A8,1PE12.4)', 'N2   =', yf(6), '  (paper:', y_f_paper(6),')'
  print '(A8,1PE12.4)', 'e-   =', yf(7), '  (paper:', y_f_paper(7),')'

end program test_cesium

  !o2- ---> 1
  !cs+ ---> 2
  !cs ----> 3
  !cso2 --> 4
  !o2 ----> 5
  !n2 ----> 6
  ! e ----> 7
 
