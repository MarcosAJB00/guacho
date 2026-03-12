program test_cesium
  
  use solver

  implicit none
  integer, parameter :: ns = 7
  real    :: ymn(ns), ti, tf
  real    :: y(ns), y_f_paper(ns), epsil(ns), yi(ns)
  real    :: epsmn, epsmx, dtmn, tnot, prt
  integer :: itermx
  integer :: i 

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

  !Numero de las especies
  !o2- ---> 1
  !cs+ ---> 2
  !cs ----> 3
  !cso2 --> 4
  !o2 ----> 5
  !n2 ----> 6
  ! e ----> 7

  !Inicializacion de las especies
  yi(1) = 5.2e2
  yi(2) = 6.2e2
  yi(3) = 1e12
  yi(4) = 0.0
  yi(5) = 3.6e14
  yi(6) = 1.4e15
  yi(7) = 1e2

  !Copy initial values to y
  do i  =1,ns
    y(i) = yi(i)
  end do

  !Solve the system
  call chemeq2solve(tf, y, ns)
  
  !Calculate final electron density from densities of other charges species
  y(7) = y(2) - y(1)

  !Values from paper for comparison
  y_f_paper(1) = 2.5913e4
  y_f_paper(2) = 7.5571e4
  y_f_paper(3) = 1.5319e3
  y_f_paper(4) = 1.00e12
  y_f_paper(5) = 3.590e14
  y_f_paper(6) = 1.400e15
  y_f_paper(7) = 4.9665e4

  !Calculate relative error
  do i=1,ns
    epsil(i) = (abs(y(i) - y_f_paper(i)))/y_f_paper(i)
  end do

  print *, 'Resultados a t =', tf, ' s:'
  ! Cabecera (solo una vez)
  print '(A6, 3X, A12, 3X, A12, 3X, A12, 3X, A12)', &
      'Especie', 'y inicial', 'yf (calc)', 'yf (paper)', 'error_relativo'
  print '(A48)', '------------------------------------------------'

  ! Filas (una por especie)
  print '(A6, 4(3X,1PE12.4))', 'O2-  ', yi(1), y(1), y_f_paper(1), epsil(1)
  print '(A6, 4(3X,1PE12.4))', 'Cs+  ', yi(2), y(2), y_f_paper(2), epsil(2)
  print '(A6, 4(3X,1PE12.4))', 'Cs   ', yi(3), y(3), y_f_paper(3), epsil(3)
  print '(A6, 4(3X,1PE12.4))', 'CsO2 ', yi(4), y(4), y_f_paper(4), epsil(4)
  print '(A6, 4(3X,1PE12.4))', 'O2   ', yi(5), y(5), y_f_paper(5), epsil(5)
  print '(A6, 4(3X,1PE12.4))', 'N2   ', yi(6), y(6), y_f_paper(6), epsil(6)
  print '(A6, 4(3X,1PE12.4))', 'e-   ', yi(7), y(7), y_f_paper(7), epsil(7)
  print '(A48)', '------------------------------------------------'

end program test_cesium

 
