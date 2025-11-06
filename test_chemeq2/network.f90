module network
  implicit none
contains    

  subroutine gsub(y, q, d, t)
    implicit none
    real, intent(in)    :: y(:)
    real, intent(out)   :: q(size(y)), d(size(y))
    real, intent(in)    :: t
    real :: o2m, csp, cs, cso2, o2, n2, ne
    real :: cr1, cr2, cr3, cr_4, cr5, cr6, cr7
    !Numero de las especies
    !o2- ---> 1
    !cs+ ---> 2
    !cs ----> 3
    !cso2 --> 4
    !o2 ----> 5
    !n2 ----> 6
    ! e ----> 7

    !utilze local storage for variables
    o2m  = y(1)
    csp  = y(2)
    cs   = y(3)
    cso2 = y(4)
    o2   = y(5)
    n2   = y(6)
  !  ne   = y(7)
    ne   = max(csp - o2m, 0.0)  ! to ensure electron density is always positive

    !calculate reaction rates
    cr1 = 5.0e-8*o2m*csp
    cr2 = 1.0e-12*csp*ne
    cr3 = 3.24e-3*cs
    cr_4 = 4.0e-1*o2m
    cr5 = 1.0e-31*o2*cs*(cs + cso2 + n2 + o2)
    cr6 = 1.24e-30*o2*o2*ne
    cr7 = 1.0e-31*o2*n2*ne 

    if (t >= 700.0) then
      cr_4 = 0.0
      cr6 = 0.0
      cr7 = 0.0
    end if
    !calculate production rates (q(i)) and loss rates (d(i))

    q(1) = cr6 + cr7
    d(1) = cr1 + cr_4

    q(2) = cr3
    d(2) = cr1 + cr2

    q(3) = cr1 + cr2
    d(3) = cr3 + cr5

    q(4) = cr5 - 1.0e-30*o2*cs*cso2
    d(4) = - 1.0e-30*o2*cs*cso2

    q(5) = cr1 + cr_4
    d(5) = cr5 + cr6 + cr7

    q(6) = 0.0  !agregado por mi
    d(6) = 0.0  !agregado por mi

    q(7) = 0.0  !agregado por mi
    d(7) = 0.0  !agregado por mi

  end subroutine gsub

end module network