module network
  implicit none
contains
  !===========================================================
  subroutine csdfe(y, q, d, t)
    implicit none
    real(8), intent(inout) :: y(:)
    real(8), intent(out)   :: q(size(y)), d(size(y))
    real(8), intent(in)    :: t
    real(8) :: o2m, csp, cs, cso2, o2, n2, ne
    real(8) :: cr1, cr2, cr3, cr4, cr5, cr6, cr7

    ! Variables locales
    o2m  = y(1)
    csp  = y(2)
    cs   = y(3)
    cso2 = y(4)
    o2   = y(5)
    n2   = y(6)

    ne = max(csp - o2m, 0.0d0)
    y(7) = ne

    ! Reacciones
    cr1 = 5.00d-08 * o2m * csp
    cr2 = 1.00d-12 * csp * ne
    cr3 = 3.24d-03 * cs
    cr4 = 4.00d-01 * o2m
    cr5 = 1.00d-31 * o2 * cs * (cs + cso2 + n2 + o2)
    cr6 = 1.24d-30 * o2 * o2 * ne
    cr7 = 1.00d-31 * o2 * n2 * ne

    ! Formación y pérdida
    q(1) = cr6 + cr7
    d(1) = cr1 + cr4

    q(2) = cr3
    d(2) = cr1 + cr2

    q(3) = cr1 + cr2
    d(3) = cr3 + cr5

    q(4) = cr5 - 1.0e-31*o2*cs*cso2
    d(4) = -1.0e31*o2*cs*cso2

    q(5) = cr1 + cr4
    d(5) = cr5 + cr6 + cr7
end subroutine csdfe


end module network

