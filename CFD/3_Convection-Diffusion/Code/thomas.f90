! ----------------------------------------------------------------------
! Subroutine to solve Ax = d using the Thomas algorithm (TDMA)
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 21st November 2025
! ----------------------------------------------------------------------

! a – lower diagonal (a(1) unused)
! b – main diagonal
! c – upper diagonal (c(n) unused)
! d – RHS (overwritten)
! x – solution

module thomas_module
    implicit none
contains

subroutine thomas(n, a, b, c, d, x)
    
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: n
    real(dp), intent(inout) :: a(:), b(:), c(:), d(:)
    real(dp), intent(out) :: x(:)
    integer :: i
    real(dp), allocatable :: p(:)

    allocate(p(n))

    ! forward elimination
    p(1) = c(1) / b(1)
    d(1)  = d(1) / b(1)

    do i = 2, n
        b(i) = b(i) - a(i) * p(i-1)
        p(i) = c(i) / b(i)
        d(i) = (d(i) - a(i) * d(i-1)) / b(i)
    end do

    ! back substitution
    x(n) = d(n)
    do i = n-1, 1, -1
        x(i) = d(i) - p(i) * x(i+1)
    end do

    deallocate(p)

end subroutine thomas

end module thomas_module