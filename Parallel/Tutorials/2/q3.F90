! Trapezoidal Rule Integration in Fortran

program trapezoidal_integration
    implicit none
    integer :: n, i
    real(8) :: a, b, h, integral, x

    ! Limits of Integration and Number of Subintervals
    a = 1
    b = 4 * atan (1.0_8)
    n = 32

    ! Width of each subinterval
    h = (b - a) / dble(n)

    ! Initialize the integral
    integral = (f(a) + f(b)) / 2.0d0

    ! Sum the function values at each subinterval
    do i = 1, n - 1
        x = a + dble(i) * h
        integral = integral + f(x)
    end do

    ! Multiply by the width of the subintervals
    integral = integral * h

    ! Print the result
    print *, 'The integral from ', a, ' to ', b, ' is approximately: ', integral
    print *, ' '
    print *, 'Error in the integral is approximately: ', abs(integral - 0.198557d0)

contains

    ! Function to be integrated
    real(8) function f(x)
        real(8), intent(in) :: x
        f = sin(x) / (2.0d0*x**3)
    end function f

end program trapezoidal_integration