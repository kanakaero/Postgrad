!Using the fourth-order Pade scheme and the third-order boundary schemes discussed in the class differntiate the following fucntion using 15 equally spaced points and plot the result obtained and compare the same with the exact derivative. The TDMA algorithm discussed in the class needs to be implmented only in serial. Compare these results with a CDS2 in the interior and FD1 and BD1 schemes at the boundaries.

!  f(x) = sin(5x), 0 ≤ x ≤ 3

program q3
    implicit none
    integer, parameter :: N = 15
    real(8) :: x(N), f(N), df_exact(N)
    real(8) :: df_pade(N), df_cds(N)
    real(8) :: a(N), b(N), c(N), d(N)
    real(8) :: dx, m
    integer :: i

    dx = 3.0d0 / (N - 1)

    ! Grid + exact solution
    do i = 1, N
        x(i) = (i-1) * dx
        f(i) = sin(5.0d0 * x(i))
        df_exact(i) = 5.0d0 * cos(5.0d0 * x(i))
    end do

    ! ---- Boundary (left) ----
    b(1) = 1.0d0
    c(1) = 2.0d0
    d(1) = (-2.5d0*f(1) + 2.0d0*f(2) + 0.5d0*f(3)) / dx

    ! ---- Interior ----
    do i = 2, N-1
        a(i) = 0.25d0
        b(i) = 1.0d0
        c(i) = 0.25d0
        d(i) = (3.0d0/(4.0d0*dx)) * (f(i+1) - f(i-1))
    end do

    ! ---- Boundary (right) ----
    a(N) = 2.0d0
    b(N) = 1.0d0
    d(N) = (2.5d0*f(N) - 2.0d0*f(N-1) - 0.5d0*f(N-2)) / dx

    ! ---- TDMA (forward sweep) ----
    do i = 2, N
        m = a(i) / b(i-1)
        b(i) = b(i) - m * c(i-1)
        d(i) = d(i) - m * d(i-1)
    end do

    ! ---- Back substitution ----
    df_pade(N) = d(N) / b(N)
    do i = N-1, 1, -1
        df_pade(i) = (d(i) - c(i) * df_pade(i+1)) / b(i)
    end do

    ! ---- CDS2 + FD1/BD1 ----
    df_cds(1) = (f(2) - f(1)) / dx
    df_cds(N) = (f(N) - f(N-1)) / dx

    do i = 2, N-1
        df_cds(i) = (f(i+1) - f(i-1)) / (2.0d0 * dx)
    end do

    ! ---- Output ----
    open(unit=10, file="output.dat", status="replace")
    write(10,*) "# x   exact   pade   cds"

    do i = 1, N
        write(10,'(4f12.6)') x(i), df_exact(i), df_pade(i), df_cds(i)
    end do

    close(10)

    print *, "Data written to output.dat"

end program q3