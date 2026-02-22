! ----------------------------------------------------------------------
! Subroutine to solve the discretized 2D Convection Diffusion equation 
! using the Gauss Seidel iterative method.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 14th November 2025
! ----------------------------------------------------------------------

subroutine gauss_seidel()
    use var_dec

    real(dp) :: delta_T, f, rho, U_A, k_by_cp  ! temporary variables
    real(dp) :: dTdx, dTdy   ! Temperature gradients in x and y directions

    U_A = 1.0_dp  ! Inlet Velocity
    rho = 1.0_dp  ! Density
    k_by_cp = 1.0_dp / 50.0_dp  ! Thermal conductivity divided by specific heat capacity
    ! since c_p is not given we will calculate heat flux per c_p

    max_iter = 10000
    tol = 1.0e-8_dp

    allocate(R(max_iter))
    allocate(qx(imax, jmax))
    allocate(qy(imax, jmax))
    allocate(qmag(imax, jmax))

    print *, 'Beginning Gauss Seidel Solver'
    print *, ' '

    do iter = 1, max_iter
        print *, 'Solver Iteration: ', iter
        R(iter) = 0.0_dp   ! Initialize the residual array
        do i = 2, imax-1
            do j = 2, jmax-1

                if (a_p(i,j) /= a_p(i,j)) then
                    print *, 'NaN found in a_p at i,j =', i, j
                    stop
                end if
                if (a_p(i,j) <= 0.0_dp) then
                    print *, 'Non-positive a_p at i,j =', i, j, ' a_p=', a_p(i,j)
                    ! Optionally set to a tiny positive floor to continue:
                    a_p(i,j) = max(a_p(i,j), 1.0e-14_dp)
                end if

                ! Interior nodes
                T(i,j) = (a_e(i,j)*T(i+1,j) + a_w(i,j)*T(i-1,j) + a_n(i,j)*T(i,j+1) + a_s(i,j)*T(i,j-1)) / a_p(i,j)
            end do
        end do

        ! Calculate the residual
        do i = 2, imax-1
            do j = 2, jmax-1
                R(iter) = R(iter) + abs(-a_p(i,j)*T(i,j) + a_e(i,j)*T(i+1,j) + a_w(i,j)*T(i-1,j) + a_n(i,j)*T(i,j+1) + a_s(i,j)*T(i,j-1))
            end do
        end do

        delta_T = abs(T(1,1) - T(imax,1)) ! Inlet - Outlet temperature difference
        f = rho * U_A * delta_T

        R(iter) = R(iter) / f

        print *, 'Residual = ', R(iter)

        ! Apply Temperature at boundaries (Neumann)
        do i = 1, imax
            T(i,1) = T(i,2)
            T(i,jmax) = T(i,jmax-1)
        end do

        do j = 1, jmax
            if (yc(j) < (1.0_dp - 0.068_dp)) then
                T(1,j) = T(2,j)
            end if
            if (yc(j) < 0.068_dp) then
                T(imax,j) = T(imax-1,j)
            end if
        end do

        if (R(iter) < tol) exit

    end do

    ! Write out the Residuals
    open(unit=20, file='residuals_gsd.dat', status='replace')
    do i = 1, iter
        write(20,*) i, R(i)
    end do
    close(20)
    print *, 'Solver Converged in ', iter, ' iterations.'

    ! Calculate heat fluxes in x and y directions with second-order central difference
    do i = 2, imax-1
    do j = 2, jmax-1
        dTdx = (T(i+1,j) - T(i-1,j)) / (2.0d0 * dx(i))
        dTdy = (T(i,j+1) - T(i,j-1)) / (2.0d0 * dy(j))

        qx(i,j) = - k_by_cp * dTdx
        qy(i,j) = - k_by_cp * dTdy
        qmag(i,j) = sqrt(qx(i,j)**2 + qy(i,j)**2)
    end do
    end do

    open(unit=30, file='output_gsd.tec', status='replace')

    ! Write Tecplot header
    write(30, '(A)') 'TITLE = "2d Diffusion Data"'
    write(30, '(A)') 'VARIABLES = "X", "Y", "T", "QX", "QY", "QMAG"'
    write(30, '(A,I0,A,I0,A)') 'ZONE T="Zone 1", I=', imax, ', J=', jmax, ', DATAPACKING=POINT, ZONETYPE=ORDERED'

    ! Write data in point format: each line is one point (X, Y, T)
    do j = 1, jmax
        do i = 1, imax
            if (abs(qx(i,j)) < 1.0d-100) qx(i,j) = 0.0d0
            if (abs(qy(i,j)) < 1.0d-100) qy(i,j) = 0.0d0
            if (abs(qmag(i,j)) < 1.0d-100) qmag(i,j) = 0.0d0

            write(30,'(ES15.6,1X,ES15.6,1X,ES15.6,1X,ES15.6,1X,ES15.6,1X,ES15.6)') &
            xc(i), yc(j), T(i,j), qx(i,j), qy(i,j), qmag(i,j)
        end do
    end do

    close(30)

    print *, 'Temperature field written to output_gsd.tec'
    print *, ' '

    print *, 'Solver Finished!'
    print *, ' '

    T_gsd = T  ! Store Gauss Seidel solution

    if (allocated(R)) deallocate(R)
    if (allocated(qx)) deallocate(qx)
    if (allocated(qy)) deallocate(qy)
    if (allocated(qmag)) deallocate(qmag)

    if (allocated(a_p)) deallocate(a_p)
    if (allocated(a_w)) deallocate(a_w)
    if (allocated(a_e)) deallocate(a_e)
    if (allocated(a_s)) deallocate(a_s)
    if (allocated(a_n)) deallocate(a_n)
    if (allocated(T)) deallocate(T)

end subroutine gauss_seidel