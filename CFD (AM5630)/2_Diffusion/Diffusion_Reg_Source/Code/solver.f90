! ----------------------------------------------------------------------
! Subroutine to solve the discretized 2D Diffusion equation using
! the Gauss Seidel iterative method.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 28th September 2025
! ----------------------------------------------------------------------

subroutine solver()
    use var_dec

    real(dp) :: total_source  ! total source term for normalization
    real(dp) :: dTdx, dTdy   ! Temperature gradients in x and y directions

    print *, 'Beginning Solver'
    print *, ' '

    do iter = 1, max_iter
        print *, 'Solver Iteration: ', iter
        R(iter) = 0.0_dp   ! Initialize the residual array
        total_source = 0.0_dp ! Initialize total source term
        do i = 2, imax-1
            do j = 2, jmax-1
                ! Interior nodes
                T(i,j) = (a_e(i,j)*T(i+1,j) + a_w(i,j)*T(i-1,j) + a_n(i,j)*T(i,j+1) + a_s(i,j)*T(i,j-1) + (s_u(i,j)*dx(i)*dy(j))) / a_p(i,j)
                total_source = total_source + abs(s_u(i,j))
            end do
        end do
        ! Calculate the residual
        do i = 2, imax-1
            do j = 2, jmax-1
                R(iter) = R(iter) + abs(a_p(i,j)*T(i,j) - &
                             (a_e(i,j)*T(i+1,j) + a_w(i,j)*T(i-1,j) + &
                              a_n(i,j)*T(i,j+1) + a_s(i,j)*T(i,j-1) + &
                              (s_u(i,j)*dx(i)*dy(j))))/(total_source*lx*ly)
            end do
        end do

        print *, 'Residual = ', R(iter)

        if (R(iter) < tol) exit

    end do

    ! Write out the Residuals
    open(unit=20, file='residuals.dat', status='replace')
    do i = 1, iter
        write(20,*) i, R(i)
    end do
    close(20)
    print *, 'Solver Converged in ', iter, ' iterations.'

    ! Calculate heat fluxes in x and y directions with second-order central difference
    do i = 2, imax-1
        do j = 2, jmax-1
            dTdx = (T(i+1,j) - T(i-1,j)) / (2.0 * dx(i))
            dTdy = (T(i,j+1) - T(i,j-1)) / (2.0 * dy(j))
            qx(i,j) = -kappa(i,j) * dTdx
            qy(i,j) = -kappa(i,j) * dTdy
            qmag(i,j) = sqrt(qx(i,j)**2 + qy(i,j)**2)
        end do
    end do



    open(unit=30, file='output.tec', status='replace')

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

    print *, 'Temperature field written to output.tec'
    print *, ' '

    print *, 'Solver Finished!'

end subroutine solver