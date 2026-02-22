! ----------------------------------------------------------------------
! Subroutine to solve the discretized 2D Convection Diffusion equation 
! using the TDMA iterative method.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 21st November 2025
! ----------------------------------------------------------------------

subroutine tdma()
    use var_dec
    use thomas_module

    real(dp) :: delta_T, f, U_A, k_by_cp  ! temporary variables
    real(dp) :: dTdx, dTdy   ! Temperature gradients in x and y directions
    real(dp) :: H  ! Height of the domain

    real(dp) :: delx_e, delx_w, dely_n, dely_s  ! Distances to corresponding cell centers
    real(dp), allocatable :: F_w(:, :), F_e(:, :), F_s(:, :), F_n(:, :)  ! Array of Mass Fluxes at cell faces
    real(dp), allocatable :: D_w(:,:), D_e(:,:), D_s(:,:), D_n(:,:)      ! Array of Diffusion Conductance at cell faces
    real(dp) :: rho, Gamma      ! Density and Diffusion Coefficient
    real(dp) :: Delta_F         ! Temporary variable for flux difference
    real(dp), allocatable :: A_x(:), B_x(:), C_x(:), D_x(:), X_x(:) ! TDMA arrays - x-direction
    real(dp), allocatable :: A_y(:), B_y(:), C_y(:), D_y(:), X_y(:) ! TDMA arrays - y-direction

    allocate(F_w(imax, jmax), F_e(imax, jmax), F_s(imax, jmax), F_n(imax, jmax))
    allocate(D_w(imax, jmax), D_e(imax, jmax), D_s(imax, jmax), D_n(imax, jmax))
    allocate(a_p(imax, jmax), a_e(imax, jmax), a_w(imax, jmax), a_n(imax, jmax), a_s(imax, jmax))
    allocate(qx(imax, jmax))
    allocate(qy(imax, jmax))
    allocate(qmag(imax, jmax))

    allocate(A_x(imax-2), B_x(imax-2), C_x(imax-2), D_x(imax-2), X_x(imax-2)) ! TDMA arrays for x-direction
    allocate(A_y(jmax-2), B_y(jmax-2), C_y(jmax-2), D_y(jmax-2), X_y(jmax-2)) ! TDMA arrays for y-direction

    rho = 1.0_dp                ! Density
    Gamma = 1.0_dp / 50.0_dp    ! Diffusion Coefficient

    H = 2 ! Height of the domain

    U_A = 1.0_dp  ! Inlet Velocity
    rho = 1.0_dp  ! Density
    k_by_cp = 1.0_dp / 50.0_dp  ! Thermal conductivity divided by specific heat capacity
    ! since c_p is not given we will calculate heat flux per c_p

    tol = 1.0e-8_dp
    max_iter = 1000

! ---------------- Distances to cell centres ---------------------------

    ! Loop over all internal cells to compute coefficients
    do j = 2, jmax-1
        do i = 2, imax-1

            ! Distances to neighbouring cell centers
            delx_w = xc(i) - xc(i-1)  ! Distance to west cell center
            delx_e = xc(i+1) - xc(i)  ! Distance to east cell center
            dely_s = yc(j) - yc(j-1)  ! Distance to south cell center
            dely_n = yc(j+1) - yc(j)  ! Distance to north cell center

            ! Mass Fluxes at cell faces
            F_w(i,j) = rho * u(i-1,j) * dy(j)   ! West face
            F_e(i,j) = rho * u(i,j) * dy(j)     ! East face
            F_s(i,j) = rho * v(i,j-1) * dx(i)   ! South face
            F_n(i,j) = rho * v(i,j) * dx(i)     ! North face

            ! Diffusion Conductance at cell faces
            D_w(i,j) = (Gamma/delx_w) * dy(j)   ! West face
            D_e(i,j) = (Gamma/delx_e) * dy(j)   ! East face
            D_s(i,j) = (Gamma/dely_s) * dx(i)   ! South face
            D_n(i,j) = (Gamma/dely_n) * dx(i)   ! North face
            
! ----------------------------------------------------------------------

! ---------------- Coefficents at cell centres -------------------------

            ! Array of Coefficients at neighbouring cell centers
            ! Using Hybrid Differencing Scheme - Max doesn't allow 3 parameters in Fortran
            a_w(i,j) = max(max(F_w(i,j), D_w(i,j) + F_w(i,j)/2.0_dp), 0.0_dp )
            a_e(i,j) = max(max(-F_e(i,j), D_e(i,j) - F_e(i,j)/2.0_dp), 0.0_dp )
            a_s(i,j) = max(max(F_s(i, j), D_s(i,j) + F_s(i,j)/2.0_dp), 0.0_dp)
            a_n(i,j) = max(max(-F_n(i, j), D_n(i,j) - F_n(i,j)/2.0_dp), 0.0_dp)

            Delta_F = (F_e(i,j) - F_w(i,j)) + (F_n(i,j) - F_s(i,j)) ! Flux difference

            a_p(i,j) = a_w(i,j) + a_e(i,j) + a_s(i,j) + a_n(i,j) + Delta_F ! Sum of coefficients at P

! ----------------------------------------------------------------------

        end do
    end do

    print *, 'Beginning TDMA Solver'
    print *, ' '

    call bc()  ! Apply Boundary Conditions

    allocate(R(max_iter))

    do iter = 1, max_iter

        ! Reinforce Dirichlet Boundary Conditions
        do j = 1, jmax
            if (yc(j) >= (1.0_dp - 0.068_dp)*H) then
                T(1,j) = T_inlet
            end if
            if (yc(j) >= 0.068_dp*H) then
                T(imax,j) = T_outlet
            end if
        end do

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
        
        ! Initialize the TDMA arrays
        A_x(:) = 0.0_dp
        B_x(:) = 0.0_dp
        C_x(:) = 0.0_dp
        D_x(:) = 0.0_dp
        X_x(:) = 0.0_dp

        ! Set up and solve TDMA in x-direction for each row
        do j = 2, jmax-1
            do i = 2, imax-1
                A_x(i-1) = -a_w(i,j)
                B_x(i-1) = a_p(i,j)
                C_x(i-1) = -a_e(i,j)
                D_x(i-1) = a_n(i,j)*T(i,j+1) + a_s(i,j)*T(i,j-1)
            end do

            D_x(1) = D_x(1) + a_w(2,j)*T(1,j)       ! Adjust for left boundary
            D_x(imax-2) = D_x(imax-2) + a_e(imax-1,j)*T(imax,j) ! Adjust for right boundary

            call thomas(imax-2, A_x, B_x, C_x, D_x, X_x)
            do i = 2, imax-1
                T(i,j) = X_x(i-1)
            end do
        end do

        ! Initialize the TDMA arrays
        A_y(:) = 0.0_dp
        B_y(:) = 0.0_dp
        C_y(:) = 0.0_dp
        D_y(:) = 0.0_dp
        X_y(:) = 0.0_dp

        ! Set up and solve TDMA in y-direction for each column
        do i = 2, imax-1
            do j = 2, jmax-1
                A_y(j-1) = -a_s(i,j)
                B_y(j-1) = a_p(i,j)
                C_y(j-1) = -a_n(i,j)
                D_y(j-1) = a_e(i,j)*T(i+1,j) + a_w(i,j)*T(i-1,j)
            end do

            D_y(1) = D_y(1) + a_s(i,2)*T(i,1)       ! Adjust for bottom boundary
            D_y(jmax-2) = D_y(jmax-2) + a_n(i,jmax-1)*T(i,jmax) ! Adjust for top boundary

            call thomas(jmax-2, A_y, B_y, C_y, D_y, X_y)
            do j = 2, jmax-1
                T(i,j) = X_y(j-1)
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

        if (R(iter) < tol) exit

    end do

    ! Write out the Residuals
    open(unit=20, file='residuals_tdma.dat', status='replace')
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

    open(unit=30, file='output_tdma.tec', status='replace')

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

    print *, 'Temperature field written to output_tdma.tec'
    print *, ' '

    print *, 'Solver Finished!'

    print *, ' '

    T_tdma = T  ! Store TDMA solution

    if (allocated(R)) deallocate(R)
    if (allocated(qx)) deallocate(qx)
    if (allocated(qy)) deallocate(qy)
    if (allocated(qmag)) deallocate(qmag)

    if (allocated(F_w)) deallocate(F_w)
    if (allocated(F_e)) deallocate(F_e)
    if (allocated(F_s)) deallocate(F_s)
    if (allocated(F_n)) deallocate(F_n)
    if (allocated(D_w)) deallocate(D_w)
    if (allocated(D_e)) deallocate(D_e)
    if (allocated(D_s)) deallocate(D_s)
    if (allocated(D_n)) deallocate(D_n)

end subroutine tdma