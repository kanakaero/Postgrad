! ----------------------------------------------------------------------
! Subroutine to discretize the 2D Convection Diffusion equation using
! the Finite Volume Method.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 12th November 2025
! ----------------------------------------------------------------------

subroutine discret()
    use var_dec
    implicit none

    real(dp) :: delx_e, delx_w, dely_n, dely_s  ! Distances to corresponding cell centers
    real(dp), allocatable :: F_w(:, :), F_e(:, :), F_s(:, :), F_n(:, :)  ! Array of Mass Fluxes at cell faces
    real(dp), allocatable :: D_w(:,:), D_e(:,:), D_s(:,:), D_n(:,:)      ! Array of Diffusion Conductance at cell faces
    real(dp) :: rho, Gamma      ! Density and Diffusion Coefficient
    real(dp) :: Delta_F         ! Temporary variable for flux difference
    real(dp), allocatable :: Pe(:,:), pec(:,:)    ! Peclet Number array    

    allocate(F_w(imax, jmax), F_e(imax, jmax), F_s(imax, jmax), F_n(imax, jmax))
    allocate(D_w(imax, jmax), D_e(imax, jmax), D_s(imax, jmax), D_n(imax, jmax))
    allocate(a_p(imax, jmax), a_e(imax, jmax), a_w(imax, jmax), a_n(imax, jmax), a_s(imax, jmax))
    allocate(Pe(imax, jmax))

    rho = 1.0_dp                ! Density
    Gamma = 1.0_dp / 50.0_dp    ! Diffusion Coefficient

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

            Pe(i,j) = abs((F_e(i,j) - F_w(i,j)) * delx_e / (Gamma * dy(j)) + (F_n(i,j) - F_s(i,j)) * dely_n / (Gamma * dx(i)))  ! Peclet Number

! ----------------------------------------------------------------------

        end do
    end do

    open(unit=30, file='Pe.tec', status='replace')

    ! Write Tecplot header
    write(30, '(A)') 'TITLE = "Convection Diffusion Initial Conditions"'
    write(30, '(A)') 'VARIABLES = "X", "Y", "Pe"'
    write(30, '(A,I0,A,I0,A)') 'ZONE T="Zone 1", I=', imax, ', J=', jmax, ', DATAPACKING=POINT, ZONETYPE=ORDERED'

    ! Write data in point format: each line is one point (X, Y, T)
    do j = 1, jmax
        do i = 1, imax
            write(30,'(ES15.6,1X,ES15.6,1X,ES15.6)') &
            xc(i), yc(j), Pe(i,j)
        end do
    end do

    close(30)

    print *, 'Peclet Number field written to Pe.tec'

    allocate(pec(imax, jmax))
    pec = 0.0_dp
    
    do i = 2, imax-1
        do j = 2, jmax-1
            if (Pe(i,j) > 2) then
                pec(i,j) = 1.0_dp
            else
                pec(i,j) = 0.0_dp
            end if
        end do
    end do

    open(unit=30, file='Pec.tec', status='replace')

    ! Write Tecplot header
    write(30, '(A)') 'TITLE = "Convection Diffusion Initial Conditions"'
    write(30, '(A)') 'VARIABLES = "X", "Y", "Pe"'
    write(30, '(A,I0,A,I0,A)') 'ZONE T="Zone 1", I=', imax, ', J=', jmax, ', DATAPACKING=POINT, ZONETYPE=ORDERED'

    ! Write data in point format: each line is one point (X, Y, T)
    do j = 1, jmax
        do i = 1, imax
            write(30,'(ES15.6,1X,ES15.6,1X,ES15.6)') &
            xc(i), yc(j), pec(i,j)
        end do
    end do

    close(30)

    print *, 'Binary Peclet Number field written to Pec.tec'

    if (allocated(Pe)) deallocate(Pe)
    if (allocated(F_w)) deallocate(F_w)
    if (allocated(F_e)) deallocate(F_e)
    if (allocated(F_s)) deallocate(F_s)
    if (allocated(F_n)) deallocate(F_n)
    if (allocated(D_w)) deallocate(D_w)
    if (allocated(D_e)) deallocate(D_e)
    if (allocated(D_s)) deallocate(D_s)
    if (allocated(D_n)) deallocate(D_n)

end subroutine discret