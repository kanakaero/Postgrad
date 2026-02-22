! ----------------------------------------------------------------------
! Subroutine to discretize the 2D Diffusion equation using
! the Finite Volume Method.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 21st September 2025
! ----------------------------------------------------------------------

subroutine discret()
    use var_dec
    implicit none

    real(dp) :: delx_e, delx_w, dely_n, dely_s  ! Distances to corresponding cell centers

! ---------------- Distances to cell centres ---------------------------

    ! Loop over all internal cells to compute coefficients
    do j = 2, jmax-1
        do i = 2, imax-1

            ! Distances to neighbouring cell centers
            delx_w = xc(i) - xc(i-1)  ! Distance to west cell center
            delx_e = xc(i+1) - xc(i)  ! Distance to east cell center
            dely_s = yc(j) - yc(j-1)  ! Distance to south cell center
            dely_n = yc(j+1) - yc(j)  ! Distance to north cell center

! ----------------------------------------------------------------------

! ---------------- Coefficents at cell centres -------------------------

            ! Array of Coefficients at neighbouring cell centers
            ! Linear interpolation of kappa at cell faces
            a_w(i,j) = ((0.5_dp * (kappa(i-1,j) + kappa(i,j))) * (dy(j)/delx_w))
            a_e(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_e))
            a_s(i,j) = ((0.5_dp * (kappa(i,j-1) + kappa(i,j))) * (dx(i)/dely_s))
            a_n(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i,j+1))) * (dx(i)/dely_n))

            a_p(i,j) = a_w(i,j) + a_e(i,j) + a_s(i,j) + a_n(i,j) ! Sum of coefficients at P
            ! Unknown part of source term is added to a_p during the Solver step

            T(1,i,j) = 1.0_dp  ! Initial guess for temperature at cell center for 1st iteration
            T(2,i,j) = 1.0_dp  ! Initial guess for temperature at cell center for 2nd iteration

! ----------------------------------------------------------------------

        end do
    end do

end subroutine discret