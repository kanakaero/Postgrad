! ----------------------------------------------------------------------
! This module contains the global variables and parameters used.
! The double precision is explicitly enforced here rather than 
! enforcing it at the compiler level.
! Any variable that is to be used globally should be declared here and 
! allocated using the allocate.f90 subroutine if it is allocatable.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 12th November 2025
! ----------------------------------------------------------------------

module var_dec

    integer, parameter :: dp = kind(1.d0)   ! guarantees precision up to that of the machine-compiler-specific double precision real

! ----------------------------------------------------------------------

    integer :: i, j, iter  ! loop counters

! ------------------ Grid and Initial Conditions -----------------------

    integer :: nx, ny       ! number of nodes in x and y directions
    integer :: imax, jmax   ! nx and ny -> + 2 for boundaries
    
    real(dp) :: T_inlet, T_outlet  ! Inlet and Outlet Temperatures

    real(dp), allocatable :: xc(:), yc(:)   ! Array of cell coordinates in x and y directions
    real(dp), allocatable :: dx(:), dy(:)   ! Array of grid spacing in x and y directions
    real(dp), allocatable :: xg(:), yg(:)   ! Array of grid lines in x and y directions

    real(dp), allocatable :: u_flat(:), v_flat(:)  ! Array of flattened velocity components in x and y directions
    real(dp), allocatable :: u(:, :), v(:, :), U_mag(:, :)  ! Array of velocity components in x and y directions and magnitude

! ----------------------- Discretization -------------------------------

    real(dp), allocatable :: a_p(:, :), a_e(:, :), a_w(:, :), a_n(:, :), a_s(:, :)   ! Coefficients of the discretized equation
    real(dp), allocatable :: T(:, :)              ! Temperature array
    real(dp), allocatable :: T_gsd(:, :)          ! Gauss seidel Solution
    real(dp), allocatable :: T_tdma(:, :)         ! TDMA Solution
    real(dp), allocatable :: qx(:,:), qy(:,:), qmag(:,:) ! Heat fluxes in x and y directions

! ----------------------------------------------------------------------

! --------------------- Gauss Seidel Solver ----------------------------

    real(dp), allocatable :: R(:)  ! Residual array
    real(dp) :: tol       ! Tolerance limit for convergence
    integer :: max_iter   ! Maximum number of iterations

! ----------------------------------------------------------------------
    
end module var_dec

