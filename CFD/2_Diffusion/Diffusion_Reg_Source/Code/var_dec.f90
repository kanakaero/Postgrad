! ----------------------------------------------------------------------
! This module contains the global variables and parameters used.
! The double precision is explicitly enforced here rather than 
! enforcing it at the compiler level.
! Any variable that is to be used globally should be declared here and 
! allocated using the allocate.f90 subroutine if it is allocatable.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 28th September 2025
! ----------------------------------------------------------------------

module var_dec

    integer, parameter :: dp = kind(1.d0)   ! guarantees precision up to that of the machine-compiler-specific double precision real
    INTEGER, PARAMETER :: rn = KIND(1.0d0)          ! Precision of real numbers
    INTEGER, PARAMETER :: is = SELECTED_INT_KIND(1) ! Data type of bytecode

! ----------------------------------------------------------------------

    integer :: i, j, iter  ! loop counters

! ------------------------- Meshing ------------------------------------
    
    real(dp) :: nx, ny      ! number of nodes in x and y directions
    integer :: imax, jmax   ! integer versions of nx and ny - + 2 for boundaries
    integer :: im, jm       ! integer versions of nx and ny
    real(dp) :: lx, ly      ! length of the domain in x and y directions
    real(dp) :: beta_x, beta_y  ! skewness factors in x and y directions

    real(dp), allocatable :: dx(:), dy(:)   ! Array of grid spacing in x and y directions
    real(dp), allocatable :: xc(:), yc(:)   ! Array of cell centres in x and y directions
    real(dp), allocatable :: xg(:), yg(:)   ! Array of grid lines in x and y directions

    real(dp) :: a_x, a_y    ! First Cell length from Geometric Progression in x and y directions

! ----------------------------------------------------------------------

! ----------------------- Discretization -------------------------------

    real(dp), allocatable :: a_p(:, :), a_e(:, :), a_w(:, :), a_n(:, :), a_s(:, :)   ! Coefficients of the discretized equation
    real(dp), allocatable :: s_u(:,:), s_p(:,:)   ! Source term - known part and unknown part
    real(dp), allocatable :: kappa(:, :)          ! Thermal Conductivity array
    real(dp), allocatable :: T(:, :)              ! Temperature array

    real(dp), allocatable :: qx(:,:), qy(:,:), qmag(:,:) ! Heat fluxes in x and y directions

! ----------------------------------------------------------------------

! Boundary Conditions Variable are defined in the bc_s_k.f90 subroutine.
! These are not required to be declared globally.

! --------------------- Gauss Seidel Solver ----------------------------

    real(dp), allocatable :: R(:)  ! Residual array
    real(dp) :: tol       ! Tolerance limit for convergence
    integer :: max_iter   ! Maximum number of iterations

! ----------------------------------------------------------------------

end module var_dec

