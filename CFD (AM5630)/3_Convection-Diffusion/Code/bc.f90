! ----------------------------------------------------------------------
! Subroutine to apply boundary conditions for the 2D Convection 
! Diffusion equation 
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 14th November 2025
! ----------------------------------------------------------------------

! -------------------------------------------------------------------------
!                                 V_D (inlet)  
!                                  |
!          +-------------------------------------------+
!          |                                           |
!          |                 h_D <---->                |
!          |                                           |
!  U_A --> |   h_A                                     | --> U_C
!          |                                           |
!          |                                           |
!          |                                           |
!          |                                           |
!          |                                           |
!          |                                           |
!          |                                           |
!  <-- U_B |   h_B                                h_C  | --> U_C
!          |                                           |
!          |                                           |
!          +-------------------------------------------+
!
!   Coordinate system:
!
!           y ^
!             |
!             |
!             +-----> x
!
!   Rectangle dimensions:
!       Height = H = 2
!       Length = L = 2
!
!   Boundary velocities:
!       U_A : inlet on left wall (upper part) : 1 
!       U_B : outlet on left wall (lower part) : 0
!       U_C : outlet on right wall : 1
!       V_D : inlet from top wall : 0
!
!   Boundary thicknesses:
!       h_A/H : height of upper left inlet : 0.068
!       h_B/H : height of lower left inlet --> doesn't matter as U_B = 0
!       h_C/H : height of right outlet : 0.068
!       h_D/L : width of top inlet --> doesn't matter as V_D = 0
!
!   Temperature Boundary Conditions:
!       Left Wall : T_A = 20 (degree C)
!       Right Wall : T = 50 (degree C) --> apart from outlet
!       Rest of the walls are Adiabatic (Neumann BC: dT/dn = 0)
!       
!    Note: When prescribing a heat flux we must divide by c_p since the equation is
!          in terms of temperature and not energy.
!
! -------------------------------------------------------------------------

subroutine bc()
    use var_dec
    implicit none
    real(dp) :: D_w, D_e, D_s, D_n  ! Diffusion Conductance at cell faces
    real(dp) :: delx_e, delx_w, dely_n, dely_s  ! Distances to corresponding cell centers
    real(dp) :: Gamma      ! Density and Diffusion Coefficient

    real(dp) :: H ! Height of the domain

    allocate(T(imax, jmax))
    
    ! Initialize Temperature field
    T = 0.0_dp

    H = 2 ! Height of the domain - given

    Gamma = 1.0_dp / 50.0_dp    ! Diffusion Coefficient

    ! Apply Dirichlet Boundary Conditions
    do j = 1, jmax
        ! Left Wall Dirichlet
        if (yc(j) >= (1.0_dp - 0.068_dp) * H) then ! upper part
            T(1,j) = T_inlet         ! inlet
        end if

        ! Right Wall Dirichlet
        if (yc(j) >= 0.068_dp * H) then ! lower part
            T(imax,j) = T_outlet ! outlet
        end if
    end do

    ! Apply Neumann Boundary Conditions (Adiabatic) - Implicit

    ! Top Wall
    do i = 2, imax-1
        dely_n = yc(jmax) - yc(jmax-1) ! Distance to north cell center
        D_n = (Gamma/dely_n) * dx(i)   ! North face
        a_n(i,jmax-1) = 0.0_dp      ! No contribution from north face
        a_p(i,jmax-1) = a_p(i, jmax-1) + D_n ! Update a_p
    end do

    ! Bottom Wall
    do i = 2, imax-1
        dely_s = yc(2) - yc(1)  ! Distance to south cell center
        D_s = (Gamma/dely_s) * dx(i)   ! South face
        a_s(i,2) = 0.0_dp          ! No contribution from south face
        a_p(i,2) = a_p(i,2) + D_s      ! Update a_p
    end do

    ! Left Wall (for j < 0.068*H)
    do j = 2, jmax-1
        if (yc(j) < (1.0_dp - 0.068_dp) * H) then
            delx_w = xc(2) - xc(1)  ! Distance to west cell center
            D_w = (Gamma/delx_w) * dy(j)   ! West face
            a_w(2,j) = 0.0_dp          ! No contribution from west face
            a_p(2,j) = a_p(2,j) + D_w      ! Update a_p
        end if
    end do

    ! Right Wall (for j < 0.068*H)
    do j = 2, jmax-1
        if (yc(j) < 0.068_dp * H) then
            delx_e = xc(imax) - xc(imax-1)  ! Distance to east cell center
            D_e = (Gamma/delx_e) * dy(j)   ! East face
            a_e(imax-1,j) = 0.0_dp      ! No contribution from east face
            a_p(imax-1,j) = a_p(imax-1,j) + D_e  ! Update a_p
        end if
    end do

end subroutine bc