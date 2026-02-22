! ----------------------------------------------------------------------
! Module to classify and compute the coefficients of boundary conditions
! that are of the Neumann type.
! Not required - will only be required if the mesh changes.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 21st September 2025
! ----------------------------------------------------------------------

module bc
    use var_dec
    use fparser
    implicit none

    character(len=256) :: line  ! contains each line read from the input file
    character(len=256), dimension(4) :: bc_values  ! value associated with the boundary condition
    character(len=3), dimension(4) :: bctype  ! value associated with the boundary condition
    integer :: ios  ! I/O status variable
    integer :: bc_count

    real(dp) :: q  ! Temporary flux for Neumann BC
    real(dp) :: delx_e, dely_n, delx_w, dely_s  ! Temporary variables for distances to corresponding cell centers

    integer :: nFuncs
    character(len=2), dimension(5) :: varnames
    real(dp), dimension(5) :: varvals

    contains

        subroutine neu_bc_up()
            
            open(unit=10, file='bc_s_k.txt', status='old', action='read')

            bc_count = 0

            do

                read(10,'(A)', iostat=ios) line
                if (ios /= 0) exit

                line = adjustl(line)  ! remove leading spaces
                if (line(1:1) == '!') cycle  ! skip comment lines
                if (index(line, '!') /= 0) line = line(:index(line,'!') - 1)  ! remove inline comments
                if (trim(line) == '') cycle  ! skip blank lines

                ! Skip irrelevant lines
                if (index(line, 'lx') /= 0) cycle
                if (index(line, 'ly') /= 0) cycle
                if (index(line, 'nx') /= 0) cycle
                if (index(line, 'ny') /= 0) cycle
                if (index(line, 'beta_x') /= 0) cycle
                if (index(line, 'beta_y') /= 0) cycle
                if (index(line, 'k') /= 0) cycle
                if (index(line, 'S_U') /= 0) cycle
                if (index(line, 'S_P') /= 0) cycle
                if (index(line, 'tol') /= 0) cycle
                if (index(line, 'max_iter') /= 0) cycle

                if (index(line, 'bc_') /= 0) then
                    bc_count = bc_count + 1
                    bctype(bc_count) = trim(adjustl(line(7:9)))
                    bc_values(bc_count) = trim(adjustl(line(12:)))
                    cycle
                end if

            end do

            ! nFunc
            nFuncs = 4
            call initf(nFuncs)

            varnames = ["x", "y", "L", "H", "p"]
            call parsef(1,bc_values(1), varnames)
            call parsef(2,bc_values(2), varnames)
            call parsef(3,bc_values(3), varnames)
            call parsef(4,bc_values(4), varnames)

            close(10)

            ! Apply Boundary Conditions
            ! Bottom Boundary (j=1)
            j = 1
            do i = 1, imax - 1

                varvals(1) = xc(i)
                varvals(2) = yc(j)
                varvals(3) = lx
                varvals(4) = ly
                varvals(5) = acos(-1.0_dp) ! pi

                if (bctype(1) .eq. 'neu') then
                    q = evalf(1, varvals) ! Evaluate the Neumann BC expression

                    delx_w = xc(i) - xc(i-1)  ! Distance to west cell center
                    delx_e = xc(i+1) - xc(i)  ! Distance to east cell center
                    dely_n = yc(j+1) - yc(j)  ! Distance to north cell center

                    ! Modify the coefficients for Neumann BC
                    a_w(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_w))
                    a_e(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_e))
                    a_s(i,j) = 0.0_dp
                    a_n(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i,j+1))) * (dx(i)/dely_n))
                    a_p(i,j) = a_e(i,j) + a_w(i,j) + a_n(i,j) + a_s(i,j) + s_p(i,j)
                    s_u(i,j) = s_u(i,j) + q * dx(i)
                end if

            end do

            ! Right Boundary (i=imax)
            i = imax
            do j = 1, jmax - 1

                varvals(1) = xc(i)
                varvals(2) = yc(j)
                varvals(3) = lx
                varvals(4) = ly
                varvals(5) = acos(-1.0_dp) ! pi
                
                if (bctype(2) .eq. 'neu') then
                    q = evalf(2, varvals) ! Evaluate the Neumann BC expression

                    delx_w = xc(i) - xc(i-1)  ! Distance to west cell center
                    dely_s = yc(j) - yc(j-1)  ! Distance to south cell center
                    dely_n = yc(j+1) - yc(j)  ! Distance to north cell center

                    ! Modify the coefficients for Neumann BC
                    a_w(i,j) = ((0.5_dp * (kappa(i-1,j) + kappa(i,j))) * (dy(j)/delx_w))
                    a_e(i,j) = 0.0_dp
                    a_s(i,j) = ((0.5_dp * (kappa(i,j-1) + kappa(i,j))) * (dx(i)/dely_s))
                    a_n(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i,j+1))) * (dx(i)/dely_n))

                    a_p(i,j) = a_e(i,j) + a_w(i,j) + a_n(i,j) + a_s(i,j) + s_p(i,j)
                    s_u(i,j) = s_u(i,j) + q * dy(j)
                end if

            end do

            ! Top Boundary (j=jmax)
            j = jmax
            do i = 1, imax - 1

                varvals(1) = xc(i)
                varvals(2) = yc(j)
                varvals(3) = lx
                varvals(4) = ly
                varvals(5) = acos(-1.0_dp) ! pi
                
                if (bctype(3) .eq. 'neu') then
                    q = evalf(3, varvals) ! Evaluate the Neumann BC expression

                    delx_w = xc(i) - xc(i-1)  ! Distance to west cell center
                    delx_e = xc(i+1) - xc(i)  ! Distance to east cell center
                    dely_s = yc(j) - yc(j-1)  ! Distance to south cell center

                    ! Modify the coefficients for Neumann BC
                    a_w(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_w))
                    a_e(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_e))
                    a_s(i,j) = ((0.5_dp * (kappa(i,j-1) + kappa(i,j))) * (dx(i)/dely_s))
                    a_n(i,j) = 0.0_dp

                    a_p(i,j) = a_e(i,j) + a_w(i,j) + a_n(i,j) + a_s(i,j) + s_p(i,j)
                    s_u(i,j) = s_u(i,j) + q * dx(i)
                end if

            end do

            ! Left Boundary (i=1)
            i = 1
            do j = 1, jmax - 1

                varvals(1) = xc(i)
                varvals(2) = yc(j)
                varvals(3) = lx
                varvals(4) = ly
                varvals(5) = acos(-1.0_dp) ! pi

                if (bctype(4) .eq. 'dir') then
                    q = evalf(4, varvals) ! Evaluate the Neumann BC expression

                    delx_e = xc(i+1) - xc(i)  ! Distance to east cell center
                    dely_s = yc(j) - yc(j-1)  ! Distance to south cell center
                    dely_n = yc(j+1) - yc(j)  ! Distance to north cell center

                    ! Modify the coefficients for Neumann BC
                    a_w(i,j) = 0.0_dp
                    a_e(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_e))
                    a_s(i,j) = ((0.5_dp * (kappa(i,j-1) + kappa(i,j))) * (dx(i)/dely_s))
                    a_n(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i,j+1))) * (dx(i)/dely_n))

                    a_p(i,j) = a_e(i,j) + a_w(i,j) + a_n(i,j) + a_s(i,j) + s_p(i,j)
                    s_u(i,j) = s_u(i,j) + q * dy(j)
                end if

            end do

        end subroutine neu_bc_up
        
end module bc