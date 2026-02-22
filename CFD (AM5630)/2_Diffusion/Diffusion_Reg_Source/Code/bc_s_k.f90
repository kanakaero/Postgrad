! ----------------------------------------------------------------------
! Subroutine to input the boundary conditions, source terms, solver parameters and 
! spatially varying thermal conductivities from an input file. Aditionally,
! it also allocates the required arrays.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 28th September 2025
! ----------------------------------------------------------------------

subroutine bc_s_k()
    use var_dec
    use fparser

    implicit none

    character(len=256) :: line  ! contains each line read from the input file
    character(len=256) :: bcval_b, bcval_r, bcval_t, bcval_l  ! value associated with the boundary condition
    character(len=3) :: bctype_b, bctype_r, bctype_t, bctype_l ! value associated with the boundary condition
    character(len=100) :: source_value_S_P, source_value_S_U  ! value associated with the source term
    character(len=100) :: k_value  ! value associated with the source term

    real(dp) :: q  ! Temporary flux for Neumann BC
    real(dp) :: delx_e, dely_n, delx_w, dely_s  ! Temporary variables for distances to corresponding cell centers

    integer :: ios ! iostat variable for error handling
    integer :: bc_count

    ! For parsing expressions
    integer :: nFuncs
    character(len=2), dimension(5) :: varnames
    real(dp), dimension(5) :: varvals

    open(unit=10, file='bc_s_k.txt', status='old', action='read')

    print *, 'Reading input file bc_s_k.txt'

    ! deallocate arrays
    if (allocated(kappa)) deallocate(kappa) 
    if (allocated(T)) deallocate(T) 
    if (allocated(R)) deallocate(R)

    ! nFunc
    nFuncs = 7
    call initf(nFuncs)

    ! allocate arrays    
    allocate(kappa(imax, jmax))  ! Thermal conductivity array
    allocate(s_u(imax, jmax))    ! Source term array - known part
    allocate(s_p(imax, jmax))    ! Source term array - unknown part
    allocate(T(imax, jmax))      ! Temperature array
    allocate(qx(imax,jmax), qy(imax,jmax), qmag(imax,jmax)) ! Heat flux arrays

    T(:,:) = 0.0_dp   ! Initialize the temperature field

    bc_count = 0 ! initialize boundary condition counter

    rewind(10) ! rewind the file to read from the beginning

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

        ! Thermal conductivity expression
        if (index(line, 'k') /= 0) then
            k_value = trim(adjustl(line(3:)))
            cycle
        end if

        ! Source terms known part
        if (index(line, 'S_U') /= 0) then
            source_value_S_U = trim(adjustl(line(5:)))
            cycle
        end if

        ! Source terms unknown part
        if (index(line, 'S_P') /= 0) then
            source_value_S_P = trim(adjustl(line(5:)))
            cycle
        end if

        ! Boundary conditions (store all lines starting with 'bc_')
        if (index(line, 'bc_') /= 0) then
            select case (trim(adjustl(line(4:5))))
            case ('b')  ! bottom
                bctype_b = trim(adjustl(line(7:9)))
                bcval_b  = trim(adjustl(line(12:)))
            case ('r')  ! right
                bctype_r = trim(adjustl(line(7:9)))
                bcval_r  = trim(adjustl(line(12:)))
            case ('t')  ! top
                bctype_t = trim(adjustl(line(7:9)))
                bcval_t  = trim(adjustl(line(12:)))
            case ('l')  ! left
                bctype_l = trim(adjustl(line(7:9)))
                bcval_l  = trim(adjustl(line(12:)))
            end select
            cycle
        end if

        if (index(line,'tol') /= 0) read(line(4:), *) tol
        if (index(line,'max_iter') /= 0) read(line(9:), *) max_iter

    end do

    varnames = ["x", "y", "L", "H", "p"]
    call parsef(1,k_value, varnames)
    call parsef(2,source_value_S_U, varnames)
    call parsef(3,source_value_S_P, varnames)
    call parsef(4, bcval_b, varnames)
    call parsef(5, bcval_r, varnames)
    call parsef(6, bcval_t, varnames)
    call parsef(7, bcval_l, varnames)

    do j = 1, jmax
        do i = 1, imax

            varvals(1) = xc(i)
            varvals(2) = yc(j)
            varvals(3) = lx
            varvals(4) = ly
            varvals(5) = acos(-1.0_dp) ! pi

            kappa(i,j) = evalf(1, varvals)
            s_u(i,j) = evalf(2, varvals)
            s_p(i,j) = evalf(3, varvals)

        end do
    end do

    ! Apply Boundary Conditions
    ! Bottom Boundary (j=1)
    j = 1

    if (bctype_b .eq. 'dir') then
        do i = 1, imax

            varvals(1) = xc(i)
            varvals(2) = yc(j)
            varvals(3) = lx
            varvals(4) = ly
            varvals(5) = acos(-1.0_dp) ! pi

            T(i,j) = evalf(4, varvals) ! Evaluate the Dirichlet BC expression

            ! Modify the coefficients for Dirichlet BC
            a_w(i,j) = 0.0_dp
            a_e(i,j) = 0.0_dp
            a_s(i,j) = 0.0_dp
            a_n(i,j) = 0.0_dp
            a_p(i,j) = 1.0_dp

        end do

    else

        do i = 1, imax - 1

            varvals(1) = xc(i)
            varvals(2) = yc(j)
            varvals(3) = lx
            varvals(4) = ly
            varvals(5) = acos(-1.0_dp) ! pi

            q = evalf(4, varvals) ! Evaluate the Neumann BC expression

            delx_w = xc(i) - xc(i-1)  ! Distance to west cell center
            delx_e = xc(i+1) - xc(i)  ! Distance to east cell center
            dely_n = yc(j+1) - yc(j)  ! Distance to north cell center

            ! Modify the coefficients for Neumann BC
            a_w(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_w))
            a_e(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_e))
            a_s(i,j) = 0.0_dp
            a_n(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i,j+1))) * (dx(i)/dely_n))
            a_p(i,j) = a_e(i,j) + a_w(i,j) + a_n(i,j) + a_s(i,j) - (s_p(i,j)*dx(i)*dy(j))

            qx(i,j) = 0.0_dp ! heat flux in x direction
            qy(i,j) = q ! heat flux normal to the x direction
            qmag(i,j) = sqrt(qx(i,j)**2 + qy(i,j)**2) ! magnitude of heat flux

        end do

    end if

    ! Right Boundary (i=imax)
    i = imax

    if (bctype_r .eq. 'dir') then

        do j = 1, jmax

            varvals(1) = xc(i)
            varvals(2) = yc(j)
            varvals(3) = lx
            varvals(4) = ly
            varvals(5) = acos(-1.0_dp) ! pi
            

            T(i,j) = evalf(5, varvals) ! Evaluate the Dirichlet BC expression

            ! Modify the coefficients for Dirichlet BC
            a_w(i,j) = 0.0_dp
            a_e(i,j) = 0.0_dp
            a_s(i,j) = 0.0_dp
            a_n(i,j) = 0.0_dp
            a_p(i,j) = 1.0_dp

        end do

    else

        do j = 1, jmax - 1

            varvals(1) = xc(i)
            varvals(2) = yc(j)
            varvals(3) = lx
            varvals(4) = ly
            varvals(5) = acos(-1.0_dp) ! pi

            q = evalf(5, varvals) ! Evaluate the Neumann BC expression

            delx_w = xc(i) - xc(i-1)  ! Distance to west cell center
            dely_s = yc(j) - yc(j-1)  ! Distance to south cell center
            dely_n = yc(j+1) - yc(j)  ! Distance to north cell center

            ! Modify the coefficients for Neumann BC
            a_w(i,j) = ((0.5_dp * (kappa(i-1,j) + kappa(i,j))) * (dy(j)/delx_w))
            a_e(i,j) = 0.0_dp
            a_s(i,j) = ((0.5_dp * (kappa(i,j-1) + kappa(i,j))) * (dx(i)/dely_s))
            a_n(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i,j+1))) * (dx(i)/dely_n))

            a_p(i,j) = a_e(i,j) + a_w(i,j) + a_n(i,j) + a_s(i,j) - (s_p(i,j)*dx(i)*dy(j))

            qx(i,j) = q ! heat flux normal to the x direction
            qy(i,j) = 0.0_dp ! heat flux in y direction
            qmag(i,j) = sqrt(qx(i,j)**2 + qy(i,j)**2) ! magnitude of heat flux

        end do

    end if

    ! Top Boundary (j=jmax)
    j = jmax

    if (bctype_t .eq. 'dir') then

        do i = 1, imax

            varvals(1) = xc(i)
            varvals(2) = yc(j)
            varvals(3) = lx
            varvals(4) = ly
            varvals(5) = acos(-1.0_dp) ! pi
            
            
            T(i,j) = evalf(6, varvals) ! Evaluate the Dirichlet BC expression

            ! Modify the coefficients for Dirichlet BC
            a_w(i,j) = 0.0_dp
            a_e(i,j) = 0.0_dp
            a_s(i,j) = 0.0_dp
            a_n(i,j) = 0.0_dp
            a_p(i,j) = 1.0_dp

        end do

    else

        do i = 1, imax - 1

            varvals(1) = xc(i)
            varvals(2) = yc(j)
            varvals(3) = lx
            varvals(4) = ly
            varvals(5) = acos(-1.0_dp) ! pi

            q = evalf(6, varvals) ! Evaluate the Neumann BC expression

            delx_w = xc(i) - xc(i-1)  ! Distance to west cell center
            delx_e = xc(i+1) - xc(i)  ! Distance to east cell center
            dely_s = yc(j) - yc(j-1)  ! Distance to south cell center

            ! Modify the coefficients for Neumann BC
            a_w(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_w))
            a_e(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_e))
            a_s(i,j) = ((0.5_dp * (kappa(i,j-1) + kappa(i,j))) * (dx(i)/dely_s))
            a_n(i,j) = 0.0_dp

            a_p(i,j) = a_e(i,j) + a_w(i,j) + a_n(i,j) + a_s(i,j) - (s_p(i,j)*dx(i)*dy(j))

            qx(i,j) = 0.0_dp ! heat flux in x direction
            qy(i,j) = q ! heat flux normal to the x direction
            qmag(i,j) = sqrt(qx(i,j)**2 + qy(i,j)**2) ! magnitude of heat flux

        end do

    end if

    ! Left Boundary (i=1)
    i = 1

    if (bctype_l .eq. 'dir') then
        
        do j = 1, jmax

        varvals(1) = xc(i)
        varvals(2) = yc(j)
        varvals(3) = lx
        varvals(4) = ly
        varvals(5) = acos(-1.0_dp) ! pi

        
        T(i,j) = evalf(7, varvals) ! Evaluate the Dirichlet BC expression

        ! Modify the coefficients for Dirichlet BC
        a_w(i,j) = 0.0_dp
        a_e(i,j) = 0.0_dp
        a_s(i,j) = 0.0_dp
        a_n(i,j) = 0.0_dp
        
        end do

    else

        do j = 1, jmax - 1

            varvals(1) = xc(i)
            varvals(2) = yc(j)
            varvals(3) = lx
            varvals(4) = ly
            varvals(5) = acos(-1.0_dp) ! pi

            q = evalf(7, varvals) ! Evaluate the Neumann BC expression

            delx_e = xc(i+1) - xc(i)  ! Distance to east cell center
            dely_s = yc(j) - yc(j-1)  ! Distance to south cell center
            dely_n = yc(j+1) - yc(j)  ! Distance to north cell center

            ! Modify the coefficients for Neumann BC
            a_w(i,j) = 0.0_dp
            a_e(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i+1,j))) * (dy(j)/delx_e))
            a_s(i,j) = ((0.5_dp * (kappa(i,j-1) + kappa(i,j))) * (dx(i)/dely_s))
            a_n(i,j) = ((0.5_dp * (kappa(i,j) + kappa(i,j+1))) * (dx(i)/dely_n))

            a_p(i,j) = a_e(i,j) + a_w(i,j) + a_n(i,j) + a_s(i,j) - (s_p(i,j)*dx(i)*dy(j))

            qx(i,j) = q ! heat flux normal to the x direction
            qy(i,j) = 0.0_dp ! heat flux in y direction
            qmag(i,j) = sqrt(qx(i,j)**2 + qy(i,j)**2) ! magnitude of heat flux

        end do

    end if

    allocate(R(max_iter))   ! Residual array allocation
    R(:) = 0.0_dp  ! Initialize the residual array
    close(10)

    print *,'Boundary Conditions, Source Terms and Thermal Conductivity assigned successfully!'
    print *, ' '

end subroutine bc_s_k
