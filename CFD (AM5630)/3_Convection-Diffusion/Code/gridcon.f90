! ----------------------------------------------------------------------
! This subroutine initializes the velocity field and the grid from the
! provided input files.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 12th November 2025
! ----------------------------------------------------------------------

subroutine gridcon()

    use var_dec
    implicit none

    integer :: ios, n ! iostat variable for error handling and temporary variable 
    real(dp) :: U_inlet, U_outlet, H

    U_inlet = 1.0_dp
    U_outlet = 1.0_dp
    H = 2.0_dp

    T_inlet = 20.0_dp
    T_outlet = 50.0_dp

    ! deallocate arrays
    if (allocated(dx)) deallocate(dx) 
    if (allocated(dy)) deallocate(dy) 
    if (allocated(xc)) deallocate(xc) 
    if (allocated(yc)) deallocate(yc)

    ! First pass - determine number of grid points in x
    open(unit=10, file='data/xc.dat', status='old', action='read', iostat=ios)
    n = 0
    do
        read(10, *, iostat=ios)
        if (ios /= 0) exit
        n = n + 1
    end do
    close(10)

    nx = n-1
    print *, "Number of divisions in the x-direction: ", nx

    ! Second pass - allocate to the array
    allocate(xg(nx+1))

    open(unit=10, file='data/xc.dat', status='old', action='read', iostat=ios)
    do i = 1, nx+1
        read(10, *, iostat=ios) xg(i)
    end do

    close(10)

    ! First pass - determine number of grid points in y
    open(unit=10, file='data/yc.dat', status='old', action='read', iostat=ios)
    n = 0
    do
        read(10, *, iostat=ios)
        if (ios /= 0) exit
        n = n + 1
    end do
    close(10)

    ny = n-1
    print *, "Number of divisions in the y-direction: ", ny

    ! Second pass - allocate to the array
    allocate(yg(ny+1))

    open(unit=10, file='data/yc.dat', status='old', action='read', iostat=ios)
    do i = 1, ny+1
        read(10, *, iostat=ios) yg(i)
    end do

    close(10)

    imax = nx + 2
    jmax = ny + 2

    ! Cell Centres - Cell nodes
    allocate(xc(imax))
    allocate(yc(jmax))

    do i = 2, nx+1
        xc(i) = 0.5d0 * (xg(i) + xg(i-1))
    end do
    xc(1) = xg(1)
    xc(imax) = xg(nx+1)

    do i = 2, ny+1
        yc(i) = 0.5d0 * (yg(i) + yg(i-1))
    end do
    yc(1) = yg(1)
    yc(jmax) = yg(ny+1)

    ! Calculate grid spacing arrays
    allocate(dx(imax-1))
    allocate(dy(jmax-1))

    do i = 1, imax-1
        dx(i) = xc(i+1) - xc(i)
    end do

    do i = 1, jmax-1
        dy(i) = yc(i+1) - yc(i)
    end do

    ! Vertical lines
    open(unit=20, file="vlines.dat", status="replace")
    do i = 1, imax+1
        do j = 1, jmax+1
            write(20,'(2f15.6)') xg(i), yg(j)
        end do
        write(20,*)   ! blank line to separate datasets
    end do
    close(20)

    ! Horizontal lines
    open(unit=21, file="hlines.dat", status="replace")
    do j = 1, jmax+1
        do i = 1, imax+1
            write(21,'(2f15.6)') xg(i), yg(j)
        end do
        write(21,*)   ! blank line to separate datasets
    end do
    close(21)

    ! Cell centres
    open(unit=22, file="ccent.dat", status="replace")
    do j = 1, jmax
        do i = 1, imax
            write(22,'(2f15.6)') xc(i), yc(j)
        end do
    end do
    close(22)

    print *, ' '
    print *, "Meshing Done - Mesh files vlines.dat, hlines.dat, ccent.dat written!"
    print *,' '

    ! Allocate velocity and flattened arrays
    allocate(u(imax, jmax))
    allocate(v(imax, jmax))

    allocate(u_flat((imax)*(jmax)))
    allocate(v_flat((imax)*(jmax)))

    ! Read velocity field from file
    open(10, file='data/u.dat', status='old', action='read', iostat=ios)
    do i = 1, (imax)*(jmax)
        read(10, *, iostat=ios) u_flat(i)
    end do
    close(10)

    open(10, file='data/v.dat', status='old', action='read', iostat=ios)
    do i = 1, (imax)*(jmax)
        read(10, *, iostat=ios) v_flat(i)
    end do
    close(10)

    ! Reshape flattened arrays to 2D arrays
    do j = 1, jmax
        do i = 1, imax
            u(i,j) = u_flat(i + (j-1)*imax)
            v(i,j) = v_flat(i + (j-1)*imax)
        end do
    end do

    ! Velocity Magnitude
    allocate(U_mag(imax, jmax))
    do j = 1, jmax
        do i = 1, imax
            U_mag(i, j) = sqrt(u(i, j)**2 + v(i, j)**2)
        end do
    end do

    do i = 1, imax
        v(i, jmax) = 0.0_dp    ! enforce impermeable north wall
    end do

    do i = 1, imax
        v(i, 1) = 0.0_dp       ! enforce impermeable south wall
    end do

    do i = 1, imax
        u(i, 1) = 0.0_dp       ! enforce no-slip at bottom wall
    end do

    do i = 1, imax
        u(i, jmax) = 0.0_dp    ! enforce no-slip at top wall
    end do

    do j = 1, jmax
        if (yc(j) >= (1.0_dp - 0.068_dp)*H) then
            u(1, j) = U_inlet      ! inlet velocity at left wall
        else
            u(1, j) = 0.0_dp
        end if
    end do

    do j = 1, jmax
        if (yc(j) <= 0.068_dp*H) then
            u(imax, j) = U_outlet   ! outlet velocity at right wall
        else
            u(imax, j) = 0.0_dp
        end if
    end do

    do j = 1, jmax
        v(1, j) = 0.0_dp       ! enforce no-slip at left wall
    end do

    do j = 1, jmax
        v(imax, j) = 0.0_dp    ! enforce no-slip at right wall
    end do

    open(unit=30, file='init.tec', status='replace')

    ! Write Tecplot header
    write(30, '(A)') 'TITLE = "Convection Diffusion Initial Conditions"'
    write(30, '(A)') 'VARIABLES = "X", "Y", "u", "v", "V"'
    write(30, '(A,I0,A,I0,A)') 'ZONE T="Zone 1", I=', imax, ', J=', jmax, ', DATAPACKING=POINT, ZONETYPE=ORDERED'

    ! Write data in point format: each line is one point (X, Y, T)
    do j = 1, jmax
        do i = 1, imax
            if (abs(u(i,j)) < 1.0d-100) u(i,j) = 0.0d0
            if (abs(v(i,j)) < 1.0d-100) v(i,j) = 0.0d0
            if (abs(V(i,j)) < 1.0d-100) V(i,j) = 0.0d0

            write(30,'(ES15.6,1X,ES15.6,1X,ES15.6,1X,ES15.6,1X,ES15.6,1X)') &
            xc(i), yc(j), u(i,j), v(i,j), U_mag(i,j)
        end do
    end do

    close(30)

    print *, 'Initial Conditions written to init.tec'
    print *, ' '

    print *, u(1,1), u(2,1), u(3,1)
    print *, u(1,2), u(2,2), u(3,2)

    print *, ' '

    print *, u_flat(1), u_flat(2), u_flat(3), u_flat(4), u_flat(5), u_flat(6), u_flat(7), u_flat(8), u_flat(9), u_flat(10)
    
end subroutine gridcon