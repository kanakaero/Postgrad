program mesh_generator
    implicit none

    ! Domain parameters
    real :: Lx, Ly
    integer :: nx, ny
    real :: stretch_factor_X, stretch_factor_Y

    ! Arrays
    real, allocatable :: dx(:), dy(:)
    real, allocatable :: grid_lines_X(:), grid_lines_Y(:)
    real, allocatable :: centre_X(:), centre_Y(:)

    ! Loop counters
    integer :: i, j
    real :: a_X, a_Y, x, y

    ! -------------------------
    ! Initialize parameters
    ! -------------------------
    Lx = 1.0
    Ly = 1.0
    nx = 10
    ny = 10

    stretch_factor_X = 1.5
    stretch_factor_Y = 1.5

    ! -------------------------
    ! Compute first cell size
    ! -------------------------
    if (abs(stretch_factor_X - 1.0) < 1e-12) then
        a_X = Lx / nx
    else if (stretch_factor_X < 1.0) then
        a_X = Lx * (1.0 - stretch_factor_X) / (1.0 - stretch_factor_X**nx)
    else
        a_X = Lx * (stretch_factor_X - 1.0) / (stretch_factor_X**nx - 1.0)
    end if

    if (abs(stretch_factor_Y - 1.0) < 1e-12) then
        a_Y = Ly / ny
    else if (stretch_factor_Y < 1.0) then
        a_Y = Ly * (1.0 - stretch_factor_Y) / (1.0 - stretch_factor_Y**ny)
    else
        a_Y = Ly * (stretch_factor_Y - 1.0) / (stretch_factor_Y**ny - 1.0)
    end if

    ! -------------------------
    ! Allocate arrays
    ! -------------------------
    allocate(dx(nx), dy(ny))
    allocate(grid_lines_X(nx+1), grid_lines_Y(ny+1))
    allocate(centre_X(nx), centre_Y(ny))

    ! -------------------------
    ! Fill dx, dy arrays
    ! -------------------------
    dx(1) = a_X
    do i = 2, nx
        dx(i) = dx(i-1) * stretch_factor_X
    end do

    dy(1) = a_Y
    do i = 2, ny
        dy(i) = dy(i-1) * stretch_factor_Y
    end do

    ! -------------------------
    ! Compute grid lines
    ! -------------------------
    grid_lines_X(1) = 0.0
    do i = 1, nx
        grid_lines_X(i+1) = grid_lines_X(i) + dx(i)
        centre_X(i) = 0.5 * (grid_lines_X(i) + grid_lines_X(i+1))
    end do

    grid_lines_Y(1) = 0.0
    do i = 1, ny
        grid_lines_Y(i+1) = grid_lines_Y(i) + dy(i)
        centre_Y(i) = 0.5 * (grid_lines_Y(i) + grid_lines_Y(i+1))
    end do

    ! ----------------------------
    ! Write Mesh to Output File
    ! ----------------------------

    ! Vertical lines
    open(unit=20, file="vertical_lines.dat", status="replace")
    do i = 1, nx+1
        do j = 1, ny+1
            write(20,'(2f15.6)') grid_lines_X(i), grid_lines_Y(j)
        end do
        write(20,*)   ! blank line to separate datasets
    end do
    close(20)

    ! Horizontal lines
    open(unit=21, file="horizontal_lines.dat", status="replace")
    do j = 1, ny+1
        do i = 1, nx+1
            write(21,'(2f15.6)') grid_lines_X(i), grid_lines_Y(j)
        end do
        write(21,*)   ! blank line to separate datasets
    end do
    close(21)

    ! Cell centres
    open(unit=22, file="centres.dat", status="replace")
    do j = 1, ny
        do i = 1, nx
            write(22,'(2f15.6)') centre_X(i), centre_Y(j)
        end do
    end do
    close(22)

    print *, "Mesh written:"
    print *, " vertical_lines.dat, horizontal_lines.dat, centres.dat"

end program mesh_generator

