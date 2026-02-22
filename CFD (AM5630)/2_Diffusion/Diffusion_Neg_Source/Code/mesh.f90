! ----------------------------------------------------------------------
! This subroutine is a general purpose routine to generate a mesh. It 
! can be used to generate both skewed and uniform 2D meshes.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 28th September 2025
! ----------------------------------------------------------------------

subroutine mesh_gen()
    use var_dec

        character(len=256) :: line  ! contains each line read from the input file

    integer :: ios ! iostat variable for error handling

    open(unit=10, file='bc_s_k.txt', status='old', action='read')

    print *, 'Reading input file bc_s_k.txt'

    do 
        read(10, '(A)', iostat=ios) line
        if (ios /= 0) exit
        
        line = adjustl(line)  ! remove leading spaces
        if (line(1:1) == '!') cycle  ! skip comment lines
        if (index(line, '!') /= 0) line = line(:index(line,'!') - 1)  ! remove inline comments
        if (trim(line) == '') cycle  ! skip blank lines

        ! Domain Parameters
        if (index(line,'lx') /= 0) read(line(4:), *) lx
        if (index(line,'ly') /= 0) read(line(4:), *) ly
        if (index(line,'nx') /= 0) read(line(4:), *) nx
        if (index(line,'ny') /= 0) read(line(4:), *) ny
        if (index(line,'beta_x') /= 0) read(line(8:), *) beta_x
        if (index(line,'beta_y') /= 0) read(line(8:), *) beta_y

    end do

    print *, 'Domain Parameters Read!'


    imax = int(nx+2)
    jmax = nint(ny+2)


    ! deallocate arrays
    if (allocated(dx)) deallocate(dx) 
    if (allocated(dy)) deallocate(dy) 
    if (allocated(xc)) deallocate(xc) 
    if (allocated(yc)) deallocate(yc) 
    if (allocated(xg)) deallocate(xg)
    if (allocated(yg)) deallocate(yg)

    if (allocated(a_p)) deallocate(a_p) 
    if (allocated(a_e)) deallocate(a_e) 
    if (allocated(a_w)) deallocate(a_w) 
    if (allocated(a_n)) deallocate(a_n) 
    if (allocated(a_s)) deallocate(a_s) 

    ! allocate arrays    
    allocate(dx(imax), dy(jmax)) ! Grid spacing arrays
    allocate(xc(imax), yc(jmax)) ! Cell centre arrays
    allocate(xg(imax+1), yg(jmax+1))     ! Grid line arrays

    allocate(a_p(imax, jmax), a_e(imax, jmax), a_w(imax, jmax), a_n(imax, jmax), a_s(imax, jmax)) ! Coefficient arrays

! --------------------- First Cell Length ---------------------------------

    if (beta_x .eq. 1.0_dp) then
        a_x = lx/nx ! Uniform mesh spacing in x direction
    elseif (beta_x .lt. 1.0_dp) then
        a_x = lx * (1 - beta_x)/(1 - beta_x**nx); ! First cell length in x direction for skewed mesh
    else
        a_x = lx * (beta_x - 1)/(beta_x**nx - 1); ! First cell length in x direction for skewed mesh
    end if

    if (beta_y .eq. 1.0_dp) then
        a_y = ly/ny ! Uniform mesh spacing in y direction
    elseif (beta_y .lt. 1.0_dp) then
        a_y = ly * (1 - beta_y)/(1 - beta_y**ny); ! First cell length in y direction for skewed mesh
    else
        a_y = ly * (beta_y - 1)/(beta_y**ny - 1); ! First cell length in y direction for skewed mesh
    end if

! --------------------------------------------------------------------------

! ----------------- Array of Cell Distances --------------------------------

    if (beta_x .le. 0.0_dp .or. beta_y .le. 0.0_dp) then

        dx(imax) = a_x     ! Ghost Cell
        dy(jmax) = a_y     ! Ghost Cell

        dx(imax - 1) = dx(imax)   ! Last Cell
        dy(jmax - 1) = dy(jmax)   ! Last Cell

        do i = imax-2, 1, -1
            dx(i) = dx(i+1) * abs(beta_x)
        end do

        dx(1) = dx(2)   ! First Cell

        do j = jmax-2, 1, -1
            dy(j) = dy(j+1) * abs(beta_y)
        end do

        dy(1) = dy(2)   ! First Cell

    else
        dx(1) = a_x     ! Ghost Cell
        dy(1) = a_y     ! Ghost Cell

        dx(2) = dx(1)   ! First Cell
        dy(2) = dy(1)   ! First Cell

        do i = 3, imax-1    ! Loop from the second element/cell
            dx(i) = dx(i-1) * beta_x      ! Array of dx's
        end do

        dx(imax) = dx(imax - 1)

        do j = 3, jmax-1    ! Loop from the second element/cell
            dy(j) = dy(j-1) * beta_y      ! Array of dy's
        end do

        dy(jmax) = dy(jmax - 1)

    end if

    

! --------------------------------------------------------------------------

! --------------- Array of Grid Lines & Cell Centers -----------------------

    xg(1) = 0.0_dp  ! First grid line in x direction
    yg(1) = 0.0_dp  ! First grid line in y direction

    x = 0.0_dp
    do i = 1, imax
        xg(i+1) = xg(i) + dx(i)     ! Array of x grid lines
        xc(i) = 0.5_dp * (xg(i) + xg(i+1)) ! Array of x cell centres
    end do

    y = 0.0_dp
    do j = 1, jmax
        yg(j+1) = yg(j) + dy(j)   ! Array of y grid lines
        yc(j) = 0.5_dp * (yg(j) + yg(j+1)) ! Array of y cell centres
    end do

! ---------------------------------------------------------------------------

! ---------------------- Write Mesh to Output File --------------------------

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

    print *, "Meshing Done - Mesh files vlines.dat, hlines.dat, ccent.dat written!"
    print *,' '

! ---------------------------------------------------------------------------

end subroutine mesh_gen