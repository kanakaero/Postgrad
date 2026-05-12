! Solve Steady State Heat Conduction using OpenMP
! delta_x = 0.01 and 0.001 for different values of nthreads = 1,2,4,8
! u(x) = 7 - x*tan(x) and x = -1 to 1
! compute du/dx using first, second and fourth order finite difference methods and 
! compare the results with the analytical solution.
! Boundary points stick to first order finite difference method

program HeatConduction
    use omp_lib
    implicit none
    integer :: i, n, nthreads, k
    double precision, allocatable :: x(:), u(:), du_dx(:)
    integer, dimension(2) :: sizes
    real(8) :: start_time, end_time
    real(8), allocatable :: times(:)
    character(100) :: argument
    real(8) :: dx
    real(8), allocatable :: du_exact(:)
    
    ! Define grid sizes
    sizes = [100, 1000]
    allocate(times(size(sizes)))

    call GET_COMMAND_ARGUMENT(1, argument)
    read(argument, *) nthreads

    print *, "Number of threads: ", nthreads

    do k = 1, size(sizes)
        n = sizes(k)
        print *, "Grid size: ", n

        ! Allocate arrays
        allocate(x(1:n))
        allocate(u(1:n))
        allocate(du_dx(1:n))
        
        ! Initialize x and u using the analytical solution
        do i = 1, n
            x(i) = -1.0d0 + (i-1)*2.0d0/(n-1) ! x from -1 to 1
            u(i) = 7.0d0 - x(i)*tan(x(i)) ! Analytical solution
        end do

        ! Record start time
        start_time = omp_get_wtime()

        ! Boundary points using first order finite difference method
        dx = 2.0d0/(n-1)

        ! Compute du/dx using finite difference methods in parallel
        !$OMP PARALLEL DO NUM_THREADS(nthreads)
        do i = 2, n-1
            du_dx(i) = (u(i+1) - u(i-1)) / (2.0d0 * dx) ! Central difference for interior points - second order accurate
        end do
        !$OMP END PARALLEL DO

        du_dx(1) = (u(2) - u(1)) / dx
        du_dx(n) = (u(n) - u(n-1)) / dx

        ! Record end time
        end_time = omp_get_wtime()
        
        ! Store execution time for this grid size
        times(k) = end_time - start_time
        print *, "Execution time: ", times(k), " seconds"

        ! write execution times to a file
        open(unit=10, file='execution_times_heat_conduction_2nd.txt', status='unknown', action='write', position='append')
        write(10,*) nthreads, sizes(k), times(k)
        close(10)

        du_exact = -tan(x) - x*(1.0d0 + tan(x)**2)

        print *, "Execution time (2nd order): ", times(k), " seconds"
        print *, "Error with second order finite difference: ", maxval(abs(du_dx(2:n-1) - du_exact(2:n-1))) ! Compare with analytical derivative

        ! Deallocate arrays
        deallocate(x)
        deallocate(u)
        deallocate(du_dx)
        deallocate(du_exact)
    end do

    do k = 1, size(sizes)
        n = sizes(k)
        print *, "Grid size: ", n

        ! Allocate arrays
        allocate(x(1:n))
        allocate(u(1:n))
        allocate(du_dx(1:n))
        
        ! Initialize x and u using the analytical solution
        do i = 1, n
            x(i) = -1.0d0 + (i-1)*2.0d0/(n-1) ! x from -1 to 1
            u(i) = 7.0d0 - x(i)*tan(x(i)) ! Analytical solution
        end do

        ! Record start time
        start_time = omp_get_wtime()

        ! Boundary points using first order finite difference method
        dx = 2.0d0/(n-1)

        ! Compute du/dx using finite difference methods in parallel
        !$OMP PARALLEL DO NUM_THREADS(nthreads)
        do i = 3, n-2
            du_dx(i) = (-u(i+2) + 8.0d0*u(i+1) - 8.0d0*u(i-1) + u(i-2)) / (12.0d0 * dx) ! Fourth order central difference for interior points
        end do
        !$OMP END PARALLEL DO

        du_dx(1) = (u(2) - u(1)) / dx
        du_dx(n) = (u(n) - u(n-1)) / dx

        ! Record end time
        end_time = omp_get_wtime()
        
        ! Store execution time for this grid size
        times(k) = end_time - start_time
        
        du_exact = -tan(x) - x*(1.0d0 + tan(x)**2)

        print *, "Execution time (4th order): ", times(k), " seconds"
        print *, "Error with fourth order finite difference: ", maxval(abs(du_dx(3:n-2) - du_exact(3:n-2))) ! Compare with analytical derivative

        ! write execution times to a file
        open(unit=10, file='execution_times_heat_conduction_4th.txt', status='unknown', action='write', position='append')
        write(10,*) nthreads, sizes(k), times(k)
        close(10)

        ! Deallocate arrays
        deallocate(x)
        deallocate(u)
        deallocate(du_dx)
        deallocate(du_exact)
    end do

end program HeatConduction