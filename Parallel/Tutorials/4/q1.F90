! Program to add two square matrices of NxN size using opnemp and test for different values of N and number of threads

program MatAdd
    use omp_lib
    implicit none
    integer :: i, j, nthreads, m, k
    double precision, allocatable :: A(:,:), B(:,:), C(:,:)
    integer, dimension(10) :: sizes
    real(8) :: start_time, end_time
    real(8), allocatable :: times(:)
    character(100) :: argument
    
    ! Define matrix size
    sizes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    allocate(times(size(sizes)))

    call GET_COMMAND_ARGUMENT(1, argument)
    read(argument, *) nthreads

    print *, "Number of threads: ", nthreads

    do k = 1, size(sizes)
        m = sizes(k)
        print *, "Matrix size: ", m, "x", m

        ! Allocate matrices
        allocate(A(1:m, 1:m))
        allocate(B(1:m, 1:m))
        allocate(C(1:m, 1:m))
        
        ! Initialize matrices A and B using Random numbers
        call random_number(A)
        call random_number(B)

        ! Record start time
        start_time = omp_get_wtime()

        ! Perform matrix addition C = A + B
        !$OMP PARALLEL DO NUM_THREADS(nthreads)
        do i = 1, m
            do j = 1, m
                C(i,j) = A(i,j) + B(i,j)
            end do
        end do
        !$OMP END PARALLEL DO

        ! Record end time
        end_time = omp_get_wtime()
        
        ! Store execution time for this matrix size
        times(k) = end_time - start_time
        print *, "Execution time: ", times(k), " seconds"

        ! Deallocate matrices
        deallocate(A)
        deallocate(B)
        deallocate(C)
    end do

    write(*,*) "Matrix Size (NxN) and Execution Time (seconds):"
    do k = 1, size(sizes)
        write(*,*) sizes(k), times(k)
    end do

    ! write execution times to a file
    open(unit=10, file='execution_times.txt', status='unknown', action='write', position='append')
    do k = 1, size(sizes)
        write(10,*) nthreads, sizes(k), times(k)
    end do
    close(10)

end program MatAdd