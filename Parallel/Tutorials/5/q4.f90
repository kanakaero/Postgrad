! MPI Program for trapezoidal rule to calculate the integral of a function f(x) = x^2 from 0 to 1

program trapezoidal
    implicit none
    include "mpif.h"
    integer :: i, myid, size, mpierror, tag, status(MPI_STATUS_SIZE)
    real :: a, b, h, local_sum, total_sum, local_sum1
    integer :: n, local_n, istart, iend

    call MPI_INIT(mpierror)       

    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, mpierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierror)

    a = 0.0 ! Start of the interval
    b = 1.0 ! End of the interval
    n = 1000 ! Number of subintervals
    h = (b - a) / n ! Width of each subinterval

    n = 1000
    local_n = n / size

    istart = myid * local_n
    iend   = istart + local_n

    ! Each process calculates its local sum
    local_sum = 0.0

    do i = istart, iend - 1
        local_sum = local_sum + 0.5 * h * (f(a + i*h) + f(a + (i+1)*h))
    end do


    !print *, 'Process', myid, 'local sum:', local_sum

    !print *, size, myid

    if (myid /= 0) then 
        call MPI_SEND(local_sum, 1, MPI_REAL, 0, myid, MPI_COMM_WORLD, mpierror)
        print *, 'Process', myid, 'sent local sum:', local_sum
    else
        print *, 'In the host'
        total_sum = local_sum ! Start with the local sum from process 0
        do i = 1, size-1
            print *, 'Waiting to receive local sum from process', i
            call MPI_RECV(local_sum1, 1, MPI_REAL, i, i, MPI_COMM_WORLD, status, mpierror)
            total_sum = total_sum + local_sum1 ! Accumulate the total sum
            print *, 'Received local sum from process', i, ':', local_sum1
        end do
        print *, ' '
        write(*, *) 'The integral of f(x) from', a, 'to', b, 'is approximately:', total_sum     
    end if


    call MPI_FINALIZE(mpierror)

contains

    real function f(x)
        real :: x
        f = x**2 ! Define the function f(x) = x^2
    end function f
    
end program trapezoidal