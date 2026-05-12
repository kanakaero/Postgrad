! Hello world program using OpenMP
program main

    use OMP_LIB
    implicit none
    integer :: i

    !$OMP PARALLEL
    i = OMP_GET_NUM_THREADS()-1
    !$OMP END PARALLEL
    
    write(*, *) 'About to enter the parallel world...'

    !  write(*, *) 'Number of threads = ', OMP_GET_NUM_THREADS()
    
    !$OMP PARALLEL

    ! Print from ranks in ascending order - deterministic output

    do while (i <= OMP_GET_NUM_THREADS())
        if (i==OMP_GET_THREAD_NUM()) then
            write(*, *) 'Hello world from thread number', OMP_GET_THREAD_NUM()
            i = i - 1
            exit
        end if
    end do

    !  write(*, *) 'Hello world from thread number', OMP_GET_THREAD_NUM()
    !  write(*, *) 'Number of threads = ', OMP_GET_NUM_THREADS()

    !$OMP END PARALLEL

    write(*, *) 'Entered the serial world...'

    !write(*, *) 'Number of threads = ', OMP_GET_NUM_THREADS() 

end program main