! Hello World in OpenACC

program hello_world
    use openacc
    implicit none

    print *, "Hello World from the host!"

    !$acc parallel
    print *, "Hello World from the device!"
    !$acc end parallel

end program hello_world