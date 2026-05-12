! Question 1 - Program to prove that A + A' is symmetric
! Kanak Agarwal - ID25M802

program symmetric_matrix
    implicit none
    integer :: n
    real, allocatable :: A(:,:), AT(:,:), S(:,:)
    integer :: i, j

    write(*,*) "Enter the number of elements:"
    read(*,*) n
    print *, " "

    allocate(A(n,n), AT(n,n), S(n,n))

    print *, "Generating a random matrix A of size", n, "x", n

    ! Initialize matrix A
    call random_seed()
    call random_number(A)

    print *, " "
    print *, "Matrix A:"
    do i = 1, n
        print *, A(i,:)
    end do
    print *, " "

    AT = transpose(A) ! Compute Transpose of A

    ! Compute S = A + A'
    S = A + AT

    ! Print the resulting matrix S
    print *, "Matrix S (A + A'):"
    do i = 1, n
        print *, S(i,:)
    end do
    print *, " "

    ! Check if S is symmetric
    print *, "Checking if S is symmetric:"
    do i = 1, n
        do j = 1, n
        if (S(i,j) /= S(j,i)) then
            print *, "Matrix S is not symmetric."
            stop
        end if
        end do
    end do
    print *, "Matrix S is symmetric."

end program symmetric_matrix