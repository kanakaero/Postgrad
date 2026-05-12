! Question 2 - Program to Multiply two matrices to perform matrix vector multiplication
! Kanak Agarwal - ID25M802

subroutine MatrixVectorGeneration(m, n, mat, vec)
    implicit none
    integer, intent(in) :: m, n
    double precision, dimension(1:m, 1:n), intent(out) :: mat
    double precision, dimension(1:n), intent(out) :: vec

    ! Initialize matrices with random values
    call random_seed()
    call random_number(mat)
    call random_number(vec)

end subroutine MatrixVectorGeneration

subroutine multiplyTwoMatrices(m, n, mat, vec, Product, sizes_index, times)
    implicit none
    integer, intent(in) :: m, n
    double precision, dimension(1:m, 1:n), intent(in) :: mat
    double precision, dimension(1:n), intent(in) :: vec
    double precision, dimension(1:m), intent(out) :: Product
    integer :: i, j, k
    double precision :: t1, t2
    integer, intent(in) :: sizes_index
    double precision, intent(inout) :: times(4)

    ! time at the beginning
    call cpu_time(t1)

    ! multiply the matrix and vector
    do i = 1, m
        Product(i) = 0.0
        do j = 1, n
            Product(i) = Product(i) + mat(i, j) * vec(j)
        end do
    end do
    
    ! time at the end 
    call cpu_time(t2)

    write(*, *) 'Multiplication took ', t2-t1, ' seconds'
    write(*, *) ' '

    times(sizes_index) = t2 - t1

end subroutine multiplyTwoMatrices

subroutine printMatrix(m, n, mat)
    implicit none
    integer, intent(in) :: m, n 
    double precision, dimension(1:m, 1:n), intent(in) :: mat
    integer :: i, j 
    
    print *, " "
    write(*, *) 'Printing Product:'
    
    do i = 1, m
        write(*, '(100f10.4)') (mat(i, j), j = 1, n)
    end do
  
end subroutine printMatrix

program MatrixVectorMultiplication
    implicit none
    integer :: m, n, sizes_index
    double precision, allocatable :: mat(:,:), vec(:), Product(:)
    integer, parameter :: sizes(4) = [64, 128, 256, 512]
    double precision :: times(4)

    do sizes_index = 1, size(sizes)
        m = sizes(sizes_index)
        n = sizes(sizes_index)

        print *, "Matrix-Vector Multiplication for size: ", m, " x ", n

        allocate(mat(1:m, 1:n))
        allocate(vec(1:n))
        allocate(Product(1:m))

        call MatrixVectorGeneration(m, n, mat, vec)

        call multiplyTwoMatrices(m, n, mat, vec, Product, sizes_index, times)

        deallocate(mat)
        deallocate(vec)
        deallocate(Product)
    end do

    print *, "Completed Matrix-Vector Multiplication for all sizes."
    print *, " "
    
    print *, "Sizes and corresponding times taken (in seconds):"
    do sizes_index = 1, size(sizes)
        write(*, '(I6, 2X, F10.6)') sizes(sizes_index), times(sizes_index)
    end do

end program MatrixVectorMultiplication