! Matrix Addition using OMP

program MatAdd
   use omp_lib
   implicit none
   integer :: i, j, m, nthreads
   double precision, allocatable :: A(:, :), B(:, :), C(:, :)
   character(100) :: argument

   ! Define matrix size
   m = 10

   call GET_COMMAND_ARGUMENT(1, argument)
   read (argument, *) nthreads

   ! Allocate matrices
   allocate (A(1:m, 1:m))
   allocate (B(1:m, 1:m))
   allocate (C(1:m, 1:m))

   ! Initialize matrices A and B using Random numbers
   call random_number(A)
   call random_number(B)

   ! Print the input matrices

   do i = 1, m
      do j = 1, m
         C(i, j) = 0.0
      end do
   end do

   print *, "Matrix A"
   do i = 1, m
      write (*, '(10F10.2)') (A(i, j), j=1, m)
   end do

   print *, " "
   print *, "Matrix B"
   do i = 1, m
      write (*, '(10F10.2)') (B(i, j), j=1, m)
   end do

   ! Perform matrix addition C = A + B
   !$OMP PARALLEL DO NUM_THREADS(nthreads)
   do i = 1, m
      do j = 1, m
         C(i, j) = A(i, j) + B(i, j)
      end do
   end do
   !$OMP END PARALLEL DO

   print *, " "
   print *, "Matrix C = A + B"
   do i = 1, m
      write (*, '(10F10.2)') (C(i, j), j=1, m)
   end do

   ! Deallocate matrices
   deallocate (A)
   deallocate (B)
   deallocate (C)

end program MatAdd
