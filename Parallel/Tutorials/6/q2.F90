! Write an MPI Program to matrix-vector multiplication using block decomposition and verify the result with a serial program. Use collective communication.

program mpi_matvec_block
   use mpi
   implicit none

   integer :: ierr, rank, size
   integer :: N, local_n, i, j
   real(8), allocatable :: A(:, :), x(:), y(:)
   real(8), allocatable :: A_local(:, :), y_local(:)
   real(8), allocatable :: y_serial(:)
   real(8) :: error

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

   ! -----------------------------
   ! Problem size (must be divisible by size)
   ! -----------------------------
   N = 8
   if (mod(N, size) /= 0) then
      if (rank == 0) print *, "N must be divisible by number of processes"
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
   end if

   local_n = N/size

   ! -----------------------------
   ! Allocate arrays
   ! -----------------------------
   allocate (A_local(local_n, N))
   allocate (y_local(local_n))
   allocate (x(N))

   if (rank == 0) then
      allocate (A(N, N))
      allocate (y(N))
      allocate (y_serial(N))

      ! Initialize matrix and vector
      do i = 1, N
         x(i) = 1.0d0
         do j = 1, N
            A(i, j) = real(i + j, 8)
         end do
      end do
   end if

   ! Broadcast vector x
   call MPI_Bcast(x, N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   ! Scatter matrix rows
   call MPI_Scatter(A, local_n*N, MPI_DOUBLE_PRECISION, &
                    A_local, local_n*N, MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)

   ! Local matrix-vector multiply
   y_local = 0.0d0
   do i = 1, local_n
      do j = 1, N
         y_local(i) = y_local(i) + A_local(i, j)*x(j)
      end do
   end do

   ! Gather result vector
   call MPI_Gather(y_local, local_n, MPI_DOUBLE_PRECISION, &
                   y, local_n, MPI_DOUBLE_PRECISION, &
                   0, MPI_COMM_WORLD, ierr)

   ! Serial
   if (rank == 0) then
      y_serial = 0.0d0
      do i = 1, N
         do j = 1, N
            y_serial(i) = y_serial(i) + A(i, j)*x(j)
         end do
      end do

      error = maxval(abs(y - y_serial))

      print *, "Parallel result y:"
      print *, y
      print *, "Serial result y_serial:"
      print *, y_serial
      print *, "Max error =", error
   end if

   ! -----------------------------
   ! Cleanup
   ! -----------------------------
   if (rank == 0) then
      deallocate (A, y, y_serial)
   end if
   deallocate (A_local, y_local, x)

   call MPI_Finalize(ierr)
end program mpi_matvec_block
