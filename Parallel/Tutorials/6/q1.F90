program mpi_collectives_demo
   use mpi
   implicit none

   integer :: ierr, rank, size
   integer :: i
   integer :: send_val, recv_sum, all_sum
   integer, allocatable :: send_array(:), recv_array(:)
   integer, allocatable :: gather_array(:), allgather_array(:)

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

   ! MPI_Reduce
   send_val = rank + 1

   call MPI_Reduce(send_val, recv_sum, 1, MPI_INTEGER, MPI_SUM, &
                   0, MPI_COMM_WORLD, ierr)

   if (rank == 0) then
      print *, "MPI_Reduce (SUM) result =", recv_sum
   end if

   ! MPI_Allreduce
   call MPI_Allreduce(send_val, all_sum, 1, MPI_INTEGER, MPI_SUM, &
                      MPI_COMM_WORLD, ierr)

   print *, "Rank", rank, ": MPI_Allreduce result =", all_sum

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   if (rank == 0) then
      print *, "----------------------------------------"
   end if
   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   ! MPI_Scatter

   if (rank == 0) then
      allocate (send_array(size))
      do i = 1, size
         send_array(i) = i*10
      end do
   end if

   call MPI_Scatter(send_array, 1, MPI_INTEGER, &
                    send_val, 1, MPI_INTEGER, &
                    0, MPI_COMM_WORLD, ierr)

   print *, "Rank", rank, ": MPI_Scatter received =", send_val

   if (rank == 0) deallocate (send_array)

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   ! MPI_Gather
   if (rank == 0) then
      allocate (gather_array(size))
   end if

   call MPI_Gather(rank, 1, MPI_INTEGER, &
                   gather_array, 1, MPI_INTEGER, &
                   0, MPI_COMM_WORLD, ierr)

   if (rank == 0) then
      print *, "MPI_Gather result:"
      print *, gather_array
      deallocate (gather_array)
   end if

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   ! MPI_Allgather
   allocate (allgather_array(size))

   call MPI_Allgather(rank, 1, MPI_INTEGER, &
                      allgather_array, 1, MPI_INTEGER, &
                      MPI_COMM_WORLD, ierr)

   print *, "Rank", rank, ": MPI_Allgather =", allgather_array

   deallocate (allgather_array)

   call MPI_Finalize(ierr)
end program mpi_collectives_demo
