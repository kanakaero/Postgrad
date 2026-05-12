! OpenMPI program to demonstrate  MPI Gatherv, MPI Scatterv, MPI Alltoall, and the timing functions MPI Wtime, MPI Wtick.

program demo
   use mpi
   implicit none
   integer :: rank, num_procs, i, ierror
   real :: start, end
   integer :: a(8), recvbuf(2), sendcounts(4), displs(4)
   integer :: sendrow(4), recvcol(4)
   integer :: j

   call mpi_init(ierror)
   call mpi_comm_rank(mpi_comm_world, rank, ierror)
   call mpi_comm_size(mpi_comm_world, num_procs, ierror)

   start = mpi_wtime()

   a = 0
   recvbuf = 0
   displs = 0
   sendcounts = 0

   if (rank == 0) then
      a = [1, 2, 3, 4, 5, 6, 7, 8]
      sendcounts = [2, 2, 2, 2]
      displs = [0, 2, 4, 6]
   end if

   call mpi_scatterv(a, sendcounts, displs, mpi_integer, recvbuf, 2, mpi_integer, 0, mpi_comm_world, ierror)

   print *, "Process ", rank, " received: ", recvbuf

   call mpi_barrier(mpi_comm_world, ierror)

   if (rank == 0) then
      a = 0
      print *, ' '
   end if

   call mpi_gatherv(recvbuf, 2, mpi_integer, a, sendcounts, displs, mpi_integer, 0, mpi_comm_world, ierror)

   if (rank == 0) then
      print *, "Gathered array: ", a
      print *, ' '
   end if

   call mpi_barrier(mpi_comm_world, ierror)

   ! each rank creates its own row
   do j = 1, 4
      sendrow(j) = rank*4 + j
   end do

   print *, "Rank", rank, " row:", sendrow

   call mpi_barrier(mpi_comm_world, ierror)

   call mpi_alltoall(sendrow, 1, mpi_integer, &
                     recvcol, 1, mpi_integer, &
                     mpi_comm_world, ierror)

   print *, "Rank", rank, " column after Alltoall:", recvcol

   end = mpi_wtime()

   if (rank == 0) then
      print *, ' '
      print *, "Time taken: ", end - start, " seconds"
   end if

   call mpi_finalize(ierror)

end program demo
