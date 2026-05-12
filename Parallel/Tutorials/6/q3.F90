! Write an MPI program that uses a derived data type to broadcasts these values to other processes in one go.

program derived_type
   use mpi
   implicit none

   integer :: n
   real :: a(3), b(2)
   integer :: rank, i
   integer :: ierr
   integer(kind=MPI_ADDRESS_KIND) :: a_addr, b_addr, n_addr, t_addr
   integer :: count
   integer :: array_of_blocks(3), array_of_types(3)
   integer :: new_type
   integer(kind=MPI_ADDRESS_KIND) :: array_of_displacements(3)

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

   if (rank == 0) then
      a = [12.67, 56.45, 25.32]
      b = [56.79, 98.26]
      n = 10
   end if

   if (rank == 0) then
      print *, "-----------------------------"
      print *, "Process ", rank, ": a = ", a
      print *, "Process ", rank, ": b = ", b
      print *, "Process ", rank, ": n = ", n
      print *, "-----------------------------"
   end if

   ! Create the derived data type
   array_of_blocks = [size(a), size(b), 1]
   count=3
   call MPI_Get_address(a(1), a_addr, ierr)
   call MPI_Get_address(b(1), b_addr, ierr)
   call MPI_Get_address(n, n_addr, ierr)
   t_addr = 0
   array_of_displacements = [t_addr, b_addr - a_addr, n_addr - a_addr]
   array_of_types = [MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_INTEGER]

   call MPI_Type_create_struct(count, array_of_blocks, array_of_displacements, array_of_types, new_type, ierr)

   call MPI_Type_commit(new_type, ierr)

   call MPI_Bcast(a, 1, new_type, 0, MPI_COMM_WORLD, ierr)

   if (rank /= 0) then
      print *, "-----------------------------"
      print *, "Process ", rank, ": a = ", a
      print *, "Process ", rank, ": b = ", b
      print *, "Process ", rank, ": n = ", n
      print *, "-----------------------------"
   end if

   call MPI_Type_free(new_type, ierr)
   call MPI_Finalize(ierr)

end program derived_type