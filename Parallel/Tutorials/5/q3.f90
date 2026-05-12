! MPI Program to Receive single real floating point from all processes and print the sum of all the received values in process 0

program msum
  implicit none
  include "mpif.h"
  integer :: i, myid, size, mpierror, tag, status(MPI_STATUS_SIZE)
  real :: value_send, value_recv, sum

  call MPI_INIT(mpierror)       

  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, mpierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierror)

  tag = 100
  value_send = real(myid) ! Each process sends its rank as a real number

  print *, 'Process', myid, 'sending value:', value_send

  if (myid /= 0) then 
     call MPI_SEND(value_send, 1, MPI_REAL, 0, tag, MPI_COMM_WORLD, mpierror)
  else
     sum = value_send ! Start with the value from process 0
     do i = 1, size-1
        call MPI_RECV(value_recv, 1, MPI_REAL, i, tag, MPI_COMM_WORLD, status, mpierror)
        sum = sum + value_recv ! Accumulate the sum
     end do
     print *, ' '
     write(*, *) 'The sum of all values is:', sum     
  end if

  call MPI_FINALIZE(mpierror)
  
end program msum