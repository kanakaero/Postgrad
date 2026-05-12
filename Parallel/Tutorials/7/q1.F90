program derivative
use mpi
implicit none

integer :: rank,num_procs,ierr
integer :: i,n_global,n_local
integer :: left,right
real :: x0,xn,h,x
real,allocatable :: u(:)
real :: du_forward,du_backward,du_central,du_central2,analytical
character(len=30) :: filename

call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,num_procs,ierr)

x0 = 0.0
xn = 3.0
h  = 0.01

n_global = int((xn-x0)/h)
n_local  = n_global/num_procs

allocate(u(-1:n_local+2))

! Fill local data
do i=1,n_local
    x = x0 + (rank*n_local + i)*h
    u(i) = f(x)
end do


! Neighbor ranks
left  = rank-1
right = rank+1

if(rank==0) left  = MPI_PROC_NULL
if(rank==num_procs-1) right = MPI_PROC_NULL


! Exchange TWO halo layers
call MPI_Sendrecv(u(1),2,MPI_REAL,left,0, &
                  u(-1),2,MPI_REAL,left,0, &
                  MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

call MPI_Sendrecv(u(n_local-1),2,MPI_REAL,right,0, &
                  u(n_local+1),2,MPI_REAL,right,0, &
                  MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)


! Physical boundaries
if(rank==0) then
    u(0)  = f(x0-h)
    u(-1) = f(x0-2*h)
endif

if(rank==num_procs-1) then
    u(n_local+1) = f(xn+h)
    u(n_local+2) = f(xn+2*h)
endif


! Output file
write(filename,'("derivative_rank",I0,".txt")') rank
open(10,file=filename,status='replace')

do i=1,n_local

    x = x0 + (rank*n_local + i)*h

    ! Forward
    if(i < n_local) then
        du_forward = (u(i+1)-u(i))/h
    else
        du_forward = 0.0
    endif

    ! Backward
    if(i > 1) then
        du_backward = (u(i)-u(i-1))/h
    else
        du_backward = 0.0
    endif

    ! Central
    if(i>1 .and. i<n_local) then
        du_central = (u(i+1)-u(i-1))/(2*h)
    else
        du_central = 0.0
    endif

    ! 4th order central
    if(i>=3 .and. i<=n_local-2) then
        du_central2 = (-u(i+2)+8*u(i+1)-8*u(i-1)+u(i-2))/(12*h)
    else
        cycle
    endif

    analytical = 3*x**2 - 5*cos(5*x)

    write(10,'(6ES16.8)') x,analytical,du_forward,du_backward,du_central,du_central2

end do

close(10)

call MPI_Finalize(ierr)

contains

function f(x) result(res)
real :: x,res
res = x**3 - sin(5*x)
end function

end program derivative