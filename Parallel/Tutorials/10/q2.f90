program conv1d_test
  implicit none
  integer, parameter :: K = 5
  integer :: N_list(3) = (/1000, 2000, 3000/)
  integer :: N, i, j, m, left, right, t
  integer, allocatable :: A(:), A_old(:), A_new(:)
  integer :: F(K)
  real(8) :: t1, t2

  F = (/3, 4, 5, 4, 3/)
  m = (K - 1) / 2

  do t = 1, 3
     N = N_list(t)

     allocate(A(N), A_old(N), A_new(N))

     ! Initialize
     do i = 1, N
        A(i) = i
     end do
     A_old = A

     !$acc data copyin(A_old, F) copyout(A_new)
     !$acc parallel loop
     do i = 1, N
        A_new(i) = 0

        left  = max(1, i - m)
        right = min(N, i + m)

        do j = left, right
           A_new(i) = A_new(i) + F(j - i + m + 1) * A_old(j)
        end do
     end do
     !$acc end parallel loop
     !$acc end data

     print *, "N =", N
     print *, "Output (center):", A_new(N/2)

     deallocate(A, A_old, A_new)
  end do

end program