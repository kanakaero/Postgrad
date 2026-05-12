! OpenACC parallelized Modified Gram Schmidt Algorithm for QR Decomposition where M = 1500 and N = 1000

program modified_gram_schmidt
   implicit none
   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: M = 15, N = 10
   real(dp), dimension(M, N) :: A, A_orig, R
   real(dp), dimension(N) :: norm
   integer :: i, j, k, p, q
   real(dp) :: tmp, rtmp
   real(dp) :: err, val, recon_err
   real(dp), parameter :: eps = 1.0d-12

   print *, "Starting Modified Gram-Schmidt QR Decomposition with M =", M, "and N =", N

   ! Initialize matrix A
   do j = 1, N
      do i = 1, M
         call random_number(A(i, j))
      end do
   end do

   print *, "Matrix A initialized with random values."

   ! Store original matrix for reconstruction check
   A_orig = A

   !$acc data copy(A, R, norm)
   ! Modified Gram-Schmidt Algorithm

   do j = 1, N

      ! ---- Compute norm ----
      tmp = 0.0
      !$acc parallel loop reduction(+:tmp)
      do i = 1, M
         tmp = tmp + A(i, j)**2
      end do
      norm(j) = sqrt(tmp)

      if (norm(j) < eps) then
         print *, "Warning: near-zero norm at column ", j
         norm(j) = eps
      end if

      R(j, j) = norm(j)

      !$acc update device(norm(j))
      ! ---- Normalize column j ----
      !$acc parallel loop
      do i = 1, M
         A(i, j) = A(i, j)/norm(j)
      end do

      !$acc wait

      ! ---- Orthogonalize remaining columns ----

      do k = j + 1, N

         rtmp = 0.0
         !$acc parallel loop reduction(+:rtmp)
         do i = 1, M
            rtmp = rtmp + A(i, j)*A(i, k)
         end do
         R(j, k) = rtmp

         !$acc update device(R(j,k))
         !$acc parallel loop
         do i = 1, M
            A(i, k) = A(i, k) - R(j, k)*A(i, j)
         end do

      end do
      
   end do

   !$acc end data   ! bring A, R, norm back to CPU

   ! ORTHOGONALITY CHECK: Q^T Q = I
   err = 0.0

   do p = 1, N
      do q = 1, N
         val = dot_product(A(:, p), A(:, q))
         if (p == q) then
            err = err + abs(val - 1.0)
         else
            err = err + abs(val)
         end if
      end do
   end do

   print *, "Orthogonality error =", err

   ! RECONSTRUCTION CHECK: A ≈ Q R
   recon_err = 0.0

   do j = 1, N
      do i = 1, M
         val = 0.0
         do k = 1, N
            val = val + A(i, k)*R(k, j)
         end do
         recon_err = recon_err + abs(val - A_orig(i, j))
      end do
   end do

   print *, "Reconstruction error =", recon_err

end program modified_gram_schmidt
