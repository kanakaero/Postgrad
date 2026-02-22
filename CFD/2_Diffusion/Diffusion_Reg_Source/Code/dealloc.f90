! ----------------------------------------------------------------------
! This subroutine deallocates all the allocated arrays at the end.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 20th September 2025
! ----------------------------------------------------------------------

subroutine dealloc()
    use var_dec

    ! Deallocate all allocated arrays
    if (allocated(dx)) deallocate(dx)
    if (allocated(dy)) deallocate(dy)
    if (allocated(xc)) deallocate(xc)
    if (allocated(yc)) deallocate(yc)
    if (allocated(xg)) deallocate(xg)
    if (allocated(yg)) deallocate(yg)
    if (allocated(a_p)) deallocate(a_p)
    if (allocated(a_e)) deallocate(a_e)
    if (allocated(a_w)) deallocate(a_w)
    if (allocated(a_n)) deallocate(a_n)
    if (allocated(a_s)) deallocate(a_s)
    if (allocated(kappa)) deallocate(kappa)
    if (allocated(T)) deallocate(T)
    if (allocated(R)) deallocate(R)

end subroutine dealloc
