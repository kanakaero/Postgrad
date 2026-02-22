! ----------------------------------------------------------------------
! This subroutine deallocates all the allocated arrays at the end.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 14th November 2025
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

    if (allocated(u_flat)) deallocate(u_flat)
    if (allocated(v_flat)) deallocate(v_flat)
    if (allocated(U_mag)) deallocate(U_mag)
    if (allocated(u)) deallocate(u)
    if (allocated(v)) deallocate(v)

    if (allocated(T)) deallocate(T)
    if (allocated(qx)) deallocate(qx)
    if (allocated(qy)) deallocate(qy)
    if (allocated(qmag)) deallocate(qmag)

    if (allocated(R)) deallocate(R)

end subroutine dealloc
