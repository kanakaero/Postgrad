! ----------------------------------------------------------------------
! Main program to solve the 2D Diffusion equation using the Finite Volume
! Method Discretization and a Gauss Seidel iterative solver.
! Supports both uniform and skewed meshes.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 21st September 2025
! ----------------------------------------------------------------------

program main
    
    use var_dec

    call mesh_gen() ! Read the Domain Parameters and Generate the mesh

    call bc_s_k() ! Input the boundary conditions, source terms & solver parameters from the input file

    call discret() ! Discretize the governing equations

    call solver() ! Solve the discretized equations using Gauss Seidel method

    call dealloc() ! Deallocate arrays

end program main
