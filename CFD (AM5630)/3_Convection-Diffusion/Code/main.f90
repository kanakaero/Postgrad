! ----------------------------------------------------------------------
! Main program to solve the 2D Diffusion equation using the Finite Volume
! Method Discretization and a Gauss Seidel iterative solver.
! Supports both uniform and skewed meshes.
! ----------------------------------------------------------------------
! Author: Kanak Agarwal (kanakaero.github.io)
! Last Modified: 14th November 2025
! ----------------------------------------------------------------------

program main
    
    use var_dec

    call gridcon() ! Read the Domain Parameters and Generate the mesh

    call discret() ! Discretize the Convection-Diffusion Equation

    call bc() ! Apply Boundary Conditions

    call gauss_seidel() ! Solve the Discretized Equations using the Gauss Seidel Method

    call tdma() ! Solve the Discretized Equations using the TDMA Method

    call postproc() ! Post Processing

    call debug()

    call dealloc() ! Deallocate arrays

end program main
