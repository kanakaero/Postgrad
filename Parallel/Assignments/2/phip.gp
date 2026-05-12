set encoding utf8
set terminal pngcairo size 1200,1000 enhanced
set output "q2p.png"

set multiplot layout 2,1 title "Jacobi Solver Comparison (Δx = 0.01)"

set grid
set key top center offset 0,-2

# ---------- Plot 1: φ vs x at y = 0 ----------
set xlabel "x"
set ylabel "𝜙"
set title "𝜙 vs x at y = 0"

plot \
    "phi_y0_p2.dat"   using 1:2 with lines lw 2 title "MPI p=2", \
    "phi_y0_p4.dat"   using 1:2 with lines lw 2 title "MPI p=4", \
    "phi_y0_p8.dat"   using 1:2 with lines lw 2 title "MPI p=8", \
    "phi_y0_pt01.dat" every 3 using 1:2 with points pt 6 ps 1.4 lc rgb "black" title "Serial"

# ---------- Plot 2: φ vs y at x = 0 ----------
set xlabel "y"
set ylabel "𝜙"
set title "𝜙 vs y at x = 0"

plot \
    "phi_x0_p2.dat"   using 1:2 with lines lw 2 title "MPI p=2", \
    "phi_x0_p4.dat"   using 1:2 with lines lw 2 title "MPI p=4", \
    "phi_x0_p8.dat"   using 1:2 with lines lw 2 title "MPI p=8", \
    "phi_x0_pt01.dat" every 3 using 1:2 with points pt 6 ps 1.4 lc rgb "black" title "Serial"

unset multiplot