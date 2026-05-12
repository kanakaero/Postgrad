set encoding utf8
set terminal pngcairo size 1200,900 enhanced
set output "q2s.png"

set grid
set key top left
set xlabel "x or y"
set ylabel "𝜙"
set title "𝜙 Profiles from Jacobi Solver (Δx = 0.1)"

plot "phi_y0.dat" using 1:2 with linespoints lw 2 pt 7 title "𝜙 vs x (y=0)", \
     "phi_x0.dat" using 1:2 with linespoints lw 2 pt 5 title "𝜙 vs y (x=0)"