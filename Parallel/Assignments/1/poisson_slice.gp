set terminal pngcairo size 900,600 enhanced font "Arial,12"
set output "phi_comparison.png"

set xlabel "x"
set ylabel "phi"
set title "Numerical vs Exact Solution"

set grid
set key top left

plot "phi_slice.dat" using 1:2 w lp lw 2 pt 7 title "Numerical", \
     "phi_slice.dat" using 1:3 w l lw 2 title "Exact"

unset output