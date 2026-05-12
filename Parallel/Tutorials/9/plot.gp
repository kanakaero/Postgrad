set terminal pngcairo size 1000,600 enhanced
set output "plot.png"

set grid
set xlabel "x"
set ylabel "Derivative"
set title "Derivative Comparison: Exact vs Pade vs CDS"

plot "output.dat" using 1:2 with lines lw 2 title "Exact", \
     "output.dat" using 1:3 with lines lw 2 title "Pade", \
     "output.dat" using 1:4 with lines lw 2 title "CDS"

unset output