set terminal pngcairo size 900,600
set output "montecarlo_error.png"

set title "Monte Carlo Error Convergence"
set xlabel "Number of Samples (N)"
set ylabel "Integral Error Estimate"

set logscale x
set logscale y
set grid

V = 24
c = 1500   # scale constant so curves overlap

plot \
"montecarlo_stats.dat" using 1:(V*$4/sqrt($1)) \
with linespoints lw 2 pt 7 title "Estimated Error", \
c/sqrt(x) with lines lw 2 dt 2 title "Theoretical O(N^{-1/2})"