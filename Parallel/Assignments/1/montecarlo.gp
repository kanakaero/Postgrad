set terminal pngcairo size 900,600
set output "montecarlo_threads.png"

set title "Monte Carlo Integration Runtime"
set xlabel "Number of Threads"
set ylabel "Execution Time (seconds)"

set grid
set key top right

plot \
"montecarlo_timing.dat" using ($1==100 ? $2 : 1/0):3 with linespoints lw 2 pt 7 title "N = 100", \
"montecarlo_timing.dat" using ($1==1000 ? $2 : 1/0):3 with linespoints lw 2 pt 9 title "N = 1000", \
"montecarlo_timing.dat" using ($1==10000 ? $2 : 1/0):3 with linespoints lw 2 pt 5 title "N = 10^4", \
"montecarlo_timing.dat" using ($1==100000 ? $2 : 1/0):3 with linespoints lw 2 pt 11 title "N = 10^5", \
"montecarlo_timing.dat" using ($1==1000000 ? $2 : 1/0):3 with linespoints lw 2 pt 13 title "N = 10^6"