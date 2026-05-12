set terminal pngcairo size 900,600
set output "runtime_vs_grid.png"

set xlabel "Grid spacing (Δ)"
set ylabel "Runtime (seconds)"
set title "Runtime vs Grid Size for 8 Threads"

set grid
set logscale x
set logscale y
set key right top

plot \
"runtime_vs_threads.dat" using ($2==8 ? $1 : 1/0):3 \
with linespoints lw 2 pt 7 title "Parallel (8 threads)"