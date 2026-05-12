set terminal pngcairo size 900,600
set output "runtime_vs_threads.png"

set xlabel "Number of Threads"
set ylabel "Runtime (seconds)"
set title "Jacobi Runtime vs Threads (Δ = 0.005)"

set grid
set key top right

plot \
"runtime_vs_threads.dat" using ($1==0.005 ? $2 : 1/0):3 \
with linespoints lw 2 pt 7 ps 1.5 title "Δ = 0.005"