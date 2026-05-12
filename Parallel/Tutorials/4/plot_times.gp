set terminal pngcairo size 1000,700 enhanced font "Sans,12"
set output "scaling_mul.png"

set title "OpenMP Matrix Multiplication Scaling"
set xlabel "Matrix size N (NxN)"
set ylabel "Execution time (seconds)"

set grid
set key outside right
set logscale y
set datafile commentschars "#"

threads = "1 2 4 8"

plot for [i=1:words(threads)] \
    'execution_times_mul.txt' using \
    ($1==word(threads,i) ? $2 : 1/0): \
    ($1==word(threads,i) ? $3 : 1/0) \
    with linespoints lw 2 pt 7 \
    title sprintf("%s threads", word(threads,i))

unset output
