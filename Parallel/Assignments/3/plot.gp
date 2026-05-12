# Output settings
set terminal pngcairo size 1000,700 enhanced
set output "cholesky_timing.png"

# General settings
set title "Cholesky Decomposition Performance"
set xlabel "Matrix Size (N)"
set ylabel "Time (seconds)"
set format y "10^{%L}"
set grid
set key left top

# Log scales (important for performance plots)
set logscale x
set logscale y

# Plot timing
plot "cholesky_results.dat" using 1:2 with linespoints lw 2 pt 7 title "Serial", \
     "cholesky_results.dat" using 1:3 with linespoints lw 2 pt 5 title "OpenACC"


# -------- Error Plot --------
set output "cholesky_error.png"
unset logscale x
set logscale y

set title "Cholesky Error"
set xlabel "Matrix Size (N)"
set ylabel "Error"

plot "cholesky_results.dat" using 1:4 with linespoints lw 2 pt 7 title "Relative Error (Serial and Parallel)", \
     "cholesky_results.dat" using 1:5 with linespoints lw 2 pt 5 title "Reconstruction Error (Parallel)"