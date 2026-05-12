set terminal pngcairo size 900,600
set output "jacobi_verification.png"

set datafile commentschars "#"

set xlabel "x"
set ylabel "phi(x,0.5)"
set title "Jacobi Solver Verification (Δ = 0.1)"
set grid

plot \
"phi_slice_all.dat" every :::0::0 u 1:2 w lp lw 2 pt 7 ps 1 title "Serial", \
"phi_slice_all.dat" every :::1::1 u 1:2 w lp lw 2 pt 5 ps 1 title "2 threads", \
"phi_slice_all.dat" every :::2::2 u 1:2 w lp lw 2 pt 9 ps 1 title "4 threads", \
"phi_slice_all.dat" every :::3::3 u 1:2 w lp lw 2 pt 11 ps 1 title "8 threads", \
"phi_slice_all.dat" every :::0::0 u 1:3 w l lw 3 dt 2 title "Exact"