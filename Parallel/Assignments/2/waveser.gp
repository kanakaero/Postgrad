set terminal pngcairo size 1200,1000 enhanced
set output "q1s.png"

set multiplot layout 3,1 title "1D Wave Equation: Serial Implementation"

set xlabel "x"
set ylabel "u(x,t)"
set grid
set xrange [0:2]

set title "t = 0"
plot \
    "t0s.dat" using 1:2 w l lw 2 dt (20,10) title "Upwind", \
    "t0s.dat" using 1:3 w l lw 2 title "QUICK", \
    "t0s.dat" using 1:4 every 25 w p pt 6 ps 1.2 lc rgb "black" title "Exact"

set title "t = 0.5"
plot \
    "t05s.dat" using 1:2 w l lw 2 dt (20,10) title "Upwind", \
    "t05s.dat" using 1:3 w l lw 2 title "QUICK", \
    "t05s.dat" using 1:4 every 25 w p pt 6 ps 1.2 lc rgb "black" title "Exact"

set title "t = 1.0"
plot \
    "t1s.dat" using 1:2 w l lw 2 dt (20,10) title "Upwind", \
    "t1s.dat" using 1:3 w l lw 2 title "QUICK", \
    "t1s.dat" using 1:4 every 25 w p pt 6 ps 1.2 lc rgb "black" title "Exact"

unset multiplot