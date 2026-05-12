set terminal pngcairo size 1600,1200 enhanced
set output "q1p.png"

set multiplot layout 3,2 title "1D Wave Equation: Parallel Implementation"

set xlabel "x"
set ylabel "u(x,t)"
set grid
set xrange [0:2]

# ---------- t = 0 ----------
set title "t = 0  (np = 2)"
plot \
    "t0p.dat" using ($1==2?$2:1/0):3 w l lw 2 dt (20,10) title "Upwind", \
    "t0p.dat" using ($1==2?$2:1/0):4 w l lw 2 title "QUICK", \
    "t0p.dat" using ($1==2?$2:1/0):5 every 25 w p pt 6 ps 1.2 lc rgb "black" title "Exact"

set title "t = 0  (np = 4)"
plot \
    "t0p.dat" using ($1==4?$2:1/0):3 w l lw 2 dt (20,10) title "Upwind", \
    "t0p.dat" using ($1==4?$2:1/0):4 w l lw 2 title "QUICK", \
    "t0p.dat" using ($1==4?$2:1/0):5 every 25 w p pt 6 ps 1.2 lc rgb "black" title "Exact"

# ---------- t = 0.5 ----------
set title "t = 0.5  (np = 2)"
plot \
    "t05p.dat" using ($1==2?$2:1/0):3 w l lw 2 dt (20,10) title "Upwind", \
    "t05p.dat" using ($1==2?$2:1/0):4 w l lw 2 title "QUICK", \
    "t05p.dat" using ($1==2?$2:1/0):5 every 25 w p pt 6 ps 1.2 lc rgb "black" title "Exact"

set title "t = 0.5  (np = 4)"
plot \
    "t05p.dat" using ($1==4?$2:1/0):3 w l lw 2 dt (20,10) title "Upwind", \
    "t05p.dat" using ($1==4?$2:1/0):4 w l lw 2 title "QUICK", \
    "t05p.dat" using ($1==4?$2:1/0):5 every 25 w p pt 6 ps 1.2 lc rgb "black" title "Exact"

# ---------- t = 1.0 ----------
set title "t = 1.0  (np = 2)"
plot \
    "t1p.dat" using ($1==2?$2:1/0):3 w l lw 2 dt (20,10) title "Upwind", \
    "t1p.dat" using ($1==2?$2:1/0):4 w l lw 2 title "QUICK", \
    "t1p.dat" using ($1==2?$2:1/0):5 every 25 w p pt 6 ps 1.2 lc rgb "black" title "Exact"

set title "t = 1.0  (np = 4)"
plot \
    "t1p.dat" using ($1==4?$2:1/0):3 w l lw 2 dt (20,10) title "Upwind", \
    "t1p.dat" using ($1==4?$2:1/0):4 w l lw 2 title "QUICK", \
    "t1p.dat" using ($1==4?$2:1/0):5 every 25 w p pt 6 ps 1.2 lc rgb "black" title "Exact"

unset multiplot