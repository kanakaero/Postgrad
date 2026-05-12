set terminal pngcairo size 1200,800
set output "derivative.png"

set xlabel "x"
set ylabel "du/dx"
set grid
set key outside

plot \
"derivative_all.txt" u 1:2 w l lw 3 title "Analytical", \
"" u 1:3 w l lw 2 title "Forward", \
"" u 1:4 w l lw 2 title "Backward", \
"" u 1:5 w l lw 2 title "Central", \
"" u 1:6 w l lw 2 title "Central 4th"