set term postscript
set out 'shuosher.eps'
set key font ",24"
set xtics font ",20"
set ytics font ",20"
set xlabel font ",24"
set ylabel font ",24"

set xlabel 'x'
set ylabel 'Density'
set xran[-5:5]
p 'sol' u 1:2 t 'DG pol' w l lw 2, \
  'avg' u 1:2 t 'DG avg' w p pt 6, \
  'shuosher.dat' u 1:2 t 'Exact' w l lw 2
