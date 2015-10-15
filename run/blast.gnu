set term postscript
set out 'blast.eps'
set key font ",24"
set xtics font ",20"
set ytics font ",20"
set xlabel font ",24"
set ylabel font ",24"

set xlabel 'x'
set ylabel 'Density'
set xran[0.5:1.0]
p 'sol' u 1:2 t 'DG pol' w l lw 2, \
  'avg' u 1:2 t 'DG avg' w p pt 6, \
  'blast.dat' u 1:2 t 'Exact' w l lw 2
