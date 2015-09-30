set xlabel 'x'
set ylabel 'Density'
p 'sol' u 1:2 t 'DG' w l lw 2, \
  'avg' u 1:2 t 'Avg' w p pt 6, \
  'lax.dat' u 1:2 t 'Exact' w l lw 2
