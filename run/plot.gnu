set term postscript enhanced color
set out 'plot.ps'

unset key

set title 'Density'
plot 'sol' u 1:2 w l,'avg' u 1:2 w p

set title 'Velocity'
plot 'sol' u 1:3 w l,'avg' u 1:3 w p

set title 'Pressure'
plot 'sol' u 1:4 w l,'avg' u 1:4 w p

set title 'Mesh size'
plot 'h' u 2:3 w lp
