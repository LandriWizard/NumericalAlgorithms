reset
set xlabel "x"
set ylabel "y"
set log x
set log y
plot "ode2_Conv.dat" using 1:2
replot "ode2_Conv.dat" using 1:3
replot "ode2_Conv.dat" using 1:4
