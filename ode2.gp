reset
set xlabel "x"
set ylabel "y"
set xrange[-10:10]
set yrange[-10:10]
plot "ode2_Euler.dat" using 2:3 with lines
replot "ode2_RK2.dat" using 2:3 with lines
replot "ode2_RK4.dat" using 2:3 with lines
