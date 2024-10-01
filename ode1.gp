set xlabel "x"
set ylabel "f(x)"
plot "ode1_Euler.dat" with lines
replot "ode1_RK2.dat" with lines
replot "ode1_RK4.dat" with lines
replot exp(-0.5*x*x)
