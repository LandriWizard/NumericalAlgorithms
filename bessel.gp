reset
set grid
set title "Bessel Function (nu = 1) with yL = 0, yR = 1" font ",18"
# 1st plot, numerical & analytical solution
set xlabel "x" font ",18"
set ylabel "y(x)" font ",18"
plot "bessel.dat" u 1:2 t "Solution of the code" w p lc "red" ps 1.5
replot besj1(x)/besj1(10.0) t "Bessel Function" lc "blue"
