reset
set title "Bessel Function (nu = 1) with yL = 0, yR = 1" font ",18"

set ylabel "delta y(x)" font ",18"
plot "bessel.dat" u 1:($2-besj1($1)/besj1(10.0)) lc "red"
 
