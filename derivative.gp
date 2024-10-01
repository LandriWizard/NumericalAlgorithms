set logscale x
set logscale y
plot "derivative.dat" using(1.0/$1):2
replot "derivative.dat" using(1.0/$1):3
replot "derivative.dat" using(1.0/$1):4
