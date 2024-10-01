reset

set autoscale
set grid

wave = "wave.dat"

plot wave i 0 u 1:2 w l t "mode 1"
replot wave i 1 u 1:2 w l t "mode 2"
replot wave i 2 u 1:2 w l t "mode 3"
replot wave i 3 u 1:2 w l t "mode 4"
replot wave i 4 u 1:2 w l t "mode 5"
replot wave i 5 u 1:2 w l t "mode 6"
