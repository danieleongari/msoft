set yrange [-1:2]
set xrange [*:*]
set xlabel "step"
set ylabel "population"
p "population.dat" u 1:2 w l ti "population C1", "population.dat" u 1:3 w l ti "population C2"
pause -1
