set terminal wxt
set xlabel "X"
set yrange [-0.015:0.015]

  p "potential_ad.dat" u 2:3 w l lc 1 ti "V_ad(x)", "potential_ad.dat" u 2:4 w l lc 1 notitle
rep  "potential_di.dat" u 2:3 w l lc 2 ti "V_di(x)", "potential_di.dat" u 2:4 w l lc 2 notitle 
 
pause 1

do for [ii=0:10] {

replot "psisq.dat" i ii u 2:($3+0.0003) w l ls 3 lc rgb 'blue' notitle
replot "psisq.dat" i ii u 2:($4+0.0003) w l ls 3 lc rgb 'blue' notitle

pause 0.5
}

