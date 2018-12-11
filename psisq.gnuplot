set xlabel "X"
#set yrange [-0.015:0.015] #for pot_type=1
set yrange [0:0.1]       #for pot_type=2

  p "potential_ad.dat" u 2:3 w l lc 1 ti "V_{adiab}(x)", "potential_ad.dat" u 2:4 w l lc 1 notitle,\
    "potential_di.dat" u 2:3 w l lc 2 ti "V_{diab}(x)", "potential_di.dat" u 2:4 w l lc 2 notitle

pause -1

do for [ii=0:900]{
  p "potential_ad.dat" u 2:3 w l lc 1 ti "V_{adiab}(x)", "potential_ad.dat" u 2:4 w l lc 1 notitle,\
    "potential_di.dat" u 2:3 w l lc 2 ti "V_{diab}(x)", "potential_di.dat" u 2:4 w l lc 2 notitle,\
    "psisq.dat" i ii u 2:($3+0.0003) w l ls 3 lc rgb 'blue' notitle,\
    "psisq.dat" i ii u 2:($4+0.0003) w l ls 3 lc rgb 'blue' notitle
pause -1
}
