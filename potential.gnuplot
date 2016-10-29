set terminal wxt
set xlabel "X"
set yrange [-0.015:0.015]

  p "potential_ad.dat" u 2:3 w l lc 1 ti "V_ad(x)", "potential_ad.dat" u 2:4 w l lc 1 notitle
rep  "potential_di.dat" u 2:3 w l lc 2 ti "V_di(x)", "potential_di.dat" u 2:4 w l lc 2 notitle 
 

pause -1
