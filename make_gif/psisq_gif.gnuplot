set xlabel "X"
#set yrange [-0.015:0.015] #for pot_type=1
set yrange [0:0.2]       #for pot_type=2

set term png size 600,300

do for [ii=0:99]{
outputfile = "./psisq_".ii.".png"
set output outputfile
  p "../potential_ad.dat" u 2:($3+0.000) w l lt rgb "green"  lw 2 ti "V_{adiab}(x)",\
    "../potential_ad.dat" u 2:($4+0.000) w l lt rgb "green"  lw 2 notitle,\
    "../potential_di.dat" u 2:($3+0.003) w l lt rgb "orange" lw 2 ti "V_{diab}(x)",\
    "../potential_di.dat" u 2:($4+0.003) w l lt rgb "orange" lw 2 notitle,\
    "../psisq.dat" i ii   u 2:($3+0.006) w l lc rgb "blue"   lw 1 ti "psisq_{1}(x)",\
    "../psisq.dat" i ii   u 2:($4+0.006) w l lc rgb "red"    lw 1 ti "psisq_{2}(x)"
#pause -1
}
