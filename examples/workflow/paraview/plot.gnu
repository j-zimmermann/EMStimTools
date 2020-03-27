set term pdf lw 4

set out 'plot_comparison_numerical_analytical.pdf'
set datafile separator ','
set xrange [0:22]
set key center
set xlabel "distance to left electrode [mm]"
set ylabel "Electric field strength [V/m]"
p 'field_tupaj.txt' u ($1*1e3):2 w l t 'analytical' ,\
'data_field.txt' u ($5*1e3):(sqrt($1*$1+$2*$2+$3*$3)) w l t 'numerical'
