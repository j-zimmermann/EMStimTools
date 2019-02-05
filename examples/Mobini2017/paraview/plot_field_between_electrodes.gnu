set term pdf 
set out 'field_cut_center_mobini.pdf'

set key center top box 
set key width 1.1 
set key height 1.1 
set xrange [-0.5:22.5]
set yrange [-25:190]
set ylabel "Electric Field [V/m]" 
set xlabel 'distance to left electrode [mm]'
set ytics -20,20,180

p 'field_between_electrodes.tsv' u (($7+0.011)*1e3):1 w l lw 3 lc rgb 'red' t 'E_x',\
'field_between_electrodes.tsv' u (($7+0.011)*1e3):2 w l lw 3 lc rgb 'blue' t 'E_y',\
'field_between_electrodes.tsv' u (($7+0.011)*1e3):3 w l lw 3 lc rgb 'green' t 'E_z',\
'field_between_electrodes.tsv' u (($7+0.011)*1e3):(sqrt($1**2+$2**2+$3**2)) w p lc rgb 'black' pt 5 ps .5 t '|E|',\
'field_tupaj.txt' u ($1*1e3):2 w l lw 3 lc rgb 'dark-violet' t 'Analytical'
