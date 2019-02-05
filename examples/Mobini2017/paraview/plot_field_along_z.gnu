set term pdf 
set out 'field_along_z_mobini.pdf'

set key center box 
set key width 1.5 
set key height 2 
set yrange [-2:*]
set ylabel "Electric Field [V/m]" 
set xlabel 'z [mm]'

p 'field_along_z.tsv' u ($8*1e3):1 w l lw 3 lc rgb 'red' t 'E_x',\
'field_along_z.tsv' u ($8*1e3):2 w l lw 3 lc rgb 'blue' t 'E_y',\
'field_along_z.tsv' u ($8*1e3):3 w l lw 3 lc rgb 'green' t 'E_z',\
'field_along_z.tsv' u ($8*1e3):(sqrt($1**2+$2**2+$3**2)) w p lc rgb 'black' pt 5 ps .5 t '|E|'
