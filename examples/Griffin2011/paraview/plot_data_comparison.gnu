set term pdf 
set out 'field_Griffin.pdf'
set key center top
set grid x y2
unset grid
set ytics nomirror
set xtics nomirror
set tics out

set yrange [*:1e-7]
set y2range [*:.1]
set ylabel "f=15 Hz, {/Symbol e}_r=2.5\nE-Field strength [V/m]" 
set y2label  "E-Field strength [V/m]\nf=60 kHz, {/Symbol e}_r=3000"
set y2tics
set xlabel 'y [mm]'
set format y "%.0e"

p 'field_griffin.tsv' u ($6*1e3):2 w l lw 4 axes x1y1 t 'f=15 Hz, {/Symbol e}_r=2.5', 'field_griffinHFeps.tsv' u ($6*1e3):2 w l axes x1y2 t 'f=60 kHz, {/Symbol e}_r=3000' lw 4
