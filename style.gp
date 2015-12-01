set encoding utf8

unset key

# Color
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green
set style line 3 lc rgb '#800000' lt 1 lw 2
set style line 4 lc rgb '#ff0000' lt 1 lw 2
set style line 5 lc rgb '#ff4500' lt 1 lw 2
set style line 6 lc rgb '#ffa500' lt 1 lw 2
set style line 7 lc rgb '#006400' lt 1 lw 2
set style line 8 lc rgb '#0000ff' lt 1 lw 2
set style line 9 lc rgb '#9400d3' lt 1 lw 2

# Grey
set style line 11 lc rgb 'gray30' lt 1 lw 2
set style line 12 lc rgb 'gray40' lt 1 lw 2
set style line 13 lc rgb 'gray70' lt 1 lw 2
set style line 14 lc rgb 'gray90' lt 1 lw 2
set style line 15 lc rgb 'black' lt 1 lw 1.5
set style line 16 lc rgb "black" lt 4 lw 2

# Borders, etc.
set style line 21 lc rgb 'black' lt 1 lw 1.5
set border 3 back ls 21
set tics nomirror
set style line 22 lc rgb 'grey20' lt 0 lw 1
set grid back ls 22

# Palette
set palette maxcolors 3
set palette defined ( 0 '#8b1a0e',\
                      1 '#5e9c36',\
                      2 '#800000' )

set xtics 1
set ytics 1
