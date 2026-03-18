set terminal pngcairo enhanced
set output "fig1.png"

set log xy
set xlabel "NX = NY"
set ylabel "L2-norm (velocity)"

set xrange [:160]

set title "Ma = 0.1, Re = 100, dt/tau = 2"

plot "tgv_fig1.txt" i 0 u (10 + $0*5):1 w lp t "Bardow (FDM)",\
     "tgv_fig1.txt" i 1 u (10 + $0*5):1 w lp t "Bardow (FEM)",\
     "tgv_fig1.txt" i 2 u (10 + $0*5):1 w lp t "DUGKS",\
     2*(1/x)**2 dt 2 t "O(dx^2)"
