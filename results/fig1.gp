set terminal pngcairo enhanced
set output "fig1.png"

set log y
set xlabel "NX = NY"
set ylabel "L2-norm (velocity)"

set xrange [:200]

set title "Ma = 0.1, Re = 100, dt/tau = 2"

plot "tgv_fig1.txt" i 0 u (10 + $0*5):1 w lp t "Bardow (FDM)",\
     "tgv_fig1.txt" i 1 u (10 + $0*5):1 w lp t "Bardow (FEM)",\
     "tgv_fig1.txt" i 2 u (10 + $0*5):1 w lp t "DUGKS",\
     "tgv_fig1.txt" i 3 u 1:2 w lp t "FDM (isotropic)",\
     "tgv_fig1.txt" i 4 u 1:2 w lp t "SL (quadratic)",\
     2*(1/x)**2 dt 2 t "O(dx^2)"
