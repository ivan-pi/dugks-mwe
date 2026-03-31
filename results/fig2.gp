set terminal pngcairo enhanced
set output "fig2.png"

set log xy
set xlabel "NX = NY"
set ylabel "L2-norm (velocity)"

set title "Ma = 0.01, Re = 100, dt/tau = 2"

plot "tgv_zhu2017.txt" i 0 u 1:2 w lp t "DUGKS (Zhu, 2017)",\
     "tgv_zhu2017.txt" i 1 u 1:2 w lp t "Bardow-FVM (Zhu, 2017)",\
     "tgv_fig2.txt" i 3 u 1:2 w lp t "Bardow-FDM (present)",\
     "tgv_fig2.txt" i 4 u 1:2 w lp t "Bardow-FEM (present)",\
     "tgv_fig2.txt" i 5 u 1:2 w lp t "DUGKS (present)",\
     0.5*(1/x)**2 dt 2 t "O(dx^2)",\
     "tgv_fig2.txt" i 6 u 1:2 w lp t "FDM v2 (present)",\
     "tgv_fig2.txt" i 9 u 1:2 w lp t "FDM (isotropic))",\
     "tgv_fig2.txt" i 10 u 1:2 w lp t "SL-Parabolic"

