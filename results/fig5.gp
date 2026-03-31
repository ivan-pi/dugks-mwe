set terminal pngcairo enhanced
set output "fig5.png"

set log xy
set xlabel "NX = NY"
set ylabel "L2-norm (velocity)"

set title "Re = 10^6, Ma = 10^{-4}"
set key top right

results="tgv_kraemer2020.txt"

plot results i 0 u 1:2 w lp t "SL (quadratic)",\
     results i 1 u 1:2 w lp t "FDM (isotropic)",\
     results i 2 u 1:2 w lp t "FDM (central)",\
     0.1*(1/x)**2 dt 2 lw 2 t "O(dx^2)"

results="tgv_kraemer2020_avg.txt"
set output "fig5_avg.png"

area = (2*pi)**2

plot results i 1 u 1:($2/area) w lp t "SL (quadratic)",\
     results i 2 u 1:($2/area) w lp t "FDM (isotropic)",\
     results i 3 u 1:($2/area) w lp t "FDM (central)",\
     100*(1.0/x)**2 dt 2 lw 2 t "O(dx^2)"
