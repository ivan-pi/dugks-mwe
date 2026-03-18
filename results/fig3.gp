set terminal pngcairo enhanced
set output "fig3.png"

set log xy
set xlabel "NX = NY"
set ylabel "L2-norm (velocity)"

plot "tgv_wu2018.txt" i 0 u 1:2 w lp t "T2S2 (Wu et al., 2018)",\
     "tgv_wu2018.txt" i 1 u 1:2 w lp t "T2S3 (Wu et al., 2018)",\
     "tgv_wu2018.txt" i 2 u 1:2 w lp t "T3S3 (Wu et al., 2018)",\
     "tgv_wu2018.txt" i 6 u 1:2 w lp t "DUGKS (GXW, 2013)",\
     13*(1/x)**2 dt 2 lw 2 t "O(dx^2)", \
     130*(1/x)**3 dt 4 lw 2 t "O(dx^3)", \
     "tgv_wu2018.txt" i 3 u 1:2 w lp t "Bardow-FEM (present)",\
     "tgv_wu2018.txt" i 4 u 1:2 w lp t "Bardow-FDM (present)",\
     "tgv_wu2018.txt" i 5 u 1:2 w lp t "DUGKS (present)"
