set terminal pngcairo enhanced
set output "fig4.png"

set log xy
set xlabel "NX = NY"
set ylabel "L2-norm (velocity)"

set key bottom right

plot "tgv_sm2024.txt" i 0 u 1:2 w lp t "SL (quadratic), CFL=0.1",\
     "tgv_sm2024.txt" i 0 u 1:3 w lp t "SL (quadratic), CFL=0.4",\
     "tgv_sm2024.txt" i 1 u 1:2 w lp t "FDM (isotropic), CFL=0.1",\
     "tgv_sm2024.txt" i 1 u 1:3 w lp t "FDM (isotropic), CFL=0.4",\
     "tgv_sm2024.txt" i 2 u 1:2 w lp t "FDM (central), CFL=0.1",\
     "tgv_sm2024.txt" i 2 u 1:3 w lp t "FDM (central), CFL=0.4",\
     "tgv_sm2024.txt" i 3 u 1:2 w lp t "SL (quadratic), CFL=0.1, umax = 0.01/2",\
     "tgv_sm2024.txt" i 4 u 1:2 w lp t "SL (quadratic), CFL=0.1, umax = 0.01/4",\
     "tgv_sm2024.txt" i 5 u 1:2 w lp t "SL (quadratic), CFL=0.1, umax = 0.01/8",\
     "tgv_sm2024.txt" i 6 u 1:2 w lp t "SL (quadratic), Re = 100, dt/tau=8/3, umax = 0.01",\
     "tgv_sm2024.txt" i 7 u 1:2 w lp t "SL (quadratic), Re = 100, dt/tau=8/3, umax = 0.01/2",\
     "tgv_sm2024.txt" i 8 u 1:2 w lp t "SL (quadratic), Re = 100, dt/tau=8/3, umax = 0.01/4",\
     13*(1/x)**2 dt 2 lw 2 t "O(dx^2)",\
     30*(1/x)**4 dt 2 lw 2 t "O(dx^4)"
