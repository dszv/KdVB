# KdVB

Todo list:
- Code for exact solution [x]
- Code for perturbative solution [x]
- Code to extract alpha and beta vectors [doing]

El codigo "kdvb_solver" da la solucion numerica total, corre hasta el tiempo T = 300 s. La info se guarda en "soliton.bin". El codigo "lin_kdvb_solver.f90" da la solucion de la parte linerizada, es decir, phi 1 y corre hasta el tiempo T = 2400 s. La info de esta parte se guarda en "phi1.bin". Ambos resultados se deben separa en arrays del tamaño correspondiente para el tiempo total.

Para encontrar $ \Delta t $ se debe comparar el resultado de kdvb con el resultado de kdv + epsilon.phi1. El soliton de kdv se puede construir en python sin problemas.
Para hallar la matriz V solo basta con usar el resultado final de lin_kdvb (phi1).

Obs: Borrar o mover los .bin antes de hacer push a github.

Obs2: ifort -r8 -O3 cheby_fft_scan.f90 -o test.out -L/home/diego/.local/src/fftw/lib -lfftw3 -I/home/diego/.local/src/fftw/include
