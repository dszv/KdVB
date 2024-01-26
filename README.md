# KdVB

Todo list:
- Code for exact solution [x]
- Code for perturbative solution [x]
- Code to extract alpha and beta vectors [x]
- Use renormalized solution as background solution [doing]

El codigo "lin_kdvb_solver.f90" da la solucion de la parte linerizada, es decir phi1, y corre hasta el tiempo necesario para el valor de v. La info de esta parte se guarda en "phi1.bin". Ambos resultados se deben separar en arrays del tamaño correspondiente para el tiempo total.

Obs: Borrar o mover los .bin antes de hacer push a github.

Obs2: ifort -r8 -O3 code.f90 -o a.out -L/home/diego/.local/src/fftw/lib -lfftw3 -I/home/diego/.local/src/fftw/include
