program fft_phi_1
use, intrinsic :: iso_c_binding; implicit none

include "fftw3.f03"

integer, parameter :: NN = 2**13, NN2 = NN/2 + 1
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0
real, parameter :: vel = 2.0

real, parameter :: x0 =100.0
real, parameter :: t0 = 0.0

real, parameter :: dt = 0.001
integer, parameter :: TT = int(2400.0/dt) ! , lapse = int(2400.0/dt)
real, parameter :: cut = real(int(real(TT)/5000.0))

! real, parameter :: epsilon = 0.1

real, parameter :: L = 2560.0, dx = L/NN, dk = 2.0*pi/L

real(8) van(NN), x(NN), v(NN)
type(C_PTR) :: plan, plan2
integer i, j

print *, cut

forall (i=1:NN) x(i) = (i-1)*dx

print *, "KdVB linearized solution"
print *, "Total time:", TT*dt
print *, "x0:", x0
print *, "v:", vel

forall (i=1:NN) v(i) = 0.0
call dump_sol(v, 0, 0.0)

do j = 1, TT
 call gl8(v, dt, j, vel)
 call kdv(van, x, vel, j*dt)
 if (mod(real(j), cut) == 0.0) then
  print *, "Progress (%): ", 100.0*real(j)/TT
  call dump_sol(v, j, j*dt)
 endif
end do

contains

subroutine evalf(v, dvdt, vel, j)
! real v(NN), dvdt(NN), rhs(NN), k(NN2), x(NN), phi(NN), dphi(NN), phi_sqr(NN)
real v(NN), dvdt(NN), rhs(NN), k(NN2), x(NN), v_nl(NN), back(NN)
complex(8) fft_kdv(NN2), fft_kdv_nl(NN2), fft_kdv_back(NN2), fft_rhs(NN2)
real vel
integer j

forall (i=1:NN) x(i) = (i-1)*dx
forall (i=1:NN2) k(i) = (i-1)*dk
! forall (i=1:NN) phi(i) = v(i)
! forall (i=1:NN) dphi(i) = dvdt(i)

! unmangle phase space state vector contents into human-readable form
! forall (i=1:NN) phi_sqr(i) = phi(i)*phi(i)
! forall (i=1:NN) v_sqr(i) = v(i)*v(i)
call kdv(back, x, vel, j*dt)
forall (i=1:NN) v_nl(i) = back(i)*v(i)

if (j == 1) plan = fftw_plan_dft_r2c_1d(NN, v, fft_kdv, FFTW_ESTIMATE)
call fftw_execute_dft_r2c(plan, v, fft_kdv)
call fftw_execute_dft_r2c(plan, v_nl, fft_kdv_nl)
call fftw_execute_dft_r2c(plan, back, fft_kdv_back)

do i = 1, NN2
 ! fft_rhs(i) = -3.0*(0,1)*k(i)*fft_kdv_sqr(i) + (0,1)*k(i)*k(i)*k(i)*fft_kdv(i) - (1,0)*epsilon*k(i)*k(i)*k(i)*k(i)*fft_kdv(i)
  fft_rhs(i) = -6.0*(0,1)*k(i)*fft_kdv_nl(i) + (0,1)*k(i)*k(i)*k(i)*fft_kdv(i) - k(i)*k(i)*k(i)*k(i)*fft_kdv_back(i)
end do

if (j==1) plan2 = fftw_plan_dft_c2r_1d(NN, fft_rhs, rhs, FFTW_ESTIMATE)
call fftw_execute_dft_c2r(plan2, fft_rhs, rhs)

! forall (i=1:NN) dphi(i) = rhs(i)/(1.0*NN)
! forall (i=1:NN) dvdt(i) = dphi(i)
forall (i=1:NN) dvdt(i) = rhs(i)/(1.0*NN)

end subroutine evalf

subroutine kdv(v, x, vel, t)
real v(NN), x(NN), vel, t

forall (i=1:NN) v(i) = vel/2.0*(1.0/cosh(sqrt(vel)/2.0*((x(i)-x0)-vel*(t-t0))))**2

end subroutine kdv

subroutine dump_sol(v, mark, time)
real v(NN), time
integer c2, mark

if (mark == 0) open(unit = 22, file = 'phi1.bin', access='stream', status='unknown')
if (mark > 0) open(unit = 22, file = 'phi1.bin', access='stream', status='old', position = 'append')

do c2 = 1, NN
    write (22) time, x(c2), v(c2)
end do
close (22)

end subroutine dump_sol

subroutine gl8(y, dt, j, vel)
integer, parameter :: s = 4, n = NN
real y(n), g(n,s), dt; integer i, k, j
real vel

! Butcher tableau for 8th order Gauss-Legendre method
real, parameter :: a(s,s) = reshape((/ &
         0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
         0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
         0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
        -0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
         0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
         0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
         0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
         0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
real, parameter ::   b(s) = (/ &
         0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
         0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)

! iterate trial steps
g = 0.0; do k = 1,16
        g = matmul(g,a)
        do i = 1,s
                call evalf(y + g(:,i)*dt, g(:,i), vel, j)
        end do
end do

! update the solution
y = y + matmul(g,b)*dt

end subroutine gl8

end program
