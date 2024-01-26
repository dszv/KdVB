program ode_solver
use, intrinsic :: iso_c_binding; implicit none

! include "fftw3.f03"

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

real, parameter :: v0 = 2.0
real, parameter :: x0 = 300.0
real, parameter :: epsilon = 0.1

real, parameter :: t0 = 0.0
real, parameter :: dt = 0.001
integer, parameter :: TT = int(900.0/dt)
real, parameter :: cut = real(int(real(TT)/5000.0))

real(8) X(2)

integer j

print *, "ODE Solver"
print *, "Total time:", TT*dt
print *, "x0:", x0
print *, "v0:", v0
print *, "epsilon:", epsilon

X(1) = v0 + epsilon * (0.18873650869279932 * v0**(1.4788388361331373)) ! initial condition for vR
X(2) = x0 + epsilon * (0.6437723750396211 * v0**2 + 0.08475603754618799 * v0 + 0.10742185262873971 ) ! initial condition for xR
call dump_sol(X, 0, 0.0)

do j = 1, TT
 call gl8(X, dt, j)
 if (mod(real(j), cut) == 0.0) then
  print *, "Progress (%): ", 100.0*real(j)/TT
  call dump_sol(X, j, j*dt)
 endif
end do

contains

subroutine evalf(X, dXdt, j)
real X(2), dXdt(2)
integer j

dXdt(1) = -epsilon * 0.1904765191012984 * X(1)**(2.999989047049034) ! eqn for vR
dXdt(2) = X(1) ! eqn for xR
end subroutine evalf

subroutine dump_sol(XR, mark, time)
real XR(2), time
integer mark

if (mark == 0) open(unit = 22, file = 'lambdaR.bin', access='stream', status='unknown')
if (mark > 0) open(unit = 22, file = 'lambdaR.bin', access='stream', status='old', position = 'append')

write (22) time, XR(1), XR(2)
close (22)

end subroutine dump_sol

subroutine gl8(y, dt, j)
integer, parameter :: s = 4, n = 2
real y(n), g(n,s), dt; integer i, k, j

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
                call evalf(y + g(:,i)*dt, g(:,i), j)
        end do
end do

! update the solution
y = y + matmul(g,b)*dt

end subroutine gl8

end program
