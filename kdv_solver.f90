program fft_phi_1
use, intrinsic :: iso_c_binding; implicit none

include "fftw3.f03"

integer, parameter :: NN = 2**13, NN2 = NN/2 + 1
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0
real, parameter :: beta = 2.0
real, parameter :: x0 = 500.0
real, parameter :: t0 = 0.0
real, parameter :: dt = 0.001
integer, parameter :: TT = int(300.0/dt)+1, lapse = int(3.0/dt)
real, parameter :: epsilon = 0.1
real, parameter :: L = 2560.0, dx = L/NN, dk = 2.0*pi/L
! real(8) in(NN), x(NN), v(2*NN)
real(8) in(NN), van(NN), x(NN), v(NN)
type(C_PTR) :: plan, plan2
integer i, j
!! integer, parameter :: x_lim_1 = 0.0, x_lim_2 = int((540.0+L/2.0)/dx)
!! in = (0,0)

! forall (i=1:NN) x(i) = (i-1)*dx - L/2.0
forall (i=1:NN) x(i) = (i-1)*dx

! call kdv(in, x, 0.0)
call kdv(v, x, 0.0)
!call dump_sol(v, 1, 0.0)
! forall (i=1:NN) v(NN+i) = 0.0

do j = 1, TT
! call kdv(in, x, (j-1)*dt); forall (i=1:NN) v(i) = in(i);
    ! if (mod(j-1,lapse) ==0) then
    !     do i = 1, NN
    !         ! write (*,'(8G25.15e3)') (j-1)*dt, x(i), v(NN+i)
    !         write (*, '(8G25.15e3)') (j-1)*dt, x(i), v(i)
    !     end do
    !     write (*,'(8G25.15e3)') '',''; write (*,'(8G25.15e3)') '',''
    ! end if
 call gl8(v, dt, j)
 call kdv(van, x, j*dt)
 call dump_sol(v, j, j*dt)
end do

contains

subroutine evalf(v, dvdt, j)
! real v(2*NN), dvdt(2*NN), rhs_pert(NN), k(NN2), x(NN), phi(NN), phi_pert(NN), phi_prod(NN)
real v(NN), dvdt(NN), rhs(NN), k(NN2), x(NN), phi(NN), dphi(NN), phi_sqr(NN)
! complex(8) fft_kdv2(NN2), fft_kdv(NN2), fft_pert(NN2), fft_kdv_pert(NN2)
complex(8) fft_kdv(NN2), fft_kdv_sqr(NN2), fft_rhs(NN2)
integer j

! forall (i=1:NN) x(i) = (i-1)*dx - L/2.0
forall (i=1:NN) x(i) = (i-1)*dx
forall (i=1:NN2) k(i) = (i-1)*dk
forall (i=1:NN) phi(i) = v(i)
forall (i=1:NN) dphi(i) = dvdt(i)

! unmangle phase space state vector contents into human-readable form
! associate (phi => v(1:NN), phi_pert => v(NN+1:2*NN), dphi => dvdt(1:NN), dphi_pert => dvdt(NN+1:2*NN))
! associate (phi => v(1:NN), dphi => dvdt(1:NN))
! phi_prod = phi*phi_pert
forall (i=1:NN) phi_sqr(i) = phi(i)*phi(i)

if (j == 1) plan = fftw_plan_dft_r2c_1d(NN, phi, fft_kdv, FFTW_ESTIMATE)
call fftw_execute_dft_r2c(plan, phi, fft_kdv)
! call fftw_execute_dft_r2c(plan, phi_pert, fft_pert)
! call fftw_execute_dft_r2c(plan, phi_prod, fft_kdv_pert)
call fftw_execute_dft_r2c(plan, phi_sqr, fft_kdv_sqr)

do i = 1, NN2
! fft_pert(i) = -(0,1)*k(i)*k(i)*k(i)*fft_pert(i)+6.0*(0,1)*k(i)*fft_kdv_pert(i) + k(i)*k(i)*fft_kdv(i)
 fft_rhs(i) = -3.0*(0,1)*k(i)*fft_kdv_sqr(i) + (0,1)*k(i)*k(i)*k(i)*fft_kdv(i) - (1,0)*epsilon*k(i)*k(i)*k(i)*k(i)*fft_kdv(i)
end do

! if (j==1) plan2 = fftw_plan_dft_c2r_1d(NN,fft_pert,rhs_pert,FFTW_ESTIMATE)
if (j==1) plan2 = fftw_plan_dft_c2r_1d(NN, fft_rhs, rhs, FFTW_ESTIMATE)
! call fftw_execute_dft_c2r(plan2, fft_pert, rhs_pert)
call fftw_execute_dft_c2r(plan2, fft_rhs, rhs)

forall (i=1:NN) dphi(i) = rhs(i)/(1.0*NN)
forall (i=1:NN) dvdt(i) = dphi(i)

! dphi = 0.0
! dphi_pert = -rhs_pert/(1.0*(NN))

! end associate
end subroutine evalf

subroutine kdv(v, x, t)
real v(NN), x(NN), t
forall (i=1:NN) v(i) = beta/2.0*(1.0/cosh(sqrt(beta)/2.0*((x(i)-x0)-beta*(t-t0))))**2
end subroutine kdv

subroutine dump_sol(state, mark, time)
real state(NN), time
integer c2, mark

if (mark == 1) open(unit = 22, file = 'soliton.bin', access='stream', status='unknown')
if (mark> 1) open(unit = 22, file = 'soliton.bin', access='stream', status='old', position = 'append')

do c2 = 1, NN
    write (22) time, x(c2), state(c2)
end do
close (22)

end subroutine dump_sol

subroutine gl8(y, dt, j)
        ! integer, parameter :: s = 4, n = 2*NN
        integer, parameter :: s = 4, n = NN
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
