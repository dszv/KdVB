
program wave
use, intrinsic :: iso_c_binding; implicit none
include "fftw3.f03"

! solver control parameters
integer, parameter :: nn = 2**11                 ! number of nodes to sample on (i.e. spectral order)
real, parameter :: ell = 200.0
real, parameter :: eps = 0.0
real, parameter :: vel = 7.8

! this is exactly what you think it is...
real, parameter :: pi = 3.1415926535897932384626433832795028841971694Q0

real r(nn), red(nn), gamma(nn)
! phase space state vector packs [phi,u,v,w] into a single array
real state(nn), F(nn), der_an(nn), der_num(nn), omega_span(50), phase_init_span(50)
integer i, j, ind(1), break
integer, parameter :: init_phase = 1
integer, parameter :: final_phase = 5

real, parameter :: dt = 0.001               ! time step size (simulated timespan is tt*dt)
integer, parameter :: tt = int(2000.0/dt)                ! total number of time steps to take (i.e. runtime)
real, parameter :: cut = real(int(real(tt)/5000.0))
integer, parameter :: out_amp = int(tt/cut)
real amp_vec(out_amp)
call radial_grid(r, gamma);
do i = 0,49
omega_span(i+1) = 10.0**(-1.0+0.02*i)
phase_init_span(i+1) = pi*i/50.0
end do

!do i = init_phase, final_phase; do j =1,50
call evol_point(vel, phase_init_span(1), amp_vec)
 !   call dump(amp_vec, 1, 1)
write (*, '(6G24.12e3)') "mode", 50*(i-1)+j, "finished..."
!end do; end do


contains

subroutine evol_point(omega, phase_init, amp_vec)
real, intent(in) :: omega
real, intent(in) :: phase_init
real, intent(out) :: amp_vec(out_amp)
real state_evol(nn)
integer count, j

call init_KdV_soliton(state_evol, omega, phase_init, r)
j=0
amp_vec = 0.0
do count = 1,tt
        call gl8(state_evol, dt, count)
        if (mod(real(count), cut) == 0.0) then
        j = j+1
        amp_vec(j) = state_evol(1)
        call dump_sol(state_evol, j, count*dt)
        end if
end do
end subroutine evol_point


function Dr(phi, time_count)
real Dr(nn), phi(nn)
real dct_phi(nn), dct_shft_phi(nn), mapping(nn)
type(C_PTR) :: plan_dct, plan_dst
integer i, time_count

if (time_count==1) then
plan_dct = fftw_plan_r2r_1d(nn, phi, dct_phi, FFTW_REDFT10, FFTW_MEASURE)
end if

call fftw_execute_r2r(plan_dct, phi, dct_phi)
forall (i=1:nn) mapping(i) = -2.0*ell/(ell*ell+r(i)*r(i))
forall (i=1:nn) dct_shft_phi(i) = -real(i-1)*dct_phi(i)
forall (i=1:nn-1) dct_shft_phi(i) = dct_shft_phi(i+1)
dct_shft_phi(nn) = 0.0

if (time_count==1) then
plan_dst = fftw_plan_r2r_1d(nn, dct_shft_phi, Dr, FFTW_RODFT01, FFTW_MEASURE)
end if
call fftw_execute_r2r(plan_dst, dct_shft_phi, Dr)

Dr = -Dr*mapping/(2.0*nn)
end function Dr

subroutine radial_grid(r, gamma)
real red(nn)
real, dimension(nn), intent(out) :: r, gamma
integer j, i
integer, parameter :: pml = 2**4
forall (j=1:nn) red(j) = (nn-j+0.5)*pi/real(nn);
forall (j=1:nn) r(j) = ell*cos(0.5*red(j))/sin(0.5*red(j))
! forall (j=1:nn) gamma(j) = exp(-real(nn-j)**2/(1.0*pml)**2)+&
! exp(-real(j-1)**2/pml**2) - 1.25e-4
forall (i=1:nn) gamma(i) = exp(-real(nn-i)**2/pml**2) - 1.25e-4
! where (gamma<0.0) gamma = 0.0; gamma = 0.5/dt/gamma(1)*gamma
where (gamma<0.0) gamma = 0.0; gamma = (1.0/dt)/gamma(nn)*gamma
end subroutine radial_grid

! initial breather
subroutine init_KdV_soliton(state, omega, phase_init, r)
real, dimension(nn), intent(out) :: state
real, intent(in) :: omega, phase_init
real KdV_soliton(nn), sech_aux(nn)
real, dimension(nn), intent(in) :: r
integer i

sech_aux = 1.0/cosh((sqrt(vel)/2)*(r - 30.0))
KdV_soliton = (vel/2)*(sech_aux**2)
do i = 1, nn
    state(i) = KdV_soliton(i)
end do
end subroutine init_KdV_soliton


! evaluate equations of motion
subroutine evalf(y, dydt, time_count)
        real, dimension(nn) :: y, dydt
        real t, phi(nn), dphi(nn)
        integer fill, time_count
! break state into pieces
do fill = 1, nn
    phi(fill) = y(fill)
end do

! set pieces in EOMs
dphi = -6.0 * phi*Dr(phi, time_count) - Dr(Dr(Dr(phi, time_count), time_count), time_count)

! build derivatives
do fill = 1, nn
    dydt(fill) = dphi(fill)
end do
end subroutine evalf


subroutine dump_sol(state, mark, time)
real state(nn), time
integer c2, mark

if (mark == 1) open(unit = 22, file = 'soliton.bin', access='stream', status='unknown')
if (mark> 1) open(unit = 22, file = 'soliton.bin', access='stream', status='old', position = 'append')

do c2 = 1, nn
    write (22) time, r(c2), state(c2)
end do
close (22)

end subroutine dump_sol

subroutine gl8(y, dt, time_count)
        integer, parameter :: s = 4, n = nn
        real y(n), g(n,s), dt; integer i, k, time_count
        
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
                        call evalf(y + g(:,i)*dt, g(:,i), time_count)
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl8
end
