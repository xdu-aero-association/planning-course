% calculate x(t) y(t)
close all;clear;
syms t T a1 a2 b1 b2 r1 r2 J dpx dpy dvx dvy dax day
syms pxf pyf vxf vyf axf ayf px0 py0 vx0 vy0 ax0 ay0 jx jy
jx = 1/2*a1*t^2+b1*t+r1;
jy = 1/2*a2*t^2+b2*t+r2;
J_inner = jx^2+jy^2+1;

A = diag(-ones(4,1),-2);
Am = expm(A*t);
B = sym(zeros(6,6));
B(1,1) = 1/120*T^5;
B(1,3) = 1/24*T^4;
B(1,5) = 1/6*T^3;
B(2,2:6) = B(1,1:5);
B(3,1) = 1/24*T^4;
B(3,3) = 1/6*T^3;
B(3,5) = 1/2*T^2;
B(4,2:6) = B(3,1:5);
B(5,1) = 1/6*T^3;
B(5,3) = 1/2*T^2;
B(5,5) = T;
B(6,2:6) = B(5,1:5);
dpx = pxf-px0-vx0*T-1/2*ax0*T^2;
dvx = vxf-vx0-ax0*T;
dax = axf-ax0;
dpy = pyf-py0-vy0*T-1/2*ay0*T^2;
dvy = vyf-vy0-ay0*T;
day = ayf-ay0;
Result = B^(-1)*[dpx dpy dvx dvy dax day].';
a1 = Result(1);
a2 = Result(2);
b1 = Result(3);
b2 = Result(4);
r1 = Result(5);
r2 = Result(6);

J = int(J_inner,0,T);
J_T = subs(J);
% solve for dJ=0
dJ = diff(J_T,T);
dJ_sim = simplify(dJ);
dJ_sim = dJ_sim * T^6;
eqn = dJ == 0;
S_an = solve(eqn,T);
S = solve(eqn, T, 'MaxDegree', 4);        
% coeffs = sym2poly(dJ_sim);
[coeffs_dj,T_set] = coeffs(dJ_sim,T);