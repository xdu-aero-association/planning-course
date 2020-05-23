close all;clear;
syms t T a1 a2 a3 b1 b2 b3 J dpx dpy dpz dvx dvy dvz 
syms pxf pyf pzf vxf vyf vzf px0 py0 pz0 vx0 vy0 vz0
% A is sym matrix for lambda, exp_At is state transition matrix
A = [0 0 0;-1 0 0;0 -1 0]; % for ppt
A = zeros(6,6); % for hw4
A(4,1) = -1;
A(5,2) = -1;
A(6,3) = -1;
exp_At = expm(A*t);
% B is sym matrix for [alpha,beta,gamma]
B = sym(zeros(6,6));
B(1,1) = -12/T^3;
B(1,4) = 6/T^2;
B(2,2) = -12/T^3;
B(2,5) = 6/T^2;
B(3,3) = -12/T^3;
B(3,6) = 6/T^2;
B(4,1) = 6/T^2;
B(4,4) = -2/T;
B(5,2) = 6/T^2;
B(5,5) = -2/T;
B(6,3) = 6/T^2;
B(6,6) = -2/T;
dpx = pxf-vx0*T-px0;
dpy = pyf-vy0*T-py0;
dpz = pzf-vz0*T-pz0;
% if we want quadcopter to stay at [px py pz]
vxf = 0;
vyf = 0;
vzf = 0;
dvx = vxf-vx0;
dvy = vyf-vy0;
dvz = vzf-vz0;
Result = B*[dpx dpy dpz dvx dvy dvz].';
a1 = Result(1);
a2 = Result(2);
a3 = Result(3);
b1 = Result(4);
b2 = Result(5);
b3 = Result(6);
J = T+(1/3*a1^2*T^3+a1*b1*T^2+b1^2*T) ...
     +(1/3*a2^2*T^3+a2*b2*T^2+b2^2*T) ...
     +(1/3*a3^2*T^3+a3*b3*T^2+b3^2*T);
dJ = diff(J,T);
dJ_sim = simplify(dJ);
dJ_sim = dJ_sim * T^4;
eqn = dJ == 0;
S_an = solve(eqn,T);
S = solve(eqn,T, 'MaxDegree', 4);        
% coeffs = sym2poly(dJ_sim);
[coeffs_dj,T_set] = coeffs(dJ_sim,T);