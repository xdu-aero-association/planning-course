% we can directly change t with x, T with xf
% to get solution for y(x) = poly(x)
% Note: we are using region around end point to get a polyfit curve
% and it's a polynomial of T(usually 3rd order)
% Note: we should use interior point method's indicator function 
% I(u) = (-1/t)*log(-u), u = Ax-b<=0 inequality constraints,
% u<0, I is valid, u>0, I is inf, t is a param and t bigger, I closer to
% real indicator
% add I(u) to cost fun to bound the T value
close all;clear;
syms t T a b r J dp dv
syms pf vf af p0 v0 a0 j
syms c0 c1 c2 c3 % for coeffs of end point polyfit curve
syms T_min T_max
syms ind % indicator parameter
A = [0 0 0; -1 0 0; 0 -1 0].*t;
Am = expm(A);
j = 1/2*a*t^2+b*t+r;
J_inner = j^2+1;
J_inner = j^2; % the optimization always choose end point at T_min(?), so we try to change J
% add indicator
% note: this will make dJ=0 has no analytical sol!!!
% So this time we have to use newton iterative method
% ind = 1; % indicator parameter
I1 = -1/ind*log(-T+T_min);
I2 = -1/ind*log(T-T_max);
% J = J+I1+I2;

B = sym(zeros(3,3));
B(1,1) = 720;
B(1,2) = -360*T;
B(1,3) = 60*T^2;
B(2,1) = -360*T;
B(2,2) = 168*T^2;
B(2,3) = -24*T^3;
B(3,1) = 60*T^2;
B(3,2) = -24*T^3;
B(3,3) = 3*T^4;
B = B./T^5;
% set pf vf af as poly of T
pf = c0+c1*T+c2*T^2+c3*T^3;
vf = c1+2*c2*T+3*c3*T^2;
af = 2*c2+6*c3*T;

dp = pf-p0-v0*T-1/2*a0*T^2;
dv = vf-v0-a0*T;
da = af-a0;
Result = B*[dp dv da].';
a = Result(1);
b = Result(2);
r = Result(3);

J = int(J_inner,t,0,T);
J_T = subs(J);
% solve for dJ=0
dJ = diff(J_T,T);
dJ_sim = simplify(dJ);
% dJ_sim = dJ_sim * T^6;

% add dI1/dT & dI2/dT to dJ
dI1dI2 = 1/ind*(2*T-T_max+T_min)/(-T+T_min)/(T-T_max);
dJ_sim = dJ_sim * T^6 * (-T+T_min)*(T-T_max);
dJ_sim = simplify(dJ_sim);
dI1dI2_sim = 1/ind*(2*T-T_max+T_min)*T^6;
dJ_sim = dJ_sim+dI1dI2_sim;

% add I1 & I2 to J
% J = J+I1+I2;

% eqn = dJ == 0;
% S_an = solve(eqn,T);
% S = solve(eqn, T, 'MaxDegree', 4);        
% coeffs = sym2poly(dJ_sim);
[coeffs_dj,T_set] = coeffs(dJ_sim,T);