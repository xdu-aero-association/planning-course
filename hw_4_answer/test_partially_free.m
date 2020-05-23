close all
clear
% partially free OBVP
syms p0 v0 a0 pf vf af
syms t T T0 dp
dp = pf-p0-v0*T-a0*T^2/2;
a = 20/T^5*dp;
j = a/2*(t-T)^2;
J_inner = 1/T*j^2;
J = int(J_inner,t,0,T);
% T0 = T;
J_T = subs(J);
% solve for dJ=0
dJ = diff(J_T,T);
dJ_sim = simplify(dJ);
dJ_sim = expand(dJ_sim);
dJ_sim = dJ_sim*T^7;
dJ_sim = expand(dJ_sim);
dJ_sim = simplify(dJ_sim);
[coeffs_dj,T_set] = coeffs(dJ_sim,T);
eqn = dJ == 0;
S = solve(eqn,T, 'MaxDegree', 4);        