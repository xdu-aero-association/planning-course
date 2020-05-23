close all
clear
% calculate cost function J = int(jerk^2+1)dx|x=[0,xf]
% for 5th degree poly y = f(x)
syms x J xf
syms c0 c1 c2 c3 c4 c5

y = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4 + c5*x^5;
j = diff(y,x,3);
J_inner = j^2+1;
J = int(J_inner,0,xf);
J = simplify(J);
