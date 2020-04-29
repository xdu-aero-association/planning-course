syms a11 a12 a21 a22
syms y1 y2 a b c
A = [a b;b c];
y = [y1;y2];
M = A*(y*y.')*A^(-1);
N = y*y.';