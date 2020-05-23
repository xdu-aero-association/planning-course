function [C] = getC(n,p,knot)
% C is bspline control point deri mtx
% n+1 is number of control points P
% p is degree of bspline
% knot is origin knot of bspline
% Q is set of 1st deri control points of P(the original control points)
% size(Q) = n = size(P)-1
% Q[i] = p/(knot[i+p+1]-knot[i+1]) * (P[i+1]-P[i]), i in [0,n-1]
% Q = C * P
C = zeros(n,n+1);
for i = 1:n
    a = -p/(knot(i+p+1)-knot(i+1));
    C(i,i) = a;
    C(i,i+1) = -a;
end

end

