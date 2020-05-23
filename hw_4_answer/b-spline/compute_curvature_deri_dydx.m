function [dkappa] = compute_curvature_deri_dydx(dydx,d2ydx2,d3ydx3)
% compute curvature derivative
a = sqrt(1+dydx^2);
dkappa = d3ydx3*a^(-3)-3*dydx*d2ydx2^2*a^(-5);
end

