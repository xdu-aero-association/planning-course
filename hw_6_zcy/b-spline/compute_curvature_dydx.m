function [kappa] = compute_curvature_dydx(dydx,d2ydx2)
% compute curvature
den = sqrt((1+dydx^2)^3);
kappa = d2ydx2 / den;
end

