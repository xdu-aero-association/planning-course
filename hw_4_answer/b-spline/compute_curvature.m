function [kappa] = compute_curvature(dx,d2x,dy,d2y)
% compute curvature
a = dx * d2y - dy * d2x;
norm = sqrt(dx^2 + dy^2);
b = norm * norm^2;
kappa = a / b;
end

