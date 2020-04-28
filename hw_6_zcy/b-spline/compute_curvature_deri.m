function [dkappa] = compute_curvature_deri(dx,d2x,d3x,dy,d2y,d3y)
% compute curvature derivative
a = dx * d2y - dy * d2x;
b = dx * d3y - dy * d3x;
c = dx * d2x + dy * d2y;
d = dx^2 + dy^2;
dkappa = (b * d - 3 * a * c) / d^3;
end

