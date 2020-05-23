function [d_c_pts, dknot] = get_deri_c_pts(degree, n, c_pts, knot)
% calculate derivative control points of original b-spline
% degree: degree of b-spline
% n for n+1 c_pts
% c_pts: control points
% knot: array of knot positions

p = degree;
% derivative knots
dknot = knot(2:max(size(knot))-1);
% derivative control points
d_c_pts = zeros(n,2);
for i = 1 : n
    d_c_pts(i,:) = (p / (knot(i+p+1) - knot(i+1))) * (c_pts(i+1,:) - c_pts(i,:));
end

end

