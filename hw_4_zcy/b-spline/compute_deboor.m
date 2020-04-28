function [pts, deri_pts] = compute_deboor(degree, c_pts, dense, type)
% type: knot type, 0 for start curve, 1 for consequent curve
% compute points on curve
n = size(c_pts,1)-1;
knot = deboor_knot(degree,n, type);
pts = deboor_to_points(degree, n, c_pts, knot, dense);
% compute derivative points
degree_deri = degree - 1;
n_deri = n-1;
[d_c_pts, dknot] = get_deri_c_pts(degree, n, c_pts, knot);
deri_pts = deboor_to_points(degree_deri, n_deri, d_c_pts, dknot, dense);
end

