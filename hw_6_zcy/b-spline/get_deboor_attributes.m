function [pts,theta_set,kappa_set,dkappa_set,d_c_pts,dd_c_pts,ddd_c_pts] = get_deboor_attributes(c_pts,degree,dense,knot_type)
n = max(size(c_pts))-1;
ns = n-degree+1; % number of segments
% compute knots
knots = deboor_knot(degree, n, knot_type);
algorithm_type = 1; % 0 for deboor algorithm, 1 for standard def.
% compute points on curve
pts = deboor_to_points(degree, n, c_pts, knots, dense, algorithm_type);
% compute derivative points
degree_d = degree - 1;
n_d = n - 1;
[d_c_pts, dknots] = get_deri_c_pts(degree, n, c_pts, knots);
d_pts = deboor_to_points(degree_d, n_d, d_c_pts, dknots, dense, algorithm_type);
assignin('base','d_pts',d_pts)
% compute 2nd deri points
degree_dd = degree_d - 1;
n_dd = n_d - 1;
[dd_c_pts, ddknots] = get_deri_c_pts(degree_d, n_d, d_c_pts, dknots);
dd_pts = deboor_to_points(degree_dd, n_dd, dd_c_pts, ddknots, dense, algorithm_type);
assignin('base','dd_pts',dd_pts)
% compute 3rd deri points
degree_ddd = degree_dd - 1;
n_ddd = n_dd - 1;
[ddd_c_pts, dddknot] = get_deri_c_pts(degree_dd, n_dd, dd_c_pts, ddknots);
ddd_pts = deboor_to_points(degree_ddd, n_ddd, ddd_c_pts, dddknot, dense, algorithm_type);
assignin('base','ddd_pts',ddd_pts)
% calculate theta
theta_set = zeros(ns * dense, 1);
theta_set = atan2(d_pts(:,2), d_pts(:,1));
% theta_set = atan2(pts(:,1), pts(:,2));
% calculate curvature
kappa_set = zeros(ns * dense, 1);
for i = 1:ns*dense
    kappa_set(i) = compute_curvature(d_pts(i,1),dd_pts(i,1),d_pts(i,2),dd_pts(i,2));
end
% calculate curvature deri
dkappa_set = zeros(ns * dense, 1);
for i = 1:ns*dense
    dkappa_set(i) = compute_curvature_deri(d_pts(i,1),dd_pts(i,1),ddd_pts(i,1),d_pts(i,2),dd_pts(i,2),ddd_pts(i,2));
end

end

