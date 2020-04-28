%% example data
clear
% spline order
k = 4;
% knot sequence
t = [0 0 0 0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1 1 1];
% control points
P = [ 0.1993 0.4965 0.6671 0.7085 0.6809 0.2938 0.1071 0.3929 0.5933 0.8099 0.8998 0.8906 ...
    ; 0.8377 0.8436 0.7617 0.6126 0.212 0.1067 0.3962 0.5249 0.5015 0.3991 0.6477 0.8553 ];
%% convert to own data
degree = k-1;
c_pts = P.';
dense = 30;
knot = t;
n = size(c_pts,1)-1;
ns = n-degree+1; % number of segments
% compute points on curve
pts = deboor_to_points(degree, n, c_pts, knot, dense);
% compute derivative points
degree_deri = degree - 1;
n_deri = n-1;
[d_c_pts, dknot] = get_deri_c_pts(degree, n, c_pts, knot);
deri_pts = deboor_to_points(degree_deri, n_deri, d_c_pts, dknot, dense);
% calculate deri vector
theta_set = zeros(ns * dense, 1);
theta_set = atan2(deri_pts(:,2), deri_pts(:,1));
% theta_set = atan2(pts(:,1), pts(:,2));
%% plot
c_x = c_pts(:,1);
c_y = c_pts(:,2);
figure
plot(c_x,c_y,'*')
axis equal
hold on
plot(pts(:,1),pts(:,2),'-')
title('deboor points with control points')

% plot deri points
figure
plot(deri_pts(:,1), deri_pts(:,2), '.-')
title('deboor derivative points')
axis equal

% plot theta set
figure
plot([1:ns*dense], theta_set*180/pi, '.-')
title('theta in degree')