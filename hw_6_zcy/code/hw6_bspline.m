clc;clear;close all
v_max = 400;
a_max = 400;
color = ['r', 'b', 'm', 'g', 'k', 'c', 'c'];

%% specify the center points of the flight corridor and the region of corridor
c_x = [0 10 20.5 30   40.5 50 60 70 80 90 100]';
c_y = [0 -4 1    6.5  8    10 6  5  10  0  -2]';
c_pts = [c_x c_y];
degree = 4;
dense = 80;
n = max(size(c_pts))-1;
n_seg = n-degree+1;
[pts,theta_set,kappa_set,dkappa_set,d_c_pts,dd_c_pts,ddd_c_pts] = get_deboor_attributes(c_pts,degree,dense,2);
%% optimization process
% add 2*degree-2 control points to make up for start & end cond
init_type = 1;
[c_pts1] = get_extended_c_pts(c_pts,degree,init_type);
% optimization
start_cond_xy = zeros(4,2);
end_cond_xy = zeros(4,2);
for i = 1:2
    start_cond_xy(:,i) = [c_pts(1,i);d_c_pts(1,i);dd_c_pts(1,i);ddd_c_pts(1,i)];
    end_cond_xy(:,i) = [c_pts(end,i);d_c_pts(end,i);dd_c_pts(end,i);ddd_c_pts(end,i)];
end
[c_pts_opt,xy_limits] = GetOptBsplineCpts(c_pts1,degree,start_cond_xy,end_cond_xy,v_max,a_max,init_type);
% trim c_pts_opt according to init_type
[c_pts_opt] = trim_c_pts_opt(c_pts_opt,degree,init_type);
% compute points on curve
n1 = max(size(c_pts_opt))-1;
n_seg1 = n1-degree+1;
dense1 = floor(dense*n_seg/n_seg1);
[pts_opt,theta_set_opt,kappa_set_opt,dkappa_set_opt,d_c_pts_opt,dd_c_pts_opt,ddd_c_pts_opt] = get_deboor_attributes(c_pts_opt,degree,dense1,2);
%% specify ts for each segment
ts = zeros(n_seg, 1);
for i = 1:n_seg
    ts(i,1) = 1;
end
%% display the trajectory and cooridor
figure
plot(pts(:,1),pts(:,2),'.-')
axis equal
hold on
plot(pts_opt(:,1),pts_opt(:,2),'.-')
hold on
plot(c_pts_opt(:,1),c_pts_opt(:,2),'*')
hold on
for i = 1:n+1
    plot_rect([c_pts(i,1);c_pts(i,2)], xy_limits(i,1), xy_limits(i,2));hold on;
end
title('points on curve comparison')
legend('origin points','opt points')

% plot comparison of theta,kappa,dkappa
figure
subplot(1,4,1);
plot(c_x,c_y,'*')
ylim([-10 10])
% axis equal
hold on
plot(pts(:,1),pts(:,2),'.-')
hold on
plot(pts_opt(:,1),pts_opt(:,2),'.-')
title('deboor points with control points')
subplot(1,4,2);
plot(theta_set*180/pi, '.-')
hold on
plot(theta_set_opt*180/pi, '.-')
title('theta in degree')
subplot(1,4,3);
plot(kappa_set, '.-')
hold on
plot(kappa_set_opt,'.-')
title('kappa')
subplot(1,4,4);
plot(dkappa_set, '.-')
hold on
plot(dkappa_set_opt,'.-')
title('dkappa')

function plot_rect(center, x_r, y_r)
    p1 = center+[-x_r;-y_r];
    p2 = center+[-x_r;y_r];
    p3 = center+[x_r;y_r];
    p4 = center+[x_r;-y_r];
    plot_line(p1,p2);
    plot_line(p2,p3);
    plot_line(p3,p4);
    plot_line(p4,p1);
end

function plot_line(p1,p2)
    a = [p1(:),p2(:)];    
    plot(a(1,:),a(2,:),'b');
end