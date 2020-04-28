clear
close all
%% prepare data
c_x = [0 10 20.5 30   40.5 50 60 70 80 90 100]';
c_y = [0 -4 1    6.5  8    10 6  5  10  0  -2]';
c_pts = [c_x c_y];
degree = 4;
dense = 80;
n = size(c_pts,1)-1;
ns = n-degree+1; % number of segments
% compute knots
knots = deboor_knot(degree, n, 0);
% compute points on curve
pts = deboor_to_points(degree, n, c_pts, knots, dense, 0);
% compute derivative points
degree_d = degree - 1;
n_d = n - 1;
[d_c_pts, dknots] = get_deri_c_pts(degree, n, c_pts, knots);
d_pts = deboor_to_points(degree_d, n_d, d_c_pts, dknots, dense, 0);
% compute 2nd deri points
degree_dd = degree_d - 1;
n_dd = n_d - 1;
[dd_c_pts, ddknots] = get_deri_c_pts(degree_d, n_d, d_c_pts, dknots);
dd_pts = deboor_to_points(degree_dd, n_dd, dd_c_pts, ddknots, dense, 0);
% compute 3rd deri points
degree_ddd = degree_dd - 1;
n_ddd = n_dd - 1;
[ddd_c_pts, dddknot] = get_deri_c_pts(degree_dd, n_dd, dd_c_pts, ddknots);
ddd_pts = deboor_to_points(degree_ddd, n_ddd, ddd_c_pts, dddknot, dense, 0);
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
% use two method to calculate transition path
flag_method = 2;
if flag_method == 0
%% generate transition path by y = f(x)
    % first get start & end points
    start_idx = 15;
    end_idx = 500;
    xs = pts(start_idx,1);
    ys = pts(start_idx,2);
    thetas = theta_set(start_idx);
    ks = kappa_set(start_idx);
    dks = dkappa_set(start_idx);
    xf = pts(end_idx,1);
    yf = pts(end_idx,2);
    thetaf = theta_set(end_idx);
    kf = kappa_set(end_idx);
    dkf = dkappa_set(end_idx);
    % select a region around end point to use for polyfit
    idx_range = 20;
    range_start = end_idx-idx_range;
    range_end = end_idx+idx_range;
    if range_end > size(pts,1)
        range_end = size(pts,1);
        range_start = range_end - 2*idx_range;
    elseif range_start < 1
        range_start = 1;
        range_end = range_start + 2*idx_range; 
    end 
    polyfit_pts = pts(range_start:range_end,:);
    % transform end point & polyfit pts to start point coord frame
    [xfl,yfl] = global2local(xs,ys,thetas,xf,yf);
    thetafl = thetaf - thetas;
    for i = 1:size(polyfit_pts,1)
        xg = polyfit_pts(i,1);
        yg = polyfit_pts(i,2);
        [xl,yl] = global2local(xs,ys,thetas,xg,yg);
        polyfit_pts(i,:) = [xl,yl];
    end
    % calculate transition traj
    py0 = 0;
    vy0 = 0;
    ay0 = ks;
    pyf = yfl;
    vyf = atan(thetafl);
    ayf = kf*(1+vyf^2)^(1.5);
    start_end_states = [py0 vy0 ay0 pyf vyf ayf];
    [optimal_T, optimal_states] = solve_obvp_x(py0,vy0,ay0,polyfit_pts);
    % transform transition pathpoints from local to global
%     t_set = linspace(0,xfl,size(optimal_states,1));
    delta_time = 0.05;
    t_size = floor(optimal_T/delta_time)+1;
    t_set = linspace(0,optimal_T,t_size);
    trans_path = zeros(size(optimal_states,1),2);
    for i = 1:size(optimal_states,1)
        xl = t_set(i);
        yl = optimal_states(i,1);
        [xg,yg] = local2global(xs,ys,thetas,xl,yl);
        trans_path(i,:) = [xg,yg];
    end
    % get real end point idx
    real_end_idx = 1;
    min_dist = sqrt((pts(1,1)-xg)^2+(pts(1,2)-yg)^2);
    for i = 2:size(pts,1)
        dist = sqrt((pts(i,1)-xg)^2+(pts(i,2)-yg)^2);
        if dist < min_dist
            real_end_idx = i;
            min_dist = dist;
        end
    end
    % calculate theta kappa dkappa for transition path
    theta_tran_set = atan(optimal_states(:,2));
    theta_tran_set = theta_tran_set + thetas;
    theta_tran_set = [theta_set(1:start_idx);theta_tran_set;theta_set(real_end_idx:end)];
    kappa_tran_set = zeros(size(optimal_states,1),1);
    for i = 1:size(optimal_states,1)
        dydx = optimal_states(i,2);
        d2ydx2 = optimal_states(i,3);
        kappa_tran_set(i) = compute_curvature_dydx(dydx,d2ydx2);
    end
    kappa_tran_set = [kappa_set(1:start_idx);kappa_tran_set;kappa_set(real_end_idx:end)];
    t_set = [1:size(theta_tran_set,1)];
elseif flag_method == 1
%% generate transition path by x = f(t) & y = f(t)
    % first get start & end points
    start_idx = 15;
    end_idx = 500;
    xs = pts(start_idx,1);
    ys = pts(start_idx,2);
    thetas = theta_set(start_idx);
    ks = kappa_set(start_idx);
    dks = dkappa_set(start_idx);
    xf = pts(end_idx,1);
    yf = pts(end_idx,2);
    thetaf = theta_set(end_idx);
    kf = kappa_set(end_idx);
    dkf = dkappa_set(end_idx);
    % calculate transition traj
    px0 = xs;
    vx0 = 10; % current speed in global x direction
    ax0 = 1; % current acc in global x direction
    pxf = xf;
    vxf = 8;
    axf = 0;
    py0 = ys;
    vy0 = vx0*tan(thetas);
    ay0 = ks*vx0^2*(1+(tan(thetas))^2)^1.5+ax0*tan(thetas);
    pyf = yf;
    vyf = vxf*tan(thetaf);
    ayf = kf*vxf^2*(1+(tan(thetaf))^2)^1.5+axf*tan(thetaf);
    [optimal_T, optimal_states] = solve_obvp_t(px0,vx0,ax0,pxf,vxf,axf,py0,vy0,ay0,pyf,vyf,ayf);
    trans_path = zeros(size(optimal_states,1),2);
    trans_path(:,1) = optimal_states(:,1);
    trans_path(:,2) = optimal_states(:,4);
    % calculate theta kappa dkappa for transition path
%     delta_time = 0.05;
%     t_size = floor(optimal_T/delta_time)+1;
%     t_set = linspace(0,optimal_T,t_size);
    theta_tran_set = atan2(optimal_states(:,5),optimal_states(:,2));
    theta_tran_set = [theta_set(1:start_idx);theta_tran_set;theta_set(end_idx:end)];
    kappa_tran_set = zeros(size(optimal_states,1),1);
    for i = 1:size(optimal_states,1)
        dx = optimal_states(i,2);
        d2x = optimal_states(i,3);
        dy = optimal_states(i,5);
        d2y = optimal_states(i,6);
        kappa_tran_set(i) = compute_curvature(dx,d2x,dy,d2y);
    end
    kappa_tran_set = [kappa_set(1:start_idx);kappa_tran_set;kappa_set(end_idx:end)];
    t_set = [1:size(theta_tran_set,1)];
else
%% generate transition path by 7th order y=f(x) for dkappa continuiity at start & end
% first get start & end points
    start_idx = 15;
    end_idx = 500;
    xs = pts(start_idx,1);
    ys = pts(start_idx,2);
    thetas = theta_set(start_idx);
    ks = kappa_set(start_idx);
    dks = dkappa_set(start_idx);
    xf = pts(end_idx,1);
    yf = pts(end_idx,2);
    thetaf = theta_set(end_idx);
    kf = kappa_set(end_idx);
    dkf = dkappa_set(end_idx);
    % transform end point fo start point coord frame
    [xfl,yfl] = global2local(xs,ys,thetas,xf,yf);
    thetafl = thetaf - thetas;
    % calculate transition traj
    py0 = 0;
    vy0 = 0;
    ay0 = ks*(1+vy0^2)^(1.5);
    jy0 = (dks*(1+vy0^2)^3+3*vy0*ay0^2*sqrt(1+vy0^2))/(sqrt(1+vy0^2))^3;
    pyf = yfl;
    vyf = tan(thetafl);
    ayf = kf*(1+vyf^2)^(1.5);    
    jyf = (dkf*(1+vyf^2)^3+3*vyf*ayf^2*sqrt(1+vyf^2))/(sqrt(1+vyf^2))^3;

    order = 7;
    P = zeros(order+1,order+1);
    % start cond & end con
    for i = 0:3
        P(i+1,i+1) = factorial(i);
        for j = i:order
           P(i+1+4,j+1) = factorial(j)/factorial(j-i)*xfl^(j-i);
        end
    end
    % solve for y = f(x) coeffs
%     coeffs_c = (P.'*P)^(-1)*P.' * [py0 vy0 ay0 jy0 pyf vyf ayf jyf].';
    coeffs_c = P^(-1) * [py0 vy0 ay0 jy0 pyf vyf ayf jyf].';
    delta_x = 0.05;
    x_size = floor(xfl/delta_x)+1;
    x_set = linspace(0,xfl,x_size);
    trans_path = zeros(x_size,2);
    for i = 1:x_size
        xl = x_set(i);
        yl = polyval(flip(coeffs_c),xl);
        [xg,yg] = local2global(xs,ys,thetas,xl,yl);
        trans_path(i,:) = [xg,yg];
    end
    Q = zeros(4,order+1);
    for i = 0:3
        for j = i:order
           Q(i+1,j+1) = factorial(j)/factorial(j-i)*coeffs_c(j+1);
        end
    end
    theta_tran_set = zeros(x_size,1);
    kappa_tran_set = zeros(x_size,1);
    dkappa_tran_set = zeros(x_size,1);
    for i = 1:x_size
        x_temp = x_set(i);
        x_poly = zeros(8,4);
        for k = 1:4
            for j = k:8
                x_poly(j,k) = x_temp^(j-k);
            end
        end
        result = Q * x_poly;
        dydx = result(2,2);
        d2ydx2 = result(3,3);
        d3ydx3 = result(4,4);
        theta_tran_set(i) = atan(dydx);
        kappa_tran_set(i) = compute_curvature_dydx(dydx,d2ydx2);
        dkappa_tran_set(i) = compute_curvature_deri_dydx(dydx,d2ydx2,d3ydx3);
    end
    theta_tran_set = theta_tran_set + thetas;
    theta_tran_set = [theta_set(1:start_idx);theta_tran_set;theta_set(end_idx:end)];
    kappa_tran_set = [kappa_set(1:start_idx);kappa_tran_set;kappa_set(end_idx:end)];
    dkappa_tran_set = [dkappa_set(1:start_idx);dkappa_tran_set;dkappa_set(end_idx:end)];
    t_set = [1:size(theta_tran_set,1)];
    % calculate cost function J = int(jerk^2+1)dx|x=[0,xf]
    % for 7th degree poly y = f(x)
    c0 = coeffs_c(1);c1 = coeffs_c(2);c2 = coeffs_c(3);c3 = coeffs_c(4);c4 = coeffs_c(5);c5 = coeffs_c(6);c6 = coeffs_c(7);c7 = coeffs_c(8);
    J = xf^5*(720*c5^2 + 504*c3*c7 + 1152*c4*c6) + xf*(36*c3^2 + 1) + xf^3*(192*c4^2 + 240*c3*c5) + xf^7*((14400*c6^2)/7 + 3600*c5*c7) + xf^4*(360*c3*c6 + 720*c4*c5) + xf^6*(1680*c4*c7 + 2400*c5*c6) + 4900*c7^2*xf^9 + 144*c3*c4*xf^2 + 6300*c6*c7*xf^8;
end
%% plot
c_x = c_pts(:,1);
c_y = c_pts(:,2);
figure
plot(c_x,c_y,'*')
axis equal
hold on
plot(pts(:,1),pts(:,2),'.-')
title('deboor points with control points')
hold on
plot(trans_path(:,1),trans_path(:,2),'.-') % transition path
hold on
plot(trans_path(1,1),trans_path(1,2),'sq') % transition start
hold on
plot(trans_path(end,1),trans_path(end,2),'o') % transition end
hold on
plot(xf,yf,'o') % real transition end

% % plot deri points
% figure
% plot(d_pts(:,1), d_pts(:,2), '.-')
% title('deboor derivative points')
% axis equal
% 
% % plot theta set
% figure
% plot([1:ns*dense], theta_set*180/pi, '.-')
% title('theta in degree')
% 
% % plot curvature
% figure
% plot([1:ns*dense], kappa_set, '.-')
% title('kappa')
% 
% % plot curvature deri
% figure
% plot([1:ns*dense], dkappa_set, '.-')
% title('dkappa')

% plot in one
figure
subplot(1,4,1);
plot(c_x,c_y,'*')
ylim([-10 10])
% axis equal
hold on
plot(pts(:,1),pts(:,2),'.-')
title('deboor points with control points')
subplot(1,4,2);
plot([1:ns*dense], theta_set*180/pi, '.-')
title('theta in degree')
subplot(1,4,3);
plot([1:ns*dense], kappa_set, '.-')
title('kappa')
subplot(1,4,4);
plot([1:ns*dense], dkappa_set, '.-')
title('dkappa')

% transition theta
figure
title('trans path theta kappa')
subplot(1,2,1);
plot(t_set, theta_tran_set/pi*180, '.-')
title('theta_tran_set in degree')
subplot(1,2,2);
plot(t_set, kappa_tran_set, '.-')
title('kappa_tran_set')
if exist('dkappa_tran_set','var')
    figure
    plot(t_set, dkappa_tran_set, '.-')
    title('dkappa_tran_set')
end

function [xl,yl] = global2local(dx,dy,dtheta,xg,yg)
    % B: local, A: global
    ATB = [cos(dtheta) -sin(dtheta) dx;sin(dtheta) cos(dtheta) dy;0 0 1];
    Ap = [xg yg 1].';
    Bp = ATB^-1 * Ap;
    xl = Bp(1);
    yl = Bp(2);
end

function [xg,yg] = local2global(dx,dy,dtheta,xl,yl)
    % B: local, A: global
    ATB = [cos(dtheta) -sin(dtheta) dx;sin(dtheta) cos(dtheta) dy;0 0 1];
    Bp = [xl yl 1].';
    Ap = ATB * Bp;
    xg = Ap(1);
    yg = Ap(2);
end