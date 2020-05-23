close all
clear
% calculate transition traj
px0 = 0;
vx0 = 10;
ax0 = 0;
pxf = 100;
vxf = 20;
axf = 5;
py0 = 0;
vy0 = 0;
ay0 = 0;
pyf = 10;
vyf = 0;
ayf = 0;
start_end_states = [py0 vy0 ay0 pyf vyf ayf];
[optimal_T, optimal_states] = solve_obvp_t(px0,vx0,ax0,pxf,vxf,axf,py0,vy0,ay0,pyf,vyf,ayf);
% plot traj
delta_time = 0.1;
t_size = floor(optimal_T/delta_time)+1;
t_set = linspace(0,optimal_T,t_size);
figure
plot(t_set,optimal_states(:,1),'.-')
title('x = f(t)')
figure
plot(t_set,optimal_states(:,4),'.-')
title('y = f(t)')
figure
plot(optimal_states(:,1),optimal_states(:,4),'.-')
title('final trajectory')
axis equal