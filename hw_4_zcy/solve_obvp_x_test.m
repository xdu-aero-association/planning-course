close all
clear
% calculate transition traj
py0 = 0;
vy0 = 0;
ay0 = 0;
pyf = 100;
vyf = 1;
ayf = 0;
start_end_states = [py0 vy0 ay0 pyf vyf ayf];
[optimal_T, optimal_states] = solve_obvp_x(py0,vy0,ay0,pyf,vyf,ayf);
% plot traj
delta_time = 0.1;
t_size = floor(optimal_T/delta_time)+1;
t_set = linspace(0,optimal_T,t_size);
figure
plot(t_set,optimal_states(:,1),'.-')
title('y = f(x) traj')