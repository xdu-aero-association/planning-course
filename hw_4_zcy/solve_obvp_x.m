function [optimal_T, optimal_states] = solve_obvp_x(py0,vy0,ay0,polyfit_pts)
% solve OBVP with xf/T as variable
% first solve for xf/T by calculating companion matrix C's eigen values
% C's char. poly. det(sI-C) is poly(s)
% C's rank can be decreased by initial & final conditions, so this method
% is not valid for all circumstances
% Update: after changing the cost function to int(1+j^2)dt, we can always
% have full rank C
p0 = py0;
v0 = vy0;
a0 = ay0;
coeffs_polyfit = polyfit(polyfit_pts(:,1),polyfit_pts(:,2),3);
% plot polyfit result
figure
f = polyval(coeffs_polyfit,polyfit_pts(:,1));
plot(polyfit_pts(:,1),polyfit_pts(:,2),'o',polyfit_pts(:,1),f,'-')
legend('data','linear fit')
coeffs_polyfit = flip(coeffs_polyfit);
c0 = coeffs_polyfit(1);
c1 = coeffs_polyfit(2);
c2 = coeffs_polyfit(3);
c3 = coeffs_polyfit(4);
ind = 16; % indcator param, larger the better approximation of indicator to real one
T_min = polyfit_pts(1,1);
T_max = polyfit_pts(end,1);
T_extend_range = (T_max-T_min)/2;
T_min = T_min + T_extend_range/2;
T_max = T_max - T_extend_range/2;
fprintf(1,'T_min: %.3f, T_max: %.3f\n',T_min,T_max)

coeffs = zeros(9,1); % [T^0 , ... , T^n];
coeff_set_select = 0;
if coeff_set_select == 0
    % set 1: for J = int(1+j^2,dt)
    coeffs(1) = T_min*T_max*(3600*c0^2 - 7200*c0*p0 + 3600*p0^2);
    coeffs(2) = T_min*T_max*(2880*c0*c1 - 2880*c1*p0 - 2880*c0*v0 + 2880*p0*v0) - (T_min + T_max)*(3600*c0^2 - 7200*c0*p0 + 3600*p0^2);
    coeffs(3) = 3600*c0^2 - (T_min + T_max)*(2880*c0*c1 - 2880*c1*p0 - 2880*c0*v0 + 2880*p0*v0) - 7200*c0*p0 + 3600*p0^2 + T_min*T_max*(576*c1^2 - 1152*c1*v0 + 576*v0^2 - 360*a0*c0 + 720*c0*c2 + 360*a0*p0 - 720*c2*p0);
    coeffs(4) = 2880*c0*c1 - (T_min + T_max)*(576*c1^2 - 1152*c1*v0 + 576*v0^2 - 360*a0*c0 + 720*c0*c2 + 360*a0*p0 - 720*c2*p0) - 2880*c1*p0 - 2880*c0*v0 + 2880*p0*v0 - T_min*T_max*(144*a0*c1 - 288*c1*c2 - 144*a0*v0 + 288*c2*v0);
    coeffs(5) = 720*c0*c2 - 360*a0*c0 + 360*a0*p0 - 720*c2*p0 - 1152*c1*v0 + (T_min + T_max)*(144*a0*c1 - 288*c1*c2 - 144*a0*v0 + 288*c2*v0) + 576*c1^2 + 576*v0^2 + T_min*T_max*(9*a0^2 - 36*a0*c2 + 36*c2^2);
    coeffs(6) = 288*c1*c2 - 144*a0*c1 - (T_min + T_max)*(9*a0^2 - 36*a0*c2 + 36*c2^2) + 144*a0*v0 - 288*c2*v0;
    coeffs(7) = 9*a0^2 - 36*a0*c2 + 36*c2^2 + (T_min - T_max)/ind - T_min*T_max*(36*c3^2 + 1);
    coeffs(8) = 2/ind + (T_min + T_max)*(36*c3^2 + 1);
    coeffs(9) = - 36*c3^2 - 1;
else
    % set 2: for J = int(j^2,dt)
    coeffs(1) = T_min*T_max*(3600*c0^2 - 7200*c0*p0 + 3600*p0^2);
    coeffs(2) = T_min*T_max*(2880*c0*c1 - 2880*c1*p0 - 2880*c0*v0 + 2880*p0*v0) - (T_min + T_max)*(3600*c0^2 - 7200*c0*p0 + 3600*p0^2);
    coeffs(3) = 3600*c0^2 - (T_min + T_max)*(2880*c0*c1 - 2880*c1*p0 - 2880*c0*v0 + 2880*p0*v0) - 7200*c0*p0 + 3600*p0^2 + T_min*T_max*(576*c1^2 - 1152*c1*v0 + 576*v0^2 - 360*a0*c0 + 720*c0*c2 + 360*a0*p0 - 720*c2*p0);
    coeffs(4) = 2880*c0*c1 - (T_min + T_max)*(576*c1^2 - 1152*c1*v0 + 576*v0^2 - 360*a0*c0 + 720*c0*c2 + 360*a0*p0 - 720*c2*p0) - 2880*c1*p0 - 2880*c0*v0 + 2880*p0*v0 - T_min*T_max*(144*a0*c1 - 288*c1*c2 - 144*a0*v0 + 288*c2*v0);
    coeffs(5) = 720*c0*c2 - 360*a0*c0 + 360*a0*p0 - 720*c2*p0 - 1152*c1*v0 + (T_min + T_max)*(144*a0*c1 - 288*c1*c2 - 144*a0*v0 + 288*c2*v0) + 576*c1^2 + 576*v0^2 + T_min*T_max*(9*a0^2 - 36*a0*c2 + 36*c2^2);
    coeffs(6) = 288*c1*c2 - 144*a0*c1 - (T_min + T_max)*(9*a0^2 - 36*a0*c2 + 36*c2^2) + 144*a0*v0 - 288*c2*v0;
    coeffs(7) = 9*a0^2 - 36*a0*c2 + 36*c2^2 - 36*T_min*T_max*c3^2 + (T_min - T_max)/ind;
    coeffs(8) = 36*c3^2*(T_min + T_max) + 2/ind;
    coeffs(9) = -36*c3^2;
end
coeffs_matlab = flip(coeffs);
assignin('base','coeffs_matlab',coeffs_matlab)
roots_matlab = roots(coeffs_matlab);
assignin('base','roots_matlab',roots_matlab)
coeffs = coeffs./coeffs(max(size(coeffs))); % we need to normalize coeff(n_order)
flag = 1;
if flag == 1
    n = max(size(coeffs))-1;
    C = diag(ones(n-1,1),-1);
    for i = 1:n
        C(i,n) = -coeffs(i);
    end
    fprintf(1,'rank of C: %d\n',rank(C));
    root_set = eig(C); % possible values of xf/T
    assignin('base','roots',root_set)
else
    root_set = roots_matlab;
end
% choose root with smallest cost J
first_flag = 0;
optimal_T = 0;
optimal_J = 0;
for i = 1:max(size(root_set))
    if abs(imag(root_set(i))) > 1e-10 || real(root_set(i)) <= 0
        fprintf(1,' % .3f + %.3fi not a valid root !!!\n',real(root_set(i)),imag(root_set(i)));
        continue;
    elseif real(root_set(i)) > T_max+T_extend_range || real(root_set(i)) < T_min-T_extend_range
        fprintf(1,' % .3f + %.3fi not a valid root(not within limits) !!!\n',real(root_set(i)),imag(root_set(i)));
        continue;
    end
    T = real(root_set(i));
    fprintf(1,' %.3f is a valid root\n',T);
    if T < T_min 
        T = T_min;
    elseif T > T_max
        T = T_max;
    end
    fprintf(1,' %.3f is T after bounding within limits\n',T);
    % calculate a b r
    [a,b,r] = calculate_abr(T,p0,v0,a0,c0,c1,c2,c3);
    % calculate cost J
    if coeff_set_select == 0
        % set 1: for J = int(1+j^2,dt)
        J = T^3*(b^2/3 + (a*r)/3) + (T^5*a^2)/20 + T*(r^2 + 1) + (T^4*a*b)/4 + T^2*b*r;
    else
        % set 2: for J = int(j^2,dt)
        J = T*r^2 + T^3*(b^2/3 + (a*r)/3) + (T^5*a^2)/20 + (T^4*a*b)/4 + T^2*b*r;
    end
    % indicator for T_min/T_max
    I1 = -1/ind*log(-T+T_min);
    I2 = -1/ind*log(T-T_max);
    J = J+I1+I2;
    % select T with min J
    if first_flag == 0
        first_flag = 1;
        optimal_T = T;
        optimal_J = J;
    elseif optimal_J > J
        optimal_T = T;
        optimal_J = J;
    end
end
% reconstruct optimal states from optimal xf/T
[optimal_states] = reconstructOptimalTrajactory(optimal_T,py0,vy0,ay0,c0,c1,c2,c3);

end

function [a,b,r] = calculate_abr(T,p0,v0,a0,c0,c1,c2,c3)
% calculate a b r
B = zeros(3,3);
B(1,1) = 720;
B(1,2) = -360*T;
B(1,3) = 60*T^2;
B(2,1) = -360*T;
B(2,2) = 168*T^2;
B(2,3) = -24*T^3;
B(3,1) = 60*T^2;
B(3,2) = -24*T^3;
B(3,3) = 3*T^4;
B = B/T^5;
% set pf vf af as poly of T
pf = c0+c1*T+c2*T^2+c3*T^3;
vf = c1+2*c2*T+3*c3*T^2;
af = 2*c2+6*c3*T;

dp = pf-p0-v0*T-1/2*a0*T^2;
dv = vf-v0-a0*T;
da = af-a0;
Result = B*[dp dv da].';
a = Result(1);
b = Result(2);
r = Result(3);
end

function [optimal_states] = reconstructOptimalTrajactory(optimal_T,p0,v0,a0,c0,c1,c2,c3)
    fprintf(1,'optimal_T: %.3f\n',optimal_T)
    delta_time = 0.05;
    t_size = floor(optimal_T/delta_time)+1;
    t_set = linspace(0,optimal_T,t_size);
    optimal_states = zeros(t_size,3);
    [a,b,r] = calculate_abr(optimal_T,p0,v0,a0,c0,c1,c2,c3);
    for i = 1:t_size
        t = t_set(i);
        optimal_states(i,1) = a/120*t^5+b/24*t^4+r/6*t^3+a0/2*t^2+v0*t+p0;
        optimal_states(i,2) = a/24*t^4+b/6*t^3+r/2*t^2+a0*t+v0;
        optimal_states(i,3) = a/6*t^3+b/2*t^2+r*t+a0;
    end
end
