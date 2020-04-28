function [optimal_T, optimal_states] = solve_obvp_t(px0,vx0,ax0,pxf,vxf,axf,py0,vy0,ay0,pyf,vyf,ayf)
% solve OBVP with xf/T as variable
% first solve for xf/T by calculating companion matrix C's eigen values
% C's char. poly. det(sI-C) is poly(s)
% C's rank will always be full, so it's always a square matrix and can have
% non-zero eigen values, because the poly of eigen values are of max
% order(max/min order coeff != 0)
coeffs = zeros(7,1); % [T^0 , ... , T^n];
coeffs(1) = - 3600*px0^2 + 7200*px0*pxf - 3600*pxf^2 - 3600*py0^2 + 7200*py0*pyf - 3600*pyf^2;
coeffs(2) = 2880*pxf*vx0 - 2880*px0*vxf - 2880*px0*vx0 + 2880*pxf*vxf - 2880*py0*vy0 - 2880*py0*vyf + 2880*pyf*vy0 + 2880*pyf*vyf;
coeffs(3) = - 576*vx0^2 - 1008*vx0*vxf - 576*vxf^2 - 576*vy0^2 - 1008*vy0*vyf - 576*vyf^2 - 360*ax0*px0 + 360*ax0*pxf + 360*axf*px0 - 360*axf*pxf - 360*ay0*py0 + 360*ay0*pyf + 360*ayf*py0 - 360*ayf*pyf;
coeffs(4) = 96*axf*vx0 - 96*ax0*vxf - 144*ax0*vx0 + 144*axf*vxf - 144*ay0*vy0 - 96*ay0*vyf + 96*ayf*vy0 + 144*ayf*vyf;
coeffs(5) = - 9*ax0^2 + 6*ax0*axf - 9*axf^2 - 9*ay0^2 + 6*ay0*ayf - 9*ayf^2;
coeffs(6) = 0;
coeffs(7) = 1;
coeffs_matlab = flip(coeffs);
assignin('base','coeffs_matlab',coeffs_matlab)
roots_matlab = roots(coeffs_matlab);
assignin('base','roots_matlab',roots_matlab)
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
    end
    T = real(root_set(i));
    % calculate a b r
    [a1,b1,r1,a2,b2,r2] = calculate_abr(T,px0,vx0,ax0,pxf,vxf,axf,py0,vy0,ay0,pyf,vyf,ayf);
    % calculate cost J
    J = T^3*(b1^2/3 + b2^2/3 + (a1*r1)/3 + (a2*r2)/3) + T^5*(a1^2/20 + a2^2/20) + T*(r1^2 + r2^2 + 1) + T^4*((a1*b1)/4 + (a2*b2)/4) + T^2*(b1*r1 + b2*r2);
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
[optimal_states] = reconstructOptimalTrajactory(optimal_T,px0,vx0,ax0,pxf,vxf,axf,py0,vy0,ay0,pyf,vyf,ayf);

end

function [a1,b1,r1,a2,b2,r2] = calculate_abr(T,px0,vx0,ax0,pxf,vxf,axf,py0,vy0,ay0,pyf,vyf,ayf)
% calculate a b r
B = sym(zeros(6,6));
B(1,1) = 1/120*T^5;
B(1,3) = 1/24*T^4;
B(1,5) = 1/6*T^3;
B(2,2:6) = B(1,1:5);
B(3,1) = 1/24*T^4;
B(3,3) = 1/6*T^3;
B(3,5) = 1/2*T^2;
B(4,2:6) = B(3,1:5);
B(5,1) = 1/6*T^3;
B(5,3) = 1/2*T^2;
B(5,5) = T;
B(6,2:6) = B(5,1:5);
dpx = pxf-px0-vx0*T-1/2*ax0*T^2;
dvx = vxf-vx0-ax0*T;
dax = axf-ax0;
dpy = pyf-py0-vy0*T-1/2*ay0*T^2;
dvy = vyf-vy0-ay0*T;
day = ayf-ay0;
Result = B^(-1)*[dpx dpy dvx dvy dax day].';
a1 = Result(1);
a2 = Result(2);
b1 = Result(3);
b2 = Result(4);
r1 = Result(5);
r2 = Result(6);
end

function [optimal_states] = reconstructOptimalTrajactory(optimal_T,px0,vx0,ax0,pxf,vxf,axf,py0,vy0,ay0,pyf,vyf,ayf)
    fprintf(1,'optimal_T: %.3f\n',optimal_T)
    delta_time = 0.05;
    t_size = floor(optimal_T/delta_time)+1;
    t_set = linspace(0,optimal_T,t_size);
    optimal_states = zeros(t_size,6);
    % calculate a b r
    [a1,b1,r1,a2,b2,r2] = calculate_abr(optimal_T,px0,vx0,ax0,pxf,vxf,axf,py0,vy0,ay0,pyf,vyf,ayf);
    for i = 1:t_size
        t = t_set(i);
        optimal_states(i,1) = a1/120*t^5+b1/24*t^4+r1/6*t^3+ax0/2*t^2+vx0*t+px0;
        optimal_states(i,2) = a1/24*t^4+b1/6*t^3+r1/2*t^2+ax0*t+vx0;
        optimal_states(i,3) = a1/6*t^3+b1/2*t^2+r1*t+ax0;
        optimal_states(i,4) = a2/120*t^5+b2/24*t^4+r2/6*t^3+ay0/2*t^2+vy0*t+py0;
        optimal_states(i,5) = a2/24*t^4+b2/6*t^3+r2/2*t^2+ay0*t+vy0;
        optimal_states(i,6) = a2/6*t^3+b2/2*t^2+r2*t+ay0;
    end
end
