close all
% clear
c_x = [0 10 20.5 30   40.5 50 60 70 80 90 100]';
c_y = [0 -4 1    6.5  8    10 6  5  10  0  -2]';
c_pts = [c_x c_y];
p = 3;
n = max(size(c_pts))-1;
knots = deboor_knot(p,n,2); % we need each curve seg to be clamped at both ends
dknots = knots(2:max(size(knots))-1);
ddknots = dknots(2:max(size(dknots))-1);
dp = p-1;
ddp = dp-1;
dn = n-1;
ddn = dn-1;
C = getC(n,p,knots);
dC = getC(dn,dp,dknots);
ddC = getC(ddn,ddp,ddknots);
assignin('base','C',C)
assignin('base','dC',dC)
assignin('base','ddC',ddC)
d_c_pts1 = C * c_pts;
dd_c_pts1 = dC * C * c_pts;
ddd_c_pts1 = ddC * dC * C * c_pts;

C_set = zeros((p-1)*n,n+1);
knots_temp = deboor_knot(p,n,2);
C_current = eye(n+1,n+1);
Bs = zeros(p,p);
Bs(1,1) = 1;
Be = zeros(p,p);
Be(1,p) = 1;
for i = 1:p-1
    p_temp = p-i+1;
    n_temp = n-i+1;
    C_temp = getC(n_temp,p_temp,knots_temp);
    C_current = C_temp * C_current;
    indices_row = 1+(i-1)*n:size(C_current,1)+(i-1)*n;
    indices_coln = 1:size(C_current,2);
    C_set(indices_row,indices_coln) = C_current;
    knots_temp = knots_temp(2:max(size(knots_temp))-1);
    Bs(i+1,1:p) = C_current(1,1:p);
    Be(i+1,1:p) = C_current(size(C_current,1),n-p+2:n+1);
end