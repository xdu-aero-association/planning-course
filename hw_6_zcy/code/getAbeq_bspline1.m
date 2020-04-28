function [Aeq, beq] = getAbeq_bspline1(n, p, start_cond, end_cond)
% Note: you don't need continuity constraints because the nature of bspline
% garentees you the continuity if you reuse p control points for each segment
% n: n+1 control points
% p: degree of bspline
% note we only need to optimize middle n+1-2*p control points
% the first & last p control points are determined by start & end
% conditions, so start & end conditions dimension are also p x 1
% so in this case, degree = order+1
% first calculate knots
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
%#####################################################
% STEP 2.1 p,v,a constraint in start & end
beq_start_end = [start_cond;end_cond];
Aeq_start = zeros(p, n+1);
Aeq_end = zeros(p, n+1);

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

Aeq_start(1:p,1:p) = Bs;
Aeq_end(1:p,n-p+2:n+1) = Be;
Aeq = [Aeq_start;Aeq_end];
beq = [beq_start_end];    

end