function [Aeq, beq] = getAbeq_bspline(n_seg, n, p, ts, start_cond, end_cond)
% Note: you don't need continuity constraints because the nature of bspline
% garentees you the continuity if you reuse p control points for each segment
% first calculate knots
knots = deboor_knot(p,n,2); % we need each curve seg to be clamped at both ends
dknots = knots(2:max(size(knots))-1);
ddknots = dknots(2:max(size(dknots))-1);
dp = p-1;
ddp = dp-1;
n_all_poly = n_seg*(n+1);
d_order = 3;
n_coef1d = n+1;
n_coef1d_d = n_coef1d-1;
n_coef1d_d2 = n_coef1d_d-1;
% compute bspline basis Nip
Ns = compute_bspline_basis(p,n_coef1d-1,knots,knots(1));
Ne = compute_bspline_basis(p,n_coef1d-1,knots,knots(end));
dNs = compute_bspline_basis(dp,n_coef1d_d-1,dknots,dknots(1));
dNe = compute_bspline_basis(dp,n_coef1d_d-1,dknots,dknots(end));
ddNs = compute_bspline_basis(ddp,n_coef1d_d2-1,ddknots,ddknots(1));
ddNe = compute_bspline_basis(ddp,n_coef1d_d2-1,ddknots,ddknots(end));
N_set = [Ns Ne dNs dNe ddNs ddNe];
%#####################################################
% STEP 2.1 p,v,a constraint in start & end
Aeq_start_end = zeros(d_order*2, n_all_poly);
% i*(i-1)*...*(i-j+1)*ts(1)^(i-j)
beq_start_end = [start_cond.';end_cond.'];
for m = 1:2
    row_offset = (m-1)*d_order;
    col_offset = (m-1)*n_coef1d*(n_seg-1);
    row_idx = row_offset+1;
    % p constraint
    N = N_set(:,m).';
    indices = 1+col_offset:n_coef1d+col_offset;
    Aeq_start_end(row_idx,indices) = N;
    row_idx = row_idx + 1;
    % v constraint
    N = N_set(:,m+2).';
    A_temp = zeros(n_coef1d_d,n_coef1d);
    for i = 1:n_coef1d_d
        c0 = -p/(knots(i+p+2)-knots(i+2));
        c1 = -c0;
        A_temp(i,i) = c0*N(i);
        A_temp(i,i+1) = c1*N(i);
    end
    indices = 1+col_offset:n_coef1d+col_offset;
    Aeq_start_end(row_idx,indices) = sum(A_temp,1);
    row_idx = row_idx + 1;
    % a constraint
    N = N_set(:,m+4).';
    A_temp = zeros(n_coef1d_d2,n_coef1d);
    for i = 1:n_coef1d_d2
        c0 = dp/(dknots(i+dp+2)-dknots(i+1));
        c1 = p/(knots(i+p+2)-knots(i+1));
        c3 = p/(knots(i+p+3)-knots(i+2));
        c2 = -c1-c3;
        A_temp(i,i) = c0*c1*N(i);
        A_temp(i,i+1) = c0*c2*N(i);
        A_temp(i,i+2) = c0*c3*N(i);
    end
    indices = 1+col_offset:n_coef1d+col_offset;
    Aeq_start_end(row_idx,indices) = sum(A_temp,1);
end
Aeq = [Aeq_start_end];
beq = [beq_start_end];    

end