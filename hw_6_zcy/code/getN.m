function [N] = getN(n,p)
% calculate matrix for transform btwn bspline cost function params P & real
% control points P_real, size (n-p+1) x (n+1)
% P = N*P_real, J = P.'*Q*P = P_real.'*N.'*Q*N*P_real, Q_real = N.'*Q*N
% basically one curve segment relates to p+1 control points (P[i] P[i+1]
% ... P[i+p]), i = 0:n_seg-1
n_seg = n-p+1;
N = zeros(n_seg*(p+1),n+1);
for i = 1:n_seg
    row_indices = 1+(i-1)*(p+1):1+p+(i-1)*(p+1);
    coln_indices = 1+(i-1):1+p+(i-1);
    N(row_indices,coln_indices) = eye(p+1);
end

end

