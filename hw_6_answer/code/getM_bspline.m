function M = getM_bspline(p)
% General matrices for UNIFORM!! bspline curves
% p: degree
% Basis = [B[i-p,p+1](u) ... B[i,p+1](u)]
% U = [1 u u^2 ... u^p]
% Basis = U * M[p+1](i), u = (t-t[i])/(t[i+1]-t[i]), u in [0,1)
% Control points P = [P[i-p] ... P[i]]
% The bspline curve seg C[i-p](u) = U * M[p+1](i) * P
% M[p+1](i) is referd to as the ith basis matrix of bspline basis functions
% of degree p, (p+1) x (p+1) matrix
if p == 0
    M = [1];
else
    M_p_1 = getM_bspline(p-1);
    A = zeros(p,p+1);
    B = zeros(p,p+1);
    for i = 1:p
        A(i,i) = i;
        A(i,i+1) = p-i;
        B(i,i) = -1;
        B(i,i+1) = 1;
    end
    assignin('base','A',A)
    assignin('base','B',B)
    Z = zeros(1,size(M_p_1,2));
%     fprintf('p: %d, Z colns: %d\n', p,size(Z,2));
    M = 1/p.*([M_p_1;Z]*A+[Z;M_p_1]*B);
end

end