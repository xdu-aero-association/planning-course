function [Q] = getQ_bspline(n, p)
% Calculate smoothness cost Js by elastic band method
% n+1 is number of control points P
% p is degree of bspline, p >= 2(must)
method = 0;
if method == 0
    % method 1 elastic band method
    Q = zeros(n+1,n+1);
    Q_temp = [1 -2 1;-2 4 -2;1 -2 1];
    for i = p:n-p+2
        Q(i-1:i+1,i-1:i+1) = Q(i-1:i+1,i-1:i+1) + Q_temp;
    end
    assignin('base','Q_elastic',Q)
else
    % method 2 conventional
    Q = [];
    M = [];
    knots = deboor_knot(p,n,2); % we need each curve seg to be clamped at both ends
    n_seg = n-p+1; % number of curve segments
    n_coeffs1d = p+1; % number of coeffs for one segment
    M_k = getM_bspline(p); % matrix for coeffs = M * c_pts[m-p,m], m=[p,n+1]
    n_order = p;
    for k = 1:n_seg
        Q_k = zeros(n_order+1,n_order+1);
        %#####################################################
        % STEP 2.1 calculate Q_k of the k-th segment 
        dt = knots(k+1+p)-knots(k+p);
        if n_order >= 4 % min snap
            % (i*(i-1)*(i-2)*(i-3)*j*(j-1)*(j-2)*(j-3))/(i+j-7)*ts(k)^(i+j-7)
            for i = 4:n_order
                for j = 4:n_order
                    Q_k(i+1,j+1) = (i*(i-1)*(i-2)*(i-3)*j*(j-1)*(j-2)*(j-3))/(i+j-7)*dt^(i+j-7);
                end
            end
        else % min jerk
            for i = 3:n_order
                for j = 3:n_order
                    Q_k(i+1,j+1) = (i*(i-1)*(i-2)*j*(j-1)*(j-2))/(i+j-5)*dt^(i+j-5);
                end
            end
        end
        Q = blkdiag(Q, Q_k);
        M = blkdiag(M, M_k);
    end
    % N for rearrange control points
    [N] = getN(n,p);
    assignin('base','Q_origin',Q)
    assignin('base','M',M)
    assignin('base','N',N)
    Q_0 = N.'*M'*Q*M*N;
    % Q_0 = nearestSPD(Q_0);
    assignin('base','Q_0',Q_0)
    Q = Q_0;
end

end