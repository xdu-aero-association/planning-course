function M = getM(n_seg, n_order, ts)
    M = [];
    for i = 1:n_seg
        M_i = [];
        %#####################################################
        % STEP 1.1: calculate M_i of the k-th segment 
        M_i = zeros(8,n_order+1);
        for k = 0:3
            M_i(k+1,k+1) = factorial(k); % start pos/vel/acc/jerk of jth seg
            for j = k:n_order
                M_i(k+5,j+1) = factorial(j)/factorial(j-k)*ts(i)^(j-k); % end pos/vel/acc/jerk of jth seg
            end
        end
        M = blkdiag(M, M_i);
    end
end