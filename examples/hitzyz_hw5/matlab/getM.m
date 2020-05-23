function M = getM(n_seg, n_order, ts)
M = [];
for k = 1:n_seg
    M_k = [];
    M_k = zeros(8,8);
    %#####################################################
    % STEP 1.1: calculate M_k of the k-th segment
    %
    %
    %
    %
    M_k(1,1) = 1;
    M_k(2,2) = 1;
    M_k(3,3) = 2;
    M_k(4,4) = 6;
    
    
    for k_end = 0:3
    for i_end=0:n_order
        if i_end>=k_end
            
            M_k(k_end+5, i_end +1) = factorial(i_end)/ factorial(i_end-k_end)*ts(n_seg)^(i_end-k_end);
            
        end
    end
end

    
    
    
    M = blkdiag(M, M_k);
end
end