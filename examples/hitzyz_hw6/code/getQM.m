function [Q, M] = getQM(n_seg, n_order, ts)
    Q = [];
    M = [];
    M_k = getM(n_order);
    for k = 1:n_seg
        %#####################################################
        % STEP 2.1 calculate Q_k of the k-th segment 
        Q_k = [];
        qa = [];
        for j=4:n_order
            for i = 4:n_order
                
             qa =[qa    i*(i-1)*(i-2)*(i-3)*j*(j-1)*(j-2)*(j-3)*ts(k)^(i+j-7)/(i+j-7)];
                
            end
        end
        Q_k= blkdiag(zeros(4,4), reshape(qa,[4,4])');
        Q = blkdiag(Q, Q_k);
        M = blkdiag(M, M_k);
    end
end