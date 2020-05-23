function Ct = getCt(n_seg, n_order)
    %#####################################################
    % STEP 2.1: finish the expression of Ct
    % df: constraint, dp: free
    Ct = zeros(8*n_seg,4*(n_seg+1));
    for i = 0:3
        Ct(i+1,i+1) = 1; % 1st point pos/vel/... of 1st seg (df)
        Ct(i+1+8*n_seg-4,i+1+4+n_seg-1) = 1; % last point pos/vel/... of last seg (df)
    end
    for i = 1:n_seg-1
        Ct(5+8*(i-1),4+i) = 1; % mid points pos (df)
        Ct(1+8*i,4+i) = 1; % mid points pos (df)
        for j = 0:2
            Ct(j+1+5+8*(i-1),j+1+3*(i-1)+8+n_seg-1) = 1; % mid points vel/acc/jerk (dp)
            Ct(j+1+1+8*i,j+1+3*(i-1)+8+n_seg-1) = 1; % mid points vel/acc/jerk (dp)
        end
    end
end