function [Aeq, beq] = getAbeq(n_seg, n_order, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    d_order = 3;
    %#####################################################
    % STEP 2.1 p,v,a constraint in start 
    Aeq_start = zeros(d_order, n_all_poly);
    beq_start = zeros(d_order, 1);
    % i*(i-1)*...*(i-j+1)*ts(1)^(i-j)
    beq_start = start_cond.';
    for j = 0:d_order-1
        Aeq_start(j+1,j+1) = factorial(j); % this time j=k
    end
    %#####################################################
    % STEP 2.2 p,v,a constraint in end
    Aeq_end = zeros(d_order, n_all_poly);
    beq_end = zeros(d_order, 1);
    % STEP 2.2: write expression of Aeq_end and beq_end
    beq_end = end_cond.';
    for i = 0:d_order-1
        for j = i:n_order
            Aeq_end(i+1,j+1+(n_seg-1)*(n_order+1)) = factorial(j)/factorial(j-i)*ts(end)^(j-i);
        end
    end
    %#####################################################
    % STEP 2.3 position continuity constrain between 2 segments
    [Aeq_con_p, beq_con_p] = get_continuity_constraint(n_seg, n_order, ts, 0);
    %#####################################################
    % STEP 2.4 velocity continuity constrain between 2 segments
    [Aeq_con_v, beq_con_v] = get_continuity_constraint(n_seg, n_order, ts, 1);
    %#####################################################
    % STEP 2.5 acceleration continuity constrain between 2 segments
    [Aeq_con_a, beq_con_a] = get_continuity_constraint(n_seg, n_order, ts, 2);
    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a];
    beq_con = [beq_con_p; beq_con_v; beq_con_a];
    Aeq = [Aeq_start; Aeq_end; Aeq_con];
    beq = [beq_start; beq_end; beq_con];    
end

function [Aeq, beq] = get_continuity_constraint(n_seg, n_order, ts, order_deri)
    k = order_deri; % order of deri
    n_all_poly = n_seg*(n_order+1);
    Aeq = zeros(n_seg-1, n_all_poly);
    beq = zeros(n_seg-1, 1);
    for i = 1:n_seg-1
        for j = k:n_order
            Aeq(i,j+1+(i-1)*(n_order+1)) = factorial(j)/factorial(j-k)*ts(i)^(j-k); % end pos/vel/... of 1st to 2nd last seg
        end
        Aeq(i,1+k+i*(n_order+1)) = -factorial(k); % start pos/vel/... of 2nd to last seg, this time j=k, dont forget the minus sign!!!!!
        beq(i) = 0;
    end
end