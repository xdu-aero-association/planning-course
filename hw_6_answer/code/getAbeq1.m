function [Aeq, beq] = getAbeq1(n_seg, n_order, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    d_order = 3;
    n_coef1d = n_order+1;
    %#####################################################
    % STEP 2.1 p,v,a constraint in start & end
    ts_set = [0, ts(end)];
    Aeq_start_end = zeros(d_order*2, n_all_poly);
    % i*(i-1)*...*(i-j+1)*ts(1)^(i-j)
    beq_start_end = [start_cond.';end_cond.'];
    for m = 1:2
        row_offset = (m-1)*d_order;
        col_offset = (m-1)*n_coef1d*(n_seg-1);
        row_idx = row_offset+1;
        % p constraint
        for i = 1:n_coef1d
            Aeq_start_end(row_idx,i+col_offset) = get_bern(n_order,i-1,ts_set(m));
        end
        row_idx = row_idx + 1;
        % v constraint
        Aeq_start_end(row_idx,1+col_offset) = -get_bern(n_order-1,0,ts_set(m));
        Aeq_start_end(row_idx,n_coef1d+col_offset) = get_bern(n_order-1,n_order-1,ts_set(m));
        for i = 2:n_coef1d-1
            Aeq_start_end(row_idx,i+col_offset) = get_bern(n_order-1,i-2,ts_set(m)) - get_bern(n_order-1,i-1,ts_set(m));
        end
        indices = 1+col_offset:n_coef1d+col_offset;
        Aeq_start_end(row_idx,indices) = Aeq_start_end(row_idx,indices) * n_order;
        row_idx = row_idx + 1;
        % a constraint
        Aeq_start_end(row_idx,1+col_offset) = get_bern(n_order-2,0,ts_set(m));
        Aeq_start_end(row_idx,2+col_offset) = -2 * get_bern(n_order-2,0,ts_set(m)) + get_bern(n_order-2,1,ts_set(m));
        Aeq_start_end(row_idx,n_coef1d-1+col_offset) = get_bern(n_order-2,n_order-3,ts_set(m)) - 2 * get_bern(n_order-2,n_order-2,ts_set(m));
        Aeq_start_end(row_idx,n_coef1d+col_offset) = get_bern(n_order-2,n_order-2,ts_set(m));
        for j = 3:n_coef1d-2
            Aeq_start_end(row_idx,j+col_offset) = get_bern(n_order-2,j-3,ts_set(m)) - 2 * get_bern(n_order-2,j-2,ts_set(m)) + get_bern(n_order-2,j-1,ts_set(m));
        end
        indices = 1+col_offset:n_coef1d+col_offset;
        Aeq_start_end(row_idx,indices) = Aeq_start_end(row_idx,indices) * n_order * (n_order-1) * (-1)^(m+1);
    end
    %#####################################################
    % STEP 2.2 p,v,a constraint in end
    %#####################################################
    % STEP 2.3 position continuity constrain between 2 segments
    [Aeq_con_p, beq_con_p] = get_continuity_constraint_p(n_seg, n_order, ts);
    %#####################################################
    % STEP 2.4 velocity continuity constrain between 2 segments
    [Aeq_con_v, beq_con_v] = get_continuity_constraint_v(n_seg, n_order, ts);
    %#####################################################
    % STEP 2.5 acceleration continuity constrain between 2 segments
    [Aeq_con_a, beq_con_a] = get_continuity_constraint_a(n_seg, n_order, ts);
    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a];
    beq_con = [beq_con_p; beq_con_v; beq_con_a];
    Aeq = [Aeq_start_end; Aeq_con];
    beq = [beq_start_end; beq_con];    
end

function [Aeq, beq] = get_continuity_constraint_p(n_seg, n_order, ts)
    n_all_poly = n_seg*(n_order+1);
    n_coef1d = n_order+1;
    Aeq = zeros(n_seg-1, n_all_poly);
    beq = zeros(n_seg-1, 1);
    for j = 1:n_seg-1
        for i = 1:n_coef1d
            Aeq(j,(j-1)*n_coef1d+i) = get_bern(n_order,i-1,ts(j));
            Aeq(j,j*n_coef1d+i) = -get_bern(n_order,i-1,0);
        end
    end
end

function [Aeq, beq] = get_continuity_constraint_v(n_seg, n_order, ts)
    n_all_poly = n_seg*(n_order+1);
    n_coef1d = n_order+1;
    Aeq = zeros(n_seg-1, n_all_poly);
    beq = zeros(n_seg-1, 1);
    for j = 1:n_seg-1
        ts_set = [ts(j) 0];
        for m = 1:2
            col_offset = (m-1)*n_coef1d+(j-1)*n_coef1d;
            Aeq(j,1+col_offset) = -get_bern(n_order-1,0,ts_set(m));
            Aeq(j,n_coef1d+col_offset) = get_bern(n_order-1,n_order-1,ts_set(m));
            for i = 2:n_coef1d-1
                Aeq(j,i+col_offset) = get_bern(n_order-1,i-2,ts_set(m)) - get_bern(n_order-1,i-1,ts_set(m));
            end
            indices = 1+col_offset:n_coef1d+col_offset;
            Aeq(j,indices) = Aeq(j,indices) * n_order * (-1)^(m+1);
        end
    end
end

function [Aeq, beq] = get_continuity_constraint_a(n_seg, n_order, ts)
    n_all_poly = n_seg*(n_order+1);
    n_coef1d = n_order+1;
    Aeq = zeros(n_seg-1, n_all_poly);
    beq = zeros(n_seg-1, 1);
    for j = 1:n_seg-1
        ts_set = [ts(j) 0];
        for m = 1:2
            col_offset = (m-1)*n_coef1d+(j-1)*n_coef1d;
            Aeq(j,1+col_offset) = get_bern(n_order-2,0,ts_set(m));
            Aeq(j,2+col_offset) = -2 * get_bern(n_order-2,0,ts_set(m)) + get_bern(n_order-2,1,ts_set(m));
            Aeq(j,n_coef1d-1+col_offset) = get_bern(n_order-2,n_order-3,ts_set(m)) - 2 * get_bern(n_order-2,n_order-2,ts_set(m));
            Aeq(j,n_coef1d+col_offset) = get_bern(n_order-2,n_order-2,ts_set(m));
            for i = 3:n_coef1d-2
                Aeq(j,i+col_offset) = get_bern(n_order-2,i-3,ts_set(m)) - 2 * get_bern(n_order-2,i-2,ts_set(m)) + get_bern(n_order-2,i-1,ts_set(m));
            end
            indices = 1+col_offset:n_coef1d+col_offset;
            Aeq(j,indices) = Aeq(j,indices) * n_order * (n_order-1) * (-1)^(m+1);
        end
    end
end

% calculate bernstein basis
function [bni] = get_bern(n,i,t)
    i = max(min(n,i),0);
    bni = 1;
    for j = n:-1:n-i+1
        bni = bni * j;
    end
    bni = bni / factorial(i);
    bni = bni * t^i * (1-t)^(n-i);
end

% function [bni] = get_bern(n,i,t)
%     i = max(min(n,i),0);
%     bni = factorial(n) / factorial(i) / factorial(n-i);
%     bni = bni * t^i * (1-t)^(n-i);
%     fprintf(1,'i:%d, n:%d, bni:%.1f\n',i,n,bni)
% end