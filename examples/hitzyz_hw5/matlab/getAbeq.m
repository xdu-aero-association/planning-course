function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
n_all_poly = n_seg*(n_order+1);
%#####################################################
% p,v,a,j constraint in start,
Aeq_start = zeros(4, n_all_poly);
beq_start = zeros(4, 1);
% STEP 2.1: write expression of Aeq_start and beq_start
%
%
%
%

Aeq_start(1,1) = factorial(0);
Aeq_start(2,2) =  factorial(1);
Aeq_start(3,3) =  factorial(2);
Aeq_start(4,4) =  factorial(3);
beq_start(1) = start_cond(1);
beq_start(2) = start_cond(2);
beq_start(3) = start_cond(3);
beq_start(4) = start_cond(4);
%#####################################################
% p,v,a constraint in end
Aeq_end = zeros(4, n_all_poly);
beq_end = zeros(4, 1);
% STEP 2.2: write expression of Aeq_end and beq_end
%
%
%
%
beq_end(1) = end_cond(1);
beq_end(2) = end_cond(2);
beq_end(3) = end_cond(3);
beq_end(4) = end_cond(4);

for k_end = 0:3
    for i_end=0:n_order
        if i_end>=k_end
            
            Aeq_end(k_end+1,n_all_poly -n_order-1+i_end +1) = factorial(i_end)/ factorial(i_end-k_end)*ts(n_seg)^(i_end-k_end);
            
        end
    end
end


%#####################################################
% position constrain in all middle waypoints
Aeq_wp = zeros(n_seg-1, n_all_poly);
beq_wp = zeros(n_seg-1, 1);
% STEP 2.3: write expression of Aeq_wp and beq_wp
%
%
%
%
for r = 1:n_seg-1
    Aeq_wp(r,r*(n_order+1)+1) = 1;
    beq_wp(r) = waypoints(r+1);
end


%#####################################################
% position continuity constrain between each 2 segments
Aeq_con_p = zeros(n_seg-1, n_all_poly);
beq_con_p = zeros(n_seg-1, 1);
% STEP 2.4: write expression of Aeq_con_p and beq_con_p
%
%
%
%
for r = 1:n_seg-1
    for i=0:n_order
        Aeq_con_p(r,(r-1)*(n_order+1)+i+1) = ts(r)^(i);
    end
    Aeq_con_p(r,r*(n_order+1)+1) =-1;
end


%#####################################################
% velocity continuity constrain between each 2 segments
Aeq_con_v = zeros(n_seg-1, n_all_poly);
beq_con_v = zeros(n_seg-1, 1);
% STEP 2.5: write expression of Aeq_con_v and beq_con_v
%
%
%
%


for r = 1:n_seg-1
    
    k_end = 1;
    for i_end=0:n_order
        if i_end>=k_end
            
            Aeq_con_v(r, (r-1)*(n_order+1)+i_end +1) = factorial(i_end)/ factorial(i_end-k_end)*ts(r)^(i_end-k_end);
            
        end
        
        if i_end==k_end
            
            Aeq_con_v(r, (r)*(n_order+1)+i_end+1 ) = -factorial(i_end)/ factorial(i_end-k_end)*ts(r)^(i_end-k_end);
            
        end
        
    end
    
end



%#####################################################
% acceleration continuity constrain between each 2 segments
Aeq_con_a = zeros(n_seg-1, n_all_poly);
beq_con_a = zeros(n_seg-1, 1);
% STEP 2.6: write expression of Aeq_con_a and beq_con_a
%
%
%
%


for r = 1:n_seg-1
    
    k_end = 2;
    for i_end=0:n_order
        if i_end>=k_end
            
            Aeq_con_a(r, (r-1)*(n_order+1)+i_end +1) = factorial(i_end)/ factorial(i_end-k_end)*ts(r)^(i_end-k_end);
            
        end
        
        if i_end==k_end
            
            Aeq_con_a(r, (r)*(n_order+1)+i_end+1 ) = -factorial(i_end)/ factorial(i_end-k_end)*ts(r)^(i_end-k_end);
            
        end
        
    end
    
end






%#####################################################
% jerk continuity constrain between each 2 segments
Aeq_con_j = zeros(n_seg-1, n_all_poly);
beq_con_j = zeros(n_seg-1, 1);
% STEP 2.7: write expression of Aeq_con_j and beq_con_j
%
%
%
%


for r = 1:n_seg-1
    
    k_end = 3;
    for i_end=0:n_order
        if i_end>=k_end
            
            Aeq_con_j(r, (r-1)*(n_order+1)+i_end +1) = factorial(i_end)/ factorial(i_end-k_end)*ts(r)^(i_end-k_end);
            
        end
        
        if i_end==k_end
            
            Aeq_con_j(r, (r)*(n_order+1)+i_end+1 ) = -factorial(i_end)/ factorial(i_end-k_end)*ts(r)^(i_end-k_end);
            
        end
        
    end
    
end


%#####################################################
% combine all components to form Aeq and beq
Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
beq = [beq_start; beq_end; beq_wp; beq_con];
end