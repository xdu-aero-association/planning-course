function [Aieq, bieq] = getAbieq(n_seg, n_order, corridor_range, ts, v_max, a_max)
n_all_poly = n_seg*(n_order+1);
%#####################################################
% STEP 3.2.1 p constraint


Aieq_p = [eye(n_all_poly); -eye(n_all_poly) ];

bieq_p1 = zeros(n_all_poly,1);
for i = 1:n_seg
    for j = 1: n_order+1
        bieq_p1( (i-1)*(n_order+1) +j )  = corridor_range(i,2);
    end
end
bieq_p2 = zeros(n_all_poly,1);
for i = 1:n_seg
    for j = 1: n_order+1
        bieq_p2( (i-1)*(n_order+1) +j )  =- corridor_range(i,1);
    end
end


bieq_p = [bieq_p1;bieq_p2];
%#####################################################
% STEP 3.2.2 v constraint
Aieq_v =  zeros(n_seg*n_order,n_seg*(n_order+1));
bieq_v =  v_max*ones( n_seg*n_order,1);

for i = 1:n_seg
    for j = 1: n_order
        Aieq_v(  (i-1)*(n_order) +j,(i-1)*(n_order) +j+1 )  =1*n_order;
        Aieq_v(  (i-1)*(n_order) +j,(i-1)*(n_order) +j )  =-1*n_order;
    end
end


%#####################################################
% STEP 3.2.3 a constraint
Aieq_a = zeros(n_seg*(n_order-1),n_seg*(n_order+1));
bieq_a = a_max*ones( n_seg*(n_order-1),1);

for i = 1:n_seg
    for j = 1: n_order-1
        Aieq_a(  (i-1)*(n_order-1) +j,(i-1)*(n_order-1) +j +2)  =1*n_order*( n_order-1);
        Aieq_a(  (i-1)*(n_order-1) +j,(i-1)*(n_order-1) +j+1 )  =-2*n_order*( n_order-1);
        Aieq_a(  (i-1)*(n_order-1) +j,(i-1)*(n_order-1) +j )  =1*n_order*( n_order-1);
    end
end



%#####################################################
% combine all components to form Aieq and bieq
Aieq = [Aieq_p; Aieq_v; Aieq_a];
bieq = [bieq_p; bieq_v; bieq_a];
Aieq = Aieq_p;
bieq = bieq_p;
end