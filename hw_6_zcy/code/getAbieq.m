function [Aieq, bieq] = getAbieq(n_seg, n_order, corridor_range, ts, v_max, a_max)
    n_all_poly = n_seg*(n_order+1);
    n_coef1d = n_order+1;
    %#####################################################
    % STEP 3.2.1 p constraint
    Aieq_p = [eye(n_all_poly);-eye(n_all_poly)];
    bieq_p = zeros(n_all_poly*2,1);
    for i = 1:n_seg
        bieq_p((i-1)*n_coef1d+1:i*n_coef1d) = ones(n_coef1d,1)*corridor_range(i,2); % max limit
        bieq_p((i-1)*n_coef1d+1+n_all_poly:i*n_coef1d+n_all_poly) = -ones(n_coef1d,1)*corridor_range(i,1); % min limit
    end
    %#####################################################
    % STEP 3.2.2 v constraint   
    n_coef1d_d = n_coef1d-1;
    n_all_poly_d = n_seg*n_coef1d_d;
    Aieq_v = zeros(n_all_poly_d*2,n_all_poly);
    bieq_v = zeros(n_all_poly_d*2,1);
    for i = 1:n_seg
        for j = 1:n_coef1d_d
            Aieq_v((i-1)*n_coef1d_d+j,j) = -1;
            Aieq_v((i-1)*n_coef1d_d+j,j+1) = 1;
            Aieq_v((i-1)*n_coef1d_d+j+n_all_poly_d,j) = 1;
            Aieq_v((i-1)*n_coef1d_d+j+n_all_poly_d,j+1) = -1;
        end
        bieq_v((i-1)*n_coef1d_d+1:i*n_coef1d_d) = ones(n_coef1d_d,1)*v_max;
        bieq_v((i-1)*n_coef1d_d+1+n_all_poly_d:i*n_coef1d_d+n_all_poly_d) = ones(n_coef1d_d,1)*v_max;
    end
    Aieq_v = Aieq_v * n_order;
    %#####################################################
    % STEP 3.2.3 a constraint   
    n_coef1d_d2 = n_coef1d_d-1;
    n_all_poly_d2 = n_seg*n_coef1d_d2;
    Aieq_a = zeros(n_all_poly_d2*2,n_all_poly);
    bieq_a = zeros(n_all_poly_d2*2,1);
    for i = 1:n_seg
        for j = 1:n_coef1d_d2
            Aieq_a((i-1)*n_coef1d_d2+j,j) = 1;
            Aieq_a((i-1)*n_coef1d_d2+j,j+1) = -2;
            Aieq_a((i-1)*n_coef1d_d2+j,j+2) = 1;
            Aieq_a((i-1)*n_coef1d_d2+j+n_all_poly_d2,j) = 1;
            Aieq_a((i-1)*n_coef1d_d2+j+n_all_poly_d2,j+1) = -2;
            Aieq_a((i-1)*n_coef1d_d2+j+n_all_poly_d2,j+2) = 1;
        end
        bieq_a((i-1)*n_coef1d_d2+1:i*n_coef1d_d2) = ones(n_coef1d_d2,1)*a_max;
        bieq_a((i-1)*n_coef1d_d2+1+n_all_poly_d2:i*n_coef1d_d2+n_all_poly_d2) = ones(n_coef1d_d2,1)*a_max;
    end
    Aieq_a = Aieq_a * n_order*(n_order-1);
    %#####################################################
    % combine all components to form Aieq and bieq   
    Aieq = [Aieq_p; Aieq_v; Aieq_a];
    bieq = [bieq_p; bieq_v; bieq_a];
    Aieq = Aieq_p;
    bieq = bieq_p;
end