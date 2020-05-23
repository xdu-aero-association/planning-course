function Ct = getCt(n_seg, n_order)
%#####################################################
% STEP 2.1: finish the expression of Ct
%
%
%
%
%
k =n_seg;


Ct=zeros(k*8  , 4*(k+1) );

Ct(1,1) = 1;
Ct(2,2) = 1;
Ct(3,3) = 1;
Ct(4,4) = 1;


Ct(k*8-3 ,k+4) = 1;
Ct(k*8 -2,k+5) = 1;
Ct(k*8-1 ,k+6) = 1;
Ct(k*8 ,k+7) = 1;

for  i = 1:(k-1)  
       Ct( 8*i-3, 4+i) = 1;
    Ct( 8*i+1, 4+i) = 1;
end



for i = 1:(k-1)  
    for j=0:2  
        Ct( 8*i-8+6+j, k+7+j+1+3*i-3) = 1;
        Ct( 8*i+2+j, k+7+j+1+3*i-3) = 1;
    end 
end





end