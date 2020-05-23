syms Delta_px Delta_py Delta_pz Delta_vx Delta_vy Delta_vz
syms px0 py0 pz0 vx0 vy0 vz0 px_F py_F pz_F
syms alpha1 alpha2 alpha3 beta1 beta2 beta3 f(T)
Delta_px = px_F - vx0 * T - px0;
Delta_py = py_F - vy0 * T - py0;
Delta_pz = pz_F - vz0 * T - pz0;
Delta_vx = - vx0;
Delta_vy = - vy0;
Delta_vz = - vz0;
matrixT = [-12/T^3 0 0 6/T^2 0 0;
    0 -12/T^3 0 0 6/T^2 0;
    0 0 -12/T^3 0 0 6/T^2;
    6/T^2 0 0 -2/T 0 0;
    0 6/T^2 0 0 -2/T 0;
    0 0 6/T^2 0 0 -2/T;];
ab = matrixT * [Delta_px Delta_py Delta_pz Delta_vx Delta_vy Delta_vz].';
alpha1 = (12*(px0 - px_F + T*vx0))/T^3 - (6*vx0)/T^2;
alpha2 = (12*(py0 - py_F + T*vy0))/T^3 - (6*vy0)/T^2;
alpha3 = (12*(pz0 - pz_F + T*vz0))/T^3 - (6*vz0)/T^2;
beta1 = (2*vx0)/T - (6*(px0 - px_F + T*vx0))/T^2;
beta2 = (2*vy0)/Tv - (6*(py0 - py_F + T*vy0))/T^2;
beta3 = (2*vz0)/T - (6*(pz0 - pz_F + T*vz0))/T^2;
f(T) = T + (1/3 * alpha1^2 *T^3 + alpha1 * beta1 * T^2 + beta1^2*T) + ...
    (1/3 * alpha2^2 *T^3 + alpha2 * beta2 * T^2 + beta2^2*T) + ...
    (1/3 * alpha3^2 *T^3 + alpha3 * beta3 * T^2 + beta3^2*T);
eqn = diff(f,T)==0;
df_sim = simplify(diff(f,T));
S = solve(eqn,T);




