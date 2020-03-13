#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}
/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment
    int d_num = d_order;  // dimension of one point for close form sol.

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);
    ROS_WARN("[TG] d_order: %d", d_order);
    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    for (int idx = 0; idx < 3; ++idx) {
        ROS_WARN("[TG] calculating %d dimension trajectory", idx);
        VectorXd start_pose = VectorXd::Zero(d_num);
        start_pose(0) = Path(0, idx);
        start_pose(1) = Vel(0, idx);
        start_pose(2) = Acc(0, idx);
        VectorXd end_pose = VectorXd::Zero(d_num);
        end_pose(0) = Path(m, idx);
        end_pose(1) = Vel(1, idx);
        end_pose(2) = Acc(1, idx);
        int n_all_poly = m * p_num1d;
        
        // calculate Aeq_start
        MatrixXd Aeq_start = MatrixXd::Zero(d_num, n_all_poly);
        VectorXd beq_start = VectorXd::Zero(d_num);
        for (int k = 0; k < d_order; ++k) {
            Aeq_start(k,k) = Factorial(k);
        }
        beq_start = start_pose;
        ROS_INFO("[TG] cal Aeq_start done, size: %ld, %ld", Aeq_start.rows(), Aeq_start.cols());
        // calculate Aeq_end
        MatrixXd Aeq_end = MatrixXd::Zero(d_num, n_all_poly);
        VectorXd beq_end = VectorXd::Zero(d_num);
        for (int k = 0; k < d_order; ++k) {
            for (int i = k; i <= p_order; ++i) {
                Aeq_end(k, (m-1)*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(m-1), i-k);
            }
        }
        beq_end = end_pose;
        ROS_INFO("[TG] cal Aeq_end done, size: %ld, %ld", Aeq_end.rows(), Aeq_end.cols());
        // calculate Aeq_wp
        MatrixXd Aeq_wp = MatrixXd::Zero(m-1, n_all_poly);
        VectorXd beq_wp = VectorXd::Zero(m-1);
        for (int j = 0; j < m-1; ++j) {
            Aeq_wp(j, (j+1)*p_num1d) = 1;
            beq_wp(j) = Path(j+1, idx);
        }
        ROS_INFO("[TG] cal Aeq_wp done, size: %ld, %ld", Aeq_wp.rows(), Aeq_wp.cols());
        // calculate pos continuity
        MatrixXd Aeq_con_p = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, p_order, Time, 0, Aeq_con_p);
        ROS_INFO("[TG] cal Aeq_con_p done, size: %ld, %ld", Aeq_con_p.rows(), Aeq_con_p.cols());
        // calculate vel continuity
        MatrixXd Aeq_con_v = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, p_order, Time, 1, Aeq_con_v);
        ROS_INFO("[TG] cal Aeq_con_v done, size: %ld, %ld", Aeq_con_v.rows(), Aeq_con_v.cols());
        // calculate acc continuity
        MatrixXd Aeq_con_a = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, p_order, Time, 2, Aeq_con_a);
        ROS_INFO("[TG] cal Aeq_con_a done, size: %ld, %ld", Aeq_con_a.rows(), Aeq_con_a.cols());
        // calculate jerk continuity
        MatrixXd Aeq_con_j = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, p_order, Time, 3, Aeq_con_j);
        ROS_INFO("[TG] cal Aeq_con_j done, size: %ld, %ld", Aeq_con_j.rows(), Aeq_con_j.cols());
        MatrixXd Aeq;
        VectorXd beq;
        int Aeq_row_num = 2*(d_num) + 5*(m-1);
        Aeq.resize(Aeq_row_num, n_all_poly);
        beq.resize(Aeq_row_num);
        Aeq << Aeq_start,
               Aeq_end,
               Aeq_wp,
               Aeq_con_p,
               Aeq_con_v,
               Aeq_con_a,
               Aeq_con_j;
        beq << beq_start,
               beq_end,
               beq_wp,
               VectorXd::Zero(4*(m-1));
        ROS_WARN("[TG] concatenated Aeq & beq");
    /*   Produce the dereivatives in X, Y and Z axis directly.  */
        MatrixXd M = MatrixXd::Zero(2*(d_num)*m, p_num1d*m);
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < d_order; ++k) {
                M(j*p_num1d+k, j*p_num1d+k) = Factorial(k);
                for (int i = k; i <= p_order; ++i) {
                    M(j*p_num1d+4+k, j*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(j), i-k);
                }
            }
        }
        ROS_INFO("[TG] cal M done, size: %ld, %ld", M.rows(), M.cols());
        std::cout << M << std::endl << std::endl;
        MatrixXd Ct = MatrixXd::Zero(2*d_num*m, d_num*(m+1));
        for (int k = 0; k < d_order; ++k) {
            Ct(k, k) = 1;
            Ct(k+d_num+2*d_num*(m-1), k+d_num+m-1) = 1;
        }
        for (int j = 0; j < m-1; ++j) {
            Ct(d_num+2*d_num*j,d_num+j) = 1;
            Ct(2*d_num*(j+1), d_num+j) = 1;
            for (int i = 0; i < d_order-1; ++i) {
                Ct(d_num+1+2*d_num*j+i, 2*d_num+m-1+j*(d_order-1)+i) = 1;
                Ct(1+2*d_num*(j+1)+i, 2*d_num+m-1+j*(d_order-1)+i) = 1;
            }
        }
        ROS_INFO("[TG] cal Ct done, size: %ld, %ld", Ct.rows(), Ct.cols());
        std::cout << Ct << std::endl << std::endl;
    /*   Produce the Minimum Snap cost function, the Hessian Matrix   */
        MatrixXd Q = MatrixXd::Zero(m*p_num1d, m*p_num1d);
        for (int j = 0; j < m; ++j) {
            for (int k = 4; k <= p_order; ++k) {
                for (int l = 4; l <= p_order; ++l) {
                    Q(k+j*p_num1d, l+j*p_num1d) = k*(k-1)*(k-1)*(k-3)*l*(l-1)*(l-2)*(l-3)*pow(Time(j), k+l-7)/(k+l-7);
                }
            }
        }
        ROS_INFO("[TG] cal Q done, size: %ld, %ld", Q.rows(), Q.cols());
        std::cout << Q << std::endl << std::endl;
        // compute R
        MatrixXd R = Ct.transpose() * (M.inverse()).transpose() * Q * M.inverse() * Ct;
        ROS_INFO("[TG] cal R done, size: %ld, %ld", R.rows(), R.cols());
        std::cout << R << std::endl << std::endl;
        // compute df
        VectorXd df = VectorXd::Zero(2*d_num+m-1);
        df.segment(0,d_num) = start_pose;
        df.segment(d_num+m-1,d_num) = end_pose;
        for (int j = 0; j < m-1; ++j) {
            df(j+d_num) = Path(j+1, idx);
        }
        ROS_INFO("[TG] cal df done, size: %ld, %ld", df.rows(), df.cols());
        std::cout << df << std::endl << std::endl;
        MatrixXd R_pp = R.block(2*d_num+m-1, 2*d_num+m-1, (d_order-1)*(m-1), (d_order-1)*(m-1));
        ROS_INFO("[TG] cal R_pp done, size: %ld, %ld", R_pp.rows(), R_pp.cols());
        MatrixXd R_fp = R.block(0, 2*d_num+m-1, 2*d_num+m-1, (d_order-1)*(m-1));
        ROS_INFO("[TG] cal R_fp done, size: %ld, %ld", R_fp.rows(), R_fp.cols());
        VectorXd dp_best = -R_pp.inverse() * R_fp.transpose() * df;
        ROS_INFO("[TG] cal dp_best done, size: %ld, %ld", dp_best.rows(), dp_best.cols());
        VectorXd dfp;
        dfp.resize(d_num * (m+1));
        dfp << df,
               dp_best;
        ROS_INFO("[TG] cal dfp done, size: %ld, %ld", dfp.rows(), dfp.cols());
        VectorXd coeff = M.inverse() * Ct * dfp;
        ROS_INFO("[TG] cal coeff done, size: %ld, %ld", coeff.rows(), coeff.cols());
        std::cout << coeff << std::endl << std::endl;
        std::cout << "printing coeffs " << std::endl;
        for (int j = 0; j < m; ++j) {
            PolyCoeff.block(j, idx*p_num1d, 1, p_num1d) = coeff.segment(j*p_num1d, p_num1d).transpose();
            std::cout << "seg [" << j << "]" << std::endl;
            std::cout << PolyCoeff.block(j, idx*p_num1d, 1, p_num1d) << std::endl << std::endl;
        }
    }
    return PolyCoeff;
}

void TrajectoryGeneratorWaypoint::GenContinuityConstraint(
        const int n_seq, 
        const int p_order, 
        const Eigen::VectorXd &Time, 
        const int order_deri, 
        Eigen::MatrixXd& Aeq_con) {
    int n_all_poly = n_seq * (p_order+1);
    int k = order_deri;
    int p_num1d = p_order + 1;
    Aeq_con = MatrixXd::Zero(n_seq-1, n_all_poly);
    for (int j = 0; j < n_seq-1; j++) {
        for (int i = k; i <= p_order; ++i) {
            Aeq_con(j, j*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(j), i-k);
        }
        Aeq_con(j,(j+1)*p_num1d+k) = - Factorial(k);
    }
}
