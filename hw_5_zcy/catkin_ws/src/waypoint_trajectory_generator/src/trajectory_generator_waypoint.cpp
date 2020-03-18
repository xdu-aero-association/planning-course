#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

// ooqp related headers
#include <ooqp/QpGenData.h>
#include <ooqp/QpGenVars.h>
#include <ooqp/QpGenResiduals.h>
#include <ooqp/GondzioSolver.h>
#include <ooqp/QpGenSparseMa27.h>

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
    int d_num = d_order;  // dimension of one point for close form sol. (>=3)

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);
    ROS_WARN("[TG] d_order: %d", d_order);
    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    for (int idx = 0; idx < Path.cols(); ++idx) {
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
        
    /*   Produce the derivatives in X, Y and Z axis directly.  */
        MatrixXd M = MatrixXd::Zero(2*(d_num)*m, p_num1d*m);
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < d_order; ++k) {
                M(j*p_num1d+k, j*p_num1d+k) = Factorial(k);
                for (int i = k; i <= p_order; ++i) {
                    M(j*p_num1d+d_num+k, j*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(j), i-k);
                }
            }
        }
        ROS_INFO("[TG] cal M done, size: %ld, %ld", M.rows(), M.cols());
        // std::cout << M << std::endl << std::endl;
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
        // std::cout << Ct << std::endl << std::endl;
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
        // std::cout << Q << std::endl << std::endl;
        // compute R
        MatrixXd R = Ct.transpose() * (M.inverse()).transpose() * Q * M.inverse() * Ct;
        ROS_INFO("[TG] cal R done, size: %ld, %ld", R.rows(), R.cols());
        // std::cout << R << std::endl << std::endl;
        // compute df
        VectorXd df = VectorXd::Zero(2*d_num+m-1);
        df.segment(0,d_num) = start_pose;
        df.segment(d_num+m-1,d_num) = end_pose;
        for (int j = 0; j < m-1; ++j) {
            df(j+d_num) = Path(j+1, idx);
        }
        ROS_INFO("[TG] cal df done, size: %ld, %ld", df.rows(), df.cols());
        // std::cout << df << std::endl << std::endl;
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
        // std::cout << coeff << std::endl << std::endl;
        std::cout << "printing coeffs " << std::endl;
        for (int j = 0; j < m; ++j) {
            PolyCoeff.block(j, idx*p_num1d, 1, p_num1d) = coeff.segment(j*p_num1d, p_num1d).transpose();
            std::cout << "seg [" << j << "]" << std::endl;
            // std::cout << PolyCoeff.block(j, idx*p_num1d, 1, p_num1d) << std::endl << std::endl;
        }
    }
    return PolyCoeff;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::SolvebyOOQP(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time) {
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment
    int d_num = d_order;  // dimension of one point for close form sol. (>=3)

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    ROS_WARN("[TG] d_order: %d", d_order);
    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    for (int idx = 0; idx < Path.cols(); ++idx) {   

        int nnzA = GetnnzA(m, d_order);
        ROS_INFO("[TG][OOQP] nnzA: %d", nnzA);
        int *irowA = new (std::nothrow) int[nnzA];
        int *jcolA = new (std::nothrow) int[nnzA];
        double *dA = new (std::nothrow) double[nnzA];
        if (!irowA || !jcolA || !dA) {
            // Handle error
            ROS_ERROR("ERROR: memory allocation failed !!!!");
            return PolyCoeff;
        }
        int k_A = 0;
        int row_idx_A = 0;

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
        GenStartConstraint(m, d_order, Time, 0, Aeq_start, &row_idx_A, irowA, jcolA, dA, &k_A);
        beq_start = start_pose;
        ROS_INFO("[TG] cal Aeq_start done, size: %ld, %ld", Aeq_start.rows(), Aeq_start.cols());
        // calculate Aeq_end
        MatrixXd Aeq_end = MatrixXd::Zero(d_num, n_all_poly);
        VectorXd beq_end = VectorXd::Zero(d_num);
        GenEndConstraint(m, d_order, Time, 0, Aeq_end, &row_idx_A, irowA, jcolA, dA, &k_A);
        beq_end = end_pose;
        ROS_INFO("[TG] cal Aeq_end done, size: %ld, %ld", Aeq_end.rows(), Aeq_end.cols());
        // calculate Aeq_wp
        MatrixXd Aeq_wp = MatrixXd::Zero(m-1, n_all_poly);
        VectorXd beq_wp = VectorXd::Zero(m-1);
        GenWPConstraint(m, d_order, Time, 0, Path.block(0, idx, Path.rows(), 1), Aeq_wp, beq_wp, &row_idx_A, irowA, jcolA, dA, &k_A);
        ROS_INFO("[TG] cal Aeq_wp done, size: %ld, %ld", Aeq_wp.rows(), Aeq_wp.cols());
        // calculate pos continuity
        MatrixXd Aeq_con_p = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 0, Aeq_con_p, &row_idx_A, irowA, jcolA, dA, &k_A);
        ROS_INFO("[TG] cal Aeq_con_p done, size: %ld, %ld", Aeq_con_p.rows(), Aeq_con_p.cols());
        // calculate vel continuity
        MatrixXd Aeq_con_v = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 1, Aeq_con_v, &row_idx_A, irowA, jcolA, dA, &k_A);
        ROS_INFO("[TG] cal Aeq_con_v done, size: %ld, %ld", Aeq_con_v.rows(), Aeq_con_v.cols());
        // calculate acc continuity
        MatrixXd Aeq_con_a = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 2, Aeq_con_a, &row_idx_A, irowA, jcolA, dA, &k_A);
        ROS_INFO("[TG] cal Aeq_con_a done, size: %ld, %ld", Aeq_con_a.rows(), Aeq_con_a.cols());
        // calculate jerk continuity
        MatrixXd Aeq_con_j = MatrixXd::Zero(m-1, n_all_poly);
        GenContinuityConstraint(m, d_order, Time, 3, Aeq_con_j, &row_idx_A, irowA, jcolA, dA, &k_A);
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
        ROS_INFO("[TG][OOQP] row_idx_A: %d", row_idx_A);

        MatrixXd Q = MatrixXd::Zero(m*p_num1d, m*p_num1d);
        for (int j = 0; j < m; ++j) {
            for (int k = 4; k <= p_order; ++k) {
                for (int l = 4; l <= p_order; ++l) {
                    Q(k+j*p_num1d, l+j*p_num1d) = k*(k-1)*(k-1)*(k-3)*l*(l-1)*(l-2)*(l-3)*pow(Time(j), k+l-7)/(k+l-7);
                }
            }
        }
        ROS_INFO("[TG] cal Q done, size: %ld, %ld", Q.rows(), Q.cols());
        // std::cout << Q << std::endl << std::endl;

        int my = beq.rows();
        double *b = new (std::nothrow) double[my];
        if (!b) {
            // Handle error
            ROS_ERROR("ERROR: memory allocation failed !!!!");
            return PolyCoeff;
        }
        for (int i = 0; i < my; ++i) {
            b[i] = beq(i);
        }

        delete [] b;
        delete [] irowA;
        delete [] jcolA;
        delete [] dA;

    }
    
    return PolyCoeff;

}
void TrajectoryGeneratorWaypoint::GenStartConstraint(
        const int n_seq, 
        const int d_order, 
        const Eigen::VectorXd &Time, 
        const int order_deri, 
        Eigen::MatrixXd& Aeq_start,
        int *p_row_idx_A,
        int *irowA,
        int *jcolA,
        double *dA,
        int *p_k_A) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int d_num      = d_order;  // dimension of one point for close form sol. (>=3)
    Aeq_start = MatrixXd::Zero(d_num, n_all_poly);        
    for (int k = 0; k < d_num; ++k) {
        Aeq_start(k,k) = Factorial(k);
        AddtodA(Aeq_start(k,k), k, k, *p_row_idx_A, irowA, jcolA, dA, p_k_A);
    }
    *p_row_idx_A += d_order;
    ROS_INFO("[TG][OOQP] row_idx_A: %d, k_A: %d", *p_row_idx_A, *p_k_A);
}

void TrajectoryGeneratorWaypoint::GenEndConstraint(
        const int n_seq, 
        const int d_order, 
        const Eigen::VectorXd &Time, 
        const int order_deri, 
        Eigen::MatrixXd& Aeq_end,
        int *p_row_idx_A,
        int *irowA,
        int *jcolA,
        double *dA,
        int *p_k_A) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int d_num      = d_order;  // dimension of one point for close form sol. (>=3)
    Aeq_end = MatrixXd::Zero(d_num, n_all_poly);
    for (int k = 0; k < d_num; ++k) {
        for (int i = k; i <= p_order; ++i) {
            Aeq_end(k, (n_seq-1)*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(n_seq-1), i-k);
            AddtodA(Aeq_end(k, (n_seq-1)*p_num1d+i), k, (n_seq-1)*p_num1d+i, *p_row_idx_A, irowA, jcolA, dA, p_k_A);
        }
    }
    *p_row_idx_A += d_order;
    ROS_INFO("[TG][OOQP] row_idx_A: %d, k_A: %d", *p_row_idx_A, *p_k_A);
}

void TrajectoryGeneratorWaypoint::GenWPConstraint(
        const int n_seq, 
        const int d_order, 
        const Eigen::VectorXd &Time, 
        const int order_deri, 
        const Eigen::VectorXd& Path,
        Eigen::MatrixXd& Aeq_wp,
        Eigen::VectorXd& beq_wp,
        int *p_row_idx_A,
        int *irowA,
        int *jcolA,
        double *dA,
        int *p_k_A) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int d_num      = d_order;  // dimension of one point for close form sol. (>=3)
    Aeq_wp = MatrixXd::Zero(n_seq-1, n_all_poly);
    beq_wp = VectorXd::Zero(n_seq-1);
    for (int j = 0; j < n_seq-1; ++j) {
        Aeq_wp(j, (j+1)*p_num1d) = 1;
        beq_wp(j) = Path(j+1);
        AddtodA(Aeq_wp(j, (j+1)*p_num1d), j, (j+1)*p_num1d, *p_row_idx_A, irowA, jcolA, dA, p_k_A);
    }
    *p_row_idx_A += n_seq-1;
    ROS_INFO("[TG][OOQP] row_idx_A: %d, k_A: %d", *p_row_idx_A, *p_k_A);
} 

void TrajectoryGeneratorWaypoint::GenContinuityConstraint(
        const int n_seq, 
        const int d_order, 
        const Eigen::VectorXd &Time, 
        const int order_deri, 
        Eigen::MatrixXd& Aeq_con,
        int *p_row_idx_A,
        int *irowA,
        int *jcolA,
        double *dA,
        int *p_k_A) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int k          = order_deri;
    Aeq_con = MatrixXd::Zero(n_seq-1, n_all_poly);
    for (int j = 0; j < n_seq-1; j++) {
        for (int i = k; i <= p_order; ++i) {
            Aeq_con(j, j*p_num1d+i) = Factorial(i) / Factorial(i-k) * pow(Time(j), i-k);
            AddtodA(Aeq_con(j, j*p_num1d+i), j, j*p_num1d+i, *p_row_idx_A, irowA, jcolA, dA, p_k_A);
        }
        Aeq_con(j,(j+1)*p_num1d+k) = - Factorial(k);
        AddtodA(Aeq_con(j,(j+1)*p_num1d+k), j, (j+1)*p_num1d+k, *p_row_idx_A, irowA, jcolA, dA, p_k_A);
    }
    *p_row_idx_A += n_seq-1;
    ROS_INFO("[TG][OOQP] row_idx_A: %d, k_A: %d", *p_row_idx_A, *p_k_A);
}

void TrajectoryGeneratorWaypoint::AddtodA(
        const double value, 
        const int irow,
        const int jcol,
        const int row_idx_A,
        int *irowA,
        int *jcolA,
        double *dA,
        int *p_k_A) {
    irowA[*p_k_A] = irow + row_idx_A;
    jcolA[*p_k_A] = jcol;
    dA[*p_k_A] = value;
    *p_k_A += 1;
}

int TrajectoryGeneratorWaypoint::GetnnzA(const int n_seq, const int d_order) {
    int p_order    = 2 * d_order - 1;              // the order of polynomial
    int p_num1d    = p_order + 1;
    int n_all_poly = n_seq * p_num1d;
    int d_num      = d_order;  // dimension of one point for close form sol. (>=3)
    int nnzA = 0;
    // start
    for (int k = 0; k < d_num; ++k) {
        nnzA += 1;
    }
    // end
    for (int k = 0; k < d_num; ++k) {
        for (int i = k; i <= p_order; ++i) {
            nnzA += 1;
        }
    }
    // wp
    for (int j = 0; j < n_seq-1; ++j) {
        nnzA += 1;
    }
    // continuity
    for (int k = 0; k < 4; ++k) {
        for (int j = 0; j < n_seq-1; j++) {
            for (int i = k; i <= p_order; ++i) {
                nnzA += 1;
            }
            nnzA += 1;
        }
    }
    return nnzA;
}
