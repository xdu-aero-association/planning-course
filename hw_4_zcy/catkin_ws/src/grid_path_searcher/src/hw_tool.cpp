#include <hw_tool.h>

using namespace std;
using namespace Eigen;

void Homeworktool::initGridMap(double _resolution, Vector3d global_xyz_l,
                               Vector3d global_xyz_u, int max_x_id,
                               int max_y_id, int max_z_id) {
  gl_xl = global_xyz_l(0);
  gl_yl = global_xyz_l(1);
  gl_zl = global_xyz_l(2);

  gl_xu = global_xyz_u(0);
  gl_yu = global_xyz_u(1);
  gl_zu = global_xyz_u(2);

  GLX_SIZE = max_x_id;
  GLY_SIZE = max_y_id;
  GLZ_SIZE = max_z_id;
  GLYZ_SIZE = GLY_SIZE * GLZ_SIZE;
  GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

  resolution = _resolution;
  inv_resolution = 1.0 / _resolution;

  data = new uint8_t[GLXYZ_SIZE];
  memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
}

void Homeworktool::setObs(const double coord_x, const double coord_y,
                          const double coord_z) {
  if (coord_x < gl_xl || coord_y < gl_yl || coord_z < gl_zl ||
      coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu)
    return;

  int idx_x = static_cast<int>((coord_x - gl_xl) * inv_resolution);
  int idx_y = static_cast<int>((coord_y - gl_yl) * inv_resolution);
  int idx_z = static_cast<int>((coord_z - gl_zl) * inv_resolution);

  data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

bool Homeworktool::isObsFree(const double coord_x, const double coord_y,
                             const double coord_z) {
  Vector3d pt;
  Vector3i idx;

  pt(0) = coord_x;
  pt(1) = coord_y;
  pt(2) = coord_z;
  idx = coord2gridIndex(pt);

  int idx_x = idx(0);
  int idx_y = idx(1);
  int idx_z = idx(2);

  return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE &&
          idx_z >= 0 && idx_z < GLZ_SIZE &&
          (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

Vector3d Homeworktool::gridIndex2coord(const Vector3i &index) {
  Vector3d pt;

  pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
  pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
  pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

  return pt;
}

Vector3i Homeworktool::coord2gridIndex(const Vector3d &pt) {
  Vector3i idx;
  idx << min(max(int((pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
      min(max(int((pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
      min(max(int((pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);

  return idx;
}

Eigen::Vector3d Homeworktool::coordRounding(const Eigen::Vector3d &coord) {
  return gridIndex2coord(coord2gridIndex(coord));
}

double Homeworktool::OptimalBVP(Eigen::Vector3d _start_position,
                                Eigen::Vector3d _start_velocity,
                                Eigen::Vector3d _target_position) {
  double optimal_cost =
      100000; // this just to initial the optimal_cost, you can delete it
  /*




  STEP 2: go to the hw_tool.cpp and finish the function Homeworktool::OptimalBVP
  the solving process has been given in the document

  because the final point of trajectory is the start point of OBVP, so we input
  the pos,vel to the OBVP

  after finish Homeworktool::OptimalBVP, the Trajctory_Cost will record the
  optimal cost of this trajectory


  */
  ROS_INFO("p0:[%.3f,%.3f,%.3f] pf:[%.3f,%.3f,%.3f]", _start_position(0),
           _start_position(1), _start_position(2), _target_position(0),
           _target_position(1), _target_position(2));
  ROS_INFO("v0:[%.3f,%.3f,%.3f", _start_velocity(0), _start_velocity(1),
           _start_velocity(2));
  // first solve the 3rd order polynomial of T
  double px0 = _start_position(0), py0 = _start_position(1),
         pz0 = _start_position(2);
  double vx0 = _start_velocity(0), vy0 = _start_velocity(1),
         vz0 = _start_velocity(2);
  double pxf = _target_position(0), pyf = _target_position(1),
         pzf = _target_position(2);
  VectorXd coeffs(5);
  coeffs[0] = -36 * pow(px0, 2) + 72 * px0 * pxf - 36 * pow(pxf, 2) -
              36 * pow(py0, 2) + 72 * py0 * pyf - 36 * pow(pyf, 2) -
              36 * pow(pz0, 2) + 72 * pz0 * pzf - 36 * pow(pzf, 2);
  coeffs[1] = 24 * pxf * vx0 - 24 * px0 * vx0 - 24 * py0 * vy0 +
              24 * pyf * vy0 - 24 * pz0 * vz0 + 24 * pzf * vz0;
  coeffs[2] = -4 * pow(vx0, 2) - 4 * pow(vy0, 2) - 4 * pow(vz0, 2);
  coeffs[3] = 0; // why?? in matlab coeffs has only 4 members
  coeffs[4] = 1;
  ROS_INFO("[OBVP] coeffs initialized ");
  // method 1: solve by compute eigen values of matrix A
  MatrixXd A(4, 4);
  for (int i = 0; i < 4; ++i) {
    if (i < 3) {
      A(i + 1, i) = 1;
    }
    A(i, 3) = -coeffs[i];
  }
  ROS_INFO("[OBVP] matrix A initialized ");
  Eigen::EigenSolver<MatrixXd> es(A);
  auto roots1 = es.eigenvalues();
  ROS_INFO("[OBVP] eigen solved");
  // method 2: solve by eigen polynomialsolver
  Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
  solver.compute(coeffs);
  const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &roots2 =
      solver.roots();
  ROS_INFO("[OBVP] polynomial solved");
  bool first_flag = 0;
  double T = 1;
  for (int r = 0; r < roots1.rows(); ++r) {
    auto root_eigen = roots1(r);
    if (abs(root_eigen.imag()) > 1e-10 || root_eigen.real() <= 0) {
      continue;
    }
    ROS_WARN("[OBVP] r:%d, root_eigen:%.3f", r, root_eigen.real());
  }
  for (int r = 0; r < roots2.rows(); ++r) {
    //   for (int r = 0; r < 4; ++r) {
    auto root_poly = roots2(r);
    // if (abs(root_eigen.imag()) > 1e-10 || root_eigen.real() <= 0) {
    if (abs(root_poly.imag()) > 1e-10 || root_poly.real() <= 0) {
      continue;
    }
    // T = root_eigen.real();
    T = root_poly.real();
    // int z = r + 1;
    // T = pow(z, 4) -
    //     pow(z, 2) * (4 * pow(vx0, 2) + 4 * pow(vy0, 2) + 4 * pow(vz0, 2)) -
    //     z * (24 * px0 * vx0 - 24 * pxf * vx0 + 24 * py0 * vy0 - 24 * pyf *
    //     vy0 +
    //          24 * pz0 * vz0 - 24 * pzf * vz0) +
    //     72 * pz0 * pzf + 72 * py0 * pyf + 72 * px0 * pxf - 36 * pow(pzf, 2) -
    //     36 * pow(pz0, 2) - 36 * pow(pyf, 2) - 36 * pow(py0, 2) -
    //     36 * pow(pxf, 2) - 36 * pow(px0, 2);
    ROS_WARN("[OBVP] r:%d, root_poly:%.3f", r, root_poly.real());
    // ROS_INFO("[OBVP] z: %d, T:%.3f",z,T);
    // if (T < 0) {
    //   continue;
    // }
    // calculate J for given T
    VectorXd delta_pose(6, 1); // [dpx,dpy,dpz,dvx,dvy,dvz].T
    for (int i = 0; i < 3; ++i) {
      delta_pose(i) =
          _target_position(i) - _start_velocity(i) * T - _start_position(i);
      delta_pose(i + 3) =
          0.0 - _start_velocity(i); // because we want to stay at final pos
    }
    MatrixXd B(6, 6);
    for (int i = 0; i < 3; ++i) {
      B(i, i) = -12 / pow(T, 3);
      B(i, i + 3) = 6 / pow(T, 2);
      B(i + 3, i) = 6 / pow(T, 2);
      B(i + 3, i + 3) = -2 / T;
    }
    // ROS_INFO("[OBVP] matrix B initialized");
    VectorXd Result = B * delta_pose;
    printf("[OBVP] Result: \n");
    for (int i = 0; i < 6; ++i) {
        printf("R[%d]:%.3f ", i, Result(i));
    }
    printf("\n");
    // ROS_INFO("[OBVP] result calculated");
    double a[3] = {0.0};
    double b[3] = {0.0};
    for (int i = 0; i < 3; ++i) {
      a[i] = Result(i);
      b[i] = Result(i + 3);
    }
    a[0] = (12*(px0-pxf+T*vx0))/(T*T*T) - (6*vx0)/(T*T);
    a[1] = (12*(py0-pyf+T*vy0))/(T*T*T) - (6*vy0)/(T*T);
    a[2] = (12*(pz0-pzf+T*vz0))/(T*T*T) - (6*vz0)/(T*T);
    b[0] = (2*vx0)/T - (6*(px0-pxf+T*vx0))/(T*T);
    b[1] = (2*vy0)/T - (6*(py0-pyf+T*vz0))/(T*T);
    b[2] = (2*vz0)/T - (6*(pz0-pzf+T*vy0))/(T*T);
    printf("[OBVP] ab: \n");
    for (int i = 0; i < 3; ++i) {
        printf("a[%d]:%.3f ", i, a[i]);
    }
    for (int i = 0; i < 3; ++i) {
        printf("b[%d]:%.3f ", i, b[i]);
    }
    printf("\n");
    double J = T;
    for (int i = 0; i < 3; ++i) {
      J += (1 / 3 * pow(a[i], 2) * pow(T, 3) + a[i] * b[i] * pow(T, 2) +
            pow(b[i], 2) * T);
    }
    if (first_flag == 0) {
      first_flag = 1;
      optimal_cost = J;
    } else if (optimal_cost > J) {
      optimal_cost = J;
    }
    ROS_WARN("[OBVP] optimal cost: %.3f", optimal_cost);
    if (optimal_cost < 0) {
      optimal_cost = 0;
    }
  }
  return optimal_cost;
}
