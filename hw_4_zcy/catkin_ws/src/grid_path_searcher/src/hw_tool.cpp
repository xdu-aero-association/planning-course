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

OBVPResult Homeworktool::OptimalBVP(Eigen::Vector3d _start_position,
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
  VectorXd coeffs = VectorXd::Zero(5, 1);
  coeffs[0] = -36 * pow(px0, 2) + 72 * px0 * pxf - 36 * pow(pxf, 2) -
              36 * pow(py0, 2) + 72 * py0 * pyf - 36 * pow(pyf, 2) -
              36 * pow(pz0, 2) + 72 * pz0 * pzf - 36 * pow(pzf, 2);
  coeffs[1] = 24 * pxf * vx0 - 24 * px0 * vx0 - 24 * py0 * vy0 +
              24 * pyf * vy0 - 24 * pz0 * vz0 + 24 * pzf * vz0;
  coeffs[2] = -4 * pow(vx0, 2) - 4 * pow(vy0, 2) - 4 * pow(vz0, 2);
  coeffs[3] = 0; // in matlab we have to use [c,t] = coeffs(f(t)) to get
                 // corresponding terms!!
  coeffs[4] = 1;
  printf("coeffs: \n");
  for (int i = 0; i < 5; ++i) {
    printf("[%d]: %.3f ", i, coeffs[i]);
  }
  printf("\n");
  // another way for coeffs
  Vector3d _target_velocity = Vector3d::Zero(3);
  Vector3d param1 = _target_position - _start_position;
  Vector3d param2 = _start_velocity + _target_velocity;
  Vector3d param3 = Vector3d::Zero(3);
  param3(0) =
      pow(_start_velocity(0), 2) + 2 * _start_velocity(0) * _target_velocity(0);
  param3(1) =
      pow(_start_velocity(1), 2) + 2 * _start_velocity(1) * _target_velocity(1);
  param3(2) =
      pow(_start_velocity(2), 2) + 2 * _start_velocity(2) * _target_velocity(2);
  double coef1 = pow(param1(0), 2) + pow(param1(1), 2) + pow(param1(2), 2);
  double coef2 =
      param1(0) * param2(0) + param1(1) * param2(1) + param1(2) * param2(2);
  double coef3 = param3(0) + param3(1) + param3(2);
  VectorXd coeffs_an = VectorXd::Zero(5, 1);
  coeffs_an[0] = -36.0 * coef1;
  coeffs_an[1] = 24.0 * coef2;
  coeffs_an[2] = -4.0 * coef3;
  coeffs_an[3] = 0.0;
  coeffs_an[4] = 1.0;
  printf("coeffs_an: \n");
  for (int i = 0; i < 5; ++i) {
    printf("[%d]: %.3f ", i, coeffs_an[i]);
  }
  printf("\n");

  ROS_INFO("[OBVP] coeffs initialized ");
  // method 1: solve by compute eigen values of matrix A
  MatrixXd A = MatrixXd::Zero(4, 4);
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
  double optimal_T = 1.0;
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
    VectorXd delta_pose = VectorXd::Zero(6, 1); // [dpx,dpy,dpz,dvx,dvy,dvz].T
    for (int i = 0; i < 3; ++i) {
      delta_pose(i) =
          _target_position(i) - _start_velocity(i) * T - _start_position(i);
      delta_pose(i + 3) =
          0.0 - _start_velocity(i); // because we want to stay at final pos
    }
    MatrixXd B = MatrixXd::Zero(6, 6);
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
    a[0] = (12 * (px0 - pxf + T * vx0)) / (T * T * T) - (6 * vx0) / (T * T);
    a[1] = (12 * (py0 - pyf + T * vy0)) / (T * T * T) - (6 * vy0) / (T * T);
    a[2] = (12 * (pz0 - pzf + T * vz0)) / (T * T * T) - (6 * vz0) / (T * T);
    b[0] = (2 * vx0) / T - (6 * (px0 - pxf + T * vx0)) / (T * T);
    b[1] = (2 * vy0) / T - (6 * (py0 - pyf + T * vy0)) / (T * T);
    b[2] = (2 * vz0) / T - (6 * (pz0 - pzf + T * vz0)) / (T * T);
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
      // J += (1 / 3 * pow(a[i], 2) * pow(T, 3) + a[i] * b[i] * pow(T, 2) +
      //       pow(b[i], 2) * T);
      double temp = (1.0 / 3.0 * pow(a[i], 2) * pow(T, 3) +
                     a[i] * b[i] * pow(T, 2) + pow(b[i], 2) * T);
      double temp_1 = 1.0 / 3.0 * a[i] * a[i] * T * T * T + a[i] * b[i] * T * T +
                      b[i] * b[i] * T;
      ROS_INFO("[OBVP] temp:%.3f, temp_1:%.3f", temp, temp_1);
      J += 1.0 / 3.0 * a[i] * a[i] * T * T * T + a[i] * b[i] * T * T +
           b[i] * b[i] * T;
    }
    double J_ex =
        T + 12 * coef1 / pow(T, 3) - 12 * coef2 / pow(T, 2) + 4 * coef3 / T;
    double J_ex1 = (T * T * T * a[0] * a[0]) / 3 +
                   (T * T * T * a[1] * a[1]) / 3 +
                   (T * T * T * a[2] * a[2]) / 3 + T * T * a[0] * b[0] +
                   T * T * a[1] * b[1] + T * T * a[2] * b[2] + T * b[0] * b[0] +
                   T * b[1] * b[1] + T * b[2] * b[2] + T;
    ROS_WARN("[OBVP] original J: %.3f, J_ex: %.3f, J_ex1: %.3f", J, J_ex,
             J_ex1);
    J = J_ex; // J shouldn't < 0 !!
    if (J < 0.0) {
      continue;
    }
    if (first_flag == 0) {
      first_flag = 1;
      optimal_cost = J;
      optimal_T = T;
    } else if (optimal_cost > J) {
      optimal_cost = J;
      optimal_T = T;
    }
    ROS_WARN("[OBVP] optimal cost: %.3f, T: %.3f", optimal_cost, optimal_T);
    // if (optimal_cost < 0) {
    //   optimal_cost = 0;
    // }
  }
  return OBVPResult(optimal_cost, optimal_T);
}

vector<Vector3d> Homeworktool::reconstructOptimalTrajactory(
    const OBVPResult &obvp_result, Eigen::Vector3d _start_position,
    Eigen::Vector3d _start_velocity, Eigen::Vector3d _target_position) {
  double T = obvp_result.optimal_T;
  vector<Vector3d> optimal_states;
  optimal_states.clear();
  const double min_delta_time = 0.01;
  const int max_num = 30;
  int num = static_cast<int>(T / min_delta_time);
  double delta_time = min_delta_time;
  if (num > max_num) {
    num = max_num;
    delta_time = static_cast<double>(T / max_num);
  }
  optimal_states.resize(num);
  ROS_INFO("[RECON] optimal_T:%.3f, optimal_cost:%.3f ", T,
           obvp_result.optimal_cost);
  ROS_INFO("[RECON] delta_time:%.3f, num:%d", delta_time, num);
  VectorXd delta_pose = VectorXd::Zero(6, 1); // [dpx,dpy,dpz,dvx,dvy,dvz].T
  for (int i = 0; i < 3; ++i) {
    delta_pose(i) =
        _target_position(i) - _start_velocity(i) * T - _start_position(i);
    delta_pose(i + 3) =
        0.0 - _start_velocity(i); // because we want to stay at final pos
  }
  MatrixXd B = MatrixXd::Zero(6, 6);
  for (int i = 0; i < 3; ++i) {
    B(i, i) = -12 / pow(T, 3);
    B(i, i + 3) = 6 / pow(T, 2);
    B(i + 3, i) = 6 / pow(T, 2);
    B(i + 3, i + 3) = -2 / T;
  }
  VectorXd Result = B * delta_pose;
  printf("[RECON] Result: \n");
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
  double px0 = _start_position(0), py0 = _start_position(1),
         pz0 = _start_position(2);
  double vx0 = _start_velocity(0), vy0 = _start_velocity(1),
         vz0 = _start_velocity(2);
  double pxf = _target_position(0), pyf = _target_position(1),
         pzf = _target_position(2);
  a[0] = (12 * (px0 - pxf + T * vx0)) / (T * T * T) - (6 * vx0) / (T * T);
  a[1] = (12 * (py0 - pyf + T * vy0)) / (T * T * T) - (6 * vy0) / (T * T);
  a[2] = (12 * (pz0 - pzf + T * vz0)) / (T * T * T) - (6 * vz0) / (T * T);
  b[0] = (2 * vx0) / T - (6 * (px0 - pxf + T * vx0)) / (T * T);
  b[1] = (2 * vy0) / T - (6 * (py0 - pyf + T * vy0)) / (T * T);
  b[2] = (2 * vz0) / T - (6 * (pz0 - pzf + T * vz0)) / (T * T);
  printf("[RECON] ab: \n");
  for (int i = 0; i < 3; ++i) {
    printf("a[%d]:%.3f ", i, a[i]);
  }
  for (int i = 0; i < 3; ++i) {
    printf("b[%d]:%.3f ", i, b[i]);
  }
  printf("\n");
  double J = T;
  for (int i = 0; i < 3; ++i) {
    J += (1.0 / 3.0 * pow(a[i], 2) * pow(T, 3) + a[i] * b[i] * pow(T, 2) +
          pow(b[i], 2) * T);
  }
  double J_ex = (T * T * T * a[0] * a[0]) / 3 + (T * T * T * a[1] * a[1]) / 3 +
                (T * T * T * a[2] * a[2]) / 3 + T * T * a[0] * b[0] +
                T * T * a[1] * b[1] + T * T * a[2] * b[2] + T * b[0] * b[0] +
                T * b[1] * b[1] + T * b[2] * b[2] + T;
  ROS_INFO("[RECON] actual cost: %.3f, ex cost: %.3f", J, J_ex);
  // calculate states
  double time = 0;
  for (auto &state : optimal_states) {
    for (int i = 0; i < 3; ++i) {
      state(i) = 1.0 / 6.0 * a[i] * pow(time, 3) + 1.0 / 2.0 * b[i] * pow(time, 2) +
                 _start_velocity(i) * time + _start_position(i);
    }
    // ROS_INFO("[RECON] origin  state: %.3f %.3f %.3f", state(0), state(1),
    //          state(2));
    // state(0) = a[0] * time * time * time / 6 + b[0] * time * time / 2 +
    //            vx0 * time + px0;
    // state(1) = a[1] * time * time * time / 6 + b[1] * time * time / 2 +
    //            vy0 * time + py0;
    // state(2) = a[2] * time * time * time / 6 + b[2] * time * time / 2 +
    //            vz0 * time + pz0;
    // ROS_INFO("[RECON] another state: %.3f %.3f %.3f", state(0), state(1),
    //          state(2));
    time += delta_time;
  }
  auto end_state = optimal_states.back();
  ROS_INFO("[RECON] end state: %.3f %.3f %.3f", end_state(0), end_state(1),
           end_state(2));
  return optimal_states;
}
