#ifndef _ASTART_SEARCHER_H
#define _ASTART_SEARCHER_H

#include "backward.hpp"
#include "node.h"
#include <Eigen/Eigen>
#include <iostream>
#include <ros/console.h>
#include <ros/ros.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <string>

class AstarPathFinder {
private:
protected:
  uint8_t *data; // obstacles
  GridNodePtr ***GridNodeMap;
  Eigen::Vector3i goalIdx;
  int GLX_SIZE, GLY_SIZE, GLZ_SIZE;
  int GLXYZ_SIZE, GLYZ_SIZE;

  double resolution, inv_resolution;
  double gl_xl, gl_yl, gl_zl;
  double gl_xu, gl_yu, gl_zu;

  // static auto key_comp = [] (double &f1, double &f2) {
  //   if (f1 == f2 && tie_breaker_ == 1) {
  //     srand (time(NULL));
  //     f1 += double(rand() % 1000) / 1000.0 / 10000.0;
  //     return f1 < f2;
  //   }
  //   return f1 < f2;
  // }

  GridNodePtr terminatePtr;
  std::multimap<double, GridNodePtr> openSet;

  double getHeu(GridNodePtr node1, GridNodePtr node2);
	double getGScore(GridNodePtr node1, GridNodePtr node2);
  void AstarGetSucc(GridNodePtr currentPtr,
                    std::vector<GridNodePtr> &neighborPtrSets,
                    std::vector<double> &edgeCostSets);

  bool isOccupied(const int &idx_x, const int &idx_y, const int &idx_z) const;
  bool isOccupied(const Eigen::Vector3i &index) const;
  bool isFree(const int &idx_x, const int &idx_y, const int &idx_z) const;
  bool isFree(const Eigen::Vector3i &index) const;

  Eigen::Vector3d gridIndex2coord(const Eigen::Vector3i &index);
  Eigen::Vector3i coord2gridIndex(const Eigen::Vector3d &pt);

  int h_selector_ = 0;
  bool tie_breaker_ = 0;
  std::string name_ = "Astar";

public:
  AstarPathFinder(){};
  ~AstarPathFinder(){};
  void AstarGraphSearch(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);
  void resetGrid(GridNodePtr ptr);
  void resetUsedGrids();

  void initGridMap(double _resolution, Eigen::Vector3d global_xyz_l,
                   Eigen::Vector3d global_xyz_u, int max_x_id, int max_y_id,
                   int max_z_id);
  void setObs(const double coord_x, const double coord_y, const double coord_z);

  Eigen::Vector3d coordRounding(const Eigen::Vector3d &coord);
  std::vector<Eigen::Vector3d> getPath();
  std::vector<Eigen::Vector3d> getVisitedNodes();

  double TieBreaker(const double fs);

  void set_h_type(const int type) {h_selector_ = type;}
  void set_tie_breaker(const bool type) {tie_breaker_ = type;}
  std::string get_name() const {return name_;}
};

#endif