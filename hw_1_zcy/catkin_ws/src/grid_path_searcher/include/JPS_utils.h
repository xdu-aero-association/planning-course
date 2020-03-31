#ifndef _JPS_UTILS_H_
#define _JPS_UTILS_H_

#include <iostream>
/// Search and prune neighbors for JPS 3D
struct JPS3DNeib {
  // for each (dx,dy,dz) these contain:
  //    ns: neighbors that are always added
  //    f1: forced neighbors to check (for black obs nodes)
  //    f2: neighbors to add if f1 is forced (for pink nodes)
  // note: the paper of 3d JPS's neighbors do not consider symmetry, and below we should consider symmetry
  //       (so the number to check is larger) !!!
  // 1st dim: there are 26+1(self) directions that may be the source node direction
  // 2nd dim: there are three dims of coords
  // 3rd dim: each explained
  // note: 1st stored as (idx = i + j * 3 + k * 3^2)
  int ns[27][3][26]; // 3rd: max 26 node always added
  int f1[27][3][12]; // 3rd: max 8 forced to check 
  int f2[27][3][12]; // 3rd: max 12 to add if forced
  // nsz contains the number of neighbors for the four different types of moves:
  // no move (norm 0):        26 neighbors always added
  //                          0 forced neighbors to check (never happens)
  //                          0 neighbors to add if forced (never happens)
  // straight (norm 1):       1 neighbor always added
  //                          8 forced neighbors to check
  //                          8 neighbors to add if forced
  // diagonal (norm sqrt(2)): 3 neighbors always added
  //                          8 forced neighbors to check
  //                          12 neighbors to add if forced
  // diagonal (norm sqrt(3)): 7 neighbors always added
  //                          6 forced neighbors to check
  //                          12 neighbors to add if forced
  static constexpr int nsz[4][2] = {{26, 0}, {1, 8}, {3, 12}, {7, 12}};
  JPS3DNeib();

private:
  void Neib(int dx, int dy, int dz, int norm1, int dev, int &tx, int &ty,
            int &tz);
  void FNeib(int dx, int dy, int dz, int norm1, int dev, int &fx, int &fy,
             int &fz, int &nx, int &ny, int &nz);
};

#endif