///////////////////////////////////////////////////////////////////////////////
// Authors: Ilgweon Kang and Lutong Wang
//          (respective Ph.D. advisors: Chung-Kuan Cheng, Andrew B. Kahng),
//          based on Dr. Jingwei Lu with ePlace and ePlace-MS
//
//          Many subsequent improvements were made by Mingyu Woo
//          leading up to the initial release.
//
// BSD 3-Clause License
//
// Copyright (c) 2018, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#include <replace_private.h>
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <ctime>

#include "replace_private.h"
#include "initPlacement.h"
#include "wlen.h"
#include "plot.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <vector>

using std::vector;
using std::string;
using std::to_string;
using std::max;
using std::min;

void initial_placement() {
  using namespace Eigen;
  printf("PROC:  Conjugate Gradient (CG) method to obtain the IP\n");

  int itmax = 100;

  prec target_tol = 0.000001;
  prec x_err = 0, y_err = 0;

  double time_s = 0;

  auto hpwl = GetUnscaledHpwl(); 

  if(isSkipIP) {
    return;
  }

  printf("INFO:  The Initial HPWL is %.6lf\n", hpwl.first + hpwl.second);

  if(hpwl.first + hpwl.second <= 0) {
    printf("ERROR: HPWL <= 0, skip initial QP\n");
    return;
  }

  printf("INFO:  The Matrix Size is %d\n", moduleCNT);
  fflush(stdout);

  // malloc to solve PCG
  //
  // Ax = b
  //
  // x : variable vector to solve
  // b : constant vector
  //

  setNbThreads(numThread);

  // BCGSTAB settings
  SMatrix eMatX(moduleCNT, moduleCNT), eMatY(moduleCNT, moduleCNT);

  VectorXf xcg_x(moduleCNT), xcg_b(moduleCNT), ycg_x(moduleCNT),
      ycg_b(moduleCNT);
  // cout<<"IP loop start"<<endl;
  for(int i = 0;; i++) {
    if(i >= numInitPlaceIter) {
      break;
    }
    // cout<<"IP loop debug 0"<<endl;
    time_start(&time_s);
    CreateSparseMatrix(xcg_x, xcg_b, ycg_x, ycg_b, eMatX, eMatY);
    // cout<<"IP loop debug 1"<<endl;
    BiCGSTAB< SMatrix, IdentityPreconditioner > solver;
    solver.setMaxIterations(itmax);

    solver.compute(eMatX);
    xcg_x = solver.solveWithGuess(xcg_b, xcg_x);
    x_err = solver.error();
    // cout<<"IP loop debug 2"<<endl;
    solver.compute(eMatY);
    ycg_x = solver.solveWithGuess(ycg_b, ycg_x);
    y_err = solver.error();
    // cout<<"IP loop debug 3"<<endl;
    update_module(xcg_x, ycg_x);
    // cout<<"IP loop debug 4"<<endl;
    update_pin_by_module();
    // cout<<"IP loop debug 5"<<endl;
    update_net_by_pin();
    // cout<<"IP loop debug 6"<<endl;
    if(isPlot && i % 5 == 0) {
      SaveCellPlotAsJPEG(string("FIP - Iter: ") + to_string(i), false,
                         string(dir_bnd) + string("/initPlace/initPlacement_") +
                             intoFourDigit(i));
      // SavePlot( string("FIP - Iter: ") + to_string(i) );
    }
    prec mean = 1.0;
    prec std = 0.5;
    std::normal_distribution<> nd(mean, std);

    std::random_device rd;
    std::mt19937 gen(rd());
    // for(int i = 0;i < 16;i++) {
    //   xcg_x[i] *= nd(gen);
    //   ycg_x[i] *= nd(gen);
    // }
    update_module(xcg_x, ycg_x);
    time_end(&time_s);

    hpwl = GetUnscaledHpwl();

    printf("INFO:  IP%3d,  CG Error %.6lf,  HPWL %.6lf,  CPUtime %.2lf\n", i,
           max(x_err, y_err), hpwl.first + hpwl.second, time_s);
    fflush(stdout);

    if(fabs(x_err) < target_tol && fabs(y_err) < target_tol && i > 4) {
      break;
    }
  }
}

void build_data_struct(bool initCoordi) {
  MODULE *mdp = NULL;
  TERM *term = NULL;
  PIN *pin = NULL;
  NET *curNet = NULL;
  FPOS pof;

  prec min_x = 0;
  prec min_y = 0;
  prec max_x = 0;
  prec max_y = 0;

  for(int i = 0; i < moduleCNT; i++) {
    mdp = &moduleInstance[i];

    if(initCoordi) {
      mdp->center = place.center;
    }

    mdp->pmin.x = mdp->center.x - 0.5 * mdp->size.x;
    mdp->pmin.y = mdp->center.y - 0.5 * mdp->size.y;

    mdp->pmax.x = mdp->center.x + 0.5 * mdp->size.x;
    mdp->pmax.y = mdp->center.y + 0.5 * mdp->size.y;

    for(int j = 0; j < mdp->pinCNTinObject; j++) {
      pof = mdp->pof[j];
      pin = mdp->pin[j];

      pin->fp.x = mdp->center.x + pof.x;
      pin->fp.y = mdp->center.y + pof.y;

      pin->X_MIN = 0;
      pin->X_MAX = 0;
      pin->Y_MIN = 0;
      pin->Y_MAX = 0;
    }
  }

  term_pmin.x = PREC_MAX;
  term_pmin.x = PREC_MAX;
  term_pmax.x = 0;
  term_pmax.y = 0;

  for(int i = 0; i < terminalCNT; i++) {
    term = &terminalInstance[i];

    term->pmin.x = term->center.x - 0.5 * term->size.x;
    term->pmin.y = term->center.y - 0.5 * term->size.y;

    term->pmax.x = term->center.x + 0.5 * term->size.x;
    term->pmax.y = term->center.y + 0.5 * term->size.y;

    if(term_pmin.x > term->pmin.x)
      term_pmin.x = term->pmin.x;
    if(term_pmin.y > term->pmin.y)
      term_pmin.y = term->pmin.y;
    if(term_pmax.x < term->pmax.x)
      term_pmax.x = term->pmax.x;
    if(term_pmax.y < term->pmax.y)
      term_pmax.y = term->pmax.y;

    for(int j = 0; j < term->pinCNTinObject; j++) {
      pof = term->pof[j];
      pin = term->pin[j];

      pin->fp.x = term->center.x + pof.x;
      pin->fp.y = term->center.y + pof.y;

      pin->X_MIN = 0;
      pin->X_MAX = 0;
      pin->Y_MIN = 0;
      pin->Y_MAX = 0;
    }
  }

  for(int i = 0; i < netCNT; i++) {
    curNet = &netInstance[i];

    min_x = PREC_MAX;
    min_y = PREC_MAX;
    max_x = PREC_MIN;
    max_y = PREC_MIN;

    PIN *pin_xmin = NULL;
    PIN *pin_ymin = NULL;
    PIN *pin_xmax = NULL;
    PIN *pin_ymax = NULL;

    for(int j = 0; j < curNet->pinCNTinObject; j++) {
      pin = curNet->pin[j];
      if(pin_xmin) {
        if(min_x > pin->fp.x) {
          min_x = pin->fp.x;
          pin_xmin->X_MIN = 0;
          pin_xmin = pin;
          pin->X_MIN = 1;
        }
      }
      else {
        min_x = pin->fp.x;  // mdp->center.x ;
        pin_xmin = pin;
        pin->X_MIN = 1;
      }

      if(pin_ymin) {
        if(min_y > pin->fp.y)  // mdp->center.y)
        {
          min_y = pin->fp.y;  // mdp->center.y ;
          pin_ymin->Y_MIN = 0;
          pin_ymin = pin;
          pin->Y_MIN = 1;
        }
      }
      else {
        min_y = pin->fp.y;  // mdp->center.y ;
        pin_ymin = pin;
        pin->Y_MIN = 1;
      }

      if(pin_xmax) {
        if(max_x < pin->fp.x)  // mdp->center.x)
        {
          max_x = pin->fp.x;  // mdp->center.x ;
          pin_xmax->X_MAX = 0;
          pin_xmax = pin;
          pin->X_MAX = 1;
        }
      }
      else {
        max_x = pin->fp.x;  // mdp->center.x ;
        pin_xmax = pin;
        pin->X_MAX = 1;
      }

      if(pin_ymax) {
        if(max_y < pin->fp.y)  // mdp->center.y)
        {
          max_y = pin->fp.y;  // mdp->center.y ;
          pin_ymax->Y_MAX = 0;
          pin_ymax = pin;
          pin->Y_MAX = 1;
        }
      }
      else {
        max_y = pin->fp.y;  // mdp->center.y ;
        pin_ymax = pin;
        pin->Y_MAX = 1;
      }
    }

    curNet->min_x = min_x;
    curNet->min_y = min_y;

    curNet->max_x = max_x;
    curNet->max_y = max_y;
  }
}

void update_module(VectorXf &xcg_x, VectorXf &ycg_x) {
  MODULE *mdp = NULL;
  for(int i = 0; i < moduleCNT; i++) {
    mdp = &moduleInstance[i];

    mdp->center.x = xcg_x(i);
    mdp->center.y = ycg_x(i);

    if((mdp->center.x + 0.5 * mdp->size.x) > place.end.x)
      mdp->center.x = place.end.x - 0.5 * mdp->size.x - Epsilon;

    if((mdp->center.y + 0.5 * mdp->size.y) > place.end.y)
      mdp->center.y = place.end.y - 0.5 * mdp->size.y - Epsilon;

    if((mdp->center.x - 0.5 * mdp->size.x) < place.org.x)
      mdp->center.x = place.org.x + 0.5 * mdp->size.x + Epsilon;

    if((mdp->center.y - 0.5 * mdp->size.y) < place.org.y)
      mdp->center.y = place.org.y + 0.5 * mdp->size.y + Epsilon;

    mdp->pmin.x = mdp->center.x - 0.5 * mdp->size.x;
    mdp->pmin.y = mdp->center.y - 0.5 * mdp->size.y;

    mdp->pmax.x = mdp->center.x + 0.5 * mdp->size.x;
    mdp->pmax.y = mdp->center.y + 0.5 * mdp->size.y;
  }
}

void update_pin_by_module(void) {
  FPOS pof;
  PIN *pin = NULL;
  MODULE *mdp = NULL;
  // cout<<"update pin by module debug 0"<<endl;
  for(int i = 0; i < moduleCNT; i++) {
    mdp = &moduleInstance[i];
    // cout<<"update pin by module debug 1"<<endl;
    for(int j = 0; j < mdp->pinCNTinObject; j++) {
      pof = mdp->pof[j];
      pin = mdp->pin[j];
      // cout<<"update pin by module debug 2"<<endl;
      if(pin->moduleID != i || pin->pinIDinModule != j || pin->term == 1) {
        cout<<"module ID = "<<pin->moduleID<<" i = "<<i<<endl;
        cout<<"pin id = "<<pin->pinIDinModule<<" j = "<<j<<endl;
        cout<<"pin->term = "<<pin->term<<endl;; 
        exit(1);
      }

      pin->fp.x = mdp->center.x + pof.x;
      pin->fp.y = mdp->center.y + pof.y;

      pin->X_MIN = 0;
      pin->X_MAX = 0;
      pin->Y_MIN = 0;
      pin->Y_MAX = 0;
      // cout<<"update pin by module debug 3"<<endl;
    }
  }
}

void update_net_by_pin() {
  for(int i = 0; i < netCNT; i++) {
    NET *curNet = &netInstance[i];

    PIN *pin_xmin = NULL, *pin_ymin = NULL;
    PIN *pin_xmax = NULL, *pin_ymax = NULL;

    prec min_x = PREC_MAX, min_y = PREC_MAX;
    prec max_x = PREC_MIN, max_y = PREC_MIN;

    for(int j = 0; j < curNet->pinCNTinObject; j++) {
      PIN *pin = curNet->pin[j];

      if(pin_xmin) {
        if(min_x > pin->fp.x)  // mdp->center.x)
        {
          min_x = pin->fp.x;  // mdp->center.x ;
          pin_xmin->X_MIN = 0;
          pin_xmin = pin;
          pin->X_MIN = 1;
        }
      }
      else {
        min_x = pin->fp.x;  // mdp->center.x ;
        pin_xmin = pin;
        pin->X_MIN = 1;
      }

      if(pin_ymin) {
        if(min_y > pin->fp.y)  // mdp->center.y)
        {
          min_y = pin->fp.y;  // mdp->center.y ;
          pin_ymin->Y_MIN = 0;
          pin_ymin = pin;
          pin->Y_MIN = 1;
        }
      }
      else {
        min_y = pin->fp.y;  // mdp->center.y ;
        pin_ymin = pin;
        pin->Y_MIN = 1;
      }

      if(pin_xmax) {
        if(max_x < pin->fp.x)  // mdp->center.x)
        {
          max_x = pin->fp.x;  // mdp->center.x ;
          pin_xmax->X_MAX = 0;
          pin_xmax = pin;
          pin->X_MAX = 1;
        }
      }
      else {
        max_x = pin->fp.x;  // mdp->center.x ;
        pin_xmax = pin;
        pin->X_MAX = 1;
      }

      if(pin_ymax) {
        if(max_y < pin->fp.y)  // mdp->center.y)
        {
          max_y = pin->fp.y;  // mdp->center.y ;
          pin_ymax->Y_MAX = 0;
          pin_ymax = pin;
          pin->Y_MAX = 1;
        }
      }
      else {
        max_y = pin->fp.y;  // mdp->center.y ;
        pin_ymax = pin;
        pin->Y_MAX = 1;
      }
    }
    curNet->min_x = min_x;
    curNet->min_y = min_y;
    curNet->max_x = max_x;
    curNet->max_y = max_y;
  }
}

//
// CreateSparseMatrix Routine
//
// using current Pin's structure,
// based on the B2B models,
// it genereates Sparsematrix into Eigen formats.
//
void CreateSparseMatrix(VectorXf &xcg_x, VectorXf &xcg_b, VectorXf &ycg_x,
                        VectorXf &ycg_b, SMatrix &eMatX, SMatrix &eMatY) {
  int pinCNTinObject = 0;
  int moduleID1 = 0;
  int moduleID2 = 0;
  int is_term1 = 0;
  int is_term2 = 0;

  prec common1;
  prec common2;

  NET *tempNet = NULL;
  PIN *pin1 = NULL;
  PIN *pin2 = NULL;

  MODULE *mdp1 = NULL;
  MODULE *mdp2 = NULL;
  TERM *term1 = NULL;
  TERM *term2 = NULL;

  FPOS center1, center2;
  FPOS fp1, fp2;

  // to easily convert (i, j, value) -> CSR sparse Matrix.
  // using Eigen library..
  vector< T > tripletListX, tripletListY;
  tripletListX.reserve(10000000);
  tripletListY.reserve(10000000);

  // xcg_x & ycg_x update
  MODULE *curModule = NULL;
  for(int i = 0; i < moduleCNT; i++) {
    curModule = &moduleInstance[i];

    // 1d prec array
    xcg_x(i) = curModule->center.x;
    ycg_x(i) = curModule->center.y;

    xcg_b(i) = ycg_b(i) = 0;
  }

  for(int i = 0; i < netCNT; i++) {
    tempNet = &netInstance[i];
    pinCNTinObject = tempNet->pinCNTinObject;
    common1 = 1.0 / ((prec)pinCNTinObject - 1.0);

    for(int j = 0; j < pinCNTinObject; j++) {
      pin1 = tempNet->pin[j];
      moduleID1 = pin1->moduleID;
      fp1 = pin1->fp;

      // if is not Terminal -> moduleInstance
      if(!pin1->term) {
        mdp1 = &moduleInstance[moduleID1];
        center1 = mdp1->center;
        is_term1 = 0;
      }
      // is terminal -> terminalInstance
      else {
        term1 = &terminalInstance[moduleID1];
        center1 = term1->center;
        is_term1 = 1;
      }

      for(int k = j + 1; k < pinCNTinObject; k++) {
        pin2 = tempNet->pin[k];
        moduleID2 = pin2->moduleID;
        fp2 = pin2->fp;

        if(!pin2->term) {
          mdp2 = &moduleInstance[moduleID2];
          center2 = mdp2->center;
          is_term2 = 0;
        }
        else {
          term2 = &terminalInstance[moduleID2];
          center2 = term2->center;
          is_term2 = 1;
        }

        // there is no need to calculate (for same nodes)
        if(moduleID1 == moduleID2 && is_term1 == is_term2) {
          continue;
        }

        if(pin1->X_MIN || pin1->X_MAX || pin2->X_MIN || pin2->X_MAX) {
          prec len_x = fabs(fp1.x - fp2.x);

          prec wt_x = 0.0f;
          if(dge(len_x, MIN_LEN)) {
            wt_x = common1 / len_x;
          }
          else {
            wt_x = common1 / MIN_LEN;
          }

          common2 = (-1.0) * wt_x;

          if(wt_x < 0)
            printf("ERROR WEIGHT\n");

          // both is module
          if(!is_term1 && !is_term2) {
            tripletListX.push_back(T(moduleID1, moduleID1, wt_x));
            tripletListX.push_back(T(moduleID2, moduleID2, wt_x));
            tripletListX.push_back(T(moduleID1, moduleID2, common2));
            tripletListX.push_back(T(moduleID2, moduleID1, common2));

            xcg_b(moduleID1) +=
                common2 * ((fp1.x - center1.x) - (fp2.x - center2.x));
            xcg_b(moduleID2) +=
                common2 * ((fp2.x - center2.x) - (fp1.x - center1.x));
          }
          // 1 is terminal, 2 is module
          else if(is_term1 && !is_term2) {
            tripletListX.push_back(T(moduleID2, moduleID2, wt_x));
            xcg_b(moduleID2) += wt_x * (fp1.x - (fp2.x - center2.x));
          }
          // 2 is terminal, 1 is module
          else if(!is_term1 && is_term2) {
            tripletListX.push_back(T(moduleID1, moduleID1, wt_x));
            xcg_b(moduleID1) += wt_x * (fp2.x - (fp1.x - center1.x));
          }
        }

        if(pin1->Y_MIN || pin1->Y_MAX || pin2->Y_MIN || pin2->Y_MAX) {
          prec len_y = fabs(fp1.y - fp2.y);

          prec wt_y = 0.0f;
          if(dge(len_y, MIN_LEN)) {
            wt_y = common1 / len_y;
          }
          else {
            wt_y = common1 / MIN_LEN;
          }
          common2 = (-1.0) * wt_y;

          // both is module
          if(!is_term1 && !is_term2) {
            tripletListY.push_back(T(moduleID1, moduleID1, wt_y));
            tripletListY.push_back(T(moduleID2, moduleID2, wt_y));
            tripletListY.push_back(T(moduleID1, moduleID2, common2));
            tripletListY.push_back(T(moduleID2, moduleID1, common2));

            ycg_b(moduleID1) +=
                common2 * ((fp1.y - center1.y) - (fp2.y - center2.y));
            ycg_b(moduleID2) +=
                common2 * ((fp2.y - center2.y) - (fp1.y - center1.y));
          }
          // 1 is terminal, 2 is module
          else if(is_term1 && !is_term2) {
            tripletListY.push_back(T(moduleID2, moduleID2, wt_y));
            ycg_b(moduleID2) += wt_y * (fp1.y - (fp2.y - center2.y));
          }
          // 2 is terminal, 1 is module
          else if(!is_term1 && is_term2) {
            tripletListY.push_back(T(moduleID1, moduleID1, wt_y));
            ycg_b(moduleID1) += wt_y * (fp2.y - (fp1.y - center1.y));
          }
        }
      }
    }
  }

  eMatX.setFromTriplets(tripletListX.begin(), tripletListX.end());
  eMatY.setFromTriplets(tripletListY.begin(), tripletListY.end());
}


///////////////////////////////////////////////////////////////////////////////////////////
void generateSearchOrder(std::vector<int> &order,int scale){
  //buffer first
  int containerCNT = scale*(scale+1);
  int scaleSquare = scale*scale;
  order.resize(containerCNT);
  /* plan 1  buffer first
  for(int i = scaleSquare;i < containerCNT;i++)
  {
    order[i-scaleSquare] = i;
  }
  //then PEs 
  for(int i = 0;i < scaleSquare;i++)
  {
    order[i+scale] = i;
  }
  */

  //蛇形方案
  int labeledNum = 0;
  // while(labeledNum<scale*scale){
  //   int bestabsolutDis = scale*scale;
  //   int bestid;
  //   int halfscale = scale/2;
  //   for(int i = 0;i < scale*scale;i++)
  //   {
  //     POS cord;
  //     cord.x = i/scale;
  //     cord.y = i%scale;
  //     int tempdistance = abs(cord.x-halfscale)+abs(cord.y-halfscale);
  //     if(tempdistance<bestabsolutDis) 
  //     {
  //       if(labeledNum == 1|| find(order.begin(),order.begin()+labeledNum,i)==order.begin()+labeledNum)
  //       {
  //         bestabsolutDis = tempdistance;
  //         bestid = i;
  //       }
  //     }
  //   }
  //   order[labeledNum] = bestid;
  //   labeledNum++;
  // }
  POS lastCord,curCord;
  lastCord.x = 7;
  lastCord.y = 7;
  order[0] = lastCord.x*scale+lastCord.y;
  int curOrderId = 1;
  for(int step = 1;step < 17;step++)
  {
    int finish = step == 16;
    
    if(step%2 == 1)
    {
      //向右
      for(int j = 0;j < step;j++)
      {
        curCord.x = lastCord.x +1;
        curCord.y = lastCord.y;
        order[curOrderId++] = curCord.x*scale+curCord.y;
        lastCord.x = curCord.x;
        lastCord.y = curCord.y;
      }
      
      //向上
      for(int j = 0;j < step;j++)
      {
        curCord.x = lastCord.x;
        curCord.y = lastCord.y+1;
        order[curOrderId++] = curCord.x*scale+curCord.y;
        lastCord.x = curCord.x;
        lastCord.y = curCord.y;
      }
    }
    else{
      if(finish)
      {
      step = 15;
      }
      //向左
      for(int j = 0;j < step;j++)
      {
        curCord.x = lastCord.x -1;
        curCord.y = lastCord.y;
        order[curOrderId++] = curCord.x*scale+curCord.y;
        lastCord.x = curCord.x;
        lastCord.y = curCord.y;
      }
      
      //向下
      if(finish)
      {
        break;
      }
      for(int j = 0;j < step;j++)
      {
        curCord.x = lastCord.x;
        curCord.y = lastCord.y-1;
        order[curOrderId++] = curCord.x*scale+curCord.y;
        lastCord.x = curCord.x;
        lastCord.y = curCord.y;
      }
    }
    if(finish)
    {
      break;;
    }
  }
  for(int i = scale*scale;i < containerCNT;i++)
  {
    order[i] = i;
  }
  
  printDebugInfo(__FUNCTION__,0,"");
  for(int i = 0;i < order.size();i++)
  {
    cout<<order[i]<<" ";
  }
  //check miss
  for(int i = 0;i < scale*scale+scale;i++){
    auto it  = find(order.begin(),order.end(),i);
    if(it == order.end())
    {
      cout<<"mising "<<i<<endl;
    }
  }
  cout<<endl;
  // exit(0);
}
void doDeterPlace(int scale){
  std::vector<int> _Grid;
    std::vector<MODULE> _Containers;
    std::vector<PIN> _TermPins;
    std::vector<NET> _Nets;
    int debugIdx = 0;
    printDebugInfo(__FUNCTION__, debugIdx++,"");
    FPOS _GridSize;
    POS _GridNum;
    
    setupGreedyPlace(_GridSize,_GridNum,_Grid);

    vector<int> searchOrder;
    generateSearchOrder(searchOrder,scale);
    for(int i = 0;i < searchOrder.size();i++)
    {
        int moduleID = searchOrder[i];
        vector<FPOS> candidates;
        int candidateIdx;
        candidates.clear();
        findCandidate(moduleInstance[moduleID],searchOrder,candidates,scale);
        if(i<searchOrder.size()-1){
          int nextModuleID = searchOrder[i+1];
          candidateIdx = lookaheadDetermine(moduleInstance[moduleID], moduleInstance[nextModuleID], searchOrder, candidates,scale);
          // candidateIdx = 0;
        }
        else{
          candidateIdx = 0;
        }
        placeContainer(moduleInstance[moduleID], candidates[candidateIdx]);
    }
    cout<<"place end = "<<place.end.x<<" "<<place.end.y<<endl;
    plot(scale);
}

void findCandidate(MODULE& curModule,vector<int> &order,vector<FPOS> &candidates,int scale){
  vector<int> activeNets;
  int debugNum = 0;
  printDebugInfo(__FUNCTION__, debugNum++, "");
  findActiveNet(curModule, order, activeNets, scale);
  vector<prec> pinXlist;
  vector<prec> pinYlist;
  printDebugInfo(__FUNCTION__, debugNum++, "");
  fillPinCordList(pinXlist, pinYlist, activeNets, curModule, order);
  printDebugInfo(__FUNCTION__, debugNum++, "pinXlist size "+to_string(pinXlist.size()));
  sort(pinXlist.begin(),pinXlist.end());
  sort(pinYlist.begin(),pinYlist.end());
  FPOS fp;
  int halfPinSize = pinYlist.size()/2;
  if(pinYlist.size() == 0){
    fp.x = curModule.center.x;
    fp.y = curModule.center.y;
  }
  else{
    fp.x = pinXlist[halfPinSize];
    fp.y = pinYlist[halfPinSize];
  }
  
  FPOS direction;
  // direction.x = 1.0;
  // direction.y = 0.0;
  fp.x = fp.x > place.org.x + curModule.half_size.x?fp.x:place.org.x + curModule.half_size.x;
  fp.x = fp.x < place.end.x - curModule.half_size.x?fp.x:place.end.x - curModule.half_size.x;
  fp.y = fp.y > place.org.y + curModule.half_size.y?fp.y:place.org.y + curModule.half_size.y;
  fp.y = fp.y < place.end.y - curModule.half_size.y?fp.y:place.end.y - curModule.half_size.y;
  
  move2nonOverlap(curModule,fp,direction,order,candidates);
  int idx = halfPinSize;
  while(candidates.size() == 0){
      // int idX,idY;
      // idX = idx>=0?idx:0;
      // idY = idx>=0?idx:0;
      // idx--;
      fp.x += curModule.size.x;
      fp.y += curModule.size.y;
      move2nonOverlap(curModule,fp,direction,order,candidates);
    
  }
  printDebugInfo(__FUNCTION__, debugNum++, "candidate size:="+ to_string(candidates.size()));
  for(int i = 0;i < candidates.size();i++)
  {
    cout<<"module "<<curModule.idx<<" cadidate "<<i<<" : ("<<candidates[i].x<<", "<<candidates[i].y<<endl;
  }

  // exit(0);
  // direction.x = -1.0;
  // direction.y = 0.0;

  // move2nonOverlapFp(curModule,fp,direction,order,candidates);
  // direction.x = 0.0;
  // direction.y = 1.0;

  // move2nonOverlapFp(curModule,fp,direction,order,candidates);
  // direction.x = 0.0;
  // direction.y = -1.0;

  // move2nonOverlapFp(curModule,fp,direction,order,candidates);
}
void fillPinCordList(vector<prec> pinXlist,vector<prec> pinYlist,vector<int> &Nets,MODULE& curModule,vector<int>& order){
  int ModuleIdx = curModule.idx;
  auto it = find(order.begin(),order.end(),ModuleIdx);
  for(int i = 0;i < Nets.size();i++)
  {
    NET* net = &netInstance[i];
    for(int j = 0;j < net->pinCNTinObject;j++)
    {
      PIN* pin = net->pin[j];
      int moduleIdOfPin = pin->moduleID;
      auto temp_it = find(order.begin(),it,moduleIdOfPin);
      if(temp_it!= it ||pin->term){
        pinXlist.push_back(pin->fp.x);
        pinYlist.push_back(pin->fp.y);
      }
    }
  }
}
void move2nonOverlap(MODULE& curModule,FPOS fp,FPOS direction,vector<int> &order,vector<FPOS> &candidates){
  auto it = find(order.begin(),order.end(),curModule.idx);
  
  if(it == order.begin()){
    candidates.push_back(fp);
  }
  // else{
  //   prec x_max,x_min,y_max,y_min;
  //   int firstX = 1;
  //   int firstY = 1;
  //   for(int i = 0;i < (it - order.begin());i++)
  //   {
  //     MODULE* tmpModule = &moduleInstance[order[i]];
  //     if(tmpModule->center.x-tmpModule->half_size.x < fp.x+curModule.half_size.x && tmpModule->center.x+tmpModule->half_size.x > fp.x-curModule.half_size.x)
  //     {
  //       if(firstY)
  //       {
  //         firstY = 0;
  //         y_min = tmpModule->center.y-tmpModule->half_size.y;
  //         y_max = tmpModule->center.y+tmpModule->half_size.y;
  //       }
  //       else{
  //         y_min = (tmpModule->center.y-tmpModule->half_size.y) < y_min?(tmpModule->center.y-tmpModule->half_size.y) : y_min;
  //         y_max = (tmpModule->center.y+tmpModule->half_size.y) > y_max?(tmpModule->center.y+tmpModule->half_size.y) : y_max;
  //       }
  //     }

  //     if(tmpModule->center.y-tmpModule->half_size.y < fp.y+curModule.half_size.y && tmpModule->center.y+tmpModule->half_size.y > fp.x-curModule.half_size.y)
  //     {
  //       if(firstX)
  //       {
  //         firstX = 0;
  //         x_min = tmpModule->center.x-tmpModule->half_size.x;
  //         x_max = tmpModule->center.x+tmpModule->half_size.x;
  //       }
  //       else{
  //         x_min = (tmpModule->center.x-tmpModule->half_size.x) < x_min?(tmpModule->center.x-tmpModule->half_size.x) : x_min;
  //         x_max = (tmpModule->center.x+tmpModule->half_size.x) > x_max?(tmpModule->center.x+tmpModule->half_size.x) : x_max;
  //       }
  //     }


  //   }
  //   FPOS candidate1, candidate2,candidate3,candidate4;
  //   candidate1.x = x_min - curModule.half_size.x;
  //   candidate1.y = fp.y;

  //   candidate2.x = x_max + curModule.half_size.x;
  //   candidate2.y = fp.y;

  //   candidate3.y = y_min - curModule.half_size.y;
  //   candidate3.x = fp.x;

  //   candidate4.y = y_max + curModule.half_size.y;
  //   candidate4.x = fp.x;
  //   printDebugInfo(__FUNCTION__, 0, "X MIN "+ to_string(x_min)+" X MAX "+to_string(x_max));
  //   if(candidate1.x>place.org.x &&candidate1.x<place.end.x && candidate1.y>place.org.y &&candidate1.y<place.end.y&&firstX==0)
  //   {
  //     candidates.push_back(candidate1);
  //   }
  //   if(candidate2.x>place.org.x &&candidate2.x<place.end.x && candidate2.y>place.org.y &&candidate2.y<place.end.y&&firstX==0)
  //   {
  //     candidates.push_back(candidate2);
  //   }
  //   if(candidate3.x>place.org.x &&candidate3.x<place.end.x &&candidate3.y>place.org.y &&candidate3.y<place.end.y&&firstY==0)
  //   {
  //     candidates.push_back(candidate3);
  //   }
  //   if(candidate4.x>place.org.x &&candidate4.x<place.end.x &&candidate4.y>place.org.y &&candidate4.y<place.end.y&&firstY==0)
  //   {
  //     candidates.push_back(candidate4);
  //   }
  // }

  else{
    vector<prec> Xcontainer,Ycontainer;
    for(int i = 0;i < (it - order.begin());i++)
    {
      MODULE* tmpModule = &moduleInstance[order[i]];
      if(tmpModule->center.x-tmpModule->half_size.x < fp.x+curModule.half_size.x && tmpModule->center.x+tmpModule->half_size.x > fp.x-curModule.half_size.x)
      {
        // if(tmpModule->center.y-tmpModule->half_size.y > fp.y+curModule.half_size.y || tmpModule->center.y+tmpModule->half_size.y < fp.x-curModule.half_size.y)
        // {
        //   Ycontainer.push_back(tmpModule->center.y);
        // }
        Ycontainer.push_back(tmpModule->center.y);
        
      }
      if(tmpModule->center.y-tmpModule->half_size.y < fp.y+curModule.half_size.y && tmpModule->center.y+tmpModule->half_size.y > fp.x-curModule.half_size.y)
      {
        // if(tmpModule->center.x-tmpModule->half_size.x > fp.x+curModule.half_size.x || tmpModule->center.x+tmpModule->half_size.x < fp.x-curModule.half_size.x)
        // {
        //    Xcontainer.push_back(tmpModule->center.x);
        // }
        Xcontainer.push_back(tmpModule->center.x);
       
      }

    }
    Xcontainer.push_back(place.org.x);
    Xcontainer.push_back(place.end.x);
    Ycontainer.push_back(place.org.y);
    Ycontainer.push_back(place.end.y);
    sort(Xcontainer.begin(),Xcontainer.end());
    sort(Ycontainer.begin(),Ycontainer.end());
    cout<<"XCONTAINER SIZE = "<<Xcontainer.size()<<endl;
    cout<<"YCONTAINER SIZE = "<<Ycontainer.size()<<endl;
    FPOS closestLocation;
    prec distance;
    distance = DBL_MAX;
    int found = 0;
    int debugIDX = 0;
    for(int i = 0;i < Xcontainer.size()-1;i++)
    {
      // printDebugInfo(__FUNCTION__,debugIDX++ , "curModule size = "+to_string(curModule.size.x));
      if((Xcontainer[i+1] - Xcontainer[i])>2* curModule.size.x){
        prec tmpDis = abs(Xcontainer[i+1]-fp.x)<abs(Xcontainer[i]-fp.x)?abs(Xcontainer[i+1]-fp.x):abs(Xcontainer[i]-fp.x);
        if(tmpDis<distance){
          distance = tmpDis;
          if(abs(Xcontainer[i+1]-fp.x)<abs(Xcontainer[i]-fp.x))
          {
            closestLocation.x = Xcontainer[i+1] - curModule.size.x;
            closestLocation.y = fp.y;
            found= 1;
          }
          else{
            closestLocation.x = Xcontainer[i] + curModule.size.x;
            closestLocation.y = fp.y;
            found = 1;
          }
        }
      }
    }
    if(found){
      if(Xcontainer.size() == 2){
        FPOS tmp;
        tmp.x = fp.x;
        tmp.y = fp.y;
        candidates.push_back(tmp);
      }
      else{
        candidates.push_back(closestLocation);
      }
      
    }
     if(found == 0){
      cout<<"not found in X"<<endl;
    }
    found = 0;
    FPOS closestLocation2;
    distance = DBL_MAX;
    for(int i = 0;i < Ycontainer.size()-1;i++)
    {
      if((Ycontainer[i+1] - Ycontainer[i])> 2*curModule.size.y){
        prec tmpDis = abs(Ycontainer[i+1]-fp.y)<abs(Ycontainer[i]-fp.y)?abs(Ycontainer[i+1]-fp.y):abs(Ycontainer[i]-fp.y);
        if(tmpDis<distance){
          distance = tmpDis;
          if(abs(Ycontainer[i+1]-fp.y)<abs(Ycontainer[i]-fp.y))
          {
            closestLocation2.y = Ycontainer[i+1] - curModule.size.y;
            closestLocation2.x = fp.x;
            found = 1;
          }
          else{
            closestLocation2.y = Ycontainer[i] + curModule.size.y;
            closestLocation2.x = fp.x;
            found = 1;
          }
        }
      }
    }
    if(found){
      if(Ycontainer.size() == 2){
        FPOS tmp;
        tmp.x = fp.x;
        tmp.y = fp.y;
        candidates.push_back(tmp);
      }
      else{
        candidates.push_back(closestLocation2);
      }
      
    }
    if(found == 0){
      cout<<"not found in Y"<<endl;
      // for(int i = 0;i < Ycontainer.size();i++)
      // {
      //   cout<<Ycontainer[i]<<endl;
      // }
    }
  }
  
}
int lookaheadDetermine(MODULE& curModule,MODULE& nextModule,vector<int> order,vector<FPOS>& candidates,int scale){
  vector<int> nets1;
  findActiveNet(curModule, order, nets1, scale);
  prec minWL = DBL_MAX;
  int bestCandidateId;
  for(int i = 0;i <  candidates.size();i++){
    vector<FPOS> tmpCandidates;
    curModule.center.x = candidates[i].x;
    curModule.center.y = candidates[i].y;
    findCandidate(nextModule,order,tmpCandidates,scale);
    vector<int> nets2;
    nets2.clear();
    findActiveNet(nextModule, order, nets2, scale);
    for(int j = 0;j < nets1.size();j++)
    {
      nets2.push_back(nets1[j]);
    }
    
    for(int j = 0;j < tmpCandidates.size();j++)
    {
      prec hpwl = getHpwl(nextModule,nets2,tmpCandidates[j]);
      if(hpwl<minWL)
      {
        minWL = hpwl;
        bestCandidateId = i;
      }
    }
    

  }

  return bestCandidateId;
}

void placeContainer(MODULE &curContainer, FPOS location){
  int debugIdx = 0;
  // printDebugInfo(__FUNCTION__,debugIdx++,"binId = "+to_string(binId)+"total bin #= "+to_string(Grid.size()));
  
  // printDebugInfo(__FUNCTION__,debugIdx++,"location = "+to_string(location.x)+" "+to_string(location.y));
  
  curContainer.center.x = location.x;
  // printDebugInfo(__FUNCTION__,debugIdx++,"");
  curContainer.center.y = location.y;
  // printDebugInfo(__FUNCTION__,debugIdx++,"");
  for(int i = 0;i < curContainer.pinCNTinObject;i++){
    curContainer.pin[i]->fp.x = location.x + curContainer.pof[i].x;
    curContainer.pin[i]->fp.y = location.y + curContainer.pof[i].y;
  }
}
void doGreedyPlace(int scale){
    std::vector<int> _Grid;
    std::vector<MODULE> _Containers;
    std::vector<PIN> _TermPins;
    std::vector<NET> _Nets;
    int debugIdx = 0;
    printDebugInfo(__FUNCTION__, debugIdx++,"");
    FPOS _GridSize;
    POS _GridNum;
    
    setupGreedyPlace(_GridSize,_GridNum,_Grid);

    vector<int> searchOrder;
    generateSearchOrder(searchOrder,scale);
    for(int i = 0;i < searchOrder.size();i++)
    {
        int moduleID = searchOrder[i];
        printDebugInfo(__FUNCTION__, debugIdx++,"moduleId = "+to_string(moduleID));
        int placedBinId = findBestBin(moduleInstance[moduleID],searchOrder,scale,_Grid,_GridSize);
        FPOS location;
        location.x  = (prec)(placedBinId/_GridNum.y) * _GridSize.x;
        location.y = (prec)(placedBinId%_GridNum.y)* _GridSize.y;
        printDebugInfo(__FUNCTION__, debugIdx++,"placed to  ");
        placeContainerToBin(moduleInstance[moduleID],placedBinId,location,_Grid);
    }
    plot(scale);
}

void setupGreedyPlace(FPOS& gridSize,POS& gridNum,std::vector<int>& Grid){
  setupGrid(gridSize,gridNum,Grid);
}

void setupGrid(FPOS& gridSize,POS& gridNum,std::vector<int>& Grid){
  
  // gridSize.x = moduleInstance[0].size.x*1.2;
  // gridSize.y = moduleInstance[0].size.y*1.2;
    gridSize.x = terminalInstance[0].size.x*1.3;
  gridSize.y = terminalInstance[0].size.y*1.3;
  gridNum.x = (int) (place.end.x - place.org.x)/gridSize.x;
  gridNum.y = (int) (place.end.y - place.org.y)/gridSize.y;
  Grid.resize(gridNum.x*gridNum.y);
  for(int i = 0;i<Grid.size();i++)
  {
    Grid[i] = 0;
  }
  printDebugInfo(__FUNCTION__, 0,"GridNum := "+to_string(gridNum.x)+" "+to_string(gridNum.y));
}

int findBestBin(MODULE& curContainer,vector<int>& order,int scale,std::vector<int> Grid,FPOS gridSize){
  vector<int> activeNetId;
  FPOS LL,UR;
  int debugIdx = 0;
  // printDebugInfo(__FUNCTION__,debugIdx++,"");
  findActiveNet(curContainer, order,activeNetId,scale);
  //  printDebugInfo(__FUNCTION__,debugIdx++,"");
  findBoundingBox(activeNetId,LL,UR,curContainer,order);
  //  printDebugInfo(__FUNCTION__,debugIdx++,"");
  vector<FPOS> candidates;
  vector<int> candidateBinId;
  findCandidate(Grid,LL,UR,candidates,candidateBinId,gridSize);
  //  printDebugInfo(__FUNCTION__,debugIdx++,"");
  prec minWL = DBL_MAX;
  FPOS bestLocation;
  int bestBinId;
  for(int i = 0;i <  candidates.size();i++)
  {
    //  printDebugInfo(__FUNCTION__,debugIdx++,"");
    prec tmp = getHpwl(curContainer,activeNetId,candidates[i]);
    // printDebugInfo(__FUNCTION__,debugIdx++,"getHpwl = "+to_string(tmp));
    if(tmp < minWL)
    {
      minWL = tmp;
      // printDebugInfo(__FUNCTION__,debugIdx++,"");
      bestLocation.x = candidates[i].x;
      bestLocation.y = candidates[i].y;
      // printDebugInfo(__FUNCTION__,debugIdx++,"candidateBin SIZE"+to_string(candidateBinId.size()));
      bestBinId = candidateBinId[i];
      // printDebugInfo(__FUNCTION__,debugIdx++,"");
    }
    // printDebugInfo(__FUNCTION__,debugIdx++,"update minWL  "+to_string(minWL));
  }
  
  // printDebugInfo(__FUNCTION__,0,"bestBinID "+to_string(bestBinId)+" ("+to_string(bestLocation.x)+", "+to_string(bestLocation.y)+")");

  // exit(0);
  return bestBinId;
  
}
void findActiveNet(MODULE& curContainer,vector<int>& order,vector<int>& activeNetId,int scale)
{
  int debugidx = 0;
  // printDebugInfo(__FUNCTION__,debugidx++,"");
  int curModuleIdx = curContainer.containerID(scale);
  auto it = find(order.begin(),order.end(),curModuleIdx);
  int searchIdx;
  
  if(it!= order.end())
  {
    searchIdx = it - order.begin();
  }
  // printDebugInfo(__FUNCTION__,debugidx++,"");
  for(int i = 0;i < curContainer.pinCNTinObject;i++)
  {
    int netid = curContainer.pin[i]->netID;
    NET* curNet = &netInstance[netid];
    int active = 1;
    // printDebugInfo(__FUNCTION__,debugidx++,"");
    for(int j = 0;j < curNet->pinCNTinObject;j++)
    {
      PIN* curPin = curNet->pin[j];
      if(curPin->term)
      {
        continue;
      }
      auto it2 = find(order.begin(),it,curPin->moduleID);
      if(it2==order.end())
      {
        active = 0;
        break;
      }
    }
    // printDebugInfo(__FUNCTION__,debugidx++,"");
    if(active)
    {
      activeNetId.push_back(netid);
    }
  }
  // printDebugInfo(__FUNCTION__,debugidx++,"");
  // for(int i = 0;i < activeNetId.size();i++)
  // {
  //   cout<<activeNetId[i]<<" ";
  // }
  // cout<<endl;
}

void doSchemPlace_TermContainer(int scale){
  std::vector<int> _Grid;
    std::vector<MODULE> _Containers;
    std::vector<PIN> _TermPins;
    std::vector<NET> _Nets;
    int debugIdx = 0;
    printDebugInfo(__FUNCTION__, debugIdx++,"");
    FPOS _GridSize;
    POS _GridNum;
    
    setupGreedyPlace(_GridSize,_GridNum,_Grid);

    vector<int> searchOrder;
    generateSearchOrder(searchOrder,scale);
    for(int i = 0;i < searchOrder.size();i++)
    {
        int moduleID = searchOrder[i];
        printDebugInfo(__FUNCTION__, debugIdx++,"moduleId = "+to_string(moduleID));
        int placedBinId = findBestBin(terminalInstance[moduleID],searchOrder,scale,_Grid,_GridSize);
        FPOS location;
        location.x  = (prec)(placedBinId/_GridNum.y+0.9) * _GridSize.x  ;
        location.y = (prec)(placedBinId%_GridNum.y+0.9)* _GridSize.y;
        printDebugInfo(__FUNCTION__, debugIdx++,"placed to  ");
        placeContainerToBin(terminalInstance[moduleID],placedBinId,location,_Grid);
    }
    plotTerm(scale);
}

int findBestBin(TERM& curContainer,vector<int>& order,int scale,std::vector<int> Grid,FPOS gridSize){

  int bestBinId;
  int offset = 0;
  int idx = curContainer.idx;
  POS containerCORD,gridCORD;
  POS gridNum;
  printDebugInfo(__FUNCTION__, 0, "IDX = "+to_string(idx));
  gridNum.x = (int) (place.end.x - place.org.x)/gridSize.x;
  gridNum.y = (int) (place.end.y - place.org.y)/gridSize.y;
  if(idx>=scale*scale){
    //buffer
    containerCORD.x = idx - scale*scale;
    gridCORD.x = 0;
    gridCORD.y = containerCORD.x + offset;
  }
  else{
    //pe
    containerCORD.x = idx/scale;
    containerCORD.y = idx%scale;
    gridCORD.x = containerCORD.x+1;
    gridCORD.y = containerCORD.y+offset;
  }
  bestBinId = gridCORD.x*gridNum.y + gridCORD.y;

  // exit(0);
  return bestBinId;
}
prec getHpwl(MODULE& curContainer,std::vector<int> activeNets,FPOS location){
  curContainer.center.x = location.x;
  curContainer.center.y = location.y;
  for(int i = 0;i < curContainer.pinCNTinObject;i++)
  {
    curContainer.pin[i]->fp.x=location.x + curContainer.pof[i].x;
    curContainer.pin[i]->fp.y=location.y + curContainer.pof[i].y;
  }
  prec hpwl = 0.0;
  for(int i = 0;i < activeNets.size();i++)
  {
    NET* curNet = &netInstance[activeNets[i]];
    prec x_max,x_min,y_max,y_min;
    for(int j = 0;j <  curNet->pinCNTinObject;j++)
    {
      if(j == 0)
      {
        x_max = curNet->pin[j]->fp.x;
        x_min = curNet->pin[j]->fp.x;
        y_max = curNet->pin[j]->fp.y;
        y_min = curNet->pin[j]->fp.y;
      }
      else{
        x_max = x_max>curNet->pin[j]->fp.x?x_max:curNet->pin[j]->fp.x;
        x_min = x_min<curNet->pin[j]->fp.x?x_min:curNet->pin[j]->fp.x;
        y_max = y_max>curNet->pin[j]->fp.y?x_max:curNet->pin[j]->fp.y;
        y_min = y_min<curNet->pin[j]->fp.y?x_min:curNet->pin[j]->fp.y;
      }
    }
    hpwl +=  x_max-x_min+y_max-y_min;
  }
  return hpwl;
}

void findBoundingBox(vector<int> &activeNetId,FPOS& LL,FPOS& UR,MODULE& curContainer,vector<int> Order){
  int curModuleID = curContainer.idx;
  auto it = find(Order.begin(),Order.end(),curModuleID);

  for(int i = 0;i < activeNetId.size();i++)
  {
    NET* curNet = &netInstance[activeNetId[i]];
    prec x_max,x_min,y_max,y_min;
    bool first = 1;
    for(int j = 0;j < curNet->pinCNTinObject;j++)
    {
      int _pinInModuleId = curNet->pin[j]->moduleID;
      if(first == 1)
      {
        
        if(curNet->pin[j]->term||find(Order.begin(),it,_pinInModuleId)!=it)
        {
          first = 0;
          x_max = curNet->pin[j]->fp.x;
          x_min = curNet->pin[j]->fp.x;
          y_max = curNet->pin[j]->fp.y;
          y_min = curNet->pin[j]->fp.y;
        }
      }
      else{
        if(curNet->pin[j]->term||find(Order.begin(),it,_pinInModuleId)!=it)
        {
          
          x_max = x_max>curNet->pin[j]->fp.x?x_max:curNet->pin[j]->fp.x;
          x_min = x_min<curNet->pin[j]->fp.x?x_min:curNet->pin[j]->fp.x;
          y_max = y_max>curNet->pin[j]->fp.y?y_max:curNet->pin[j]->fp.y;
          y_min = y_min<curNet->pin[j]->fp.y?y_min:curNet->pin[j]->fp.y;
        }
      }
    }
    if(i == 1){
      LL.x = x_min;
      LL.y = y_min;
      UR.x = x_max;
      UR.y = y_max;
    }
    else{
      LL.x = LL.x<x_min?LL.x:x_min;
      LL.y = LL.y<y_max?LL.y:y_max;
      UR.x = UR.x<x_min?UR.x:x_min;
      UR.y = UR.y<y_max?UR.y:y_max;
    }
  }
}

void findCandidate(vector<int>&Grid,FPOS LL,FPOS UR,vector<FPOS> &candidates,vector<int>&candidateBinId,FPOS gridSize){
  POS gridNum;
  gridNum.x = (int) (place.end.x - place.org.x)/gridSize.x;
  gridNum.y = (int) (place.end.y - place.org.y)/gridSize.y;
  for(int i = 0;i < gridNum.x;i++)
  {
    int gridIDBase = i * gridNum.y;
    for(int j = 0;j < gridNum.y;j++)
    {
      int gridID = gridIDBase+ j;
      FPOS location;
      location.x = (i+0.5)*gridSize.x;
      location.y = (j+0.5)*gridSize.y;
      if(Grid[gridID] == 0)
      {
        // if(location.x <= UR.x+gridSize.x+2)
        // {
        //   if(location.y <= UR.y+gridSize.y+2)
        //   {
        //     if(location.x >= LL.x-gridSize.x-2)
        //     {
        //       if(location.y >= LL.y-gridSize.y-2){
        //         candidates.push_back(location);
        //         candidateBinId.push_back(gridID);
        //       }
        //     }
        //   }
        // }
        candidates.push_back(location);
        candidateBinId.push_back(gridID);
      }
    }
  }
 
  // for(int i = 0;i < candidateBinId.size();i++)
  // {
  //   cout<<"cadidate "<<i<<" id = "<<candidateBinId[i]<<" ("<<candidates[i].x<<", "<<candidates[i].y<<")"<<endl;
  // }
}


void placeContainerToBin(MODULE& curContainer,int binId,FPOS location,vector<int> &Grid){
  int debugIdx = 0;
  // printDebugInfo(__FUNCTION__,debugIdx++,"binId = "+to_string(binId)+"total bin #= "+to_string(Grid.size()));
  Grid[binId] = 1;
  // printDebugInfo(__FUNCTION__,debugIdx++,"location = "+to_string(location.x)+" "+to_string(location.y));
  
  curContainer.center.x = location.x;
  // printDebugInfo(__FUNCTION__,debugIdx++,"");
  curContainer.center.y = location.y;
  // printDebugInfo(__FUNCTION__,debugIdx++,"");
  for(int i = 0;i < curContainer.pinCNTinObject;i++){
    curContainer.pin[i]->fp.x = location.x + curContainer.pof[i].x;
    curContainer.pin[i]->fp.y = location.y + curContainer.pof[i].y;
  }
  printDebugInfo(__FUNCTION__,0,"placed to" +to_string(location.x)+to_string(location.y));

}

void placeContainerToBin(TERM& curContainer,int binId,FPOS location,vector<int> &Grid){
  int debugIdx = 0;
  printDebugInfo(__FUNCTION__,debugIdx++,"binId = "+to_string(binId)+"total bin #= "+to_string(Grid.size()));
  Grid[binId] = 1;
  printDebugInfo(__FUNCTION__,debugIdx++,"location = "+to_string(location.x)+" "+to_string(location.y));
  
  curContainer.center.x = location.x;
  printDebugInfo(__FUNCTION__,debugIdx++,"");
  curContainer.center.y = location.y;
  printDebugInfo(__FUNCTION__,debugIdx++,"");
  for(int i = 0;i < curContainer.pinCNTinObject;i++){
    curContainer.pin[i]->fp.x = location.x + curContainer.pof[i].x;
    curContainer.pin[i]->fp.y = location.y + curContainer.pof[i].y;
  }
  printDebugInfo(__FUNCTION__,0,"placed to" +to_string(location.x)+to_string(location.y));

}

void plot(int scale){
  int containerNum = scale*(scale+1);
  int PENum = scale*scale;
  vector<int> color1,color2;
  // color1.push_back(255);
  // color1.push_back(255);
  // color1.push_back(255);
  // color2.push_back(100);
  // color2.push_back(100);
  // color2.push_back(100);
  int color;
  string filename = "greedyPlot.gnu";
  std::ofstream outfile;

  outfile.open(filename, std::ios_base::trunc);
  outfile<<"set terminal png size 5500,5500 crop font \"Helvetica,30\""<<endl;
  outfile<<"set output \'output.png\'"<<endl;
  outfile<<"set xrange [-500:5500]"<<endl;
  outfile<<"set yrange [-500:5500]"<<endl;
  outfile<<"set style rect back fs empty border lc rgb \'#008800\'"<<endl;
  outfile.close();
  for(int i =0;i < containerNum;i++){
    FPOS LL,UR;
    LL.x = moduleInstance[i].center.x - moduleInstance[i].half_size.x;
    LL.y = moduleInstance[i].center.y - moduleInstance[i].half_size.y;
    UR.x = moduleInstance[i].center.x + moduleInstance[i].half_size.x;
    UR.y = moduleInstance[i].center.y + moduleInstance[i].half_size.y;
    string label = moduleInstance[i].type+to_string(i/scale)+","+to_string(i%scale);
    if(moduleInstance[i].type == "buffer")
    {
      label = moduleInstance[i].type+to_string(i-scale*scale);
    }
    if(i< PENum){
      color = 1;
      drawRectangle(LL,  UR, color,filename,i,label);
    }
    else{
      color = 2;
      drawRectangle(LL,  UR, color,filename,i,label);
    }
  }
  outfile.open(filename, std::ios_base::app);
  outfile<<"plot x"<<endl;
  outfile.close();
}

void plotTerm(int scale){
  // int containerNum = scale*(scale+1);
  int containerNum = scale;
  int PENum = scale*scale;
  vector<int> color1,color2;
  // color1.push_back(255);
  // color1.push_back(255);
  // color1.push_back(255);
  // color2.push_back(100);
  // color2.push_back(100);
  // color2.push_back(100);
  int color;
  string filename = "greedyPlot.gnu";
  std::ofstream outfile;

  outfile.open(filename, std::ios_base::trunc);
  outfile<<"set terminal png size 5500,5500 crop font \"Helvetica,30\""<<endl;
  outfile<<"set output \'output.png\'"<<endl;
  outfile<<"set xrange [-500:10000]"<<endl;
  outfile<<"set yrange [-500:10000]"<<endl;
  outfile<<"set style rect back fs empty border lc rgb \'#008800\'"<<endl;
  outfile.close();
  for(int i =0;i < containerNum;i++){
    FPOS LL,UR;
    LL.x = terminalInstance[i].center.x - terminalInstance[i].size.x*0.5;
    LL.y = terminalInstance[i].center.y - terminalInstance[i].size.y*0.5;
    UR.x = terminalInstance[i].center.x + terminalInstance[i].size.x*0.5;
    UR.y = terminalInstance[i].center.y + terminalInstance[i].size.y*0.5;
    // string label = "PE"+to_string(i/scale)+","+to_string(i%scale);
    string label = "buffer";
    // if(terminalInstance[i].idx >= scale*scale)
    // {
    //   label = "buffer"+to_string(i-scale*scale);
    // }
    if(i< PENum){
      color = 1;
      drawRectangle(LL,  UR, color,filename,i,label);
    }
    else{
      color = 2;
      drawRectangle(LL,  UR, color,filename,i,label);
    }
  }
  outfile.open(filename, std::ios_base::app);
  outfile<<"plot x"<<endl;
  outfile.close();
}
void drawRectangle(FPOS LL,FPOS UR,int colors,string filename,int objectID,string label){
  std::ofstream outfile;
  outfile.open(filename, std::ios_base::app); // append instead of overwrite
  outfile<<"LABEL"<<objectID+1<< "= \""<< label<<"\""<<endl;
  if(colors ==1)
  {
    
    outfile<<"set object "<<objectID+1<<" rect from "<<(int)LL.x<<","<<(int)LL.y<<" to "<< (int)UR.x<<","<<(int)UR.y<<" lw 5 fs empty border lc rgb \'#880088\'"<<endl;
  }
  else{
    outfile<<"set object "<<objectID+1<<" rect from "<<(int)LL.x<<","<<(int)LL.y<<" to "<< (int)UR.x<<","<<(int)UR.y<<" lw 5 fs empty border lc rgb \'#008800\'"<<endl;
  }
  outfile<<"set label "<<objectID+1<<" at "<<(int)LL.x<<","<<(int)LL.y<<" LABEL"<<objectID+1<<" front center"<<endl;
  
  outfile.close();
  
}
void printDebugInfo(string funcName,int debugIndex,string info){
  
  cout<<"[FUNC]"<<funcName<<" debug "<<debugIndex<<" "<<info<<endl;
    
}


void doBoundingBoxPlaceForBuffer(int scale){
    printDebugInfo(__FUNCTION__, 0, "start");
    FPOS ioCenter;
    findIOCenter(ioCenter);
    
    FPOS peArraySize;
    findPEArraySize(peArraySize);

    
    int bufferNum = scale;

    FPOS bufferSize;
    findBufferSize(bufferSize);
    printDebugInfo(__FUNCTION__, 0, "before exit");
    // exit(0);
    prec PEarrayCenterMax = 2*bufferSize.x;
    std::vector<PIN> _TermPins;
    std::vector<NET> _Nets;
    int debugIdx = 0;
    printDebugInfo(__FUNCTION__, debugIdx++,"");
    prec bestWL = 1e15;
    FPOS bestPEArrayCenter;
    std::vector<FPOS> bestBufferCenters;
    bestPEArrayCenter.y = 2*bufferSize.y;
    for(int i = 0;i  < PEarrayCenterMax;i++){
      prec WL = 0;
      int arrangedNum = 0;
      prec width = (prec) i*bufferSize.x;
      int dis = width;
      vector<FPOS> BufferCenters;
      // BufferCenters.resize(bufferNum);
      int j = 0;
      while (arrangedNum<bufferNum+8)
      {
        prec height =  j*bufferSize.y;
        
        if(j >0){
          arrangedNum += i*2+2*j -1;
          WL += dis/bufferSize.x*(width*2+2*height - 1);
          for(int k = 0;k <= i;k++){
            FPOS center1, center2;
            center1.x = k*bufferSize.x;
            center1.y = height;
            center2.x = k*bufferSize.x;
            center2.y = -height;
            BufferCenters.push_back(center1);
            BufferCenters.push_back(center2);
          }
         
          FPOS center;
          center.x = width+bufferSize.x*j;
          center.y = 0.0;
          BufferCenters.push_back(center);
          for(int k = 1;k < j;k++){
            FPOS center1, center2;
            center1.x = width+bufferSize.x*k;
            center1.y = height - k*bufferSize.y;
            center2.x = width+bufferSize.x*k;
            center2.y = -height + k*bufferSize.y;
            BufferCenters.push_back(center1);
            BufferCenters.push_back(center2);
          }
          
        }

        if(j == 0){
          arrangedNum += i;
          WL += dis/bufferSize.x*(width-1);
          for(int k = 0;k < i;k++){
            FPOS center1;
            center1.x = k*bufferSize.x;
            center1.y = height;
            BufferCenters.push_back(center1);
            
          }
          // printDebugInfo(__FUNCTION__, 11, "Height < bufferSize.y");
        }
        j++;
      }
      if(WL < bestWL){
        bestWL = WL;
        bestPEArrayCenter.x = width;
        bestPEArrayCenter.y = 4.0*bufferSize.y;
        bestBufferCenters.swap(BufferCenters);
        BufferCenters.clear();
      }
      
    }
    printDebugInfo(__FUNCTION__, 12, "bestPEArrayCenter"+to_string(bestPEArrayCenter.x)+" "+to_string(bestPEArrayCenter.y));
    for(int i = 7;i < bufferNum+7;i++){
      cout<<"buffer "<<i<<" "<<bestBufferCenters[i].x<<" "<<bestBufferCenters[i].y<<endl;
      FPOS center;
      center.x = bestBufferCenters[i].x+bufferSize.x/2;
      center.y = bestBufferCenters[i].y +4.0*bufferSize.y;
      placeContainer(terminalInstance[i-7], center);
    }
    // exit(0);
    for(int i = 0;i < moduleCNT;i++){
      
        moduleInstance[i].center.x = bestPEArrayCenter.x;
        moduleInstance[i].center.y = bestPEArrayCenter.y;
      
    }
    // spaceCenter.x = bestPEArrayCenter.x;
    // spaceCenter.y = bestPEArrayCenter.y;
    plotTerm(scale);
    // exit(0);
    // place.end.x *= 30;
    // place.end.y *= 30;
}


void findIOCenter(FPOS &ioCenter){
  for(int i= 0;i <  terminalCNT;i++){
    std::string str(terminalInstance[i].Name());
    if(str.compare("w[3]") == 0){
      ioCenter.x = terminalInstance[i].center.x;
      ioCenter.y = terminalInstance[i].center.y;
      
      std::string str(terminalInstance[i].Name());
      printDebugInfo(__FUNCTION__, 0, "w[3] found");
      printDebugInfo(__FUNCTION__, 1, "ioCenter.x = "+to_string(ioCenter.x)+" ioCenter.y = "+to_string(ioCenter.y));
      break;
    }
    
    
  }
  // printDebugInfo(__FUNCTION__, 0, "not find w[3]");
}
void findPEArraySize(FPOS &peArraySize){
  prec area = 0;
  for(int i = 0;i < moduleCNT;i++){
    if(moduleInstance[i].type == "PE"){
      area += moduleInstance[i].area;
    }
  }
  peArraySize.x = sqrt(area);
  peArraySize.y = peArraySize.x;
  printDebugInfo(__FUNCTION__, 0, "peArraySize.x = "+to_string(peArraySize.x)+" peArraySize.y = "+to_string(peArraySize.y));
}
void findBufferSize(FPOS &bufferSize){
  bufferSize.x = terminalInstance[1].size.x;
  bufferSize.y = terminalInstance[1].size.y;
  printDebugInfo(__FUNCTION__, 0, "bufferSize.x = "+to_string(bufferSize.x)+" bufferSize.y = "+to_string(bufferSize.y));
}

void placeContainer(TERM &curContainer, FPOS location){
  curContainer.center.x = location.x;
  curContainer.center.y = location.y;
  // curContainer.LL.x = location.x - curContainer.size.x/2;
  // curContainer.LL.y = location.y - curContainer.size.y/2;
  // curContainer.UR.x = location.x + curContainer.size.x/2;
  // curContainer.UR.y = location.y + curContainer.size.y/2;
  printDebugInfo(__FUNCTION__, 0, "curContainer.center.x = "+to_string(curContainer.center.x)+" curContainer.center.y = "+to_string(curContainer.center.y));
  // printDebugInfo(__FUNCTION__, 1, "curContainer.LL.x = "+to_string(curContainer.LL.x)+" curContainer.LL.y = "+to_string(curContainer.LL.y));
  // printDebugInfo(__FUNCTION__, 2, "curContainer.UR.x = "+to_string(curContainer.UR.x)+" curContainer.UR.y = "+to_string(curContainer.UR.y));
}


void roughLegalization(int PENum){
  std::vector<prec> leg_grim;
  
  std::vector<int> grimIDX_list;
  grimIDX_list.resize(moduleCNT);
  int maxGrimDim = 1024;
  prec maxGrimDimprec = 1024.0;
  leg_grim.resize(maxGrimDim*maxGrimDim);
  FPOS dxy_grim;
  dxy_grim.x = place.end.x/maxGrimDimprec;
  dxy_grim.y = place.end.x/maxGrimDimprec;
  prec grim_area = dxy_grim.x*dxy_grim.y;
  //register the dust cell
  for(int i = 0;i < moduleCNT;i++){
    FPOS center;
    center.x = moduleInstance[i].center.GetX(); 
    center.y = moduleInstance[i].center.GetX(); 
    prec area = moduleInstance[i].area;
    int grimIDX = (int)(center.x/dxy_grim.GetX()) + (int)1024*(center.y/dxy_grim.GetY());
    grimIDX_list[i] = grimIDX;

    leg_grim[grimIDX] += area;

  }
  //register the macro cell
  for(int i = 0;i < PENum;i++)
  {
    FPOS center,LL,UR;
    
    LL.x = terminalInstance[i].pmin.GetX(); 
    LL.y = terminalInstance[i].pmin.GetY(); 
    UR.x = terminalInstance[i].pmax.GetX(); 
    UR.y = terminalInstance[i].pmax.GetY(); 
    
    POS occupyLength;
    occupyLength.x = (int)(UR.x - LL.x)/(dxy_grim.GetX());
    occupyLength.y = (int)(UR.y - LL.y)/(dxy_grim.GetY());
    int LLgrimIDX = (int)(LL.x/dxy_grim.GetX()) + (int)1024*(LL.y/dxy_grim.GetY());
    int base = LLgrimIDX;
    for(int dy=0;dy < occupyLength.y;dy++)
    {
      
      for(int dx = 0;dx < occupyLength.x;dx++)
      {
        int grimIDX = base + dx;
        leg_grim[grimIDX] = grim_area;
      }
      base += maxGrimDim;
    }
  }
  //move dust cells
  prec threhold = 0.9 * grim_area;
  int illegalNum = 0;
  int illegalNum_origin = 0;
  for(int i = 0;i < moduleCNT;i++){
    int grimIDX = grimIDX_list[i] ;
    prec area = moduleInstance[i].area;
    if(leg_grim[grimIDX] > threhold){
      illegalNum_origin++;
      if(leg_grim[grimIDX] >= grim_area)//overlap by macro
      {
        int length = 1;
        bool found = false;
        bool legalSearch = true;
        int foundIDX = -1;
        while(legalSearch && !found)
        {
          
          int LL,LR,UL,UR;
          LL = grimIDX - maxGrimDim*length - length;
          LR = grimIDX - maxGrimDim*length + length;
          UL = grimIDX + maxGrimDim*length - length;
          UR = grimIDX + maxGrimDim*length + length;
          if(LL<0 ||UR>leg_grim.size())
          {
            legalSearch = false;
            break;
          }
          for(int idx = LL;idx < LR;idx++){
            if(leg_grim[idx]+area<threhold){
              
              found = true;
              foundIDX = idx;
              break;
            }
          }
          for(int idx = LL;idx < UL;idx+=maxGrimDim){
            if(leg_grim[idx]+area<threhold){
              
              found = true;
              foundIDX = idx;
              break;
            }
          }
          for(int idx = UL;idx < UR;idx++){
            if(leg_grim[idx]+area<threhold){
              
              found = true;
              foundIDX = idx;
              break;
            }
          }
          for(int idx = LR;idx < UR;idx+=maxGrimDim){
            if(leg_grim[idx]+area<threhold){
              
              found = true;
              foundIDX = idx;
              break;
            }
          }
          length++;

        }
        if(legalSearch&&found){
          leg_grim[grimIDX] -= area;
          leg_grim[foundIDX] += area;
          grimIDX_list[i] = foundIDX;
          moduleInstance[i].center.x = dxy_grim.x*(prec)(foundIDX%maxGrimDim) ;
          moduleInstance[i].center.y = dxy_grim.y*(prec)(foundIDX/maxGrimDim) ;
        }
        if(!found){
          cout<<"illegal dust cell id: "<<i<<endl;
          illegalNum++;
        }
      }
    }
  }
  cout<<"original illegal rate  = "<<illegalNum_origin<<" / "<<moduleCNT<<endl;
  cout<<"illegal rate  = "<<illegalNum<<" / "<<moduleCNT<<endl;
}

void stitchingLegalization(int PENum)
{
  cout<<"debug 1 stitch"<<endl;
  FPOS PE_gravecenter;
  PE_gravecenter.x = 0;
  PE_gravecenter.y = 0;
  for(int i = 0; i < moduleCNT;i++)
  {
    PE_gravecenter.x += moduleInstance[i].center.x;
    PE_gravecenter.y += moduleInstance[i].center.y;
  }
  PE_gravecenter.x /= (prec)moduleCNT;
  PE_gravecenter.y /= (prec)moduleCNT;
  cout<<"debug 2 stitch"<<endl;
  FPOS bufferSize;
  findBufferSize(bufferSize);
  for(int i = 0; i < moduleCNT;i++)
  {
    
    moduleInstance[i].center.y += moduleInstance[i].center.y+bufferSize.y*3;
  }

  vector<vector<int>> grids_dust_cell_cnt;//-1 means no dust cell but terminal
  cout<<"debug 3 stitch"<<endl;

  int grid_size = 100;
  vector<vector<prec>> grid_used_area;
  for(int i = 0;i < grid_size;i++){
    vector<int> tmp;
    tmp.resize(grid_size);
    grids_dust_cell_cnt.push_back(tmp);
    vector<prec> tmp2;
    tmp2.resize(grid_size);
    grid_used_area.push_back(tmp2);
    
    for(int j = 0;j< grid_size;j++){
      grids_dust_cell_cnt[i][j] = 0;
      grid_used_area[i][j] = 0.0;
    }
  }
  cout<<"debug 4 stitch"<<endl;
  FPOS grid_interval;
  vector<POS> grid_center_of_module;
  grid_center_of_module.resize(moduleCNT);
  grid_interval.x = (prec)(place.end.x/(prec)grid_size);
  grid_interval.y = (prec)(place.end.y/(prec)grid_size);
  prec grid_area = grid_interval.x*grid_interval.y;
  cout<<"debug 5 stitch"<<endl;
  for(int i = 0;i < 16;i++){
    FPOS LL,UR;
    LL.x = terminalInstance[i].pmin.GetX();
    LL.y = terminalInstance[i].pmin.GetY();
    UR.x = terminalInstance[i].pmax.GetX();
    UR.y = terminalInstance[i].pmax.GetY();

    int LL_grid_x = (int)(LL.x/grid_interval.x);
    LL_grid_x = max(0,LL_grid_x);
    int LL_grid_y = (int)(LL.y/grid_interval.y);
    LL_grid_y = max(0,LL_grid_y);
    int UR_grid_x = (int)(UR.x/grid_interval.x);
    UR_grid_x = min(grid_size-1,UR_grid_x);
    int UR_grid_y = (int)(UR.y/grid_interval.y);
    UR_grid_y = min(grid_size-1,UR_grid_y);

    for(int x = LL_grid_x;x <= UR_grid_x;x++){
      for(int y = LL_grid_y;y <= UR_grid_y;y++){
        grids_dust_cell_cnt[x][y] = -1;
        grid_used_area[x][y] = grid_area;
      }
    }
  }
  cout<<"debug 6 stitch"<<endl;

  for(int i = 0;i < moduleCNT;i++){
    FPOS center;
    center.x = moduleInstance[i].center.x;
    center.y = moduleInstance[i].center.y;
    int grid_x = (int)(center.x/grid_interval.x);
    grid_x = max(0,grid_x);
    grid_x = min(grid_size-1,grid_x);

    int grid_y = (int)(center.y/grid_interval.y);
    grid_x = max(0,grid_y);
    grid_x = min(grid_size-1,grid_y);
    cout<<"debug 61 stitch"<<endl;
    if(grids_dust_cell_cnt[grid_x][grid_y] >= 0){
      grids_dust_cell_cnt[grid_x][grid_y]++;
      grid_used_area[grid_x][grid_y] += moduleInstance[i].area;
    }
    else{
      grids_dust_cell_cnt[grid_x][grid_y]--;
      moduleInstance[i].center.x = 0;
      moduleInstance[i].center.y = 0;
      moduleInstance[i].pmin.x = moduleInstance[i].center.x;
      moduleInstance[i].pmin.y = moduleInstance[i].center.y;
      moduleInstance[i].pmax.x = moduleInstance[i].center.x;
      moduleInstance[i].pmax.y = moduleInstance[i].center.y;
      moduleInstance[i].half_size.x = 0;
      moduleInstance[i].half_size.y = 0;
      moduleInstance[i].area = 0;
    }
    // cout<<"debug 62 stitch"<<endl;
    // POS id;
    // id.x = grid_x;
    // id.y = grid_y;
    // grid_center_of_module[i].x = grid_x;
    // grid_center_of_module[i].y = grid_y;
    // cout<<"debug 63 stitch"<<endl;
  }
  //find grid not covered by terminal from LL to UR
  // cout<<"debug 7 stitch"<<endl;
  // for(int i = 0;i < grid_size;i++){
  //   for(int j = 0;j< grid_size;j++){
  //     if(grids_dust_cell_cnt[i][j] >= 0){
  //       cout<<"debug 70 stitch"<<endl;
  //       if(grid_used_area[i][j] <=0.9*grid_area){
  //         //find the closest grid with dust cell
  //         int closest_grid_x = i;
  //         int closest_grid_y = j;
  //         int closest_grid_dist = 100000000;
  //         cout<<"debug 71 stitch"<<endl;
  //         for(int x=closest_grid_x;x< grid_size;x++){
  //           for(int y=closest_grid_y;y< grid_size;y++){
  //             if(grids_dust_cell_cnt[x][y] > 0){
  //               int dist = abs(x-i)+abs(y-j);
  //               if(dist < closest_grid_dist){
  //                 closest_grid_dist = dist;
  //                 closest_grid_x = x;
  //                 closest_grid_y = y;
  //               }
  //             }
  //           }
  //         }
  //         if(closest_grid_dist>200){
  //           continue;
  //         }
  //         cout<<"debug 72 stitch"<<endl;
  //         // move the module from the closest grid
  //         for(int x = 0;x < moduleCNT;x++)
  //         {
  //           if(grid_center_of_module[x].x == closest_grid_x && grid_center_of_module[x].y == closest_grid_y&&grid_used_area[i][j] <=0.9*grid_area){
  //             moduleInstance[x].center.x -= (prec)(closest_grid_x - i)*grid_interval.x;
  //             moduleInstance[x].center.y -= (prec)(closest_grid_y - j)*grid_interval.y;
  //             grid_center_of_module[x].x = i;
  //             grid_center_of_module[x].y = j;
  //             grid_used_area[i][j] += moduleInstance[x].area;
  //             grid_used_area[closest_grid_x][closest_grid_y] -= moduleInstance[x].area;
  //             break;
  //           }
  //         }
  //         cout<<"debug 73 stitch"<<endl;
  //       }
  //     }
      
  //   }
  // }



}