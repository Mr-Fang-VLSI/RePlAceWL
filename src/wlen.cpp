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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <omp.h>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <vector>

#include "replace_private.h"
#include "macro.h"
#include "opt.h"
#include "wlen.h"

#include "lefdefIO.h"
using namespace std;
using std::fstream;
using std::make_pair;
using std::max;
using std::min;
using std::stringstream;
#define debug_2pinnet 0
#define compare 1

int wcof_flg;
int MAX_EXP;
int NEG_MAX_EXP;
prec hpwl_mGP3D;
prec hpwl_mGP2D;
prec hpwl_cGP3D;
prec hpwl_cGP2D;
prec TSV_WEIGHT;

FPOS wcof00;
FPOS wcof00_org;
FPOS total_hpwl;
FPOS total_stnwl;  // lutong
FPOS total_wlen;
FPOS gp_wlen_weight;
FPOS dp_wlen_weight;
FPOS base_wcof;
FPOS wlen_cof;
FPOS wlen_cof_inv;
prec fastWLCPU;
EXP_ST *exp_st;

// extern vector<prec> ctrl_pts;
// extern vector<FPOS> ctrl_pt_grad;
void SetMAX_EXP_wlen() {
  MAX_EXP = 300;
  NEG_MAX_EXP = -300;
}

FPOS get_wlen_cof(prec ovf) {
  switch(wcof_flg) {
    case 1:
      return get_wlen_cof1(ovf);
      break;
    case 2:
      return get_wlen_cof2(ovf);
      break;
    default:
      return get_wlen_cof3(ovf);
      break;
  }
}

FPOS get_wlen_cof3(prec ovf) {
  FPOS cof;
  prec tmp = 0;
  if(ovf > 1.0) {
    cof.x = 1.0 / 10.0;
    cof.y = 1.0 / 10.0;
  }
  else if(ovf < 0.1) {
    cof.x = 1.0 / 0.01;
    cof.y = 1.0 / 0.01;
  }
  else {
    tmp = 1.0 / pow(10.0, (ovf - 0.1) * 10.0 / 3.0 - 2.0);
    cof.x = cof.y = tmp;
  }
  return cof;
}

FPOS get_wlen_cof1(prec ovf) {
  FPOS cof;

  if(ovf > 1.0)
    cof.x = cof.y = 1.0 / 10.0;
  else if(ovf < 0.1)
    cof.x = cof.y = 1.0 / 0.01;
  else
    cof.x = cof.y = 1.0 / pow(10.0, (ovf - 0.1) * 10.0 / 3.0 - 2.0);

  return cof;
}

//
// Described in ePlace-MS paper !!!!
//
//
FPOS get_wlen_cof2(prec ovf) {
  FPOS cof;
  prec tmp = 0.0;
  if(ovf > 1.0) {
    cof.x = 0.1;
    cof.y = 0.1;
  }
  else if(ovf < 0.1) {
    cof.x = 10.0;
    cof.y = 10.0;
  }
  else {
    tmp = 1.0 / pow(10.0, (ovf - 0.1) * 20 / 9.0 - 1.0);
    cof.x = cof.y = tmp;
  }
  //  cout << "ovfl: " << ovf << endl;
  //  cof.Dump("current_wlen_cof:");
  return cof;
}

void wlen_init() {
  //  int i = 0 /* ,cnt=exp_st_cnt */;
  /* prec interval = exp_interval ; */
  //  exp_st = (EXP_ST *)malloc(sizeof(EXP_ST) * exp_st_cnt);
  //  for(i = 0; i < exp_st_cnt; i++) {
  //    exp_st[i].x = (prec)i * exp_interval - MAX_EXP;
  //    exp_st[i].val = exp(exp_st[i].x);
  //    if(i > 0) {
  //      exp_st[i - 1].y_h = (exp_st[i].val - exp_st[i - 1].val) /
  //      exp_interval;
  //    }
  //  }

  gp_wlen_weight.x = gp_wlen_weight.y = 1.0;
  dp_wlen_weight.x = dp_wlen_weight.y = 1.0;
}

//
void wcof_init(FPOS bstp) {
  // 0.5*(~) = binSize;
  //
  base_wcof.x = wcof00.x / (0.5 * (bstp.x + bstp.y));
  base_wcof.y = wcof00.y / (0.5 * (bstp.x + bstp.y));

  wlen_cof = fp_scal(0.1, base_wcof);
  //  wlen_cof = base_wcof;
  wlen_cof_inv = fp_inv(wlen_cof);
}

prec get_wlen() {
  ////////////////////////////

  switch(WLEN_MODEL) {
    case WA:
      return get_wlen_wa();
      break;

    case LSE:
      return get_wlen_lse();
      break;
  }

  ////////////////////////////
  // return 0;
}

prec get_wlen_wa() {
  FPOS net_wlen;
  prec tot_wlen = 0;

  total_wlen.SetZero();

  for(int i = 0; i < netCNT; i++) {
    net_wlen = get_net_wlen_wa(&netInstance[i]);

    total_wlen.x += net_wlen.x;
    total_wlen.y += net_wlen.y;
  }

  tot_wlen = total_wlen.x + total_wlen.y;  // +

  return tot_wlen;
}

prec get_wlen_lse(void) {
  NET *net = NULL;
  FPOS net_wlen;
  prec tot_wlen = 0;

  total_wlen.SetZero();

  for(int i = 0; i < netCNT; i++) {
    net = &netInstance[i];
    if(net->pinCNTinObject <= 1)
      continue;
    net_wlen = get_net_wlen_lse(net);
    total_wlen.x += net_wlen.x;
    total_wlen.y += net_wlen.y;
  }

  tot_wlen = total_wlen.x + total_wlen.y;

  return tot_wlen;

  /* return  */
  /*   total_wlen.x * wlen_cof_inv.x * gp_wlen_weight.x + */
  /*   total_wlen.y * wlen_cof_inv.y * gp_wlen_weight.y +  */
}

FPOS get_net_wlen_wa(NET *net) {
  FPOS net_wlen;
  FPOS sum_num1 = net->sum_num1;
  FPOS sum_num2 = net->sum_num2;
  FPOS sum_denom1 = net->sum_denom1;
  FPOS sum_denom2 = net->sum_denom2;

  if(net->pinCNTinObject <= 1)
    return FPOS(0, 0);

  if(fastWL && net->pinCNTinObject == 2) {
    // 2-pin
    net_wlen.x = net->max_x - net->min_x;
    net_wlen.y = net->max_y - net->min_y;
  }
  else {
    // multiple pin
    net_wlen.x = sum_num1.x / sum_denom1.x - sum_num2.x / sum_denom2.x;
    net_wlen.y = sum_num1.y / sum_denom1.y - sum_num2.y / sum_denom2.y;
  }

  return net_wlen;
}

FPOS get_net_wlen_lse(NET *net) {
  FPOS sum1, sum2;
  FPOS fp, wlen;
  PIN *pin = NULL;

  for(int i = 0; i < net->pinCNTinObject; i++) {
    pin = net->pin[i];
    fp = pin->fp;

    sum1.x += get_exp(wlen_cof.x * fp.x);
    sum1.y += get_exp(wlen_cof.y * fp.y);

    sum2.x += get_exp(-1.0 * wlen_cof.x * fp.x);
    sum2.y += get_exp(-1.0 * wlen_cof.y * fp.y);
  }

  wlen.x = log(sum1.x) + log(sum2.x);
  wlen.y = log(sum1.y) + log(sum2.y);

  return wlen;
}

//
// this only calculate HPWL based on the stored value in NET's ure
prec GetHpwl() {
  total_hpwl.SetZero();
  total_stnwl.SetZero();

  for(int i = 0; i < netCNT; i++) {
    NET *curNet = &netInstance[i];
    if(curNet->pinCNTinObject <= 1)
      continue;

    total_hpwl.x += curNet->max_x - curNet->min_x;
    total_hpwl.y += curNet->max_y - curNet->min_y;

    //// lutong
    // total_stnwl.x += curNet->stn_cof * (curNet->max_x - curNet->min_x) ;
    // total_stnwl.y += curNet->stn_cof * (curNet->max_y - curNet->min_y) ;
  }

  return total_hpwl.x + total_hpwl.y;
}

//
// this calculate current NET's informations (min_xyz/max_xyz) & calculate HPWL
// simultaneously
prec UpdateNetAndGetHpwl() {
  FPOS pof, fp;
  total_hpwl.SetZero();

  for(int i = 0; i < netCNT; i++) {
    NET *curNet = &netInstance[i];

    curNet->min_x = curNet->terminalMin.x;
    curNet->min_y = curNet->terminalMin.y;

    curNet->max_x = curNet->terminalMax.x;
    curNet->max_y = curNet->terminalMax.y;

    // update min_xyz, max_xyz in NET's info
    for(int j = 0; j < curNet->pinCNTinObject; j++) {
      PIN *pin = curNet->pin[j];

      // only for modules
      if(pin->term) {
        continue;
      }

      MODULE *curModule = &moduleInstance[pin->moduleID];
      pin->fp.SetAdd(curModule->center, curModule->pof[pin->pinIDinModule]);

      curNet->min_x = min(curNet->min_x, pin->fp.x);
      curNet->min_y = min(curNet->min_y, pin->fp.y);

      curNet->max_x = max(curNet->max_x, pin->fp.x);
      curNet->max_y = max(curNet->max_y, pin->fp.y);
    }

    if(curNet->pinCNTinObject <= 1) {
      continue;
    }

    // calculate HPWL
    total_hpwl.x += curNet->max_x - curNet->min_x;
    total_hpwl.y += curNet->max_y - curNet->min_y;
  }

  return total_hpwl.x + total_hpwl.y;
}

void wlen_grad2(int cell_idx, FPOS *grad2) {
  switch(WLEN_MODEL) {
    case LSE:
      wlen_grad2_lse(cell_idx, grad2);
      break;

    case WA:
      wlen_grad2_wa(grad2);
      break;
  }
  return;
}

void wlen_grad(int cell_idx, FPOS *grad) {
  grad->SetZero();
#ifdef NO_WLEN
  return;
#endif

  switch(WLEN_MODEL) {
    case LSE:
      wlen_grad_lse(cell_idx, grad);
      break;

    case WA:
      wlen_grad_wa(cell_idx, grad);
      break;
  }
  grad->x *= -1.0 * gp_wlen_weight.x;
  grad->y *= -1.0 * gp_wlen_weight.y;
  // cout<<"grad->x == "<<grad->x<<endl;
  // cout<<"grad->y == "<<grad->y<<endl;
}

void wlen_grad2_wa(FPOS *grad) {
  grad->x = 1.0;
  grad->y = 1.0;
  return;
}

void wlen_grad2_lse(int cell_idx, FPOS *grad2) {
  FPOS net_grad2;
  CELL *cell = &gcell_st[cell_idx];
  NET *net = NULL;
  PIN *pin = NULL;

  grad2->SetZero();

  for(int i = 0; i < cell->pinCNTinObject; i++) {
    pin = cell->pin[i];
    net = &netInstance[pin->netID];

    if(net->pinCNTinObject <= 1)
      continue;

    get_net_wlen_grad2_lse(net, pin, &net_grad2);

    grad2->x += net_grad2.x;
    grad2->y += net_grad2.y;
  }
}

void wlen_grad_lse(int cell_idx, FPOS *grad) {
  FPOS net_grad;
  CELL *cell = &gcell_st[cell_idx];
  NET *net = NULL;
  PIN *pin = NULL;

  grad->SetZero();

  for(int i = 0; i < cell->pinCNTinObject; i++) {
    pin = cell->pin[i];
    net = &netInstance[pin->netID];

    if(net->pinCNTinObject <= 1)
      continue;

    get_net_wlen_grad_lse(net, pin, &net_grad);

    grad->x += net_grad.x;
    grad->y += net_grad.y;
  }
}

void wlen_grad_wa(int cell_idx, FPOS *grad) {
  CELL *cell = &gcell_st[cell_idx];
  PIN *pin = NULL;
  NET *net = NULL;
  FPOS net_grad;

  grad->SetZero();

  for(int i = 0; i < cell->pinCNTinObject; i++) {
    pin = cell->pin[i];
    net = &netInstance[pin->netID];

    if(net->pinCNTinObject <= 1)
      continue;
    bool timeon = true;
    double time = 0.0f;
    if(timeon) {
      // cout << "GetCost Grad: " << time << endl;
      time_start(&time);
    }
    
    
    get_net_wlen_grad_wa(pin->fp, net, pin, &net_grad);
    if(timeon && net->pinCNTinObject == 2) {
      time_end(&time);
      // cout << "GetCost Grad: " << time << endl;
      // grad_update_runtime += time;
      // lc_update_runtime += time;
      runtime_2pinnet += time;
      time_start(&time);
    }
    if(timeon && net->pinCNTinObject == 3) {
      time_end(&time);
      // cout << "GetCost Grad: " << time << endl;
      // grad_update_runtime += time;
      // lc_update_runtime += time;
      runtime_3pinnet += time;
      time_start(&time);
    }
    float curTimingWeight = netInstance[pin->netID].timingWeight;
    // Timing Control Parts
    // curTimingWeight = netWeightBase + min(max(0.0f, netWeightBound -
    // netWeightBase), curTimingWeight / netWeightScale);
    //
    // cout << "calNetWeight: " << curTimingWeight << endl;
    if(hasUnitNetWeight) {
      net_grad.x *= netWeight;
      net_grad.y *= netWeight;
    }
    else if(hasCustomNetWeight) {
      net_grad.x *= net->customWeight;
      net_grad.y *= net->customWeight;
    }
    else if(isTiming && netWeightApply && curTimingWeight > 0) {
      net_grad.x *= curTimingWeight;
      net_grad.y *= curTimingWeight;
    }

    grad->x += net_grad.x;
    grad->y += net_grad.y;
  }
}

// customWeight update functions:
void initCustomNetWeight(string netWeightFile) {
  PrintProcBegin("CustomNetWeightInit");
  // save netName / weight pair
  std::unordered_map< string, prec > tempMap;

  fstream fin(netWeightFile.c_str());
  string line;
  if(fin.is_open()) {
    while(fin.good()) {
      getline(fin, line);
      if(line.empty() || line[0] != '#') {
        char delimiter = ' ';
        int pos = line.find(delimiter);
        string field = line.substr(0, pos);
        string value = line.substr(pos + 1);
        stringstream ss(value);
        if(line == "")
          continue;
        tempMap[field] = atof(value.c_str());
        //        cout <<"Net " <<field <<" has weight " <<tempMap[field]
        //        <<endl;
      }
    }
    fin.close();
  }
  else {
    cout << "ERROR: Cannot open " << netWeightFile << endl;
    exit(1);
  }

  // fill in net->customWeight
  int customWeightCnt = 0;
  for(int i = 0; i < netCNT; i++) {
    if(tempMap.find(string(netInstance[i].Name())) != tempMap.end()) {
      netInstance[i].customWeight = tempMap[string(netInstance[i].Name())];
      customWeightCnt++;
      // cout << netInstance[i].Name()
      //  <<" weight " <<netInstance[i].customWeight
      //  <<" assigned" <<endl;
    }
  }
  PrintInfoInt("CustomNetWeightCount", customWeightCnt);
  PrintProcBegin("CustomNetWeightEnd");
}

void get_net_wlen_grad2_lse(NET *net, PIN *pin, FPOS *grad2) {
  POS flg1 = pin->flg1, flg2 = pin->flg2;
  FPOS e1 = pin->e1, e2 = pin->e2;
  FPOS grad2_1, grad2_2;
  FPOS sum_denom1 = net->sum_denom1;
  FPOS sum_denom2 = net->sum_denom2;

  if(flg1.x) {
    grad2_1.x = (e1.x) * (sum_denom1.x - e1.x) / (sum_denom1.x * sum_denom1.x);
  }

  if(flg2.x) {
    grad2_2.x = (e2.x) * (sum_denom2.x - e2.x) / (sum_denom2.x * sum_denom2.x);
  }

  grad2->x = grad2_1.x + grad2_2.x;

  if(flg1.y) {
    grad2_1.y = (e1.y) * (sum_denom1.y - e1.y) / (sum_denom1.y * sum_denom1.y);
  }

  if(flg2.y) {
    grad2_2.y = (e2.y) * (sum_denom2.y - e2.y) / (sum_denom2.y * sum_denom2.y);
  }

  grad2->y = grad2_1.y + grad2_2.y;

  return;
}

void get_net_wlen_grad_lse(NET *net, PIN *pin, FPOS *grad) {
  POS flg1 = pin->flg1, flg2 = pin->flg2;
  FPOS grad1, grad2;
  FPOS e1 = pin->e1, e2 = pin->e2;
  FPOS sum_denom1 = net->sum_denom1;
  FPOS sum_denom2 = net->sum_denom2;

  if(flg1.x) {
    grad1.x = e1.x / sum_denom1.x;
  }

  if(flg2.x) {
    grad2.x = e2.x / sum_denom2.x;
  }

  grad->x = grad1.x - grad2.x;

  if(flg1.y) {
    grad1.y = e1.y / sum_denom1.y;
  }

  if(flg2.y) {
    grad2.y = e2.y / sum_denom2.y;
  }

  grad->y = grad1.y - grad2.y;

  return;
}

// wlen_cof
// obj?
//
void get_net_wlen_grad_wa(FPOS obj, NET *net, PIN *pin, FPOS *grad) {
  if(fastWL && net->pinCNTinObject == 2) {
    // cout<<"fast gradient debug 0"<<endl;
    PIN *pin1;

    pin1 = net->pin[1 - pin->pinIDinNet];
    // cout<<"fast gradient debug 0.1"<<endl;

    FPOS fp0, fp1;
    // if(!pin->term)
    // {
    // fp0.x = moduleInstance[pin->moduleID].center.x +
    // moduleInstance[pin->moduleID].pof[pin->pinIDinModule].x; fp0.y =
    // moduleInstance[pin->moduleID].center.y +
    // moduleInstance[pin->moduleID].pof[pin->pinIDinModule].y;
    // }
    // else
    // {
    fp0 = pin->fp;
    // }
    // if(!pin1->term)
    // {
    // fp1.x = moduleInstance[pin1->moduleID].center.x +
    // moduleInstance[pin1->moduleID].pof[pin->pinIDinModule].x; fp1.y =
    // moduleInstance[pin1->moduleID].center.y +
    // moduleInstance[pin1->moduleID].pof[pin->pinIDinModule].y;
    // }
    // else
    // {
    fp1 = pin1->fp;
    // }

    // cout<<"Pin module id is"<<pin->moduleID;
    // cout<<" module id is at"<<moduleInstance[pin->moduleID].center.x<<"
    // "<<moduleInstance[pin->moduleID].center.y; cout<<" fp 0 "<<fp0.x<<"
    // "<<fp0.y<<endl; cout<<"Pin1 module id is"<<pin1->moduleID; cout<<" module
    // id is at"<<moduleInstance[pin1->moduleID].center.x<<"
    // "<<moduleInstance[pin1->moduleID].center.y; cout<<"fast gradient debug
    // 0.2"<<endl; cout<<" fp 1 "<<fp1.x<<" "<<fp1.y<<endl;
    // if(isnan(fp1.x)||isnan(fp0.x)||isnan(fp1.y)||isnan(fp0.y))
    // {
    //   exit(0);
    // }
    FPOS dist;
    dist.x = fp0.x - fp1.x;
    dist.y = fp0.y - fp1.y;
    int i = 0;
    // cout<<"ctrl pt num is "<<ctrl_pt_num<<endl;
    // if(std::isnan(dist.x))
    // {
    //   // cout<<"it is nan"<<endl;
    //   // cout<<dist.x<<endl;
    // }
    // else
    // prec pts_start = ctrl_pts[0].x;
    // prec pts_end = ctrl_pts[0].x;
    if(dist.x <= ctrl_pts[0].x) {
      grad->x = ctrl_pt_grad.at(0).x;
      // grad->x = -1.0007;
      // cout<<" < min grad x = "<<grad->x<<endl;
    }
    else if(dist.x >= ctrl_pts[ctrl_pt_num - 1].x) {
      grad->x = ctrl_pt_grad[ctrl_pt_num - 1].x;
      // grad->x = 1.0007;
      // cout<<" > max grad x = "<<grad->x<<endl;
    }
    else {
      // cout<<"fast gradient debug 1"<<endl;
      prec tot = FASTWL_HALF.x;

      // prec interval = FASTWL_INTERVAL.x;

      i = (int)((dist.x + tot) * FASTWL_INTERVAL.x);
      // cout<<"fast gradient debug 1.1 i = "<<i<<endl;
      // cout<<"fast gradient debug 1.1 dist = "<<dist.x<<endl;
      // cout<<"i = "<<i<<endl;

      // cout<<"module id"<<pin->moduleID<<endl;
      // grad->x = ctrl_pt_grad[i].x + (dist.x -
      // ctrl_pts[i])*(ctrl_pt_grad[i+1].x -
      // ctrl_pt_grad[i].x)/(ctrl_pts[i+1]-ctrl_pts[i]);
      grad->x = linearFuncX[i].x * dist.x + linearFuncX[i].y;
      // cout<<"fast gradient debug 1.2"<<endl;
      // if(i ==0 ||i ==20)
      // {
      //   //  cout<<"module id is" << pin->moduleID<< "at "<<
      //   moduleInstance[pin->moduleID].center.x<<"
      //   "<<moduleInstance[pin->moduleID].center.y;
      //   //   cout<<" dist x = "<<dist.x<<endl;
      //   //   cout<<"i = "<<i<<endl;
      //   //   cout<<"grad x = "<<grad->x<<endl;
      // }
      // if(grad->x > 3||grad->x <-3)
      // {
      //   cout<<"fast gradient debug 1.1 dist = "<<dist.x<<endl;
      //   cout<<"i = "<<i<<endl;
      //   cout<<"tot = "<<tot<<endl;
      //   cout<<"grad x  = "<<grad->x<<endl;
      //   cout<<"linearFuncX[i].x  = "<<linearFuncX[i].x<<endl;
      //   cout<<"linearFuncX[i].y  = "<<linearFuncX[i].y<<endl;
      //   cout<<"ctrl_pts[i].x  = "<<ctrl_pts[i].x <<endl;
      //   cout<<"ctrl_pts[i+1].x   = "<<ctrl_pts[i+1].x <<endl;
      // }
    }
    // cout<<"fast gradient debug 0.3"<<endl;
    // cout<<"fast gradient debug 2"<<endl;
    // if(std::isnan(dist.y))
    // {

    // }
    // else
    if(dist.y <= ctrl_pts[0].y) {
      grad->y = ctrl_pt_grad[0].y;
    }
    else if(dist.y >= ctrl_pts[ctrl_pt_num - 1].y) {
      grad->y = ctrl_pt_grad[ctrl_pt_num - 1].y;
    }
    else {
      // cout<<"fast gradient debug 3"<<endl;

      prec tot = FASTWL_HALF.y;

      // prec interval = FASTWL_INTERVAL.y;
      i = (int)((dist.y + tot) * FASTWL_INTERVAL.y);
      // grad->y = ctrl_pt_grad[i].y + (dist.y -
      // ctrl_pts[i])*(ctrl_pt_grad[i+1].y -
      // ctrl_pt_grad[i].y)/(ctrl_pts[i+1]-ctrl_pts[i]);
      grad->y = linearFuncY[i].x * dist.y + linearFuncY[i].y;
      // if(grad->y > 2||grad->y <-2)
      // {
      //   cout<<"fast gradient debug 3.1 dist = "<<dist.y<<endl;
      //   cout<<"i = "<<i<<endl;
      //   cout<<"tot = "<<tot<<endl;
      //   cout<<"grad y  = "<<grad->y<<endl;
      //   cout<<"k = "<<linearFuncY[i].x<<endl;
      //   cout<<"b = "<<linearFuncY[i].y<<endl;
      // }
    }
    /* disgard
    // for(int i = 0;i < ctrl_pt_num-1;++i)// can be improved : calculate the
    interval instead of for+if
    // {
    //   if(dist.x > ctrl_pts[i]&&dist.x > ctrl_pts[i+1])
    //   {
    //     grad->x = ctrl_pt_grad[i] + dist.x*(ctrl_pt_grad[i+1] -
    ctrl_pt_grad[i])/(ctrl_pts[i+1]-ctrl_pts[i]);
    //     break;rad
    //   }
    // }

    // for(int i = 0;i < ctrl_pt_num-1;++i)
    // {
    //   if(dist.y > ctrl_pts[i]&&dist.y > ctrl_pts[i+1])
    //   {
    //     grad->y = ctrl_pt_grad[i] + dist.y*(ctrl_pt_grad[i+1] -
    ctrl_pt_grad[i])/(ctrl_pts[i+1]-ctrl_pts[i]);
    //     break;
    //   }
    // }
    */

    // cout<<"fast gradient debug 1"<<endl;
  }
  else if(fastWL && net->pinCNTinObject == 3) {
    prec min_x = net->min_x;
    prec max_x = net->max_x;
    prec min_y = net->min_y;
    prec max_y = net->max_y;
    FPOS dist;
    dist.x = max_x - min_x;
    dist.y = max_y - min_y;

    if(dist.x < FASTWL_HALF.x) {
      int i = (int)((dist.x + FASTWL_HALF.x) * FASTWL_INTERVAL.x);
      // cout<<"3 pin i = "<<i<<endl;
      // cout<<"3 pin debug 0"<<endl;
      prec mid;
      for(int i = 0; i < net->pinCNTinObject; i++) {
        if(net->pin[i]->fp.x < max_x && net->pin[i]->fp.x > min_x) {
          mid = net->pin[i]->fp.x;
          break;
        }
      }
      mid = mid - min_x;
      // cout<<"3 pin debug 0"<<endl;
      if(pin->fp.x == max_x) {
        // cout<<"3 pin debug 1"<<endl;
        grad->x = (-linearFuncX_3pin[i].second.y * mid +
                     linearFuncX_3pin[i].first.x);
        // if(i < ctrl_pt_num - 2) {
        //   prec ratio = 1;
        //   // prec ratio = 1 - (dist.x - ctrl_pts[i].x) * FASTWL_INTERVAL.x;
        //   grad->x = (-linearFuncX_3pin[i].second.y * mid +
        //              linearFuncX_3pin[i].first.x) *
        //             (ratio);
        //   grad->x += (-linearFuncX_3pin[i + 1].second.y * mid +
        //               linearFuncX_3pin[i + 1].first.x) *
        //              (1 - ratio);
        // }
        // else {
        //   grad->x = (-linearFuncX_3pin[i].second.y * mid +
        //              linearFuncX_3pin[i].first.x);
        // }
      }
      else if(pin->fp.x == min_x) {
        // cout<<"3 pin debug 2"<<endl;
        grad->x = -(linearFuncX_3pin[i].second.y * mid +
                      linearFuncX_3pin[i].first.y);
        // if(i < ctrl_pt_num - 2) {
        //   prec ratio = 1;
        //   // prec ratio = 1 - (dist.x - ctrl_pts[i].x) * FASTWL_INTERVAL.x;
        //   grad->x = -(linearFuncX_3pin[i].second.y * mid +
        //               linearFuncX_3pin[i].first.y) *
        //             (ratio);
        //   grad->x += -(linearFuncX_3pin[i + 1].second.y * mid +
        //                linearFuncX_3pin[i + 1].first.y) *
        //              (1 - ratio);
        // }
        // else {
        //   grad->x = -(linearFuncX_3pin[i].second.y * mid +
        //               linearFuncX_3pin[i].first.y);
        // }
      }
      else {
        // cout<<"3 pin debug 3"<<endl;
        grad->x = (linearFuncX_3pin[i].second.x * mid -
                      linearFuncX_3pin[i].first.y);
        // if(i < ctrl_pt_num - 2) {
        //   prec ratio = 1;
        //   // prec ratio = 1 - (dist.x - ctrl_pts[i].x) * FASTWL_INTERVAL.x;
        //   grad->x = (linearFuncX_3pin[i].second.x * mid -
        //               linearFuncX_3pin[i].first.y) *
        //             (ratio);
        //   grad->x += (linearFuncX_3pin[i + 1].second.x * mid -
        //                linearFuncX_3pin[i + 1].first.y) *
        //              (1 - ratio);
        // }
        // else {
        //   grad->x = (linearFuncX_3pin[i].second.x * mid -
        //               linearFuncX_3pin[i].first.y);
        // }
      }
    }

    else {
      // cout<<"3 pin debug 4"<<endl;
      FPOS grad_sum_num1, grad_sum_num2;
      FPOS grad_sum_denom1, grad_sum_denom2;
      FPOS grad1;
      FPOS grad2;
      FPOS e1 = pin->e1;
      FPOS e2 = pin->e2;
      POS flg1 = pin->flg1;
      POS flg2 = pin->flg2;
      FPOS sum_num1 = net->sum_num1;
      FPOS sum_num2 = net->sum_num2;
      FPOS sum_denom1 = net->sum_denom1;
      FPOS sum_denom2 = net->sum_denom2;
      if(flg1.x) {
        grad_sum_denom1.x = wlen_cof.x * e1.x;
        grad_sum_num1.x = e1.x + obj.x * grad_sum_denom1.x;
        grad1.x =
            (grad_sum_num1.x * sum_denom1.x - grad_sum_denom1.x * sum_num1.x) /
            (sum_denom1.x * sum_denom1.x);
      }

      

      if(flg2.x) {
        grad_sum_denom2.x = wlen_cof.x * e2.x;
        grad_sum_num2.x = e2.x - obj.x * grad_sum_denom2.x;
        grad2.x =
            (grad_sum_num2.x * sum_denom2.x + grad_sum_denom2.x * sum_num2.x) /
            (sum_denom2.x * sum_denom2.x);
      }

      

      grad->x = grad1.x - grad2.x;
    }

    if(dist.y < FASTWL_HALF.y) {
      // cout<<"3 pin debug 5"<<endl;
      int i = (int)((dist.y + FASTWL_HALF.y) * FASTWL_INTERVAL.y);
      // cout<<"3 pin i = "<<i<<endl;
      prec mid;
      for(int i = 0; i < net->pinCNTinObject; i++) {
        if(net->pin[i]->fp.y < max_y && net->pin[i]->fp.y > min_y) {
          mid = net->pin[i]->fp.y;
          break;
        }
      }
      mid = mid - min_y;
      // cout<<"3 pin debug 6"<<endl;
      // cout<<"mid = "<<mid;
      if(pin->fp.y == max_y) {
        // cout<<"3 pin debug 6.1"<<endl;
        grad->y = (-linearFuncY_3pin[i].second.y * mid +
                     linearFuncY_3pin[i].first.x);
        // if(i < ctrl_pt_num - 2) {
        //   // cout<<"3 pin debug 6.2"<<endl;
        //   prec ratio = 1;
        //   // prec ratio = 1 - (dist.x - ctrl_pts[i].y) * FASTWL_INTERVAL.y;
        //   // cout<<"3 pin debug 6.21"<<endl;
        //   grad->y = (-linearFuncY_3pin[i].second.y * mid +
        //              linearFuncY_3pin[i].first.x) *
        //             (ratio);
        //   // cout<<"3 pin debug 6.22"<<endl;
        //   grad->y += (-linearFuncY_3pin[i + 1].second.y * mid +
        //               linearFuncY_3pin[i + 1].first.x) *
        //              (1 - ratio);
        //   // cout<<"3 pin debug 6.23"<<endl;
        //   // cout<<"grad y = "<<grad->y<<endl;;
        // }
        // else {
        //   // cout<<"3 pin debug 6.3"<<endl;
        //   grad->y = (-linearFuncY_3pin[i].second.y * mid +
        //              linearFuncY_3pin[i].first.x);
        // }
      }
      
      else if(pin->fp.y == min_y) {
        // cout<<"3 pin debug 7"<<endl;
        grad->y = -(linearFuncY_3pin[i].second.y * mid +
                      linearFuncY_3pin[i].first.y);
        // if(i < ctrl_pt_num - 2) {
        //   prec ratio = 1;
        //   // prec ratio = 1 - (dist.x - ctrl_pts[i].y) * FASTWL_INTERVAL.y;
        //   grad->y = -(linearFuncY_3pin[i].second.y * mid +
        //               linearFuncY_3pin[i].first.y) *
        //             (ratio);
        //   grad->y += -(linearFuncY_3pin[i + 1].second.y * mid +
        //                linearFuncY_3pin[i + 1].first.y) *
        //              (1 - ratio);
        // }
        // else {
        //   grad->y = -(linearFuncY_3pin[i].second.y * mid +
        //               linearFuncY_3pin[i].first.y);
        // }
      }
      else {
        // cout<<"3 pin debug 8"<<endl;
        grad->y = (linearFuncY_3pin[i].second.x * mid -
                      linearFuncY_3pin[i].first.y);
        // if(i < ctrl_pt_num - 2) {
        //   prec ratio = 1;
        //   // prec ratio = 1 - (dist.x - ctrl_pts[i].y) * FASTWL_INTERVAL.y;
        //   grad->y = (linearFuncY_3pin[i].second.x * mid -
        //               linearFuncY_3pin[i].first.y) *
        //             (ratio);
        //   grad->y += (linearFuncY_3pin[i + 1].second.x * mid -
        //                linearFuncY_3pin[i + 1].first.y) *
        //              (1 - ratio);
        // }
        // else {
        //   grad->y = (linearFuncY_3pin[i].second.x * mid -
        //               linearFuncY_3pin[i].first.y);
        // }
      }
    //   if(1)
    // {
    //   FPOS grad_ref;
    //   FPOS grad_sum_num1, grad_sum_num2;
    // FPOS grad_sum_denom1, grad_sum_denom2;
    // FPOS grad1;
    // FPOS grad2;
    // FPOS e1 = pin->e1;
    // FPOS e2 = pin->e2;
    // POS flg1 = pin->flg1;
    // POS flg2 = pin->flg2;
    // FPOS sum_num1 = net->sum_num1;
    // FPOS sum_num2 = net->sum_num2;
    // FPOS sum_denom1 = net->sum_denom1;
    // FPOS sum_denom2 = net->sum_denom2;
    

    // if(flg1.y) {
    //   grad_sum_denom1.y = wlen_cof.y * e1.y;
    //   grad_sum_num1.y = e1.y + obj.y * grad_sum_denom1.y;
    //   grad1.y =
    //       (grad_sum_num1.y * sum_denom1.y - grad_sum_denom1.y * sum_num1.y) /
    //       (sum_denom1.y * sum_denom1.y);
    // }

    

    // if(flg2.y) {
    //   grad_sum_denom2.y = wlen_cof.y * e2.y;
    //   grad_sum_num2.y = e2.y - obj.y * grad_sum_denom2.y;
    //   grad2.y =
    //       (grad_sum_num2.y * sum_denom2.y + grad_sum_denom2.y * sum_num2.y) /
    //       (sum_denom2.y * sum_denom2.y);
    // }

    
    // grad_ref.y = grad1.y - grad2.y;
    // cout<<"diff = "<<(grad_ref.y - grad->y)/(grad_ref.y+0.001)<<endl;
    // cout<<"ref grad = "<<grad_ref.y<<", grad = "<<grad->y<<endl;
    // }
    }
    else {
      // cout<<"3 pin debug 9"<<endl;
    FPOS grad_sum_num1, grad_sum_num2;
    FPOS grad_sum_denom1, grad_sum_denom2;
    FPOS grad1;
    FPOS grad2;
    FPOS e1 = pin->e1;
    FPOS e2 = pin->e2;
    POS flg1 = pin->flg1;
    POS flg2 = pin->flg2;
    FPOS sum_num1 = net->sum_num1;
    FPOS sum_num2 = net->sum_num2;
    FPOS sum_denom1 = net->sum_denom1;
    FPOS sum_denom2 = net->sum_denom2;
    

    if(flg1.y) {
      grad_sum_denom1.y = wlen_cof.y * e1.y;
      grad_sum_num1.y = e1.y + obj.y * grad_sum_denom1.y;
      grad1.y =
          (grad_sum_num1.y * sum_denom1.y - grad_sum_denom1.y * sum_num1.y) /
          (sum_denom1.y * sum_denom1.y);
    }

    

    if(flg2.y) {
      grad_sum_denom2.y = wlen_cof.y * e2.y;
      grad_sum_num2.y = e2.y - obj.y * grad_sum_denom2.y;
      grad2.y =
          (grad_sum_num2.y * sum_denom2.y + grad_sum_denom2.y * sum_num2.y) /
          (sum_denom2.y * sum_denom2.y);
    }

    
    grad->y = grad1.y - grad2.y;

    }
    
  }
  else {
    FPOS grad_sum_num1, grad_sum_num2;
    FPOS grad_sum_denom1, grad_sum_denom2;
    FPOS grad1;
    FPOS grad2;
    FPOS e1 = pin->e1;
    FPOS e2 = pin->e2;
    POS flg1 = pin->flg1;
    POS flg2 = pin->flg2;
    FPOS sum_num1 = net->sum_num1;
    FPOS sum_num2 = net->sum_num2;
    FPOS sum_denom1 = net->sum_denom1;
    FPOS sum_denom2 = net->sum_denom2;
    if(flg1.x) {
      grad_sum_denom1.x = wlen_cof.x * e1.x;
      grad_sum_num1.x = e1.x + obj.x * grad_sum_denom1.x;
      grad1.x =
          (grad_sum_num1.x * sum_denom1.x - grad_sum_denom1.x * sum_num1.x) /
          (sum_denom1.x * sum_denom1.x);
    }

    if(flg1.y) {
      grad_sum_denom1.y = wlen_cof.y * e1.y;
      grad_sum_num1.y = e1.y + obj.y * grad_sum_denom1.y;
      grad1.y =
          (grad_sum_num1.y * sum_denom1.y - grad_sum_denom1.y * sum_num1.y) /
          (sum_denom1.y * sum_denom1.y);
    }

    if(flg2.x) {
      grad_sum_denom2.x = wlen_cof.x * e2.x;
      grad_sum_num2.x = e2.x - obj.x * grad_sum_denom2.x;
      grad2.x =
          (grad_sum_num2.x * sum_denom2.x + grad_sum_denom2.x * sum_num2.x) /
          (sum_denom2.x * sum_denom2.x);
    }

    if(flg2.y) {
      grad_sum_denom2.y = wlen_cof.y * e2.y;
      grad_sum_num2.y = e2.y - obj.y * grad_sum_denom2.y;
      grad2.y =
          (grad_sum_num2.y * sum_denom2.y + grad_sum_denom2.y * sum_num2.y) /
          (sum_denom2.y * sum_denom2.y);
    }

    grad->x = grad1.x - grad2.x;
    grad->y = grad1.y - grad2.y;
  }
}

void net_update_init(void) {
  if(!fastWL) {
    for(int i = 0; i < netCNT; i++) {
      NET *net = &netInstance[i];

      net->terminalMin.Set(place.end);
      net->terminalMax.Set(place.org);

      bool first_term = true;

      for(int j = 0; j < net->pinCNTinObject; j++) {
        PIN *pin = net->pin[j];
        if(pin->term) {
          if(first_term) {
            first_term = false;

            net->terminalMin.Set(pin->fp);
            net->terminalMax.Set(pin->fp);
          }
          else {
            net->terminalMin.Min(pin->fp);
            net->terminalMax.Max(pin->fp);
          }
        }
      }
    }
  }
  else {
    for(int i = 0; i < netCNT; i++) {
      NET *net = &netInstance[i];
      if(net->pinCNTinObject > 2) {
        net->terminalMin.Set(place.end);
        net->terminalMax.Set(place.org);
        net->terminal2ndMin.Set(place.end);
        net->terminal2ndMax.Set(place.org);
        bool first_term = true;

        for(int j = 0; j < net->pinCNTinObject; j++) {
          PIN *pin = net->pin[j];
          if(pin->term) {
            if(first_term) {
              first_term = false;

              net->terminalMin.Set(pin->fp);
              net->terminalMax.Set(pin->fp);
            }
            else {
              // net->terminalMin.Min(pin->fp);
              // net->terminalMax.Max(pin->fp);
              // compare to max and min
              if(pin->fp.x > net->terminalMax.x) {
                net->terminal2ndMax.x = net->terminalMax.x;
                net->terminalMax.x = pin->fp.x;
              }
              else if(pin->fp.x > net->terminal2ndMax.x) {
                net->terminal2ndMax.x = pin->fp.x;
              }
              if(pin->fp.y > net->terminalMax.y) {
                net->terminal2ndMax.y = net->terminalMax.y;
                net->terminalMax.y = pin->fp.y;
              }
              else if(pin->fp.y > net->terminal2ndMax.y) {
                net->terminal2ndMax.y = pin->fp.y;
              }

              // for min
              if(pin->fp.x < net->terminalMin.x) {
                net->terminal2ndMin.x = net->terminalMin.x;
                net->terminalMin.x = pin->fp.x;
              }
              else if(pin->fp.x < net->terminal2ndMin.x) {
                net->terminal2ndMin.x = pin->fp.x;
              }
              if(pin->fp.y < net->terminalMin.y) {
                net->terminal2ndMin.y = net->terminalMin.y;
                net->terminalMin.y = pin->fp.y;
              }
              else if(pin->fp.y < net->terminal2ndMin.y) {
                net->terminal2ndMin.y = pin->fp.y;
              }
            }
          }
        }
        // cout<<"2nd max terminal:"<<net->terminal2ndMax.x<<",
        // "<<net->terminal2ndMax.y<<endl; cout<<"2nd min
        // terminal:"<<net->terminal2ndMin.x<<", "<<net->terminal2ndMin.y<<endl;
      }
      else {
        net->terminalMin.Set(place.end);
        net->terminalMax.Set(place.org);

        bool first_term = true;

        for(int j = 0; j < net->pinCNTinObject; j++) {
          PIN *pin = net->pin[j];
          if(pin->term) {
            if(first_term) {
              first_term = false;

              net->terminalMin.Set(pin->fp);
              net->terminalMax.Set(pin->fp);
            }
            else {
              net->terminalMin.Min(pin->fp);
              net->terminalMax.Max(pin->fp);
            }
          }
        }
      }
    }
  }
}
void fastWL_init(int controlNum)  // controlNum must be odd
{
  // std::cout<<" size "<< ctrl_pt_grad.size()<<endl;
  // std::cout<<"fastwl_init debug 0 "<<controlNum<<" size "<<
  // ctrl_pt_grad.size()<<endl;
  ctrl_pt_grad.clear();
  ctrl_pts.clear();
  ctrl_pt_num = controlNum;
  // free(ctrl_pts);
  // free(ctrl_pt_grad);
  // ctrl_pt_grad = (struct FPOS*)malloc(controlNum*sizeof(FPOS));
  // ctrl_pts = (prec*)malloc(controlNum*sizeof(prec));

  // ctrl_pt_grad.resize(ctrl_pt_num);
  // ctrl_pts.resize(ctrl_pt_num);
  // std::cout<<" size "<< ctrl_pt_grad.size()<<endl;
  // set ctrl_pts,
  prec intervalX = (prec)2.56 / wlen_cof.GetX() / ((prec)controlNum - 1.0);
  prec totX = intervalX * ((prec)controlNum - 1.0);
  prec intervalY = (prec)2.56 / wlen_cof.GetX() / ((prec)controlNum - 1.0);
  prec totY = intervalY * ((prec)controlNum - 1.0);

  // prec intervalX = (prec)25.6/wlen_cof.GetX() /((prec)controlNum-1.0) ;
  // prec totX = intervalX*((prec)controlNum-1.0);
  // prec intervalY = (prec)25.6/wlen_cof.GetX() /((prec)controlNum-1.0) ;
  // prec totY = intervalY*((prec)controlNum-1.0);
  FASTWL_TOT.x = totX;
  FASTWL_TOT.y = totY;
  FASTWL_HALF.x = totX * 0.5;
  FASTWL_HALF.y = totY * 0.5;
  FASTWL_INTERVAL.x = 1.0 / intervalX;
  FASTWL_INTERVAL.y = 1.0 / intervalY;
  // std::cout<<"fastwl_init debug 0.1 "<<endl;
  FPOS tmp;
  tmp.x = -totX * 0.5;
  tmp.y = -totY * 0.5;
  ctrl_pts.push_back(tmp);
  // std::cout<<"fastwl_init debug 0.2 "<< ctrl_pts.at(0)<<endl;
  FPOS zerofp;
  zerofp.SetZero();
  ctrl_pt_grad.push_back(zerofp);
  // std::cout<<"fastwl_init debug 0.3 "<<endl;
  // cout<<"fastwl_init debug 1"<<endl;
  // cout<<ctrl_pts[0]<<endl;
  // std::cout<<" size "<< ctrl_pt_grad.size()<<endl;
  for(int i = 1; i < controlNum; i++) {
    FPOS tmp2;
    tmp2.x = ctrl_pts[i - 1].x + intervalX;
    tmp2.y = ctrl_pts[i - 1].y + intervalY;
    ctrl_pts.push_back(tmp2);
    ctrl_pt_grad.push_back(zerofp);
    ctrl_pt_grad[i].x = ctrl_pt_grad[i].y = 0.0;
    // linearFunc.push_back(zerofp);
    // cout<<ctrl_pts.at(i).x<<endl;
    // cout<<ctrl_pts.at(i).y<<endl;
  }
  for(int i = 1; i < controlNum - 1; i++) {
    linearFuncX.push_back(zerofp);
    linearFuncY.push_back(zerofp);
  }
  // cout<<"fastwl_init debug 2"<< ctrl_pt_grad.size()<<":"<<endl;
  for(int i = 0; i < ctrl_pt_num; i++) {
    // cout<<"fastwl_update debug 0.1"<<" i = "<<i<<endl;
    if(i != (ctrl_pt_num - 1) / 2) {
      // cout<<"fastwl_update debug 0.2"<<" i = "<<i<<endl;
      // cout<<"fastwl_update debug 0.2"<<" ctrl pt "<<ctrl_pts[i]<<endl;
      // cout<<"pts size = "<<ctrl_pts.size()<<endl;
      // cout<<"fastwl_update debug 0.2"<<" ctrl pt grad
      // "<<ctrl_pt_grad[i].x<<endl;
      cal_2pin_WA_grads(ctrl_pts[i], ctrl_pt_grad[i]);
      // cout<<"fastwl_update debug 0.2"<<i<<endl;
    }
    else {
      ctrl_pt_grad[i].SetZero();
    }
    // cout<<"grad for 2-pin dist = "<<ctrl_pts[i].x<<" "<<ctrl_pts[i].y <<" =
    // "<<ctrl_pt_grad[i].x<<" "<<ctrl_pt_grad[i].y<<endl; cout<<"begin linear
    // Func update"<<endl;

    // cout<<"linear Func updated"<<endl;
  }
  fastWL_init_3pin(controlNum);
  linearFunc_update();
}
void fastWL_init_3pin(int controlNum) {
  FPOS zerofp;
  pair< FPOS, FPOS > zeroPair;
  zeroPair.first = zerofp;
  zeroPair.second = zerofp;
  zerofp.SetZero();
  ctrl_pt_grad_3pin.clear();
  ctrl_pt_grad_3pin.push_back(zerofp);

  for(int i = 1; i < controlNum; i++) {
    ctrl_pt_grad_3pin.push_back(zerofp);
    linearFuncX_3pin.push_back(zeroPair);
    linearFuncY_3pin.push_back(zeroPair);
  }
  for(int i = 0; i < controlNum; i++) {
    if(i != (ctrl_pt_num - 1) / 2) {
      cal_3pin_WA_grads(ctrl_pts[i], ctrl_pt_grad_3pin[i]);
    }
    else {
      ctrl_pt_grad_3pin[i].SetZero();
    }
    cout << "3 pin grad at " << ctrl_pts[i].x << " = " << ctrl_pt_grad_3pin[i].x
         << " " << ctrl_pt_grad_3pin[i].y << endl;
  }
  linearFunc_update_3pin();
}
void cal_3pin_WA_grads(FPOS dist, FPOS &grad) {
  prec max_x = dist.x > 0.0 ? dist.x : 0.0;
  prec min_x = dist.x < 0.0 ? dist.x : 0.0;
  prec max_y = dist.y > 0.0 ? dist.y : 0.0;
  prec min_y = dist.y < 0.0 ? dist.y : 0.0;
  PIN pin1, pin2, pin3;
  FPOS sum_num1, sum_denom1, sum_num2, sum_denom2;
  sum_num1.SetZero();
  sum_num2.SetZero();
  sum_denom1.SetZero();
  sum_denom2.SetZero();

  for(int j = 0; j < 3; j++) {
    // PIN *pin = net->pin[j];
    FPOS fp;
    PIN *pin;

    if(j == 0) {
      pin = &pin1;
    }
    else if(j == 1) {
      pin = &pin2;
    }
    else {
      pin = &pin3;
    }
    if(j < 2) {
      fp.SetZero();
    }
    else {
      fp.Set(dist);
      cout << "cal 3pin fp :=" << fp.x << " " << fp.y << endl;
    }
    prec exp_max_x = (fp.x - max_x) * wlen_cof.x;
    prec exp_min_x = (min_x - fp.x) * wlen_cof.x;
    prec exp_max_y = (fp.y - max_y) * wlen_cof.y;
    prec exp_min_y = (min_y - fp.y) * wlen_cof.y;

    if(exp_max_x > NEG_MAX_EXP) {
      pin->e1.x = get_exp(exp_max_x);
      sum_num1.x += fp.x * pin->e1.x;
      sum_denom1.x += pin->e1.x;
      pin->flg1.x = 1;
    }
    else {
      pin->flg1.x = 0;
    }

    if(exp_min_x > NEG_MAX_EXP) {
      pin->e2.x = get_exp(exp_min_x);
      sum_num2.x += fp.x * pin->e2.x;
      sum_denom2.x += pin->e2.x;
      pin->flg2.x = 1;
    }
    else {
      pin->flg2.x = 0;
    }

    if(exp_max_y > NEG_MAX_EXP) {
      pin->e1.y = get_exp(exp_max_y);
      sum_num1.y += fp.y * pin->e1.y;
      sum_denom1.y += pin->e1.y;
      pin->flg1.y = 1;
    }
    else {
      pin->flg1.y = 0;
    }

    if(exp_min_y > NEG_MAX_EXP) {
      pin->e2.y = get_exp(exp_min_y);
      sum_num2.y += fp.y * pin->e2.y;
      sum_denom2.y += pin->e2.y;
      pin->flg2.y = 1;
    }
    else {
      pin->flg2.y = 0;
    }
  }

  FPOS grad_sum_num1, grad_sum_num2;
  FPOS grad_sum_denom1, grad_sum_denom2;
  FPOS grad1;
  FPOS grad2;
  FPOS e1 = pin3.e1;
  FPOS e2 = pin3.e2;
  POS flg1 = pin3.flg1;
  POS flg2 = pin3.flg2;

  if(flg1.x) {
    grad_sum_denom1.x = wlen_cof.x * e1.x;
    grad_sum_num1.x = e1.x + dist.x * grad_sum_denom1.x;
    grad1.x =
        (grad_sum_num1.x * sum_denom1.x - grad_sum_denom1.x * sum_num1.x) /
        (sum_denom1.x * sum_denom1.x);
  }

  if(flg1.y) {
    grad_sum_denom1.y = wlen_cof.y * e1.y;
    grad_sum_num1.y = e1.y + dist.y * grad_sum_denom1.y;
    grad1.y =
        (grad_sum_num1.y * sum_denom1.y - grad_sum_denom1.y * sum_num1.y) /
        (sum_denom1.y * sum_denom1.y);
  }

  if(flg2.x) {
    grad_sum_denom2.x = wlen_cof.x * e2.x;
    grad_sum_num2.x = e2.x - dist.x * grad_sum_denom2.x;
    grad2.x =
        (grad_sum_num2.x * sum_denom2.x + grad_sum_denom2.x * sum_num2.x) /
        (sum_denom2.x * sum_denom2.x);
  }

  if(flg2.y) {
    grad_sum_denom2.y = wlen_cof.y * e2.y;
    grad_sum_num2.y = e2.y - dist.y * grad_sum_denom2.y;
    grad2.y =
        (grad_sum_num2.y * sum_denom2.y + grad_sum_denom2.y * sum_num2.y) /
        (sum_denom2.y * sum_denom2.y);
  }

  grad.x = grad1.x - grad2.x;
  grad.y = grad1.y - grad2.y;
}

void linearFunc_update_3pin() {
  for(int i = 0; i < ctrl_pt_num - 1; i++) {
    //截距
    if(i != (ctrl_pt_num - 1) / 2)
    {
      linearFuncX_3pin[i].first.x = ctrl_pt_grad_3pin[i + 1].x;
      linearFuncX_3pin[i].first.y = ctrl_pt_grad_3pin[i + 1].x * 0.5;

      linearFuncY_3pin[i].first.x = ctrl_pt_grad_3pin[i + 1].y;
      linearFuncY_3pin[i].first.y = ctrl_pt_grad_3pin[i + 1].y * 0.5;
    //斜率
      linearFuncX_3pin[i].second.x = ctrl_pt_grad_3pin[i + 1].x / ctrl_pts[i].x;
      linearFuncX_3pin[i].second.y = linearFuncX_3pin[i].second.x * 0.5;

      linearFuncY_3pin[i].second.x = ctrl_pt_grad_3pin[i + 1].y / ctrl_pts[i].y;
      linearFuncY_3pin[i].second.y = linearFuncY_3pin[i].second.x * 0.5;
    }
    else{
      linearFuncX_3pin[i].first.x = 0.0;
      linearFuncX_3pin[i].first.y = 0.0;

      linearFuncY_3pin[i].first.x = 0.0;
      linearFuncY_3pin[i].first.y = 0.0;
    //斜率
      linearFuncX_3pin[i].second.x = 0.0;
      linearFuncX_3pin[i].second.y = 0.0;

      linearFuncY_3pin[i].second.x = 0.0;
      linearFuncY_3pin[i].second.y = 0.0;
    }
    
  }
}
void fastWL_update() {
  // update ctgrl_pt_grad with new wcof
  //  cout<<"fastwl_update debug 0"<<endl;
  //  fastWL_init(ctrl_pt_num);
  //  cout<<"fastwl_update debug 0.1 size = "<<ctrl_pts.size()<<endl;
  //  cout<<"fastwl_update debug 0.1 size = "<<ctrl_pt_grad.size()<<endl;
  //  for(int i = 0;i < ctrl_pts.size();i++)
  //  {
  //    cout<<ctrl_pts.at(i)<<endl;
  //    // cout<<ctrl_pt_grad[i].x<<endl;
  //  }
  //  for(int i = 0;i < ctrl_pt_num;i++)
  //  {
  //    // cout<<"fastwl_update debug 0.1"<<" i = "<<i<<endl;
  //    if(i!=(ctrl_pt_num-1)/2)
  //    {
  //      // cout<<"fastwl_update debug 0.2"<<" i = "<<i<<endl;
  //      // cout<<"fastwl_update debug 0.2"<<" ctrl pt "<<ctrl_pts[i]<<endl;
  //      // cout<<"pts size = "<<ctrl_pts.size()<<endl;
  //      // cout<<"fastwl_update debug 0.2"<<" ctrl pt grad
  //      "<<ctrl_pt_grad[i].x<<endl;
  //      cal_2pin_WA_grads(ctrl_pts[i],ctrl_pt_grad[i]);
  //      // cout<<"fastwl_update debug 0.2"<<i<<endl;

  //   }
  //   else{
  //     ctrl_pt_grad[i].SetZero();
  //   }
  //   // cout<<"grad for 2-pin dist = "<<ctrl_pts[i].x<<" "<<ctrl_pts[i].y <<"
  //   = "<<ctrl_pt_grad[i].x<<" "<<ctrl_pt_grad[i].y<<endl;
  //   // cout<<"begin linear Func update"<<endl;

  //   // cout<<"linear Func updated"<<endl;
  // }
  // linearFunc_update();
  // // cout<<"fastwl_update debug 1"<<endl;

  ctrl_pts.clear();

  // free(ctrl_pts);
  // free(ctrl_pt_grad);
  // ctrl_pt_grad = (struct FPOS*)malloc(controlNum*sizeof(FPOS));
  // ctrl_pts = (prec*)malloc(controlNum*sizeof(prec));

  // ctrl_pt_grad.resize(ctrl_pt_num);
  // ctrl_pts.resize(ctrl_pt_num);
  // std::cout<<" size "<< ctrl_pt_grad.size()<<endl;
  // set ctrl_pts,
  prec intervalX = (prec)2.56 / wlen_cof.GetX() / ((prec)ctrl_pt_num - 1.0);
  prec totX = intervalX * ((prec)ctrl_pt_num - 1.0);
  prec intervalY = (prec)2.56 / wlen_cof.GetX() / ((prec)ctrl_pt_num - 1.0);
  prec totY = intervalY * ((prec)ctrl_pt_num - 1.0);
  // prec intervalX = (prec)25.6/wlen_cof.GetX() /((prec)ctrl_pt_num -1.0) ;
  // prec totX = intervalX*((prec)ctrl_pt_num -1.0);
  // prec intervalY = (prec)25.6/wlen_cof.GetX() /((prec)ctrl_pt_num -1.0) ;
  // prec totY = intervalY*((prec)ctrl_pt_num -1.0);
  FASTWL_TOT.x = totX;
  FASTWL_TOT.y = totY;
  FASTWL_HALF.x = totX * 0.5;
  FASTWL_HALF.y = totY * 0.5;
  FASTWL_INTERVAL.x = 1.0 / intervalX;
  FASTWL_INTERVAL.y = 1.0 / intervalY;
  // std::cout<<"fastwl_init debug 0.1 "<<endl;
  FPOS tmp;
  tmp.x = -totX * 0.5;
  tmp.y = -totY * 0.5;
  ctrl_pts.push_back(tmp);
  // std::cout<<"fastwl_init debug 0.2 "<< ctrl_pts.at(0)<<endl;
  FPOS zerofp;
  zerofp.SetZero();
  // ctrl_pt_grad.push_back(zerofp);
  // std::cout<<"fastwl_init debug 0.3 "<<endl;
  // cout<<"fastwl_init debug 1"<<endl;
  // cout<<ctrl_pts[0]<<endl;
  // std::cout<<" size "<< ctrl_pt_grad.size()<<endl;
  for(int i = 1; i < ctrl_pt_num; i++) {
    FPOS tmp2;
    tmp2.x = ctrl_pts[i - 1].x + intervalX;
    tmp2.y = ctrl_pts[i - 1].y + intervalY;
    ctrl_pts.push_back(tmp2);
    cal_2pin_WA_grads(ctrl_pts[i], ctrl_pt_grad[i]);
    // ctrl_pt_grad.push_back(zerofp);
    // ctrl_pt_grad[i].x = ctrl_pt_grad[i].y = 0.0;
    // linearFunc.push_back(zerofp);
    // cout<<ctrl_pts.at(i).x<<endl;
    // cout<<ctrl_pts.at(i).y<<endl;
  }
  linearFunc_update();
  linearFunc_update_3pin();
}
void linearFunc_update() {
  FPOS inteval;
  // inteval.x = 1.0/FASTWL_INTERVAL.x;
  // inteval.y = 1.0/FASTWL_INTERVAL.y;
  // cout<<"linear Func_update interval is"<<inteval<<endl;
  for(int i = 0; i < ctrl_pt_num - 1; i++) {
    linearFuncX[i].x =
        (ctrl_pt_grad[i + 1].x - ctrl_pt_grad[i].x) * FASTWL_INTERVAL.x;
    linearFuncX[i].y = ctrl_pt_grad[i].x - (linearFuncX[i].x * ctrl_pts[i].x);
    linearFuncY[i].x =
        (ctrl_pt_grad[i + 1].y - ctrl_pt_grad[i].y) * FASTWL_INTERVAL.y;
    linearFuncY[i].y = ctrl_pt_grad[i].y - (linearFuncY[i].x * ctrl_pts[i].y);
    // cout<<"grad is"<<ctrl_pt_grad[i+1].x<<" "<<ctrl_pt_grad[i].x<<endl;
    // cout<<"ipos="<<ctrl_pts[i].x<<", xk = "<<linearFuncX[i].x<<",xb
    // ="<<linearFuncX[i].y<<" interval x"<<inteval.x<<endl; cout<<"ipos
    // ="<<ctrl_pts[i].y<<", yk = "<<linearFuncY[i].x<<",yb
    // ="<<linearFuncY[i].y<<" interval y"<<inteval.y<<endl;
  }
}
void cal_2pin_WA_grads(FPOS dist, FPOS &grad) {
  // cout<<"cal_2pin_wa debug 0"<<endl;
  FPOS sum_num1, sum_num2;
  FPOS sum_denom1, sum_denom2;
  prec max_x, min_x, max_y, min_y;
  max_x = dist.x > 0.0 ? dist.x : 0.0;
  min_x = dist.x < 0.0 ? dist.x : 0.0;
  max_y = dist.y > 0.0 ? dist.y : 0.0;
  min_y = dist.y < 0.0 ? dist.y : 0.0;
  prec exp_max_x = (dist.x - max_x) * wlen_cof.x;
  prec exp_min_x = (min_x - dist.x) * wlen_cof.x;
  prec exp_max_y = (dist.y - max_y) * wlen_cof.y;
  prec exp_min_y = (min_y - dist.y) * wlen_cof.y;

  prec exp_max_x_0 = (0 - max_x) * wlen_cof.x;
  prec exp_min_x_0 = (min_x - 0) * wlen_cof.x;
  prec exp_max_y_0 = (0 - max_y) * wlen_cof.y;
  prec exp_min_y_0 = (min_y - 0) * wlen_cof.y;
  // cout<<"cal 2-pin gradient exp max and min"<<exp_max_x<<"
  // "<<exp_min_x<<endl;
  FPOS e1, e2;
  POS flg1, flg2;
  // cout<<"cal_2pin_wa debug 1"<<endl;
  if(exp_max_x > NEG_MAX_EXP) {
    e1.x = get_exp(exp_max_x);
    sum_num1.x += dist.x * e1.x;
    sum_denom1.x += e1.x;
    flg1.x = 1;
  }
  else {
    flg1.x = 0;
  }
  if(exp_max_x_0 > NEG_MAX_EXP) {
    sum_denom1.x += get_exp(exp_max_x_0);
  }
  // cout<<"cal 2-pin e1x"<<e1.x<<endl;
  // cout<<"cal_2pin_wa debug 2"<<endl;
  if(exp_min_x > NEG_MAX_EXP) {
    e2.x = get_exp(exp_min_x);
    sum_num2.x += dist.x * e2.x;
    sum_denom2.x += e2.x;
    flg2.x = 1;
  }
  else {
    flg2.x = 0;
  }

  if(exp_min_x_0 > NEG_MAX_EXP) {
    sum_denom2.x += get_exp(exp_min_x_0);
  }
  // cout<<"cal 2-pin e2x"<<e2.x<<endl;
  // cout<<"cal_2pin_wa debug 3"<<endl;
  if(exp_max_y > NEG_MAX_EXP) {
    e1.y = get_exp(exp_max_y);
    sum_num1.y += dist.y * e1.y;
    sum_denom1.y += e1.y;
    flg1.y = 1;
  }
  else {
    flg1.y = 0;
  }

  if(exp_max_y_0 > NEG_MAX_EXP) {
    sum_denom1.y += get_exp(exp_max_y_0);
    ;
  }
  // cout<<"cal_2pin_wa debug 4"<<endl;
  if(exp_min_y > NEG_MAX_EXP) {
    e2.y = get_exp(exp_min_y);
    sum_num2.y += dist.y * e2.y;
    sum_denom2.y += e2.y;
    flg2.y = 1;
  }
  else {
    flg2.y = 0;
  }

  if(exp_min_y_0 > NEG_MAX_EXP) {
    sum_denom2.y += get_exp(exp_min_y_0);
  }

  // cout<<"cal_2pin_wa debug 5"<<endl;
  FPOS grad_sum_denom1, grad_sum_denom2, grad_sum_num1, grad_sum_num2;
  FPOS grad1, grad2;
  if(flg1.x) {
    grad_sum_denom1.x = wlen_cof.x * e1.x;
    grad_sum_num1.x = e1.x + dist.x * grad_sum_denom1.x;
    grad1.x =
        (grad_sum_num1.x * sum_denom1.x - grad_sum_denom1.x * sum_num1.x) /
        (sum_denom1.x * sum_denom1.x);
  }
  // cout<<"cal 2-pin grad1x"<<grad1.x<<endl;
  // cout<<"cal 2-pin wlencof"<<wlen_cof.x<<endl;
  // cout<<"cal 2-pin dist"<<dist<<endl;
  if(flg2.x) {
    grad_sum_denom2.x = wlen_cof.x * e2.x;
    grad_sum_num2.x = e2.x - dist.x * grad_sum_denom2.x;
    grad2.x = (grad_sum_num2.x * sum_denom2.x +
               grad_sum_denom2.x *
                   sum_num2.x) /  //((e2-200*0.01e2)*e2 + 0.01e2*200e2)/(e2^2)
              (sum_denom2.x * sum_denom2.x);
  }
  // cout<<"cal 2-pin grad2x"<<grad2.x<<endl;
  // cout<<"cal_2pin_wa debug 6"<<endl;
  if(flg1.y) {
    grad_sum_denom1.y = wlen_cof.y * e1.y;
    grad_sum_num1.y = e1.y + dist.y * grad_sum_denom1.y;
    grad1.y =
        (grad_sum_num1.y * sum_denom1.y - grad_sum_denom1.y * sum_num1.y) /
        (sum_denom1.y * sum_denom1.y);
  }

  // cout<<"cal_2pin_wa debug 7"<<endl;
  if(flg2.y) {
    grad_sum_denom2.y = wlen_cof.y * e2.y;
    grad_sum_num2.y = e2.y - dist.y * grad_sum_denom2.y;
    grad2.y =
        (grad_sum_num2.y * sum_denom2.y + grad_sum_denom2.y * sum_num2.y) /
        (sum_denom2.y * sum_denom2.y);
  }
  // cout<<"cal_2pin_wa debug 8"<<endl;
  grad.x = (prec)grad1.x - grad2.x;
  grad.y = (prec)grad1.y - grad2.y;
  // cout<<"cal_2pin_wa debug 9 "<<grad.x<<" "<< grad.y<<endl;
}
void net_update(FPOS *st) {
  switch(WLEN_MODEL) {
    case LSE:
      return net_update_lse(st);
      break;
    case WA:
      if(fastWL) {
        // cout<<"fast wl net update"<<endl;
        return net_update_wa_fast(st);
      }
      else {
        // cout<<"slow wl net update"<<endl;
        return net_update_wa(st);
      }

      break;
  }
}

void net_update_lse(FPOS *st) {
  int i = 0, j = 0;
  CELL *cell = NULL;
  NET *net = NULL;
  PIN *pin = NULL;
  MODULE *curModule = NULL;
  FPOS fp, pof, center;
  FPOS sum_denom1, sum_denom2;
  prec exp_val = 0;
  prec min_x = 0, min_y = 0;
  prec max_x = 0, max_y = 0;
  prec exp_min_x = 0, exp_min_y = 0;
  prec exp_max_x = 0, exp_max_y = 0;

  for(i = 0; i < gcell_cnt; i++) {
    cell = &gcell_st[i];

    cell->center = st[i];

    cell->den_pmin.x = cell->center.x - cell->half_den_size.x;
    cell->den_pmin.y = cell->center.y - cell->half_den_size.y;
    // cell->den_pmin.z = cell->center.z - cell->half_den_size.z;

    cell->den_pmax.x = cell->center.x + cell->half_den_size.x;
    cell->den_pmax.y = cell->center.y + cell->half_den_size.y;
    // cell->den_pmax.z = cell->center.z + cell->half_den_size.z;
  }

  for(i = 0; i < netCNT; i++) {
    net = &netInstance[i];

    net->min_x = net->terminalMin.x;
    net->min_y = net->terminalMin.y;
    // net->min_z = net->terminalMin.z;

    net->max_x = net->terminalMax.x;
    net->max_y = net->terminalMax.y;
    // net->max_z = net->terminalMax.z;

    for(j = 0; j < net->pinCNTinObject; j++) {
      pin = net->pin[j];

      if(!pin->term) {
        curModule = &moduleInstance[pin->moduleID];
        pof = curModule->pof[pin->pinIDinModule];
        center = st[pin->moduleID];
        fp.x = center.x + pof.x;
        fp.y = center.y + pof.y;
        // fp.z = center.z + pof.z ;
        pin->fp = fp;

        net->min_x = min(net->min_x, fp.x);
        net->min_y = min(net->min_y, fp.y);
        // net->min_z = min ( net->min_z , fp.z ) ;

        net->max_x = max(net->max_x, fp.x);
        net->max_y = max(net->max_y, fp.y);
        // net->max_z = max ( net->max_z , fp.z ) ;
      }
      else {
        continue;
      }
    }

    min_x = net->min_x;
    min_y = net->min_y;
    // min_z = net->min_z ;

    max_x = net->max_x;
    max_y = net->max_y;
    // max_z = net->max_z;

    sum_denom1.x = sum_denom1.y = 0;
    sum_denom2.x = sum_denom2.y = 0;

    for(j = 0; j < net->pinCNTinObject; j++) {
      pin = net->pin[j];

      if(!pin->term) {
#ifdef CELL_CENTER_WLEN_GRAD
        fp = st[pin->moduleID];
#else
        fp = pin->fp;
#endif
      }
      else {
        fp = pin->fp;
      }

      exp_max_x = (fp.x - max_x) * wlen_cof.x;

      if(fabs(exp_max_x) < MAX_EXP) {
        exp_val = get_exp(exp_max_x);
        sum_denom1.x += exp_val;
        pin->flg1.x = 1;
        pin->e1.x = exp_val;
      }
      else {
        exp_val = 0;
        pin->flg1.x = 0;
      }

      exp_min_x = (min_x - fp.x) * wlen_cof.x;

      if(fabs(exp_min_x) < MAX_EXP) {
        exp_val = get_exp(exp_min_x);
        sum_denom2.x += exp_val;
        pin->flg2.x = 1;
        pin->e2.x = exp_val;
      }
      else {
        pin->flg2.x = 0;
      }

      exp_max_y = (fp.y - max_y) * wlen_cof.y;

      if(fabs(exp_max_y) < MAX_EXP) {
        exp_val = get_exp(exp_max_y);
        sum_denom1.y += exp_val;
        pin->flg1.y = 1;
        pin->e1.y = exp_val;
      }
      else {
        pin->flg1.y = 0;
      }

      exp_min_y = (min_y - fp.y) * wlen_cof.y;

      if(fabs(exp_min_y) < MAX_EXP) {
        exp_val = get_exp(exp_min_y);
        sum_denom2.y += exp_val;
        pin->flg2.y = 1;
        pin->e2.y = exp_val;
      }
      else {
        pin->flg2.y = 0;
      }
    }

    net->sum_denom1 = sum_denom1;
    net->sum_denom2 = sum_denom2;
  }
}

prec net_update_hpwl_mac(void) {
  int i = 0, j = 0;
  prec hpwl = 0;
  NET *net = NULL;
  PIN *pin = NULL;
  MODULE *curModule = NULL;
  FPOS fp, pof, p0;

  total_hpwl.x = total_hpwl.y = 0;

  for(i = 0; i < netCNT; i++) {
    net = &netInstance[i];

    net->min_x = net->terminalMin.x;
    net->min_y = net->terminalMin.y;

    net->max_x = net->terminalMax.x;
    net->max_y = net->terminalMax.y;

    for(j = 0; j < net->pinCNTinObject; j++) {
      pin = net->pin[j];
      if(!pin->term) {
        curModule = &moduleInstance[pin->moduleID];
        // cell = & gcell_st [pin->moduleID];
        p0 = curModule->center;
        pof = curModule->pof[pin->pinIDinModule];
        fp.x = p0.x + pof.x;
        fp.y = p0.y + pof.y;
        pin->fp = fp;

        net->min_x = min(net->min_x, fp.x);
        net->min_y = min(net->min_y, fp.y);

        net->max_x = max(net->max_x, fp.x);
        net->max_y = max(net->max_y, fp.y);
      }
    }

    if(net->pinCNTinObject <= 1)
      continue;

    total_hpwl.x += (net->max_x - net->min_x);
    total_hpwl.y += (net->max_y - net->min_y);
  }

  /* total_hpwl_xy = total_hpwl.x + total_hpwl.y; */
  //  total_hpwl_xyz = total_hpwl.x + total_hpwl.y + total_hpwl.z;

  hpwl = total_hpwl.x * dp_wlen_weight.x + total_hpwl.y * dp_wlen_weight.y;

  return hpwl;
}

// WA
//
void net_update_wa(FPOS *st) {
  int i = 0;

  bool timeon = false;
  double time = 0.0f;
  if(timeon)
    time_start(&time);

  omp_set_num_threads(numThread);
  // cout<<"num thread := "<<numThread<<endl;
#pragma omp parallel default(none) shared(gcell_cnt, gcell_st, st) private(i)
  {
//        CELL* cell = NULL;
#pragma omp for
    for(i = 0; i < gcell_cnt; i++) {
      CELL *cell = &gcell_st[i];
      cell->center = st[i];
      cell->den_pmin.x = cell->center.x - cell->half_den_size.x;
      cell->den_pmin.y = cell->center.y - cell->half_den_size.y;
      cell->den_pmax.x = cell->center.x + cell->half_den_size.x;
      cell->den_pmax.y = cell->center.y + cell->half_den_size.y;
    }
  }
  if(timeon) {
    time_end(&time);
    cout << "parallelTime : " << time << endl;
  }
  //  wlen_cof.Dump("current_wlen_cof");
  //  cout << "NEG_MAX_EXP: " << NEG_MAX_EXP << endl;

#pragma omp parallel default(none) shared( \
    netInstance, moduleInstance, st, netCNT, NEG_MAX_EXP, wlen_cof) private(i)
  {
    //        NET *net = NULL;
    //        PIN *pin = NULL;
    //        MODULE *curModule = NULL;
    //        prec exp_min_x = 0, exp_min_y = 0;
    //        prec exp_max_x = 0, exp_max_y = 0;
    //        prec min_x = 0, min_y = 0;
    //        prec max_x = 0, max_y = 0;

    //        FPOS fp, pof, center;
    //        FPOS sum_num1, sum_num2, sum_denom1, sum_denom2;

#pragma omp for
    for(i = 0; i < netCNT; i++) {
      bool timeon = true;
      double time = 0.0f;
      if(timeon) {
        // cout << "GetCost Grad: " << time << endl;
        time_start(&time);
      }
      NET *net = &netInstance[i];
      net->min_x = net->terminalMin.x;
      net->min_y = net->terminalMin.y;
      net->max_x = net->terminalMax.x;
      net->max_y = net->terminalMax.y;

      //        cout << "size: " << sizeof(PIN) << endl;
      //        cout << net->pinCNTinObject << endl;
      for(int j = 0; j < net->pinCNTinObject; j++) {
        PIN *pin = net->pin[j];

        //            cout << j << " " << pin << endl;

        if(!pin->term) {
          MODULE *curModule = &moduleInstance[pin->moduleID];
          FPOS pof = curModule->pof[pin->pinIDinModule];
          FPOS center = st[pin->moduleID];
          FPOS fp;
          fp.x = center.x + pof.x;
          fp.y = center.y + pof.y;
          pin->fp = fp;

          net->min_x = min(net->min_x, fp.x);
          net->min_y = min(net->min_y, fp.y);
          net->max_x = max(net->max_x, fp.x);
          net->max_y = max(net->max_y, fp.y);
        }
        else {
          continue;
        }
      }
      //        if( i >=2 )
      //        exit(1);

      prec min_x = net->min_x;
      prec min_y = net->min_y;
      prec max_x = net->max_x;
      prec max_y = net->max_y;

      FPOS sum_num1, sum_num2;
      FPOS sum_denom1, sum_denom2;

      // UPDATE
      // pin->e1 (MAX)
      // net->sum_num1 (MAX)
      // net->sum_denom1 (MAX)
      // pin->flg1 (MAX)
      //
      // pin->e2 (MIN)
      // net->sum_num2 (MIN)
      // net->sum_denom2 (MIN)
      // pin->flg2 (MIN)
      //
      //
      // Note that NEG_MAX_EXP is -300
      // The meaning NEG_MAX_EXP is. not to have weird out of range values
      // in floating vars.
      //
      // we know that wlen_cof is 1/ gamma.
      // See main.cpp wcof00 and wlen.cpp: wcof_init.
      //

      for(int j = 0; j < net->pinCNTinObject; j++) {
        PIN *pin = net->pin[j];
        FPOS fp = pin->fp;
        prec exp_max_x = (fp.x - max_x) * wlen_cof.x;
        prec exp_min_x = (min_x - fp.x) * wlen_cof.x;
        prec exp_max_y = (fp.y - max_y) * wlen_cof.y;
        prec exp_min_y = (min_y - fp.y) * wlen_cof.y;

        if(exp_max_x > NEG_MAX_EXP) {
          pin->e1.x = get_exp(exp_max_x);
          sum_num1.x += fp.x * pin->e1.x;
          sum_denom1.x += pin->e1.x;
          pin->flg1.x = 1;
        }
        else {
          pin->flg1.x = 0;
        }

        if(exp_min_x > NEG_MAX_EXP) {
          pin->e2.x = get_exp(exp_min_x);
          sum_num2.x += fp.x * pin->e2.x;
          sum_denom2.x += pin->e2.x;
          pin->flg2.x = 1;
        }
        else {
          pin->flg2.x = 0;
        }

        if(exp_max_y > NEG_MAX_EXP) {
          pin->e1.y = get_exp(exp_max_y);
          sum_num1.y += fp.y * pin->e1.y;
          sum_denom1.y += pin->e1.y;
          pin->flg1.y = 1;
        }
        else {
          pin->flg1.y = 0;
        }

        if(exp_min_y > NEG_MAX_EXP) {
          pin->e2.y = get_exp(exp_min_y);
          sum_num2.y += fp.y * pin->e2.y;
          sum_denom2.y += pin->e2.y;
          pin->flg2.y = 1;
        }
        else {
          pin->flg2.y = 0;
        }
      }

      net->sum_num1 = sum_num1;
      net->sum_num2 = sum_num2;
      net->sum_denom1 = sum_denom1;
      net->sum_denom2 = sum_denom2;
      if(timeon && net->pinCNTinObject == 2) {
        time_end(&time);
        // cout << "GetCost Grad: " << time << endl;
        // grad_update_runtime += time;
        // lc_update_runtime += time;

        netupdate_runtime_pinnet += time;
        time_start(&time);
      }
      if(timeon && net->pinCNTinObject == 3) {
        time_end(&time);
        // cout << "GetCost Grad: " << time << endl;
        // grad_update_runtime += time;
        // lc_update_runtime += time;

        netupdate_runtime_3pinnet += time;
        time_start(&time);
      }
    }
  }
}

void net_update_wa_fast(FPOS *st) {
  int i = 0;

  bool timeon = false;
  double time = 0.0f;
  if(timeon)
    time_start(&time);

  omp_set_num_threads(numThread);
#pragma omp parallel default(none) shared(gcell_cnt, gcell_st, st) private(i)
  {
//        CELL* cell = NULL;
#pragma omp for
    for(i = 0; i < gcell_cnt; i++) {
      CELL *cell = &gcell_st[i];
      cell->center = st[i];
      cell->den_pmin.x = cell->center.x - cell->half_den_size.x;
      cell->den_pmin.y = cell->center.y - cell->half_den_size.y;
      cell->den_pmax.x = cell->center.x + cell->half_den_size.x;
      cell->den_pmax.y = cell->center.y + cell->half_den_size.y;
    }
  }
  if(timeon) {
    time_end(&time);
    cout << "parallelTime : " << time << endl;
  }
  //  wlen_cof.Dump("current_wlen_cof");
  //  cout << "NEG_MAX_EXP: " << NEG_MAX_EXP << endl;

#pragma omp parallel default(none) shared( \
    netInstance, moduleInstance, st, netCNT, NEG_MAX_EXP, wlen_cof) private(i)
  {
    //        NET *net = NULL;
    //        PIN *pin = NULL;
    //        MODULE *curModule = NULL;
    //        prec exp_min_x = 0, exp_min_y = 0;
    //        prec exp_max_x = 0, exp_max_y = 0;
    //        prec min_x = 0, min_y = 0;
    //        prec max_x = 0, max_y = 0;

    //        FPOS fp, pof, center;
    //        FPOS sum_num1, sum_num2, sum_denom1, sum_denom2;
    FPOS MinAdmitLength;
    MinAdmitLength.x = 30 * wlen_cof_inv.x;
    MinAdmitLength.y = 30 * wlen_cof_inv.y;
    FPOS AdmitInterval;
    AdmitInterval.x = MinAdmitLength.x * 0.05;
    AdmitInterval.y = MinAdmitLength.y * 0.05;
    // cout << "MinAdmitLength = " << MinAdmitLength.x << ", " <<
    // MinAdmitLength.y
    //      << endl;
    // cout << "AdmitInterval = " << AdmitInterval.x << ", " << AdmitInterval.y
    //      << endl;
#pragma omp for
    for(i = 0; i < netCNT; i++) {
      NET *net = &netInstance[i];

      net->min_x = net->terminalMin.x;
      net->min_y = net->terminalMin.y;
      net->max_x = net->terminalMax.x;
      net->max_y = net->terminalMax.y;
      // if(net->pinCNTinObject < 4) {
      if(1) {
        for(int j = 0; j < net->pinCNTinObject; j++) {
          PIN *pin = net->pin[j];

          //            cout << j << " " << pin << endl;

          if(!pin->term) {
            MODULE *curModule = &moduleInstance[pin->moduleID];
            FPOS pof = curModule->pof[pin->pinIDinModule];
            FPOS center = st[pin->moduleID];
            FPOS fp;
            fp.x = center.x + pof.x;
            fp.y = center.y + pof.y;
            pin->fp = fp;

            net->min_x = min(net->min_x, fp.x);
            net->min_y = min(net->min_y, fp.y);
            net->max_x = max(net->max_x, fp.x);
            net->max_y = max(net->max_y, fp.y);
          }
          else {
            continue;
          }
        }
      }
      else {
        net->min2nd_x = net->terminal2ndMin.x;
        net->min2nd_y = net->terminal2ndMin.y;
        net->max2nd_x = net->terminal2ndMax.x;
        net->max2nd_y = net->terminal2ndMax.y;

        for(int j = 0; j < net->pinCNTinObject; j++) {
          PIN *pin = net->pin[j];

          //            cout << j << " " << pin << endl;

          if(!pin->term) {
            MODULE *curModule = &moduleInstance[pin->moduleID];
            FPOS pof = curModule->pof[pin->pinIDinModule];
            FPOS center = st[pin->moduleID];
            FPOS fp;
            fp.x = center.x + pof.x;
            fp.y = center.y + pof.y;
            pin->fp = fp;

            if(fp.x > net->max_x) {
              net->max2nd_x = net->max_x;
              net->max_x = fp.x;
            }
            else if(fp.x > net->max2nd_x) {
              net->max2nd_x = fp.x;
            }
            if(fp.y > net->max_y) {
              net->max2nd_y = net->max_y;
              net->max_y = fp.y;
            }
            else if(fp.y > net->max2nd_y) {
              net->max2nd_y = fp.y;
            }
            // for min
            if(fp.x < net->min_x) {
              net->min2nd_x = net->min_x;
              net->min_x = fp.x;
            }
            else if(fp.x < net->min2nd_x) {
              net->min2nd_x = fp.x;
            }
            if(fp.y < net->min_y) {
              net->min2nd_y = net->min_y;
              net->min_y = fp.y;
            }
            else if(fp.y < net->min2nd_y) {
              net->min2nd_y = fp.y;
            }
          }
          else {
            continue;
          }
        }
        // cout<<"net  min ="<<net->min_x<<", "<<net->min_y<<endl;
        // cout<<"net  max ="<<net->max_x<<", "<<net->max_y<<endl;
        // cout<<"net 2nd min ="<<net->min2nd_x<<", "<<net->min2nd_y<<endl;
        // cout<<"net 2nd max ="<<net->max2nd_x<<", "<<net->max2nd_y<<endl;
        // if((net->max_x - net->min_x) > MinAdmitLength.x) {
        //   cout << "find long net x" << endl;
        // }
        // if((net->max_x - net->min_x) > MinAdmitLength.x &&
        //    (net->min2nd_x - net->min_x) > AdmitInterval.x &&
        //    (net->max_x > net->max2nd_x) > AdmitInterval.x) {
        //   net->simple_multinet.x = 1;
        //   cout << "find simple net x" << endl;
        // }
        // else {
        //   net->simple_multinet.x = 0;
        // }

        // if((net->max_y - net->min_y) > MinAdmitLength.y &&
        //    (net->min2nd_y - net->min_y) > AdmitInterval.y &&
        //    (net->max_y > net->max2nd_y) > AdmitInterval.y) {
        //   net->simple_multinet.y = 1;
        //   cout << "find simple net y" << endl;
        // }
        // else {
        //   net->simple_multinet.y = 0;
        // }
      }

      //        cout << "size: " << sizeof(PIN) << endl;
      //        cout << net->pinCNTinObject << endl;

      //        if( i >=2 )
      //        exit(1);

      prec min_x = net->min_x;
      prec min_y = net->min_y;
      prec max_x = net->max_x;
      prec max_y = net->max_y;

      FPOS sum_num1, sum_num2;
      FPOS sum_denom1, sum_denom2;

      // UPDATE
      // pin->e1 (MAX)
      // net->sum_num1 (MAX)
      // net->sum_denom1 (MAX)
      // pin->flg1 (MAX)
      //
      // pin->e2 (MIN)
      // net->sum_num2 (MIN)
      // net->sum_denom2 (MIN)
      // pin->flg2 (MIN)
      //
      //
      // Note that NEG_MAX_EXP is -300
      // The meaning NEG_MAX_EXP is. not to have weird out of range values
      // in floating vars.
      //
      // we know that wlen_cof is 1/ gamma.
      // See main.cpp wcof00 and wlen.cpp: wcof_init.
      //
      if(net->pinCNTinObject < 4)  // apply fast model in 2 pin net
      {
        continue;
      }
      for(int j = 0; j < net->pinCNTinObject; j++) {
        PIN *pin = net->pin[j];
        FPOS fp = pin->fp;
        prec exp_max_x = (fp.x - max_x) * wlen_cof.x;
        prec exp_min_x = (min_x - fp.x) * wlen_cof.x;
        prec exp_max_y = (fp.y - max_y) * wlen_cof.y;
        prec exp_min_y = (min_y - fp.y) * wlen_cof.y;

        if(exp_max_x > NEG_MAX_EXP) {
          pin->e1.x = get_exp(exp_max_x);
          sum_num1.x += fp.x * pin->e1.x;
          sum_denom1.x += pin->e1.x;
          pin->flg1.x = 1;
        }
        else {
          pin->flg1.x = 0;
        }

        if(exp_min_x > NEG_MAX_EXP) {
          pin->e2.x = get_exp(exp_min_x);
          sum_num2.x += fp.x * pin->e2.x;
          sum_denom2.x += pin->e2.x;
          pin->flg2.x = 1;
        }
        else {
          pin->flg2.x = 0;
        }

        if(exp_max_y > NEG_MAX_EXP) {
          pin->e1.y = get_exp(exp_max_y);
          sum_num1.y += fp.y * pin->e1.y;
          sum_denom1.y += pin->e1.y;
          pin->flg1.y = 1;
        }
        else {
          pin->flg1.y = 0;
        }

        if(exp_min_y > NEG_MAX_EXP) {
          pin->e2.y = get_exp(exp_min_y);
          sum_num2.y += fp.y * pin->e2.y;
          sum_denom2.y += pin->e2.y;
          pin->flg2.y = 1;
        }
        else {
          pin->flg2.y = 0;
        }
      }

      net->sum_num1 = sum_num1;
      net->sum_num2 = sum_num2;
      net->sum_denom1 = sum_denom1;
      net->sum_denom2 = sum_denom2;
    }
  }
}

prec get_mac_hpwl(int idx) {
  MODULE *mac = macro_st[idx];
  PIN *pin = NULL;
  NET *net = NULL;
  int i = 0, j = 0;
  int moduleID = mac->idx;
  CELL *cell = &gcell_st[moduleID];
  FPOS fp;

  mac_hpwl.x = mac_hpwl.y = 0;
  mac_hpwl_xyz = 0;

  for(i = 0; i < cell->pinCNTinObject; i++) {
    pin = cell->pin[i];
    net = &netInstance[pin->netID];

    net->min_x = net->terminalMin.x;
    net->min_y = net->terminalMin.y;
    net->max_x = net->terminalMax.x;
    net->max_y = net->terminalMax.y;

    for(j = 0; j < net->pinCNTinObject; j++) {
      fp = net->pin[j]->fp;
      net->min_x = min(net->min_x, fp.x);
      net->min_y = min(net->min_y, fp.y);
      net->max_x = max(net->max_x, fp.x);
      net->max_y = max(net->max_y, fp.y);
    }

    if(net->pinCNTinObject <= 1)
      continue;

    mac_hpwl.x += (net->max_x - net->min_x);
    mac_hpwl.y += (net->max_y - net->min_y);
  }

  mac_hpwl_xyz = mac_hpwl.x * dp_wlen_weight.x + mac_hpwl.y * dp_wlen_weight.y;

  return mac_hpwl_xyz;
}

// update net->pin2->fp as updated version (modules' center + modules' offset)
void update_pin2(void) {
  //    cout << "called update_pin2" << endl;
  for(int i = 0; i < netCNT; i++) {
    NET *net = &netInstance[i];

    //        cout << net->pinCNTinObject2 << endl;
    for(int j = 0; j < net->pinCNTinObject2; j++) {
      PIN *pin = net->pin2[j];

      // if pin is terminal pin, then skip this procedure.
      if(pin->term) {
        continue;
      }

      MODULE *curModule = &moduleInstance[pin->moduleID];
      pin->fp.SetAdd(curModule->center, curModule->pof[pin->pinIDinModule]);
    }
  }
}

// Update HPWL based on net->pin2 definition.
void UpdateNetMinMaxPin2(void) {
  NET *net = NULL;
  PIN *pin = NULL;
  FPOS fp;

  for(int i = 0; i < netCNT; i++) {
    net = &netInstance[i];

    net->min_x = net->terminalMin.x;
    net->min_y = net->terminalMin.y;

    net->max_x = net->terminalMax.x;
    net->max_y = net->terminalMax.y;

    for(int j = 0; j < net->pinCNTinObject2; j++) {
      pin = net->pin2[j];

      if(pin->term) {
        continue;
      }
      fp = pin->fp;

      net->min_x = min(net->min_x, fp.x);
      net->min_y = min(net->min_y, fp.y);

      net->max_x = max(net->max_x, fp.x);
      net->max_y = max(net->max_y, fp.y);
    }

    if(net->pinCNTinObject2 <= 1)
      continue;
  }
}

// Get HPWL as micron units
pair< double, double > GetUnscaledHpwl() {
  double x = 0.0f, y = 0.0f;

  NET *curNet = NULL;
  for(int i = 0; i < netCNT; i++) {
    curNet = &netInstance[i];
    x += (curNet->max_x - curNet->min_x) * GetUnitX() / GetDefDbu();
    y += (curNet->max_y - curNet->min_y) * GetUnitY() / GetDefDbu();

    if(curNet->max_x - curNet->min_x < 0) {
      cout << "NEGATIVE HPWL ERROR! " << curNet->Name() << " " << curNet->max_x
           << " " << curNet->min_x << endl;
    }
    if(curNet->max_y - curNet->min_y < 0) {
      cout << "NEGATIVE HPWL ERROR! " << curNet->Name() << " " << curNet->max_y
           << " " << curNet->min_y << endl;
    }

    if(x < 0 || y < 0) {
      printf("NEGATIVE HPWL ERROR! \n");
      cout << curNet->Name() << x << " " << y << endl;
      exit(1);
    }
  }
  return make_pair(x, y);
}

void PrintUnscaledHpwl(string mode) {
  pair< double, double > hpwl = GetUnscaledHpwl();
  cout << "===HPWL(MICRON)====================================" << endl;
  cout << "  Mode  : " << mode << endl;
  cout << "  HPWL  : " << std::fixed << std::setprecision(4)
       << hpwl.first + hpwl.second << endl
       << "          x= " << hpwl.first << " y= " << hpwl.second << endl;
  cout << "===================================================" << endl;
}
