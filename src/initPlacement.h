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

#ifndef __PL_IP__
#define __PL_IP__
#include <vector>
#include <Eigen/SparseCore>
#include "replace_private.h"
#include <fstream>
#include <random>

using Eigen::VectorXf;
typedef Eigen::SparseMatrix< prec, Eigen::RowMajor > SMatrix;
typedef Eigen::Triplet< prec > T;

void initial_placement();

void CreateSparseMatrix(VectorXf &xcg_x, VectorXf &xcg_b, VectorXf &ycg_x,
                        VectorXf &ycg_b, SMatrix &eMatX, SMatrix &eMatY);

void update_module(VectorXf &xcg_x, VectorXf &ycg_x);
void build_data_struct(bool initCoordi = true);
void update_pin_by_module();


///////////////////////////////////////////
// std::vector<int> _Grid;
// std::vector<MODULE> _Containers;
// std::vector<PIN> _TermPins;
// std::vector<NET> _Nets;

void doGreedyPlace(int scale);
void doBoundingBoxPlaceForBuffer(int scale);
void doSchemPlace_TermContainer(int scale);
void doDeterPlace(int scale);
void setupGreedyPlace(FPOS& gridSize,POS& gridNum,std::vector<int>& Grid);
void setupGrid(FPOS& gridSize,POS& gridNum,std::vector<int>& Grid);
void setupContainers();
void setupTermPins();
void setupNets();
void generateSearchOrder(std::vector<int> &order,int scale);
int findBestBin(MODULE& curContainer,vector<int>& order,int scale,std::vector<int> Grid,FPOS gridSize);
int findBestBin(TERM& curContainer,vector<int>& order,int scale,std::vector<int> Grid,FPOS gridSize);
void findActiveNet(MODULE& curContainer,vector<int>& order,vector<int>& activeNetId,int scale);
void placeContainerToBin(MODULE& curContainer,int binId,FPOS location,vector<int> &Grid);
void placeContainerToBin(TERM& curContainer,int binId,FPOS location,vector<int> &Grid);
prec getHpwl(MODULE& curContainer,std::vector<int> activeNets,FPOS location);
FPOS getHpwl(NET& curNet);
void findCandidate(vector<int>&Grid,FPOS LL,FPOS UR,vector<FPOS> &candidates,vector<int>&candidateBinId,FPOS gridSize);
void findCandidate(MODULE& curModule,vector<int> &order,vector<FPOS> &candidates,int scale);
void findBoundingBox(vector<int> &activeNetId,FPOS& LL,FPOS& UR,MODULE& curContainer,vector<int> Order);
void plot(int scale);
void plotTerm(int scale);
void drawRectangle(FPOS LL,FPOS UR,int colors,string filename, int objectID,string label);

void printDebugInfo(string funcName,int debugIndex,string info);

void fillPinCordList(vector<prec> pinXlist,vector<prec> pinYlist,vector<int> &Nets,MODULE& curModule,vector<int>& order);
void move2nonOverlap(MODULE& curModule,FPOS fp,FPOS direction,vector<int> &order,vector<FPOS> &candidates);
int lookaheadDetermine(MODULE& curModule,MODULE& nextModule,vector<int> order,vector<FPOS>& candidates,int scale);
void placeContainer(MODULE &curContainer, FPOS location);
void placeContainer(TERM &curContainer, FPOS location);

void findIOCenter(FPOS &ioCenter);
void findPEArraySize(FPOS &peArraySize);
void findBufferSize(FPOS &bufferSize);

void roughLegalization(int PENum);
// void placeContainer_Term(TERM &curContainer, FPOS location);
void stitchingLegalization(int PENum);
// FPOS spaceCenter;
#endif
