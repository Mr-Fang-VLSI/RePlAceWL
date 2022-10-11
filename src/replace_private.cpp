#include "replace_private.h"
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>


using std::string;
using std::cout;
using std::endl;

using std::vector;
using std::pair;
using std::tuple;
using std::max;
using std::min;

using std::fixed;
using std::scientific;
using std::map;
using std::to_string;

using std::numeric_limits;
using std::make_pair;

#ifdef USE_GOOGLE_HASH
using google::dense_hash_map;
using google::dense_hash_set;
#else
using std::unordered_map;
using std::unordered_set;
#endif


// Functions for FPOS
FPOS::FPOS() {
  x = y = 0; 
};

FPOS::FPOS(prec xloc, prec yloc) : x(xloc), y(yloc) {};

void FPOS::Set(prec a) {
  x = y = a;
}

void FPOS::Set(FPOS a) {
  x = a.x;
  y = a.y;
}

void FPOS::Set(prec xloc, prec yloc) {
  x = xloc;
  y = yloc;
}

void FPOS::SetZero() {
  x = y = 0;
}

prec FPOS::GetX(){ return x; };
prec FPOS::GetY(){ return y; }; 

void FPOS::Add(FPOS a) {
  x += a.x;
  y += a.y;
}

void FPOS::SetAdd(FPOS a, FPOS b) {
  x = a.x + b.x;
  y = a.y + b.y;
}
void FPOS::Min(FPOS a) {
  x = min(x, a.x);
  y = min(y, a.y);
}
void FPOS::SetMin(FPOS a, FPOS b) {
  x = min(a.x, b.x);
  y = min(a.y, b.y);
}
void FPOS::Max(FPOS a) {
  x = max(x, a.x);
  y = max(y, a.y);
}
void FPOS::SetMax(FPOS a, FPOS b) {
  x = max(a.x, b.x);
  y = max(a.y, b.y);
}

prec FPOS::GetProduct() {
  return x * y;
}
void FPOS::Dump() {
  cout << "(" << x << " " << y << ")" << endl;
}
void FPOS::Dump(string a) {
  cout << a << ": (" << x << " " << y << ")" << endl;
}
void FPOS::Set(POS p) {
  x = p.x;
  y = p.y;
}



POS::POS() {
  x = y = 0;
};

POS::POS(int xloc, int yloc) : x(xloc), y(yloc) {};

void POS::Set(FPOS fp) {
  x = INT_CONVERT(fp.x);
  y = INT_CONVERT(fp.y);
}


void POS::Set(int a) {
  x = y = a;
}

void POS::Set(POS a) {
  x = a.x;
  y = a.y;
}

void POS::Set(int xloc, int yloc){
  x = xloc;
  y = yloc;
}
void POS::SetZero() {
  x = y = 0;
}
void POS::Add(POS a) {
  x += a.x;
  y += a.y;
}
void POS::SetAdd(POS a, POS b) {
  x = a.x + b.x;
  y = a.y + b.y;
}
void POS::Min(POS a) {
  x = min(x, a.x);
  y = min(y, a.y);
}
void POS::SetMin(POS a, POS b) {
  x = min(a.x, b.x);
  y = min(a.y, b.y);
}
void POS::Max(POS a) {
  x = max(x, a.x);
  y = max(y, a.y);
}
void POS::SetMax(POS a, POS b) {
  x = max(a.x, b.x);
  y = max(a.y, b.y);
}
int POS::GetProduct() {
  return x * y ;
}
void POS::SetXProjection(int a, int b) {
  x = (x < a) ? a : (x > b) ? b : x;
}
void POS::SetYProjection(int a, int b) {
  y = (y < a) ? a : (y > b) ? b : y;
}
void POS::SetProjection(POS a, POS b) {
  SetXProjection(a.x, b.x);
  SetYProjection(a.y, b.y);
}
void POS::SetXYProjection(POS a, POS b) {
  SetXProjection(a.x, b.x);
  SetYProjection(a.y, b.y);
}
void POS::Dump() {
  cout << "(" << x << " " << y << ")" << endl;
}

void POS::Dump(std::string a) {
  cout << a << ": (" << x << " " << y << ")" << endl;
}



RECT::RECT() {
  pmin.SetZero();
  pmax.SetZero();
};
void RECT::Dump() {
  pmin.Dump("RECT: pmin");
  pmax.Dump("RECT: pmax");
  cout << endl;
}

PIN::PIN() {
  fp.SetZero();
  e1.SetZero();
  e2.SetZero();
  flg1.SetZero();
  flg2.SetZero();
}



const char* MODULE::Name() { 
  return moduleNameStor[idx].c_str(); 
}

void MODULE::SetContainerType(){
  std::string Name = moduleNameStor[idx];
  std::vector<std::string> clause;
  char breaksign = '_';
  // cout<<"name is "<<Name<<endl;
  
  BreakDownName(Name, breaksign, clause);
  // for(int i = 0;i<clause.size();i++)
  // {
  //   cout<<clause.at(i)<<" ";
  // }
  
  // cout<<Name<<endl;
  std::vector<std::string>::iterator it;
  it = find(clause.begin(),clause.end(),"genper");
  type = "";
  if(it!=clause.end())
  {
    type = "PE";
    containerCORD.x = std::stoi(*(it+1));
    it = find(clause.begin(),clause.end(),"genpec");
    containerCORD.y = std::stoi(*(it+1));
    // cout<<"PE cord = "<<containerCORD.x<<" "<<containerCORD.y<<endl;
    // cout<<Name<<endl;
  }
  else 
  {
    it = find(clause.begin(),clause.end(),"genbc");
    if(it != clause.end())
    {
      type = "buffer";
      containerCORD.x = std::stoi(*(it+1));
    }
    else{
      type = "";
    }
    
    
  }
  
  std::vector<std::string>().swap(clause);
  
}
void MODULE::SetContainerCORD(std::vector<std::string>& clause){
  if(type ==  "PE")
  {
    std::vector<std::string>::iterator it;
    // it = find(clause.begin().)
  }
  else{

  }
}
void BreakDownName(std::string Name,char breaksign,std::vector<std::string>& clause){
  
  int pos = 0;
  for(int i = 0;i < Name.length();i++)
  {
    if(Name.at(i) == breaksign)
    {
      clause.push_back(Name.substr(pos,i-pos));
      pos = (i+1)<Name.length()?(i+1):i;
    }
  }
  clause.push_back(Name.substr(pos,Name.length()-1-pos));
 
}

void testNameBreak(){
  MODULE* module;
  for(int i = 0;i < moduleCNT;i++)
  {
    
    module = &moduleInstance[i];
    
    string type;
    
    module->SetContainerType();
    
  }
}
int MODULE::containerID(int scale){
  int id = 0;
  if(type == "PE")
  {
    id = containerCORD.x*scale+containerCORD.y;
  }
  else if(type == "buffer"){
    id = scale*scale + containerCORD.x;
  }
  return id;
}
void ClusterModuleAndNet(int scale){
  cout<<"start cluster"<<endl;
  MODULE* module;
  int clusterModuleCNT = scale*(scale+1);
  int j = clusterModuleCNT;
  int PEpinCNT= 0;
  int BuffpinCNT = 0;
  for(int i = 0;i < moduleCNT;i++)
  {
    
    module = &moduleInstance[i];
    
    
    
    module->SetContainerType();
    if(module->type == "")
    {
      clusterModuleCNT++;
    }
    if(module->type =="PE"&& module->containerCORD.x == 0&&module->containerCORD.y == 0)
    {
      PEpinCNT +=  module->pinCNTinObject;
      
    }
    if(module->type =="buffer"&& module->containerCORD.x == 0)
    {
      BuffpinCNT +=  module->pinCNTinObject;
      
    }
    
  }
  // cout<<"pinCNT in PE"<<PEpinCNT<<endl;
  // cout<<"clusterCNT "<<clusterModuleCNT<<endl;
  // moduleInstance_origin = (struct MODULE *)malloc(sizeof(struct MODULE) * moduleCNT);
  // moduleInstance_origin = moduleInstance;
  
  // free(moduleInstance);
 
  // moduleInstance_origin = (struct MODULE *)malloc(sizeof(struct MODULE) * clusterModuleCNT);
  cout<<" cluster debug 0"<<endl;
 moduleInstance_origin.resize(clusterModuleCNT);
  int PEcnt = scale*scale;
  for(int i = 0;i < j;i++)
  {
    MODULE *tmpModule = &moduleInstance_origin[i];
    
    
    if(j<PEcnt)
    {
      // tmpModule->pin = (struct PIN **)malloc(sizeof(struct PIN *) * PEpinCNT);
      // tmpModule->pof = (struct FPOS*)malloc(sizeof(struct FPOS) * PEpinCNT);
      tmpModule->pin.resize(PEpinCNT);
      tmpModule->pof.resize(PEpinCNT);
    }
    else{
      // tmpModule->pin = (struct PIN **)malloc(sizeof(struct PIN *) * BuffpinCNT);
      // tmpModule->pof = (struct FPOS*)malloc(sizeof(struct FPOS) * BuffpinCNT);
      tmpModule->pin.resize(BuffpinCNT);
      tmpModule->pof.resize(BuffpinCNT);
    }
    
  }
  cout<<"cluster debug 4"<<endl;
  MODULE* curModule;
  for(int i = 0;i < moduleCNT;i++)
  {
    cout<<"cluster debug 4.4"<<endl;
    module = &moduleInstance[i];
    cout<<"cluster debug 4.5 type = "<<module->type<<endl;
    if(module->type == "")
    {
      cout<<"cluster debug 4.7"<<endl;
      moduleInstance_origin[j].cpy(module);
      cout<<"cluster debug 5"<<endl;
      j++;
    }
    else
    {
      int PEid;
      cout<<"cluster debug 6"<<endl;
      curModule = &moduleInstance_origin[PEid];
      PEid = module->containerID(scale);
      curModule->area += module->area;
      for(int k = 0;k <  module->pinCNTinObject;k++)
      {
        curModule->pin[curModule->pinCNTinObject] = module->pin[k];
        curModule->pinCNTinObject++;
      }
    }
    
  }
  cout<<"cluster debug 11"<<endl;

}

void MODULE::cpy(MODULE* COPY)
{
  area = COPY->area;
  pinCNTinObject = COPY->pinCNTinObject;
  netCNTinObject = 0;
  flg = COPY->flg;
  tier = COPY->tier;
  mac_idx = COPY->mac_idx;
  ovlp_flg = COPY->ovlp_flg;
  pmin = COPY->pmin;
  pmax= COPY->pmax;
  size= COPY->size;
  half_size= COPY->pmin;
  center= COPY->center;
  pmin_lg= COPY->pmin_lg;
  pmax_lg= COPY->pmax_lg;
  cout<<"cluster debug 8"<<endl;
  // pof = (struct FPOS*)malloc(sizeof(struct FPOS) * pinCNTinObject);
  // pin= (struct PIN**)malloc(sizeof(struct PIN*) * pinCNTinObject);
  pof.resize(pinCNTinObject);
  pin.resize(pinCNTinObject);
  cout<<"cluster debug 9"<<endl;
  // type = COPY->type;
  cout<<"cluster debug 9.5"<<endl;
  containerCORD.x = COPY->containerCORD.x;
  containerCORD.y = COPY->containerCORD.y;
  
  
  for(int i = 0;i < pinCNTinObject;i++)
  {
    pin[i] = COPY->pin[i];
    pof[i].x = COPY->pof[i].x;
    pof[i].y = COPY->pof[i].y;
  }
  cout<<"cluster debug 10"<<endl;
}
MODULE::MODULE()
  : pof(0),
  pin(0),
  area(0.0f),
  idx(0),
  netCNTinObject(0),
  pinCNTinObject(0),
  flg(0),
  tier(0),
  mac_idx(0),
  ovlp_flg(0) {
    pmin.SetZero();
    pmax.SetZero();
    size.SetZero();
    half_size.SetZero();
    center.SetZero();
    pmin_lg.SetZero();
    pmax_lg.SetZero();
    containerCORD.SetZero();
    type = "";
  }

void MODULE::Dump(string a) {
  cout << a << endl;
  cout << "tier: " << tier << endl;
  cout << "mac_idx: " << mac_idx << endl;
  cout << "ovlp_flg: " << ovlp_flg << endl;
  pmin.Dump("pmin");
  pmax.Dump("pmax");
  size.Dump("size");
  half_size.Dump("half_size");
  center.Dump("center");
  cout << "area: " << area << endl;
  cout << "netCNTinObject: " << netCNTinObject << endl;
  cout << "pinCNTinObject: " << pinCNTinObject << endl;
  pmin_lg.Dump("pmin_lg");
  pmax_lg.Dump("pmax_lg");
  cout << endl;
}

const char* TERM::Name(){ 
  return terminalNameStor[idx].c_str(); 
}


TERM::TERM()
  : area(0.0f),
  pof(0),
  pin(0),
  idx(0),
  netCNTinObject(0),
  pinCNTinObject(0),
  IO(0),
  isTerminalNI(0),
  PL_area(0.0f) {
    pmin.SetZero();
    pmax.SetZero();
    size.SetZero();
    center.SetZero();
  }

void TERM::Dump() {
  printf("terminal[%d]: name: %s \n", idx, Name());
  fflush(stdout);
  cout << "isTerminalNI: " << (isTerminalNI ? "YES" : "NO") << endl;
  cout << "IO: " << ((IO == 0) ? "Input" : "Output") << endl;
  //        cout << "tier: " << tier << endl;
  pmin.Dump("pmin");
  pmax.Dump("pmax");
  cout << "area: " << area << endl;
  size.Dump("size");
  center.Dump("center");
  cout << "netCNTinObject: " << netCNTinObject << endl;
  cout << "pinCNTinObject: " << pinCNTinObject << endl;
  cout << "PL_area: " << PL_area << endl;
  cout << endl;
}
  
const char* CELL::Name() { 
  return cellNameStor[idx].c_str(); 
}

void CELL::Dump() {
  printf("gcell[%d]:name: %s\n", idx, Name());
  fflush(stdout);
  pmin.Dump("pmin");
  pmax.Dump("pmax");
  den_pmin.Dump("den_pmin");
  den_pmax.Dump("den_pmax");
  cout << "grad: x: " << x_grad << " y: " << y_grad << endl;
  cout << "den_scal: " << den_scal << endl << endl;
}
  
SHAPE::SHAPE(std::string _name, std::string _instName, int _idx, prec _llx, prec _lly,
        prec _width, prec _height)
      : name(_name),
        instName(_instName),
        idx(_idx),
        llx(_llx),
        lly(_lly),
        width(_width),
        height(_height){};

  void SHAPE::Dump() {
    printf(
        "shape[%d]: name: %s, instName: %s, llx: %lf, lly: %lf, width: %lf, "
        "height: %lf\n",
        idx, name.c_str(), instName.c_str(), llx, lly, width, height);
    fflush(stdout);
  }

UFPin::UFPin() : parent(0), rank(0), modu(0) {};
UFPin::UFPin(int moduleID) {
  // id = 0;
  parent = moduleID;
  rank = 0;
  modu = moduleID;
}

TwoPinNets::TwoPinNets() : selected(false), 
  start_modu(0), end_modu(0), rect_dist(0), 
  i(0), j(0) {};

TwoPinNets::TwoPinNets(int start_modu, int end_modu, 
      prec rect_dist, int i, int j) {
    selected = false;
    this->start_modu = start_modu;
    this->end_modu = end_modu;
    this->rect_dist = rect_dist;
    this->i = i;
    this->j = j;
}

ROUTRACK::ROUTRACK() : layer(INT_MAX), netIdx(INT_MAX) { 
  from.SetZero(); 
  to.SetZero(); 
};

ROUTRACK::ROUTRACK(struct FPOS _from, struct FPOS _to, 
    int _layer, int _netIdx) {
  from.Set(_from);
  to.Set(_to);
  layer = _layer;
  netIdx = _netIdx;
};

void ROUTRACK::Dump() {
  from.Dump("from"); 
  to.Dump("to"); 
  cout << "layer: " << layer << endl;
  cout << "netIdx: " << netIdx << endl << endl;
}
  
const char* NET::Name() { 
  return netNameStor[idx].c_str(); 
}

NET::NET() : min_x(PREC_MAX), min_y(PREC_MAX), 
  max_x(PREC_MIN), max_y(PREC_MIN),
  pin(0), pin2(0), hpwl_x(PREC_MIN), hpwl_y(PREC_MIN),  
  outPinIdx(INT_MAX), pinCNTinObject(INT_MAX), pinCNTinObject2(INT_MAX),
  pinCNTinObject_tier(INT_MAX), idx(INT_MAX), mod_idx(INT_MAX), 
  timingWeight(1.0f), 
  customWeight(1.0f),
  wl_rsmt(0.0f) { 
    sum_num1.SetZero();
    sum_num2.SetZero();
    sum_denom1.SetZero();
    sum_denom2.SetZero();
    terminalMin.SetZero();
    terminalMax.SetZero();
  }

ROW::ROW() : site_wid(0),
  site_spa(0),
  ori(""),
  isXSymmetry(false),
  isYSymmetry(false),
  isR90Symmetry(false),
  x_cnt(0) {
    pmin.SetZero();
    pmax.SetZero();
    size.SetZero();
  }

void ROW::Dump(std::string a) {
    cout << a << endl;
    cout << "site_wid: " << site_wid << endl;
    cout << "site_spa: " << site_spa << endl;
    cout << "ori: " << ori << endl;
    cout << "x_cnt: " << x_cnt << endl;
    pmin.Dump("pmin");
    pmax.Dump("pmax");
    size.Dump("size");
    cout << endl;
  }


void PLACE::Dump(std::string a) {
    cout << a << endl;
    org.Dump("origin");
    end.Dump("end");
    center.Dump("center");
    stp.Dump("stp");
    cnt.Dump("cnt");
    cout << endl;
  }
