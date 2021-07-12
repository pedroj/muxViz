#ifndef GREEDY_H
#define GREEDY_H

#include "MersenneTwister.h"
#include "GreedyBase.h"
#include "Node.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <stack>
#include <map>
#include <algorithm>
using namespace std;


class Greedy : public GreedyBase{
 public:
  Greedy(MTRand *RR,int nnode,Node **node,int nmembers,double Alpha);
  virtual ~Greedy();
  virtual void initiate(void);
  virtual void calibrate(void);
  virtual void tune(void);
  virtual void prepare(bool sort);
  virtual void level(Node ***,bool sort);
  virtual void move(int &moved);
  virtual void determMove(vector<int> &moveTo);
  virtual void eigenvector(void);
  virtual void eigenfactor(void);


  int Nempty;
  vector<int> mod_empty;

  vector<double> mod_enter;
  vector<double> mod_exit;
  vector<double> mod_size;
  vector<int> mod_members;
  vector<map<int,pair<int,double> > > mod_nodeSize;

 protected:

  vector<int> modSnode;
};

#endif
