#include "Greedy.h"
#define plogp( x ) ( (x) > 0.0 ? (x)*log(x) : 0.0 )

Greedy::~Greedy(){
  vector<int>().swap(modSnode);
}

Greedy::Greedy(MTRand *RR,int nnode,Node **ah,int nmember,double Alpha){

  R = RR;
  Nnode = nnode;
  NphysicalNode = nmember;
  node = ah;
  Nmod = Nnode;

  alpha = Alpha; // teleportation probability
  beta = 1.0-alpha; // probability to take normal step

  Ndanglings = 0;

}

void Greedy::move(int &moved){

  moved = 0;

  // Generate random enumeration of nodes
  vector<int> randomOrder(Nnode);
  for(int i=0;i<Nnode;i++)
    randomOrder[i] = i;
  for(int i=0;i<Nnode-1;i++){
    int randPos = i + R->randInt(Nnode-i-1);
    int tmp = randomOrder[i];
    randomOrder[i] = randomOrder[randPos];
    randomOrder[randPos] = tmp;
  }

  unsigned int offset = 1;
  vector<unsigned int> redirect(Nnode,0);
  vector<pair<int,pair<double,double> > > flowNtoM(Nnode);
  vector<pair<double,double> > overlapNtoM(Nnode);
  pair<double,double> overlapNtoOldM;

  for(int k=0;k<Nnode;k++){

    // Pick nodes in random order
    int flip = randomOrder[k];
    int oldM = node[flip]->index;

    // Reset offset when int overflows
    if(offset > INT_MAX){
      for(int j=0;j<Nnode;j++)
        redirect[j] = 0;
      offset = 1;
    }

    // Size of vector with module links
    int NmodLinks = 0;

    // For all outLinks
    int NoutLinks = node[flip]->outLinks.size();
    if(NoutLinks == 0){ //dangling node, add node to calculate flow below
      redirect[oldM] = offset + NmodLinks;
      flowNtoM[NmodLinks].first = oldM;
      flowNtoM[NmodLinks].second.first = 0.0;
      flowNtoM[NmodLinks].second.second = 0.0;
      overlapNtoM[NmodLinks].first = 0.0;
      overlapNtoM[NmodLinks].second = 0.0;
      NmodLinks++;
    }
    else{
      for(int j=0; j<NoutLinks; j++){
        int nb_M = node[node[flip]->outLinks[j].first]->index;
        double nb_flow = node[flip]->outLinks[j].second;
        if(redirect[nb_M] >= offset){
          flowNtoM[redirect[nb_M] - offset].second.first += nb_flow;
        }
        else{
          redirect[nb_M] = offset + NmodLinks;
          flowNtoM[NmodLinks].first = nb_M;
          flowNtoM[NmodLinks].second.first = nb_flow;
          flowNtoM[NmodLinks].second.second = 0.0;
          overlapNtoM[NmodLinks].first = 0.0;
          overlapNtoM[NmodLinks].second = 0.0;
          NmodLinks++;
        }
      }
    }

    // For all inLinks
    int NinLinks = node[flip]->inLinks.size();
    for(int j=0; j<NinLinks; j++){
      int nb_M = node[node[flip]->inLinks[j].first]->index;
      double nb_flow = node[flip]->inLinks[j].second;

      if(redirect[nb_M] >= offset){
        flowNtoM[redirect[nb_M] - offset].second.second += nb_flow;
      }
      else{
        redirect[nb_M] = offset + NmodLinks;
        flowNtoM[NmodLinks].first = nb_M;
        flowNtoM[NmodLinks].second.first = 0.0;
        flowNtoM[NmodLinks].second.second = nb_flow;
        overlapNtoM[NmodLinks].first = 0.0;
        overlapNtoM[NmodLinks].second = 0.0;
        NmodLinks++;
      }
    }

    // Calculate flow to/from own module (default value if no link to own module)
    double outFlowOldM = 0.0;
    double inFlowOldM = 0.0;
    if(redirect[oldM] >= offset){
      outFlowOldM = flowNtoM[redirect[oldM] - offset].second.first;
      inFlowOldM = flowNtoM[redirect[oldM] - offset].second.second;
    }

    // Option to move to empty module (if node not already alone)
    if(mod_members[oldM] > static_cast<int>(node[flip]->members.size())){
      if(Nempty > 0){
        flowNtoM[NmodLinks].first = mod_empty[Nempty-1];
        flowNtoM[NmodLinks].second.first = 0.0;
        flowNtoM[NmodLinks].second.second = 0.0;
        overlapNtoM[NmodLinks].first = 0.0;
        overlapNtoM[NmodLinks].second = 0.0;
        NmodLinks++;
      }
    }

    overlapNtoOldM.first = 0.0;
    overlapNtoOldM.second = 0.0;
    // For all multiple assigned nodes
    int flipNphysicalNode = node[flip]->physicalNodes.size();
    for(int j=0;j<flipNphysicalNode;j++){
      int physicalNode = node[flip]->physicalNodes[j].first;
      double physicalNodeSize = node[flip]->physicalNodes[j].second;
      for(map<int,pair<int,double> >::iterator overlap_it = mod_nodeSize[physicalNode].begin(); overlap_it != mod_nodeSize[physicalNode].end(); ++overlap_it){
        int nb_M = overlap_it->first;
            if(nb_M == oldM){ // From where the multiple assigned node is moved
              double oldP = overlap_it->second.second;
              double newP = overlap_it->second.second - physicalNodeSize;
              overlapNtoOldM.first += plogp(newP) - plogp(oldP);
              overlapNtoOldM.second += plogp(physicalNodeSize);
            }
            else{ // To where the multiple assigned node is moved
              double oldP = overlap_it->second.second;
              double newP = overlap_it->second.second + physicalNodeSize;
              if(redirect[nb_M] >= offset){
                overlapNtoM[redirect[nb_M] - offset].first += plogp(newP) - plogp(oldP);
                overlapNtoM[redirect[nb_M] - offset].second += plogp(physicalNodeSize);
                // cout << "(" << flip << ">" << nb_M << " " << (plogp(newP) - plogp(oldP))/log(2.0) << ") ";
              }
              else{
                redirect[nb_M] = offset + NmodLinks;
                flowNtoM[NmodLinks].first = nb_M;
                flowNtoM[NmodLinks].second.first = 0.0;
                flowNtoM[NmodLinks].second.second = 0.0;
                overlapNtoM[NmodLinks].first = plogp(newP) - plogp(oldP);
                overlapNtoM[NmodLinks].second = plogp(physicalNodeSize);
                NmodLinks++;
              }
            }
          }
        }

    // Randomize link order for optimized search
        for(int j=0;j<NmodLinks-1;j++){
          int randPos = j + R->randInt(NmodLinks-j-1);
          pair<int,pair<double,double> > tmp_flow = flowNtoM[j];
          pair<double,double> tmp_overlap = overlapNtoM[j];
          flowNtoM[j] = flowNtoM[randPos];
          overlapNtoM[j] = overlapNtoM[randPos];
          flowNtoM[randPos] = tmp_flow;
          overlapNtoM[randPos] = tmp_overlap;
          // int randPos = j + R->randInt(NmodLinks-j-1);
          // int tmp_M = flowNtoM[j].first;
          // double tmp_outFlow = flowNtoM[j].second.first;
          // double tmp_inFlow = flowNtoM[j].second.second;
          // pair<double,double> tmp_overlap = overlapNtoM[j];
          // flowNtoM[j].first = flowNtoM[randPos].first;
          // flowNtoM[j].second.first = flowNtoM[randPos].second.first;
          // flowNtoM[j].second.second = flowNtoM[randPos].second.second;
          // overlapNtoM[j] = overlapNtoM[randPos];
          // flowNtoM[randPos].first = tmp_M;
          // flowNtoM[randPos].second.first = tmp_outFlow;
          // flowNtoM[randPos].second.second = tmp_inFlow;
          // overlapNtoM[randPos] = tmp_overlap;
        }

        int bestM = oldM;
        double best_outFlow = 0.0;
        double best_inFlow = 0.0;
        double best_delta = 0.0;
        double best_delta_nodeSize_log_nodeSize = 0.0;

    // Find the move that minimizes the description length
        for (int j=0; j<NmodLinks; j++) {

          int newM = flowNtoM[j].first;
          double outFlowNewM = flowNtoM[j].second.first;
          double inFlowNewM = flowNtoM[j].second.second;

          if(newM != oldM){

            double delta_enter = plogp(enterFlow + outFlowOldM + inFlowOldM - outFlowNewM - inFlowNewM) - enter;

            double delta_enter_log_enter = - plogp(mod_enter[oldM]) - plogp(mod_enter[newM]) \
            + plogp(mod_enter[oldM] - node[flip]->enter + outFlowOldM + inFlowOldM) + plogp(mod_enter[newM] + node[flip]->enter - outFlowNewM - inFlowNewM);

            double delta_exit_log_exit = - plogp(mod_exit[oldM]) - plogp(mod_exit[newM]) \
            + plogp(mod_exit[oldM] - node[flip]->exit + outFlowOldM + inFlowOldM) + plogp(mod_exit[newM] + node[flip]->exit - outFlowNewM - inFlowNewM);

            double delta_size_log_size = - plogp(mod_exit[oldM] + mod_size[oldM]) - plogp(mod_exit[newM] + mod_size[newM]) \
            + plogp(mod_exit[oldM] + mod_size[oldM] - node[flip]->exit - node[flip]->size + outFlowOldM + inFlowOldM) \
            + plogp(mod_exit[newM] + mod_size[newM] + node[flip]->exit + node[flip]->size - outFlowNewM - inFlowNewM);

            double delta_nodeSize_log_nodeSize = overlapNtoOldM.first + overlapNtoM[j].first + overlapNtoOldM.second - overlapNtoM[j].second;

            double deltaL = delta_enter - delta_enter_log_enter - delta_exit_log_exit + delta_size_log_size - delta_nodeSize_log_nodeSize;

            // cout << oldM+1 << " --> " << newM+1 << " " << deltaL << endl;

            if(deltaL < best_delta){
              bestM = newM;
              best_outFlow = outFlowNewM;
              best_inFlow = inFlowNewM;
              best_delta = deltaL;
              best_delta_nodeSize_log_nodeSize = delta_nodeSize_log_nodeSize;
            }
          }
        }

    // Make best possible move
        if(bestM != oldM){

      //Update empty module vector
          if(mod_members[bestM] == 0){
            Nempty--;
          }
          if(mod_members[oldM] == static_cast<int>(node[flip]->members.size())){
            mod_empty[Nempty] = oldM;
            Nempty++;
          }

          enterFlow -= mod_enter[oldM] + mod_enter[bestM];
          enter_log_enter -= plogp(mod_enter[oldM]) + plogp(mod_enter[bestM]);
          exit_log_exit -= plogp(mod_exit[oldM]) + plogp(mod_exit[bestM]);
          size_log_size -= plogp(mod_exit[oldM] + mod_size[oldM]) + plogp(mod_exit[bestM] + mod_size[bestM]);

          mod_enter[oldM] -= node[flip]->enter - outFlowOldM - inFlowOldM;
          mod_exit[oldM] -= node[flip]->exit - outFlowOldM - inFlowOldM;
          mod_size[oldM] -= node[flip]->size;
          mod_members[oldM] -= node[flip]->members.size();
          mod_enter[bestM] += node[flip]->enter - best_outFlow - best_inFlow;
          mod_exit[bestM] += node[flip]->exit - best_outFlow - best_inFlow;
          mod_size[bestM] += node[flip]->size;
          mod_members[bestM] += node[flip]->members.size();

      // For all multiple assigned nodes
          int flipNphysicalNode = node[flip]->physicalNodes.size();
          for(int j=0;j<flipNphysicalNode;j++){
            int physicalNode = node[flip]->physicalNodes[j].first;
            double physicalNodeSize = node[flip]->physicalNodes[j].second;

        // Remove contribution to old module
            map<int,pair<int,double> >::iterator overlap_it = mod_nodeSize[physicalNode].find(oldM);
            if(overlap_it == mod_nodeSize[physicalNode].end() ){
              cout << "ERROR: Could not find module " << oldM << " among physical node " << physicalNode << "'s assignments." << endl;
              abort();
            }
            else{
              overlap_it->second.second -= physicalNodeSize;
              overlap_it->second.first--;
              if(overlap_it->second.first == 0)
                mod_nodeSize[physicalNode].erase(overlap_it);
            }

        // Add contribution to new module
            overlap_it = mod_nodeSize[physicalNode].find(bestM);
            if(overlap_it == mod_nodeSize[physicalNode].end() ){
              mod_nodeSize[physicalNode].insert(make_pair(bestM,make_pair(1,physicalNodeSize)));
            }
            else{
              overlap_it->second.second += physicalNodeSize;
              overlap_it->second.first++;
            }
          }

          enterFlow += mod_enter[oldM] + mod_enter[bestM];

      // Update terms in map equation
          enter_log_enter += plogp(mod_enter[oldM]) + plogp(mod_enter[bestM]);
          exit_log_exit += plogp(mod_exit[oldM]) + plogp(mod_exit[bestM]);
          size_log_size += plogp(mod_exit[oldM] + mod_size[oldM]) + plogp(mod_exit[bestM] + mod_size[bestM]);

          nodeSize_log_nodeSize += best_delta_nodeSize_log_nodeSize;

          enter = plogp(enterFlow);

      // Update code length

          codeLength = enter - enter_log_enter - exit_log_exit + size_log_size - nodeSize_log_nodeSize;

          node[flip]->index = bestM;
          moved++;
        }

        offset += Nnode;

      }
    }

    void Greedy::initiate(void){

  // Take care of dangling nodes, normalize outLinks, and calculate total teleport weight
      for(int i=0;i<Nnode;i++){

        if(!node[i]->outLinks.empty() || (node[i]->selfLink > 0.0)){
          int NoutLinks = node[i]->outLinks.size();
      double sum = node[i]->selfLink; // Take care of self-links
      for(int j=0;j < NoutLinks; j++)
        sum += node[i]->outLinks[j].second;
      node[i]->selfLink /= sum;
      for(int j=0;j < NoutLinks; j++){
        node[i]->outLinks[j].second /= sum;
      }
    }
  }

  //   if(Ndanglings > 0){
  //     cout << "Found " << Ndanglings << " dangling node(s)." << endl;
  //   }


  // Calculate steady state matrix

  eigenvector();

  // Store the PageRank sizes
  vector<double> pr_size(Nnode);
  for(int i=0;i<Nnode;i++)
    pr_size[i] = node[i]->size;

  // Take eigenfactor step to update sizes
  eigenfactor();

  // Normalize such that code length measure average description length of encoded steps (and not encoded and non-encoded steps)
  double normalizeFactor = 0.0;
  for(int i=0;i<Nnode;i++)
    normalizeFactor += node[i]->size;

  for(int i=0;i<Nnode;i++){
    pr_size[i] /= normalizeFactor;
    node[i]->size /= normalizeFactor;
    // cout << i+1 << " " << pr_size[i] << " " << node[i]->size << endl;
  }

  // Set sizes of state nodes
  for(int i=0;i<Nnode;i++){
    if(node[i]->physicalNodes.size() > 0){
      if(node[i]->physicalNodes.size() > 1)
        cout << "Strange, the number of members is more than 1." << endl;

      node[i]->physicalNodes[0].second = node[i]->size;
    }
  }

  // Update links to represent flow (based on PageRank frequencies)
  for(int i=0;i<Nnode;i++){
    node[i]->selfLink = pr_size[i]*node[i]->selfLink;
    if(!node[i]->outLinks.empty()){
      int NoutLinks = node[i]->outLinks.size();
      for(int j=0;j < NoutLinks; j++){
        node[i]->outLinks[j].second = pr_size[i]*node[i]->outLinks[j].second;
      }

      // Update values for corresponding inlinks
      for(int j=0; j < NoutLinks; j++){
        int NinLinks = node[node[i]->outLinks[j].first]->inLinks.size();
        for(int k=0; k < NinLinks; k++){
          if(node[node[i]->outLinks[j].first]->inLinks[k].first == i){
            node[node[i]->outLinks[j].first]->inLinks[k].second = node[i]->outLinks[j].second;
            break;
          }
        }
      }
    }
  }


    // Calculate enter- and exitflow
  for(int i=0;i<Nnode;i++){
    node[i]->exit = 0.0;
    node[i]->enter = 0.0;
  }
  for(int i=0;i<Nnode;i++){
    int NoutLinks = node[i]->outLinks.size();
    for(int j=0;j<NoutLinks;j++){
      node[i]->exit += node[i]->outLinks[j].second;
      node[node[i]->outLinks[j].first]->enter += node[i]->outLinks[j].second;
    }
  }

  // Calculate p log p over all nodes
  nodeSize_log_nodeSize = 0.0;
  for(int i=0;i<Nnode;i++)
    nodeSize_log_nodeSize += plogp(node[i]->size);

  calibrate();

  //  cout << "Initial bit rate is " << codeLength/log(2.0) << ", starting merging " << Nnode << " nodes..." << endl;

}

void Greedy::tune(void){

  enter_log_enter = 0.0;
  exit_log_exit = 0.0;
  size_log_size = 0.0;
  enterFlow = 0.0;

  for(int i=0;i<Nmod;i++){
    mod_enter[i] = 0.0;
    mod_exit[i] = 0.0;
    mod_size[i] = 0.0;
    mod_members[i] = 0;
  }

  // Update all values except contribution from teleportation
  for(int i=0;i<Nnode;i++){
    int i_M = node[i]->index;
    int Nlinks = node[i]->outLinks.size();
    mod_size[i_M] += node[i]->size;
    mod_members[i_M]++;
    for(int j=0;j<Nlinks;j++){
      int nb = node[i]->outLinks[j].first;
      double nb_w = node[i]->outLinks[j].second;
      int nb_M = node[nb]->index;
      if(i_M != nb_M){
       mod_exit[i_M] += nb_w;
       mod_enter[nb_M] += nb_w;
     }
   }
 }

 for(int i=0;i<Nmod;i++){
  enter_log_enter += plogp(mod_enter[i]);
  exit_log_exit += plogp(mod_exit[i]);
  size_log_size += plogp(mod_exit[i] + mod_size[i]);
  enterFlow += mod_enter[i];
}

enter = plogp(enterFlow);

codeLength = enter - enter_log_enter - exit_log_exit + size_log_size - nodeSize_log_nodeSize;

}

void Greedy::calibrate(void){

  vector<int>(Nmod).swap(mod_empty);
  Nempty = 0;

  vector<double>(Nmod).swap(mod_enter);
  vector<double>(Nmod).swap(mod_exit);
  vector<double>(Nmod).swap(mod_size);
  vector<int>(Nmod).swap(mod_members);

    vector<map<int,pair<int,double> > >(NphysicalNode).swap(mod_nodeSize); // What happens with the vectors within the vector?

    enter_log_enter = 0.0;
    exit_log_exit = 0.0;
    size_log_size = 0.0;
    enterFlow = 0.0;

    for(int i=0;i<Nmod;i++){

      enter_log_enter += plogp(node[i]->enter);
      exit_log_exit += plogp(node[i]->exit);
      size_log_size += plogp(node[i]->exit + node[i]->size);
      enterFlow += node[i]->enter;

      mod_enter[i] = node[i]->enter;
      mod_exit[i] = node[i]->exit;
      mod_size[i] = node[i]->size;
      mod_members[i] = node[i]->members.size();
      node[i]->index = i;

      int Nmem = node[i]->physicalNodes.size();
      for(int j=0;j<Nmem;j++){
        int overlapModule = node[i]->physicalNodes[j].first;
        mod_nodeSize[overlapModule].insert(mod_nodeSize[overlapModule].end(),make_pair(i,make_pair(1,node[i]->physicalNodes[j].second)));
      }
    }

    enter = plogp(enterFlow);

    nodeSize_log_nodeSize = 0.0;
    for(int i=0;i<NphysicalNode;i++)
      for(map<int,pair<int,double> >::iterator overlap_it = mod_nodeSize[i].begin(); overlap_it != mod_nodeSize[i].end(); ++overlap_it)
        nodeSize_log_nodeSize += plogp(overlap_it->second.second);

      codeLength = enter - enter_log_enter - exit_log_exit + size_log_size - nodeSize_log_nodeSize;

    }

    void Greedy::prepare(bool sort){

      Nmod = 0;
      vector<int>().swap(modSnode);

      if(sort){

        multimap<double,int> Msize;
        for(int i=0;i<Nnode;i++){
          if(mod_members[i] > 0){
            Nmod++;
            Msize.insert(make_pair(mod_size[i],i));
          }
        }

        for(multimap<double,int>::reverse_iterator it = Msize.rbegin(); it != Msize.rend(); it++)
          modSnode.push_back(it->second);

      }
      else{

        for(int i=0;i<Nnode;i++){
          if(mod_members[i] > 0){
            Nmod++;
            modSnode.push_back(i);
          }
        }

      }

    }

    void Greedy::level(Node ***node_tmp, bool sort){

 // cout << "Level! " << endl;

      prepare(sort);

  //Node ***ntmp = node_tmp;

      (*node_tmp) = new Node*[Nmod];

      vector<int> nodeInMod = vector<int>(Nnode);

      for(int i=0;i<Nmod;i++){
        (*node_tmp)[i] = new Node();
        (*node_tmp)[i]->index = i;
        (*node_tmp)[i]->enter = mod_enter[modSnode[i]];
        (*node_tmp)[i]->exit = mod_exit[modSnode[i]];
        (*node_tmp)[i]->size = mod_size[modSnode[i]];
        nodeInMod[modSnode[i]] = i;
      }

  // Calculate outflow of links to different modules
      vector<map<int,double> > outFlowNtoM(Nmod);
      map<int,double>::iterator it_M;

      for(int i=0;i<Nnode;i++){

        int i_M = nodeInMod[node[i]->index];

        (*node_tmp)[i_M]->members.insert((*node_tmp)[i_M]->members.end(),node[i]->members.begin(),node[i]->members.end());

    // copy(node[i]->members.begin(),node[i]->members.end(),back_inserter((*node_tmp)[i_M]->members));

        int NoutLinks = node[i]->outLinks.size();
        for(int j=0;j<NoutLinks;j++){
          int nb = node[i]->outLinks[j].first;
          int nb_M = nodeInMod[node[nb]->index];
          double nb_flow = node[i]->outLinks[j].second;
          if (nb != i) {
            it_M = outFlowNtoM[i_M].find(nb_M);
            if (it_M != outFlowNtoM[i_M].end())
              it_M->second += nb_flow;
            else
              outFlowNtoM[i_M].insert(make_pair(nb_M,nb_flow));
          }
        }
      }

  // Create outLinks at new level
      for(int i=0;i<Nmod;i++){
        for(it_M = outFlowNtoM[i].begin(); it_M != outFlowNtoM[i].end(); it_M++){
          if(it_M->first != i){
            (*node_tmp)[i]->outLinks.push_back(make_pair(it_M->first,it_M->second));
          }
        }
      }


  // Calculate inflow of links from different modules
      vector<map<int,double> > inFlowNtoM(Nmod);

      for(int i=0;i<Nnode;i++){

        int i_M = nodeInMod[node[i]->index];

        int NinLinks = node[i]->inLinks.size();
        for(int j=0;j<NinLinks;j++){
          int nb = node[i]->inLinks[j].first;
          int nb_M = nodeInMod[node[nb]->index];
          double nb_flow = node[i]->inLinks[j].second;
          if (nb != i) {
            it_M = inFlowNtoM[i_M].find(nb_M);
            if (it_M != inFlowNtoM[i_M].end())
              it_M->second += nb_flow;
            else
              inFlowNtoM[i_M].insert(make_pair(nb_M,nb_flow));
          }
        }
      }

  // Create inLinks at new level
      for(int i=0;i<Nmod;i++){
        for(it_M = inFlowNtoM[i].begin(); it_M != inFlowNtoM[i].end(); it_M++){
          if(it_M->first != i){
            (*node_tmp)[i]->inLinks.push_back(make_pair(it_M->first,it_M->second));
          }
        }
      }

//   set<int> validModules;
// for(int i=0;i<Nmod;i++)
// validModules.insert(modSnode[i]);

//   // Check physicalNodes
//   if(mod_nodeSize.size() > 0){
//   for(int i=0;i<NphysicalNode;i++){
//     int Noverlaps = mod_nodeSize[i].size();
//     map<int,int> count;
//     int duplicate;
//         set<int>::iterator valid_it;
//     for(int j=0;j<Noverlaps;j++){
//       valid_it = validModules.find(mod_nodeSize[i][j].first);
//       if(valid_it == validModules.end())
//         cout << validModules.size() << "Physical node " << i << " is assigned to empty module " << mod_nodeSize[i][j].first << " ";

//       duplicate = count[mod_nodeSize[i][j].first] += 1;
//       if(duplicate > 1){
//           cout << "ERROR IN LEVEL! Physical node: " << i << endl;
//           for(int k=0;k<Noverlaps;k++)
//             cout << mod_nodeSize[i][k].first << " ";
//           abort();
//       }
//     }
//   }
// }

  // for(int i=0;i<Nnode;i++){
  //   int Noverlaps = node[i]->physicalNodes.size();
  //   map<int,int> count;
  //   int duplicate;
  //   for(int j=0;j<Noverlaps;j++){
  //     duplicate = count[node[i]->physicalNodes[j].first] += 1;
  //     if(duplicate > 1){
  //         cout << "ERROR BEFORE UPDATE IN NODES! Physical node: " << i << endl;
  //         for(int k=0;k<Noverlaps;k++)
  //           cout << node[i]->physicalNodes[k].first << " ";
  //         abort();
  //     }
  //   }
  // }

  // Update physicalNodes
      map<int,map<int,int> > validate;

      for(int i=0;i<NphysicalNode;i++){
        for(map<int,pair<int,double> >::iterator overlap_it = mod_nodeSize[i].begin(); overlap_it != mod_nodeSize[i].end(); ++overlap_it){
          if(++validate[overlap_it->first][i] > 1){
            cout << "DUPLICATION ERROR" << endl;
            abort();

          }

          (*node_tmp)[nodeInMod[overlap_it->first]]->physicalNodes.push_back(make_pair(i,overlap_it->second.second));

        }

      }

  //   for(int i=0;i<Nmod;i++){
  //   int Noverlaps = (*node_tmp)[i]->physicalNodes.size();
  //   map<int,int> count;
  //   int duplicate;
  //   for(int j=0;j<Noverlaps;j++){
  //     duplicate = count[(*node_tmp)[i]->physicalNodes[j].first] += 1;
  //     if(duplicate > 1){
  //         cout << "ERROR AFTER UPDATE IN NODES! Physical node: " << i << "-" << (*node_tmp)[i]->physicalNodes[j].first << endl;
  //         for(int k=0;k<Noverlaps;k++)
  //           cout << (*node_tmp)[i]->physicalNodes[k].first << " ";
  //         abort();
  //     }
  //   }
  // }



  // Option to move to empty module
      vector<int>().swap(mod_empty);
      Nempty = 0;
      for(int i=0;i<Nnode;i++){
        delete node[i];
      }
      delete [] node;

      Nnode = Nmod;
      node = (*node_tmp);

//   for(int i=0;i<Nnode;i++){
//    int Noverlaps = node[i]->physicalNodes.size();
//       for(int j=0;j<Noverlaps;j++){
//      cout << i << " " << node[i]->physicalNodes[j].first << " " << node[i]->physicalNodes[j].second << endl;
//       }
//   }

      calibrate();

    }


    void Greedy::determMove(vector<int> &moveTo){

      for(int i=0;i<Nnode;i++){
        int oldM = i;
        int newM = moveTo[i];

        if(newM != oldM){

          double outFlowOldM = 0.0;
          double inFlowOldM = 0.0;
          double outFlowNewM = 0.0;
          double inFlowNewM = 0.0;

      // For all outLinks
          int NoutLinks = node[i]->outLinks.size();
          for(int j=0; j<NoutLinks; j++){
            int nb_M = node[node[i]->outLinks[j].first]->index;
            double nb_flow = node[i]->outLinks[j].second;
            if(nb_M == oldM){
              outFlowOldM += nb_flow;
            }
            else if(nb_M == newM){
              outFlowNewM += nb_flow;
            }
          }

      // For all inLinks
          int NinLinks = node[i]->inLinks.size();
          for(int j=0; j<NinLinks; j++){
            int nb_M = node[node[i]->inLinks[j].first]->index;
            double nb_flow = node[i]->inLinks[j].second;
            if(nb_M == oldM){
              inFlowOldM += nb_flow;
            }
            else if(nb_M == newM){
              inFlowNewM += nb_flow;
            }
          }

      //Update empty module vector
          if(mod_members[newM] == 0){
            Nempty--;
          }
          if(mod_members[oldM] == static_cast<int>(node[i]->members.size())){
            mod_empty[Nempty] = oldM;
            Nempty++;
          }

      // For all multiple assigned nodes
          double delta_nodeSize_log_nodeSize = 0.0;

          int iNphysicalNode = node[i]->physicalNodes.size();
          for(int j=0;j<iNphysicalNode;j++){
            int physicalNode = node[i]->physicalNodes[j].first;
            int Nassignments = mod_nodeSize[physicalNode].size();
        // Measure changes in nodeSize_log_nodeSize
            double physicalNodeSize = node[i]->physicalNodes[j].second;

        // Remove contribution to old module
            map<int,pair<int,double> >::iterator overlap_it = mod_nodeSize[physicalNode].find(oldM);
            if(overlap_it == mod_nodeSize[physicalNode].end() ){
              cout << "ERROR: Could not find module " << oldM << " among physical node " << physicalNode << "'s assignments." << endl;
              abort();
            }
            else{
              double oldP = overlap_it->second.second;
              double newP = overlap_it->second.second - physicalNodeSize;
              delta_nodeSize_log_nodeSize += plogp(newP) - plogp(oldP);
              overlap_it->second.second -= physicalNodeSize;
              overlap_it->second.first--;
              if(overlap_it->second.first == 0)
                mod_nodeSize[physicalNode].erase(overlap_it);
            }

        // Add contribution to new module
            overlap_it = mod_nodeSize[physicalNode].find(newM);
            if(overlap_it == mod_nodeSize[physicalNode].end() ){

              delta_nodeSize_log_nodeSize += plogp(physicalNodeSize);
              mod_nodeSize[physicalNode].insert(make_pair(newM,make_pair(1,physicalNodeSize)));

            }
            else{

              double oldP = overlap_it->second.second;
              double newP = overlap_it->second.second + physicalNodeSize;
              delta_nodeSize_log_nodeSize += plogp(newP) - plogp(oldP);
              overlap_it->second.second += physicalNodeSize;
              overlap_it->second.first++;

            }

          }

          enterFlow -= mod_enter[oldM] + mod_enter[newM];
          enter_log_enter -= plogp(mod_enter[oldM]) + plogp(mod_enter[newM]);
          exit_log_exit -= plogp(mod_exit[oldM]) + plogp(mod_exit[newM]);
          size_log_size -= plogp(mod_exit[oldM] + mod_size[oldM]) + plogp(mod_exit[newM] + mod_size[newM]);

          mod_enter[oldM] -= node[i]->enter - outFlowOldM - inFlowOldM;
          mod_exit[oldM] -= node[i]->exit - outFlowOldM - inFlowOldM;
          mod_size[oldM] -= node[i]->size;
          mod_members[oldM] -= node[i]->members.size();
          mod_enter[newM] += node[i]->enter - outFlowNewM - inFlowNewM;
          mod_exit[newM] += node[i]->exit - outFlowNewM - inFlowNewM;
          mod_size[newM] += node[i]->size;
          mod_members[newM] += node[i]->members.size();

          enterFlow += mod_enter[oldM] + mod_enter[newM];
          enter_log_enter += plogp(mod_enter[oldM]) + plogp(mod_enter[newM]);
          exit_log_exit += plogp(mod_exit[oldM]) + plogp(mod_exit[newM]);
          size_log_size += plogp(mod_exit[oldM] + mod_size[oldM]) + plogp(mod_exit[newM] + mod_size[newM]);
          nodeSize_log_nodeSize += delta_nodeSize_log_nodeSize;
          enter = plogp(enterFlow);

          codeLength = enter - enter_log_enter - exit_log_exit + size_log_size - nodeSize_log_nodeSize;

          node[i]->index = newM;

        }

      }

    };

    void Greedy::eigenvector(void){

  // cout << "Calculating steady state distribution of flow...";

      vector<double> size_tmp = vector<double>(Nnode);
      for(int i=0;i<Nnode;i++){
        size_tmp[i] = node[i]->teleportWeight;
      }

      int Niterations = 0;
      double sqdiff = 1.0;
      double sqdiff_old;
      double sum;

      double danglingSize;
      int Ndanglings = 0;
      vector<int> danglings;
      for(int i=0;i<Nnode;i++){
        if(node[i]->outLinks.empty() && (node[i]->selfLink <= 0.0)){
          danglings.push_back(i);
          Ndanglings++;
        }
      }

      do{

    // Calculate dangling size
        danglingSize = 0.0;
        for(int i=0;i<Ndanglings;i++){
          danglingSize += size_tmp[danglings[i]];
        }

    // Flow from teleportation
        for(int i=0;i<Nnode;i++)
          node[i]->size = (alpha + beta*danglingSize)*node[i]->teleportWeight;

    // Flow from network steps
        for(int i=0;i<Nnode;i++){
          node[i]->size += beta*node[i]->selfLink*size_tmp[i];
          int Nlinks = node[i]->outLinks.size();
          for(int j=0; j < Nlinks; j++)
            node[node[i]->outLinks[j].first]->size += beta*node[i]->outLinks[j].second*size_tmp[i];
        }

    // Normalize
        sum = 0.0;
        for(int i=0;i<Nnode;i++){
          sum += node[i]->size;
        }
        sqdiff_old = sqdiff;
        sqdiff = 0.0;
        for(int i=0;i<Nnode;i++){
          node[i]->size /= sum;
          sqdiff += fabs(node[i]->size - size_tmp[i]);
          size_tmp[i] = node[i]->size;
        }
        Niterations++;

        if(sqdiff == sqdiff_old){
      //    fprintf(stderr,"\n1.0e-10 added to alpha for convergence (precision error)\n");
          alpha += 1.0e-10;
          beta = 1.0-alpha;
        }

    //fprintf(stderr,"\rCalculating steady state distribution of flow...the error is %e after %d iterations\n",sqdiff,Niterations);

      }  while((Niterations < 200) && (sqdiff > 1.0e-15 || Niterations < 50));

      danglingSize = 0.0;
      for(int i=0;i<Ndanglings;i++){
        danglingSize += size_tmp[danglings[i]];
      }

  //cout << "done! (the error is " << sqdiff << " after " << Niterations << " iterations)" << endl;

    }

    void Greedy::eigenfactor(void){

      vector<double> size_tmp = vector<double>(Nnode);
      for(int i=0;i<Nnode;i++){
        size_tmp[i] = node[i]->size;
        node[i]->size = 0.0;
      }

  // Flow from network steps
      for(int i=0;i<Nnode;i++){
        node[i]->size += node[i]->selfLink*size_tmp[i];
        int Nlinks = node[i]->outLinks.size();
        for(int j=0; j < Nlinks; j++)
          node[node[i]->outLinks[j].first]->size += node[i]->outLinks[j].second*size_tmp[i];
      }

    }
