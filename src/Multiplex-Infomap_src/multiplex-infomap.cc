#include "multiplex-infomap.h"

using namespace std;
using std::cout;
using std::cin;
using std::endl;

unsigned stou(char *s){
  return strtoul(s,(char **)NULL,10);
}

void printTree(string s,multimap<double,treeNode,greater<double> >::iterator it_tM,ofstream *outfile,bool flip);
void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent,int Ntrials,vector<int> &nodeassignments);
void repeated_partition_init(MTRand *R, Node ***node, GreedyBase *greedy, Node ***node_init, bool silent,int Ntrials,vector<int> &nodeassignments,map<int,int> &layerNetworkNodes);

void partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent);

  // Call: trade <seed> <Ntries>
int main(int argc,char *argv[]){

  cout << "Version: June 20, 2014." << endl;
  cout << "Command: ";
  cout << argv[0];
  for(int i=1;i<argc; i++)
    cout << " " << argv[i];
  cout << endl;

  if( argc == 1 ){
    cout << "Call: ./multiplex-infomap [-s <seed>] [-N <attempts>] [-teleportrate <rate>] [-switchrate <rate>] [-proportionalswitch|-uniformswitch] [-noselflinks|-selflinks] [-multiplex|-aggregated|-expanded] [-classical|-physical] [-clusters <nodeassignments.txt>|-guess <nodeassignments.txt>] [-smartinit] [-nodenames <nodenames.txt>] network_level1.net network_level2.net [network_level3.net ... network_levelN.net] level_network.net" << endl;
    exit(-1);
  }

  unsigned int seed = 1234;
  int Ntrials = 1;

  bool selflink = false;
  bool guess = false;
  bool smartinit = false;
  double alpha = 0.15;
  double switchrate = 0.15;
  NetworkMode netMode = MULTIPLEX;
  DynamicMode dynMode = PHYSICAL;
  SwitchMode switchMode = LINKS;

  int Nlayers = 0;
  vector<string> networkLayerFiles;
  string layerNetworkFile;

  string nodenamesFile = "nonames";
  string nodeassignmentsFile = "noassignments";

  int argNr = 1;
  while(argNr < argc){
    if(to_string(argv[argNr]) == "-h"){
      cout << "Call: ./multiplex-infomap [-s <seed>] [-N <attempts>] [-teleportrate <rate>] [-switchrate <rate>] [-proportionalswitch|-uniformswitch] [-noselflinks|-selflinks] [-multiplex|-aggregate|-expand] [-classical|-physical] [-clusters <nodeassignments.txt>|-guess <nodeassignments.txt>] [-smartinit] [-nodenames <nodenames.txt>] network_level1.net network_level2.net [network_level3.net ... network_levelN.net] level_network.net" << endl;
      exit(-1);
    }
    else if(to_string(argv[argNr]) == "-s"){
      argNr++;
      seed = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-N"){
      argNr++;
      Ntrials = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-teleportrate"){
      argNr++;
      alpha = atof(argv[argNr]);
      if(alpha < 0.0 || alpha > 0.999){
        cout << "Input error: " << alpha << " is illegal value of -teleportrate, use value in range (0,1). Aborting..." << endl;
        exit(-1);
      }
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-switchrate"){
      switchMode = PROPORTIONAL;
      argNr++;
      switchrate = atof(argv[argNr]);
      if(switchrate < 0.0 || switchrate > 1){
        cout << "Input error: " << alpha << " is illegal value of -switchrate, use value in range [0,1]. Aborting..." << endl;
        exit(-1);
      }
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-proportionalswitch"){
      switchMode = PROPORTIONAL;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-uniformswitch"){
      switchMode = UNIFORM;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-aggregated"){
      netMode = AGGREGATED;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-expanded"){
      netMode = EXPANDED;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-multiplex"){
      netMode = MULTIPLEX;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-physical"){
      dynMode = PHYSICAL;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-classical"){
      dynMode = CLASSICAL;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-noselflinks"){
      selflink = false;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-selflinks"){
      selflink = true;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-smartinit"){
      smartinit = true;
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-nodenames"){
      argNr++;
      nodenamesFile = string(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-clusters"){
      argNr++;
      nodeassignmentsFile = string(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-guess"){
      guess = true;
      argNr++;
      nodeassignmentsFile = string(argv[argNr]);
      argNr++;
    }
    else{

      if(argv[argNr][0] == '-'){
        cout << "Unknown command: " << to_string(argv[argNr]) << endl;
        cout << "Call: ./multiplex-infomap [-s <seed>] [-N <attempts>] [-teleportrate <rate>] [-switchrate <rate>] [-proportionalswitch|-uniformswitch] [-noselflinks|-selflinks] [-multiplex|-aggregate|-expand] [-classical|-physical] [-clusters <nodeassignments.txt>|-guess <nodeassignments.txt>] [-nodenames <nodenames.txt>] network_level1.net network_level2.net [network_level3.net ... network_levelN.net] level_network.net" << endl;
        exit(-1);
      }

      while(argNr < (argc - 1)){
        networkLayerFiles.push_back(string(argv[argNr]));
        Nlayers++;
        argNr++;
      }
      layerNetworkFile = string(argv[argNr]);
      argNr++;
    }
  }

  cout << "Setup:" << endl;
  cout << "-->Using seed: " << seed << endl;
  cout << "-->Using teleportation rate: " << alpha << " (teleportatin steps are not encoded)" << endl;
  cout << "-->Number of trials: " << Ntrials << endl;
  if(netMode == AGGREGATED){
    cout << "-->Analysis with aggregated layers and ignored inter-layer links, using standard random walker, " << endl;
  }
  else{
    if(netMode == EXPANDED){
      cout << "-->Analysis with expanded network, using standard random walker, " << endl;
    }
    else if(netMode == MULTIPLEX){
      cout << "-->Analysis with multiplex network structure in multiple layers, ";
      if(dynMode == CLASSICAL){
        cout << "using classical random walker, ";
      }
      else if(dynMode == PHYSICAL){
        cout << "using physical random walker, ";
      }
      if(switchMode == PROPORTIONAL){
        cout << "and inter-layer switch rate " << switchrate << " porportional to link weights, ";
      }
      else if(switchMode == UNIFORM){
        cout << "and inter-layer switch rate " << switchrate << " uniformly over layers, ";
      }
      else{
        cout << "and inter-layer link weights from link list, ";
      }
    }
  }
  if(smartinit){
    cout << "with additional smart initialization from undirected expanded solution, ";
  }
  else{
   cout << "without smart initialization, "; 
 }
 if(!selflink){
  cout << "without self-links." << endl;
}
else{
  cout << "with self-links." << endl;
}
if(nodenamesFile == "nonames"){
  cout << "-->No node-names file provided. Using node id." << endl;
}
else{
  cout << "-->Using node names from file: " << nodenamesFile << endl;
}
if(nodeassignmentsFile == "noassignments"){
  cout << "-->No node assignmnets provided. Performing community detection." << endl;
}
else{
  cout << "-->Using node assigments from file: " << nodeassignmentsFile;
  if(guess){
    cout << " as an initial guess." << endl;
  }
  else{
    cout << " as the final solution." << endl;
  }
}

for(int i=0;i<Nlayers;i++)
  cout << "-->Network layer " << i+1 << " from file: " << networkLayerFiles[i].c_str() << endl;
if(Nlayers == 0)
  cout << "-->Multiplex network from file: " << layerNetworkFile.c_str() << endl;
else
  cout << "-->Layer network from file: " << layerNetworkFile.c_str() << endl;

MTRand *R = new MTRand(seed);

  // Extract layer network name without extension to be used as stem in names of outfiles
string networkName = layerNetworkFile;
size_t lastPeriod = networkName.find_last_of(".");
if(lastPeriod != string::npos)
 networkName = string(networkName.begin(),networkName.begin() + lastPeriod);

Multiplex multiplex(networkLayerFiles,layerNetworkFile,nodenamesFile,nodeassignmentsFile);

loadMnetwork(multiplex,netMode,selflink,switchMode,switchrate);

vector<pair<int,int> > reverseMNodeMap(multiplex.MNnode);

 // Partition network
Node **node = new Node*[multiplex.MNnode];
for(map<pair<int,int>,int>::iterator it = multiplex.MnodeMap.begin(); it != multiplex.MnodeMap.end(); it++){
  node[it->second] = new Node(it->second,multiplex.MnodeWeights[it->second]/multiplex.totMNodeWeights);
  if(netMode == EXPANDED)
      node[it->second]->physicalNodes.push_back(make_pair(it->second,0.0)); // Each state node represents its own fake physical node
    else
      node[it->second]->physicalNodes.push_back(make_pair(it->first.second,0.0)); // State nodes represent their real physical nodes

    reverseMNodeMap[it->second] = it->first;
  }

  if(dynMode == PHYSICAL && netMode == MULTIPLEX){
    cout << "-->Modifying links for physical random walkers..." << flush;
    map<pair<pair<int,int>,pair<int,int> >,double> physicalMLinks;
    map<pair<int,int>,double > MWeights;
    map<pair<int,int>,vector<pair<pair<int,int>,double> > > SortedMLinks;
    for(map<pair<pair<int,int>,pair<int,int> >,double>::iterator it = multiplex.MLinks.begin(); it != multiplex.MLinks.end(); it++){
      if(it->first.first.first == it->first.second.first){ // Only consider links within the same layer
        if(it->first.first.second != it->first.second.second || (it->first.first.second == it->first.second.second&& selflink)){
          MWeights[it->first.first] += it->second;
          SortedMLinks[it->first.first].push_back(make_pair(it->first.second,it->second));
        }
      }
    }
    for(map<pair<pair<int,int>,pair<int,int> >,double>::iterator it = multiplex.MLinks.begin(); it != multiplex.MLinks.end(); it++){
      if((it->first.first.first != it->first.second.first) && (it->first.first.second == it->first.second.second)){ //Identify links between layers

        map<pair<int,int>,vector<pair<pair<int,int>,double> > >::iterator it_m = SortedMLinks.find(it->first.second);
        map<pair<int,int>,double >::iterator it_w = MWeights.find(it->first.second);
        if(it_m == SortedMLinks.end()){
          cout << "Error, could not find state node" << endl;
          exit(-1);
        }
        else{
          for(vector<pair<pair<int,int>,double> >::iterator it_l = it_m->second.begin(); it_l != it_m->second.end(); it_l++){
            physicalMLinks[make_pair(it->first.first,it_l->first)] += it->second*it_l->second/it_w->second;
          }
        }

      }
      else{
        physicalMLinks[it->first] += it->second;
      }
    }

    cout << "adding " << physicalMLinks.size() - multiplex.MLinks.size() << " links to a total of " << physicalMLinks.size() << " links..." << flush;
    physicalMLinks.swap(multiplex.MLinks);
    map<pair<pair<int,int>,pair<int,int> >,double>().swap(physicalMLinks);

    cout << "done!" << endl;

    // for(map<pair<pair<int,int>,pair<int,int> >,double>::iterator it = multiplex.MLinks.begin(); it != multiplex.MLinks.end(); it++){
    //   if(it->first.first == make_pair(0,0))
    //     cout << it->first.first.first << " " << it->first.first.second << " --> " << it->first.second.first << " " << it->first.second.second << " " << it->second << endl;
    // }

    // for(map<pair<pair<int,int>,pair<int,int> >,double>::iterator it = multiplex.MLinks.begin(); it != multiplex.MLinks.end(); it++){
    //   if(it->first.first == make_pair(1,0))
    //     cout << it->first.first.first << " " << it->first.first.second << " --> " << it->first.second.first << " " << it->first.second.second << " " << it->second << endl;
    // }

  }

  // for(map<pair<pair<int,int>,pair<int,int> >,double>::iterator it = multiplex.MLinks.begin(); it != multiplex.MLinks.end(); it++){
  //   if(it->first.first.first == 0)
  //     cout << it->first.first.second+1 << " " << it->first.second.second+1 << " " << it->second << endl;
  // }



  int NselfLinks = 0;
  for(map<pair<pair<int,int>,pair<int,int> >,double>::iterator it = multiplex.MLinks.begin(); it != multiplex.MLinks.end(); it++){

    int from = multiplex.MnodeMap[it->first.first];
    int to = multiplex.MnodeMap[it->first.second];

    double weight = it->second;

    if(weight > 0.0){
      if(from == to){
        if(selflink)
          node[from]->selfLink += weight;
        NselfLinks++;
      }
      else{
        node[from]->outLinks.push_back(make_pair(to,weight));
        node[to]->inLinks.push_back(make_pair(from,weight));
      }
    }
  }

  if(selflink)
    cout << "-->Included " <<  NselfLinks << " within-layer self-link(s)." << endl;
  else
    cout << "-->Ignored " <<  NselfLinks << " within-layer self-link(s)." << endl;

  //Swap vector to free memory
  map<pair<int,int>,double>().swap(multiplex.Links);
  map<pair<pair<int,int>,pair<int,int> >,double>().swap(multiplex.MLinks);

   // Initiation
  GreedyBase* greedy;
  if(netMode == EXPANDED)
    greedy = new Greedy(R,multiplex.MNnode,node,multiplex.MNnode,alpha);
  else
    greedy = new Greedy(R,multiplex.MNnode,node,multiplex.Nnode,alpha);

  greedy->initiate();

  // vector<double> nodeSize(multiplex.Nnode,0.0);
  // for(int i=0;i<multiplex.MNnode;i++)
  // nodeSize[reverseMNodeMap[i].second] += node[i]->size;

  //  //Print flow in node order
  // ofstream outfile;
  // ostringstream oss;
  // oss.str("");
  // if(netMode == MULTIPLEX)
  //   oss << networkName << "_Multiplex_flow.txt";
  // else if(netMode == AGGREGATED)
  //   oss << networkName << "_Aggregated_flow.txt";
  // else if(netMode == EXPANDED)
  //   oss << networkName << "_Expanded_flow.txt";
  // outfile.open(oss.str().c_str());
  // for(int i=0;i<multiplex.Nnode;i++)
  // outfile << "\"" << multiplex.nodeNames[i] << "\" " << nodeSize[i] << endl;
  // outfile.close();

  double withinModuleFlow = 0.0;
  double betweenModuleFlow = 0.0;

  cout << "Community detection:" << endl;
  if(multiplex.nodeassignments.size() == multiplex.MNnode && !guess){
    cout << "-->Using nodes assignments from: " << nodeassignmentsFile << endl;
    greedy->determMove(multiplex.nodeassignments);

    // Print out flow across modules
    for(int i=0;i<greedy->Nnode;i++){
      int NoutLinks = node[i]->outLinks.size();
      for(int j=0;j<NoutLinks;j++){
        if(node[i]->index == node[node[i]->outLinks[j].first]->index)
          withinModuleFlow += node[i]->outLinks[j].second;
        else
          betweenModuleFlow += node[i]->outLinks[j].second;
      }
    }
    cout << "-->Within/between-module flow: " << withinModuleFlow << "/" << betweenModuleFlow << endl;

    greedy->level(&node,true);
  }
  else{
    // cout << "-->Running community-detection algorithm." << endl;
    vector<int> noAssignments(0);

    // -------------- Start of alternative initialization ---------------------- //

    if(smartinit){

      cout << "-->Preparing smart initiation." << endl;

      string nodenamesFile_init = "nonames";
      string nodeassignmentsFile_init = "noassignments";
      Multiplex multiplex_init(networkLayerFiles,layerNetworkFile,nodenamesFile_init,nodeassignmentsFile_init);
      NetworkMode netMode_init = INIT;
      if(netMode == AGGREGATED)
        netMode_init = INITAGG;


      loadMnetwork(multiplex_init,netMode_init,selflink,switchMode,switchrate);

      // For creating separate networks of layers
      map<int,int> layerNetworkNodes;

    // Partition network
      Node **node_init = new Node*[multiplex_init.MNnode];
      for(map<pair<int,int>,int>::iterator it = multiplex_init.MnodeMap.begin(); it != multiplex_init.MnodeMap.end(); it++){
       node_init[it->second] = new Node(it->second,multiplex_init.MnodeWeights[it->second]/multiplex_init.totMNodeWeights);
     // As if expanded
       layerNetworkNodes[it->first.first] = it->second;
         node_init[it->second]->physicalNodes.push_back(make_pair(it->second,0.0)); // Each state node represents its own fake physical node
       }

       int NselfLinks = 0;
       for(map<pair<pair<int,int>,pair<int,int> >,double>::iterator it = multiplex_init.MLinks.begin(); it != multiplex_init.MLinks.end(); it++){

         int from = multiplex_init.MnodeMap[it->first.first];
         int to = multiplex_init.MnodeMap[it->first.second];

         double weight = it->second;

         if(weight > 0.0){
           if(from == to){
             if(selflink)
               node_init[from]->selfLink += weight;
             NselfLinks++;
           }
           else{
             node_init[from]->outLinks.push_back(make_pair(to,weight));
             node_init[to]->inLinks.push_back(make_pair(from,weight));
           }
         }
       }

       // if(selflink)
       //   cout << "-->Included " <<  NselfLinks << " within-layer self-link(s)." << endl;
       // else
       //   cout << "-->Ignored " <<  NselfLinks << " within-layer self-link(s)." << endl;

     //Swap vector to free memory
       map<pair<int,int>,double>().swap(multiplex_init.Links);
       map<pair<pair<int,int>,pair<int,int> >,double>().swap(multiplex_init.MLinks);

        cout << "-->Running community-detection algorithm." << endl;
       // repeated_partition(R,&node,greedy,false,Ntrials,multiplex.nodeassignments);
       if(multiplex.nodeassignments.size() == multiplex.MNnode){
         cout << "-->Using nodes assignments from: " << nodeassignmentsFile << endl;
         repeated_partition_init(R,&node,greedy,&node_init,false,Ntrials,multiplex.nodeassignments,layerNetworkNodes);
       }
       else{
        repeated_partition_init(R,&node,greedy,&node_init,false,Ntrials,noAssignments,layerNetworkNodes);
      }

      for(int i=0;i<multiplex_init.MNnode;i++){
       delete node_init[i];
     }
     delete [] node_init;

   }
   else{
    cout << "-->Running community-detection algorithm." << endl;

    // multiplex.nodeassignments = vector<int>(multiplex.MNnode,0);

    if(multiplex.nodeassignments.size() == multiplex.MNnode){
      cout << "-->Using nodes assignments from: " << nodeassignmentsFile << endl;
      repeated_partition(R,&node,greedy,false,Ntrials,multiplex.nodeassignments);
    }
    else{
     repeated_partition(R,&node,greedy,false,Ntrials,noAssignments);
   }


 }

       // -------------- End of alternative initialization ---------------------- //
}

int Nmod = greedy->Nnode;

  if(netMode == EXPANDED){ // Remap fake physical nodes to real physical nodes
    for(int i=0;i<Nmod;i++){
      map<int,double> realPhysicalNodesMap;
      int Nmembers = node[i]->physicalNodes.size();
      for(int j=0;j<Nmembers;j++){
        realPhysicalNodesMap[reverseMNodeMap[node[i]->physicalNodes[j].first].second] += node[i]->physicalNodes[j].second;
      }
      node[i]->physicalNodes.resize(realPhysicalNodesMap.size());
      int k=0;
      for(map<int,double>::iterator it = realPhysicalNodesMap.begin(); it != realPhysicalNodesMap.end(); it++){
        node[i]->physicalNodes[k] = make_pair(it->first,it->second);
        k++;
      }
    }
  }

   // Multiplex clustering
  vector<int> Mclustering(multiplex.MNnode);
  for(int i=0;i<Nmod;i++){
    int Nmembers = node[i]->members.size();
    for(int j=0;j<Nmembers;j++){
     Mclustering[node[i]->members[j]] = i;
   }
 }

 ofstream outfile;
 ostringstream oss;
 oss.str("");
 if(netMode == MULTIPLEX){
   oss << networkName << "_Multiplex";
   if(dynMode == CLASSICAL)
     oss << "_Classical.clu";
   else
     oss << "_Physical.clu";
 }
 else if(netMode == AGGREGATED)
   oss << networkName << "_Aggregated.clu";
 else if(netMode == EXPANDED)
   oss << networkName << "_Expanded.clu";

 cout << "Results:"  << endl;
 cout << "-->Best code length " << greedy->codeLength/log(2.0) << " with " << Nmod << " modules." << endl;
 cout << "-->Writing state node assignments to \"" << oss.str() << "\"..." << flush;
 outfile.open(oss.str().c_str());
 for(int i=0;i<multiplex.MNnode;i++)
   outfile << reverseMNodeMap[i].first+1 << " " << reverseMNodeMap[i].second+1 << " " << Mclustering[i]+1 << endl;
 outfile.close();
 cout << "done!" << endl;

 int Noverlaps = 0;
 int Nassignments = 0;
 vector<int> assignments = vector<int>(multiplex.Nnode,0);
 for(int i=0;i<Nmod;i++){
  int NNodeAssignments = node[i]->physicalNodes.size();
  Nassignments += NNodeAssignments;
  for(int j=0;j<NNodeAssignments;j++){
    assignments[node[i]->physicalNodes[j].first]++;
    
    if(assignments[node[i]->physicalNodes[j].first] > 1)
      Noverlaps++;
  }
}



   // Order links by size
vector<double> exit(Nmod,0.0);
multimap<double,pair<int,int>,greater<double> > sortedLinks;
for(int i=0;i<Nmod;i++){
  int NoutLinks = node[i]->outLinks.size();
  for(int j=0;j<NoutLinks;j++){
    double linkFlow = node[i]->outLinks[j].second/greedy->beta;
    sortedLinks.insert(make_pair(linkFlow,make_pair(i+1,node[i]->outLinks[j].first+1)));
    exit[i] += linkFlow;
  }
}

   // Order modules by size
double typModSize = 0.0;
double effModSize = 0.0;
int NSig01Mod = 0;
int NSig005Mod = 0;
multimap<double,treeNode,greater<double> > treeMap;
multimap<double,treeNode,greater<double> >::iterator it_tM;
for(int i=0;i<Nmod;i++){
  typModSize += node[i]->size*node[i]->size;
  effModSize -= node[i]->size*log(node[i]->size)/log(2.0);

  if(node[i]->size > 0.01)
    NSig01Mod++;
  if(node[i]->size > 0.005)
    NSig005Mod++;

  int Nmembers = node[i]->physicalNodes.size();
  treeNode tmp_tN;
  it_tM = treeMap.insert(make_pair(node[i]->size,tmp_tN));
  it_tM->second.exit = exit[i];
  for(int j=0;j<Nmembers;j++){
    it_tM->second.members.insert(make_pair(node[i]->physicalNodes[j].second,make_pair(node[i]->physicalNodes[j].first,multiplex.nodeNames[node[i]->physicalNodes[j].first])));
  }
}
effModSize = pow(2.0,effModSize);

   //Print partition in format "module:rank size name"
oss.str("");
if(netMode == MULTIPLEX){
  oss << networkName << "_Multiplex";
  if(dynMode == CLASSICAL)
    oss << "_Classical.tree";
  else
    oss << "_Physical.tree";
}
else if(netMode == AGGREGATED)
  oss << networkName << "_Aggregated.tree";
else if(netMode == EXPANDED)
  oss << networkName << "_Expanded.tree";

cout << "-->Writing physical node assignments to \"" << oss.str() << "\"..." << flush;
outfile.open(oss.str().c_str());
outfile << "# Code length " << greedy->codeLength/log(2.0) << " in " << Nmod << " modules." << endl;
int k = 1;
for(multimap<double,treeNode,greater<double> >::iterator it = treeMap.begin(); it != treeMap.end(); it++){
  string s;
  s.append(to_string(k));
  s.append(":");
  printTree(s,it,&outfile,false);
  k++;
}
outfile.close();
cout << "done!" << endl;

cout << "-->" << Noverlaps << " overlapping module assignments of all physical nodes." << endl;
cout << "-->" << 1.0*Nassignments/multiplex.NConnectedNodes << " assignments per connected physical node." << endl;
cout << "-->Number of modules larger than 1\%: " << NSig01Mod << "." << endl;
cout << "-->Number of modules larger than 0.5\%: " << NSig005Mod << "." << endl;
cout << "-->Typical number of modules: " << 1.0/typModSize << "." << endl;
cout << "-->Typical module size: " << typModSize*multiplex.NConnectedNodes << "." << endl;
cout << "-->Effective number of modules: " << effModSize << "." << endl;
cout << "-->Effective module size: " << multiplex.NConnectedNodes/effModSize << "." << endl;

// // Print results for paper to cerr
// cerr << multiplex.Nnode/effModSize << " " << 1.0*Nassignments/multiplex.Nnode << " " << withinModuleFlow << " " << greedy->codeLength/log(2.0) << endl;
// cerr << switchrate << " " << multiplex.NConnectedNodes/effModSize << " " << 1.0*Nassignments/multiplex.NConnectedNodes << " " << greedy->codeLength/log(2.0) << endl;

for(int i=0;i<Nmod;i++){
  delete node[i];
}
delete [] node;

delete greedy;
delete R;
}

void partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent){

  int Nnode = greedy->Nnode;
  Node **cpy_node = new Node*[Nnode];
  for(int i=0;i<Nnode;i++){
    cpy_node[i] = new Node();
    cpyNode(cpy_node[i],(*node)[i]);
  }

  int iteration = 0;
  double outer_oldCodeLength;
  do{
    outer_oldCodeLength = greedy->codeLength;


    if((iteration > 0) && (iteration % 2 == 0) && (greedy->Nnode > 1)){  // Partition the partition

      if(!silent)
        cout << "-->Iteration " << iteration+1 << ", moving " << flush;

      Node **rpt_node = new Node*[Nnode];
      for(int i=0;i<Nnode;i++){
        rpt_node[i] = new Node();
        cpyNode(rpt_node[i],cpy_node[i]);
      }
      vector<int> subMoveTo(Nnode);
      vector<int> moveTo(Nnode);
      int subModIndex = 0;

      for(int i=0;i<greedy->Nnode;i++){

        int sub_Nnode = (*node)[i]->members.size();

        if(sub_Nnode > 1){
          Node **sub_node = new Node*[sub_Nnode];
          set<int> sub_mem;
          set<int> sub_physicalNodes;
          double totNodeWeights = 0.0;
          for(int j=0;j<sub_Nnode;j++){
            totNodeWeights += cpy_node[(*node)[i]->members[j]]->teleportWeight;
            sub_mem.insert((*node)[i]->members[j]);
          }
          set<int>::iterator it_mem = sub_mem.begin();
          vector<int> sub_renumber = vector<int>(Nnode);
          vector<int> sub_rev_renumber = vector<int>(sub_Nnode);
          for(int j=0;j<sub_Nnode;j++){
            int orig_nr = (*it_mem);
            int orig_NoutLinks = cpy_node[orig_nr]->outLinks.size();
            int orig_NinLinks = cpy_node[orig_nr]->inLinks.size();
            sub_renumber[orig_nr] = j;
            sub_rev_renumber[j] = orig_nr;
            sub_node[j] = new Node(j,cpy_node[orig_nr]->teleportWeight/totNodeWeights);
            sub_node[j]->selfLink =  cpy_node[orig_nr]->selfLink; // Take care of self-link
            for(int k=0;k<orig_NoutLinks;k++){
              int orig_link = cpy_node[orig_nr]->outLinks[k].first;
              int orig_link_newnr = sub_renumber[orig_link];
              double orig_weight = cpy_node[orig_nr]->outLinks[k].second;
              if(orig_link < orig_nr){
                if(sub_mem.find(orig_link) != sub_mem.end()){
                  sub_node[j]->outLinks.push_back(make_pair(orig_link_newnr,orig_weight));
                  sub_node[orig_link_newnr]->inLinks.push_back(make_pair(j,orig_weight));
                }
              }
            }
            for(int k=0;k<orig_NinLinks;k++){
              int orig_link = cpy_node[orig_nr]->inLinks[k].first;
              int orig_link_newnr = sub_renumber[orig_link];
              double orig_weight = cpy_node[orig_nr]->inLinks[k].second;
              if(orig_link < orig_nr){
                if(sub_mem.find(orig_link) != sub_mem.end()){
                  sub_node[j]->inLinks.push_back(make_pair(orig_link_newnr,orig_weight));
                  sub_node[orig_link_newnr]->outLinks.push_back(make_pair(j,orig_weight));
                }
              }
            }

            // Overlap members
            int orig_NphysicalNodes = cpy_node[orig_nr]->physicalNodes.size();
            for(int k=0;k<orig_NphysicalNodes;k++){
              sub_physicalNodes.insert(cpy_node[orig_nr]->physicalNodes[k].first);
              sub_node[j]->physicalNodes.push_back(cpy_node[orig_nr]->physicalNodes[k]);
            }
            //sub_node[j]->physicalNodes = cpy_node[orig_nr]->physicalNodes;

            it_mem++;
          }

          // Renumber overlap members
          //int sub_NphysicalNodes = 0;
          int sub_NphysicalNodes = sub_physicalNodes.size();
          map<int,int> sub_NphysicalNodesMap;
          int physicalNode = 0;
          for(set<int>::iterator it = sub_physicalNodes.begin(); it != sub_physicalNodes.end(); it++){
            sub_NphysicalNodesMap.insert(make_pair((*it),physicalNode));
            physicalNode++;
          }

          //         cout << i << "/" << greedy->Nnode << ": " << sub_Nnode << " " << sub_NphysicalNodes << ":" << flush;

          //          for(set<int>::iterator it = sub_physicalNodes.begin(); it != sub_physicalNodes.end(); it++)
          //            cout << (*it) << " ";

          for(int j=0;j<sub_Nnode;j++){
            int sub_jNphysicalNodes = sub_node[j]->physicalNodes.size();
            //           cout << sub_node[j]->physicalNodes.size() << " " << sub_node[j]->physicalNodes[0].first << endl;
            for(int k=0;k<sub_jNphysicalNodes;k++){
              //             cout << sub_node[j]->physicalNodes[k].first << " --> " << sub_NphysicalNodesMap[sub_node[j]->physicalNodes[k].first] << " " << flush;
              sub_node[j]->physicalNodes[k].first = sub_NphysicalNodesMap[sub_node[j]->physicalNodes[k].first];
              //cout << sub_Nnode << " " << " " << sub_node[j]->physicalNodes[k].first << " " << sub_renumber[sub_node[j]->physicalNodes[k].first] << endl;
            }
            //         cout << endl;
          }

          //sub_NphysicalNodes = sub_Nnode;


          GreedyBase* sub_greedy;
          sub_greedy = new Greedy(R,sub_Nnode,sub_node,sub_NphysicalNodes,greedy->alpha);
          sub_greedy->initiate();
          partition(R,&sub_node,sub_greedy,true);
          for(int j=0;j<sub_greedy->Nnode;j++){
            int Nmembers = sub_node[j]->members.size();
            for(int k=0;k<Nmembers;k++){
              subMoveTo[sub_rev_renumber[sub_node[j]->members[k]]] = subModIndex;
            }
            moveTo[subModIndex] = i;
            subModIndex++;
            delete sub_node[j];
          }

          delete [] sub_node;
          delete sub_greedy;
        }
        else{

          subMoveTo[(*node)[i]->members[0]] = subModIndex;
          moveTo[subModIndex] = i;

          subModIndex++;

        }
      }

      for(int i=0;i<greedy->Nnode;i++)
        delete (*node)[i];
      delete [] (*node);

      greedy->Nnode = Nnode;
      greedy->Nmod = Nnode;
      greedy->Ndanglings = 0;
      greedy->node = rpt_node;
      greedy->calibrate();
      greedy->determMove(subMoveTo);
      greedy->level(node,false);
      greedy->determMove(moveTo);
      (*node) = rpt_node;

      outer_oldCodeLength = greedy->codeLength;

      if(!silent)
        cout << greedy->Nnode << " modules, looping ";

      //abort();

    }
    else if(iteration > 0){

      //if(iteration > 0){

      if(!silent)
        cout << "-->Iteration " << iteration+1 << ", moving " << Nnode << " nodes, looping ";

      Node **rpt_node = new Node*[Nnode];
      for(int i=0;i<Nnode;i++){
        rpt_node[i] = new Node();
        cpyNode(rpt_node[i],cpy_node[i]);
      }

      vector<int>moveTo(Nnode);
      for(int i=0;i<greedy->Nnode;i++){
        int Nmembers = (*node)[i]->members.size();
        for(int j=0;j<Nmembers;j++){
          moveTo[(*node)[i]->members[j]] = i;
        }
      }
      // for(int i=0;i<18;i++)
      // 	moveTo[i] = 0;
      // for(int i=18;i<32;i++)
      // 	moveTo[i] = 1;
      //
      // cout << "Setting nodes in two modules" << endl;

      for(int i=0;i<greedy->Nnode;i++)
        delete (*node)[i];
      delete [] (*node);

      greedy->Nnode = Nnode;
      greedy->Nmod = Nnode;
      greedy->Ndanglings = 0;
      greedy->node = rpt_node;
      greedy->calibrate();
      greedy->determMove(moveTo);
      (*node) = rpt_node;

    }
    else{

      if(!silent)
        cout << "-->Iteration " << iteration+1 << ", moving " << Nnode << " nodes, looping ";

    }


    double oldCodeLength;
    bool aggregated; // To run one more extra time if no node is moved but nodes are aggregated
    do{
      oldCodeLength = greedy->codeLength;
      aggregated = false;
      int moved = 1;
      int Nloops = 0;
      int count = 0;
      while(moved > 0){
        // moved = false;
        double inner_oldCodeLength = greedy->codeLength;
        greedy->move(moved);

        Nloops++;
        count++;
        // cout << setprecision(12) << greedy->codeLength << "-" << inner_oldCodeLength << " " << flush;
        if(fabs(greedy->codeLength - inner_oldCodeLength) < 1.0e-10)
          moved = 0;

        // if(count == 10){
        //   greedy->tune();
        //   count = 0;
        // }
      }

      int NnodeBeforeLevel = greedy->Nnode;

      greedy->level(node,true);
      
      int NnodeAfterLevel = greedy->Nnode;

      if(NnodeBeforeLevel != NnodeAfterLevel)
        aggregated = true;

      if(!silent)
        cout << Nloops << " " << flush;

    } while(oldCodeLength - greedy->codeLength >  1.0e-10 || aggregated);

    iteration++;
    if(!silent)
      cout << "times between mergings to code length " <<  greedy->codeLength/log(2.0) << " with " << greedy->Nmod << " modules." << endl;

  } while(outer_oldCodeLength - greedy->codeLength > 1.0e-10);

  for(int i=0;i<Nnode;i++)
    delete cpy_node[i];
  delete [] cpy_node;

}

void repeated_partition_init(MTRand *R, Node ***node, GreedyBase *greedy, Node ***node_init, bool silent,int Ntrials,vector<int> &nodeassignments,map<int,int> &layerNetworkNodes){

  double shortestCodeLength = 1000.0;
  int Nnode = greedy->Nnode;
  vector<int> cluster(Nnode);
  vector<int> cluster_init(Nnode);

  for(int trial = 0; trial<Ntrials;trial++){

    if(!silent)
      cout << "-->Attempt " << trial+1 << "/" << Ntrials << endl;

    Node **cpy_node = new Node*[Nnode];
    for(int i=0;i<Nnode;i++){
      cpy_node[i] = new Node();
      cpyNode(cpy_node[i],(*node)[i]);
    }

    greedy->Nnode = Nnode;
    greedy->Nmod = Nnode;
    greedy->Ndanglings = 0;
    greedy->node = cpy_node;
    greedy->calibrate();

    if(nodeassignments.size() == Nnode){
      greedy->determMove(nodeassignments);
    }
    cout << "-->Intial code length " << greedy->codeLength/log(2.0) << "." << endl;

    partition(R,&cpy_node,greedy,silent);

    if(greedy->codeLength < shortestCodeLength){

      shortestCodeLength = greedy->codeLength;

      // Store best partition
      for(int i=0;i<greedy->Nnode;i++){
        for(vector<int>::iterator mem = cpy_node[i]->members.begin(); mem != cpy_node[i]->members.end(); mem++){
          cluster[(*mem)] = i;
        }
      }
      
    }

    for(int i=0;i<greedy->Nnode;i++){
      delete cpy_node[i];
    }
    delete [] cpy_node;

    if(!silent)
      cout << "-->Alternative initialization..." << flush;

    int minNodeNr = 0;
    int maxClusterNr = 0;

    for(map<int,int>::iterator it = layerNetworkNodes.begin(); it != layerNetworkNodes.end(); it++){

      int maxNodeNr = it->second;
      int NnodeLayer = maxNodeNr - minNodeNr + 1;
      int maxClusterNrTmp = 0;

      Node **cpy_node_init = new Node*[NnodeLayer];
      double totNodeWeights = 0.0;
      for(int i=0;i<NnodeLayer;i++)
        totNodeWeights += (*node_init)[minNodeNr+i]->teleportWeight;

      for(int i=0;i<NnodeLayer;i++){
        cpy_node_init[i] = new Node();
        cpyNode(cpy_node_init[i],(*node_init)[minNodeNr+i]);
        // Renumber
        if(minNodeNr > 0){
          cpy_node_init[i]->index -= minNodeNr;
          cpy_node_init[i]->teleportWeight /= totNodeWeights;
          cpy_node_init[i]->members[0] -= minNodeNr;
          cpy_node_init[i]->physicalNodes[0].first -= minNodeNr;
          int NoutLinks = cpy_node_init[i]->outLinks.size();
          for(int j=0;j<NoutLinks;j++)
            cpy_node_init[i]->outLinks[j].first -= minNodeNr;
          int NinLinks = cpy_node_init[i]->inLinks.size();
          for(int j=0;j<NinLinks;j++)
            cpy_node_init[i]->inLinks[j].first -= minNodeNr;
        }
      }

      // Initiation
      GreedyBase* greedy_init;
      greedy_init = new Greedy(R,NnodeLayer,cpy_node_init,NnodeLayer,greedy->alpha);
      
      greedy_init->initiate();

      partition(R,&cpy_node_init,greedy_init,true);

      for(int i=0;i<greedy_init->Nnode;i++){
        for(vector<int>::iterator mem = cpy_node_init[i]->members.begin(); mem != cpy_node_init[i]->members.end(); mem++){
          int clusterNr = maxClusterNr + i;
          maxClusterNrTmp = max(maxClusterNrTmp,clusterNr);

          cluster_init[minNodeNr + (*mem)] = clusterNr;
        }
      }

      for(int i=0;i<greedy_init->Nnode;i++){
        delete cpy_node_init[i];
      }
      delete [] cpy_node_init;
      delete greedy_init;

      maxClusterNr = maxClusterNrTmp + 1;
      minNodeNr = maxNodeNr + 1;

    }

    cout << "with " << maxClusterNr << " module(s)." << endl;

    cpy_node = new Node*[Nnode];
    for(int i=0;i<Nnode;i++){
      cpy_node[i] = new Node();
      cpyNode(cpy_node[i],(*node)[i]);
    }

    greedy->Nnode = Nnode;
    greedy->Nmod = Nnode;
    greedy->Ndanglings = 0;
    greedy->node = cpy_node;
    greedy->calibrate();

    greedy->determMove(cluster_init);

    cout << "-->Intial code length " << greedy->codeLength/log(2.0) << "." << endl;

    partition(R,&cpy_node,greedy,silent);

    if(greedy->codeLength < shortestCodeLength){

      shortestCodeLength = greedy->codeLength;

      // Store best partition
      for(int i=0;i<greedy->Nnode;i++){
        for(vector<int>::iterator mem = cpy_node[i]->members.begin(); mem != cpy_node[i]->members.end(); mem++){
          cluster[(*mem)] = i;
        }
      }
      
    }

    for(int i=0;i<greedy->Nnode;i++){
      delete cpy_node[i];
    }
    delete [] cpy_node;

  }

  // Commit best partition
  greedy->Nnode = Nnode;
  greedy->Nmod = Nnode;
  greedy->Ndanglings = 0;
  greedy->node = (*node);
  greedy->calibrate();
  greedy->determMove(cluster);
  greedy->level(node,true);

}

void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent,int Ntrials,vector<int> &nodeassignments){

  double shortestCodeLength = 1000.0;
  int Nnode = greedy->Nnode;
  vector<int> cluster(Nnode);

  for(int trial = 0; trial<Ntrials;trial++){

    if(!silent)
      cout << "-->Attempt " << trial+1 << "/" << Ntrials << endl;

    Node **cpy_node = new Node*[Nnode];
    for(int i=0;i<Nnode;i++){
      cpy_node[i] = new Node();
      cpyNode(cpy_node[i],(*node)[i]);
    }

    greedy->Nnode = Nnode;
    greedy->Nmod = Nnode;
    greedy->Ndanglings = 0;
    greedy->node = cpy_node;
    greedy->calibrate();

    if(nodeassignments.size() == Nnode){
      greedy->determMove(nodeassignments);
    }
    cout << "-->Intial code length " << greedy->codeLength/log(2.0) << "." << endl;

    partition(R,&cpy_node,greedy,silent);

    if(greedy->codeLength < shortestCodeLength){

      shortestCodeLength = greedy->codeLength;

      // Store best partition
      for(int i=0;i<greedy->Nnode;i++){
        for(vector<int>::iterator mem = cpy_node[i]->members.begin(); mem != cpy_node[i]->members.end(); mem++){
          cluster[(*mem)] = i;
        }
      }
      
    }

    for(int i=0;i<greedy->Nnode;i++){
      delete cpy_node[i];
    }
    delete [] cpy_node;

  }

  // Commit best partition
  greedy->Nnode = Nnode;
  greedy->Nmod = Nnode;
  greedy->Ndanglings = 0;
  greedy->node = (*node);
  greedy->calibrate();
  greedy->determMove(cluster);
  greedy->level(node,true);

}


void printTree(string s,multimap<double,treeNode,greater<double> >::iterator it_tM,ofstream *outfile,bool flip){

  multimap<double,treeNode,greater<double> >::iterator it;
  if(it_tM->second.nextLevel.size() > 0){
    int i=1;
    for(it = it_tM->second.nextLevel.begin(); it != it_tM->second.nextLevel.end(); it++){
      string cpy_s(s + to_string(i) + ":");
      printTree(cpy_s,it,outfile,flip);
      i++;
    }
  }
  else{
    int i = 1;
    for(multimap<double,pair<int,string>,greater<double> >::iterator mem = it_tM->second.members.begin(); mem != it_tM->second.members.end(); mem++){
      if(flip){
        string cpy_s(s + to_string(i) + " \"" + mem->second.second + "\" " + to_string(mem->first));
        (*outfile) << cpy_s << endl;
      }
      else{
        string cpy_s(s + to_string(i) + " " + to_string(mem->first) + " \"" + mem->second.second + "\"");
        (*outfile) << cpy_s << endl;
      }
      i++;
    }
  }
}

