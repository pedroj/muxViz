#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "MersenneTwister.h"
#include "GreedyBase.h"
#include "Greedy.h"
#include "Node.h"
#define PI 3.14159265
using namespace std;

enum NetworkMode { MULTIPLEX, AGGREGATED, EXPANDED, INIT, INITAGG };
enum DynamicMode { CLASSICAL, PHYSICAL };
enum SwitchMode { UNIFORM, PROPORTIONAL, LINKS};


unsigned stou(char *s);

vector<string> tokenize(const string& str,string& delimiters){

	vector<string> tokens;

      // skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while(string::npos != pos || string::npos != lastPos){

        // found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));

        // skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);

        // find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}

	return tokens;

}

class Multiplex{

public:

	Multiplex(vector<string> networkLayerFiles,string layerNetworkFile,string nodenamesFile,string nodeassignmentsFile);
	vector<string> networklayernames;
	string layernetworkname;
	string nodename;
	string assignmentname;
	int Nnode;
	int NConnectedNodes;
	int Nlinks;
	vector<string> nodeNames;
	map<pair<int,int>,double> Links;

	int Nlayers;
	int MNnode;
	int MNlinks;
	double totMNodeWeights;
	vector<int> nodeassignments;
	map<pair<int,int>,int> MnodeMap;
	vector<double> MnodeWeights;
	map<pair<pair<int,int>,pair<int,int> >,double> MLinks;
	map<pair<int,pair<int,int> >,double> interLinks;
};

Multiplex::Multiplex(vector<string> networkLayerFiles,string layerNetworkFile,string nodenamesFile,string nodeassignmentsFile){
	networklayernames = networkLayerFiles;
	layernetworkname = layerNetworkFile;
	nodename = nodenamesFile;
	assignmentname = nodeassignmentsFile;
}

class treeNode{
public:
	double exit;
	multimap<double,pair<int,string>,greater<double> > members;
	multimap<double,treeNode,greater<double> > nextLevel;
};

template <class T>
inline std::string to_string (const T& t){
	std::stringstream ss;
	ss << t;
	return ss.str();
}

void cpyNode(Node *newNode,Node *oldNode){

	newNode->index = oldNode->index;
	newNode->enter = oldNode->enter;
	newNode->exit = oldNode->exit;
	newNode->size = oldNode->size;
	newNode->teleportWeight = oldNode->teleportWeight;

	newNode->members = oldNode->members;

	newNode->physicalNodes = oldNode->physicalNodes;

	newNode->selfLink = oldNode->selfLink;

	newNode->outLinks = oldNode->outLinks;

	newNode->inLinks = oldNode->inLinks;

}


void loadMnetwork(Multiplex &multiplex,NetworkMode netMode,bool selflink,SwitchMode switchMode,double switchrate){

	if(!(netMode == INIT || netMode == INITAGG))
		cout << "Inputs:" << endl;

	set<int> ConnectedNodes;

	string line;
	string buf;
	map<pair<int,int>,double > Mnodes;
	multiplex.Nnode = 0;
	multiplex.Nlayers = multiplex.networklayernames.size();

    // Read intra-layer networks
	for(int layer = 0;layer < multiplex.Nlayers; layer++){

		cout << "-->Reading network layer " << layer+1 << ": " << multiplex.networklayernames[layer] << "..." << flush;

		int activeLayer = layer;
		if(netMode == AGGREGATED || netMode == INITAGG){
			activeLayer = 0;
			cout << "collapsing layer..." << flush;
		}

		ifstream net(multiplex.networklayernames[layer].c_str());
		if(!net){
			cout << "failed to open \"" << multiplex.networklayernames[layer] << "\" exiting…" << endl;
			exit(-1);
		}
		istringstream ss;
		while(getline(net,line) != 0){

			ss.clear();
			ss.str(line);
			ss >> buf;
			if(buf[0] != '#'){

				int n1 = atoi(buf.c_str()) - 1;
				ss >> buf;
				int n2 = atoi(buf.c_str()) - 1;
				buf.clear();
				ss >> buf;
				double linkWeight;
     if( buf.empty() ) // If no information
     	linkWeight = 1.0;
     else
     	linkWeight = atof(buf.c_str());

     multiplex.Nnode = max(multiplex.Nnode,n1+1);
     multiplex.Nnode = max(multiplex.Nnode,n2+1);

     multiplex.MLinks[make_pair(make_pair(activeLayer,n1),make_pair(activeLayer,n2))] += linkWeight;
     multiplex.Links[make_pair(n1,n2)] += linkWeight;
     if(selflink || (n1 != n2) ){
     	Mnodes[make_pair(activeLayer,n1)] += linkWeight;
     }
     Mnodes[make_pair(activeLayer,n2)] += 0.0;
   
     if(netMode == INIT || netMode == INITAGG){ // Add link in other direction for undirected network
     	multiplex.MLinks[make_pair(make_pair(activeLayer,n2),make_pair(activeLayer,n1))] += linkWeight;
     	multiplex.Links[make_pair(n2,n1)] += linkWeight;
     	if(selflink || (n1 != n2) ){
     		Mnodes[make_pair(activeLayer,n2)] += linkWeight;
     	}
     	Mnodes[make_pair(activeLayer,n1)] += 0.0;
     }

 }
}

net.close();
cout << "done!" << endl;

}

  // Read intra- and inter-layer network
if(multiplex.Nlayers == 0)
	cout << "-->Reading multiplex network: " << multiplex.layernetworkname << "..." << flush;
else
	cout << "-->Reading layer network: " << multiplex.layernetworkname << "..." << flush;
ifstream net(multiplex.layernetworkname.c_str());
if(!net){
	cout << "failed to open \"" << multiplex.layernetworkname << "\" exiting…" << endl;
	exit(-1);
}
bool intra = true;
istringstream ss;
while(getline(net,line) != 0){

	ss.clear();
	ss.str(line);
	ss >> buf;
	if(buf[0] != '#'){
		if(buf == "*Intra" || buf == "*intra"){
			intra = true;
		}
		else if(buf == "*Inter" || buf == "*inter"){
			intra = false;
		}
		else{

			int layer1,layer2,n1,n2;

			if(intra){
				layer1 = atoi(buf.c_str()) - 1;
				layer2 = layer1;
				ss >> buf;
				n1 = atoi(buf.c_str()) - 1;
				ss >> buf;
				n2 = atoi(buf.c_str()) - 1;
			}
			else{//inter
				n1 = atoi(buf.c_str()) - 1;
				n2 = n1;
				ss >> buf;
				layer1 = atoi(buf.c_str()) - 1;
				ss >> buf;
				layer2 = atoi(buf.c_str()) - 1;
			}
			buf.clear();
			ss >> buf;
			double linkWeight;
            if( buf.empty() )// If no information
            	linkWeight = 1.0;
            else
            	linkWeight = atof(buf.c_str());

            multiplex.Nlayers = max(multiplex.Nlayers,layer1+1);
            multiplex.Nlayers = max(multiplex.Nlayers,layer2+1);

            // Collapse layers but don't create self-links from inter-layer links
            if(intra && (netMode == AGGREGATED || netMode == INITAGG)){
            	layer1 = 0;
            	layer2 = 0;
            }

            ConnectedNodes.insert(n1);
            ConnectedNodes.insert(n2);

            multiplex.Nnode = max(multiplex.Nnode,n1+1);
            multiplex.Nnode = max(multiplex.Nnode,n2+1);

            // if(netMode == INIT){ // Add link in both directions for undirected network 
            //     multiplex.MLinks[make_pair(make_pair(layer1,n1),make_pair(layer2,n2))] += linkWeight;
            //     multiplex.Links[make_pair(n1,n2)] += linkWeight;
            //     Mnodes[make_pair(layer1,n1)] += linkWeight;
            //     Mnodes[make_pair(layer2,n2)] += 0.0;
            //     multiplex.MLinks[make_pair(make_pair(layer2,n2),make_pair(layer1,n1))] += linkWeight;
            //     multiplex.Links[make_pair(n2,n1)] += linkWeight;
            //     Mnodes[make_pair(layer2,n2)] += linkWeight;
            //     Mnodes[make_pair(layer1,n1)] += 0.0;
            // }
            if(intra){
                 if(netMode == INIT || netMode == INITAGG){ // Add link in both directions for undirected network if not a between-layer link
                 	if(layer1 == layer2){
                 		multiplex.MLinks[make_pair(make_pair(layer1,n1),make_pair(layer2,n2))] += linkWeight;
                 		multiplex.Links[make_pair(n1,n2)] += linkWeight;
                 		if(selflink || (n1 != n2) ){
                 			Mnodes[make_pair(layer1,n1)] += linkWeight;
                 		}
                 		Mnodes[make_pair(layer2,n2)] += 0.0;

                 		multiplex.MLinks[make_pair(make_pair(layer2,n2),make_pair(layer1,n1))] += linkWeight;
                 		multiplex.Links[make_pair(n2,n1)] += linkWeight;
                 		if(selflink || (n1 != n2) ){
                 			Mnodes[make_pair(layer2,n2)] += linkWeight;
                 		}
                 		Mnodes[make_pair(layer1,n1)] += 0.0;

                 	}
                 	else{
                 		Mnodes[make_pair(layer1,n1)] += 0.0;
                 		Mnodes[make_pair(layer2,n2)] += 0.0;
                 	}
                 }
                 else{

                 	multiplex.MLinks[make_pair(make_pair(layer1,n1),make_pair(layer2,n2))] += linkWeight;
                 	multiplex.Links[make_pair(n1,n2)] += linkWeight;
                 	if(selflink || (n1 != n2) ){
                 	// if(selflink || !((n1 == n2) && (layer1 == layer2)) ){
                 		Mnodes[make_pair(layer1,n1)] += linkWeight;
                 	}
                 	Mnodes[make_pair(layer2,n2)] += 0.0;

                 }
             }
             else if(switchMode == LINKS){

             	multiplex.interLinks[make_pair(n1,make_pair(layer1,layer2))] += linkWeight;

             }        

         }
     }
 }
 net.close();
 cout << "done!" << endl;

 multiplex.NConnectedNodes = ConnectedNodes.size();

 if(netMode == AGGREGATED || netMode == INITAGG || netMode == INIT){
 	cout << "-->Ignoring inter-layer link weights." << endl;
 }
 else if(switchMode == LINKS){

	//Rescale inter-layer link weights for one-step dynamics.
 	cout << "-->Rescaling " << multiplex.interLinks.size() << " inter-layer link weights..." << flush;
 	for(map<pair<int,pair<int,int> >,double>::iterator it = multiplex.interLinks.begin(); it != multiplex.interLinks.end(); it++){
 		int n1 = it->first.first;
 		int layer1 = it->first.second.first;
 		int layer2 = it->first.second.second;
 		if(layer1 != layer2){

 			map<pair<int,pair<int,int> >,double>::iterator selfLayer_it = multiplex.interLinks.find(make_pair(n1,make_pair(layer1,layer1)));
 			if(selfLayer_it != multiplex.interLinks.end()){

 				map<pair<int,int>,double>::iterator strength_it = Mnodes.find(make_pair(layer1,n1));
 				if(strength_it != Mnodes.end()){
 					double linkWeight = strength_it->second/selfLayer_it->second*it->second;
 					multiplex.MLinks[make_pair(make_pair(layer1,n1),make_pair(layer2,n1))] += linkWeight;
 					Mnodes[make_pair(layer1,n1)] += 0;
 					Mnodes[make_pair(layer2,n1)] += 0;
 				}
 				else{
 					double linkWeight = it->second;
 					multiplex.MLinks[make_pair(make_pair(layer1,n1),make_pair(layer2,n1))] += linkWeight;
 					Mnodes[make_pair(layer1,n1)] += 0;
 					Mnodes[make_pair(layer2,n1)] += 0;
 				}
 			}
 			else{

 				cout << "The self-layer link for node " << n1+1 << " at layer " << layer1+1 << " is not provided, using the node strength instead." << endl;
 				double linkWeight = it->second;
 				multiplex.MLinks[make_pair(make_pair(layer1,n1),make_pair(layer2,n1))] += linkWeight;
 				Mnodes[make_pair(layer1,n1)] += 0;
 				Mnodes[make_pair(layer2,n1)] += 0;

 			}
 		}
 	}
 	cout << "done!" << endl;
 }
 else if(multiplex.Nlayers > 1){

 	//Generate inter-layer link weights from parameters.
 	cout << "-->Generating parametrized inter-layer link weights..." << flush;

 	if(switchMode == PROPORTIONAL){

 		vector<double> physicalNodeStrength(multiplex.Nnode,0.0);

 		for(map<pair<int,int>,double >::iterator it = Mnodes.begin(); it != Mnodes.end(); it++){
 			physicalNodeStrength[it->first.second] += it->second;
 		}

 		for(map<pair<int,int>,double >::iterator it = Mnodes.begin(); it != Mnodes.end(); it++){

 			if(it->second > 0.0){

 				for(int layer = 0; layer < multiplex.Nlayers; layer++){

 					if(layer != it->first.first){

 						map<pair<int,int>,double >::iterator layer_it = Mnodes.find(make_pair(layer,it->first.second));
 						if(layer_it != Mnodes.end()){
 							if(layer_it->second > 0.0){
 								double linkWeight = it->second*switchrate*layer_it->second/(switchrate*it->second + (1.0-switchrate)*physicalNodeStrength[it->first.second]);  //Switching includes current layer
 								// double linkWeight = it->second*switchrate/(1.0-switchrate)*layer_it->second/(physicalNodeStrength[it->first.second]-it->second); //Switching excludes current layer
 								// cout << it->first.first << " " << it->first.second << " " << linkWeight << endl;
 								multiplex.MLinks[make_pair(it->first,layer_it->first)] += linkWeight;
 							}
 						}
 					}
 				}
 			}
 		}
 	}
 	else{

 		vector<int> physicalNodeStrength(multiplex.Nnode,0);

 		for(map<pair<int,int>,double >::iterator it = Mnodes.begin(); it != Mnodes.end(); it++){
 			if(it->second > 0)
 				physicalNodeStrength[it->first.second]++;
 		}

 		for(map<pair<int,int>,double >::iterator it = Mnodes.begin(); it != Mnodes.end(); it++){

 			if(it->second > 0.0){

 				for(int layer = 0; layer < multiplex.Nlayers; layer++){

 					if(layer != it->first.first){

 						map<pair<int,int>,double >::iterator layer_it = Mnodes.find(make_pair(layer,it->first.second));
 						if(layer_it != Mnodes.end()){
 							if(layer_it->second > 0.0){
 								double linkWeight = it->second*switchrate*layer_it->second/(switchrate*it->second + (1.0-switchrate)*physicalNodeStrength[it->first.second]); //Switching includes current layer
 								// double linkWeight = it->second*switchrate/(1.0-switchrate)/(physicalNodeStrength[it->first.second]-1); //Switching excludes current layer
 								multiplex.MLinks[make_pair(it->first,layer_it->first)] += linkWeight;
 							}
 						}
 					}
 				}
 			}
 		}
 	}
 	cout << "done!" << endl;
 	
 }

 // for(map<pair<pair<int,int>,pair<int,int> >,double>::iterator it = multiplex.MLinks.begin(); it != multiplex.MLinks.end(); it++)
 // 	cout << it->first.first.first+1 << " " << it->first.first.second+1 << " " << it->first.second.first+1 << " " << it->first.second.second+1 << " " << it->second << endl;

  // Read node names
 if(multiplex.nodename == "nonames"){
 	multiplex.nodeNames = vector<string>(multiplex.Nnode);
 	for(int i=0;i<multiplex.Nnode;i++){
 		multiplex.nodeNames[i] = "Node " + to_string(i+1);
 	}
 }
 else{
 	cout << "-->Reading node names from: " << multiplex.nodename << "..." << flush;
 	ifstream netnames(multiplex.nodename.c_str());
 	if(!netnames){
 		cout << "failed to open \"" << multiplex.nodename << "\" exiting…" << endl;
 		exit(-1);
 	}

 	while(getline(netnames,line) != 0){
 		int nameStart = line.find_first_of("\"");
 		int nameEnd = line.find_last_of("\"");
 		if(nameStart < nameEnd){
 			multiplex.nodeNames.push_back(string(line.begin() + nameStart + 1,line.begin() + nameEnd));
 		}
 		else{
 			multiplex.nodeNames.push_back(line);
 		}
 		multiplex.Nnode = multiplex.nodeNames.size();
 	}
 	cout << "done!" << endl;
 }


 multiplex.Nlinks = multiplex.Links.size();
 multiplex.MNnode = Mnodes.size();
 multiplex.MNlinks = multiplex.MLinks.size();

 int MnodeNr = 0;
 for(map<pair<int,int>,double >::iterator it = Mnodes.begin(); it != Mnodes.end(); it++){
 	multiplex.MnodeMap[it->first] = MnodeNr;
 	MnodeNr++;
 }

 multiplex.MnodeWeights = vector<double>(multiplex.MNnode);
 multiplex.totMNodeWeights = 0.0;
 MnodeNr = 0;
 for(map<pair<int,int>,double >::iterator it = Mnodes.begin(); it != Mnodes.end(); it++){
 	multiplex.MnodeWeights[MnodeNr] += it->second;
 	multiplex.totMNodeWeights += it->second;
 	MnodeNr++;
 }

 if(netMode == AGGREGATED)
 	cout << "-->Found " << multiplex.Nlinks << " links between " << multiplex.Nnode << " nodes in " << multiplex.Nlayers << " layer(s), collapsed them into a single layer, and generated " << multiplex.MNnode << " nodes and " << multiplex.MNlinks << " links." << endl;
 else if(netMode == EXPANDED)
 	cout << "-->Found " << multiplex.Nlinks << " links between " << multiplex.Nnode << " nodes in " << multiplex.Nlayers << " layer(s), expanded them into a single layer, and generated " << multiplex.MNnode << " fake physical nodes and " << multiplex.MNlinks << " links." << endl;
 else if(netMode == INIT)
 	cout << "-->Found " << multiplex.Nlinks << " links between " << multiplex.Nnode << " nodes in " << multiplex.Nlayers << " disconnected layer(s), expanded them into a single layer of disconnected networks, and generated " << multiplex.MNnode << " fake physical nodes and " << multiplex.MNlinks/2 << " undirected links." << endl;
// else if(netMode == INIT)
//   cout << "-->Found " << multiplex.Nlinks << " links between " << multiplex.Nnode << " nodes in " << multiplex.Nlayers << " layer(s), expanded them into a single layer, and generated " << multiplex.MNnode << " fake physical nodes and " << multiplex.MNlinks/2 << " undirected links." << endl;
 else if(netMode == INITAGG)
 	cout << "-->Found " << multiplex.Nlinks << " links between " << multiplex.Nnode << " nodes in " << multiplex.Nlayers << " layer(s), collapsed them into a single layer, and generated " << multiplex.MNnode << " nodes and " << multiplex.MNlinks/2 << " undirected links." << endl;
 else
 	cout << "-->Found " << multiplex.Nlinks << " links between " << multiplex.Nnode << " nodes in " << multiplex.Nlayers << " layer(s), and generated " << multiplex.MNnode << " state nodes and " << multiplex.MNlinks << " state links." << endl;

  // Read node assignments
 if(multiplex.assignmentname != "noassignments"){

 	cout << "-->Reading node assignments from: " << multiplex.assignmentname << "..." << endl;
 	ifstream nodeass(multiplex.assignmentname.c_str());
 	if(!nodeass){
 		cout << "-->Failed to open \"" << multiplex.assignmentname << "\" exiting..." << endl;
 		exit(-1);
 	}

 	multiplex.nodeassignments = vector<int>(multiplex.MNnode);
 	set<pair<int,int> > checkOffAssignments;
 	for(map<pair<int,int>,int>::iterator it = multiplex.MnodeMap.begin(); it != multiplex.MnodeMap.end(); it++)
 		checkOffAssignments.insert(checkOffAssignments.end(),it->first);

 	map<int,int> physicalAssignment;

 	while(getline(nodeass,line) != 0){
 		ss.clear();
 		ss.str(line);
 		ss >> buf;
 		if(buf[0] != '#'){
 			int layer = atoi(buf.c_str()) - 1;
 			ss >> buf;
 			int node = atoi(buf.c_str()) - 1;
 			ss >> buf;
 			int assignment = atoi(buf.c_str()) - 1;

 			physicalAssignment[node] = assignment;

 			map<pair<int,int>,int>::iterator it = multiplex.MnodeMap.find(make_pair(layer,node));
 			if(it != multiplex.MnodeMap.end()){
 				multiplex.nodeassignments[it->second] = assignment;
 				checkOffAssignments.erase(it->first);
 			}
 			else{
 				cout  << "-->Warning, state node " << layer+1 << " " << node+1 << " never used." << endl;
 			}
 		}
 	}

 	int Nmissing = checkOffAssignments.size();
 	int Nfound = 0;
 	if(Nmissing > 0){
 		cout << "-->Warning, missing assignments for state nodes:" << endl;
 		for(set<pair<int,int> >::iterator it = checkOffAssignments.begin(); it != checkOffAssignments.end(); it++){
 			map<int,int>::iterator it_p = physicalAssignment.find(it->second);
 			if(it_p != physicalAssignment.end()){
 				map<pair<int,int>,int>::iterator it_a = multiplex.MnodeMap.find(*(it));
 				if(it_a != multiplex.MnodeMap.end()){
 					multiplex.nodeassignments[it_a->second] = it_p->second;
 					Nfound++;
 					cout << "-->Warning, missing assignment for state node " << it->first+1 << " " << it->second+1 << ", using physical node assignment " << it_p->second+1 << endl;	
 				}
 				else{
 					cout << "-->Warning, missing assignment for state node " << it->first+1 << " " << it->second+1 << endl;
 				}
 			}
 			else{
 				cout << "-->Warning, missing assignment for state node " << it->first+1 << " " << it->second+1 << endl;
 			}
 		}
 		
 		if(Nfound != Nmissing){
 			cout << "-->Could not recover assignments for all state nodes from physical nodes, please provide complete assignment data. Exiting..." << endl;
 			exit(-1);
 		}
 	}
 	
 	cout << "-->Completed reading node assignments!" << endl;

 }

	}





