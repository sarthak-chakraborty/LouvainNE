// File: main_louvain.cpp
// -- community detection, sample main file
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// And based on the article
// Copyright (C) 2013 R. Campigotto, P. Conde Céspedes, J.-L. Guillaume
//
// This file is part of Louvain algorithm.
// 
// Louvain algorithm is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Louvain algorithm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time     : July 2013
//-----------------------------------------------------------------------------
// see README.txt for more details


#include "graph_binary.h"
#include "louvain.h"
#include <unistd.h>
#include <set>
#include <string.h>
#include <fstream>
#include <unordered_map>
#include "math.h"

#include "modularity.h"
#include "zahn.h"
#include "owzad.h"
#include "goldberg.h"
#include "condora.h"
#include "devind.h"
#include "devuni.h"
#include "dp.h"
#include "shimalik.h"
#include "balmod.h"

#include "snap-adv/n2v.h"


using namespace std;


char *filename = NULL;
char *filename_w = NULL;
char *filename_part = NULL;
int type = UNWEIGHTED;

int nb_pass = 0;
long double precision = 0.000001L;
int display_level = -2;

long double exponent = 0.05L;
unsigned int dimension = 10;

unsigned short id_qual = 0;

long double alpha = 0.05L;
int kmin = 1;

long double sum_se = 0.0L;
long double sum_sq = 0.0L;

long double max_w = 1.0L;

Quality *q;

bool verbose = false;

double ParamP = 1.;
double ParamQ = 1.;
int walklen = 80;
int numwalks = 10;
int winsize = 10;
int iter = 1;
bool outputwalks = false;


void
usage(char *prog_name, const char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " input_file [-q id_qual] [-c alpha] [-k min] [-w weight_file] [-p part_file] [-e epsilon] [-l display_level] [-v] [-h]" << endl << endl;
  cerr << "input_file: file containing the graph to decompose in communities" << endl;

  cerr << "-a exponent\ta the exponent parameter for the generation of random vectors" << endl;
  cerr << "-n dimension\ta the dimension of the vector" << endl;

  cerr << "-q id\tthe quality function used to compute partition of the graph (modularity is chosen by default):" << endl << endl;

  cerr << "\tid = 0\t -> the classical Newman-Girvan criterion (also called \"Modularity\")" << endl;
  cerr << "\tid = 1\t -> the Zahn-Condorcet criterion" << endl;
  cerr << "\tid = 2\t -> the Owsinski-Zadrozny criterion (you should specify the value of the parameter with option -c)" << endl;
  cerr << "\tid = 3\t -> the Goldberg Density criterion" << endl;
  cerr << "\tid = 4\t -> the A-weighted Condorcet criterion" << endl;
  cerr << "\tid = 5\t -> the Deviation to Indetermination criterion" << endl;
  cerr << "\tid = 6\t -> the Deviation to Uniformity criterion" << endl;
  cerr << "\tid = 7\t -> the Profile Difference criterion" << endl;
  cerr << "\tid = 8\t -> the Shi-Malik criterion (you should specify the value of kappa_min with option -k)" << endl;
  cerr << "\tid = 9\t -> the Balanced Modularity criterion" << endl;

  cerr << endl;

  cerr << "-c al\tthe parameter for the Owsinski-Zadrozny quality function (between 0.0 and 1.0: 0.5 is chosen by default)" << endl;
  cerr << "-k min\tthe kappa_min value (for Shi-Malik quality function) (it must be > 0: 1 is chosen by default)" << endl;

  cerr << endl;

  cerr << "-w file\tread the graph as a weighted one (weights are set to 1 otherwise)" << endl;
  cerr << "-p file\tstart the computation with a given partition instead of the trivial partition" << endl;
  cerr << "\tfile must contain lines \"node community\"" << endl;
  cerr << "-e eps\ta given pass stops when the quality is increased by less than epsilon" << endl;
  cerr << "-l k\tdisplays the graph of level k rather than the hierachical structure" << endl;
  cerr << "\tif k=-1 then displays the hierarchical structure rather than the graph at a given level" << endl;
  cerr << "-v\tverbose mode: gives computation time, information about the hierarchy and quality" << endl;
  cerr << "-h\tshow this usage message" << endl;

  exit(0);
}

void
parse_args(int argc, char **argv) {
  if (argc<2)
    usage(argv[0], "Bad arguments number\n");

  for (int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
      case 'w':
        type = WEIGHTED;
        filename_w = argv[i+1];
  i++;
  break;
      case 'q':
  id_qual = (unsigned short)atoi(argv[i+1]);
  i++;
  break;
      case 'c':
  alpha = atof(argv[i+1]);
  i++;
  break;
      case 'k':
  kmin = atoi(argv[i+1]);
  i++;
  break;
      case 'p':
        filename_part = argv[i+1];
  i++;
  break;
      case 'e':
  precision = atof(argv[i+1]);
  i++;
  break;
      case 'a':
        exponent = atof(argv[i+1]);
  i++;
  break;
      case 'n':
  dimension = atoi(argv[i+1]);
  i++;
  break;
      case 'l':
  display_level = atoi(argv[i+1]);
  i++;
  break;
      case 'v':
  verbose = true;
  break;
      case 'h':
  usage(argv[0], "");
  break;
      default:
  usage(argv[0], "Unknown option\n");
      }
    } else {
      if (filename==NULL)
        filename = argv[i];
      else
        usage(argv[0], "More than one filename\n");
    }
  }
  if (filename == NULL)
    usage(argv[0], "No input file has been provided\n");
}

void
display_time(const char *str) {
  time_t rawtime;
  time ( &rawtime );
  cerr << str << ": " << ctime (&rawtime);
}

void
init_quality(Graph *g, unsigned short nbc) {

  if (nbc > 0)
    delete q;

  switch (id_qual) {
  case 0:
    q = new Modularity(*g);
    break;
  case 1:
    if (nbc == 0)
      max_w = g->max_weight();
    q = new Zahn(*g, max_w);
    break;
  case 2:
    if (nbc == 0)
      max_w = g->max_weight();
    if (alpha <= 0. || alpha >= 1.0)
      alpha = 0.5;
    q = new OwZad(*g, alpha, max_w);
    break;
  case 3:
    if (nbc == 0)
      max_w = g->max_weight();
    q = new Goldberg(*g, max_w);
    break;
  case 4:
    if (nbc == 0) {
      g->add_selfloops();
      sum_se = CondorA::graph_weighting(g);
    }
    q = new CondorA(*g, sum_se);
    break;
  case 5:
    q = new DevInd(*g);
    break;
  case 6:
    q = new DevUni(*g);
    break;
  case 7:
    if (nbc == 0) {
      max_w = g->max_weight();
      sum_sq = DP::graph_weighting(g);
    }
    q = new DP(*g, sum_sq, max_w);
    break;
  case 8:
    if (kmin < 1)
      kmin = 1;
    q = new ShiMalik(*g, kmin);
    break;
  case 9:
    if (nbc == 0)
      max_w = g->max_weight();
    q = new BalMod(*g, max_w);
    break;
  default:
    q = new Modularity(*g);
    break;
  }
}

// generate a random vector
vector<double>
random_vector() {
  vector<double> res(dimension);
  for (unsigned int i = 0; i < dimension; i++)
    res[i] = (double)rand() / RAND_MAX;
  return res;
}



/* NODE2VEC EMBEDDING */
// given a vector, multiply it by alpha and add another random vector to it
void
multiply_and_add(Graph g, vector<int> pnodes, vector<int> cnodes, vector<vector<double>>& embedding, unordered_map<int, vector<double>> embed2, unsigned int step) {
  cout << "...EMBED " << step <<endl;
  vector<double> res(dimension);
  TIntFltVH new_embed;
  PWNet graph = PWNet::New();
  TVVec <TInt, int64> walksvv;

  for(int i=0; i < g.nb_nodes; i++){
    int src = i;
    if(!graph->IsNode(src))
        graph->AddNode(src);

    unsigned long long nb_conn_node;
    if(i > 0)
      nb_conn_node = g.out_degrees[i] - g.out_degrees[i-1];
    else
      nb_conn_node = g.out_degrees[i];

    for(unsigned long long j = g.out_degrees[i] - nb_conn_node; j < g.out_degrees[i]; j++){
      int dst = g.unlinks[j];
      int src_node = cnodes[src];
      int dst_node = cnodes[dst];
  
      bool cond = (find(pnodes.begin(), pnodes.end(), src_node) == pnodes.end()) && (find(pnodes.begin(), pnodes.end(), dst_node) == pnodes.end());

      if(!cond){
        if(!graph->IsNode(dst))
          graph->AddNode(dst);
        double weight = 1.;

        // cout << src_node << " " << dst_node << endl;
        graph->AddEdge(src, dst, weight);
        graph->AddEdge(dst, src, weight);
      }     
    }
  }

  int walklength = walklen;
  int numwalk = numwalks;

  if(walklen > g.nb_nodes){
    walklength = g.nb_nodes/4 + 1;
    numwalk = 5;
  }

  cout << "... n2v start " << g.nb_nodes << endl;
  node2vec(graph, ParamP, ParamQ, dimension, walklength, numwalk, winsize, iter, verbose, outputwalks, walksvv, new_embed);
  cout << "... n2v end" << endl;

  double al = pow(exponent, (double)step);

  for(int i = new_embed.FFirstKeyId(); new_embed.FNextKeyId(i);){
    int node = cnodes[new_embed.GetKey(i)];

    for(int j=0; j < new_embed[i].Len(); j++){
      vector<int>::iterator it = find(pnodes.begin(), pnodes.end(), node);
      if(it != pnodes.end()){
        if( embed2.find(node) != embed2.end() )
          embedding[node][j] += al * (0.8 * new_embed[i][j] + 0.2 * embed2[node][j]);
        else
          embedding[node][j] += al * new_embed[i][j];
      }
      else{
        if( embed2.find(node) != embed2.end() )
          embedding[node][j] += al * (0.8 * al * new_embed[i][j] + 0.2 * embed2[node][j]);
        else
          embedding[node][j] += al * al * new_embed[i][j];
      }
    }
  }
}



void
h_louvain(Graph g, unsigned int step, vector<int> nodes, vector<vector<double>> & embedding) {
  unsigned short nb_calls = 0;

  Graph g_new = g;
  init_quality(&g_new, nb_calls);
  nb_calls++;

  cout << "quality: " << q->quality() <<  endl;
  
  Louvain c(-1, precision, q);
  
  bool improvement = true;
  long double quality = (c.qual)->quality();
  long double new_qual;
  
  vector<int> partition(g.nb_nodes);

  for (unsigned int i = 0; i < partition.size(); i++)
    partition[i] = i;

  do {
    improvement = c.one_level();
    new_qual = (c.qual)->quality();

    // cout << "quality: " << q->quality() <<  endl;
    // cout << "neigh_last: " << c.neigh_last << endl;

    c.update_partition(partition);

    g_new = c.partition2graph_binary();
    init_quality(&g_new, nb_calls);
    nb_calls++;

    c = Louvain(-1, precision, q);

    quality = new_qual;
  } while(improvement);


  vector<vector<int> > parts;
  for (unsigned int i = 0; i<partition.size(); i++) {
    if (partition[i] >= parts.size())
      parts.resize(partition[i]+1);
    parts[partition[i]].push_back(i);
  }
  
  cout << "Size: " << parts.size() << endl;
  // cout << "PARTS" << endl;
  // for(unsigned int i=0; i<parts.size(); i++){
    // cout << parts[i].size() << " ";
    // for(unsigned int j=0; j<parts[i].size(); j++){
    //   cout << nodes[parts[i][j]] << " ";
    // }
    // cout << endl;
  // } 
  // cout <<  endl;


  // cout << "COMMUNITY" << endl;
  vector<Graph> com_graph(parts.size());
  vector<vector<int>> com_nodes(parts.size());
  vector<unordered_map<int, int>> com_edge(parts.size());

  for(unsigned int i=0; i<parts.size(); i++){
    set<int> node;
    set<int> cnodes;
    unordered_map<int, int> edge;

    for(unsigned int j=0; j<parts[i].size(); j++){
      node.insert(parts[i][j]);
      cnodes.insert(nodes[parts[i][j]]);

      unsigned long long nb_conn_node;
      if(parts[i][j] > 0)
        nb_conn_node = g.degrees[parts[i][j]] - g.degrees[parts[i][j] - 1];
      else
        nb_conn_node = g.degrees[parts[i][j]];

      for(unsigned long long k = g.degrees[parts[i][j]] - nb_conn_node; k < g.degrees[parts[i][j]]; k++){
        if( ((double)rand() / RAND_MAX) < 0.3){
          node.insert(g.links[k]);
          cnodes.insert(nodes[g.links[k]]);

          // if(std::find(parts[i].begin(), parts[i].end(), g.unlinks[k]) == parts[i].end()){
          //   if(edge.find(parts[i][j]) == edge.end())
          //     edge[parts[i][j]] = g.unlinks[k];
          //   else
          //     edge[g.unlinks[k]] = parts[i][j];
          // }
        }
      }
    }
    vector<int> a(cnodes.begin(), cnodes.end());
    com_nodes[i] = a;

    vector<int> v(node.begin(), node.end());
    com_graph[i] = g.subgraph(v);

    // cout << v.size() << endl;

    // for(unsigned int i1=0; i1<v.size(); i1++)
    //   cout << nodes[v[i1]] << " ";
    // cout << endl;

    // cout << "Nb Links: " << com_graph[i].nb_links << endl;



    // for(unsigned int j=0; j<com_graph[i].nb_nodes; j++){
    //   unsigned long long nb_conn_node;
    //   if(j > 0)
    //     nb_conn_node = com_graph[i].out_degrees[j] - com_graph[i].out_degrees[j - 1];
    //   else
    //     nb_conn_node = com_graph[i].out_degrees[j];

    //   for(unsigned long long k = com_graph[i].out_degrees[j] - nb_conn_node; k < com_graph[i].out_degrees[j]; k++){
    //     cout << j << " " << com_graph[i].unlinks[k] << endl;
    //   }
    // }


    for(unsigned int j = 0 ; j < v.size(); j++){
      unsigned long long nb_conn_node;
      if(v[j] > 0)
        nb_conn_node = g.out_degrees[v[j]] - g.out_degrees[v[j] - 1];
      else
        nb_conn_node = g.out_degrees[v[j]];

      for(unsigned long long k = g.out_degrees[v[j]] - nb_conn_node; k < g.out_degrees[v[j]]; k++){
        bool from = (std::find(parts[i].begin(), parts[i].end(), v[j]) != parts[i].end()) && (std::find(parts[i].begin(), parts[i].end(), g.unlinks[k]) == parts[i].end());
        bool to = (std::find(parts[i].begin(), parts[i].end(), v[j]) == parts[i].end()) && (std::find(parts[i].begin(), parts[i].end(), g.unlinks[k]) != parts[i].end());
        if(from || to){
          if(edge.find(v[j]) == edge.end())
            edge[v[j]] = g.unlinks[k];
          else
            edge[g.unlinks[k]] = v[j];
        }
      }
    }
    com_edge[i] = edge;  
  }

  cout << "Hello" << endl;


   // if there are several parts, then call louvain recursively on the induced subgraph of each part
  if (parts.size() > 1) {
    // cout << "size > 1" << endl;
    unordered_map<int, vector<double>> embed2;
    for (unsigned int i = 0; i < parts.size(); i++) {
      vector<int> new_nodes;
      for (unsigned int j = 0; j < parts[i].size(); j++) {
        new_nodes.push_back(nodes[parts[i][j]]);
      }

      // Second component calc.
      for(auto x: com_edge[i]){
        int node1 = nodes[x.first];
        int node2 = nodes[x.second];

        if((embed2.find(node1) == embed2.end()) && (embed2.find(node2) == embed2.end())){
          vector<double> res(dimension);
          for (unsigned int d = 0; d < dimension; d++)
            res[d] = (double)rand() / RAND_MAX;
          embed2[node1] = res;
          embed2[node2] = res;
        }
        else if((embed2.find(node1) != embed2.end()) && (embed2.find(node2) == embed2.end())){
          embed2[node2] = embed2[node1];
        }
        else if((embed2.find(node1) == embed2.end()) && (embed2.find(node2) != embed2.end())){
          embed2[node1] = embed2[node2];
        }
        else{
          vector<double> res(dimension);
          for(unsigned int d = 0; d < dimension; d++)
            res[d] = (embed2[node1][d] + embed2[node2][d]) / 2;
          embed2[node1] = res;
          embed2[node2] = res;
        }
      }

      multiply_and_add(com_graph[i], new_nodes, com_nodes[i], embedding, embed2, step);

      cout << step << endl;
      Graph s = g.subgraph(parts[i]);
      h_louvain(s, step+1, new_nodes, embedding);
      
    }
  }
  // else generate random vectors for each node in the part
  // vector are generated by multiplying 
  else {
    unordered_map<int, vector<double>> embed2;
    vector<int> new_nodes;
    for (unsigned int i = 0; i < parts[0].size(); i++) {
      new_nodes.push_back(nodes[parts[0][i]]);
      }

      multiply_and_add(com_graph[0], new_nodes, com_nodes[0], embedding, embed2, step);
      
  }


  // // if there are several parts, then call louvain recursively on the induced subgraph of each part
  // if (parts.size() > 1) {
  //   depth++;
  //   for (unsigned int i = 0; i < parts.size(); i++) {
  //     vector<int> new_nodes(parts[i].size());
  //     for (unsigned int j = 0; j < parts[i].size(); j++) {
  //       new_nodes[j] = nodes[parts[i][j]];
  //     }

  //     vector<double> embedding_next = multiply_and_add(embedding);

  //     Graph s = g.subgraph(parts[i]);
  //     h_louvain(s, step+1, new_nodes, embedding_next);
  //   }
  //   depth--;
  // }
  // // else generate random vectors for each node in the part
  // // vector are generated by multiplying 
  // else {
  //   for (unsigned int i = 0; i < parts[0].size(); i++) {
  //     vector<double> embedding_next = multiply_and_add(embedding);
      
  //     cout << nodes[parts[0][i]] ;
  //     for (unsigned int j = 0; j < embedding_next.size(); j++) {
  //       cout << " " << embedding_next[j] / pow(exponent, depth) ;
  //     }
  //     cout << endl;

  //   }
  // }
}


int
main(int argc, char **argv) {
  srand(time(NULL)+getpid());
  
  parse_args(argc, argv);
  
  time_t time_begin, time_end;
  time(&time_begin);
  display_time("Begin");

  Graph g(filename, filename_w, type);

  vector<int> nodes(g.nb_nodes);
  vector<vector<double> > embedding(g.nb_nodes, vector<double> (dimension, 0.0L));  

  for (int i = 0; i < g.nb_nodes; i++) {
    nodes[i] = i;
  }

  h_louvain(g, 0, nodes, embedding);

  string out_file = filename;
  int pos = out_file.find(".");
  out_file = out_file.substr(0, pos) + "-emb_rand.txt";

  ofstream f(out_file);
  if(f.is_open()){
    f << g.nb_nodes << " " << dimension << endl;
    for(unsigned int i = 0; i < embedding.size(); i++){
      f << i;
      for(unsigned int j = 0; j < embedding[i].size(); j++)
        f << " " << embedding[i][j];
      f << endl;
    }
    f.close();
  }
  else{
    cout << "Write not possible" << endl;
  }

  
  time(&time_end);
  display_time("End");
  cerr << "Total duration: " << (time_end-time_begin) << " sec" << endl;

  /*
  g.display();

  vector<int> sub(5);
  sub[0]=0;
  sub[1]=1;
  sub[2]=2;
  sub[3]=32;
  sub[4]=33;

  Graph s = g.subgraph(sub);
  s.display();
  exit(0);
  
  */
}
