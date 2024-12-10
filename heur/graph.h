#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <set>

using namespace std;

struct EDGE {
    int v1, v2;
};

extern int v_num, e_num;
extern int** neighbor;
extern int* degree;
extern double Density;
extern int Max_degree;

// Function declarations
int readGraph(char* file_name);
void freeGraph();
set<int> vndQuasiClique(int maxTime, double threshold);  // Now takes threshold as a parameter

bool isQuasiClique(const set<int>& subset, double threshold);
double calculateDensity(const set<int>& subset);

#endif // GRAPH_H
