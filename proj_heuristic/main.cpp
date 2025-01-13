#include "graph.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <DIMACS_file> <maxTimeInSeconds> <densityThreshold>" << endl;
        return EXIT_FAILURE;
    }

    int maxTime = atoi(argv[2]);
    double threshold = atof(argv[3]);

    if (readGraph(argv[1]) != 0) {
        cerr << "Error reading the graph!" << endl;
        return EXIT_FAILURE;
    }
    int seed = atoi(argv[4]);
    srand(seed);
  
    
    set<int> quasiClique = vndQuasiClique(maxTime, threshold);

   

    freeGraph();
    return EXIT_SUCCESS;
}
