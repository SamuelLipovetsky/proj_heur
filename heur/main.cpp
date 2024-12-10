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

    cout << "Running VND with threshold: " << threshold << endl;
    set<int> quasiClique = vndQuasiClique(maxTime, threshold);

    // cout << "Found Quasi-Clique of size: " << quasiClique.size() << "\nNodes: ";
    // for (int node : quasiClique) {
    //     cout << (node + 1) << " ";  // Convert back to 1-based indexing
    // }
    // cout << endl;

    freeGraph();
    return EXIT_SUCCESS;
}
