
#include <set>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <mutex>
#include "graph.h"
#include <ctime>
#include <set>
#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <unordered_set>

using namespace std;

// Global variables for graph representation
extern int v_num;      // Number of vertices
extern int *degree;    // Degree of each vertex
extern int **neighbor; // Adjacency list representation of the graph

// Function declarations
bool isAdjacent(int u, int v);
set<int> greedyMaximalClique();
set<int> randomNodeClique();
set<int> largestDegreeClique();
set<int> vndQuasiCliqueThread(int maxTime, double threshold, set<int> &initialClique);

// Global variables
int v_num = 0, e_num = 0;
int **neighbor = nullptr;
int *degree = nullptr;
double Density = 0.0;
int Max_degree = 0;

int readGraph(char *file_name)
{
    int u, v, count = 0;
    char temp1[10], temp2[10];
    ifstream ffs;

    ffs.open(file_name);
    if (ffs.fail())
    {
        printf("### Error Open, File Name:%s\n", file_name);
        return 1;
    }

    ffs >> temp1 >> temp2 >> v_num >> e_num;
    neighbor = (int **)malloc(v_num * sizeof(int *));
    degree = (int *)malloc(v_num * sizeof(int));
    EDGE *edge = (EDGE *)malloc(e_num * sizeof(EDGE));

    memset(degree, 0, v_num * sizeof(int));

    while (ffs >> temp1 >> u >> v)
    {
        u--;
        v--;
        degree[u]++;
        degree[v]++;
        edge[count].v1 = u;
        edge[count].v2 = v;
        count++;
    }

    Density = 2.0 * e_num / (v_num - 1) / v_num;
    Max_degree = 0;

    for (int i = 0; i < v_num; ++i)
    {
        neighbor[i] = (int *)malloc(degree[i] * sizeof(int));
        if (degree[i] > Max_degree)
            Max_degree = degree[i];
        degree[i] = 0;
    }

    for (int i = 0; i < e_num; ++i)
    {
        u = edge[i].v1;
        v = edge[i].v2;
        neighbor[u][degree[u]++] = v;
        neighbor[v][degree[v]++] = u;
    }

    free(edge);
    ffs.close();
    return 0;
}

void freeGraph()
{
    for (int i = 0; i < v_num; ++i)
    {
        free(neighbor[i]);
    }
    free(neighbor);
    free(degree);
}

double calculateDensity(const set<int> &subset)
{
    int edges = 0;
    for (int u : subset)
    {
        for (int i = 0; i < degree[u]; ++i)
        {
            int v = neighbor[u][i];
            if (subset.count(v))
                edges++;
        }
    }
    int n = subset.size();
    return (n <= 1) ? 0.0 : (double)edges / (n * (n - 1));
}

bool isQuasiClique(const set<int> &subset, double threshold)
{
    return calculateDensity(subset) >= threshold;
}
bool isAlmostQuasiClique(const set<int> &subset, double threshold)
{
    return calculateDensity(subset) >= 0.9 * threshold;
}

// Check if two vertices u and v are adjacent
bool isAdjacent(int u, int v)
{
    for (int i = 0; i < degree[u]; ++i)
    {
        if (neighbor[u][i] == v)
        {
            return true;
        }
    }
    return false;
}

// Generate a greedy maximal clique
set<int> greedyMaximalClique()
{
    set<int> clique;
    bool *included = new bool[v_num]();

    for (int i = 0; i < v_num; ++i)
    {
        if (clique.empty())
        {
            clique.insert(rand() % v_num); // Start with any vertex
            included[i] = true;
        }
        else
        {
            bool canAdd = true;
            for (int j : clique)
            {
                if (!isAdjacent(i, j))
                { // Not adjacent to all clique members
                    canAdd = false;
                    break;
                }
            }
            if (canAdd)
            {
                clique.insert(i);
                included[i] = true;
            }
        }
    }
    delete[] included;
    return clique;
}

// Generate a random initial clique with a single random vertex
set<int> randomNodeClique()
{
    set<int> clique;
    clique.insert(rand() % v_num);
    // clique.insert(rand() % v_num);  // Pick a random vertex to start
    return clique;
}

// Generate an initial clique starting with the largest degree vertex
set<int> largestDegreeClique()
{
    set<int> clique;
    int maxDegree = -1;
    int vertexWithMaxDegree = -1;

    // Find the vertex with the largest degree
    for (int i = 0; i < v_num; ++i)
    {
        if (degree[i] > maxDegree)
        {
            maxDegree = degree[i];
            vertexWithMaxDegree = i;
        }
    }

    clique.insert(vertexWithMaxDegree);
    return clique;
}

// set<int> vndQuasiCliqueThread(int maxTime, double threshold, set<int>& initialClique)
// {
//     set<int> bestClique = initialClique, currentClique = initialClique;
//     clock_t startTime = clock();

//     int iterations = 0;
//     int tabu = 0;
//     int noImprovementIterations = 0; // Track iterations without improvement
   

//     // Main VND loop for the given thread
//     while ((clock() - startTime) / CLOCKS_PER_SEC < maxTime)
//     {
//         bool improved = false;

//         // Try to improve the current clique by adding vertices
//         for (int i = 0; i < v_num; i++)
//         {
//             if (currentClique.count(i))  // Skip if already in the clique or the node has been removed
//                 continue;

//             // Calculate the number of neighbors the node has inside the current clique
//             int neighborCount = 0;
//             for (int neighborNode : currentClique)
//             {
//                 if (neighbor[i][neighborNode])  // Check if there's an edge between i and neighborNode
//                 {
//                     neighborCount++;
//                 }
//             }

//             // If there are no neighbors inside the clique, skip
//             if (neighborCount == 0)
//                 continue;

//             // Calculate the probability to add this node
//             double probability = static_cast<double>(currentClique.size()) / neighborCount;

//             // Generate a random number to decide whether to add the node
//             if (rand() / static_cast<double>(RAND_MAX) < probability)
//             {
//                 currentClique.insert(i); // Add the node with the calculated probability

//                 if (isQuasiClique(currentClique, threshold))
//                 {
//                     if (currentClique.size() > bestClique.size())
//                     {
//                         bestClique = currentClique;
//                         improved = true;
//                         noImprovementIterations = 0;  // Reset no improvement counter
//                     }
//                 }
//                 else
//                 {
//                     currentClique.erase(i);  // Remove the node if it doesn't fit the clique
//                 }
//             }
//         }

//         // If no improvement in the last 10 iterations, remove the node with the least connections
//         if (!improved)
//         {
//             noImprovementIterations++;

//             if (noImprovementIterations >= 10)
//             {
//                 // Find the node with the least number of connections inside the clique
//                 int nodeToRemove = -1;
//                 int minConnections = INT_MAX;

//                 // Check each node in the current clique
//                 for (int node : currentClique)
//                 {
//                     int connectionCount = 0;
//                     for (int neighborNode : currentClique)
//                     {
//                         if (neighbor[node][neighborNode])  // Check if there is an edge between `node` and `neighborNode`
//                         {
//                             connectionCount++;
//                         }
//                     }

//                     // If this node has fewer connections, mark it for removal
//                     // If connection counts are equal, choose randomly between the nodes
//                     if (connectionCount < minConnections)
//                     {
//                         minConnections = connectionCount;
//                         nodeToRemove = node;
//                     }
//                     else if (connectionCount == minConnections)
//                     {
//                         // Randomly choose between the nodes with the same connection count
//                         if (rand() % 2 == 0)
//                         {
//                             nodeToRemove = node;
//                         }
//                     }
//                 }

//                 // Remove the node with the least connections (or randomly chosen one in case of tie)
//                 if (nodeToRemove != -1)
//                 {
//                     currentClique.erase(nodeToRemove);
                    
//                 }

//                 // Reset the counter after a removal
//                 noImprovementIterations = 0;
//             }
//         }

//         iterations++;
//     }

//     cout << "Iterations: " << iterations << endl;
//     return bestClique; // Return the best clique found by this thread
// }

// VND procedure executed by each thread (this will run independently)
set<int> vndQuasiCliqueThread(int maxTime, double threshold, set<int> &initialClique)
{
    set<int> bestClique = initialClique, currentClique = initialClique;
    clock_t startTime = clock();

    int iterations = 0;
    int tabu = 0;
    // Main VND loop for the given thread
    while ((clock() - startTime) / CLOCKS_PER_SEC < maxTime)
    {
        bool improved = false;

        // Try to improve the current clique by adding vertices

        for (int i = 0; i < v_num; i++)
        {
            if (currentClique.count(i))
                continue; // Skip if already in the clique
            currentClique.insert(i);

            if (isQuasiClique(currentClique, threshold))
            {
                if (currentClique.size() > bestClique.size())
                {
                    bestClique = currentClique;
                    improved = true;
                }
            }
            else
            {
                currentClique.erase(i); // Remove the vertex if it doesn't fit the clique
                if (iterations % 10 == 0)
                {

                    currentClique.erase(rand() % v_num);
                }
                 iterations+=1;
            }
        }
    }
    cout << iterations << endl;
    return bestClique; // Return the best clique found by this thread
}
double calculateDensity(const set<int> &clique, int **neighbor, int *degree)
{
    int edges = 0;
    for (int u : clique)
    {
        for (int i = 0; i < degree[u]; ++i)
        {
            int v = neighbor[u][i];
            if (clique.count(v))
                edges++;
        }
    }
    int n = clique.size();
    return (n <= 1) ? 0.0 : (double)edges / (n * (n - 1));
}

// Check if the set is a quasi-clique (density above threshold)
bool isQuasiClique(const set<int> &clique, double threshold, int **neighbor, int *degree)
{
    return calculateDensity(clique, neighbor, degree) >= threshold;
}
// Function to merge three cliques and remove nodes with the fewest neighbors until density is above threshold
set<int> mergeCliques(set<int> &clique1, set<int> &clique2, set<int> &clique3, double threshold, int **neighbor, int *degree)
{
    // Step 1: Union all three cliques
    set<int> mergedClique = clique1;
    mergedClique.insert(clique2.begin(), clique2.end());
    mergedClique.insert(clique3.begin(), clique3.end());

    set<int> currentClique = mergedClique;

    while (!isQuasiClique(currentClique, threshold, neighbor, degree))
    {

        vector<pair<int, int>> nodeNeighbors;

        for (int node : currentClique)
        {
            int neighborCount = 0;
            for (int i = 0; i < degree[node]; ++i)
            {
                int neighborNode = neighbor[node][i];
                if (currentClique.count(neighborNode))
                {
                    neighborCount++;
                }
            }
            nodeNeighbors.push_back({node, neighborCount});
        }

        // Sort the nodes based on their neighbor count in the current clique (ascending)
        sort(nodeNeighbors.begin(), nodeNeighbors.end(), [](const pair<int, int> &a, const pair<int, int> &b)
             { return a.second < b.second; });

        // Remove the node with the fewest neighbors
        currentClique.erase(nodeNeighbors.front().first);
    }

    return currentClique;
}
// Modify the multithreaded VND Quasi-Clique Search to include merging step
set<int> vndQuasiClique(int maxTime, double threshold)
{
    set<int> bestClique1;
    set<int> bestClique2;
    set<int> bestClique3;
    thread greedyThread, randomThread, largestDegreeThread;

    // Create the initial cliques for each thread
    set<int> greedyClique = greedyMaximalClique();
    set<int> randomClique = randomNodeClique();
    set<int> largestDegreeCliqueSet = largestDegreeClique();

    // Start the threads with their respective initial cliques
    greedyThread = thread([&]()
                          {
        set<int> result = vndQuasiCliqueThread(maxTime, threshold, greedyClique);
        if (result.size() > bestClique1.size()) {
            bestClique1 = result;
        } });

    randomThread = thread([&]()
                          {
        set<int> result = vndQuasiCliqueThread(maxTime, threshold, randomClique);
        if (result.size() > bestClique2.size()) {
            bestClique2 = result;
        } });

    largestDegreeThread = thread([&]()
                                 {
        set<int> result = vndQuasiCliqueThread(maxTime, threshold, largestDegreeCliqueSet);
        if (result.size() > bestClique3.size()) {
            bestClique3 = result;
        } });

    // Join all threads
    greedyThread.join();
    randomThread.join();
    largestDegreeThread.join();

    // Merge the best cliques found by the threads and see if there's a better one
    // set<int> mergedClique = mergeCliques(bestClique1, bestClique2, bestClique3, threshold);
    set<int> mergedClique = mergeCliques(bestClique1, bestClique2, bestClique3, threshold, neighbor, degree);

    // Optionally, print the sizes of the cliques to compare
    cout << "Best Clique from Greedy: " << bestClique1.size() << "d =" << calculateDensity(bestClique1) << endl;
    cout << "Best Clique from Random: " << bestClique2.size() << "d =" << calculateDensity(bestClique2) << endl;
    cout << "Best Clique from Largest Degree: " << bestClique3.size() << "d =" << calculateDensity(bestClique3) << endl;
    cout << "Merged Clique: " << mergedClique.size() << "d =" << calculateDensity(mergedClique) << endl;

    // Return the merged clique as the final result
    return bestClique1;
}
