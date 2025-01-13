
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
#include <chrono>
using namespace std;

extern int v_num;      // Number of vertices
extern int *degree;    // Degree of each vertex
extern int **neighbor; // Adjacency list representation of the graph

bool isAdjacent(int u, int v);
set<int> greedyMaximalClique();
set<int> randomNodeClique();
set<int> largestDegreeClique();
set<int> quasiCliqueThread(int maxTime, double threshold, set<int> &initialClique);

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
        // for (int i : subset)
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

set<int> greedyMaximalClique()
{
    set<int> clique;
    bool *included = new bool[v_num]();

    for (int i = 0; i < v_num; ++i)
    {
        if (clique.empty())
        {
            clique.insert(rand() % v_num);
            included[i] = true;
        }
        else
        {
            bool canAdd = true;
            for (int j : clique)
            {
                if (!isAdjacent(i, j))
                {
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

set<int> randomNodeClique()
{
    set<int> clique;
    clique.insert(rand() % v_num);

    return clique;
}

set<int> largestDegreeClique()
{
    set<int> clique;
    int maxDegree = -1;
    int vertexWithMaxDegree = -1;

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
set<int> quasiCliqueThread(int maxTime, double threshold, set<int> &initialClique)
{
    using namespace std::chrono;

    set<int> bestClique = initialClique, currentClique = initialClique;
    auto startTime = steady_clock::now();

   
    int currentEdges = calculateDensity(currentClique); // Initial edge count in the current clique

    while (duration_cast<seconds>(steady_clock::now() - startTime).count() < maxTime)
    {
        bool improved = false;

        for (int i = 0; i < v_num; i++)
        {
            if (currentClique.count(i))
                continue; 

         
            int newEdges = 0;
            for (int neighborIdx = 0; neighborIdx < degree[i]; ++neighborIdx)
            {
                if (currentClique.count(neighbor[i][neighborIdx]))
                {
                    newEdges++; 
                }
            }

           
            int totalEdges = currentEdges + newEdges;
            int totalVertices = currentClique.size() + 1;
            double maxPossibleEdges = totalVertices * (totalVertices - 1) / 2.0;
            double density = (maxPossibleEdges == 0) ? 0.0 : totalEdges / maxPossibleEdges;

            
            if (density >= threshold)
            {
                currentClique.insert(i);
                currentEdges = totalEdges;  
                if (currentClique.size() > bestClique.size())
                {
                    bestClique = currentClique; 
                    improved = true;
                }
            }
        }

       
        if (!improved && !currentClique.empty())
        {
            // Randomly select a node from the current clique to remove
            auto it = currentClique.begin();
            advance(it, rand() % currentClique.size());  // Randomly select a node
            currentClique.erase(it);  // Remove the selected node
        }

        
    }

    return bestClique;
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

bool isQuasiClique(const set<int> &clique, double threshold, int **neighbor, int *degree)
{
    return calculateDensity(clique, neighbor, degree) >= threshold;
}

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

set<int> vndQuasiClique(int maxTime, double threshold)
{
    set<int> bestClique1;
    set<int> bestClique2;
    set<int> bestClique3;
    thread greedyThread, randomThread, largestDegreeThread;

    set<int> greedyClique = greedyMaximalClique();
    set<int> randomClique = randomNodeClique();
    set<int> largestDegreeCliqueSet = largestDegreeClique();

    // Start the threads with their respective initial cliques
    greedyThread = thread([&]()
                          {
        set<int> result = quasiCliqueThread(maxTime, threshold, greedyClique);
        if (result.size() > bestClique1.size()) {
            bestClique1 = result;
        } });

    randomThread = thread([&]()
                          {
        set<int> result = quasiCliqueThread(maxTime, threshold, randomClique);
        if (result.size() > bestClique2.size()) {
            bestClique2 = result;
        } });

    largestDegreeThread = thread([&]()
                                 {
        set<int> result = quasiCliqueThread(maxTime, threshold, largestDegreeCliqueSet);
        if (result.size() > bestClique3.size()) {
            bestClique3 = result;
        } });

    // join all threads
    greedyThread.join();
    randomThread.join();
    largestDegreeThread.join();

    set<int> mergedClique = mergeCliques(bestClique1, bestClique2, bestClique3, threshold, neighbor, degree);

    std::vector<std::string> cliqueNames = {"Greedy", "Random", "Largest Degree", "Merged"};
    std::vector<std::pair<size_t, double>> cliqueData = {
        {bestClique1.size(), calculateDensity(bestClique1)},
        {bestClique2.size(), calculateDensity(bestClique2)},
        {bestClique3.size(), calculateDensity(bestClique3)},
        {mergedClique.size(), calculateDensity(mergedClique)}};

    size_t largestIndex = 0;
    for (size_t i = 1; i < cliqueData.size(); ++i)
    {
        if (cliqueData[i].first > cliqueData[largestIndex].first)
        {
            largestIndex = i;
        }
    }

    // Print the largest clique
    std::cout << cliqueData[largestIndex].first;

    return bestClique1;
}
