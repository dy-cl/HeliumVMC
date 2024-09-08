// main.cpp
#include <iostream>
#include <fstream>
#include <math.h> 
#include <list>
#include <utility>
#include "vectorUtilities.h"

// Returns value of trial wavefunction for two given electron positions
double trial(std::vector<double> pos1, std::vector<double> pos2, double a)
{
    std::vector<double> zeroVector = {0, 0, 0}; // Initialise zero vector as the origin
    double r1 = vectorDistance(pos1, zeroVector); // Distance of electron 1 from origin
    double r2 = vectorDistance(pos2, zeroVector); // Distance of electron 2 from origin
    double r12 = vectorDistance(pos1, pos2); // Distance of electrons from one another

    double trial_psi = exp(-2*r1)*exp(-2*r2)*exp(r12/(2*(1 + a*r12))); // Trial wavefunction
    return trial_psi;
}

// Returns value of local energy for two given electron positions
double localEnergy(std::vector<double> pos1, std::vector<double> pos2, double a)
{
    double r12 = vectorDistance(pos1, pos2); // Distance of electrons from one another
    std::vector<double> unitpos1 = unitVector(pos1); // UnitVector of electron 1
    std::vector<double> unitpos2 = unitVector(pos2); // UnitVector of electron 2

    // Assemble expression for local energy
    double term1 = (vectorDotProduct(subtractVectors(unitpos1, unitpos2), subtractVectors(pos1, pos2)))*(1 / (r12*std::pow(1 + a*r12, 2)));
    double term2 = (1 / (r12*std::pow(1 + a*r12, 3)));
    double term3 = (1 / (4*std::pow(1 + a*r12, 4)));
    double localEnergy = -4 + term1 - term2 - term3 + 1/r12;

    return localEnergy;
}

// Initialize walkers with random starting positions
std::list<std::pair<std::vector<double>, std::vector<double>>> initializeWalkers(int nWalkers)
{
    std::list<std::pair<std::vector<double>, std::vector<double>>> walkerPositions;
    for (int i = 0; i < nWalkers; i++) {
        walkerPositions.push_back({generateRandomVector(), generateRandomVector()});
    }
    return walkerPositions;
}

// Propose a move for both walkers and check acceptance
bool isMoveAccepted(const std::pair<std::vector<double>, std::vector<double>>& oldPositions, 
                    const std::pair<std::vector<double>, std::vector<double>>& proposedPositions, 
                    double alpha)
{
    double p = std::pow(trial(proposedPositions.first, proposedPositions.second, alpha) / 
                        trial(oldPositions.first, oldPositions.second, alpha), 2);
    
    if (p >= 1) return true; // Accept move if p greater than 1 or
    return p > generateRandomNum(0, 1); // Accept with probability p
}

// Perform walker simulation and store visited positions
std::list<std::list<std::pair<std::vector<double>, std::vector<double>>>> simulateWalkers(double variationalParameter)
{   
    int nSteps = 34000;
    int nWalkers = 400;
    auto walkerPositions = initializeWalkers(nWalkers);
    std::list<std::list<std::pair<std::vector<double>, std::vector<double>>>> visitedPositions; //Create list to store all visited positions
    visitedPositions.push_back(walkerPositions); // Add initial positions to visitedPositions

    int acceptedMoves = 0, totalProposals = 0;

    for (int step = 0; step < nSteps; step++) { // nSteps iterations
        for (auto & walker : walkerPositions) { // For each walker
            totalProposals++;
            std::pair<std::vector<double>, std::vector<double>> proposedPositions = { // Make proposal
                moveWalker(walker.first),
                moveWalker(walker.second)
            };

            if (isMoveAccepted(walker, proposedPositions, variationalParameter)) { // Accept or reject proposal
                walker = proposedPositions;
                acceptedMoves++;
            }
        }
        visitedPositions.push_back(walkerPositions); // Record the positions at each step
    }

    std::cout << "Acceptance Rate: " << static_cast<double>(acceptedMoves) / totalProposals << std::endl;

    // Discard the first 10% of visited positions for thermalisation
    for (int i = 0; i < nSteps * 0.1; ++i) {
        visitedPositions.pop_front();
    }

    return visitedPositions;
}

int main()
{   
    std::ofstream outFile("results.txt"); // Open file for writing

    for (double alpha = 0.05; alpha <= 0.25; alpha += 0.05) { // For each value of the variational parameter
        std::cout << "Variational Parameter (alpha): " << alpha << std::endl;

        auto visitedPositions = simulateWalkers(alpha); // Generate the visited positions
        double totalLocalEnergy = 0.0;
        int totalPairs = 0;

        for (const auto & step : visitedPositions) { // For each step
            for (const auto & walker : step) { // For each pair of walkers at that step
                totalLocalEnergy += localEnergy(walker.first, walker.second, alpha);
                ++totalPairs;
            }
        }

        double averageLocalEnergy = totalLocalEnergy / totalPairs; // Calculate final energy averaged over all pairs of walkers and steps
        std::cout << "Average Local Energy: " << averageLocalEnergy << std::endl;
        outFile << "Average Local Energy for alpha: " << alpha << ": " << averageLocalEnergy << std::endl;
        std::cout << std::endl;
        outFile << std::endl;
    }

    outFile.close(); // Close the file
    return 0;
}


