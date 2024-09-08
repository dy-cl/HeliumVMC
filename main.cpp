// main.cpp
#include <iostream>
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

// Returns a list of all electron positions
std::list<std::list<std::pair<std::vector<double>, std::vector<double>>>> generatePositions(int a)
{
    std::list<std::pair<std::vector<double>, std::vector<double>>> currentWalkerPos; // Initialse walker starting positions (list of pairs of vectors)
    std::list<std::list<std::pair<std::vector<double>, std::vector<double>>>> visitedPositions; // Initialise storage of positions (list of lists of pairs of vectors)
    int nWalkers = 100; // Number of walkers
    int nSteps = 10000; // Number of proposed moves

    // Tracking for acceptance rate
    int totalProposals = 0;
    int acceptedMoves = 0;

    for (int i = 0; i < nWalkers; i++)
    {   
        std::vector<double> vector1 = generateRandomVector();
        std::vector<double> vector2 = generateRandomVector();
        currentWalkerPos.push_back(std::make_pair(vector1, vector2)); // Generate initial walker positions

    }

    visitedPositions.push_back(currentWalkerPos); // Add initial posistions to visitedPositions

    for (int i = 0; i < nSteps; i++) // Number of steps
    {
        for (auto& pair : currentWalkerPos) // For each pair of walkers
        {
            totalProposals++;

            // Propose a move for both walkers
            std::vector<double> old1 = pair.first;
            std::vector<double> old2 = pair.second;
            std::vector<double> proposed1 = moveWalker(pair.first);
            std::vector<double> proposed2 = moveWalker(pair.second);

            // Calculate the ratio p
            double p = std::pow((trial(proposed1, proposed2, a)/trial(old1, old2, a)), 2);

            // Acceptance criteria
            if (p >= 1) // If p is greater than 1 accept
            {
                pair.first = proposed1;
                pair.second = proposed2;
                acceptedMoves++;
            }
            if (p < 1){
                double x = generateRandomNum(0,1);
                if (p > x){ // If p is less than 1, accept with probability p
                    pair.first = proposed1;
                    pair.second = proposed2;
                    acceptedMoves++;
                }
            }
            // Else nothing is changed
        }
        visitedPositions.push_back(currentWalkerPos); // Add current positions to visitedPositions
    }

    std::cout << "Acceptance Rate: " << (double)acceptedMoves / totalProposals << std::endl;   

    // Discard the first 10% of visited positions (walkers are thermalising here)
    for (int i = 0; i < nSteps*0.1; i++)
    {
        visitedPositions.pop_front();
    }
    return visitedPositions;
}

int main()
{   
    // vectorUtilitiesTest(); // Test vector utilities
    int a = 0.05; // Example variational parameter, will vary later

    for (double a = 0.05; a <= 0.25; a += 0.05) // For each variational parameter
    {
        std::cout << "a: " << a << std::endl;
        std::list<std::list<std::pair<std::vector<double>, std::vector<double>>>> visitedPositions = generatePositions(a); // Generate positions
        double totalLocalEnergy = 0.0;
        int totalPairs = 0;

        for (const auto& step : visitedPositions) // For each specific step in visitedPositions 
        {
            for (const auto& pair : step) // For each pair of walkers at that step
            {
                double localEnergyValue = localEnergy(pair.first, pair.second, a);
                totalLocalEnergy += localEnergyValue;
                ++totalPairs;
            }
        }

        double localenergy = totalLocalEnergy / totalPairs;
        std::cout << "Energy for a: " << a << ": " << localenergy << std::endl;
        std::cout << std::endl;
    }

    return 0;
}


