// vectorUtilities.cpp
#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <list> 

// Provides a random float between upper and lower bounds
double generateRandomNum(double lower, double upper)
{

    std::random_device rd; 
    std::default_random_engine generator(rd()); 
    std::uniform_real_distribution<double> distribution (lower, upper); 

    return distribution(generator);
}

// Populates a vector with random floats
std::vector<double> generateRandomVector()
{
    double upper = 5;
    double lower = -5;
    std::vector<double> vector(3);

    for (double value : vector) 
    {
        value = generateRandomNum(lower, upper); // Add randomly generated float to element i of vector in bounds [upper, lower]
    }
    return vector;
}

// Prints a vector
void printVector(const std::vector<double> & vector)
{
    for (const double value : vector) 
    {
        std::cout << value << " ";
    }
}

// Finds unit vector rhat along r
std::vector<double> unitVector(const std::vector<double>& vector)
{
    double vectorSum = 0; 
    for (const double value : vector) 
    {
        vectorSum += value*value; 
    }

    double magnitude = std::sqrt(vectorSum); // Find magnitude 
    std::vector<double> unitVector(vector.size()); // Initialize vector with the same size as the input vector
    
    for (int i = 0; i < vector.size(); i++) 
    {
        unitVector[i] = vector[i] / magnitude; 
    }
    return unitVector;
}

// Subtracts two vectors
std::vector<double> subtractVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) 
{
    std::vector<double> newVector(vector1.size());
    for (int i = 0; i < vector1.size(); ++i) {
        newVector[i] = vector1[i] - vector2[i];
    }
    return newVector;
}

// Finds Euclidean distance between two vectors  r12 = |r1−r2|
double vectorDistance(const std::vector<double> & vector1, const std::vector<double> & vector2)
{
    std::vector<double> vector = subtractVectors(vector1, vector2);
    double vectorSum = 0; 

    for (int i = 0; i < vector.size(); i++)
    {
       vectorSum += vector.at(i)*vector.at(i); 
    }

    double distance = std::sqrt(vectorSum); 
    return distance;
}

// Calculates dot product of two vectors
double vectorDotProduct(const std::vector<double> & vector1, const std::vector<double> & vector2)
{   
    double vectorSum = 0;
    for (int i = 0; i < vector1.size(); i++)
    {
        vectorSum += vector1.at(i)*vector2.at(i); 
    }
    return vectorSum;
}

// Move a vector within a cube of values -stepSize, stepSize
std::vector<double> moveWalker(std::vector<double> oldPosition) {
    std::vector<double> newPosition = oldPosition;
    double stepSize = 0.6;
    
    for (size_t i = 0; i < oldPosition.size(); ++i) {
        newPosition[i] += (generateRandomNum(-stepSize, stepSize)); 
    }
    
    return newPosition;
}

// Displays result of above functions for checking
void vectorUtilitiesTest()
{
    std::cout << "TESTING VECTOR UTILS" << std::endl;

    std::vector<double> vector1 = generateRandomVector();
    std::vector<double> vector2 = generateRandomVector();

    std::cout << "Vector1: ";
    printVector(vector1); // Print vector1
    std::cout << std::endl;

    std::cout << "Vector2: ";
    printVector(vector2); // Print vector2
    std::cout << std::endl;

    std::cout << "Unit-Vector1: ";
    printVector(unitVector(vector1)); // Print unit vector of vector1
    std::cout << std::endl;

    std::cout << "Subtract Vectors: ";
    printVector(subtractVectors(vector1, vector2)); // Subtract vector 1 and 2
    std::cout << std::endl;

    std::cout << "Vector Distance: ";
    std::cout << vectorDistance(vector1, vector2); // Euclidean distance between vector 1 and 2
    std::cout << std::endl;

    std::cout << "Vector Dot Product: ";
    std::cout << vectorDotProduct(vector1, vector2); // Dot product between vector 1 and 2
    std::cout << std::endl;
    
}