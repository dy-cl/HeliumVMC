// vectorUtilities.h
#ifndef VECTORUTILITIES_H
#define VECTORUTILITIES_H

#include <vector>

// Provides a random float between upper and lower bounds
double generateRandomNum(double lower, double upper);

// Populates a vector with random floats
std::vector<double> generateRandomVector();

// Prints a vector
int printVector(const std::vector<double> &vector);

// Finds unit vector rhat along r
std::vector<double> unitVector(const std::vector<double> &vector);

// Subtracts two vectors
std::vector<double> subtractVectors(const std::vector<double> & vector1, const std::vector<double> & vector2);

// Finds Euclidean distance between two vectors  r12 = |r1−r2|
double vectorDistance(const std::vector<double> &vector1, const std::vector<double> &vector2);

// Calculates dot product of two vectors
double vectorDotProduct(const std::vector<double> & vector1, const std::vector<double> & vector2);

// Move a vector within a cube of values -stepSize, stepSize
std::vector<double> moveWalker(std::vector<double> oldPosition);

// Displays result of above functions for checking
void vectorUtilitiesTest();

#endif