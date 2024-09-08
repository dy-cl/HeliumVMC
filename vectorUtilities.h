// vectorUtilities.h
#ifndef VECTORUTILITIES_H
#define VECTORUTILITIES_H

#include <vector>

// Function declarations
double generateRandomNum(double lower, double upper);
std::vector<double> generateRandomVector();
int printVector(const std::vector<double> &vector);
std::vector<double> unitVector(const std::vector<double> &vector);
std::vector<double> subtractVectors(const std::vector<double> & vector1, const std::vector<double> & vector2);
double vectorDistance(const std::vector<double> &vector1, const std::vector<double> &vector2);
double vectorDotProduct(const std::vector<double> & vector1, const std::vector<double> & vector2);
std::vector<double> moveWalker(std::vector<double> oldPosition);
void vectorUtilitiesTest();

#endif