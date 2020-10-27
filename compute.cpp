#include <fstream>
#include <iostream>
#include <vector>
#include "dual.hpp"
#include "fadoopp.hpp"

// typedef void* ADFcn;
// This is our function to simulate a given DE
// void runSimulation(int N, ADFcn AD_FCN) {}
/*
Our main function to run.
We *need* to:
    -Define a DE/DAE
    -Compute DE/DAE
    -Output Result of DE/DAE
We *want* to:
    -Time the computation?
    -Loop many times to check memory leaks
    -Call the computation for numerous D(A)Es
*/

int main(void) {
  // Set our initial values for our variables
  std::vector<double> points;
  points.push_back(0.78);
  points.push_back(4);
  points.push_back(-2.3);
  std::vector<adoopp::Dual> vars;
  // Push our points into the variable vector
  for (int i = 0; i < points.size(); i++)
    vars.push_back(adoopp::Dual(points[i], 1));

  // Create our forward AD object, print out grad, print out jacobian
  adoopp::Fadoopp ff1(points, vars);
  printf("FADOOPP GRAD: %f \n", ff1.grad());
  for (int i = 0; i < ff1.size(); i++)
    printf("Using FADOOPP: %f \n", ff1.jacElem(i));
}