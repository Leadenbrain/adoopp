#ifndef INCLUDED_PLOTTEST
#define INCLUDED_PLOTTEST

#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "dual.hpp"

/* I really like visuals, I find they are a good way to communicate a lot of
 * data very quickly in a way that is easy to parse for humans
 *
 * I'm going to make a couple of 2D surfaces and demonstrate that they produce a
 * proper gradient field by plotting the value of the divergence over the
 * surface. This will give us an equation that we can verify by calculating it
 * by hand and plotting it in python.
 *
 * This is just a quick way to check some of the more complicated AD functions.
 * */

void printPointToFile(const double& x1, const double& x2, std::string output) {
  std::ofstream outFile;
  outFile.open(output, std::ios_base::app);
  adoopp::Dual func;
  adoopp::Dual x(x1, 1), y(x2, 1);
  /*plottest*/
  func = x * sin(y);
  outFile.precision(20);
  outFile << x.real() << "\t" << y.real() << "\t" << func.real() << "\t"
          << func.dual() << "\n";
  /*plottest*/
  outFile.close();
  // For the size of our problem, add our variables to an array so we can
  // access them easily
}

void plotFunc(const int& N, std::string output) {
  double xlim = 2;
  double ylim = 2;
  double h = 2 * xlim / N;
  for (double i = -xlim; i < xlim; i += h)
    for (double j = -ylim; j < ylim; j += h)
      printPointToFile(i, j, output);
  printf("Plot Test Completed, please plot with 'python3 plot.py' \n");
}
#endif