#ifndef INCLUDED_JACTEST
#define INCLUDED_JACTEST
#include <omp.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "dual.hpp"

/*
This function is to test that the jacobian for a scalar function is computed
correctly and can handle problems of large size

I devised a simple test to accomplish this.
We will use a function: f = x0^2 + x1^2 + ... + xn^2
Each element of the jacobian/gradient for this function will simply be: 2*xi

We will define our function in a for loop over all the independent variables
Each iteration we will add to it the square of the next independent variables
->Due to the AD we can simply call the dual component of the function to get the
  derivative
->Looping and printing out these duals will give us our jacobian.
->We can simply change this from a void to array of doubles and return an array
  containing all the dual components to have a function to return jacobians.
*/
void runJacTest1(const int& N) {
  // allocate our points
  double* points = new double[N];
  // For this test we just want an easy way to check our jacobian is correct
  // We'll accomplish this by knowing the answer for each element should be
  // 2*(i+1) Note: i+1 to avoid awkwardness of 0.
  for (int i = 0; i < N; i++)
    points[i] = i + 1;

// Open an omp parallel region, I want the dual pointer to be private
#pragma omp parallel
  {
    // Push our points into the variable vector
    adoopp::Dual* vars = new adoopp::Dual[N];
// Everything needed is declared in scope, except what we want shared
// anyways
// -> A simple, guided omp loop will suffice
#pragma omp for schedule(guided)
    // Loop over the independent variables, setting their index (wrt)
    for (int wrt = 0; wrt < N; wrt++) {
      // Our function of interest
      adoopp::Dual func;
      // For the size of our problem, add our variables to an array so we can
      // access them easily
      for (int i = N - 1; i > 0; i--)
        // If our variable is the one to be differentiated, set its dual to 0
        // (dx/dx=1)
        if (i == wrt)
          vars[i] = adoopp::Dual(points[i], 1);
        // Otherwise, we treat it as a real number (dc/dx=0)
        else
          vars[i] = adoopp::Dual(points[i], 0);
      double check_deriv{0};
      // The function of interest is just f = \Sigma_{i} x_{i}^2
      for (int i = 0; i < N; i++) {
        func += (vars[i] * vars[i]);
        /*check_deriv*/
        if (i == wrt)
          check_deriv = 2 * vars[i].real();
      }
      if (func.dual() != check_deriv) {
        printf("Jacobian Test 1 Error: derivative wrong at: %d \n", wrt);
        /*check_deriv*/
        // Can't break out of OMP parallel region; will require user attention
        // if parallel; the print statement will have to suffice
#ifndef _OPENMP
        break;
#endif
      }
      if (N <= 10 || (wrt + 1) % 1000 == 0) {
// I have not figured out how to get cout to play nice with omp, so printf it is
#ifdef _OPENMP
        printf("Jac(%d): %0.0f\t by: %d\n", wrt + 1, func.dual(),
               omp_get_thread_num());
#else
        std::cout << "Jac(" << wrt + 1 << "):" << func.dual() << "\n";
#endif
      }
    }
    // Delete our variable array
    delete[] vars;
  }
  delete[] points;
  printf("Jacobian Test 1 Completed with size: %d\n", N);
}

/*
This test is very similar to the last one, but now with a vector func instead

We will use the same idea as JacTest1, except split it up into different fi's
This time we slightly alter the func to give a new jacobian with unique elements

We will assign fj = xj*xi^2. This gives us a nice and easy result we can check:
jac(i,j) = 2*xi*xj + xi*kronecker_delta_ij
*/
void runJacTest2(int N, int M) {
  // allocate our points

  double* points = new double[N];
  // For this test we just want an easy way to check our jacobian is correct
  // We'll accomplish this by knowing the answer for each element should be
  // 2*(i+1) Note: i+1 to avoid awkwardness of 0.
  for (int i = 0; i < N; i++)
    points[i] = i + 1;

#pragma omp parallel
#pragma omp for schedule(guided)
  for (int j = 0; j < M; j++)
  // Open an omp parallel region, I want the dual pointer to be private
  {
    // Push our points into the variable vector
    adoopp::Dual* vars = new adoopp::Dual[N];
    // Everything needed is declared in scope, except what we want shared
    // anyways
    // -> A simple, guided omp loop will suffice
    // Loop over the independent variables, setting their index (wrt)
    for (int wrt = 0; wrt < N; wrt++) {
      // Our function of interest
      adoopp::Dual func;
      // For the size of our problem, add our variables to an array so we can
      // access them easily
      for (int i = 0; i < N; i++)
        // If our variable is the one to be differentiated, set its dual to 0
        // (dx/dx=1)
        if (i == wrt)
          vars[i] = adoopp::Dual(points[i], 1);
        // Otherwise, we treat it as a real number (dc/dx=0)
        else
          vars[i] = adoopp::Dual(points[i], 0);

      // The function of interest is just fij = x_{i}^2*x_j
      // for (int i = 0; i < N; i++) {
      func = (vars[wrt] * vars[wrt] * vars[j]);
      // }
      // Check the dual against the actual derivative, if not equal there's a
      // problem)
      /*deriv_check2*/
      // If (dual != deriv); i!=j
      if (func.dual() != 2 * vars[wrt].real() * vars[j].real())
        // If (dual != deriv); i==j
        if (func.dual() != 2 * vars[wrt].real() * vars[j].real() +
                               vars[wrt].real() * vars[wrt].real()) {
          printf("Jacobian Test 2 Error: derivative wrong at: %d %d\n", j, wrt);
          /*deriv_check2*/
#ifndef _OPENMP
          break;
#endif
        }
      if (N <= 10 || (((j + 1) % 100 == 0) && ((wrt + 1) % 100 == 0))) {
        // I have not figured out how to get cout to play nice with omp, so
        // printf it is
#ifdef _OPENMP
        printf("Jac(%d, %d): %0.0f\t\n", j + 1, wrt + 1, func.dual());
#else
        std::cout << "Jac(" << j + 1 << "," << wrt + 1 << "):" << func.dual()
                  << "\n";
#endif
      }
    }  // Delete our variable array
    delete[] vars;
  }
  delete[] points;
  printf("Jacobian Test 2 Completed with size: %d %d \n", N, M);
}

// If we want a square matrix, just assign both entries the same
// This is likely what I will actually use in the project
void runJacTest2(int N) {
  runJacTest2(N, N);
}

#endif