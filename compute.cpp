#include "jactest.hpp"
#include "plottest.hpp"
// I like to keep my main relatively empty so I can see what is going on
// The purpose of this main function is just to run through the performative
// tests needed to demonstrate the AD capacity.

adoopp::Dual func(adoopp::Dual a) {
  return pow(a, 3);
}
/*mainf*/
int main(void) {
  const int sizeTest1 = 100000;  // This problem is O(n*cost(f)), can be higher
  const int sizeTest2 = 5000;  // This problem is O(n$^2$*cost(f)); takes longer
  const int sizeTest3 = 30;  // Too many points will make our plot crowded (also
                             // python is very slow)
  // Expected behavior: "Jacobian Test 1 Completed: N" printed to terminal
  // runJacTest1(sizeTest1);
#define OMP_NUM_THREAS = 4;  // This test is poorly optimized, 4 is maximal gain
  // Expected behavior: "Jacobian Test 2 Completed: N N" printed to terminal
  runJacTest2(sizeTest2);
  std::string outputFile1 = "file1.out";
  // Expected bavior: tsv of x*sin(y) and surface tangent to (sin(y) +x(cos*y))
  // Plot with python3 plot.py (need numpy and matplotlib)
  // plotFunc(sizeTest3, outputFile1);
  /*mainf*/
}
