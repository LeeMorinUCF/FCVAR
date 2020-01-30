#include <Rcpp.h>
using namespace Rcpp;

// Based on example of cumsum() function from: 
// http://gallery.rcpp.org/articles/vector-cumulative-sum/



// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
NumericVector X_c_gen_cpp_v2(NumericVector dx_t, double c){
  // initialize the increment variables.
  // double dx_t = 0;
  double dx_c_t = 0;
  // initialize the result vector. 
  int T_x = dx_t.size();
  NumericVector x_c(dx_t.size() + 1);
  x_c[0] = 0;
  for(int i = 1; i < dx_t.size() + 1; i++){
    // dx_t = x[i] - x[i-1];
    dx_c_t = x_c[i-1]*c/T_x + dx_t[i-1];
    x_c[i] = x_c[i-1] + dx_c_t;
  }
  return x_c;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
  # Testing with dX = rep(1.0, 4) and c = 0.0:
  X_c_gen_cpp_v2(rep(1.0, 4), 0.0)
  # Result should be:
  # [1] 0 1 2 3 4
  # 
  # Testing with dX = rep(1.0, 4) and c = 1.0:
    X_c_gen_cpp_v2(rep(1.0, 4), 1.0)
  # Result should be:
  # [1] 0.000000 1.000000 2.250000 3.812500 5.765625
  # 
  # Testing with dX = rep(1.0, 4) and c = - 1.0:
  X_c_gen_cpp_v2(rep(1.0, 4), - 1.0)
  # Result should be:
  # [1] 0.000000 1.000000 1.750000 2.312500 2.734375
  # 
  # Testing with dX = c(0, 1, -1, 1) and c = - 1.0:
  X_c_gen_cpp_v2(c(0, 1, -1, 1), - 1.0)
  # Result should be:
  # [1] 0.000000 0.000000 1.000000 -0.250000 0.812500
*/



// Test it with the following R code:
// 
// library(Rcpp)
// sourceCpp('Rcpp_tests/X_c_gen_v2.cpp')
// X_c_gen_cpp_v2(rep(1.0, 4), - 1.0)


