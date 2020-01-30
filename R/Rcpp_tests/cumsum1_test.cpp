#include <Rcpp.h>
using namespace Rcpp;

// Example of cumsum() function from: 
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
NumericVector cumsum1(NumericVector x){
  // initialize an accumulator variable
  double acc = 0;
  // initialize the result vector
  NumericVector res(x.size());
  for(int i = 0; i < x.size(); i++){
    acc += x[i];
    res[i] = acc;
  }
  return res;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
print('Testing with c(1, 2, 3, 4, 5):')
cumsum1(c(1.0, 2.0, 3.0, 4.0, 5.0))
  */



// Test it with the following R code:
// 
// library(Rcpp)
// sourceCpp('Rcpp_tests/cumsum1_test.cpp')
// cumsum1(c(1.0, 2.0, 3.0, 4.0, 5.0))


