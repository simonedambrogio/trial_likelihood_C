Rcpp::Rcpp.package.skeleton("likelihood",cpp_files = "get_trial_likelihood_C.cpp")
Rcpp::compileAttributes("likelihood")
install.packages("likelihood", repos=NULL, type="source")
