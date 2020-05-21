#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector likelihood_C(NumericVector media, NumericVector correctedFixTime, NumericVector tim, int sum_correctedFixTime,
                           double stateStep, NumericMatrix changeMatrix, NumericMatrix prStates, double sigma, 
                           NumericMatrix changeUp, NumericMatrix changeDown) {
  
  int i = 0;
  double mean=0;
  NumericVector probUpCrossing(sum_correctedFixTime);
  NumericVector probDownCrossing(sum_correctedFixTime);
  NumericMatrix prStatesNew(21, sum_correctedFixTime);
  NumericVector tempUpCross(sum_correctedFixTime);
  NumericVector tempDownCross(sum_correctedFixTime);
  NumericVector sumIn(sum_correctedFixTime);
  NumericVector sumCurrent(sum_correctedFixTime);
  NumericMatrix prStatesNew_f(21, sum_correctedFixTime);
  NumericVector tempUpCross_f(sum_correctedFixTime);
  NumericVector tempDownCross_f(sum_correctedFixTime);
  NumericMatrix pdf(21, sum_correctedFixTime);
  NumericMatrix pdf_up(21,sum_correctedFixTime);
  NumericMatrix pdf_down(21,sum_correctedFixTime);
  
  for(int time = 1; time < sum_correctedFixTime; ++time) {
    if ( time == (tim[i]) ){
      ++i;
      mean = media[i];
    } else if(time == 1) {
      mean = media[0];
    }
    
    for(int i = 0; i < 21; ++i){
      pdf(_,time)= dnorm(changeMatrix(_,i), mean, sigma);
      prStatesNew(i, time) = sum(pdf(_,time) * prStates(_,time-1)) * stateStep;
    }
    
    pdf_up(_,time)= 1-pnorm(changeUp(time, _ ), mean, sigma);
    tempUpCross[time] = sum(pdf_up(_,time) * prStates( _ , time-1));
    
    pdf_down(_,time)= pnorm(changeDown(time,_), mean, sigma);
    tempDownCross[time] = sum(pdf_down(_,time) * prStates(_, time-1));
    
    sumIn[time] = sum(prStates(_, time-1));
    sumCurrent[time] = sum(prStatesNew(_,time)) + tempUpCross[time] + tempDownCross[time];
    
    prStatesNew_f(_, time) =  (prStatesNew(_,time) * sumIn[time])/ sumCurrent[time];
    tempUpCross_f[time] = (tempUpCross[time] * sumIn[time]) / sumCurrent[time];
    tempDownCross_f[time] = (tempDownCross[time] * sumIn[time]) / sumCurrent[time];
    
    prStates(_, time) = prStatesNew_f(_,time);
    probUpCrossing[time] = tempUpCross_f[time];
    probDownCrossing[time] = tempDownCross_f[time];
  }
  
  NumericVector result = {probUpCrossing[sum_correctedFixTime-1], probDownCrossing[sum_correctedFixTime-1]};
  
  // if (result[1] < 1e-10) result[1] = 1e-10;
  // if (result[2] < 1e-10) result[2] = 1e-10;
  
  return(result);
}  