#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double derivs_x_wr_t_Rcpp(const int GrowthCrvChoice, const int function_call,
                          double x, double y, NumericVector params) {

  // Calculate the derivatives of x w.r. t at
  // the current values of x, returning these as dxdt

  // The variable x refers to the size of the animal and t refers to the age of
  // the abalone (relative to some assumed age)

  double InitialLength;
  double dxdt;
  double L50_1;
  double L95_1;
  double L50_2;
  double L95_2;
  double Max_increment;
  double Prop1;
  double Prop2;
  double Gaussian_A;
  double Gaussian_u;
  double Gaussian_sd;
  double vb_Linf;
  double vb_K;
  double Gomp_Linf;
  double Gomp_G;
  dxdt = 0;

  if (function_call == 1) {
    InitialLength = x;
  }
  else {
    InitialLength = y;
  }

  // Calculate the terms that are used to produce the estimate of the derivative
  // of length
  if (GrowthCrvChoice == 1)   { // double logistic

    // L50_1, L95_1, L50_2, L95_2, Max_increment (double logistic model)
    L50_1 = exp(params(0));
    L95_1 = exp(params(1));
    L50_2 = exp(params(2));
    L95_2 = exp(params(3));
    Max_increment = exp(params(4));
    Prop1 = 1 / (1 + exp(-log(19) * (InitialLength - L50_1) / (L95_1 - L50_1)));
    Prop2 = 1 / (1 + exp(log(19) * (InitialLength - L50_2) / (L95_2 - L50_2)));

    // Calculate the derivative
    dxdt = Max_increment * Prop1 * Prop2;
  }

  if (GrowthCrvChoice == 2)   { // Gaussian function

    // Gaussian_A, Gaussian_u, Gaussian_sd (Gaussian function)
    Gaussian_A = exp(params(0));
    Gaussian_u = exp(params(1));
    Gaussian_sd = exp(params(2));

    dxdt = Gaussian_A * exp(-((InitialLength - Gaussian_u) *
      (InitialLength - Gaussian_u)) / (2 * Gaussian_sd * Gaussian_sd));
  }

  if (GrowthCrvChoice == 3)   { // von Bertalanffy growth curve. delta_t is
    // converted t to decimal years (i.e., 1/52)

    // vb_Linf, vb_K (von Bertalanffy)
    vb_Linf = exp(params(0));
    vb_K = exp(params(1));
    dxdt = (vb_Linf - InitialLength) * (1 - exp(-vb_K * (1./365)));
  }

  if (GrowthCrvChoice == 4)   { // Gompertz growth curve. delta_t is
    // converted to decimal years (i.e., 1/52)

    // Gomp_Linf, Gomp_G
    Gomp_Linf = exp(params(0));
    Gomp_G = exp(params(1));
    dxdt = Gomp_Linf * pow((InitialLength / Gomp_Linf),exp(-Gomp_G * (1./365))) - InitialLength;
  }

  return(dxdt);
}

// [[Rcpp::export]]
double RK4SYS_Rcpp(const int nstep, const int GrowthCrvChoice, const double h,
                   double t, double x, NumericVector params)  {

  // Runge-Kutta 4'th order integration of
  //  a system of n ordinary differential equations
  // of x() w.r. t, using a step size of h
  // over nstep steps.
  // number of derivatives - for runga kutta method of integration
  int k;
  int function_call;
  double h2;
  double Start;
  double dxdt;
  double F1;
  double F2;
  double F3;
  double F4;
  double y;
  y = 0;
  h2 = 0.5 * h;
  Start = t;

  for (k=0; k<nstep; k++) {

    function_call = 1;
    dxdt = derivs_x_wr_t_Rcpp(GrowthCrvChoice, function_call, x, y, params);
    F1 = dxdt;
    y = x + h2 * F1;

    function_call = 2;
    dxdt = derivs_x_wr_t_Rcpp(GrowthCrvChoice, function_call, x, y, params);
    F2 = dxdt;
    y = x + h2 * F2;

    function_call = 3;
    dxdt = derivs_x_wr_t_Rcpp(GrowthCrvChoice, function_call, x, y, params);
    F3 = dxdt;
    y = x + h * F3;

    function_call = 4;
    dxdt = derivs_x_wr_t_Rcpp(GrowthCrvChoice, function_call, x, y, params);
    F4 = dxdt;
    x = x + h * (F1 + 2.0 * (F2 + F3) + F4) / 6;
    t = Start + k * h;
  }

  return(x);
}


// [[Rcpp::export]]
double LenAtAge_Rcpp(const int j, NumericVector params, const int GrowthCrvChoice, const int nstep,
                     const int CalculationStage, const double LenPrevIntAge, const double StartAge,
                     NumericVector Obs_delta_t, NumericVector Obs_Initlen) {

  double t;
  double tAtEnd;
  double h;
  double EstLen;
  double x;
  x = 0;
  tAtEnd = 0;

  // Initial time
  t = StartAge;

  if (CalculationStage == 1) { // calculate final length, given time at liberty and initial length
    tAtEnd = Obs_delta_t(j-1);     // Final time
    x = Obs_Initlen(j-1);    // Initial size
  }

  if (CalculationStage == 2) { // calculate annual growth increment, given initial length
    tAtEnd = 365;
    x = j;    // Initial size
  }

  if (CalculationStage == 3) { // (generating growth curve) calculate annual growth increment, given initial length
    tAtEnd = 365;
    x = LenPrevIntAge;    // Initial size
  }

  // Determine the value of x() at t = tAtEnd
  h = (tAtEnd - StartAge) / nstep;
  x = RK4SYS_Rcpp(nstep, GrowthCrvChoice, h, t, x, params);
  EstLen = x;

  return(EstLen);

}

// [[Rcpp::export]]
NumericVector TaggingGrowthModelNLLCalcs_Rcpp(NumericVector params, const int nobs, const int GrowthCrvChoice,
                                              const int nstep, NumericVector Obs_delta_t, NumericVector Obs_Initlen) {

  double StartAge;
  double CalculationStage;
  double LenPrevIntAge;
  int j;
  NumericVector EstLenAtRelAge(nobs);
  NumericVector EstDeltaL(nobs);
  StartAge = 0;
  CalculationStage = 1;
  LenPrevIntAge = 0;
  for (j=1; j<=nobs; j++) {
    EstLenAtRelAge(j-1) = LenAtAge_Rcpp(j, params, GrowthCrvChoice, nstep, CalculationStage,
                                     LenPrevIntAge, StartAge, Obs_delta_t, Obs_Initlen);
  }
  // estimated growth, given parameters, and length of time at liberty
  EstDeltaL = EstLenAtRelAge - Obs_Initlen;

  return(EstDeltaL);
}






