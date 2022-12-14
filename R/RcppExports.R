# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

derivs_x_wr_t_Rcpp <- function(GrowthCrvChoice, function_call, x, y, params) {
    .Call(`_WAFishBiology_derivs_x_wr_t_Rcpp`, GrowthCrvChoice, function_call, x, y, params)
}

RK4SYS_Rcpp <- function(nstep, GrowthCrvChoice, h, t, x, params) {
    .Call(`_WAFishBiology_RK4SYS_Rcpp`, nstep, GrowthCrvChoice, h, t, x, params)
}

LenAtAge_Rcpp <- function(j, params, GrowthCrvChoice, nstep, CalculationStage, LenPrevIntAge, StartAge, Obs_delta_t, Obs_Initlen) {
    .Call(`_WAFishBiology_LenAtAge_Rcpp`, j, params, GrowthCrvChoice, nstep, CalculationStage, LenPrevIntAge, StartAge, Obs_delta_t, Obs_Initlen)
}

TaggingGrowthModelNLLCalcs_Rcpp <- function(params, nobs, GrowthCrvChoice, nstep, Obs_delta_t, Obs_Initlen) {
    .Call(`_WAFishBiology_TaggingGrowthModelNLLCalcs_Rcpp`, params, nobs, GrowthCrvChoice, nstep, Obs_delta_t, Obs_Initlen)
}

