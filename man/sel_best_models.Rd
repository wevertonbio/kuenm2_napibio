% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sel_best_models.R
\name{sel_best_models}
\alias{sel_best_models}
\title{Select best models}
\usage{
sel_best_models(
  calibration_results = NULL,
  cand_models = NULL,
  model_type = c("glmnet", "glm"),
  test_concave = TRUE,
  omrat_threshold = 10,
  allow_tolerance = TRUE,
  tolerance = 0.01,
  AIC = "ws",
  significance = 0.05,
  delta_aic = 2,
  verbose = TRUE
)
}
\arguments{
\item{calibration_results}{an object of class \code{calibration_results} returned
by the \code{\link{calibration}}() function. Default is NULL.}

\item{cand_models}{(data.frame) a summary of the evaluation metrics for each
candidate model. In the output of the \code{\link{calibration}}(), this
data.frame is located in \verb{$calibration_results$Summary}. Default is NULL.}

\item{model_type}{(character) model type, either "glm" or "glmnet".}

\item{test_concave}{(logical) whether to remove candidate models presenting
concave curves. Default is TRUE.}

\item{omrat_threshold}{(numeric) the maximum omission rate a candidate model
can have to be considered a best model. Default is 10. This value must match
one of the values specified in \code{omrat} in \code{\link{calibration}}().}

\item{allow_tolerance}{(logical) whether to allow selection of models with
minimum values of omission rates even if their omission rate surpasses the
\code{omrat_threshold}. This is only applicable if all candidate models have
omission rates higher than the \code{omrat_threshold}. Default is TRUE.}

\item{tolerance}{(numeric) The value added to the minimum omission rate if it
exceeds the \code{omrat_threshold}. If \code{allow_tolerance = TRUE}, selected models
will have an omission rate equal to or less than the minimum rate plus this
tolerance. Default is 0.01.}

\item{AIC}{(character) the type of AIC to be calculated: "ws" for AIC
proposed by Warren and Seifert (2011), or "nk" for AIC proposed by Ninomiya
and Kawano (2016). This is only applicable if model_type = "glmnet".
Default is "ws". See References for details.}

\item{significance}{(numeric) the significance level to select models
based on the partial ROC (pROC). Default is 0.05. See Details.}

\item{delta_aic}{(numeric) the value of delta AIC used as a threshold to
select models. Default is 2.}

\item{verbose}{(logical) whether to display messages during processing.
Default is TRUE.}
}
\value{
If calibration_results is provided, it returns a new calibration_results with
the new selected models and summary. If calibration_results is NULL, it
returns a list containing the following elements:
\itemize{
\item selected_models: data frame with the ID and the summary of evaluation
metrics for the selected models.
\item summary: A list containing the delta AIC values for model selection, and
the ID values of models that failed to fit, had concave curves,
non-significant pROC values, omission rates above the threshold, delta AIC
values above the threshold, and the selected models.
}
}
\description{
This function selects the best models according to user-defined criteria, evaluating statistical significance (partial ROC), predictive ability (omission rates), and model complexity (AIC).
}
\details{
Partial ROC is calculated following Peterson et al.
(2008; http://dx.doi.org/10.1016/j.ecolmodel.2007.11.008).
}
\examples{
# Import example of calibration results (output of calibration function)
## GLM
data("calib_results_glm", package = "kuenm2")

#Select new best models based on another value of omission rate
new_best_model <- sel_best_models(cand_models = calib_results_glm$calibration_results$Summary,
                                  model_type = "glm",
                                  test_concave = TRUE,
                                  omrat_threshold = 5,
                                  allow_tolerance = TRUE,
                                  tolerance = 0.01,
                                  AIC = "ws",
                                  significance = 0.05,
                                  delta_aic = 10, #Higher value of delta AIC
                                  verbose = TRUE)
#Compare with best models selected previously, with omission rate of 10 and delta AIC of 2
calib_results_glm$summary$Selected #Models 1, 2 and 5 selected
new_best_model$summary$Selected #Models 1 and 5 selected
#Replace selected models in calib_results
calib_results_glm$selected_models <- new_best_model$cand_final
calib_results_glm$summary <- new_best_model$summary

}
\references{
Ninomiya, Yoshiyuki, and Shuichi Kawano. "AIC for the Lasso in generalized
linear models." (2016): 2537-2560.

Warren, D. L., & Seifert, S. N. (2011). Ecological niche modeling in Maxent:
the importance of model complexity and the performance of model selection
criteria. Ecological applications, 21(2), 335-342.
}
