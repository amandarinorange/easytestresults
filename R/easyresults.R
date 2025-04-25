####################################################
#### easytestresults
####################################################
#' easytestresults package
#'
#' This package includes two functions for outputting A/B/n test results by experimental condition, as well as effects of condition within levels of a given user segment
#'
#'
#' @section function easy_results()
#'
#' @describeIn easy_results Reports pairwise results from an A/B/n test where each treatment condition is contrasted with the reference group
#'
#'
#' @param outcome The name of the column holding the dependent variable (aka y variable/outcome), specified inside quotation marks
#' @param treatment_name The name of the column holding the experimental treatment name, specified inside quotation marks
#' @param df The name of the dataframe, unquoted
#' @param user_segment The name of the column holding the user segment, specified inside quotation marks
#' @param family The type of test to be performed. Currently, support is only provided for the following: lm (OLS regression), logit (logistic regression), negbin (negative binomial regression), and zinb (zero-inflated negative binomial). However, confidence intervals are not output for the latter family type.
#' @import lmtest
#' @import tibble
#' @import dplyr
#' @importFrom MASS "glm.nb"
#' @import stats
#' @import emmeans
#' @import interactions
#' @import tidyr
#' @import purrr
#' @importFrom emdbook "dzinbinom"
#' @return A table of test statistics reporting the results of pairwise comparisons for each
#' treatment group versus the reference group. When user_segment is specified, the simple effects of treatment condition
#' are reported within each level of the user segment variable. Test statistics include point estimate, effective sample size (using pairwise deletion for missing values), confidence interval,
#' relative lift (percent), z- or t-statistic, and p-value. Alpha is assumed to be set at the conventional 0.05.
#'
#'
#' @export
easy_results <- function(outcome, treatment_name, df, family) {
  treatment_col_ind <- which(colnames(df) == treatment_name)
  outcome_col_ind <- which(colnames(df) == outcome)

  #if (!is.factor(df[,treatment_col_ind])) {
  #  stop("x variables must be coded as a factor")
  #}

  if (family == "logit") {
    if(!is.numeric(df[,outcome_col_ind])) {
      stop("Outcome variable must be numeric (0 or 1)")
    }
    fit_stats <- data.frame(treatment_name = levels(as.factor(df[,treatment_name])),est = NA, test_stat = NA, p = NA, lift = NA)
    f <- as.formula(paste(outcome, treatment_name, sep="~"))
    fit <- do.call("glm", list(formula=f, data=df, family="binomial"))
    fit_summary <- summary(fit)
    fit_coeffs <- rownames_to_column(as.data.frame(fit_summary$coefficients), var = "term")
    #fit_coeffs$`Pr(>|z|)`[fit_coeffs$`Pr(>|z|)` < .001] <- '<.001'

    fit_stats$est[1] <- round(plogis(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"]), 5)

    for (i in 2:nrow(fit_stats)) {
      fit_stats$est[i] <- round(plogis(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"] + fit_coeffs$Estimate[i]), 5)
      fit_stats$test_stat[i] <- round(fit_coeffs$`z value`[i], 3)
      fit_stats$p[i] <- fit_coeffs$`Pr(>|z|)`[i]
      fit_stats$lift[i] <- paste(round(100* (fit_stats$est[i] - fit_stats$est[1]) / fit_stats$est[1], 3),'%')
    }
  }

  if (family == "lm") {
    fit_stats <- data.frame(treatment_name = levels(as.factor(df[,treatment_name])),est = NA, test_stat = NA, p = NA, lift = NA)
    f <- as.formula(paste(outcome, treatment_name, sep="~"))
    fit <- do.call("lm", list(formula=f, data=df))
    fit_summary <- summary(fit)
    fit_coeffs <- rownames_to_column(as.data.frame(fit_summary$coefficients), var = "term")
    #fit_coeffs$`Pr(>|z|)`[fit_coeffs$`Pr(>|z|)` < .001] <- '<.001'

    fit_stats$est[1] <- round(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"], 3)

    for (i in 2:nrow(fit_stats)) {
      fit_stats$est[i] <- round(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"] + fit_coeffs$Estimate[i], 3)
      fit_stats$test_stat[i] <- round(fit_coeffs$`t value`[i],3)
      fit_stats$p[i] <- fit_coeffs$`Pr(>|t|)`[i]
      fit_stats$lift[i] <- paste(round(100* (fit_stats$est[i] - fit_stats$est[1]) / fit_stats$est[1], 3),'%')
    }
  }


  if (family == "negbin") {
    fit_stats <- data.frame(treatment_name = levels(as.factor(df[,treatment_name])),est = NA, test_stat = NA, p = NA, lift = NA)
    f <- as.formula(paste(outcome, treatment_name, sep="~"))
    fit <- do.call("glm.nb", list(formula=f, data=df))
    fit_summary <- summary(fit)
    fit_coeffs <- rownames_to_column(as.data.frame(fit_summary$coefficients), var = "term")
    #fit_coeffs$`Pr(>|z|)`[fit_coeffs$`Pr(>|z|)` < .001] <- '<.001'

    fit_stats$est[1] <- round(exp(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"]), 3)

    for (i in 2:nrow(fit_stats)) {
      fit_stats$est[i] <- round(exp(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"] + fit_coeffs$Estimate[i]), 3)
      fit_stats$test_stat[i] <- round(fit_coeffs$`z value`[i], 3)
      fit_stats$p[i] <- fit_coeffs$`Pr(>|z|)`[i]
      fit_stats$lift[i] <- paste(round(100* (fit_stats$est[i] - fit_stats$est[1]) / fit_stats$est[1], 3),'%')
    }
  }

  if (family == "zinb") {
    fit_stats <- data.frame(treatment_name = levels(as.factor(df[,treatment_name])),est = NA, test_stat = NA, p = NA, lift = NA)
    f <- as.formula(paste(outcome, treatment_name, sep="~"))
    fit <- do.call("zeroinfl", list(formula=f, data=df, dist = "negbin"))
    fit_summary <- summary(fit)
    fit_coeffs <- rownames_to_column(as.data.frame(fit_summary$coefficients$count), var = "term")
    #fit_coeffs$`Pr(>|z|)`[fit_coeffs$`Pr(>|z|)` < .001] <- '<.001'

    fit_stats$est[1] <- round(exp(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"]), 3)

    for (i in 2:nrow(fit_stats)) {
      fit_stats$est[i] <- round(exp(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"] + fit_coeffs$Estimate[i]), 3)
      fit_stats$test_stat[i] <- round(fit_coeffs$`z value`[i], 3)
      fit_stats$p[i] <- fit_coeffs$`Pr(>|z|)`[i]
      fit_stats$lift[i] <- paste(round(100* (fit_stats$est[i] - fit_stats$est[1]) / fit_stats$est[1], 3),'%')
    }
  }

  `%notin%` <- Negate(`%in%`)
  #if (family %notin% c("logit","lm","negbin")) {
  #  stop("Test family should be one of: logit, lm, negbin")
  #}

  fit_stats <- fit_stats %>%
    pivot_wider(names_from = 1, values_from = c(2,3,4,5)) %>%
    setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))

  # Get confidence intervals
  if (family %in% c("logit","lm","negbin")) {
    fit_data <- suppressWarnings(do.call("cat_plot", list(model = fit, pred = as.name(treatment_name)))$data)

    fit_cis <- fit_data %>%
      dplyr::mutate(est_ci = gsub(x = (paste('[',round(ymin,4),', ',round(ymax,4),']')), pattern = " ", replacement = "")) %>%
      dplyr::select(c(2,"est_ci")) %>%
      pivot_wider(names_from = 1, values_from = c("est_ci"))

    colnames(fit_cis) <- paste(levels(df[,treatment_col_ind]),"estci",sep="_")


    fit_final <- fit_stats %>%
      cbind(fit_cis)
    fit_final <- fit_final %>%
      dplyr::select(c(1,order(colnames(fit_final)))) %>%
      discard(~all(is.na(.) | . ==""))


    colnames(fit_final) <- gsub(" ", "", colnames(fit_final))
  }

  if (family %in% c("zinb")) {
    fit_data <- suppressWarnings(do.call("confint", list(object = fit)))


    fit_cis <- fit_data %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      dplyr::filter(str_starts(rowname, "count")) %>%
      dplyr::mutate(term = rowname %>%
                      str_remove_all("count_") %>%
                      str_remove_all(treatment_name),
                    est_ci = paste0("[", round(`2.5 %`, 4), ",", round(`97.5 %`, 4), "]") %>%
                      gsub(" ", "", .)) %>%
      dplyr::select(term, est_ci) %>%
      pivot_wider(names_from = term, values_from = est_ci)

    colnames(fit_cis) <- paste(levels(df[,treatment_col_ind]),"estci",sep="_")


    fit_final <- fit_stats %>%
      cbind(fit_cis)
    fit_final <- fit_final %>%
      dplyr::select(c(1,order(colnames(fit_final)))) %>%
      discard(~all(is.na(.) | . ==""))


    colnames(fit_final) <- gsub(" ", "", colnames(fit_final))
  }
  # sample sizes per cell
  sample_sizes <- table(df[,treatment_col_ind][!is.na(df[outcome])]) %>%
    data.frame() %>%
    dplyr::rename(treatment_name = Var1,
                  n = Freq) %>%
    pivot_wider(names_from = treatment_name,
                values_from = n,
                names_glue = "{treatment_name}_n")

  # Final object is one row with all results from the model
  return(as.data.frame(cbind(outcome_variable = outcome, model_type = family, sample_sizes, fit_final)))
}

#' @section function easy_results_segmented()
#'
#' @describeIn easy_results_segmented Reports simple effects of the treatment condition within each level of a categorical third variable, typically a user segment value
#'
#' @export
easy_results_segmented <- function(outcome, treatment_name, user_segment, df, family) {

  treatment_col_ind <- which(colnames(df) == treatment_name)
  segment_col_ind <- which(colnames(df) == user_segment)
  outcome_col_ind <- which(colnames(df) == outcome)

  if (!is.factor(df[,segment_col_ind])) {
    stop("User segment variable must be coded as factor")
  }
  if (!is.numeric(df[,outcome_col_ind]) & family == "logit") {
    stop("Outcome variable must be numeric (0 or 1) for logistic regression")
  }
  if (sum(df[,outcome_col_ind] - floor(df[,outcome_col_ind])!=0)!=0 & family == "negbin") {
    stop("Outcome variable must be integer for negative binomial regression")
  }
  intx <- paste(treatment_name, user_segment, sep = "*")
  f_intx <- as.formula(paste(outcome, intx, sep="~"))
  cov <- paste(treatment_name, user_segment, sep = "+")
  f_cov <- as.formula(paste(outcome, cov, sep="~"))

  if (family == "logit") {
    fit_cov <- do.call("glm", list(formula=f_cov, data=df, family="binomial"))
    fit_intx <- do.call("glm", list(formula=f_intx, data=df, family="binomial"))
  }

  if (family == "lm") {
    fit_cov <- do.call("lm", list(formula=f_cov, data=df))
    fit_intx <- do.call("lm", list(formula=f_intx, data=df))
  }

  if (family == "negbin") {
    fit_cov <- do.call("glm.nb", list(formula=f_cov, data=df))
    fit_intx <- do.call("glm.nb", list(formula=f_intx, data=df))
  }

  # if (family == "zip") {
  #   fit_cov <- do.call("zeroinfl", list(formula=f_cov, data=df, dist = "poisson"))
  #   fit_intx <- do.call("zeroinfl", list(formula=f_intx, data=df, dist = "poisson"))
  # }
  #
  #
  if (family == "zinb") {
    fit_cov <- do.call("zeroinfl", list(formula=f_cov, data=df, dist = "negbin"))
    fit_intx <- do.call("zeroinfl", list(formula=f_intx, data=df, dist = "negbin"))
  }
  `%notin%` <- Negate(`%in%`)

  # Check for significance of interaction effect
  lrtest_res <- lrtest(fit_cov,fit_intx)

  if (family %in% c("logit","lm","negbin")) {
    intx_plot <- do.call("cat_plot", list(model = fit_intx,
                                          pred = as.name(treatment_name),
                                          modx = as.name(user_segment),
                                          colors = Polychrome::createPalette(N = length(levels(df[,segment_col_ind])) * length(levels(df[,treatment_col_ind])), seedcolors = c("#ff0000", "#00ff00", "#0000ff"))))

    # Get confidence intervals
    intx_cis <- intx_plot$data %>%
      dplyr::mutate(est_ci = gsub(x = (paste('[',round(ymin,4),', ',round(ymax,4),']')), pattern = " ", replacement = "")) %>%
      dplyr::select(c(2,3,"est_ci")) %>%
      pivot_wider(names_from = 2, values_from = c("est_ci")) %>%
      dplyr::select(-1)

    colnames(intx_cis) <- paste(colnames(intx_cis),"estci",sep="_")

  }

  if (family %in% c("zinb")) {
    # confidence intervals not yet supported for zinb models, so we'll just leave them blank

    treatments <- unique(df[[treatment_col_ind]][!is.na(df[[outcome]])])
    segments <- unique(df[[segment_col_ind]][!is.na(df[[outcome]])])

    # Create grid of all combinations
    grid_df <- expand.grid(
      segment = segments,
      treatment = treatments,
      stringsAsFactors = FALSE
    )

    # Add an NA column for each treatment level
    intx_cis <- tidyr::pivot_wider(
      grid_df,
      names_from = treatment,
      values_from = treatment,  # temporary value
      values_fn = function(x) NA,
      values_fill = NA) %>%
      dplyr::select(-1)

    colnames(intx_cis) <- paste(levels(df[,treatment_col_ind]),"estci",sep="_")

    # Point estimates -- cat_plot doesn't support zeroinfl models, so I'm running glm.nb on the same data just to get point estimates and relative lifts.
    ## reset fit_intx into a glm.nb object instead of zeroinfl
    fit_intx <- do.call("glm.nb", list(formula=f_intx, data=df))
    intx_plot <- do.call("cat_plot", list(model = fit_intx,
                                          pred = as.name(treatment_name),
                                          modx = as.name(user_segment),
                                          colors = Polychrome::createPalette(N = length(levels(df[,segment_col_ind])) * length(levels(df[,treatment_col_ind])), seedcolors = c("#ff0000", "#00ff00", "#0000ff"))))
  }

  # Relative lift
  intx_data <- intx_plot$data %>%
    dplyr::select(c(1,2,3)) %>%
    pivot_wider(names_from = 3, values_from = c(1))

  intx_lift <- vector(length = nrow(intx_data))

  for (i in 3:ncol(intx_data)) {
    lift <- round((100 * (intx_data[,i] - intx_data[,2]) / intx_data[,2]), 3)
    intx_lift <- cbind(intx_lift, lift)
  }

  intx_lift <- intx_lift %>%
    dplyr::select(-c("intx_lift"))

  colnames(intx_lift) <- paste(colnames(intx_lift),"lift_percent",sep="_")

  intx_plot$data[,1] <- round(intx_plot$data[,1], 5)

  # Get coefficients, z-scores, and p-values
  em <- emmeans(fit_intx, as.formula(paste("~", intx)))
  em_c <- as.data.frame(contrast(em, method = "revpairwise",
                                 by = user_segment,
                                 adjust = "none",
                                 type = "response",
                                 ratios = F)) %>%
    dplyr::select(-c("estimate","SE","df"))

  em_c <- separate(data = em_c, col = contrast, into = c("reference", "comparison"), sep = " - ", remove = F)

  em_c_comparison <- em_c[em_c$comparison == levels(df[,treatment_name])[1],] %>%
    dplyr::rename_with(~ "test_stat", .cols = dplyr::ends_with("ratio"))

  intx_stats <- intx_plot$data %>%
    dplyr::select(c(1,2,3)) %>%
    dplyr::rename(est = 1) %>%
    left_join(y = em_c_comparison,
              by = setNames(c(user_segment,"reference"),
                            c(user_segment,treatment_name))) %>%
    dplyr::select(-c("contrast","comparison")) %>%
    pivot_wider(names_from = 3, values_from = c(1,4,5)) %>%
    setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.))) %>%
    cbind(intx_lift,intx_cis)

  colnames(intx_stats) <- gsub(" ", "", colnames(intx_stats))

  intx_final <- intx_stats %>%
    dplyr::select(c(1,order(colnames(intx_stats)))) #%>%
  #discard(~all(is.na(.) | . =="")) #discard breaks the output for zinb models?

  # sample sizes per cell

  sample_sizes <- table(df[,treatment_col_ind][!is.na(df[outcome])], df[,segment_col_ind][!is.na(df[outcome])]) %>%
    data.frame() %>%
    dplyr::rename(treatment_name = Var1,
                  user_segment = Var2,
                  n = Freq) %>%
    pivot_wider(names_from = c(treatment_name),
                values_from = n,
                names_glue = "{treatment_name}_n")


  if (lrtest_res$`Pr(>Chisq)`[2] >= .05) {
    interaction_significance <- paste("Warning! The interaction effect '",
                                      intx,
                                      "' is NOT SIGNIFICANT at p = ",
                                      round(lrtest_res$`Pr(>Chisq)`[2], 3),
                                      ". Beware with interpreting results of treatment by ",
                                      user_segment)
    return(cbind(outcome_variable = outcome, model_type = family, sample_sizes, intx_final, interaction_significance))
  }

  if (lrtest_res$`Pr(>Chisq)`[2] < .05) {
    interaction_significance <- paste("The interaction effect '",
                                      intx,
                                      "' is significant at p = ",
                                      ifelse(test = lrtest_res$`Pr(>Chisq)`[2] < .001, yes = '<.001', no = round(lrtest_res$`Pr(>Chisq)`[2], 3)))


    return(cbind(outcome_variable = outcome, model_type = family, sample_sizes, intx_final, interaction_significance))
  }

}

#' @section function easy_main_effects()
#'
#' @describeIn easy_main_effects Reports main effects of two factors on a binary, normal, or negative binomial outcome variable
#'
#' @export
easy_main_effects <- function(outcome, me1, me2, df, family) {
  me1_col_ind <- which(colnames(df) == me1)
  me2_col_ind <- which(colnames(df) == me2)
  outcome_col_ind <- which(colnames(df) == outcome)

  cov <- paste(me1, me2, sep = "+")
  f <- as.formula(paste(outcome, cov, sep="~"))
  intx <- paste(me1, me2, sep = "*")
  f_intx <- as.formula(paste(outcome, intx, sep="~"))

  #if (!is.factor(df[,treatment_col_ind])) {
  #  stop("x variables must be coded as a factor")
  #}

  if (family == "logit") {
    if(!is.numeric(df[,outcome_col_ind])) {
      stop("Outcome variable must be numeric (0 or 1)")
    }
    fit_stats <- rbind(data.frame(effect = 'control',est = NA, est_ci = NA, test_stat = NA, p = NA, lift = NA),
                       data.frame(effect = me1,est = NA,est_ci = NA, test_stat = NA, p = NA, lift = NA),
                       data.frame(effect = me2,est = NA, est_ci = NA,test_stat = NA, p = NA, lift = NA))

    fit <- do.call("glm", list(formula=f, data=df, family="binomial"))
    fit_intx <- do.call("glm", list(formula=f_intx, data=df, family="binomial"))
    fit_summary <- summary(fit)
    fit_coeffs <- rownames_to_column(as.data.frame(fit_summary$coefficients), var = "term")
    fit_cis <- as.data.frame(confint(fit))
    #fit_coeffs$`Pr(>|z|)`[fit_coeffs$`Pr(>|z|)` < .001] <- '<.001'

    fit_stats$est[1] <- round(plogis(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"]), 5)
    fit_stats$est_ci[1] <- paste('[',round(plogis(fit_cis$`2.5 %`[1]), 5), ',', round(plogis(fit_cis$`97.5 %`[1]), 5), ']')

    for (i in 2:nrow(fit_stats)) {
      fit_stats$est[i] <- round(plogis(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"] + fit_coeffs$Estimate[i]), 5)
      fit_stats$est_ci[i] <- paste('[',round(plogis(fit_cis$`2.5 %`[1] + fit_cis$`2.5 %`[i]), 5), ',', round(plogis(fit_cis$`2.5 %`[1] + fit_cis$`97.5 %`[i]), 5), ']')
      fit_stats$test_stat[i] <- round(fit_coeffs$`z value`[i], 3)
      fit_stats$p[i] <- fit_coeffs$`Pr(>|z|)`[i]
      fit_stats$lift[i] <- paste(round(100* (fit_stats$est[i] - fit_stats$est[1]) / fit_stats$est[1], 3),'%')
    }
  }

  if (family == "lm") {
    fit_stats <- rbind(data.frame(effect = 'control',est = NA,est_ci=NA, test_stat = NA, p = NA, lift = NA),
                       data.frame(effect = me1,est = NA,est_ci = NA, test_stat = NA, p = NA, lift = NA),
                       data.frame(effect = me2,est = NA,est_ci = NA, test_stat = NA, p = NA, lift = NA))

    fit <- do.call("lm", list(formula=f, data=df))
    fit_intx <- do.call("lm", list(formula=f_intx, data=df))
    fit_summary <- summary(fit)
    fit_coeffs <- rownames_to_column(as.data.frame(fit_summary$coefficients), var = "term")
    fit_cis <- as.data.frame(confint(fit))
    #fit_coeffs$`Pr(>|z|)`[fit_coeffs$`Pr(>|z|)` < .001] <- '<.001'

    fit_stats$est[1] <- round(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"], 5)
    fit_stats$est_ci[1] <- paste('[',round((fit_cis$`2.5 %`[1]), 5), ',', round((fit_cis$`97.5 %`[1]), 5), ']')

    for (i in 2:nrow(fit_stats)) {
      fit_stats$est[i] <- round(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"] + fit_coeffs$Estimate[i], 5)
      fit_stats$est_ci[i] <- paste('[',round((fit_cis$`2.5 %`[1] + fit_cis$`2.5 %`[i]), 5), ',', round((fit_cis$`2.5 %`[1] + fit_cis$`97.5 %`[i]), 5), ']')
      fit_stats$test_stat[i] <- round(fit_coeffs$`t value`[i],5)
      fit_stats$p[i] <- fit_coeffs$`Pr(>|t|)`[i]
      fit_stats$lift[i] <- paste(round(100* (fit_stats$est[i] - fit_stats$est[1]) / fit_stats$est[1], 5),'%')
    }
  }


  if (family == "negbin") {
    fit_stats <- rbind(data.frame(effect = 'control',est = NA,est_ci = NA, test_stat = NA, p = NA, lift = NA),
                       data.frame(effect = me1,est = NA, est_ci = NA, test_stat = NA, p = NA, lift = NA),
                       data.frame(effect = me2,est = NA, est_ci = NA, test_stat = NA, p = NA, lift = NA))


    fit <- do.call("glm.nb", list(formula=f, data=df))
    fit_intx <- do.call("glm.nb", list(formula=f_intx, data=df))
    fit_summary <- summary(fit)
    fit_coeffs <- rownames_to_column(as.data.frame(fit_summary$coefficients), var = "term")
    fit_cis <- as.data.frame(confint(fit))
    #fit_coeffs$`Pr(>|z|)`[fit_coeffs$`Pr(>|z|)` < .001] <- '<.001'

    fit_stats$est[1] <- round(exp(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"]), 3)
    fit_stats$est_ci[1] <- paste('[',round(exp(fit_cis$`2.5 %`[1]), 5), ',', round(exp(fit_cis$`97.5 %`[1]), 5), ']')

    for (i in 2:nrow(fit_stats)) {
      fit_stats$est[i] <- round(exp(fit_coeffs$Estimate[fit_coeffs$term == "(Intercept)"] + fit_coeffs$Estimate[i]), 3)
      fit_stats$est_ci[i] <- paste('[',round(exp(fit_cis$`2.5 %`[1] + fit_cis$`2.5 %`[i]), 5), ',', round(exp(fit_cis$`2.5 %`[1] + fit_cis$`97.5 %`[i]), 5), ']')
      fit_stats$test_stat[i] <- round(fit_coeffs$`z value`[i], 3)
      fit_stats$p[i] <- fit_coeffs$`Pr(>|z|)`[i]
      fit_stats$lift[i] <- paste(round(100* (fit_stats$est[i] - fit_stats$est[1]) / fit_stats$est[1], 3),'%')
    }
  }

  `%notin%` <- Negate(`%in%`)
  if (family %notin% c("logit","lm","negbin")) {
    stop("Test family should be one of: logit, lm, negbin")
  }


  fit_final <- fit_stats %>%
    dplyr::select(c(1,order(colnames(fit_stats)))) %>%
    discard(~all(is.na(.) | . ==""))


  lrtest_res <- lrtest(fit,fit_intx)

  colnames(fit_final) <- gsub(" ", "", colnames(fit_final))

  return(as.data.frame(cbind(outcome_variable = outcome, model_type = paste(family, ' main effects'), fit_final, intx_p = lrtest_res$`Pr(>Chisq)`[2])))
}
