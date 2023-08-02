#' Get probabilities for intervals
#'
#' @export
#'
get_prob_intervals <- function(post_smps, interval) {
    interval <- c(-Inf, interval, Inf)
    rst      <- sapply(interval,
                       function(x) mean(post_smps <= x))

    sapply(2:length(rst), function(x) rst[x] - rst[x - 1])
}


#' Simulate patients
#'
#' @param cv Coefficient of variation. If it is scalar, the CV is assumed to be
#'     the same across dose levels. If it is a vector, it has to be the same
#'     length as mean_raw
#'
#' @export
#'
#'
generate_patients_cv <- function(mean_raw, cv, ...) {

    if (1 == length(cv)) {
        cv <- rep(cv, length(mean_raw))
    } else {
        stopifnot(length(cv) == length(mean_raw))
    }

    mean_log <- log(mean_raw / sqrt(1 + cv^2))
    sd_log   <- sqrt(log(1 + cv^2))

    generate_patients(mean_log,
                      rep(sd_log, length(mean_log)),
                      ...)
}

#' Simulate patients
#'
#' the lengths of dose_level, mean_log, sd_log and dlt_rate must be the same
#'
#' right now let n_patients to be the same across dose levels
#'
#' @export
#'
generate_patients <- function(mean_log   = c(2.88, 3.8, 4.27, 4.49),
                              sd_log     = rep(0.47, 4),
                              dlt_rate   = c(0.01, 0.01, 0.02, 0.02),
                              n_patients = 6) {

    stopifnot(length(mean_log) == length(sd_log) &
              length(mean_log) == length(dlt_rate))

    func_act <- mapply(rlnorm,
                       n       = n_patients,
                       meanlog = mean_log,
                       sdlog   = sd_log)

    dlt      <- mapply(rbinom,
                       n    = n_patients,
                       size = 1,
                       prob = dlt_rate)

    data     <- as.data.frame(cbind(reshape2::melt(func_act)[, 3],
                                    reshape2::melt(dlt)[, 3],
                                    rep(seq_len(length(mean_log)),
                                        each = n_patients)))

    colnames(data) <- c("functional_activity", "DLT", "dose_level")
    data
}

#' Design point
#'
#' Whether the dose is safe
#'
#' @export
#'
get_guarded <- function(dta_fix, dta_dlt,
                        fix_interval_ind = c(5, 150),
                        dlt_thresh_r     = 0.15,
                        dlt_thresh_c     = 0.8,
                        dlt_prior        = c(0.1, 0.9)) {


    ## bayesian DLT threshold
    post_a   <- dlt_prior[1] + sum(dta_dlt)     + 1
    post_b   <- dlt_prior[1] + sum(1 - dta_dlt) + 1
    p_dlt_gt <- 1 - pbeta(dlt_thresh_r, post_a, post_b)

    ## fix activity and dlt
    rst <- c(all(dta_fix < fix_interval_ind[2]),
             p_dlt_gt < dlt_thresh_c)

    rst
}


#' Dose escalation algorithm
#'
#' 0:stay
#'
#' -12:stop after 6 patients
#' -13:stop before 6 patients
#' -14:stop for both DLT and Over FIX
#' -15:stop for Over FIX
#' -16:stop for DLT
#'
#' 11: escalate after the 1st patient;
#' 12: escalate after 6 patients;
#' 13: escalate before 6 patients
#'
#' 2:selected
#'
#' @export
#'
dose_algorithm_1 <- function(cur_data, lst_design, ...) {

    n_pt_max          <- lst_design$n_pt_max
    n_pt_min          <- lst_design$n_pt_min

    dlt_thresh_r      <- lst_design$dlt_thresh_r
    dlt_thresh_c      <- lst_design$dlt_thresh_c
    dlt_prior         <- lst_design$dlt_prior

    fix_thresh        <- lst_design$fix_thresh
    fix_prior_meanraw <- lst_design$fix_prior_meanraw
    fix_prior_cv      <- lst_design$fix_prior_cv

    fix_interval_ind       <- lst_design$fix_interval_ind
    fix_interval_meanraw   <- lst_design$fix_interval_meanraw


    dta_fix <- cur_data[, 1]
    dta_dlt <- cur_data[, 2]
    npt      <- nrow(cur_data)
    res      <- rep(NA, npt)

    ## check the first patient
    if (dta_fix[1] < fix_interval_ind[1]) {
        res[1] <- 11
        return(res)
    } else {
        res[1] <- 0
    }

    ## enroll and treat the rest patients
    for (i in 2 : npt) {

        ## dont check until the first cohort is finished
        if (i < n_pt_min) {
            res[i] <- 0
            next
        }

        cumu_fix   <- dta_fix[1:i]
        cumu_dlt   <- dta_dlt[1:i]

        ## check dlt and over fix
        is_guarded <- get_guarded(cumu_fix, cumu_dlt,
                                  fix_interval_ind,
                                  dlt_thresh_r,
                                  dlt_thresh_c,
                                  dlt_prior)


        if (all(0 == is_guarded)) {
            res[i] <- -14
            break
        } else if (0 == is_guarded[1]) {
            res[i] <- -15
            break
        } else if (0 == is_guarded[2]) {
            res[i] <- -16
            break
        }

        ## check posterior fix
        bayes_samples <- df_bayes_fix(y    = cumu_fix,
                                      L_m  = fix_prior_meanraw[1],
                                      U_m  = fix_prior_meanraw[2],
                                      L_cv = fix_prior_cv[1],
                                      U_cv = fix_prior_cv[2],
                                      ...)

        pmu  <- get_prob_intervals(bayes_samples$post_m,
                                   fix_interval_meanraw)
        py   <- get_prob_intervals(bayes_samples$post_y_tilde,
                                   fix_interval_ind)

        max_mu      <- which.max(pmu)
        is_over_fix <- py[3] > fix_thresh

        if (n_pt_max == i) {
            if (!is_over_fix & 1 == max_mu) {
                res[i] <- 12
            } else if (!is_over_fix & 2 == max_mu) {
                res[i] <- 2
            } else {
                res[i] <- -12
            }
        } else {
            if (!is_over_fix & 1 == max_mu) {
                res[i] <- 13
                break
            } else if (is_over_fix & 3 == max_mu) {
                res[i] <- -13
                break
            } else {
                res[i] <- 0
            }
        }
    }

    return(res)
}


#' Enroll patients
#'
#'
#' @export
#'
enroll_patients <- function(dta_all,
                            fun_alg = dose_algorithm_1,
                            ...) {

    n_dose  <- max(dta_all$dose_level)
    rst     <- NULL
    for (i in 1:n_dose) {
        cur_data <- dta_all %>%
            dplyr::filter(dose_level == i)

        cur_rst  <- fun_alg(cur_data, ...)
        rst      <- c(rst, cur_rst)

        if (!any(10 < cur_rst, na.rm = TRUE))
            break
    }

    to_add <- nrow(dta_all) - length(rst)
    if (to_add > 0)
        rst <- c(rst, rep(NA, to_add))

    dta_all$decision <- rst
    dta_all
}
