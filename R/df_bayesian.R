#' Bayesian log normal models
#'
#'
#' @export
#'
df_bayes_fix <- function(y, L_m, U_m, L_cv, U_cv, ...){

    if (!is.vector(y))
        stop("y must be a vector.")

    if (length(L_m)  != 1 |
        length(U_m)  != 1 |
        length(L_cv) != 1 |
        length(U_cv) != 1)
        stop("L or U must be a single number.")

    lst_data <- list(N    = length(y),
                     y    = y   / 100,
                     L_m  = L_m / 100,
                     U_m  = U_m / 100,
                     L_cv = L_cv,
                     U_cv = U_cv)

    fit <- df_stan(lst_data,
                   stan_mdl = "logn",
                   ...)

    post_par     = rstan::extract(fit)
    post_m       = post_par$m * 100
    post_cv      = post_par$cv
    post_y_tilde = post_par$y_tilde * 100

    rst <- list(post_m       = post_m,
                post_cv      = post_cv,
                post_y_tilde = post_y_tilde)

    return(rst)
}
