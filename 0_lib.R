approx_ecdf <- function(df, grp_col, ptile_col) {
  
  ptile_tbl <- df %>%
    group_by({{grp_col}}) %>%
    summarise(
      tot = sum({{ptile_col}}),
      .groups = "drop"
    ) %>%
    arrange({{grp_col}}) %>%
    mutate(
      cumu_tot = cumsum(tot)
    ) %>%
    ungroup() %>%
    mutate(
      pct_tot = cumu_tot / sum(tot)
    )
  
  ptile_x_vals <- ptile_tbl %>% pull({{grp_col}})
  
  find_p_val_argmin <- function(p_val) {
    p_val_diff <- abs((p_val - ptile_tbl$pct_tot))
    ptile_x_vals[which((p_val_diff == min(p_val_diff)))]
  }
  
  function(ptile_cuts) {
    vec_vals <- vapply(ptile_cuts, find_p_val_argmin, c(ptile_x_vals[1]))
    names(vec_vals) <- ptile_cuts
    vec_vals
  }
}

make_bit_str_to_formula <- function(formula_template, spline_locations) {
  
  function(bit_string) {
    
    bit_str_iter <- bit_string
    full_formula <- formula_template
    
    for (formula_i in 1:length(spline_locations)) {
      
      curr_part <- names(spline_locations)[formula_i]
      all_knots <- spline_locations[[curr_part]]
      curr_knot_locs <- which(bit_str_iter[1:length(all_knots)] == 1)
      curr_knots <- all_knots[curr_knot_locs]
      
      curr_knot_str <- paste0(curr_knots, collapse=",")
      curr_formula <- stringr::str_replace(curr_part, "%KNOTS%", curr_knot_str)
      full_formula <- stringr::str_replace_all(full_formula, fixed(curr_part), curr_formula)
      
      # shift bit string
      bit_str_iter <- bit_str_iter[-(1:length(all_knots))]
    }
    
    as.formula(full_formula)
  }
}

calc_aic_bic_glmnet <- carrier::crate(
  function(glmnet_fit) {

    tLL <- glmnet_fit$nulldev - stats::deviance(glmnet_fit)
    k <- glmnet_fit$df
    n <- glmnet_fit$nobs
    AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
    
    BIC <- log(n)*k - tLL
    
    list(
      "AIC" = AICc,
      "BIC" = BIC
    )
  })

log_interp_spline <- function(duration, interp_start = 10, interp_end=15, interp_pwr=0.33) {
  
  duration <- pmin(interp_end, duration)
  
  c1 <- ((interp_end - duration) / (interp_end - interp_start))^interp_pwr
  c1 <- pmin(1, pmax(0, c1))
                     
  c1 * log(duration) + (1 - c1) * log(interp_end)
}
