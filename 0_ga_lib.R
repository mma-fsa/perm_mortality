ga_knot_search <- function(base_formula, knot_formula, offset=NULL, 
                           eval_metric="BIC", weight_callback=NULL, ...) {
  
  df = as.data.frame(as.list(parent.frame()))
  
  knot_definitions <- list(...)
  
  base_formula_c <- as.character(base_formula)
  
  if (length(base_formula_c) != 3) {
    stop("Expected base formula to be 'y ~ x', got '~ x' instead" )
  } 
  
  if (is.null(weight_callback)) {
    weight_callback <- carrier::crate(
      function(param_names) {
        data.frame(
          term_name = param_names,
          penalty = rep(1, length(param_names)),
          upper = rep(Inf, length(param_names)),
          lower = rep(-Inf, length(param_names))
        )
      })
  }
  
  bitstr_formula_base <- base_formula_c[[3]]
  
  total_bits <- 0
  for (knot_name in names(knot_definitions)) {
    if (stringr::str_trim(knot_name) != "") {
      if (!stringr::str_detect(bitstr_formula_base, pattern=fixed(knot_name))) {
        stop(sprintf("No spline placeholder called '%s'", knot_name))
      }
      this_spline <- knot_definitions[[knot_name]]$formula_def
      total_bits <- total_bits + knot_definitions[[knot_name]]$n_bits
      this_spline <- stringr::str_replace(
        this_spline,
        pattern = fixed("%KNOTS%"),
        replacement = paste0("%", knot_name, "_KNOTS%")
      )
      bitstr_formula_base <- stringr::str_replace_all(
        bitstr_formula_base, knot_name, this_spline
      )
    }  
  }
  
  bitstr_to_formula <- carrier::crate(
    function(bitstr) {
      
      fixed <- stringr::fixed
      
      bitstr_proc <- bitstr
      bitstr_formula <- bitstr_formula_base
      
      for (knot_name in names(knot_definitions)) {
        
        iter_knot_def <- knot_definitions[[knot_name]]
        iter_n_bits <- iter_knot_def$n_bits
        iter_knot_locs <- iter_knot_def$knots
        iter_bits <- bitstr_proc[1:iter_n_bits]
        iter_bit_knots_selected <- iter_knot_locs[which(iter_bits == 1)]
        iter_knot_name_pattern <- fixed(paste0("%", knot_name, "_KNOTS%"))
        
        bitstr_formula <- stringr::str_replace_all(
          bitstr_formula,
          iter_knot_name_pattern,
          paste(iter_bit_knots_selected, collapse = ",")
        )
        
        bitstr_proc <- bitstr_proc[-c(1:iter_n_bits)]
      }
      stats::as.formula(paste0(base_formula_c[[2]], "~", bitstr_formula))
    },
    knot_definitions = knot_definitions,
    bitstr_formula_base = bitstr_formula_base,
    base_formula_c=base_formula_c
  )
    
  eval_bitstr <- carrier::crate(
    
    function(bitstr) {
      
      model_formula <- bitstr_to_formula(bitstr)
      
      X_mat <- stats::model.matrix(
        model_formula,
        data = df
      )
      
      df_wts <- weight_callback(colnames(X_mat))
      
      y_vec <- df[[base_formula_c[[2]]]]
      if (!is.null(offset)) {
        ga_fit <- glmnet::glmnet(
          X_mat,
          df$number_of_deaths,
          offset = offset,
          family = "poisson",
          intercept = T,
          lower.limits = df_wts$lower,
          upper.limits = df_wts$upper,
          penalty.factor = df_wts$penalty
        )  
      } else {
        ga_fit <- glmnet::glmnet(
          X_mat,
          df$number_of_deaths,
          family = "poisson",
          intercept = T,
          lower.limits = df_wts$lower,
          upper.limits = df_wts$upper,
          penalty.factor = df_wts$penalty
        )
      }
      
      aic_bic <- calc_aic_bic_glmnet(ga_fit)
      -1*min(aic_bic[[eval_metric]])
    },
    bitstr_to_formula = bitstr_to_formula,
    df=df,
    weight_callback=weight_callback,
    base_formula_c=base_formula_c,
    offset=offset,
    calc_aic_bic_glmnet=calc_aic_bic_glmnet,
    eval_metric=eval_metric
  )
    
  # provide two initial population members:
  # an all linear model (no knots) and one
  # that is overfit (all knots)
  initial_suggestions <- rbind(
    rep(0, total_bits),
    rep(1, total_bits)
  )
  
  list(
    "knot_definitions" = knot_definitions,
    "eval_bitstr" = memoise::memoise(eval_bitstr),
    "bitstr_to_formula" = bitstr_to_formula,
    "initial_suggestions" = initial_suggestions,
    "total_bits" = total_bits
  )
  
}

ga_knot_def_inner <- function(col_name, col_v, knot_locs, fmt_str) {
  
  boundary_low <- min(col_v)
  boundary_high <- max(col_v)
  
  boundary_locs <- c(boundary_low, boundary_high)
  knot_locs <- lubridate::setdiff(unique(knot_locs), boundary_locs)
  
  formula_def <- sprintf(
    fmt_str,
    deparse(col_name),
    boundary_low,
    boundary_high
  )
  
  list(
    "formula_def" = formula_def,
    "end_knots" = c(boundary_low, boundary_high),
    "knots" = knot_locs,
    "n_bits" = length(knot_locs)
  )
  
}

ga_knot_def_ns <- function(colname, knot_locs) {
  colname <- substitute(colname)
  col_v <- eval(colname, env = parent.frame())
  fmt_str <- "splines::ns(%s, Boundary.knots = c(%.0f, %.0f), knots=c(%%KNOTS%%))"
  ga_knot_def_inner(colname, col_v, knot_locs, fmt_str)
}

ga_knot_def_bs <- function(colname, knot_locs, degree=1) {
  colname <- substitute(colname)
  col_v <- eval(colname, env = parent.frame())
  fmt_str <- paste0(
    "splines::bs(%s, Boundary.knots = c(%.0f, %.0f), knots=c(%%KNOTS%%), degree=", 
    as.character(degree), ")")
  ga_knot_def_inner(colname, col_v, knot_locs, fmt_str)
}


