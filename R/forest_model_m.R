#'
#'  Creates data for forestplot, 99 pct of code from the forestmodel package
#'
#' @description This is used to make data for plots
#'
#' @details changed some things in the broom coef where it won't correctly give values
#' 
#' @param model model
#' 
#' @return data for plot
#' 
#' @examples 
#' # something
#' 
#' #forest_model_m(model = fit_mul, return_data = T, factor_separate_line = T, exponentiate = T)
#' 
#' @export
#' 
#' @importFrom stats var runif rexp
#' @importFrom Rdpack reprompt
#' @import forestmodel


forest_model_m <- function(model,
                           panels = default_forest_panels(model, factor_separate_line = factor_separate_line),
                           covariates = NULL, exponentiate = NULL, funcs = NULL,
                           factor_separate_line = FALSE,
                           format_options = forest_model_format_options(),
                           theme = theme_forest(),
                           limits = NULL, breaks = NULL, return_data = FALSE,
                           recalculate_width = TRUE, recalculate_height = TRUE,
                           model_list = NULL, merge_models = FALSE, exclude_infinite_cis = TRUE) {
  
  mapping <- aes(estimate, xmin = conf.low, xmax = conf.high)
  if (!is.null(model_list)) {
    if (!is.list(model_list)) {
      stop("`model_list` must be a list if provided.")
    }
    if (is.null(names(model_list))) {
      model_names <- rep("", length(model_list))
    } else {
      model_names <- names(model_list)
    }
    if (any(model_names == "")) {
      need_names <- which(model_names == "")
      model_names_needed <- vapply(model_list[need_names], function(x) quo_name(x$call), character(1))
      model_names[need_names] <- model_names_needed
    }
    if (!merge_models) {
      mapping <- c(mapping, aes(section = model_name))
    }
    if (is.null(exponentiate)) {
      exponentiate <- inherits(model_list[[1]], "coxph") ||
        (inherits(model_list[[1]], "glm") && model_list[[1]]$family$link == "logit")
    }
    if (missing(panels)) {
      panels <- default_forest_panels(model_list[[1]], factor_separate_line = factor_separate_line)
    }
  } else {
    if (is.null(exponentiate)) {
      exponentiate <- inherits(model, "coxph") ||
        (inherits(model, "glm") && model$family$link == "logit")
    }
  }
  
  if (exponentiate) trans <- exp else trans <- I
  
  stopifnot(is.list(panels))
  
  remove_backticks <- function(x) {
    gsub("^`|`$|\\\\(?=`)|`(?=:)|(?<=:)`", "", x, perl = TRUE)
  }
  
  make_forest_terms <- function(model) {
    tidy_model <- broom::tidy(model, conf.int = FALSE,conf.type="Wald")
    tidy_model <-  tidy_model %>% 
      mutate(conf.low=estimate-(1.96*std.error),
             conf.high=estimate+(1.96*std.error),
             estimate = if_else(estimate>5 | -5 > estimate,NA_real_,estimate),
             conf.low = if_else(conf.low>5 | -5 > conf.low,NA_real_,conf.low),
             conf.high = if_else(conf.high>5 | -5 > conf.high,NA_real_,conf.high))
    data <- stats::model.frame(model)
    
    forest_terms <- tibble::tibble(
      term_label = attr(model$terms, "term.labels"),
      variable = remove_backticks(term_label)
    ) %>%
      inner_join(
        tibble::tibble(
          variable = names(attr(model$terms, "dataClasses"))[-1],
          class = attr(model$terms, "dataClasses")[-1]
        ),
        by = "variable"
      )
    
    forest_labels <- tibble::tibble(
      variable = names(data),
      label = vapply(
        data,
        function(x) attr(x, "label", exact = TRUE) %||% NA_character_,
        character(1)
      ) %>%
        coalesce(variable)
    )
    
    create_term_data <- function(term_row) {
      if (!is.na(term_row$class)) {
        var <- term_row$variable
        if (term_row$class %in% c("factor", "character")) {
          tab <- table(data[, var])
          if (!any(paste0(term_row$term_label, names(tab)) %in% tidy_model$term)) {
            # Filter out terms not in final model summary (e.g. strata)
            out <- tibble::tibble(variable = NA)
          } else {
            out <- data.frame(
              term_row,
              level = names(tab),
              level_no = 1:length(tab),
              n = as.integer(tab),
              stringsAsFactors = FALSE
            )
            if (factor_separate_line) {
              out <- bind_rows(tibble::as_tibble(term_row), out)
            }
            if (inherits(model, "coxph")) {
              data_event <- bind_cols(data[, -1, drop = FALSE],
                                      .event_time = data[, 1][, "time"],
                                      .event_status = data[, 1][, "status"]
              )
              event_detail_tab <- data_event %>%
                group_by(!!as.name(var)) %>%
                summarise(
                  person_time = sum(.event_time),
                  n_events = sum(.event_status)
                )
              colnames(event_detail_tab)[1] <- "level"
              event_detail_tab$level <- as.character(event_detail_tab$level)
              out <- out %>% left_join(event_detail_tab, by = "level")
            }
          }
        } else {
          out <- data.frame(term_row,
                            level = NA, level_no = NA, n = sum(!is.na(data[, var])),
                            stringsAsFactors = FALSE
          )
          if (term_row$class == "logical") {
            out$term_label <- paste0(term_row$term_label, "TRUE")
          }
        }
      } else {
        out <- data.frame(term_row, level = NA, level_no = NA, n = NA, stringsAsFactors = FALSE)
      }
      out
    }
    forest_terms <- forest_terms %>%
      rowwise() %>%
      do(create_term_data(.)) %>%
      ungroup() %>%
      filter(!is.na(variable)) %>%
      mutate(term = paste0(term_label, replace(level, is.na(level), ""))) %>%
      left_join(tidy_model, by = "term") %>%
      mutate(
        reference = ifelse(is.na(level_no), FALSE, level_no == 1),
        estimate = ifelse(reference, 0, estimate),
        variable = ifelse(is.na(variable), remove_backticks(term), variable)
      ) %>%
      mutate(
        variable = ifelse(is.na(level_no) | (level_no == 1 & !factor_separate_line), variable, NA)
      ) %>%
      left_join(
        forest_labels,
        by = "variable"
      ) %>%
      mutate(
        variable = coalesce(label, variable)
      )
    if (!is.null(covariates)) {
      forest_terms <- filter(forest_terms, term_label %in% covariates)
    }
    
    forest_terms
  }
  
  if (!is.null(model_list)) {
    forest_terms <- lapply(seq_along(model_list), function(i) {
      make_forest_terms(model_list[[i]]) %>%
        mutate(model_name = model_names[i])
    }) %>%
      bind_rows()
    if (merge_models) {
      forest_terms$model_name <- NULL
    }
  } else {
    forest_terms <- make_forest_terms(model)
  }
  
  # #use_exp <- grepl("exp", deparse(trans))
  if (!is.null(limits)) {
    forest_terms <- forest_terms %>%
      mutate(
        arrow_tag.l = limits[1],
        arrow_tag.r = limits[2],
        arrow_tag.l = ifelse(conf.low < .data$arrow_tag.l, TRUE, FALSE),
        arrow_tag.r = ifelse(conf.high > .data$arrow_tag.r, TRUE, FALSE)
      ) %>%
      mutate(
        plot_range.low = ifelse(.data$arrow_tag.l, limits[1], conf.low),
        plot_range.high = ifelse(.data$arrow_tag.r, limits[2], conf.high)
      )
  }
  
  
  
  plot_data <- list(
    forest_data = forest_terms,
    mapping = mapping,
    panels = panels, trans = trans,
    funcs = funcs, format_options = format_options, theme = theme,
    limits = limits, breaks = breaks, recalculate_width = recalculate_width,
    recalculate_height = recalculate_height, exclude_infinite_cis = exclude_infinite_cis
  )
  main_plot <- do.call("panel_forest_plot", plot_data)
  if (return_data) {
    list(plot_data = plot_data, plot = main_plot)
  } else {
    main_plot
  }
}