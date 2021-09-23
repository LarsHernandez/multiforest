forest_model_lars <- function(model,
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
    tidy_model <- broom::tidy(model, conf.int = FALSE)
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











mforestmodel <- function(data, pala="grey15", palb="grey75", lim=c(-2.6,4), dependent="delta", 
                         legend_position="none", header = NULL, spaces = NULL) {

  tgt_uni <- paste0(dependent,"~", names(data)[-grep(dependent, names(data))])
  tgt_mul <- paste0(dependent,"~", paste(names(data)[-grep(dependent, names(data))], collapse="+"))
  
  fit_uni <- Map(function(x) glm(as.formula(x), family = "poisson", data = data), tgt_uni)
  fit_mul <- glm(tgt_mul, family = "poisson", data = data)
  
  
trans <- exp
  

fd2 <- forest_model_lars(model = fit_mul,      return_data = T, factor_separate_line = T, exponentiate = T)$plot_data
fd1 <- forest_model_lars(model_list = fit_uni, return_data = T, factor_separate_line = T, exponentiate = T)$plot_data

fd1$panels[[11]] <- fd2$panels[[7]]
fd1$panels[[12]] <- fd2$panels[[8]]
fd1$panels[[13]] <- fd2$panels[[9]]
fd1$panels[[14]] <- fd2$panels[[10]]

if (is.null(header)) {
  head <- c("Variable","N","Relative Risk","Univariate","P","Multivariate","P")
} else {head <- header}

if (is.null(spaces)) {
  space <- c(0.015,0.27,0.2,0.005,0.2,0.02)
} else {space <- spaces}


fd1$panels[[2]]$heading <- head[1]
fd1$panels[[4]]$heading <- head[2]
fd1$panels[[6]]$heading <- head[3]
fd1$panels[[8]]$heading <- head[4]
fd1$panels[[9]]$heading <- head[5]
fd1$panels[[12]]$heading <-head[6]
fd1$panels[[13]]$heading <-head[7]

fd1$panels[[1]]$width <- space[1]
fd1$panels[[3]]$width <- space[2]
fd1$panels[[8]]$width <- space[3]
fd1$panels[[10]]$width <-space[4]
fd1$panels[[12]]$width <-space[5]
fd1$panels[[14]]$width <-space[6]

#fd1$panels[[13]]$color <- "Multivariate"
#fd1$panels[[9]]$color <- "Univariate"

fd1$panels[[12]]$display <- new_quosure(
  quote(if_else(reference, "Reference", sprintf("%0.2f (%0.2f, %0.2f)",
                                                trans(estimate_2), trans(conf.low_2), trans(conf.high_2)))))

fd1$panels[[13]]$display<- new_quosure(
  quote(if_else(reference, "", format.pval(p.value_2, digits = 1,eps = 0.001))))


# Function ----------------------------------------------------------------

forest_data <- fd1$forest_data
forest_data_2 <- fd2$forest_data
mapping <- aes(estimate, xmin = conf.low, xmax = conf.high)
panels <- fd1$panels

funcs <- NULL
format_options <- list(colour = "black", shape = 15, banded = F, text_size = 5, point_size = 5)
theme <- theme_forest()
limits <- lim
breaks <- NULL
recalculate_width <- FALSE
recalculate_height <- TRUE
exclude_infinite_cis <- TRUE


# Update older style panels to quosures
for (i in seq_along(panels)) {
  if (!is.null(panels[[i]]$display)) {
    panels[[i]]$display <- as_quosure(panels[[i]]$display)
  }
  if (!is.null(panels[[i]]$display_na)) {
    panels[[i]]$display_na <- as_quosure(panels[[i]]$display_na)
  } else {
    panels[[i]]$display_na <- panels[[i]]$display
  }
}

format_options$colour <- format_options$colour %||% "black"
format_options$shape <- format_options$shape %||% 15
format_options$banded <- format_options$banded %||% TRUE
format_options$text_size <- format_options$text_size %||% 5

## Modify the square size
## mapping$size <- mapping$size %||% 5
mapping$whole_row <- mapping$whole_row %||% FALSE

if (!is.null(mapping$section) && !all(is.na(eval_tidy(mapping$section, forest_data)))) {
  forest_data <- forest_data %>%
    mutate(
      .section = !!mapping$section,
      .whole_row = !!mapping$whole_row
    ) %>%
    group_by(.section) %>%
    do({
      bind_rows(
        tibble(
          .section = first(.$.section),
          .subheading = first(.$.section),
          .whole_row = TRUE
        ),
        as_tibble(.)
      )
    })
  mapping$subheading <- quo(.subheading)
  mapping$whole_row <- quo(.whole_row)
}

default_cols <- function(.data, ...) {
  dots <- list2(...)
  if (is.null(names(dots)) || any(names(dots) == "")) {
    stop("All parameters (except `.data`) must be named")
  }
  for (cn in names(dots)) {
    if (!has_name(.data, cn)) {
      .data[[cn]] <- eval_tidy(dots[[cn]], .data)
    }
  }
  .data
}



fd_for_eval <- c(as.list(forest_data), trans = trans, funcs)

mapped_data <- lapply(mapping, function(el) {
  eval_tidy(el, fd_for_eval)
}) %>%
  as_tibble() %>%
  default_cols(
    band = TRUE,
    diamond = FALSE,
    section = 1
  ) %>%
  mutate(
    diamond = if_else(is.na(diamond), FALSE, diamond),
    whole_row = if_else(is.na(whole_row), FALSE, whole_row),
    y = n():1
  )

fd_for_eval_2 <- c(as.list(forest_data_2), trans = trans, funcs)
mapped_data_2 <- lapply(mapping, function(el) {
  eval_tidy(el, fd_for_eval_2)
}) %>%
  as_tibble() %>%
  default_cols(
    band = TRUE,
    diamond = FALSE,
    section = 1
  ) %>%
  mutate(
    diamond = if_else(is.na(diamond), FALSE, diamond),
    whole_row = if_else(is.na(whole_row), FALSE, whole_row),
    y = n():1
  )

fd_for_eval$estimate_2 <- fd_for_eval_2$estimate
fd_for_eval$conf.low_2 <- fd_for_eval_2$conf.low
fd_for_eval$conf.high_2 <- fd_for_eval_2$conf.high
fd_for_eval$p.value_2 <- fd_for_eval_2$p.value

mapped_text <- lapply(seq(panels), function(i) {
  if (!is.null(panels[[i]]$display)) {
    as.character(ifelse(!is.na(mapped_data$x),
                        eval_tidy(panels[[i]]$display, fd_for_eval),
                        eval_tidy(panels[[i]]$display_na, fd_for_eval)
    ))
  } else {
    NULL
  }
})

max_y <- max(mapped_data$y)

panel_positions <- lapply(panels, function(panel) {
  tibble(
    width = panel$width %||% NA,
    item = panel$item %||% {
      if (!is.null(panel$display)) "text" else NA
    },
    display = quo_text(panel$display),
    hjust = as.numeric(panel$hjust %||% 0),
    heading = panel$heading %||% NA,
    fontface = panel$fontface %||% "plain",
    linetype = panel$linetype %||% {
      if (!is.na(item) && item == "vline") "solid" else NA
    },
    line_x = as.numeric(panel$line_x %||% NA),
    parse = as.logical(panel$parse %||% FALSE),
    width_group = panel$width_group %||% NA
  )
}) %>%
  bind_rows()

if (any(panel_positions$parse & panel_positions$fontface != "plain")) {
  warning("Fontface cannot be applied to parsed text; please use the plotmath functions (e.g. bold())")
}

if (any(is.na(panel_positions$width)) && !recalculate_width) {
  recalculate_width <- TRUE
  message("Some widths are undefined; defaulting to recalculate_width = TRUE")
}

if (sum(panel_positions$item == "forest", na.rm = TRUE) != 1) {
  stop("One panel must include item \"forest\".")
}

forest_panel <- which(panel_positions$item == "forest")

if (is.null(limits)) {
  x_min_max <- c(mapped_data$xmin, mapped_data$xmax)
  x_min_max <- x_min_max[is.finite(x_min_max)]
  forest_min_max <- range(x_min_max)
  if (!is.na(forest_line_x <- panel_positions$line_x[forest_panel])) {
    if (forest_min_max[2] < forest_line_x) {
      forest_min_max[2] <- forest_line_x
      message("Resized limits to included dashed line in forest panel")
    }
    if (forest_min_max[1] > forest_line_x) {
      forest_min_max[1] <- forest_line_x
      message("Resized limits to included dashed line in forest panel")
    }
  }
} else {
  forest_min_max <- limits
}

if (!is.null(recalculate_height) && !(identical(recalculate_height, FALSE))) {
  if (identical(recalculate_height, TRUE)) {
    recalculate_height <- graphics::par("din")[2]
  }
  max_text_size <- recalculate_height / (max_y + 1) / 1.3 * 25.4
  if (format_options$text_size > max_text_size) {
    format_options$text_size <- max_text_size
  }
}

if (!is.null(recalculate_width) && !(identical(recalculate_width, FALSE))) {
  panel_positions <-
    forestmodel:::recalculate_width_panels(panel_positions,
                             mapped_text = mapped_text,
                             mapped_data = mapped_data,
                             recalculate_width = recalculate_width,
                             format_options = format_options,
                             theme = theme
    )
}

panel_positions <- panel_positions %>%
  mutate(
    rel_width = width / width[forest_panel],
    rel_x = cumsum(c(0, width[-n()])),
    rel_x = (rel_x - rel_x[forest_panel]) / width[forest_panel],
    abs_x = rel_x * diff(forest_min_max) + forest_min_max[1],
    abs_width = rel_width * diff(forest_min_max),
    abs_end_x = abs_x + abs_width,
    text_x = ifelse(hjust == 0, abs_x,
                    ifelse(hjust == 0.5, abs_x + abs_width / 2, abs_end_x)
    )
  )

forest_vlines <- panel_positions %>%
  filter(item == "vline" | !is.na(linetype)) %>%
  transmute(
    x = case_when(
      !is.na(line_x) ~ line_x,
      hjust == 1 ~ abs_end_x,
      hjust == 0.5 ~ abs_x + abs_width / 2
    ),
    xend = x,
    y = 0.5,
    yend = if_else(is.na(line_x), max_y + 1.5, max_y + 0.5),
    linetype = linetype
  )

forest_hlines <- tibble(
  x = c(min(panel_positions$abs_x), max(panel_positions$abs_end_x)),
  y = max_y + 0.5,
  linetype = "solid"
)

forest_headings <- panel_positions %>%
  filter(!is.na(heading)) %>%
  transmute(
    x = text_x,
    y = max_y + 1,
    hjust = hjust,
    label = heading,
    fontface = "bold",
    parse = FALSE
  )

if (has_name(mapped_data, "subheading")) {
  forest_subheadings <- mapped_data %>%
    filter(!is.na(subheading)) %>%
    transmute(
      x = mean(forest_min_max),
      y,
      hjust = 0.5,
      label = subheading,
      fontface = "bold",
      parse = FALSE
    )
  forest_headings <- bind_rows(forest_headings, forest_subheadings)
}

if (any(mapped_data$whole_row)) {
  forest_whole_row_back <- mapped_data %>%
    filter(whole_row) %>%
    transmute(
      y,
      xmin = min(panel_positions$abs_x),
      xmax = max(panel_positions$abs_end_x),
      ymin = y - 0.5,
      ymax = y + 0.5
    )
  forest_hlines <- bind_rows(
    forest_hlines,
    forest_whole_row_back %>%
      rowwise() %>%
      do({
        tibble(
          x = rep(range(c(panel_positions$abs_x, panel_positions$abs_end_x)), 2),
          y = rep(.$y + c(-0.5, 0.5), each = 2),
          linetype = "solid"
        )
      })
  )
}

forest_text <- lapply(seq(panels), function(i) {
  if (!is.null(mapped_text[[i]])) {
    with(
      panel_positions[i, ],
      tibble(
        x = text_x,
        y = mapped_data$y,
        hjust = hjust,
        label = mapped_text[[i]],
        fontface = fontface,
        parse = parse
      )
    )
  }
}) %>%
  bind_rows() %>%
  bind_rows(., forest_headings)

if (format_options$banded) {
  forest_rectangles <- mapped_data %>%
    filter(band) %>%
    group_by(section) %>%
    filter(row_number() %% 2 == 1) %>%
    transmute(
      xmin = min(panel_positions$abs_x),
      xmax = max(panel_positions$abs_end_x),
      y = y,
      ymin = y - 0.5,
      ymax = y + 0.5
    )
}

if (any(mapped_data$diamond)) {
  forest_diamonds <- mapped_data %>%
    filter(diamond == TRUE) %>%
    rowwise() %>%
    do({
      tibble(
        x = c(.$xmin, .$x, .$xmax, .$x, .$xmin),
        y = .$y + c(0, 0.15, 0, -0.15, 0)
      )
    }) %>%
    ungroup() %>%
    mutate(group = (row_number() + 4) %/% 5)
  
  forest_hlines <- bind_rows(
    forest_hlines,
    mapped_data %>%
      filter(diamond) %>%
      rowwise() %>%
      do({
        tibble(
          x = rep(range(c(panel_positions$abs_x, panel_positions$abs_end_x)), 2),
          y = rep(.$y + c(-0.5, 0.5), each = 2),
          linetype = "solid"
        )
      })
  )
}

forest_hlines <- mutate(forest_hlines, group = (row_number() + 1) %/% 2)

if (is.null(breaks)) {
  if (identical(trans, exp)) {
    breaks <- log(grDevices::axisTicks(log10(exp(forest_min_max)), TRUE))
  } else {
    breaks <- grDevices::axisTicks(forest_min_max, FALSE)
  }
  breaks <- breaks[breaks >= forest_min_max[1] & breaks <= forest_min_max[2]]
}

if (exclude_infinite_cis) {
  mapped_data <- mapped_data %>%
    mutate_at(
      vars(x, xmin, xmax),
      list(~if_else(
        is.infinite(x) | is.infinite(xmin) | is.infinite(xmax),
        NA_real_,
        .
      ))
    )
}




mapped_data <- rbind(mapped_data %>% mutate(type="Univariate"),
                     mapped_data_2 %>% mutate(type="Multivariate")) %>% 
  mutate(type = fct_relevel(type,"Multivariate")) %>% 
  mutate(xmin = if_else(xmin< -2.4,-2.5,xmin))

main_plot <- ggplot()

if (format_options$banded) {
  main_plot <- main_plot +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              forest_rectangles, na.rm = TRUE,
              color = "grey90", fill="white", size=0.3)
}

main_plot <- main_plot + 
  geom_line(
  aes(x, y, linetype = linetype, group = group),
  forest_hlines, na.rm = TRUE) 

main_plot <- main_plot +
  geom_segment(
    aes(x, y, xend = xend, yend = yend, linetype = linetype),
    forest_vlines, na.rm = TRUE)

if (any(mapped_data$diamond)) {
  main_plot <- main_plot +
    geom_polygon(aes(x, y, group = group), forest_diamonds, fill = format_options$colour, na.rm = TRUE)
}

main_plot <- main_plot +
  geom_errorbar(aes(y = y, x = x, xmin = xmin, xmax = xmax, color=(type)), position=position_dodge(0.6),
                width=0.4,size=0.4, na.rm = TRUE, show.legend = F,
                filter(mapped_data, !diamond & !(is.na(xmin) & is.na(xmax))))


main_plot <- main_plot +
  geom_pointrange(aes(y=y, x=x, xmin = xmin, xmax = xmax, color=(type)), position=position_dodge(0.6),
                  size = 0.5,fatten = 6, filter(mapped_data, !diamond) %>% mutate(xmin=if_else(is.na(xmin),x,xmin),
                                                                                  xmax=if_else(is.na(xmax),x,xmax)),
                  shape = format_options$shape, na.rm = TRUE)



if (any(mapped_data$whole_row)) {
  main_plot <- main_plot +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              forest_whole_row_back,
              fill = "#FFFFFF"
    )
}

for (parse_type in unique(forest_text$parse)) {
  main_plot <- main_plot +
    geom_text(aes(x, y, label = label, hjust = hjust, fontface = fontface),
              filter(forest_text, parse == parse_type),
              na.rm = TRUE, parse = parse_type,
              size = 3.5)
}






main_plot <- main_plot +
  scale_color_manual(values=c(pala, palb))+
  labs(color=NULL) +
  scale_linetype_identity() +
  scale_alpha_identity() +
  scale_x_continuous(breaks = breaks, labels = sprintf("%g", trans(breaks)), expand = c(0,0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme + 
  theme(legend.position = legend_position, 
        legend.background = element_rect(fill=NA),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.text =element_text(color="black"), 
        axis.ticks = element_line(color="black"))

main_plot
}

