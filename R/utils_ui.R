# ---------------------------------------------------------------------------
# Small UI helpers shared across modules
# ---------------------------------------------------------------------------

#' A label followed by a hover tooltip help icon
#' @noRd
ff_tooltip <- function(label, text) {
  shiny::tagList(
    label,
    bslib::tooltip(
      bsicons::bs_icon("info-circle-fill", class = "ff-help"),
      text,
      placement = "right"
    )
  )
}

#' A compact metric tile
#' @noRd
ff_metric <- function(label, value) {
  shiny::div(
    class = "ff-metric",
    shiny::div(class = "ff-metric-label", label),
    shiny::div(class = "ff-metric-value", value)
  )
}

#' The wizard step sidebar, highlighting the current step
#' @noRd
ff_stepper <- function(current) {
  steps <- c("Upload data", "Review & QC", "Filter", "Traits (optional)",
             "Build & export")
  items <- lapply(seq_along(steps), function(i) {
    cls <- if (i == current) "ff-step active" else if (i < current) "ff-step done" else "ff-step"
    icon <- if (i < current) bsicons::bs_icon("check-lg") else as.character(i)
    shiny::tags$li(
      class = cls,
      shiny::span(class = "ff-step-num", icon),
      shiny::span(steps[i])
    )
  })
  shiny::tags$ul(class = "ff-stepper", items)
}

#' A standard Back / Next button row for a wizard step
#' @noRd
ff_nav_buttons <- function(ns, back = TRUE, next_label = "Continue",
                           next_icon = "arrow-right", back_id = "back",
                           next_id = "next") {
  shiny::div(
    class = "d-flex justify-content-between mt-4",
    if (back) {
      shiny::actionButton(ns(back_id), label = shiny::tagList(
        bsicons::bs_icon("arrow-left"), " Back"
      ), class = "btn-outline-secondary")
    } else {
      shiny::span()
    },
    shiny::actionButton(ns(next_id), label = shiny::tagList(
      next_label, " ", bsicons::bs_icon(next_icon)
    ), class = "btn-primary")
  )
}
