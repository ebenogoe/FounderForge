# Wizard step 4: optional trait definition ----------------------------------

#' @noRd
mod_traits_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header(shiny::h3("4. Define traits (optional)")),
    bslib::card_body(
      shiny::p(class = "text-muted",
        "Optionally add additive traits to build an AlphaSimR SimParam. Leave at zero to export the founder population only."),
      shiny::numericInput(ns("n_traits"),
        ff_tooltip("Number of traits",
          "Each additive trait needs a mean, variance, QTL count and heritability."),
        value = 0, min = 0, max = 6, step = 1),
      shiny::uiOutput(ns("trait_rows")),
      ff_nav_buttons(ns, back = TRUE, next_label = "Continue to export")
    )
  )
}

#' @noRd
mod_traits_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$trait_rows <- shiny::renderUI({
      n <- input$n_traits %||% 0
      if (n < 1) return(NULL)
      rows <- lapply(seq_len(n), function(i) {
        bslib::card(class = "mb-2", bslib::card_body(
          shiny::strong(sprintf("Trait %d", i)),
          bslib::layout_columns(
            col_widths = c(4, 4, 4),
            shiny::textInput(ns(paste0("name_", i)), "Name", sprintf("Trait%d", i)),
            shiny::numericInput(ns(paste0("nqtl_", i)),
              ff_tooltip("QTL per chr", "Number of QTL placed on each chromosome."),
              value = 10, min = 1, step = 1),
            shiny::numericInput(ns(paste0("h2_", i)),
              ff_tooltip("Heritability", "Narrow/broad-sense h2, between 0 and 1."),
              value = 0.5, min = 0, max = 1, step = 0.05)
          ),
          bslib::layout_columns(
            col_widths = c(6, 6),
            shiny::numericInput(ns(paste0("mean_", i)), "Mean", value = 100),
            shiny::numericInput(ns(paste0("var_", i)), "Variance", value = 25, min = 0)
          )
        ))
      })
      shiny::tagList(rows)
    })

    collect_specs <- function() {
      n <- input$n_traits %||% 0
      if (n < 1) return(list())
      lapply(seq_len(n), function(i) {
        list(
          name = input[[paste0("name_", i)]] %||% sprintf("Trait%d", i),
          nQtlPerChr = input[[paste0("nqtl_", i)]] %||% 10,
          mean = input[[paste0("mean_", i)]] %||% 0,
          var  = input[[paste0("var_", i)]]  %||% 1,
          h2   = input[[paste0("h2_", i)]]   %||% 0.5
        )
      })
    }

    shiny::observeEvent(input$back, rv$step <- 3L)
    shiny::observeEvent(input$`next`, {
      rv$trait_specs <- collect_specs()
      rv$step <- 5L
    })
  })
}
