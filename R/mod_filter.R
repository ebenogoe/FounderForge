# Wizard step 3: filtering --------------------------------------------------

#' @noRd
mod_filter_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header(shiny::h3("3. Filter markers")),
    bslib::card_body(
      bslib::layout_columns(
        col_widths = c(6, 6),
        shiny::sliderInput(ns("maf"),
          ff_tooltip("Minimum MAF",
            "Markers with a minor-allele frequency below this are removed. 0.05 is a common starting point."),
          min = 0, max = 0.5, value = 0.05, step = 0.01),
        shiny::sliderInput(ns("max_miss_snp"),
          ff_tooltip("Max SNP missingness",
            "Drop markers missing in more than this fraction of individuals. 1 keeps all."),
          min = 0, max = 1, value = 1, step = 0.05)
      ),
      bslib::layout_columns(
        col_widths = c(6, 6),
        shiny::sliderInput(ns("max_miss_ind"),
          ff_tooltip("Max individual missingness",
            "Drop individuals missing more than this fraction of markers."),
          min = 0, max = 1, value = 1, step = 0.05),
        shiny::numericInput(ns("min_snp_per_chr"),
          ff_tooltip("Min SNPs per chromosome",
            "Chromosomes with fewer markers than this are dropped. AlphaSimR needs at least 2."),
          value = 2, min = 2, step = 1)
      ),
      shiny::checkboxInput(ns("drop_mono"), "Remove monomorphic markers", TRUE),
      shiny::textInput(ns("drop_ind"),
        ff_tooltip("Individuals to drop (optional)",
          "Comma-separated list of sample ids to exclude (e.g. flagged samples)."), ""),
      shiny::actionButton(ns("preview"), shiny::tagList(
        bsicons::bs_icon("funnel"), " Preview retained markers"),
        class = "btn-outline-secondary"),
      shiny::uiOutput(ns("preview_out")),
      ff_nav_buttons(ns, back = TRUE, next_label = "Apply & continue")
    )
  )
}

#' @noRd
mod_filter_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {

    do_filter <- function() {
      drops <- trimws(strsplit(input$drop_ind %||% "", ",")[[1]])
      drops <- drops[nzchar(drops)]
      filter_genotypes(
        rv$geno, rv$snp_map,
        maf = input$maf,
        max_miss_snp = input$max_miss_snp,
        max_miss_ind = input$max_miss_ind,
        drop_mono = input$drop_mono,
        min_snp_per_chr = input$min_snp_per_chr,
        drop_individuals = drops,
        logger = rv$logger
      )
    }

    shiny::observeEvent(input$preview, {
      shinybusy::show_modal_spinner(spin = "fading-circle", color = "#1D9E75",
        text = "Previewing retained markers...")
      on.exit(shinybusy::remove_modal_spinner(), add = TRUE)
      res <- tryCatch(do_filter(), error = function(e) e)
      output$preview_out <- shiny::renderUI({
        if (inherits(res, "error")) {
          bslib::card(class = "border-danger mt-3", bslib::card_body(
            shiny::span(class = "text-danger",
              bsicons::bs_icon("exclamation-triangle-fill"), " ", conditionMessage(res))))
        } else {
          warn <- if (nrow(res$geno) < 5000)
            shiny::p(class = "text-warning mb-0",
              bsicons::bs_icon("exclamation-triangle"),
              " Fewer than 5,000 markers retained - consider relaxing the MAF threshold.")
          bslib::card(class = "border-info mt-3", bslib::card_body(
            shiny::strong(sprintf("%s SNPs x %s individuals retained.",
              format(nrow(res$geno), big.mark = ","),
              format(ncol(res$geno), big.mark = ","))),
            warn))
        }
      })
    })

    shiny::observeEvent(input$back, rv$step <- 2L)
    shiny::observeEvent(input$`next`, {
      shinybusy::show_modal_spinner(spin = "fading-circle", color = "#1D9E75",
        text = "Applying filters...")
      on.exit(shinybusy::remove_modal_spinner(), add = TRUE)
      res <- tryCatch(do_filter(), error = function(e) e)
      if (inherits(res, "error")) {
        output$preview_out <- shiny::renderUI(
          bslib::card(class = "border-danger mt-3", bslib::card_body(
            shiny::span(class = "text-danger",
              bsicons::bs_icon("exclamation-triangle-fill"), " ", conditionMessage(res)))))
        return(invisible())
      }
      rv$filtered <- res
      rv$step <- 4L
    })
  })
}
