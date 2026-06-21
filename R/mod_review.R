# Wizard step 2: review & QC ------------------------------------------------

#' @noRd
mod_review_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header(shiny::h3("2. Review & quality control")),
    bslib::card_body(
      shiny::p(class = "text-muted", "A summary of your uploaded data before filtering."),
      shiny::uiOutput(ns("metrics")),
      shiny::uiOutput(ns("banner")),
      bslib::layout_columns(
        col_widths = c(6, 6),
        bslib::card(bslib::card_header("Minor allele frequency"),
                    shiny::plotOutput(ns("maf"), height = "240px")),
        bslib::card(bslib::card_header("Per-SNP missingness"),
                    shiny::plotOutput(ns("miss"), height = "240px"))
      ),
      bslib::card(bslib::card_header("SNPs per chromosome"),
                  shiny::plotOutput(ns("spc"), height = "240px")),
      ff_nav_buttons(ns, back = TRUE, next_label = "Continue to filtering")
    )
  )
}

#' @noRd
mod_review_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(rv$step, {
      if (identical(rv$step, 2L) && !is.null(rv$geno) && is.null(rv$qc)) {
        shinybusy::show_modal_spinner(spin = "fading-circle", color = "#1D9E75",
          text = "Computing QC summary...")
        on.exit(shinybusy::remove_modal_spinner(), add = TRUE)
        rv$qc <- compute_qc_summary(rv$geno, rv$snp_map, rv$logger)
      }
    })

    output$metrics <- shiny::renderUI({
      shiny::req(rv$qc)
      q <- rv$qc
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3), class = "mb-3",
        ff_metric("SNPs", format(q$n_snp, big.mark = ",")),
        ff_metric("Individuals", format(q$n_ind, big.mark = ",")),
        ff_metric("Chromosomes", q$n_chr),
        ff_metric("Missing", sprintf("%.2f%%", q$pct_missing))
      )
    })

    output$banner <- shiny::renderUI({
      shiny::req(rv$qc)
      q <- rv$qc
      msg <- sprintf("Mean MAF %.3f \u00b7 %s monomorphic marker(s) \u00b7 mean heterozygosity %.3f",
                     q$mean_maf, format(q$n_monomorphic, big.mark = ","), q$mean_het)
      bslib::card(class = "border-success mb-3",
        bslib::card_body(shiny::span(class = "text-success",
          bsicons::bs_icon("check-circle-fill"), " ", msg)))
    })

    output$maf  <- shiny::renderPlot({ shiny::req(rv$qc); hist_plot(rv$qc$maf, "", "MAF") })
    output$miss <- shiny::renderPlot({ shiny::req(rv$qc); hist_plot(rv$qc$miss_snp, "", "Missing rate") })
    output$spc  <- shiny::renderPlot({ shiny::req(rv$qc); snps_per_chr_plot(rv$qc$snps_per_chr) })

    shiny::observeEvent(input$back, rv$step <- 1L)
    shiny::observeEvent(input$`next`, rv$step <- 3L)
  })
}
