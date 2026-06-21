# Wizard step 5: build & export ---------------------------------------------

#' @noRd
mod_build_export_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header(shiny::h3("5. Build & export")),
    bslib::card_body(
      shiny::p(class = "text-muted",
        "Choose where to save the results, then build your AlphaSimR founder population."),

      shinyFiles::shinyDirButton(
        ns("outdir"), "Choose output folder", "Select an output folder",
        class = "btn-outline-primary"),
      shiny::div(class = "mt-2 text-muted", shiny::textOutput(ns("outdir_path"))),

      shiny::hr(),
      shinyWidgets::materialSwitch(
        ns("do_pca"),
        ff_tooltip("Run before/after PCA validation",
          "Compares PCA of your data before and after the AlphaSimR import. If the two agree, the import is faithful. Adds time on large datasets."),
        value = TRUE, status = "success"),

      shiny::actionButton(ns("run"), shiny::tagList(
        bsicons::bs_icon("play-fill"), " Build founder population"),
        class = "btn-primary btn-lg mt-2"),

      shiny::uiOutput(ns("result")),
      shiny::uiOutput(ns("pca_block")),

      shiny::tags$details(class = "mt-3",
        shiny::tags$summary("Run log"),
        shiny::verbatimTextOutput(ns("log"))),

      shiny::div(class = "d-flex justify-content-between mt-4",
        shiny::actionButton(ns("back"), shiny::tagList(
          bsicons::bs_icon("arrow-left"), " Back"), class = "btn-outline-secondary"),
        shiny::actionButton(ns("restart"), shiny::tagList(
          bsicons::bs_icon("arrow-repeat"), " Start a new run"),
          class = "btn-outline-secondary"))
    )
  )
}

#' @noRd
mod_build_export_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    roots <- c(Home = path.expand("~"), Root = "/")
    shinyFiles::shinyDirChoose(input, "outdir", roots = roots, session = session)
    outdir <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$outdir, {
      p <- shinyFiles::parseDirPath(roots, input$outdir)
      if (length(p) > 0) outdir(as.character(p))
    })
    output$outdir_path <- shiny::renderText({
      if (is.null(outdir())) "No folder selected yet." else paste("Saving to:", outdir())
    })

    shiny::observeEvent(input$run, {
      if (is.null(outdir())) {
        shiny::showNotification("Please choose an output folder first.", type = "warning")
        return(invisible())
      }
      shiny::req(rv$filtered)
      logger <- rv$logger

      shinybusy::show_modal_spinner(spin = "fading-circle", color = "#1D9E75",
        text = "Building founder population...")
      on.exit(shinybusy::remove_modal_spinner(), add = TRUE)
      step <- function(msg) shinybusy::update_modal_spinner(msg)

      res <- tryCatch({
        step("Imputing missing genotypes...")
        imp <- impute_missing(rv$filtered$geno, logger)

        step("Constructing haplotypes & genetic map...")
        br <- build_founder_pop(imp$geno, rv$filtered$snp_map,
                                ploidy = rv$ploidy, strategy = rv$strategy,
                                cm_per_mb = rv$cm_per_mb, logger = logger)

        SP <- NULL
        if (length(rv$trait_specs %||% list()) > 0) {
          step("Defining traits...")
          SP <- define_traits(br$founderPop, rv$trait_specs, logger)
        }

        step("Verifying heterozygosity...")
        ver <- verify_founder(br$founderPop, br$strategy, logger = logger)

        pca <- NULL
        if (isTRUE(input$do_pca)) {
          step("Running before/after PCA...")
          pca <- tryCatch(pca_validation(imp$geno, br$founderPop, logger = logger),
                          error = function(e) { logger$log(conditionMessage(e), "WARN"); NULL })
        }

        step("Writing outputs...")
        run_dir <- write_outputs(outdir(), br, SP = SP, qc_summary = rv$qc,
                                 filter_report = rv$filtered$report,
                                 verify_result = ver, pca_result = pca,
                                 source_info = rv$source_info, logger = logger)
        list(br = br, SP = SP, ver = ver, pca = pca, run_dir = run_dir)
      }, error = function(e) e)

      if (inherits(res, "error")) {
        output$result <- shiny::renderUI(
          bslib::card(class = "border-danger mt-3", bslib::card_body(
            shiny::span(class = "text-danger",
              bsicons::bs_icon("exclamation-triangle-fill"),
              " Build failed: ", conditionMessage(res)))))
        output$log <- shiny::renderText(paste(logger$lines(), collapse = "\n"))
        return(invisible())
      }

      rv$run_dir <- res$run_dir
      rv$pca <- res$pca

      output$result <- shiny::renderUI({
        v <- res$ver
        het_cls <- if (v$pass) "text-success" else "text-warning"
        het_ic  <- if (v$pass) "check-circle-fill" else "exclamation-triangle-fill"
        bslib::card(class = "border-success mt-3", bslib::card_body(
          shiny::h4(class = "text-success",
            bsicons::bs_icon("check-circle-fill"), " Founder population built"),
          shiny::p(sprintf("%d individuals \u00b7 %d chromosome(s) \u00b7 %s segregating loci",
            res$br$founderPop@nInd, res$br$founderPop@nChr,
            format(v$n_seg, big.mark = ","))),
          shiny::p(class = het_cls, bsicons::bs_icon(het_ic),
            sprintf(" Heterozygosity check: %s (observed %s; expected %s)",
              if (v$pass) "passed" else "review",
              paste(v$values, collapse = ", "), v$expected)),
          shiny::p(shiny::tags$strong("Saved to: "), res$run_dir)
        ))
      })

      output$pca_block <- shiny::renderUI({
        if (is.null(res$pca)) return(NULL)
        bslib::card(class = "mt-3",
          bslib::card_header(shiny::tagList(
            "Before/after PCA validation ",
            shiny::downloadButton(session$ns("dl_pca"), "Download plot",
                                  class = "btn-sm btn-outline-secondary float-end"))),
          shiny::plotOutput(session$ns("pca_plot"), height = "360px"))
      })

      output$pca_plot <- shiny::renderPlot({ shiny::req(rv$pca); rv$pca$plot })
      output$log <- shiny::renderText(paste(logger$lines(), collapse = "\n"))
    })

    output$dl_pca <- shiny::downloadHandler(
      filename = function() "pca_before_after.png",
      content = function(file) {
        ggplot2::ggsave(file, plot = rv$pca$plot, width = 11, height = 5, dpi = 150)
      }
    )

    shiny::observeEvent(input$back, rv$step <- 4L)
    shiny::observeEvent(input$restart, {
      for (nm in c("geno", "snp_map", "qc", "filtered", "build_result",
                   "SP", "run_dir", "pca", "trait_specs", "logger"))
        rv[[nm]] <- NULL
      rv$screen <- "welcome"
    })
  })
}
