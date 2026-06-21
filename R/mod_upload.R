# Wizard step 1: upload & configure import ----------------------------------

#' @noRd
mod_upload_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header(shiny::h3("1. Upload your genetic data")),
    bslib::card_body(
      shiny::p(class = "text-muted",
               "Choose where your data comes from, then load it."),

      shiny::radioButtons(
        ns("source"), "Data source",
        choices = c("Upload my own data" = "own",
                    "Try a sample dataset" = "sample"),
        selected = "own", inline = TRUE
      ),

      shiny::conditionalPanel(
        sprintf("input['%s'] == 'own'", ns("source")),
        shinyFiles::shinyFilesButton(
          ns("file"), "Choose genotype file", "Select a genotype file",
          multiple = FALSE, class = "btn-outline-primary"
        ),
        shiny::div(class = "mt-2 text-muted", shiny::textOutput(ns("file_path")))
      ),

      shiny::conditionalPanel(
        sprintf("input['%s'] == 'sample'", ns("source")),
        shiny::selectInput(ns("sample"),
          ff_tooltip("Sample dataset",
            "Small example datasets (300 markers x 30 individuals, 3 chromosomes) bundled with the app."),
          choices = c("Numeric matrix (.rds)" = "sample_numeric.rds",
                      "Numeric matrix (.csv)" = "sample_numeric.csv",
                      "VCF (.vcf)" = "sample.vcf",
                      "HapMap (.hmp.txt)" = "sample.hmp.txt"))
      ),

      shiny::selectInput(
        ns("format"), ff_tooltip("File format",
          "Auto-detect uses the file extension. Override if your file has an unusual name."),
        choices = c("Auto-detect" = "auto", "VCF (.vcf/.vcf.gz)" = "vcf",
                    "HapMap (.hmp.txt)" = "hapmap",
                    "Numeric matrix (CSV/RDS/RData)" = "numeric",
                    "PLINK (.raw or .ped/.map)" = "plink"),
        selected = "auto"
      ),

      shiny::uiOutput(ns("detected")),

      shiny::conditionalPanel(
        sprintf("input['%s'] == 'numeric'", ns("format")),
        shiny::hr(),
        shiny::strong(ff_tooltip("Numeric matrix options",
          "Usually optional: FounderForge auto-detects 'chrom'/'pos'/'rs#' columns and ignores standard HapMap metadata. Use these only to override.")),
        shiny::textInput(ns("chr_col"), "Chromosome column name", ""),
        shiny::textInput(ns("pos_col"), "Position (bp) column name", ""),
        shiny::textInput(ns("id_col"), "Marker id column name (optional)", ""),
        shiny::checkboxInput(ns("clean_names"),
          "Clean sample names (keep text before first ':')", FALSE)
      ),

      shiny::hr(),
      shiny::numericInput(ns("ploidy"),
        ff_tooltip("Ploidy", "Number of chromosome copies. 2 for diploid species."),
        value = 2, min = 1, max = 8, step = 1),

      shiny::radioButtons(
        ns("strategy"),
        ff_tooltip("Inbreeding strategy",
          "Outbred: keeps heterozygosity by splitting each genotype into two haplotypes (landraces, outbred panels). Inbred (DH): one haplotype per line, heterozygotes coerced to homozygous (doubled-haploid / fully inbred lines). importInbredGeno: AlphaSimR's convenience importer for fully inbred panels. Picking 'inbred' for truly outbred data destroys heterozygosity; picking 'outbred' for true inbred lines invents fake phasing."),
        choices = c("Outbred / heterozygous (newMapPop, inbred = FALSE)" = "outbred",
                    "Inbred / doubled-haploid (newMapPop, inbred = TRUE)" = "inbred",
                    "Fully inbred panel (importInbredGeno)" = "import_inbred"),
        selected = "outbred"
      ),

      shiny::uiOutput(ns("status")),
      shiny::uiOutput(ns("preview")),

      shiny::div(
        class = "d-flex justify-content-between mt-4",
        shiny::actionButton(ns("back"), shiny::tagList(
          bsicons::bs_icon("arrow-left"), " Back"), class = "btn-outline-secondary"),
        shiny::actionButton(ns("load"), shiny::tagList(
          bsicons::bs_icon("search"), " Load & validate"), class = "btn-primary")
      ),
      shiny::uiOutput(ns("continue_ui"))
    )
  )
}

#' @noRd
mod_upload_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    roots <- c(Home = path.expand("~"), Root = "/")
    shinyFiles::shinyFileChoose(input, "file", roots = roots, session = session)

    file_sel <- shiny::reactiveVal(NULL)  # path from the user's own file picker
    chosen <- shiny::reactiveVal(NULL)    # the path that will actually be read
    loaded <- shiny::reactiveVal(NULL)

    # Clear any previous validation card / preview / continue button so a new
    # file, format or source can't be confused with the last one.
    clear_results <- function() {
      loaded(NULL)
      rv$geno <- NULL; rv$snp_map <- NULL; rv$qc <- NULL; rv$filtered <- NULL
      output$status <- shiny::renderUI(NULL)
      output$preview <- shiny::renderUI(NULL)
      output$continue_ui <- shiny::renderUI(NULL)
    }
    set_chosen <- function(p) { chosen(p); clear_results() }
    sample_path <- function() app_sys("extdata", input$sample)

    # Changing the declared format also invalidates a prior validation.
    shiny::observeEvent(input$format, clear_results(), ignoreInit = TRUE)

    shiny::observeEvent(input$file, {
      info <- shinyFiles::parseFilePaths(roots, input$file)
      if (nrow(info) > 0) {
        file_sel(as.character(info$datapath[1]))
        if (identical(input$source, "own")) set_chosen(file_sel())
      }
    })

    # Switching source, or changing the sample, updates the single chosen path.
    shiny::observeEvent(input$source, {
      if (identical(input$source, "sample")) set_chosen(sample_path())
      else set_chosen(file_sel())
    })
    shiny::observeEvent(input$sample, {
      if (identical(input$source, "sample")) set_chosen(sample_path())
    })

    output$file_path <- shiny::renderText({
      if (is.null(file_sel())) "No file selected yet." else paste("Selected:", file_sel())
    })

    # Live detected-filetype feedback (before loading).
    output$detected <- shiny::renderUI({
      shiny::req(chosen())
      ext_fmt <- detect_format(chosen())
      fmt <- if (input$format == "auto") ext_fmt else input$format
      known <- fmt != "unknown"
      icon <- if (known) "check-circle-fill" else "exclamation-triangle-fill"
      cls  <- if (known) "text-info" else "text-warning"
      note <- if (input$format == "auto") {
        if (known) sprintf("Auto-detected from the file extension as %s.", toupper(fmt))
        else "Could not recognise this extension - please pick a format above."
      } else sprintf("Format set manually to %s.", toupper(fmt))
      shiny::div(class = "mt-2", shiny::span(class = cls,
        bsicons::bs_icon(icon), " Detected file type: ",
        shiny::strong(toupper(fmt)), " - ", note))
    })

    shiny::observeEvent(input$back, rv$screen <- "welcome")

    shiny::observeEvent(input$load, {
      shiny::req(chosen())
      fmt <- input$format
      opts <- list(
        chr_col = if (nzchar(input$chr_col %||% "")) input$chr_col else NULL,
        pos_col = if (nzchar(input$pos_col %||% "")) input$pos_col else NULL,
        id_col  = if (nzchar(input$id_col %||% "")) input$id_col else NULL,
        clean_sample_names = isTRUE(input$clean_names)
      )
      logger <- new_run_logger()
      logger$section("Import")
      logger$log(sprintf("User selected file: %s", chosen()))

      shinybusy::show_modal_spinner(
        spin = "fading-circle", color = "#1D9E75",
        text = "Reading and validating data...")
      on.exit(shinybusy::remove_modal_spinner(), add = TRUE)

      result <- tryCatch(
        read_genotypes(chosen(), format = fmt, opts = opts, logger = logger),
        error = function(e) e
      )

      if (inherits(result, "error")) {
        loaded(NULL)
        output$preview <- shiny::renderUI(NULL)
        output$continue_ui <- shiny::renderUI(NULL)
        output$status <- shiny::renderUI(
          bslib::card(class = "border-danger mt-3", bslib::card_body(
            shiny::tags$strong(class = "text-danger",
              bsicons::bs_icon("exclamation-triangle-fill"), " Could not read the file"),
            shiny::p(conditionMessage(result)),
            shiny::p(class = "text-muted mb-0",
              "Check that the format above matches your file, or try a sample dataset."))))
        return(invisible())
      }

      # Stash everything for downstream steps.
      rv$geno <- result$geno
      rv$snp_map <- result$snp_map
      rv$ploidy <- input$ploidy
      rv$strategy <- input$strategy
      rv$cm_per_mb <- 1
      rv$logger <- logger
      rv$source_info <- list(path = chosen(),
                             format = if (fmt == "auto") detect_format(chosen()) else fmt)
      rv$qc <- NULL; rv$filtered <- NULL; rv$build_result <- NULL
      rv$SP <- NULL; rv$run_dir <- NULL
      loaded(result)

      output$status <- shiny::renderUI(validity_card(result, rv$source_info$format))
      output$preview <- shiny::renderUI(
        bslib::card(class = "mt-3",
          bslib::card_header(ff_tooltip("Data preview",
            "First 20 markers and 20 individuals only - scroll to see more. The full file is used for analysis.")),
          shiny::div(style = "overflow:auto;", DT::DTOutput(ns("preview_table")))))
      output$preview_table <- DT::renderDT(preview_table(result), server = FALSE)
      output$continue_ui <- shiny::renderUI(
        shiny::div(class = "d-flex justify-content-end mt-3",
          shiny::actionButton(ns("continue"), shiny::tagList(
            "Continue to review ", bsicons::bs_icon("arrow-right")),
            class = "btn-success btn-lg")))
    })

    shiny::observeEvent(input$continue, rv$step <- 2L)
  })
}

#' Validity summary card shown after a successful load
#' @noRd
validity_card <- function(result, fmt) {
  geno <- result$geno
  n_chr <- length(unique(result$snp_map$chr))
  n_na <- sum(is.na(geno))
  bslib::card(class = "border-success mt-3", bslib::card_body(
    shiny::tags$strong(class = "text-success",
      bsicons::bs_icon("check-circle-fill"), " File loaded and validated"),
    shiny::tags$ul(class = "mb-0 mt-2",
      shiny::tags$li(sprintf("Detected file type: %s", toupper(fmt))),
      shiny::tags$li(sprintf("Genotypes: %s markers x %s individuals",
        format(nrow(geno), big.mark = ","), format(ncol(geno), big.mark = ","))),
      shiny::tags$li(sprintf("Chromosomes: %d", n_chr)),
      shiny::tags$li(sprintf("Coding check: values in {0, 1, 2%s} - OK",
        if (n_na > 0) ", NA" else "")),
      if (n_na > 0) shiny::tags$li(sprintf(
        "Missing genotypes: %s (will be imputed before building)",
        format(n_na, big.mark = ","))))))
}

#' Build a capped, scrollable preview table (first 20 markers x 20 individuals)
#' @noRd
preview_table <- function(result, n_row = 20L, n_col = 20L) {
  geno <- result$geno
  nr <- min(n_row, nrow(geno))
  nc <- min(n_col, ncol(geno))
  sub <- geno[seq_len(nr), seq_len(nc), drop = FALSE]
  df <- data.frame(
    marker = result$snp_map$id[seq_len(nr)],
    chr = result$snp_map$chr[seq_len(nr)],
    pos = result$snp_map$pos[seq_len(nr)],
    as.data.frame(sub, check.names = FALSE),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  cap <- sprintf("Showing %d of %s markers and %d of %s individuals.",
                 nr, format(nrow(geno), big.mark = ","),
                 nc, format(ncol(geno), big.mark = ","))
  DT::datatable(
    df, rownames = FALSE, caption = cap,
    options = list(scrollX = TRUE, scrollY = "320px", paging = FALSE,
                   searching = FALSE, info = FALSE, ordering = FALSE)
  )
}
