# Help & docs screen --------------------------------------------------------

#' @noRd
mod_help_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header(shiny::h3(shiny::tagList(
      bsicons::bs_icon("question-circle-fill"), " Help & docs"))),
    bslib::card_body(
      bslib::accordion(
        open = "What does FounderForge do?",
        bslib::accordion_panel("What does FounderForge do?",
          shiny::p("FounderForge prepares your marker data as an AlphaSimR founder population - the object you use to seed a breeding simulation - without writing any R code.")),
        bslib::accordion_panel("Which file formats are supported?",
          shiny::tags$ul(
            shiny::tags$li(shiny::strong("VCF"), " (.vcf, .vcf.gz) - biallelic SNPs; the GT field is converted to 0/1/2."),
            shiny::tags$li(shiny::strong("HapMap"), " (.hmp.txt) - TASSEL-style nucleotide calls."),
            shiny::tags$li(shiny::strong("Numeric matrix"), " (CSV/RDS/RData) - a 0/1/2 dosage table with a marker map, or a saved list(geno, snp_map)."),
            shiny::tags$li(shiny::strong("PLINK"), " - .raw (from plink --recodeA) or .ped/.map text files."))),
        bslib::accordion_panel("Outbred vs inbred - which do I pick?",
          shiny::p("Pick ", shiny::strong("Outbred"), " for landraces and heterozygous panels; FounderForge splits each genotype into two haplotypes so heterozygosity is preserved."),
          shiny::p("Pick ", shiny::strong("Inbred / doubled-haploid"), " (or importInbredGeno) only when your lines are fully homozygous. Choosing the wrong option distorts the genetic structure.")),
        bslib::accordion_panel("What comes out at the end?",
          shiny::tags$ul(
            shiny::tags$li(shiny::strong("genmap.rds & haplotypes.rds"), " - the exact inputs to AlphaSimR::newMapPop()."),
            shiny::tags$li(shiny::strong("founder_pop.rds"), " (+ sim_param.rds if you defined traits)."),
            shiny::tags$li(shiny::strong("rebuild_snippet.R"), " - the precise call used, so you can reproduce it."),
            shiny::tags$li(shiny::strong("qc/"), ", ", shiny::strong("validation/"), " and ", shiny::strong("logs/"), " folders with plots and a full text log."))),
        bslib::accordion_panel("Numeric matrix / RData column names",
          shiny::p("For a wide numeric table (CSV/RData) FounderForge auto-detects the chromosome, position and marker-id columns (e.g. ", shiny::code("chrom"), ", ", shiny::code("pos"), ", ", shiny::code("rs#"), ") and ignores standard HapMap metadata columns (", shiny::code("alleles, strand, assembly#, REFERENCE_GENOME"), ", etc.)."),
          shiny::p("If your sample names carry suffixes like ", shiny::code("LINE1:FLOWCELL:LANE"), ", tick ", shiny::strong("Clean sample names"), " on the upload page. You can also type exact column names there to override auto-detection.")),
        bslib::accordion_panel("Try a sample dataset",
          shiny::p("Download a small example (300 markers x 30 individuals, 3 chromosomes) in any format, or load it directly from the upload page:"),
          shiny::div(class = "d-flex flex-wrap gap-2",
            shiny::downloadButton(ns("dl_rds"), "Numeric .rds", class = "btn-sm btn-outline-secondary"),
            shiny::downloadButton(ns("dl_csv"), "Numeric .csv", class = "btn-sm btn-outline-secondary"),
            shiny::downloadButton(ns("dl_vcf"), "VCF", class = "btn-sm btn-outline-secondary"),
            shiny::downloadButton(ns("dl_hmp"), "HapMap", class = "btn-sm btn-outline-secondary"))),
        bslib::accordion_panel("A note on approximations",
          shiny::p("The genetic map assumes a uniform recombination rate (1 cM/Mb). QTL positions are placed at random by AlphaSimR, and any heritability you enter is your assumption. These are listed in the run log."))
      ),
      shiny::actionButton(ns("home"), shiny::tagList(
        bsicons::bs_icon("house"), " Back to home"), class = "btn-outline-secondary mt-2")
    )
  )
}

#' @noRd
mod_help_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$home, rv$screen <- "welcome")

    sample_dl <- function(fname) {
      shiny::downloadHandler(
        filename = function() fname,
        content = function(file) file.copy(app_sys("extdata", fname), file)
      )
    }
    output$dl_rds <- sample_dl("sample_numeric.rds")
    output$dl_csv <- sample_dl("sample_numeric.csv")
    output$dl_vcf <- sample_dl("sample.vcf")
    output$dl_hmp <- sample_dl("sample.hmp.txt")
  })
}
