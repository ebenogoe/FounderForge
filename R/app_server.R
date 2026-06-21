#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  rv <- shiny::reactiveValues(
    screen = "welcome",
    step = 1L,
    geno = NULL, snp_map = NULL, qc = NULL, filtered = NULL,
    build_result = NULL, SP = NULL, trait_specs = list(),
    logger = NULL, source_info = list(), run_dir = NULL, pca = NULL,
    ploidy = 2L, strategy = "outbred", cm_per_mb = 1
  )

  # Screen modules
  mod_welcome_server("welcome", rv)
  mod_upload_server("upload", rv)
  mod_review_server("review", rv)
  mod_filter_server("filter", rv)
  mod_traits_server("traits", rv)
  mod_build_export_server("build", rv)
  mod_help_server("help", rv)
  mod_about_server("about", rv)

  # Top-bar navigation
  shiny::observeEvent(input$brand_home, rv$screen <- "welcome")
  shiny::observeEvent(input$nav_home,  rv$screen <- "welcome")
  shiny::observeEvent(input$nav_help,  rv$screen <- "help")
  shiny::observeEvent(input$nav_about, rv$screen <- "about")

  # Drive the top-level screen navset
  shiny::observeEvent(rv$screen, {
    bslib::nav_select("main_nav", rv$screen, session = session)
  })

  # Drive the wizard step navset + stepper highlight
  shiny::observeEvent(rv$step, {
    bslib::nav_select("wizard_nav", paste0("s", rv$step), session = session)
  })
  output$wizard_stepper <- shiny::renderUI(ff_stepper(rv$step))
}
