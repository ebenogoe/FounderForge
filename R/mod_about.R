# About screen --------------------------------------------------------------

#' @noRd
mod_about_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header(shiny::h3(shiny::tagList(
      bsicons::bs_icon("info-circle-fill"), " About"))),
    bslib::card_body(
      shiny::p("FounderForge is a point-and-click companion for ",
               shiny::tags$a(href = "https://gaynorr.github.io/AlphaSimR/", "AlphaSimR"),
               ", built so breeders and students can prepare real marker data as a founder population for simulations without coding."),
      shiny::h5("How it works"),
      shiny::p("Your data is read into a 0/1/2 genotype matrix, quality-checked and filtered, then converted to AlphaSimR haplotypes and a genetic map and passed to newMapPop() (or importInbredGeno()). An optional PCA confirms the import preserved population structure."),
      shiny::h5("Found a bug?"),
      shiny::p("Please report issues with a copy of your run log (logs/ folder) to ",
               shiny::tags$a(href = "mailto:ebenezerogoe@gmail.com", "ebenezerogoe@gmail.com"), "."),
      shiny::p(class = "text-muted", sprintf("Version %s", utils::packageVersion("FounderForge"))),
      shiny::actionButton(ns("home"), shiny::tagList(
        bsicons::bs_icon("house"), " Back to home"), class = "btn-outline-secondary mt-2")
    )
  )
}

#' @noRd
mod_about_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$home, rv$screen <- "welcome")
  })
}
