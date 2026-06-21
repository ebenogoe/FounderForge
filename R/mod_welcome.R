# Welcome screen ------------------------------------------------------------

#' @noRd
mod_welcome_ui <- function(id) {
  ns <- shiny::NS(id)
  card_action <- function(action_id, icon, title, subtitle, colour) {
    shiny::actionLink(
      ns(action_id), class = "text-decoration-none",
      bslib::card(
        class = "ff-card-action",
        bslib::card_body(
          shiny::div(bsicons::bs_icon(icon), style = paste0("color:", colour)),
          shiny::h4(title, class = "mt-2"),
          shiny::p(subtitle, class = "text-muted mb-0")
        )
      )
    )
  }
  shiny::div(
    class = "ff-welcome-wrap",
    shiny::div(class = "ff-hex-stage", shiny::div(class = "ff-hex-plane")),
    shiny::div(
      class = "ff-hero",
      shiny::div(class = "ff-hero-icon", bsicons::bs_icon("diagram-3-fill")),
      shiny::h1("FounderForge"),
      shiny::p(class = "lead",
               "Prepare your marker data as a founder population for AlphaSimR simulations - no coding required.")
    ),
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      card_action("get_started", "rocket-takeoff-fill", "Get started",
                  "Open the import wizard", "#1D9E75"),
      card_action("go_help", "question-circle-fill", "Help & docs",
                  "Guides for first-time users", "#185FA5"),
      card_action("go_about", "info-circle-fill", "About",
                  "Credits & reporting bugs", "#BA7517")
    )
  )
}

#' @noRd
mod_welcome_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$get_started, {
      rv$step <- 1L
      rv$screen <- "wizard"
    })
    shiny::observeEvent(input$go_help,  rv$screen <- "help")
    shiny::observeEvent(input$go_about, rv$screen <- "about")
  })
}
