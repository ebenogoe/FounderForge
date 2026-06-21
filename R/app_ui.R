#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    golem_add_external_resources(),
    bslib::page_fluid(
      theme = ff_theme(),
      # An animated top bar shows whenever the server is working. Heavy steps
      # additionally raise a centred blocking modal (see the modules).
      shinybusy::add_busy_bar(color = "#1D9E75", height = "5px"),
      ff_topbar(),
      bslib::navset_hidden(
        id = "main_nav",
        bslib::nav_panel_hidden("welcome", mod_welcome_ui("welcome")),
        bslib::nav_panel_hidden("wizard", ff_wizard_ui()),
        bslib::nav_panel_hidden("help", mod_help_ui("help")),
        bslib::nav_panel_hidden("about", mod_about_ui("about"))
      )
    )
  )
}

#' Brand + navigation top bar
#' @noRd
ff_topbar <- function() {
  shiny::div(
    class = "ff-topbar",
    shiny::actionLink("brand_home", class = "ff-brand",
      shiny::tagList(bsicons::bs_icon("diagram-3-fill"), "FounderForge")),
    shiny::span(class = "spacer"),
    shiny::actionLink("nav_home", "Home", class = "btn btn-link"),
    shiny::actionLink("nav_help", "Help", class = "btn btn-link"),
    shiny::actionLink("nav_about", "About", class = "btn btn-link")
  )
}

#' The wizard layout: stepper sidebar + hidden step navset
#' @noRd
ff_wizard_ui <- function() {
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 240, open = "always",
      shiny::h6("Import wizard", class = "text-muted"),
      shiny::uiOutput("wizard_stepper")
    ),
    bslib::navset_hidden(
      id = "wizard_nav",
      bslib::nav_panel_hidden("s1", mod_upload_ui("upload")),
      bslib::nav_panel_hidden("s2", mod_review_ui("review")),
      bslib::nav_panel_hidden("s3", mod_filter_ui("filter")),
      bslib::nav_panel_hidden("s4", mod_traits_ui("traits")),
      bslib::nav_panel_hidden("s5", mod_build_export_ui("build"))
    )
  )
}

#' The bslib theme
#' @noRd
ff_theme <- function() {
  bslib::bs_theme(
    version = 5,
    primary = "#1D9E75",
    secondary = "#185FA5",
    success = "#1D9E75",
    info = "#185FA5",
    warning = "#BA7517",
    font_scale = 1.05,
    base_font = bslib::font_collection(
      bslib::font_google("Inter", local = FALSE), "system-ui", "sans-serif")
  )
}

#' Add external Resources to the Application
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path("www", app_sys("app/www"))
  tags$head(
    tags$link(rel = "icon", type = "image/svg+xml", href = "www/favicon.svg"),
    bundle_resources(path = app_sys("app/www"), app_title = "FounderForge")
  )
}
