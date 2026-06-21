# Generates a gently-waved honeycomb SVG used as the homepage background.
# It is drawn flat here; the app tilts it in 3D (CSS rotateX) for perspective.
s <- 22                 # hex radius
W <- 1400; H <- 900
w <- sqrt(3) * s        # horizontal centre spacing
vstep <- 1.5 * s        # vertical centre spacing
amp <- 13; wl <- 540    # gentle horizontal undulation
ang <- (c(30, 90, 150, 210, 270, 330)) * pi / 180
dx <- s * cos(ang); dy <- s * sin(ang)

paths <- character()
ys <- seq(-vstep, H + vstep, by = vstep)
for (ri in seq_along(ys)) {
  cy <- ys[ri]
  xoff <- if (ri %% 2 == 0) w / 2 else 0
  for (cx in seq(-w + xoff, W + w, by = w)) {
    vx <- cx + dx
    vy <- cy + dy + amp * sin(2 * pi * (cx + dx) / wl)  # wave shared across neighbours
    paths <- c(paths, paste0("M", paste(sprintf("%.1f,%.1f", vx, vy), collapse = "L"), "Z"))
  }
}
svg <- paste0(
  '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ', W, ' ', H, '">',
  '<path d="', paste(paths, collapse = ""),
  '" fill="none" stroke="#6cc0b0" stroke-width="1" stroke-linejoin="round"/></svg>'
)
writeLines(svg, "inst/app/www/hex_bg.svg")
cat("Wrote inst/app/www/hex_bg.svg :", round(file.size("inst/app/www/hex_bg.svg")/1024), "KB,",
    length(paths), "hexagons\n")
