### Block without UI and SERVER ###
### To launch from console:  R -e "shiny::runApp('.')"

# Small helper functions
plot_circ <- function(x, plt_rose, ...){
  suppressWarnings(x <- circular::as.circular(x))
  plot(x, cex = 1., bin = 720, stack = TRUE, sep = 0.035, shrink = 1.3, ...)
  if(plt_rose){
    circular::rose.diag(x, bins = 16, col = "darkgrey", cex = 1., prop = 1.3, add = TRUE)
  }
  lines(circular::density.circular(x, bw = 1), lty = 2)
}

create_mus <- function(type, n, mu0, alpha0, beta0, beta1, fcst.alpha0, fcst.beta0, fcst.beta1){
  if(type == "linear"){ 
    x1 <- seq(mu0 -1, mu0 + 1, length.out = n)
    mu.obs <- 2 * atan(alpha0) + 2 * atan(beta0 + beta1 %*% x1)
    mu.fcst <- 2 * atan(fcst.alpha0) + 2 * atan(fcst.beta0 + fcst.beta1 %*% x1)

    x1.plot <- seq(-15, 15, 0.1)
    mu.obs.plt <- 2 * atan(alpha0) + 2 * atan(beta0 + beta1 %*% x1.plot)
    mu.fcst.plt <- 2 * atan(fcst.alpha0) + 2 * atan(fcst.beta0 + fcst.beta1 %*% x1.plot)
  } else if (type == "circular"){
    x1 <- seq(-pi, pi, length.out = n)
    mu.obs <- 2 * atan(alpha0) + 2 * atan(beta0 + beta1 %*% x1)
    mu.fcst <- 2 * atan(fcst.alpha0) + 2 * atan(fcst.beta0 + fcst.beta1 %*% x1)

    x1.plot <- seq(-pi, pi, 0.1)
    mu.obs.plt <- 2 * atan(alpha0) + 2 * atan(beta0 + beta1 %*% x1.plot)
    mu.fcst.plt <- 2 * atan(fcst.alpha0) + 2 * atan(fcst.beta0 + fcst.beta1 %*% x1.plot)
  } else if (type == "circularlink"){
    x1 <- seq(-pi, pi, length.out = n)
    mu.obs <- 2 * atan(alpha0) + 2 * atan(beta0 + beta1 %*% tan(x1/2))
    mu.fcst <- 2 * atan(fcst.alpha0) + 2 * atan(fcst.beta0 + fcst.beta1 %*% tan(x1/2))

    x1.plot <- seq(-pi, pi, 0.1)
    mu.obs.plt <- 2 * atan(alpha0) + 2 * atan(beta0 + beta1 %*% tan(x1.plot/2))
    mu.fcst.plt <- 2 * atan(fcst.alpha0) + 2 * atan(fcst.beta0 + fcst.beta1 %*% tan(x1.plot/2))
  }
  rval <- list(mu.obs, mu.fcst, mu.obs.plt, mu.fcst.plt, x1, x1.plot)
  return(rval)
}

# Set values for all script
col_hcl <- c(gray(0, alpha = 0.8), colorspace::rainbow_hcl(2, alpha =0.8))

### UI ###

# Define UI
ui <- shiny::fixedPage(
  shiny::tags$head(shiny::tags$script(shiny::HTML('
      Shiny.addCustomMessageHandler("jsCode",
        function(message) {
          eval(message.code);
        }
      );
    '))),

  # App title
  shiny::titlePanel("Von Mises Distribution"),
  
  # Sidebar layout with input and output definitions ----
  shiny::fixedRow(

    # Sidebar panel for inputs ----
    shiny::column(4,

      shiny::wellPanel(

        #textOutput(outputId = "formula"),
        shiny::withMathJax(shiny::textOutput(outputId = "formula")),

        shiny::hr(),

        shiny::radioButtons(inputId = "type", "",
              c("Linear Covariable" = "linear",
                "Circular Covariable w/o link" = "circular",
                "Circular Covariable w/ link" = "circularlink")),

        shiny::fixedRow(
          shiny::column(4,
            shiny::numericInput(inputId = "n", 
                         label = "samples", 
                         value = 100,
                         width = "100px")
          ),
          shiny::column(4,
            shiny::numericInput(inputId = "mu0", 
                         label = "mean(x1)", 
                         min = -10,
                         max = +10,
                         value = 0,
                         width = "100px")
          )
        ),

        shiny::fixedRow(
         shiny::column(4,
            shiny::numericInput(inputId = "alpha0", 
                         label = "alpha0", 
                         value = 0,
                         min = -10, 
                         step = 0.5,
                         max = 10,
                         width = "60px")
          ),
          shiny::column(4,
            shiny::numericInput(inputId = "beta0", 
                         label = "beta0", 
                         value = 0,
                         min = -10, 
                         step = 0.5,
                         max = 10,
                         width = "60px")
          ),
          shiny::column(4,
            shiny::numericInput(inputId = "beta1", 
                         label = "beta1", 
                         value = 1,
                         min = -20, 
                         max = 20,
                         width = "60px")
          )
        ),

        shiny::actionButton(inputId = "startval",
                     label = "Reset to default parameters"),

        shiny::hr(),

        shiny::sliderInput(inputId = "fcst.alpha0",
                    label = "Misspecified alpha0:",
                    min = -10,
                    max = 10,
                    step = 0.5,
                    value = 0),

        shiny::sliderInput(inputId = "fcst.beta0",
                    label = "Misspecified beta0:",
                    min = -10,
                    max = 10,
                    step = 0.5,
                    value = 0),

        shiny::sliderInput(inputId = "fcst.beta1",
                    label = "Misspecified beta1:",
                    min = -20,
                    max = 20,
                    value = 1),

        shiny::actionButton(inputId = "reset",
                     label = "Reset to true parameters"),

        shiny::hr(),
        
        shiny::checkboxInput("rosediag", "Plot rose diagram", TRUE)

      )

    ),
    shiny::column(8,
      shiny::fixedRow(
        shiny::column(6,
        # Main panel for displaying outputs ----

          # Output: plot obs
          shiny::plotOutput(outputId = "obsPlot")

        ), 

        shiny::column(6,

          # Output: plot fcst 
          shiny::plotOutput(outputId = "fcstPlot")

        )
      ),

      shiny::hr(),

      shiny::plotOutput(outputId = "transPlot")
    )
  )

  

)

### SERVER ###

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  shiny::observe({
    if(input$type != "linear") {
      session$sendCustomMessage(type="jsCode",
                                list(code= "$('#mu0').prop('disabled',true)"))
    } else {
      session$sendCustomMessage(type="jsCode",
                                list(code= "$('#mu0').prop('disabled',false)"))
    }
  })

  shiny::observeEvent(input$reset, {
    #showNotification("")
    shiny::updateSliderInput(session, "fcst.alpha0", value = input$alpha0)
    shiny::updateSliderInput(session, "fcst.beta0", value = input$beta0)
    shiny::updateSliderInput(session, "fcst.beta1", value = input$beta1)

  })

  shiny::observeEvent(input$startval, {
    #showNotification("")
    shiny::updateNumericInput(session, "n", value = 100)
    shiny::updateNumericInput(session, "mu0", value = 0)
    shiny::updateNumericInput(session, "alpha0", value = 0)
    shiny::updateNumericInput(session, "beta0", value = 0)
    shiny::updateNumericInput(session, "beta1", value = 1)

  })

  # Reactive expression 
  current_mu <- shiny::reactive({ 
    
    create_mus(input$type, input$n, input$mu0, input$alpha0, input$beta0, input$beta1, input$fcst.alpha0, input$fcst.beta0, input$fcst.beta1) 

    }) 

  # Generate output
  output$formula <- shiny::renderText({
    #type <- input$type
    #if(type == "linear"){
      print(paste0("VM (location part): $$\\mu = 2 \\cdot atan(\\alpha_0) + 2 \\cdot atan(\\beta_0 + \\beta_1 \\cdot x_1)$$"))
    #} else {
    #  print(paste0("VM (location part): $$\\mu = 2 \\cdot atan(\\alpha_0) + 2 \\cdot atan(\\beta_0 + \\beta_1 \\cdot atan(\\frac{x_1}{2}))$$"))
    #}
  })

  output$obsPlot <- shiny::renderPlot({

    plot_circ(as.numeric(current_mu()[[1]]), plt_rose = input$rosediag, 
      main = sprintf("VM(obs | beta0 = %.1f, beta1 = %.1f, beta2 = %.1f)", 
      input$alpha0, input$beta0, input$beta1), col = col_hcl[1])

    })

  output$fcstPlot <- shiny::renderPlot({

    plot_circ(as.numeric(current_mu()[[2]]), plt_rose = input$rosediag, 
      main = sprintf("VM(fcst | beta0 = %.1f, beta1 = %.1f, beta2 = %.1f)", 
      input$fcst.alpha0, input$fcst.beta0, input$fcst.beta1), col = col_hcl[3])

    })

  output$transPlot <- shiny::renderPlot({

    plot(current_mu()[[6]],  as.numeric(current_mu()[[3]]),
      type = "line", xlab = "x1", ylab = "mu", ylim = c(-2*pi, 2*pi))
    lines(current_mu()[[6]],  as.numeric(current_mu()[[4]])
      , lty = 2, type = "line", xlab = "x1", ylab = "mu", ylim = c(-2*pi, 2*pi))
    points(current_mu()[[5]], as.numeric(current_mu()[[1]]), pch = 20, col = col_hcl[1], cex = 0.5)
    points(current_mu()[[5]], as.numeric(current_mu()[[2]]), pch = 20, col = col_hcl[3], cex = 0.5) 
    rug(as.numeric(current_mu()[[5]]), side = 1, col = col_hcl[1])
    rug(as.numeric(current_mu()[[1]]), side = 2, col = col_hcl[1])
    rug(as.numeric(current_mu()[[2]]), side = 2, col = col_hcl[3])
    legend("topright", c("observation", "misspecified fcst"), pch = 20, col = col_hcl[c(1,3)])

  })  

}

# Complete app with UI and server components
shiny::shinyApp(ui, server)


