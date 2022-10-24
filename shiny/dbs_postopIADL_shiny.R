# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# ---- Load packages ----

# list packages needed
pkgs <- c(
  "shiny", # to run the app
  "dplyr", "tidyverse", # for data wrangling
  "brms", "bayestestR", # for prediction calculations and summaries 
  "ggplot2", "patchwork" # for plotting
)

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}


# ---- Read objects needed ----

# read model and scaling values
m <- readRDS("stat_model.rds")$m
scl <- readRDS("stat_model.rds")$scl


# ---- Set-up ggplot colors ----

# list a colorblind-friently pallete (cbPal)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )


# ---- Define UI ----

ui <- fluidPage(
  
  sidebarLayout(
    
    # Inputs: select LEDD to predict at
    sidebarPanel(
      
      # Select input for the pre-surgery LEDD
      sliderInput( inputId = "led_0", label = "Pre-surgery LEDD (mg)", min = 0, max = 5000, value = 1820, step = 1),
      numericInput( inputId = "led_n0", label = NULL, min = 0, max = 5000, value = 1820),
      
      # Select input for the post-surgery LEDD
      sliderInput( inputId = "led_1", label = "Post-surgery LEDD (mg)", min = 0, max = 5000, value = 833, step = 1),
      numericInput( inputId = "led_n1", label = NULL, min = 0, max = 5000, value = 833 )
      
    ),
  
    # Output: Show the plot
    mainPanel(
      plotOutput( outputId = "pred_plt", height = 900 )
      )
    
  )
  
)


# ---- Define server ----

server <- function(input, output) {
  
  # update slider according to the number
  observe( { updateSliderInput( inputId = "led_0", value = input$led_n0 ) } )
  observe( { updateSliderInput( inputId = "led_1", value = input$led_n1 ) } )
  
  # update number according to the slider
  observe( { updateSliderInput( inputId = "led_n0", value = input$led_0 ) } )
  observe( { updateSliderInput( inputId = "led_n1", value = input$led_1 ) } )
  
  # compute predictions
  ppred <- reactive( { posterior_epred( object = m, 
                                        newdata = data.frame( post = factor(c(0,1)),
                                                              led = ( c(input$led_0,input$led_1)-scl["led","M"] )/scl["led","SD"] ,
                                                              id = NA, item = NA ),
                                        re_formula = NA ) } )
  
  # plot predictions
  output$pred_plt <- renderPlot(
    { data.frame( resp = rep(1:5,2), post = c( rep(1,5), rep(2,5) ) ) %>% # prepare data
        # add median and boundary (95% HDI) predictions
        mutate( Md = 100 * sapply( 1:nrow(.), function(i) median( ppred()[ , post[i], resp[i] ] ) ),
                CI.low = 100 * sapply( 1:nrow(.), function(i) ci( ppred()[ , post[i], resp[i] ], method = "HDI", ci = .95 ) %>% as.data.frame() %>% select(CI_low) %>% as.numeric() ),
                CI.upp = 100 * sapply( 1:nrow(.), function(i) ci( ppred()[ , post[i], resp[i] ], method = "HDI", ci = .95 ) %>% as.data.frame() %>% select(CI_high) %>% as.numeric() ),
                `Session:` = ifelse( post == 1, "pre-surgery", "post-surgery"), resp = as.factor(resp) ) %>%
        # plot it
        ggplot( aes( x = resp, y = Md, ymin = CI.low, ymax = CI.upp, color = `Session:`, fill = `Session:` ) ) +
        geom_pointrange( shape = 21, size = 3.5, fatten = 2, position = position_dodge( width = .25 ) ) +
        labs( x = element_blank(), y = "Probability of response (%)" ) +
        scale_x_discrete( labels = c('"cannot do"','"a lot"', '"somewhat"','"a little"','"none"') ) +
        scale_y_continuous( limits = c(0,100), breaks = seq(0,100,25), labels = seq(0,100,25) ) +
        scale_color_manual( values = cbPal[c(1,3)] ) +
        scale_fill_manual( values = c("black",cbPal[6]) ) +
        coord_flip() +
        theme_minimal( base_size = 18 ) +
        theme( legend.position = "bottom" )
    }
  )
}


# ---- Create a shiny app object ----

shinyApp(ui = ui, server = server)
