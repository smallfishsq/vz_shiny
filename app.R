# app.R

# Load required libraries
library(shiny)

# Source the R code file
source("functions.R", local = TRUE)

nhanes_comb_feat <- read.csv('nhanes_comb_feat.csv')

cox_model <- nhanesdesignmodel <- readRDS("nhanesdesignmodel.rds") # load cox regression model
nhane_xgb <- read.csv("nhane_xgb.csv") # this data is used to train cox reg model (nhanesdesignmodel), but it's still needed for survfit function

# Define UI
ui <- fluidPage(
  titlePanel("Cox Proportional Hazards Model Predictor"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("input_age", "Input Age:",
                  min = min(25), max = max(85), value = 50),
      selectInput("input_gender", "Select Gender:", choices = c("Male", "Female")),
      
      # actionButton("plot_button", "Plot")
    ),
    mainPanel(
      plotOutput("survival_plot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  # Plot survival curve
  output$survival_plot <- renderPlot({
    
    new_data <- nhanes_comb_feat[1088, ]
    # print(new_data)
    new_data$Age <- input$input_age
    new_data$c_Age <- input$input_age

    new_data$Gender <-  input$input_gender
    new_data$c_Gender <- ifelse(input$input_gender == "Male", 1, 0)

    # newdataxgb <- xgboostcox(new_data, test='Test')
    newdataxgb <- xgboost_test(new_data)
    print(newdataxgb)
    survfit(nhanesdesignmodel, newdata = newdataxgb) %>%
      ggsurvplot(data = nhanes_comb_feat, palette = "Dark2", conf.int = TRUE
                 ) +
      ggtitle("Survival Curve with Confidence Interval")
  })
}

# Run the application
shinyApp(ui = ui, server = server)