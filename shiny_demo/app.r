library(shiny)

# Define UI for data upload app ----
ui <- fluidPage(

  # App title ----
  titlePanel("Upload Files for Analysis"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Select a file ----
      fileInput("file1", "Choose counts file",
                multiple = FALSE,
                accept = c("text/tsv",
                         "text/csv",
                         "text/tab-separated-values,text/plain",
                         ".csv",
                         ".tsv")),

      fileInput("file2", "Choose coldata file",
                multiple = FALSE,
                accept = c("text/tsv",
                         "text/csv",
                         "text/tab-separated-values,text/plain",
                         ".csv",
                         ".tsv")),
    
      fileInput("file3", "Choose filter file",
                multiple = FALSE,
                accept = c("text/txt",
                         "text/plain",
                         ".txt")),

      # Horizontal line ----
      tags$hr(),

      # Input: Select Analysis Type ----
      radioButtons("Analysis-Type", "Analysis-Type",
                   choices = c("RNA-Seq" = "RNA-Seq",
                               "Microarray" = "Microarray"),
                   selected = "RNA-Seq"),

      # Horizontal line ----
      tags$hr()

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Data file ----
      tableOutput("contents")

    )

  )
)

# Define server logic to read selected files ----
server <- function(input, output) {

}

# Create Shiny app ----
shinyApp(ui, server)
