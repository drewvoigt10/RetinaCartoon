library(shiny)
library(RetinaCartoon)
# Define UI ----
ui <- fluidPage(
  fluidRow(
    column(9,
         plotOutput("cartoon", height = "750px"),
         selectInput(inputId = "dataset",
                     label = "Select a dataset for expression visualization",
                     choices = c("", "all_retina_rpe_chor", "retina_fov_vs_periph", "retina_AIR_vs_control", "RPE_choroid_unselected",
                                 "RPE_choroid_CD31_selected", "retina_fovea_perifovea",  "CD31_choroid_infant_adult"),
                     selected = ""),
         selectizeInput(inputId = "gene",
                        label = "Select a gene to plot",
                        choices = sort(c("", "RHO", "PDE6H", "PDE6A",
                                         "ABCA4", "THY1", "ONECUT2",
                                         "CA4", "VWF", "GRM6", "C1QA",
                                         "CD3D")),
                        selected = "",
                        multiple = FALSE,
                        options= list(maxOptions = 100)
                        ),
         p("Only a subset of genes are included in this basic hosted shiny application. For exploration of all genes, visit our website ",
           a("Spectacle",
             href = "https://singlecell-eye.org")))
  )
)


# Define server logic ----
server <- function(session, input, output) {


  output$cartoon <- renderPlot({
    req(input$gene != "")
    req(input$dataset != "")

    p <- RetinaCartoon::generate_cartoon_data(gene = input$gene, dataset = input$dataset) %>% RetinaCartoon::plot_cartoon()
    p
  })




}

# Run the app ----
shinyApp(ui = ui, server = server)
