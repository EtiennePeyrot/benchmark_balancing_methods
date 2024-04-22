library(shiny)
library(ggplot2)
if (!"df" %in% ls()) load("improved data.Rdata")

# See above for the definitions of ui and server
ui <- fluidPage(
  
  # Sidebar layout with input and output definitions ----
  column(
    width = 3,
    h3("La statistique"),
    fluidRow(
      column(
        width = 6,
        radioButtons(
          "stat",
          label = "Population cible :",
          choices = list("ATE" = "ATE",
                         "ATT" = "ATT"),
          
          selected = "ATE"
        )
      ),
      column(
        width = 6,
        checkboxGroupInput(
          "approche",
          label = "Approche :",
          choices = list(
            "Simple" = "",
            "Double robuste" = ".dr",
            "Alternative" = ".alt"
          ),
          selected = ""
        )
      )
    ),
    h3("Le jeu de donné"),
    fluidRow(
      column(
        width = 6,
        radioButtons(
          "trt_rarity",
          label = "Rareté du traitement :",
          choices = list(
            "Commun (~50%)" = "common",
            "Rare (~20%)" = "rare",
            "Très rare (~5%)" = "very rare"
          )
        ),
        radioButtons(
          "obs",
          label = "Taille de l'échantillon :",
          choices = list(
            "250" = "250",
            "500" = "500",
            "1000" = "1000",
            "2000" = "2000"
          )
        )
      ),
      column(
        width = 6,
        radioButtons(
          "confdg_lvl",
          label = "Niveau de confusion :",
          choices = list(
            "Faible" = "low",
            "Moyen" = "moderate",
            "Fort" = "high"
          )
        )
      )
    ),
    h3("Les méthodes"),
    fluidRow(
      column(
        width = 6,
        checkboxGroupInput(
          "method",
          label = "Méthode :",
          choices = list(
            "IPW" = "IPW",
            "KOM" = "KOM",
            "EB" = "EB",
            "TLF" = "TLF"
          ),
          selected = c("AIPW", "KOM", "EB", "TLF")
        )
      ),
      column(
        width = 6,
        checkboxGroupInput(
          "learner",
          label = "Learners :",
          choices = list(
            "Log reg" = "true",
            "Misspecified" = "misspecified",
            "Random forest" = "rf"
          ),
          selected = "true"
        )
      )
    ),
    h3("Parametre affichage"),
    fluidRow(
      sliderInput("ylim", label = "Limites", min = -2, 
                  max = 2, value = c(-.3,.3), step = .1)
    )
    
  ),
  column(width = 8,
         plotOutput(outputId = "distPlot",height = "800px"))
)
alpha = .6
meth_to_col_alpha = c("IPW" = rgb(231/255, 076/255, 060/255, alpha), # "#E74C3C",
                      "KOM" = rgb(244/255, 208/255, 063/255, alpha), # "#F4D03F",
                      "EB"  = rgb(052/255, 152/255, 219/255, alpha), # "#3498DB",
                      "TLF" = rgb(146/255, 204/255, 113/255, alpha)) # "#2ECC71"

modl_to_pch = c("true" = 0, "misspecified" = 2, "rf" = 1)

server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    n = sapply(c("method","approche","learner"),
               function(s) length(input[[s]])) |> prod()
    data = list()
    itr = 0
    col = rep("", n)
    pch = rep(NA, n)
    for (lrnr in input$learner) {
      for (met in input$method) {
        for (appr in input$approche) {
          itr = itr + 1
          data[[paste(lrnr,met,appr)]] <- df[
            input$obs,input$confdg_lvl,input$trt_rarity,met,lrnr,,paste0(input$stat,appr)]
          col[itr] = meth_to_col_alpha[met]
          pch[itr] = modl_to_pch[lrnr]
        }
      }
    }
    
    par(mar=c(2,10,5,0))
    boxplot(data, ylab = "", ylim = input$ylim,
         xlab = "", horizontal=TRUE, las = 1, col = col, pch =  pch,
         main = paste("Boxplot of balancing method for\nsample size =", input$obs,
                      input$trt_rarity,"treatment &",input$confdg_lvl,"confounding"))
    abline(v = 0, col = "red")
    
  })
  
}

shinyApp(ui = ui, server = server)

