library(shiny)
library(ggplot2)

setwd("path to folder 'display app'")
if (!"df" %in% ls()) {
  # load each data file and combine them in a single array
  # because max file size is 25MB on github
  df = array(NA,c(4,3,3,4,3,5000,10),
             list(obs=c("250","500","1000","2000"),
                  confdg_lvl=c("low","moderate","high"),
                  trt_rarity=c("common","rare","very rare"),
                  method=c("IPW","KOM","EB","TLF"),
                  model=c("true","misspecified","rf"),
                  NULL,
                  c("ATE", "ATE.dr", "ATE.var", "ATE.alt", "ATE.alt.var", "ATT",
                    "ATT.dr", "ATT.var", "ATT.alt", "ATT.alt.var")))
  
  for(i in 0:4) {
    load(paste0("data",i+1,".Rdata"))
    df[,,,,,1000*i+1:1000,] = data
    rm(data)
  }
}

# See above for the definitions of ui and server
ui <- fluidPage(
  
  # Sidebar layout with input and output definitions ----
  column(
    width = 3,
    h3("Scenario"),
    fluidRow(
      column(
        width = 6,
        radioButtons(
          "trt_rarity",
          label = "Treatment Rarity:",
          choices = list(
            "Common (~40%)" = "common",
            "Rare (~15%)" = "rare",
            "Very rare (~5%)" = "very rare"
          )
        ),
        radioButtons(
          "confdg_lvl",
          label = "Confounding Level:",
          choices = list(
            "Low" = "low",
            "Average" = "moderate",
            "High" = "high"
          )
        )
      ),
      column(
        width = 6,
        radioButtons(
          "obs",
          label = "Sample Size:",
          choices = list(
            "n = 250" = "250",
            "n = 500" = "500",
            "n = 1000" = "1000",
            "n = 2000" = "2000"
          )
        )
      )
    ),
    h3("Estimation Method"),
    fluidRow(
      column(
        width = 6,
        radioButtons(
          "stat",
          label = "Estimand:",
          choices = list("ATE" = "ATE",
                         "ATT" = "ATT"),
          
          selected = "ATE"
        )
      ),
      column(
        width = 6,
        checkboxGroupInput(
          "approche",
          label = "Estimator:",
          choices = list(
            "Weighted average (WA)" = "",
            "Augmented weighted average (AWA)" = ".dr",
            "Weighted linear regression (WLR)" = ".alt"
          ),
          selected = ""
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        checkboxGroupInput(
          "method",
          label = "Balancing Mehtod:",
          choices = list(
            "IPTW" = "IPW",
            "KOM" = "KOM",
            "EB" = "EB",
            "TLF" = "TLF"
          ),
          selected = c("IPW", "KOM", "EB", "TLF")
        )
      ),
      column(
        width = 6,
        checkboxGroupInput(
          "learner",
          label = "Regressor:",
          choices = list(
            "Well-specified logistic reg." = "true",
            "Misspecified logistic reg." = "misspecified",
            "Random forest" = "rf"
          ),
          selected = "true"
        )
      )
    ),
    h3("Display parameter:"),
    fluidRow(
      sliderInput("ylim", label = "Window's Size", min = -2, 
                  max = 2, value = c(-.3,.3), step = .1)
    )
    
  ),
  column(width = 8,
         plotOutput(outputId = "distPlot",height = "800px"))
)
alpha = .6
meth_to_col_alpha = c("IPW" = "#E74C3C", "KOM" = "#F1C40F",
                      "EB"  = "#2E86C1", "TLF" = "#229954")

modl_to_pch = c("true" = 0, "misspecified" = 2, "rf" = 1)

server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    n = sapply(c("method","approche","learner"),
               function(s) length(input[[s]])) |> prod()
    data = list()
    itr = 0
    col = rep("", n)
    pch = rep(NA, n)
    name = rep(NA,n)
    for (met in input$method) {
      for (appr in input$approche) {
        for (lrnr in input$learner) {
          itr = itr + 1
          name[itr] = paste(ifelse(met=="IPW","IPTW",met),
                       switch(appr,"WA",".dr"="AWA", ".alt"="WLR"))
          data[[itr]] <- df[input$obs,input$confdg_lvl,input$trt_rarity,met,
                             lrnr,,paste0(input$stat,appr)]
          col[itr] = meth_to_col_alpha[met]
          pch[itr] = modl_to_pch[lrnr]
        }
      }
    }
    par(mar=c(2,5,2,0))
    boxplot(data, ylab = "", ylim = input$ylim,
         xlab = "", horizontal=TRUE, las = 1, col = col, pch =  pch,
         names = name)
    legend("topleft",c("Well-specified log. reg.","Misspecified log. reg."," Random forest"),pch=modl_to_pch)
    abline(v = 0, col = "red")
    
  })
  
}

shinyApp(ui = ui, server = server)