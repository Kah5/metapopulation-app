library(shiny)
library(ggplot2)
library(tidyverse)
library(primer)
ui <- fluidPage(
  
  titlePanel("Metapopulation simulation for Sky Island black bears"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Select values each run of the model"),
      
      sliderInput("Npop", h3("Total Number of simulations to run"),
                  min = 0, max = 100, value = 25),
      
      sliderInput("pe", h3("Pe Mean Extinction Rate"),
                  min = 0, max = 1, value = 0.5),
      
      sliderInput("Ci", h3("Pi Mean Colonization Rate"),
                  min = 0, max = 1, value = 0.5),
      
      sliderInput("npatchs", h3("Total Number of habitat patches"),
                  min = 0, max = 20, value = 10),
      
      sliderInput("Nyears", h3("Total Number of years to run the simulation"),
                  min = 0, max = 100, value = 50),
      
      
      sliderInput("p0", h3("Proportion of initial patches populated"),
                  min = 0, max = 1, value = 0.5),
     
      
      
      sliderInput("D", h3("D Habitat Destruction rate"),
                  min = 0, max = 1, value = 0)),
    
    
    
      mainPanel(
               
                plotOutput("plotproprain"),
                verbatimTextOutput("equilfraction"),
                verbatimTextOutput("Px"),
                
                plotOutput("coresatmodel"),
                
                plotOutput("plot"),
                plotOutput("landemodel" ),
                plotOutput("immigrationextinctionplot")
                
                ))
  
  # selectInput("dataset", label = "Dataset", choices = ls("package:datasets")),
  # verbatimTextOutput("summary"),
  # tableOutput("table")
)


server <- function(input, output, session) {
  output$plot <- renderPlot({
    
    
   # run the levins metapopulation model
    out <- MetaSim(NSims=input$Npop, method = "levins", e = input$pe, p0 = input$p0, ci = input$Ci, Time = input$Nyears)
    pops <- out$Ns
    matplot(out$t, pops, type='l', xlab = "Years", ylab= "Proportion of patches occupied", ylim = c(0,1))
    title(main=expression(paste("p"[i]~"f(1-f) - p"[e]~"f   Levin's internal colonization model")))
    
    if(input$Npop == 1){
      # create a block of 50 grid cells:
      
      
      landscape.mat <- matrix(data = rbinom(n = 25, size = 1, out$Parameters["p0"]), ncol = 5, nrow = 5)
      
      
      landscape.list <- list()
      
      landscape.list[[1]] <- landscape.mat
      
      
      colonize.or.extinct<- function(x){
        if(x == 0){
          # if the plot is not colonized, and the random draw is < the colonization probability, then colonize it
          ifelse(runif(n = 1) < as.vector(out$Parameters["ci"]), 1,0)
        }else{
          # if the plot is colonized, and the random draw is < the extinction probability, then make that population extinct
          ifelse(runif(n = 1) < as.vector(out$Parameters["e"]), 0,1)
        }
      }
      
      
      for(t in 2:input$Nyears){
        landscape.list[[t]]<- apply(landscape.list[[t-1]], 1:2, FUN = colonize.or.extinct)
      }
      
      test <- reshape2::melt(landscape.list)
      ggplot(test, aes(x = Var1, y = Var2, fill = value))+geom_tile()+facet_wrap(~L1)+theme_bw()+theme(axis.title = element_blank())
      
    }
  })
  
  
  output$equilfraction <- renderPrint({
    # f = pi/(pi + pe).
    f <- input$Ci/ (input$Ci + input$pe)
  cat("equlibrium fraction of patches occupied (fhat) = ", f)
  })
  
  output$fractionpatches <- renderPrint({
    # f = pi/(pi + pe).
    f <- input$Ci/ (input$Ci + input$pe)
    cat("equlibrium fraction of patches occupied (fhat) = ", f)
  })
  
  output$Px <- renderPrint({
    # f = pi/(pi + pe).
    Px <- 1- as.vector(input$pe)^as.vector(input$npatchs)
    cat("Probability of global persistance = ", Px)
  })
  
  output$plotproprain <- renderPlot({
   
    
    # run the levins metapopulation model
    out <- MetaSim(NSims=input$Npop, method = "gotelli", e = input$pe, p0 = input$p0, ci = input$Ci, Time = input$Nyears)
    pops <- out$Ns
    matplot(out$t, pops, type='l', ylab= "Proportion of patches occupied", ylim = c(0,1))
    #title(main=paste(out$method, "propagule rain model (Pi*(1-f) - Pe*f*(1-f))"))
    title(main=expression(paste("p"[i]~"(1-f) - p"[e]~"f  Propagule Rain model")))
    
    
  })
  
  output$coresatmodel <- renderPlot({
    
    
    # run the levins metapopulation model
    out <- MetaSim(NSims=input$Npop, method = "hanski", e = input$pe, p0 = input$p0, ci = input$Ci, Time = input$Nyears)
    pops <- out$Ns
    matplot(out$t, pops, type='l', xlab = "Years", ylab= "Proportion of patches occupied", ylim = c(0,1))
    #title(main=paste(out$method, "core-satellite metapopulation model (Pi*f*(1-f) - Pe*f*(1-f))"))
    title(main=expression(paste("p"[i]~"f(1-f) - p"[e]~"f(1-f)   Hanski's Internal colonization + rescue effect model")))
    
    
  })
  
  output$immigrationextinctionplot <- renderPlot({
    
    input.df <- data.frame(Ci = as.vector(input$Ci), 
               Pe = as.vector(input$pe))
    # run the levins metapopulation model
    patchf <- data.frame(f = seq(from=0, to = 1,by = 0.1))
    I.table <- patchf %>% group_by(f) %>% mutate(LevinsI =   input.df$Ci*f*(1-f),
                                                 hanskisI =   input.df$Ci*f*(1-f), 
                                                 proprainI =   input.df$Ci*(1-f))
    
    
    
    
    E.table <- patchf  %>% group_by(f) %>% mutate(levinsE = input.df$Pe*f, 
                                                 hanskisE =   input.df$Pe*f*(1-f), 
                                                 proprainE =   input.df$Pe*f*(1-f))
    
    
    E.df <- reshape2::melt(E.table, id.vars = "f")
    colnames(E.df) <- c("f", "model", "E")
    
    I.df <- reshape2::melt(I.table, id.vars = "f")
    colnames(I.df) <- c("f", "model", "I")
    
    a <- ggplot(E.df, aes(x = f, y = E, group = model, color = model))+geom_line()+ylab("Extinction Rate") + xlab("Fraction of patches occupied")+
         scale_colour_manual(values = c('levinsE' = 'red','hanskisE' = 'green', "proprainE" = 'blue'),
                             name = '', 
                          labels = expression(paste("p"[e]~"(f)"), "p"[e]~"f(1-f)", "p"[e]~"f(1-f)"))
    
    b <- ggplot(I.df, aes(x = f, y = I, group = model, color = model))+geom_line()+ylab("Immigration Rate") + xlab("Fraction of patches occupied")+
      scale_colour_manual(values = c('LevinsI' = 'red','hanskisI' = 'green', "proprainI" = 'blue'),
                          name = '', 
                          labels = expression(paste("p"[i]~"f(1-f)"), "p"[i]~"f(1-f)", "p"[i]~"(1-f)"))
    
    
    cowplot::plot_grid(b, a)
  })
  
  output$landemodel <- renderPlot({
    
    
    # run the levins metapopulation model
    out <- MetaSim(NSims=input$Npop, method = "lande", e = input$pe, p0 = input$p0, ci = input$Ci, Time = input$Nyears, D =input$D)
    pops <- out$Ns
    matplot(out$t, pops, type='l', xlab = "Years", ylab= "Proportion of patches occupied", ylim = c(0,1))
    #title(main =paste( "Levins model with habitat destruction \n (Pi*(f)*(1-f - D) - Pe*(f))"))
    title(main=expression(paste("p"[i]~"f(1-f-D) - p"[e]~"f   Levin's model with habitat destruction")))
    
    
  })
  
}

shinyApp(ui, server)