library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("Age-structured Pacific salmon population model"),
  
  sidebarPanel(
    radioButtons("model_type","Choose model type:",
                 list("Model I" ="model1", 
                      "Model II"="model2")),
    conditionalPanel(
      condition="input.model_type=='model2'",
      numericInput("c_value", "c: Choose the strength of the density dependence:", 0.7)      
    ),
    
    HTML("<hr />"),
    
    numericInput("sj", "Sj: Probability of a juvenile salmon surviving migration to ocean:", 0.1),
    numericInput("sa", "Sa: Probability of survival of 2,3,4 yr old (adult) salmon returning to natal stream:", 0.15),    
    numericInput("h",  "h: Probability of returning 3 yr old salmon to its natal stream:", 1),
    numericInput("ff", "f: No of juveniles yielded from eggs produced by survived migration of adults:", 75),
    
    sliderInput("ngens", "Time scale:",
                min = 1, max=1000, value=500),
    
    HTML("<hr />"),
    radioButtons("scale_type","Choose scale of the plot:",
                 list("Semi log y" ="logar", 
                      "Linear"="linear")), 

    HTML("<hr />"),
    helpText(HTML("All source available on <a href = \"https://github.com/AnotherContext/lesliemat.git\" target=_blank>Github</a>")) 
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Description", h5("Short Intro"), "There are 7 species of Pacific salmon and several of these species have been recently listed as threatened under the US endangered species act. We will study coho salmon which has shown a decline in numbers and chinook salmon which have maintained there numbers. The model is composed of 4 age classes: (N1(t)) 1 year olds (juveniles) just before they migrate from the natal stream to the ocean, 2, 3 and 4 year olds in the ocean just before they migrate back to their natal stream to spawn (denoted by N2(t), N3(t), N4(t) respectively).",  
               br(), br(), h5("Model I"), 
               h6("Survival"), ("Individuals experience mortality primarily during migration from and to the natal stream. Mortality of individuals in the ocean is negligible. sj is the probability that juveniles survive migration to the ocean to become 2 year olds. Adults (2,3 and 4 year olds) have a probability sa of surviving the migration back to the stream to spawn."),
               h6("Reproduction"), ("Adults (2,3 and 4 year olds) return to their natal stream to spawn. With probability h they return when they are 3 years old. With probability (1 - h)/2 they return to spawn at either age 2 or 4. If the adults survive the migration to the stream they produce a number of eggs which eventually yield f 1-year old individuals. Fecundity is assumed to be independent of the age at which adults spawn."), 
               h6("Spawning age"), ("Coho salmon only spawn as 3 year olds (h = 1), while chinook salmon from the southern portion of their range spawn mainly in year 3 also (h = 0:8) and those from the northern portion of the range spawn in most years 1"), 
               h6("Other"), ("The fecundities are F2 = Sa* f* (1 - h)/2; F3 = Sa* f* 2h/(1 + h), F4 = Sa*f; and the survival probabilities are P1 = Sj , P2 = (1 + h)/2, P3 = (1 - h)/(1 + h). The default values of the parameters be Sj = 0.1, f = 75 and Sa = 0.15."), 
               br(), br(), h5("Model II"),
               ("2 year olds must compete with the larger salmon for food, so if the total population of adults is high the fecundity of the 2 years olds is reduced. This can be modelled by a Ricker function, giving the number of 1 year olds produced by 2 year old adults at a particular time by"), 
               br(), br(), ("Sa* f* (1 - h)/2*N2(t)*exp(-c(N2(t) + N3(t) + N4(t))"), 
               br(), br(), ("in which c determines the strength of the density dependence. This results in a modification to the matrix A, and entries will now depend on N2(t)+N3(t)+N4(t), giving a density-dependence Leslie matrix model. As density-dependent leslie matrix models can have richer dynamics than a density-independent model, we examine the behaviour of the model by looking at how the solution N(t) changes as a function time.")
      ), 
      tabPanel("Model", h5("Fig1. Projection of population size for age-structured populations"), plotOutput("plotmodel"), br(), h5("Tab2. Leslie matrix"), tableOutput("Lesliemat")),
      tabPanel("Stable Age Distribution", h5("Stable Age Vector"), verbatimTextOutput("stableage"), h5("Fig2. Stable Age Distribution plot"), plotOutput("plotstableage")),
      tabPanel("Elasticity", h5("Tab3. Sensitivity matrix"), tableOutput("elastTab"), h5("Final check with automatic output from eigen.analysis() fuction from popbio library"), verbatimTextOutput("elastfunc"))
    )
  )
))
