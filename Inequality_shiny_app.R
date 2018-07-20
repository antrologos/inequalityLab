rm(list=ls())
options(scipen = 999)
library(shiny)
library(shinydashboard)
library(tidyverse)
library(scales)
library(plotly)
library(Hmisc)
library(IC2)
library(data.table)

#setwd("")
#===========================================================================
# Global parameters

input = list()

input$quantile_range_add <- c(0, 1)
input$intercept_add      <- 0

input$quantile_range_mult <- c(0,1)
input$intercept_mult      <- 1

input$n_quantiles  = 100
input$poverty_line = 450

input$auto_sort = "yes"

data <- fread("data_pnad2015.csv") %>% 
        as_tibble() 

#===========================================================================
# Functions

calcShare <- function(x, wgt = NULL, quantile_lb, quantile_ub){
        
        na_x <- is.na(x)
        x <- x[!na_x] 
        
        if(is.null(wgt)){
                wgt <- rep(1, length(x))
        }else{
                wgt <- wgt[!na_x]
        }
        
        lower_bound_quantile = wtd.quantile(x, weights = wgt, probs = quantile_lb)
        upper_bound_quantile = wtd.quantile(x, weights = wgt, probs = quantile_ub)
        
        sum(x[x >= lower_bound_quantile & x <= upper_bound_quantile] * wgt[x >= lower_bound_quantile & x <= upper_bound_quantile])/sum(x * wgt)
}

quantilize <- function(x, n_quantiles = 100, wgt = NULL){
        na_x <- is.na(x)
        x <- x[!na_x] 
        
        if(is.null(wgt)){
                wgt <- rep(1, length(x))
        }else{
                wgt <- wgt[!na_x]
        }
        
        x_categories <- cut(x = x, 
                            breaks = wtd.quantile(x = x, weights = wgt,
                                                  probs = seq(0, 1, 1/n_quantiles)),include.lowest = T)
        
        tibble(x, x_categories, wgt) %>%
                group_by(x_categories) %>%
                summarise(quantile_mean_income = wtd.mean(x, wgt)) %>%
                ungroup() %>%
                select(-x_categories) %>%
                mutate(q = (1:n())/n_quantiles) 
        
}



make_modifier <- function(n_quantiles, 
                          quantile_range_mult, intercept_mult, gradual_mult,
                          quantile_range_add,  intercept_add, gradual_add){
        
        q = seq(1/n_quantiles, 1, 1/n_quantiles)
        
        multiplier = rep(1, length(q))
        adder      = rep(0, length(q))
        
        ls = quantile_range_add[[1]]
        us = quantile_range_add[[2]]
        
        lm = quantile_range_mult[[1]]
        um = quantile_range_mult[[2]]
        
        
        if(gradual_add == "descending"){
                adder[q >=  ls & q <=  us]      = intercept_add * ((q[q >= ls & q <=  us] - q[q == ls])/(q[q == us] - q[q == ls])) %>% sort(decreasing = T)
                
        } else if(gradual_add == "ascending"){
                adder[q >=  ls & q <=  us]      = intercept_add * ((q[q >= ls & q <=  us] - q[q == ls])/(q[q == us] - q[q == ls]))
        } else {
                adder[q >=  ls & q <=  us]      = intercept_add 
        }
        
        if(gradual_mult == "descending"){
                multiplier[q >=  lm & q <=  um] = 1 + (intercept_mult - 1)* ((q[q >= lm & q <=  um] - q[q == lm])/(q[q == um] - q[q == lm])) %>% sort(decreasing = T)
        } else if(gradual_mult == "ascending"){
                multiplier[q >=  lm & q <=  um] = 1 + (intercept_mult - 1)* ((q[q >= lm & q <=  um] - q[q == lm])/(q[q == um] - q[q == lm]))
        } else {
                multiplier[q >=  lm & q <=  um] = intercept_mult
        }
       
        quantiles_modifier <- tibble(q, adder, multiplier)
        
        quantiles_modifier
        
}


modify_quantiles <- function(n_quantiles , quantiles_modifier, auto_sort = "yes"){
        
        quantile_mean_income <- quantilize(data$income, n_quantiles = n_quantiles)$quantile_mean_income
        
        modified_quantile_mean_income = quantiles_modifier$adder + quantiles_modifier$multiplier * quantile_mean_income
        
        modified_quantiles <- tibble(modified_quantile_mean_income)
        
        
        if(auto_sort == "yes"){
                modified_quantiles <- modified_quantiles %>%
                        arrange(modified_quantile_mean_income) %>%
                        mutate(q = (1:n_quantiles)/n_quantiles)
        } else {
                modified_quantiles <- modified_quantiles %>%
                        mutate(q = quantiles_modifier$q)
        }
        
        modified_quantiles
} 

modify_incomes <- function(x , quantiles_modifier, wgt = NULL, auto_sort = "yes"){
        na_x <- is.na(x)
        x <- x[!na_x] 
        
        if(is.null(wgt)){
                wgt <- rep(1, length(x))
        }else{
                wgt <- wgt[!na_x]
        }
        
        n_quantiles <- length(unique(quantiles_modifier$q))
        
        x_categories <- cut(x = x, 
                            breaks = wtd.quantile(x = x, weights = wgt,
                                                  probs = seq(0, 1, 1/n_quantiles)),
                            include.lowest = T)
        
        x_categories <- as.numeric(x_categories)/n_quantiles %>% round(14)
        quantiles_modifier$q <- quantiles_modifier$q %>% round(14)
        
        
        modified_income_data <- tibble(income = x, 
                                       q      = x_categories, 
                                       wgt    = wgt) %>%
                left_join(y = quantiles_modifier, by = "q") %>% 
                mutate(modified_income = adder + income * multiplier) %>%
                select(-income, -q, -multiplier, -adder)
        
        if(auto_sort == "yes"){
                modified_income_data <- modified_income_data %>%
                        arrange(modified_income) %>%
                        mutate(person =  1:length(modified_income))
        } else {
                
                modified_income_data <- modified_income_data %>%
                        mutate(person =  1:length(modified_income))
        }
        
        
        modified_income_data
}



make_lorenz <- function(x, wgt = NULL){
        na_x <- is.na(x)
        x <- x[!na_x] 
        
        if(is.null(wgt)){
                wgt <- rep(1, length(x))
        }else{
                wgt <- wgt[!na_x]
        }
        
        
        lorenz_data <- tibble(pop_cum    = c(0, cumsum(wgt)/sum(wgt)),
                              income_cum = c(0, cumsum(x*wgt)/sum(x*wgt)))
        
        lorenz_data
}




#===========================================================================
# Shiny app



ui <- shinyUI(
        dashboardPage(title = "Income Inequality Lab",
                
                dashboardHeader(title = "Income Inequality Lab"),
                
                
                dashboardSidebar(
                        collapsed = TRUE,
                        
                        sidebarMenu(
                                menuItem("Play with an income distribution", tabName = "playTab", icon = icon("dashboard"), selected = TRUE),
                                menuItem("About me and this app", tabName = "about", icon = icon("dashboard"))
                        )
                        
                        
                        
                ),
                dashboardBody(
                        
                        tabItems(
                                tabItem(tabName = "playTab",
                                        #h2("Play"),
                                        
                                        fluidPage(
                                                
                                                fluidRow(
                                                        
                                                        column(3,
                                                               box(title = "Give some amount of money to a quantile bracket", width = 12,
                                                                   plotlyOutput('adder', width  = "90%", height = "180px"),
                                                                   sliderInput("quantile_range_add", "Quantile range:", min = 0, max = 1, value = c(0, 1), step = .01),
                                                                   radioButtons(inputId = "gradual_add", 
                                                                                label = "Make it gradual?", 
                                                                                choices = c("descending", "no", "ascending"), 
                                                                                selected = "no", choiceNames = c("Yes, ascending", "No", "Yes, descending"),
                                                                                inline = T),
                                                                   column(12,
                                                                          numericInput(inputId = 'intercept_add', label = 'Amount to give:', value = input$intercept_add)
                                                                   ),
                                                                   background = "navy"
                                                               ),
                                                               
                                                               box(collapsible = T, collapsed = T, title = "More options", width = 12,
                                                                       numericInput(inputId = 'n_quantiles',  label = 'Number of quantiles', value = input$n_quantiles),
                                                                       numericInput(inputId = 'poverty_line', label = 'Poverty Line', value = input$poverty_line),
                                                                       radioButtons(inputId = "auto_sort", "Automatically resort incomes?", choices = c("yes", "no"), selected = "yes", inline=T)
                                                               )
                                                        ),
                                                        column(6,
                                                               tabBox(width = 12,
                                                                      height = "35vw",
                                                                      
                                                                      tabPanel(title = "Plots",
                                                                               fluidRow(
                                                                                       box(
                                                                                               plotlyOutput('pen_parade', 
                                                                                                            width = "100%",
                                                                                                            height = "13vw"),
                                                                                               title = "Pen's Parade (individuals)",
                                                                                               solidHeader = TRUE, status = "primary"
                                                                                       ),
                                                                                       box(
                                                                                               plotlyOutput('quantiles_parade', 
                                                                                                            width = "100%",
                                                                                                            height = "13vw"),
                                                                                               title = "Pen's Parade (quantiles)",
                                                                                               solidHeader = TRUE, status = "primary"
                                                                                       )
                                                                               ), 
                                                                               fluidRow(
                                                                                       box(
                                                                                               plotlyOutput('lorenz', 
                                                                                                            width = "100%",
                                                                                                            height = "13vw"),
                                                                                               title = "Lorenz's Curve",
                                                                                               solidHeader = TRUE, status = "primary"
                                                                                       ),
                                                                                       
                                                                                       box(
                                                                                               plotlyOutput('density', 
                                                                                                            width = "100%",
                                                                                                            height = "13vw"),
                                                                                               title = "Probability density",
                                                                                               solidHeader = TRUE, status = "primary"
                                                                                       )                                                        
                                                                               )
                                                                      ),
                                                                      
                                                                      tabPanel(title = "Measures",
                                                                               tableOutput("inequalityTable")
                                                                      ),
                                                                      
                                                                      id="tabBox"
                                                               ) #tabbox
                                                        ), #column
                                                        column(3,
                                                               box(title = "Multiply the income of a quantile bracket",  width = 12,
                                                                   plotlyOutput('multiplier', width  = "90%", height = "180px"),
                                                                   
                                                                   sliderInput("quantile_range_mult", "Quantile Range:", min = 0, max = 1, value = c(0,1), step = .01),
                                                                   radioButtons(inputId = "gradual_mult", 
                                                                                label = "Make it gradual?", 
                                                                                choices = c("descending", "no", "ascending"), 
                                                                                selected = "no", choiceNames = c("Yes, ascending", "No", "Yes, descending"),
                                                                                inline = T),
                                                                   column(12,
                                                                          numericInput(inputId = 'intercept_mult', label = 'Amount to multiply by:', value = input$intercept_mult)
                                                                   ),
                                                                   background = "navy"
                                                               ) # box
                                                        ) #column
                                                )
                                        ) #fluidpage
                                        
                                ),
                                
                                tabItem(tabName = "about",
                                        box(title = "About me", width = 4,
                                            p("My name is", strong(a(href ="https://scholar.google.com.br/citations?user=_GGcZ8PnOA0C",
                                                                      "Rogerio Barbosa")),
                                               ". I' have got my PhD in Sociolgy at the University of Sao Paulo (USP, Brazil). My research interests are income inequality, education, Social Sciences Methology, and Econometrics."),
                                            p(),
                                            p("Found an error? Have a suggestion? Want to send me a message? Easy:"),
                                            p("antrologos [at] usp [dot] br")
                                            ),
                                        
                                        box(title = "About this app", width = 4,
                                            p("This web application was built for educational purposes. It uses a random sample (20K individuals) of household income per capita, drawn from the National Household Sample Survey (PNAD-IBGE, Brazil). Every calculation -- including the plots and inequality indexes -- makes use of sample weights. Once this is a sub-sample, the original weights had to be multiplied by a constant in order to make the number of cases add up to the estimated population in that year. Income values are expressed are local currency ('reais', R$) They are nominal values -- no deflation procedure was applied. A little uniform random noise (ranging from -R$ 0.1 to R$ 0.1) was also added to all incomes in order to avoid repeated values in the quantiles' cut (which could cause error in some procedures here). Once the original incomes were all integer values, this does not cause a reordering of cases but forces people with equal incomes to be ordered. In order to not produce people with negative values, those with zero income received a random uniform noise ranging from R$ 0.00 to R$ 0.00. It is important to stress that this strategy of adding noises must not be applied in a actual analysis of income data.")
                                            )
                                )# tabItem
                                
                        )
                                
                        ) #dashboardbody
        )#dashboardpage
)#shinyUI




server <- function(input, output, session) {
        
        values <- reactiveValues()
        
        observeEvent({
                input$n_quantiles
                
                input$quantile_range_add
                input$intercept_add
                input$gradual_add
                
                input$quantile_range_mult
                input$intercept_mult
                input$gradual_mult
                
                input$auto_sort
        }
        , {
                
                values$quantile_range_add <- input$quantile_range_add
                values$quantile_range_mult <- input$quantile_range_mult
                
                if(values$quantile_range_add[1] == 0)  values$quantile_range_add[1]  <- 1/input$n_quantiles
                if(values$quantile_range_mult[1] == 0) values$quantile_range_mult[1] <- 1/input$n_quantiles
                
                values$quantiles_modifier <- make_modifier(n_quantiles        = input$n_quantiles, 
                                                           
                                                           quantile_range_add = input$quantile_range_add,
                                                           intercept_add      = input$intercept_add,
                                                           gradual_add        = input$gradual_add,
                                                           
                                                           quantile_range_mult = input$quantile_range_mult,
                                                           intercept_mult      = input$intercept_mult,
                                                           gradual_mult        = input$gradual_mult
                                                           )
                
                values$modified_quantiles <- modify_quantiles(n_quantiles          = input$n_quantiles,
                                                              quantiles_modifier   = values$quantiles_modifier, 
                                                              auto_sort            = input$auto_sort)
                
                
                values$quantilized_original_income <- quantilize(data$income, 
                                                                 wgt = data$wgt, 
                                                                 n_quantiles = input$n_quantiles)
                
                values$modified_income_data <- modify_incomes(x   = data$income,
                                                              wgt = data$wgt,
                                                              quantiles_modifier = values$quantiles_modifier,
                                                              auto_sort = input$auto_sort)
        })
        
        observe({
                new_step <- 1/input$n_quantiles
                updateSliderInput(session, inputId = "quantile_range_add", step = new_step, 
                                  label = NULL, min = NULL, max = NULL, value = NULL)
                updateSliderInput(session, inputId = "quantile_range_mult", step = new_step, 
                                  label = NULL, min = NULL, max = NULL, value = NULL)
        })
        
        output$multiplier <- renderPlotly({
                p <- values$quantiles_modifier %>%
                        ggplot() + 
                        geom_abline(slope = 0, intercept = 0, color = "black", lwd = 1.5) +
                        geom_line(aes(x = q, y = multiplier), color = "red") +
                        ylim(-max(abs(values$quantiles_modifier$multiplier)), max(abs(values$quantiles_modifier$multiplier))) +
                        xlab("Quantile") +
                        theme(axis.text.x = element_text(size=6),
                              axis.text.y = element_text(size=6))
                
                p <- ggplotly(p) %>% layout(margin = list(l = 75))
                p$elementId <- NULL
                p
        })
        
        
        output$adder <- renderPlotly({
                p <- values$quantiles_modifier %>%
                        ggplot() + 
                        geom_abline(slope = 0, intercept = 0, color = "black", lwd = 1.5) +
                        geom_line(aes(x = q, y = adder), color = "red") +
                        ylim(-max(abs(values$quantiles_modifier$adder)), max(abs(values$quantiles_modifier$adder))) +
                        xlab("Quantile") +
                        theme(axis.text.x = element_text(size=6),
                              axis.text.y = element_text(size=6))
                
                p <- ggplotly(p) %>% layout(margin = list(l = 75))
                p$elementId <- NULL
                p
        })
        
        
        output$pen_parade <- renderPlotly({
                                p <- ggplot() + 
                        geom_line(data = data, aes(x = cumsum(wgt),
                                                   y = income),  
                                  lwd = 1.3, color = "dark blue") +
                        geom_line(data = values$modified_income_data, aes(x = cumsum(wgt), 
                                                                          y = modified_income),
                                                                          lwd = 1.3, lty = 2, color = "red") +
                        geom_abline(data = NULL, aes(intercept = input$poverty_line, slope = 0), 
                                    lwd = .5, lty = 2, color = "black") +
                        xlab("Person") + ylab("Income") +
                        theme(axis.text.x = element_text(size=6),
                              axis.text.y = element_text(size=6)) +
                        scale_y_continuous(label=dollar_format(prefix = "R$")) +
                        scale_x_continuous(labels=function(x) paste0(format(trunc(x/1000000), scientific = FALSE), "mi"))
                
                p <- ggplotly(p, tooltip="label") %>% layout(margin = list(l = 75))
                p$elementId <- NULL
                p
                
        })
        
        output$quantiles_parade <- renderPlotly({
                p <- ggplot() + 
                        geom_line(data = values$quantilized_original_income, 
                                  aes(x = q, 
                                      y = quantile_mean_income),
                                  lwd = 1.3, color = "dark blue") +
                        geom_line(data = values$modified_quantiles,
                                  aes(x = q, 
                                      y = modified_quantile_mean_income), 
                                  lwd = 1.3, lty = 2, color = "red") +
                        geom_abline(data = NULL, aes(intercept = input$poverty_line, slope = 0), 
                                    lwd = .5, lty = 2, color = "black") +
                        xlab("Quantile") + ylab("Mean income") +
                        theme(axis.text.x = element_text(size=6),
                              axis.text.y = element_text(size=6)) +
                        scale_y_continuous(label=dollar_format(prefix = "R$"))
                
                p <- ggplotly(p) %>% layout(margin = list(l = 75))
                p$elementId <- NULL
                p
        })
        
        
        output$lorenz <- renderPlotly({
                p <- ggplot() +
                        geom_line(data = make_lorenz(data$income, 
                                                     wgt = data$wgt), 
                                  aes(x = pop_cum, 
                                      y = income_cum),
                                  lwd = 1.3, color = "dark blue") + 
                        geom_line(data = make_lorenz(values$modified_income_data$modified_income,
                                                     wgt = values$modified_income_data$wgt),
                                  aes(x = pop_cum, 
                                      y = income_cum),
                                  lwd = 1.3, lty = 2, color = "red") +
                        geom_abline(slope = 1, intercept = 0) +
                        xlim(0,1) + ylim(0,1) +
                        xlab("Population (Cumulative)") + ylab("Income (Cumulative)") +
                        theme(axis.text.x = element_text(size=6),
                              axis.text.y = element_text(size=6))
                
                p <- ggplotly(p) %>% layout(margin = list(l = 75))
                p$elementId <- NULL
                p
        })
        
        
        output$density <- renderPlotly({
                p <- ggplot() + 
                        geom_density(data = data, 
                                     aes(x = income, weight = wgt/sum(wgt)),
                                     fill = "dark blue", alpha = .3) +
                        geom_density(data = values$modified_income_data, 
                                     aes(x = modified_income, weight = wgt/sum(wgt)), 
                                     fill = "red", alpha = .3) +
                        xlab("Income") + ylab("Probability Density") +
                        theme(axis.text.x = element_text(size=6),
                              axis.text.y = element_text(size=6)) +
                        scale_x_continuous(label=dollar_format(prefix = "R$"))
                
                p <- ggplotly(p) %>% layout(margin = list(l = 75))
                p$elementId <- NULL
                p
        })
        
        
        
        inequalityMeasures <- reactive({
                
                tibble(Measure = c("Mean",
                                   "Median",
                                   "Minimum",
                                   "Maximum",
                                   "Gini",
                                   "Theil-T",
                                   "Theil-L",
                                   "Atkinson (epsilon = 0.5)",
                                   "Atkinson (epsilon = 1)",
                                   "Atkinson (epsilon = 2)",
                                   "Var-log",
                                   "Top 10% Share",
                                   "Top 1% Share",
                                   "Top 0.1% Share",
                                   "Bottom 50% Share",
                                   "Poverty rate"),
                       
                       `Original Measures` = with(data, {
                               c(
                                               wtd.mean(income, weights = wgt),
                                               wtd.quantile(income, probs = .5, weights = wgt),  
                                               min(income, na.rm = T),
                                               max(income, na.rm = T),
                                               calcSGini(income, w = wgt)$ineq$index,
                                               calcGEI(income, w = wgt, alpha = 1)$ineq$index,
                                               calcGEI(income, w = wgt, alpha = 0)$ineq$index,
                                               calcAtkinson(income, w = wgt, epsilon = .5)$ineq$index,
                                               calcAtkinson(income, w = wgt, epsilon = 1)$ineq$index,
                                               calcAtkinson(income, w = wgt, epsilon = 2)$ineq$index,
                                               wtd.var(log(income), weights = wgt),
                                               calcShare(income, wgt = wgt, quantile_lb = .90, quantile_ub = 1),
                                               calcShare(income, wgt = wgt, quantile_lb = .99, quantile_ub = 1),
                                               calcShare(income, wgt = wgt, quantile_lb = .999, quantile_ub = 1),
                                               calcShare(income, wgt = wgt, quantile_lb = 0, quantile_ub = .5),
                                               wtd.mean(income <= input$poverty_line, weights = wgt)
                               )
                       }),
                       
                       `Counterfactual Measures` = with(values$modified_income_data, {
                               c(
                                       wtd.mean(modified_income, weights = wgt),
                                       wtd.quantile(modified_income, probs = .5, weights = wgt),  
                                       min(modified_income, na.rm = T),
                                       max(modified_income, na.rm = T),
                                       calcSGini(modified_income, w = wgt)$ineq$index,
                                       calcGEI(modified_income, w = wgt, alpha = 1)$ineq$index,
                                       calcGEI(modified_income, w = wgt, alpha = 0)$ineq$index,
                                       calcAtkinson(modified_income, w = wgt, epsilon = .5)$ineq$index,
                                       calcAtkinson(modified_income, w = wgt, epsilon = 1)$ineq$index,
                                       calcAtkinson(modified_income, w = wgt, epsilon = 2)$ineq$index,
                                       wtd.var(log(modified_income), weights = wgt),
                                       calcShare(modified_income, wgt = wgt, quantile_lb = .90, quantile_ub = 1),
                                       calcShare(modified_income, wgt = wgt, quantile_lb = .99, quantile_ub = 1),
                                       calcShare(modified_income, wgt = wgt, quantile_lb = .999, quantile_ub = 1),
                                       calcShare(modified_income, wgt = wgt, quantile_lb = 0, quantile_ub = .5),
                                       wtd.mean(modified_income <= input$poverty_line, weights = wgt)
                                       )
                               })
                       ) %>%
                        mutate(`Change (%)` = ((`Counterfactual Measures` - `Original Measures`)/`Original Measures`)*100)
                
                
        })
        
        
        output$inequalityTable <- renderTable({
                inequalityMeasures()   
        },digits = 3)
        
        
        
        
        
}


#===========================================================================
# Shiny app execution

shinyApp(ui = ui, server = server) 









