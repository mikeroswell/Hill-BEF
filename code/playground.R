#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
require(shiny)
require(magrittr, quietly = T)

# diversity functions (in q)
Ln = function(x) {ifelse(x != 0, log(x), 0)}
Exp = function(x, y) {ifelse(x != 0, x^y, 0)}
get.Dq = function(vec, q) {
    if(sum(vec) == 0) {return(0)}
    ifelse(q != 1,
           sum(Exp(vec/sum(vec), q))^(1/(1 - q)),
           exp(-sum((vec/sum(vec))*Ln(vec/sum(vec))))) }






# Define UI 
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            
            # stationary random input
            numericInput("sites", label = "sites", 20, width='75%'),
            numericInput("s", label = "species number:", 200, width='75%'),
            numericInput("abund_mu", label = "Abundance - mean:", 150, width='75%'),
            numericInput("abund_size", label = "Abundance - dispersion:", 1.5, width='75%'),
            numericInput("fun_beta", label = "Function - slope:", -1.1, width='75%'),
            numericInput("fun_shape", label = "Function - error:", 50, width='75%')
        ), # end sidebar panel
        
        mainPanel(fluidRow(
            verticalLayout(
                # cellHeights = c("200px", "200px"),
                splitLayout(cellWidths = c("50%", "50%"),
                            plotOutput(outputId = "plotgraph1"),
                            plotOutput(outputId = "plotgraph2")
                ), # end upper splitLayout
                splitLayout(cellWidths = c("50%", "50%"),
                            plotOutput(outputId = "plotgraph3"),
                            plotOutput(outputId = "plotgraph4")
                ) # end lower splitLayout
                
            ) # end vertical layout
        ))# end mainPanel
        
    ) # end sidebar layout
) # end fluid page

# Define server logic
server <- function(input, output) {
    
# old abundance: for each site, choose abundances across species using a fixed  
    a <- reactive({
        site.abund_mu = rgamma(input$sites, shape = input$abund_mu)
        site.abund_size = rgamma(input$sites, shape = input$abund_size)
        
        a <- sapply(1:as.numeric(input$s), function(i) {
            sapply(1:as.numeric(input$sites), function(k) {
                rnbinom(1,  
                mu = site.abund_mu[k], 
                size = site.abund_size[k])})
        }) %>% t
        a[which.max(rowSums(a)), which(colSums(a) == 0)] <- 1
        a
    })
    
    z <- reactive({sapply(1:as.numeric(input$sites), function(i) {
        y_true = exp( input$fun_beta*Ln(a()[,i])) # exp link glm
        shape = input$fun_shape # larger shape values = lower error
        rgamma(input$s, rate = shape/y_true, shape = shape)
    }) # close sapply
        }) # close reactive
    
    q <- seq(-2, 4, by = 0.2)
    # derived function and diversity
    f <- reactive({colSums(a()*z())})
    Dq <- reactive({sapply(q, function(q) apply(a(), 2, get.Dq, q))}) 
    ell <- 1-q
    # sites in rows, q' in columns
    # x.ax.val = seq(min(Dq), max(Dq), by = 1)
    
    output$plotgraph1 = renderPlot({
        plot(a(), z(), # col=alpha("gray", 0.2),
             main = "site selection function")
        for(i in 1:as.numeric(input$sites)) {
            X = a()[,i]
            mod = glm(z()[,i] ~ X, family=Gamma(link="log"))
            x = seq(1, max(a()), by=1)
            y = predict(mod, newdata=data.frame(X=x), type="response")
            points(x, y, type="l")
            } 
        })
    
    output$plotgraph2 = renderPlot({
        plot(range(ell), range(Dq()), type="n",
             xlab = "ell", ylab = "D",
             main = "site diversity profiles")
        apply(Dq(), 1, function(x) points(ell, x, 
                                          type="l", lwd=2))
        abline(v=0, lty=2)
        evenness = log(Dq()[,q==2])/log(Dq()[,q==0])
        mtext(paste("mean even =", round(mean(evenness), 2)), side=3, line=-2)
        mtext(paste("sd even =", round(sd(evenness), 2)), side=3, line=-3)
        
    })
    
    output$plotgraph3 = renderPlot({
        plot(f() ~ Dq()[,1], ylim = range(f()), xlim = range(Dq()), 
             type = "n", xlab = "D", ylab = "EF"
             , main = "D vs EF")
        Col=rainbow(ncol(Dq()))
        sapply(1:length(ell), function(x) points(Dq()[, x], f(), col=Col[x]))
        ### linear mod
        # plot(f ~ Dq[,1], ylim = range(f), xlim=range(Dq), type = "n")
        sapply(1:length(ell), function(x) abline(lm(f() ~ Dq()[, x]), col=Col[x], lwd=2))
    }) 
    
    output$plotgraph4 = renderPlot({
        plot(ell, apply(Dq(), 2, function(x) {
            cor(f(), x)
        }), ylab = "corelation coefficient", col=rainbow(ncol(Dq())), pch=19, cex=2,
        main = "ell-BEF raw cor", ylim = c(-1,1))
        abline(v=0, lty=2)

        
    })
} # close server

# Run the application 
shinyApp(ui = ui, server = server)

