library(shiny)
library(reshape2)
library(ggplot2)
library(plyr)

theme_set(theme_bw(base_size=20))
update_geom_defaults("line", list(size = 1.75))

dichroics_df <- read.csv('dichroics.csv')
dichroics <- melt(dichroics_df, id.vars='wavelength')

spectrum <- c('#0046ff', '#00c0ff', '#00ff92', '#4aff00', '#a3ff00', '#f0ff00', '#ffbe00', '#ff6300', '#ff0000')

# Longpass dichroics
# D1: 550
# D2: 600
# D3: 625
# D4: 500
# D5: 475
# D6: 525
# D7: 575

# Camera tree paths
c1 <- '(RD1) * (RD4) * (RD5)'
c2 <- '(RD1) * (RD4) * (TD5)'
c3 <- '(RD1) * (TD4) * (RD6)'
c4 <- '(RD1) * (TD4) * (TD6)'
c5 <- '(TD1) * (RD2) * (RD7)'
c6 <- '(TD1) * (RD2) * (TD7)'
c7 <- '(TD1) * (TD2) * (RD3)'
c8 <- '(TD1) * (TD2) * (TD3)'

c1e <- gsub('RD', '1 - dichroics_df$D', gsub('TD', 'dichroics_df$D', c1))
c2e <- gsub('RD', '1 - dichroics_df$D', gsub('TD', 'dichroics_df$D', c2))
c3e <- gsub('RD', '1 - dichroics_df$D', gsub('TD', 'dichroics_df$D', c3))
c4e <- gsub('RD', '1 - dichroics_df$D', gsub('TD', 'dichroics_df$D', c4))
c5e <- gsub('RD', '1 - dichroics_df$D', gsub('TD', 'dichroics_df$D', c5))
c6e <- gsub('RD', '1 - dichroics_df$D', gsub('TD', 'dichroics_df$D', c6))
c7e <- gsub('RD', '1 - dichroics_df$D', gsub('TD', 'dichroics_df$D', c7))
c8e <- gsub('RD', '1 - dichroics_df$D', gsub('TD', 'dichroics_df$D', c8))

c1d <- eval(parse(text=c1e))
c2d <- eval(parse(text=c2e))
c3d <- eval(parse(text=c3e))
c4d <- eval(parse(text=c4e))
c5d <- eval(parse(text=c5e))
c6d <- eval(parse(text=c6e))
c7d <- eval(parse(text=c7e))
c8d <- eval(parse(text=c8e))

channels_df <- data.frame(cbind(c1d, c2d, c3d, c4d, c5d, c6d, c7d, c8d))
channels_df$wavelength <- dichroics_df$wavelength
channels <- melt(channels_df, id.vars='wavelength')


# # Fake fluorophore data
lambda <- dichroics_df$wavelength
fcs <- seq(470, by=25, length.out=8)
f <- sapply(fcs, function(x) exp(-(lambda - x)^2 / (2 * 25^2)))
f <- data.frame(f)
f$wavelength <- lambda
fluorophores <- melt(f, id.vars='wavelength')

# FP data
fl <- read.csv('fluorophores_fps_1.csv')
fl[is.na(fl)] <- 0
fl <- fl[, c(2, 3, 4, 5, 6, 7, 8, 9, 1)]
fl <- fl[fl$wavelength >= 450 & fl$wavelength <= 750, ]
f <- fl
fluorophores <- melt(fl, id.vars='wavelength')

# Mixing matrix calcs


ui <- fluidPage(
    titlePanel("8 camera multispectral system"),

    sidebarLayout(
        sidebarPanel(
            sliderInput("f1b", "Fluorophore 1 brightness", min = 0, max = 5, value = 1),
            sliderInput("f2b", "Fluorophore 2 brightness", min = 0, max = 5, value = 1),
            sliderInput("f3b", "Fluorophore 3 brightness", min = 0, max = 5, value = 1),
            sliderInput("f4b", "Fluorophore 4 brightness", min = 0, max = 5, value = 1),
            sliderInput("f5b", "Fluorophore 5 brightness", min = 0, max = 5, value = 1),
            sliderInput("f6b", "Fluorophore 6 brightness", min = 0, max = 5, value = 1),
            sliderInput("f7b", "Fluorophore 7 brightness", min = 0, max = 5, value = 1),
            sliderInput("f8b", "Fluorophore 8 brightness", min = 0, max = 5, value = 1),
        ),

        mainPanel(
            h3("Dichroic spectra"),
            plotOutput("dichroicPlot", height='200px'),
            h3("Channel spectra"),
            plotOutput("channelPlot", height='200px'),
            h3("Fluorophore spectra"),
            plotOutput("fluorophorePlot", height='200px'),
            h3("Camera spectra"),
            plotOutput("cameraPlot", height='200px'),
            h3("Mixing matrix"),
            textOutput("conditionText"),
            plotOutput("mixingPlot", height='400px')
        )
    )
)



server <- function(input, output) {
    output$dichroicPlot <- renderPlot({
        # bins <- input$bins
        ggplot(dichroics, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Transmission') + theme(legend.position='none') + scale_color_brewer(palette='Set2')
    })
    
    output$channelPlot <- renderPlot({
        ggplot(channels, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Transmission') + theme(legend.position='none') + scale_color_manual(values=spectrum)
    })
    

    output$fluorophorePlot <- renderPlot({
        f[, 1] <- f[, 1] * input$f1b
        f[, 2] <- f[, 2] * input$f2b
        f[, 3] <- f[, 3] * input$f3b
        f[, 4] <- f[, 4] * input$f4b
        f[, 5] <- f[, 5] * input$f5b
        f[, 6] <- f[, 6] * input$f6b
        f[, 7] <- f[, 7] * input$f7b
        f[, 8] <- f[, 8] * input$f8b
        fluorophores <- melt(f, id.vars='wavelength')
        
        ggplot(fluorophores, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Fluorescence') + theme(legend.position='none') + scale_color_brewer(palette='Set2')
    })
    
    output$cameraPlot <- renderPlot({
        f[, 1] <- f[, 1] * input$f1b
        f[, 2] <- f[, 2] * input$f2b
        f[, 3] <- f[, 3] * input$f3b
        f[, 4] <- f[, 4] * input$f4b
        f[, 5] <- f[, 5] * input$f5b
        f[, 6] <- f[, 6] * input$f6b
        f[, 7] <- f[, 7] * input$f7b
        f[, 8] <- f[, 8] * input$f8b
        fluorophores <- melt(f, id.vars='wavelength')
        
        emission <- ddply(fluorophores, .(wavelength), summarise, value=sum(value))
        cs <- channels_df
        cs[, 1:8] <- cs[, 1:8] * emission$value
        cs <- melt(cs, id.vars='wavelength')
        
        ggplot(cs, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Signal') + theme(legend.position='none') + scale_color_manual(values=spectrum)
    })
    
    output$mixingPlot <- renderPlot({
        cmat <- t(as.matrix(channels_df[, 1:8]))
        cmat[1, ] <- cmat[1, ] * input$f1b
        cmat[2, ] <- cmat[2, ] * input$f2b
        cmat[3, ] <- cmat[3, ] * input$f3b
        cmat[4, ] <- cmat[4, ] * input$f4b
        cmat[5, ] <- cmat[5, ] * input$f5b
        cmat[6, ] <- cmat[6, ] * input$f6b
        cmat[7, ] <- cmat[7, ] * input$f7b
        cmat[8, ] <- cmat[8, ] * input$f8b
        
        fmat <- as.matrix(f[, 1:8])
        mmat <- cmat %*% fmat
        mmat <- mmat / max(mmat)
        
        ggplot(melt(mmat), aes(x=Var2, y=Var1, fill=value)) + geom_tile() + labs(x='Fluorophore', y='Channel') + theme(axis.text.x = element_text(angle=90))
    })
    
    output$conditionText <- renderText({
        cmat <- t(as.matrix(channels_df[, 1:8]))
        cmat[1, ] <- cmat[1, ] * input$f1b
        cmat[2, ] <- cmat[2, ] * input$f2b
        cmat[3, ] <- cmat[3, ] * input$f3b
        cmat[4, ] <- cmat[4, ] * input$f4b
        cmat[5, ] <- cmat[5, ] * input$f5b
        cmat[6, ] <- cmat[6, ] * input$f6b
        cmat[7, ] <- cmat[7, ] * input$f7b
        cmat[8, ] <- cmat[8, ] * input$f8b
        
        fmat <- as.matrix(f[, 1:8])
        mmat <- cmat %*% fmat
        mmat <- mmat / max(mmat)
        
        paste0("Condition number = ", round(kappa(mmat), 0))
    })
}



shinyApp(ui = ui, server = server)
