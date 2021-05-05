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
fl <- read.csv('live-cell_dye_spectra.csv', check.names = F)
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
            selectInput("f1n", "Fluorophore 1", choices=head(colnames(f), -1), selected=colnames(f)[1]),
            selectInput("f2n", "Fluorophore 2", choices=head(colnames(f), -1), selected=colnames(f)[2]),
            selectInput("f3n", "Fluorophore 3", choices=head(colnames(f), -1), selected=colnames(f)[3]),
            selectInput("f4n", "Fluorophore 4", choices=head(colnames(f), -1), selected=colnames(f)[4]),
            selectInput("f5n", "Fluorophore 5", choices=head(colnames(f), -1), selected=colnames(f)[5]),
            selectInput("f6n", "Fluorophore 6", choices=head(colnames(f), -1), selected=colnames(f)[6]),
            selectInput("f7n", "Fluorophore 7", choices=head(colnames(f), -1), selected=colnames(f)[7]),
            selectInput("f8n", "Fluorophore 8", choices=head(colnames(f), -1), selected=colnames(f)[8]),
        ),

        mainPanel(
            h3("Dichroic spectra"),
            plotOutput("dichroicPlot", height='200px'),
            h3("Channel spectra"),
            plotOutput("channelPlot", height='200px'),
            h3("Fluorophore spectra"),
            plotOutput("fluorophorePlot", height='300px'),
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
        ggplot(fluorophores, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Fluorescence', col='') + theme(legend.position='top') + scale_color_brewer(palette='Set2') + theme(legend.text=element_text(size=10))
    })
    
    output$cameraPlot <- renderPlot({
        emission <- ddply(fluorophores, .(wavelength), summarise, value=sum(value))
        cs <- channels_df
        cs[, 1:8] <- cs[, 1:8] * emission$value
        cs <- melt(cs, id.vars='wavelength')
        
        ggplot(cs, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Signal') + theme(legend.position='none') + scale_color_manual(values=spectrum)
    })
    
    output$mixingPlot <- renderPlot({
        cmat <- t(as.matrix(channels_df[, 1:8]))
        fmat <- as.matrix(f[, 1:8])
        mmat <- cmat %*% fmat
        mmat <- mmat / max(mmat)
        
        ggplot(melt(mmat), aes(x=Var2, y=Var1, fill=value)) + geom_tile() + labs(x='Fluorophore', y='Channel') + theme(axis.text.x = element_text(angle=90))
    })
    
    output$conditionText <- renderText({
        cmat <- t(as.matrix(channels_df[, 1:8]))
        fmat <- as.matrix(f[, 1:8])
        fmat <- scale(fmat, center=F, scale=colSums(fmat))
        mmat <- cmat %*% fmat
        mmat <- mmat / max(mmat)
        
        paste0("Condition number = ", round(kappa(mmat), 0))
    })
}



shinyApp(ui = ui, server = server)
