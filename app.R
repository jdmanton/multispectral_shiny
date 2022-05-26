library(shiny)
library(reshape2)
library(ggplot2)
library(plyr)

theme_set(theme_bw(base_size=20))
update_geom_defaults("line", list(size = 1.75))

dichroics_df <- read.csv('dichroics.csv')
colnames(dichroics_df) <- c('wavelength', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7')
dichroics_df[is.na(dichroics_df)] <- 0
dichroics_df <- dichroics_df[dichroics_df$wavelength <= 750 & dichroics_df$wavelength >= 450, ]
dichroics <- melt(dichroics_df, id.vars='wavelength')

spectrum <- c('#0046ff', '#00c0ff', '#00ff92', '#4aff00', '#a3ff00', '#f0ff00', '#ffbe00', '#ff6300', '#ff0000')


# Longpass dichroics
# D1: 575
# D2: 625
# D3: 650
# D4: 525
# D5: 500
# D6: 550
# D7: 600

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


# Load fpbase.org data
fpbase_data <- readRDS('fpbase_data.rds')
fpbase_data <- fpbase_data[fpbase_data$Wavelength >= 450 & fpbase_data$Wavelength <= 750, ]

# Add blank spectrum
w <- unique(fpbase_data$Wavelength)
v <- 0 * w
n <- rep("**None**", length(w))
blank_spectrum <- data.frame(Wavelength=w, Emission=v, FP=n)
fpbase_data <- rbind(fpbase_data, blank_spectrum)

fp_names <- unique(fpbase_data$FP)


ui <- fluidPage(
    titlePanel("8 camera multispectral system"),
    
    sidebarLayout(
        sidebarPanel(
            selectInput("f1n", "Fluorophore 1", choices=fp_names, selected='EGFP'),
            selectInput("f2n", "Fluorophore 2", choices=fp_names, selected='mCherry'),
            selectInput("f3n", "Fluorophore 3", choices=fp_names, selected='**None**'),
            selectInput("f4n", "Fluorophore 4", choices=fp_names, selected='**None**'),
            selectInput("f5n", "Fluorophore 5", choices=fp_names, selected='**None**'),
            selectInput("f6n", "Fluorophore 6", choices=fp_names, selected='**None**'),
            selectInput("f7n", "Fluorophore 7", choices=fp_names, selected='**None**'),
            selectInput("f8n", "Fluorophore 8", choices=fp_names, selected='**None**'),
            div(style = "height:50px"),
            downloadButton("downloadData", "Download mixing matrix")
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
    
    # Reactive function for checking selected fluorophores
    fluorophores <- reactive({
        fluorophores <- fpbase_data[fpbase_data$FP == input$f1n, ]
        if (input$f2n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$FP == input$f2n, ])
        if (input$f3n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$FP == input$f3n, ])
        if (input$f4n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$FP == input$f4n, ])
        if (input$f5n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$FP == input$f5n, ])
        if (input$f6n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$FP == input$f6n, ])
        if (input$f7n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$FP == input$f7n, ])
        if (input$f8n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$FP == input$f8n, ])
        fluorophores
    })
    
    
    # Reactive function for calculating mixing matrix
    mixing_matrix <- reactive({
        f <- dcast(fluorophores(), Wavelength ~ FP, value.var='Emission')
        f[is.na(f)] <- 0
        
        cmat <- t(as.matrix(channels_df[, 1:8]))
        
        if (nrow(f) < ncol(cmat)) {
            # We've dropped some wavelengths, so put them back with zero
            if (min(f$Wavelength) != 450) {
                low_wavelengths <- 450:(min(f$Wavelength) - 1)
                empty_data <- matrix(0, ncol=ncol(f), nrow=length(low_wavelengths))
                empty_data[, 1] <- low_wavelengths
                colnames(empty_data) <- colnames(f)
                f <- rbind(as.data.frame(empty_data), as.matrix(f))
            }
            if (max(f$Wavelength) != 750) {
                high_wavelengths <- (max(f$Wavelength) + 1):750
                empty_data <- matrix(0, ncol=ncol(f), nrow=length(high_wavelengths))
                empty_data[, 1] <- high_wavelengths
                f <- rbind(as.matrix(f), empty_data)
            }
        }
        
        fmat <- as.matrix(f[, 2:ncol(f)])
        mmat <- cmat %*% fmat
        mmat <- mmat / sum(mmat)
        mmat
    })
    
    
    # Display of dichroic spectra
    output$dichroicPlot <- renderPlot({
        ggplot(dichroics, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Transmission') + theme(legend.position='none') + scale_color_brewer(palette='Set2')
    })
    
    
    # Display of camera channel spectra
    output$channelPlot <- renderPlot({
        ggplot(channels, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Transmission') + theme(legend.position='none') + scale_color_manual(values=spectrum)
    })
    
    
    # Display of fluorophore spectra
    output$fluorophorePlot <- renderPlot({
        ggplot(fluorophores(), aes(x=Wavelength, y=Emission, col=FP)) + geom_line() + labs(x='Wavelength / nm', y='Fluorescence', col='') + theme(legend.position='top') + scale_color_brewer(palette='Set2') + theme(legend.text=element_text(size=10))
    })
    
    
    # Display of camera-frame fluorophore spectra
    output$cameraPlot <- renderPlot({
        emission <- ddply(fluorophores(), .(Wavelength), summarise, value=sum(Emission))
        cs <- channels_df
        cs[, 1:8] <- cs[, 1:8] * emission$value
        cs <- melt(cs, id.vars='wavelength')
        
        ggplot(cs, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Signal') + theme(legend.position='none') + scale_color_manual(values=spectrum)
    })
    
    
    # Display of mixing matrix
    output$mixingPlot <- renderPlot({
        ggplot(melt(mixing_matrix()), aes(x=Var2, y=Var1, fill=value)) + geom_tile() + labs(x='Fluorophore', y='Channel') + theme(axis.text.x = element_text(angle=90))
    })
    
    
    # Display of mixing matrix condition number
    output$conditionText <- renderText({
        paste0("Condition number = ", round(kappa(mixing_matrix()), 0))
    })
    
    
    # Download button for mixing matrix as CSV
    output$downloadData <- downloadHandler(
        filename = function() {
            selected_names <- c(input$f1n, input$f2n, input$f3n, input$f4n, input$f5n, input$f6n, input$f7n, input$f8n)
            selected_names <- selected_names[selected_names != "**None**"]
            selected_names <- paste(selected_names, collapse="_")
            paste('Mixing_matrix_', selected_names, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(mixing_matrix(), file, row.names = FALSE)
        }
    )
}



shinyApp(ui = ui, server = server)
