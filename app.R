library(shiny)
library(reshape2)
library(ggplot2)
library(plyr)

theme_set(theme_bw(base_size=20))
update_geom_defaults("line", list(size = 1.75))

dichroic_choices <- c('491/514/532/561/594/633/670', '458/491/514/532/561/594/633')

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

# c1d <- eval(parse(text=c1e))
# c2d <- eval(parse(text=c2e))
# c3d <- eval(parse(text=c3e))
# c4d <- eval(parse(text=c4e))
# c5d <- eval(parse(text=c5e))
# c6d <- eval(parse(text=c6e))
# c7d <- eval(parse(text=c7e))
# c8d <- eval(parse(text=c8e))

# channels_df <- data.frame(cbind(c1d, c2d, c3d, c4d, c5d, c6d, c7d, c8d))
# channels_df$wavelength <- dichroics_df$wavelength
# channels <- melt(channels_df, id.vars='wavelength')


# Load fpbase.org data
fpbase_data <- readRDS('fpbase_full_emission.rds')
fpbase_data <- fpbase_data[fpbase_data$Wavelength >= 450 & fpbase_data$Wavelength <= 750, ]
fpbase_data <- fpbase_data[order(fpbase_data$Fluor), ]
fp_names <- unique(fpbase_data$Fluor)


# Add blank spectrum
w <- unique(fpbase_data$Wavelength)
v <- 0 * w
n <- rep("**None**", length(w))
blank_spectrum <- data.frame(Wavelength=w, Emission=v, Fluor=n)
fpbase_data <- rbind(fpbase_data, blank_spectrum)
fp_names <- unique(fpbase_data$Fluor)


ui <- fluidPage(
    titlePanel("8 camera multispectral system"),
    
    sidebarLayout(
        sidebarPanel(
            selectInput("dset", "Dichroic set", choices=dichroic_choices, selected=dichroic_choices[1]),
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

    # Reactive functions    
    dichroics_df <- reactive({
        dset <- input$dset
        dichroics_filename <- ifelse(input$dset == dichroic_choices[1], 'dichroics.csv', 'dichroics_bluer.csv')
        
        dichroics_df <- read.csv(dichroics_filename)
        colnames(dichroics_df) <- c('wavelength', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7')
        dichroics_df[is.na(dichroics_df)] <- 0
        dichroics_df <- dichroics_df[dichroics_df$wavelength <= 750 & dichroics_df$wavelength >= 450, ]
        dichroics_df
    })
    
    
    dichroics <- reactive({
        dichroics_df <- dichroics_df()
        dichroics <- melt(dichroics_df, id.vars='wavelength')
        dichroics
    })
    
    
    channels_df <- reactive({
        dichroics <- dichroics()
        dichroics_df <- dichroics_df()

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
        channels_df
    })
    
    
    channels <- reactive({
        channels <- melt(channels_df(), id.vars='wavelength')
        channels
    })
    
    
    fluorophores <- reactive({
        fluorophores <- fpbase_data[fpbase_data$Fluor == input$f1n, ]
        if (input$f2n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$Fluor == input$f2n, ])
        if (input$f3n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$Fluor == input$f3n, ])
        if (input$f4n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$Fluor == input$f4n, ])
        if (input$f5n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$Fluor == input$f5n, ])
        if (input$f6n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$Fluor == input$f6n, ])
        if (input$f7n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$Fluor == input$f7n, ])
        if (input$f8n != "**None**") fluorophores <- rbind(fluorophores, fpbase_data[fpbase_data$Fluor == input$f8n, ])
        fluorophores
    })
    
    
    mixing_matrix <- reactive({
        f <- dcast(fluorophores(), Wavelength ~ Fluor, value.var='Emission')
        f[is.na(f)] <- 0
        channels_df <- channels_df()
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
        
        # This is the silly R way of dividing a matrix by the sum of its columns...
        fmat <- sweep(fmat, 2, colSums(fmat), '/')
        
        mmat <- cmat %*% fmat
        mmat
    })
    
    
    cs <- reactive({
        emission <- ddply(fluorophores(), .(Wavelength), summarise, value=sum(Emission))
        cs <- channels_df()
        
        if (nrow(emission) < nrow(cs)) {
            # We've dropped some wavelengths, so put them back with zero
            if (min(emission$Wavelength) != 450) {
                low_wavelengths <- 450:(min(emission$Wavelength) - 1)
                empty_data <- matrix(0, ncol=ncol(emission), nrow=length(low_wavelengths))
                empty_data[, 1] <- low_wavelengths
                colnames(empty_data) <- colnames(emission)
                emission <- rbind(as.data.frame(empty_data), as.matrix(emission))
            }
            if (max(emission$Wavelength) != 750) {
                high_wavelengths <- (max(emission$Wavelength) + 1):750
                empty_data <- matrix(0, ncol=ncol(emission), nrow=length(high_wavelengths))
                empty_data[, 1] <- high_wavelengths
                emission <- rbind(as.matrix(emission), empty_data)
            }
        }
        
        cs[, 1:8] <- cs[, 1:8] * emission$value
        cs <- melt(cs, id.vars='wavelength')
        cs
    })
    
    
    # Display of dichroic spectra
    output$dichroicPlot <- renderPlot({
        ggplot(dichroics(), aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Transmission') + theme(legend.position='none') + scale_color_brewer(palette='Set2') + coord_cartesian(xlim=c(450, 750))
    })
    
    
    # Display of camera channel spectra
    output$channelPlot <- renderPlot({
        ggplot(channels(), aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Transmission') + theme(legend.position='none') + scale_color_manual(values=spectrum) + coord_cartesian(xlim=c(450, 750))
    })
    
    
    # Display of fluorophore spectra
    output$fluorophorePlot <- renderPlot({
        ggplot(fluorophores(), aes(x=Wavelength, y=Emission, col=Fluor)) + geom_line() + labs(x='Wavelength / nm', y='Fluorescence', col='') + theme(legend.position='top') + scale_color_brewer(palette='Set2') + theme(legend.text=element_text(size=10)) + coord_cartesian(xlim=c(450, 750))
    })
    
    
    # Display of camera-frame fluorophore spectra
    output$cameraPlot <- renderPlot({
        ggplot(cs(), aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Signal') + theme(legend.position='none') + scale_color_manual(values=spectrum) + coord_cartesian(xlim=c(450, 750))
    })
    
    
    # Display of mixing matrix
    output$mixingPlot <- renderPlot({
        mix_mat <- mixing_matrix()
        rownames(mix_mat) <- c('Ch 1', 'Ch 2', 'Ch 3', 'Ch 4', 'Ch 5', 'Ch 6', 'Ch 7', 'Ch 8')
        ggplot(melt(mix_mat), aes(x=Var2, y=Var1, fill=value)) + geom_tile() + labs(x='Fluorophore', y='Channel', fill='') + theme(axis.text.x = element_text(angle=90)) + geom_text(aes(label=round(value, 3)), col='white')
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
            paste('PRISM_mixing_matrix_', selected_names, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(mixing_matrix(), file, row.names = FALSE)
        }
    )
}



shinyApp(ui = ui, server = server)
