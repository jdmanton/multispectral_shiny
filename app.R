library(shiny)
library(reshape2)
library(ggplot2)
library(plyr)
library(colorspace)
library(ggforce)
library(tiff)
library(png)
library(MASS)

# ggplot theme settings
theme_set(theme_bw(base_size=18))
update_geom_defaults("line", list(size = 1.75))


# List of PRISM dichroic / detector choices
dichroic_choices <- c('PRISM 491/514/532/561/594/633/670', 'PRISM 458/491/514/532/561/594/633', 'Zeiss QUASAR detector (reflection mode)')


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


# Load fpbase.org data (FPs + small molecule dyes)
fpbase_data <- readRDS('fpbase_full_emission.rds')
fpbase_data <- fpbase_data[fpbase_data$Wavelength >= 400 & fpbase_data$Wavelength <= 750, ]
fpbase_data <- fpbase_data[order(fpbase_data$Fluor), ]
fp_names <- unique(fpbase_data$Fluor)


# Add blank spectrum (needed for **None** option in UI)
w <- unique(fpbase_data$Wavelength)
v <- 0 * w
n <- rep("**None**", length(w))
qy <- rep(0, length(w))
extCoeff <- rep(0, length(w))
brightness <- rep(0, length(w))
blank_spectrum <- data.frame(Wavelength=w, Emission=v, Fluor=n, QY=qy, extCoeff=extCoeff, Brightness=brightness)
fpbase_data <- rbind(fpbase_data, blank_spectrum)
fp_names <- unique(fpbase_data$Fluor)

# Load SPECTRUM letters image
spectrum_letters <- readTIFF('spectrum_letters.tif', all=TRUE)
letters_vec <- t(sapply(spectrum_letters, c))


##################
# User interface #
##################
ui <- function(request) {
	fluidPage(
		titlePanel("Spectral unmixing explorer"),
		
		div(style="height:25px"),
		
		fluidRow(
			column(1, bookmarkButton(label="Copy link", title="Share current settings via URL")),
			column(2, downloadButton("downloadData", "Download mixing matrix")),
			column(3, checkboxInput("brightnessScale", "Scale fluorophore spectra by brightness", value = FALSE, width = NULL))
		),
		
		div(style="height:25px"),
		
		fluidRow(
			column(4, selectInput("dset", "Detector set", choices=dichroic_choices, selected=dichroic_choices[2]))
		),
		
		fluidRow(
			column(3, selectInput("f1n", "Fluorophore 1", choices=fp_names, selected='TagBFP')),
			column(3, selectInput("f2n", "Fluorophore 2", choices=fp_names, selected='Cerulean')),
			column(3, selectInput("f3n", "Fluorophore 3", choices=fp_names, selected='mAzamiGreen')),
			column(3, selectInput("f4n", "Fluorophore 4", choices=fp_names, selected='Citrine'))
		),
		
		fluidRow(
			column(3, selectInput("f5n", "Fluorophore 5", choices=fp_names, selected='mCherry')),
			column(3, selectInput("f6n", "Fluorophore 6", choices=fp_names, selected='iRFP670')),
			column(3, selectInput("f7n", "Fluorophore 7", choices=fp_names, selected='**None**')),
			column(3, selectInput("f8n", "Fluorophore 8", choices=fp_names, selected='**None**'))
		),
			
		fluidRow(
			column(6, 
				h3("Dichroic spectra"),
				plotOutput("dichroicPlot", height='200px'),
				h3("Fluorophore spectra"),
				plotOutput("fluorophorePlot", height='300px'),
				h3("Mixing matrix"),
				textOutput("conditionText"),
				plotOutput("mixingPlot", height='800px')
			),
			column(6,
				h3("Channel spectra"),
				plotOutput("channelPlot", height='200px'),
				h3("Detector spectra"),
				plotOutput("cameraPlot", height='200px'),
				h3("Ground truth image"),
				imageOutput("letters", height='64px'),
				h3("Mixed image"),
				imageOutput("mixed_letters", height='64px'),
				h3("Linearly unmixed image"),
				imageOutput("unmixed_letters", height='64px'),
				h3("Phasor plot"),
				plotOutput("phasorPlot", height='400px')
			),
		)
	)
}



####################
# Application code #
####################
server <- function(input, output) {
	
	# Reactive functions    
	dichroics_df <- reactive({
		dset <- input$dset
		dichroics_filename <- ifelse(input$dset == dichroic_choices[1], 'dichroics.csv', 'dichroics_bluer.csv')
		
		dichroics_df <- read.csv(dichroics_filename)
		colnames(dichroics_df) <- c('wavelength', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7')
		dichroics_df[is.na(dichroics_df)] <- 0
		dichroics_df <- dichroics_df[dichroics_df$wavelength <= 750 & dichroics_df$wavelength >= 400, ]
		
		if(input$dset == dichroic_choices[3]) {
			dichroics_df[, 2:ncol(dichroics_df)] <- 0
		}
		
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
		
		if(input$dset == 'Zeiss QUASAR detector (reflection mode)') {
			quasar <- read.csv('quasar_detector.csv')
			quasar <- cbind(quasar[, 2:ncol(quasar)], quasar$Wavelength)
			colnames(quasar)[ncol(quasar)] <- 'wavelength'
			channels_df <- quasar
		}
		
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
		cmat <- t(as.matrix(channels_df[, 1:(ncol(channels_df) - 1)]))
		
		if (nrow(f) < ncol(cmat)) {
			# We've dropped some wavelengths, so put them back with zero
			if (min(f$Wavelength) != 400) {
				low_wavelengths <- 400:(min(f$Wavelength) - 1)
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
			if (min(emission$Wavelength) != 400) {
				low_wavelengths <- 400:(min(emission$Wavelength) - 1)
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
		
		cs[, 1:(ncol(cs) - 1)] <- cs[, 1:(ncol(cs) - 1)] * emission$value
		cs <- melt(cs, id.vars='wavelength')
		cs
	})
	
	
	calc_s <- function(wavelengths, intensities, harmonic) {
		sine_function <- sin(2 * pi * harmonic * (wavelengths - 400) / (750 - 400))
		s <- sum(sine_function * intensities) / sum(intensities)
		s
	}
	
	calc_g <- function(wavelengths, intensities, harmonic) {
		cosine_function <- cos(2 * pi * harmonic * (wavelengths - 400) / (750 - 400))
		g <- sum(cosine_function * intensities) / sum(intensities)
		g
	}
	
	
	phasors <- reactive({
		fluors <- fluorophores()
		
		harmonics <- seq(1, 4)
		phasor_harmonics <- lapply(harmonics, function(y) {
			phasors_s <- sapply(unique(fluors$Fluor), function (x) {
				ws <- fluors[fluors$Fluor == x, 'Wavelength']
				intensities <- fluors[fluors$Fluor == x, 'Emission']
				s <- calc_s(ws, intensities, y)
			})
			phasors_g <- sapply(unique(fluors$Fluor), function (x) {
				ws <- fluors[fluors$Fluor == x, 'Wavelength']
				intensities <- fluors[fluors$Fluor == x, 'Emission']
				g <- calc_g(ws, intensities, y)
			})
			phasors <- data.frame(s=phasors_s, g=phasors_g, fluor=names(phasors_s), harmonic=paste0("Harmonic ", y))
		})
		
		phasor_harmonics <- do.call(rbind, phasor_harmonics)
	})
	
	
	mixed_letters <- reactive({
		mmat <- mixing_matrix()
		letters_vec_subset <- letters_vec[1:ncol(mmat), ]
		mixed_letters <- mmat %*% letters_vec_subset
	})
	
	unmixed_letters <- reactive({
		mmat <- mixing_matrix()
		mixed <- mixed_letters()
		unmixed_letters <- ginv(mmat) %*% mixed
	})
	
	
	# Display of dichroic spectra
	output$dichroicPlot <- renderPlot({
		ggplot(dichroics(), aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Transmission') + theme(legend.position='none') + scale_color_brewer(palette='Set2') + coord_cartesian(xlim=c(400, 750))
	})
	
	
	# Display of camera channel spectra
	output$channelPlot <- renderPlot({
		channels <- channels()
		ggplot(channels, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Transmission') + theme(legend.position='none') + scale_color_manual(values=rev(rainbow_hcl(length(levels(channels$variable))))) + coord_cartesian(xlim=c(400, 750))
	})
	
	
	# Display of fluorophore spectra
	output$fluorophorePlot <- renderPlot({
		fluors <- fluorophores()
		if(input$brightnessScale) {
			fluors$Emission <- fluors$Emission * fluors$Brightness
		}
		ggplot(fluors, aes(x=Wavelength, y=Emission, col=Fluor)) + geom_line() + labs(x='Wavelength / nm', y='Fluorescence', col='') + theme(legend.position='top') + scale_color_brewer(palette='Set2') + coord_cartesian(xlim=c(400, 750))
	})
	
	
	# Display of camera-frame fluorophore spectra
	output$cameraPlot <- renderPlot({
		cs <- cs()
		ggplot(cs, aes(x=wavelength, y=value, col=variable)) + geom_line() + labs(x='Wavelength / nm', y='Signal') + theme(legend.position='none') + scale_color_manual(values=rev(rainbow_hcl(length(levels(cs$variable))))) + coord_cartesian(xlim=c(400, 750))
	})
	
	
	# Display of mixing matrix
	output$mixingPlot <- renderPlot({
		mix_mat <- mixing_matrix()
		rownames(mix_mat) <- paste0("Ch ", seq(1, nrow(mix_mat)))
		ggplot(melt(mix_mat), aes(x=Var2, y=Var1, fill=value)) + geom_tile(show.legend=FALSE) + labs(x='Fluorophore', y='Channel', fill='') + theme(axis.text.x = element_text(angle=90)) + geom_text(aes(label=round(value, 3)), col='white')
	})
	
	
	# Display of mixing matrix condition number
	output$conditionText <- renderText({
		paste0("Condition number = ", round(kappa(mixing_matrix()), 0))
	})
	
	
	# Display of phasor plot
	output$phasorPlot <- renderPlot({
		p <- phasors()
		ggplot(p, aes(x=g, y=s, col=fluor)) + geom_point(size=5) + lims(x=c(-1, 1), y=c(-1, 1)) + labs(x='G', y='S', col='Fluorophore') + scale_color_brewer(palette='Set2') + geom_circle(inherit.aes=FALSE, aes(x0=0, y0=0, r=1), linetype='dashed', lwd=0.75) + theme(aspect.ratio=1, legend.position='top') + facet_grid(. ~ harmonic)
	})
	
	
	# Display of ground truth letters
	output$letters <- renderImage({
		letters <- letters_vec
		unmixed_letters <- unmixed_letters()
		ims <- lapply(1:nrow(unmixed_letters), function(x) {
			im <- matrix(letters[x, ], 64, 64)
		})
		ims <- do.call(cbind, ims)
		outfile <- tempfile(fileext=".tif")
		writePNG(ims / max(ims), target=outfile)
		list(src=outfile, contentType = "image/png", width=64*nrow(unmixed_letters), height=64)
	}, deleteFile=TRUE)
	
	
	# Display of mixed letters
	output$mixed_letters <- renderImage({
		mixed_letters <- mixed_letters()
		ims <- lapply(1:nrow(mixed_letters), function(x) {
			im <- matrix(mixed_letters[x, ], 64, 64)
		})
		ims <- do.call(cbind, ims)
		outfile <- tempfile(fileext=".tif")
		writePNG(ims / max(ims), target=outfile)
		list(src=outfile, contentType = "image/png", width=64*nrow(mixed_letters), height=64)
	}, deleteFile=TRUE)
	
	
	# Display of linearly unmixed letters
	output$unmixed_letters <- renderImage({
		unmixed_letters <- unmixed_letters()
		ims <- lapply(1:nrow(unmixed_letters), function(x) {
			im <- matrix(unmixed_letters[x, ], 64, 64)
		})
		ims <- do.call(cbind, ims)
		outfile <- tempfile(fileext=".tif")
		writePNG(ims / max(ims), target=outfile)
		list(src=outfile, contentType = "image/png", width=64*nrow(unmixed_letters), height=64)
	}, deleteFile=TRUE)
	
	
	# Download button for mixing matrix as CSV
	output$downloadData <- downloadHandler(
		filename = function() {
			selected_names <- c(input$f1n, input$f2n, input$f3n, input$f4n, input$f5n, input$f6n, input$f7n, input$f8n)
			selected_names <- selected_names[selected_names != "**None**"]
			selected_names <- paste(selected_names, collapse="_")
			filename <- ifelse(input$dset == 'Zeiss QUASAR detector (reflection mode)', paste('QUASAR_mixing_matrix_', selected_names, ".csv", sep = ""), paste('PRISM_mixing_matrix_', selected_names, ".csv", sep = ""))
			filename
		},
		content = function(file) {
			write.csv(mixing_matrix(), file, row.names = FALSE)
		}
	)
}



shinyApp(ui = ui, server = server, enableBookmarking = "url")
