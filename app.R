packages <- c("shiny", "dplyr", "Seurat", "ggplot2")

for (package in packages) {
	if (!require(package, character.only = TRUE)) {
		install.packages(package, dependencies = TRUE)
		library(package, character.only = TRUE)
	}
}

# Stuff that is run once, when the Shiny app launches.
integrated.0h.1h_4h_16h_27h_36h_48h <- readRDS("data/integrated.0h.1h_4h_16h_27h_36h_48h.rds")
DefaultAssay(integrated.0h.1h_4h_16h_27h_36h_48h) <- "RNA"

#RENAME IDENTS TO figure 1 cluster names correlating with UMAP
cluster_ids <- c('0' = 'S1', '1' = 'S3', '2' = 'S3','3' = 'S1', '4' = 'DCT',
	'5' = 'Loop', '6' = 'S2', '7' = 'S2',
	'8' = 'S3', '9' = 'Endothelial',
	'10' = 'Mp/DC',
	'11' = 'S3', '12' = 'S1', '13' = 'Lymphocytes','14' = 'CNT', '15' = 'CD-PC',
	'16' = 'CD-IC', '17' = 'S3T2', '18'= 'Other', '19'= 'Neutrophil',
	'20'= 'Peri_Stromal', '21'= 'Prolif', '22'= 'Other', '23'= 'Other',
	'24'= 'Other', '25'= 'Other')
integrated.0h.1h_4h_16h_27h_36h_48h <- RenameIdents(integrated.0h.1h_4h_16h_27h_36h_48h, cluster_ids)
# Reorder the category names for display purposes.
integrated.0h.1h_4h_16h_27h_36h_48h@active.ident <- factor(integrated.0h.1h_4h_16h_27h_36h_48h@active.ident, levels = c("S1", "S2", "S3", "S3T2", "Loop", "CNT", "DCT", "CD-IC", "CD-PC", "Peri_Stromal", "Prolif", "Endothelial", "Mp/DC", "Lymphocytes", "Neutrophil", "Other"))
# Store the categories in metadata as well.
integrated.0h.1h_4h_16h_27h_36h_48h$category <- Idents(integrated.0h.1h_4h_16h_27h_36h_48h)



# Define UI for application
ui <- fluidPage(

	# Application title
	titlePanel("Mouse kidney single-cell expression"),

	# Sidebar with a text input for gene
	sidebarLayout(
		sidebarPanel(
			selectizeInput("gene_name",
				label = "Gene of interest:",
				choices = as.list(rownames(integrated.0h.1h_4h_16h_27h_36h_48h)),
				multiple = FALSE,
				selected = "Cd74"
			),
			actionButton(
				inputId = "submit_loc",
				label = "Submit gene"
			),
			checkboxInput(
				inputId = "timesplit",
				label = "Split by timepoint",
				value = FALSE
			),
			sliderInput(
				inputId = "pt_size",
				label = "Point size",
				min = 0.0,
				max = 0.5,
				value = 0.01
			),
			img(src = "umap_img.jpg",
				width = "100%",
				height = "auto",
				alt = "Cell categories"
			)
		),

		# Show one or more FeaturePlots
		mainPanel(
			uiOutput("plots") # This is the dynamic UI for the plots
		)
	),
	# Footer.
	hr(),
	fluidRow(
		column(width = 1, align = "center", img(src="iu_tab_web.png", width="100%")),
		column(width = 11,
			p( "Thomas McCarthy",
				br(),
				"Takashi Hato",
				a("thato@iu.edu", href="mailto:thato@iu.edu")
			)
		)
	)
)


server <- function(input, output) {
	# Read the gene name from the textInput when the Submit actionButton is used or the time split checkbox is changed.
	gene <- eventReactive(
		eventExpr = {input[["submit_loc"]] | input$timesplit},
		valueExpr = input$gene_name
	)

	# Generate a reactive list of plots that changes when the gene changes.
	plot_list <- reactive({
		FeaturePlot(integrated.0h.1h_4h_16h_27h_36h_48h,
			features = gene(),
			#features = input$gene_name,
			cols = c("lightgrey", "blue"),
			split.by = {if (input$timesplit) "hrs" else NULL},
			pt.size = input$pt_size,
			by.col = FALSE,
			combine = FALSE
		)
	})


	# Insert the right number of plot output objects into the web page
	output$plots <- renderUI({
		plot_output_list <- lapply(1:length(plot_list()), function(i) {
			plotname <- paste0("plot", i)
			plotOutput(plotname, height = 400, width = 500)
		})

		# Convert the list to a tagList - this is necessary for the list of items to display properly.
		do.call(tagList, plot_output_list)
	})

	observe({
		# Call renderPlot for each plot in the list.
		for (i in 1:length(plot_list())) {
			# Need local so that each item gets its own number. Without it, the value of i in the renderPlot() will be the same across all instances, because of when the expression is evaluated.
			local({
				my_i <- i
				# Give them the names the UI expects.
				plotname <- paste0("plot", my_i)
				output[[plotname]] <- renderPlot({plot_list()[[my_i]]})
			})
		}
	})
}

# Run the application
shinyApp(ui = ui, server = server)
