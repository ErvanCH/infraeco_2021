library(shiny)
library(rgdal)
library(DT)
library(dygraphs)
library(xts)
library(leaflet)

# ui object
ui <- fluidPage(
    titlePanel(p("Spatial app", style = "color:#3474A7")),
    sidebarLayout(
        sidebarPanel(
            fileInput(
                inputId = "filemap",
                label = "Upload map. Choose shapefile",
                multiple = TRUE,
                accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj")
            ),
            selectInput(
                inputId = "variableselected",
                label = "Select variable",
                choices = c("cases", "population")
            ),
            selectInput(
                inputId = "yearselected",
                label = "Select year",
                choices = 1968:1988
            ),
            img(
                src = "Logo info flora-rgb.jpg",
                width = "70px", height = "70px"
            )
        ),
        
        mainPanel(
            leafletOutput(outputId = "map"),
            DTOutput(outputId = "table")
        )
    )
)

# server()
server <- function(input, output) {
    data <- reactive({
        req(input$filedata)
        read.csv(input$filedata$datapath)
    })
    
    map <- reactive({
        req(input$filemap)
        
        # shpdf is a data.frame with the name, size, type and
        # datapath of the uploaded files
        shpdf <- input$filemap
        
        # The files are uploaded with names
        # 0.dbf, 1.prj, 2.shp, 3.xml, 4.shx
        # (path/names are in column datapath)
        # We need to rename the files with the actual names:
        # fe_2007_39_county.dbf, etc.
        # (these are in column name)
        
        # Name of the temporary directory where files are uploaded
        tempdirname <- dirname(shpdf$datapath[1])
        
        # Rename files
        for (i in 1:nrow(shpdf)) {
            file.rename(
                shpdf$datapath[i],
                paste0(tempdirname, "/", shpdf$name[i])
            )
        }
        
        # Now we read the shapefile with readOGR() of rgdal package
        # passing the name of the file with .shp extension.
        
        # We use the function grep() to search the pattern "*.shp$"
        # within each element of the character vector shpdf$name.
        # grep(pattern="*.shp$", shpdf$name)
        # ($ at the end denote files that finish with .shp,
        # not only that contain .shp)
        map <- readOGR(paste(tempdirname,
                             shpdf$name[grep(pattern = "*.shp$", shpdf$name)],
                             sep = "/"
        ))
        map
    })
    
    output$table <- renderDT(data())
    
    output$timetrend <- renderDygraph({
        data <- data()
        dataxts <- NULL
        counties <- unique(data$county)
        for (l in 1:length(counties)) {
            datacounty <- data[data$county == counties[l], ]
            dd <- xts(
                datacounty[, input$variableselected],
                as.Date(paste0(datacounty$year, "-01-01"))
            )
            dataxts <- cbind(dataxts, dd)
        }
        colnames(dataxts) <- counties
        dygraph(dataxts) %>%
            dyHighlight(highlightSeriesBackgroundAlpha = 0.2) -> d1
        d1$x$css <- "
 .dygraph-legend > span {display:none;}
 .dygraph-legend > span.highlight { display: inline; }
 "
        d1
    })
    
    output$map <- renderLeaflet({
        if (is.null(data()) | is.null(map())) {
            return(NULL)
        }
        
        map <- map()
        data <- data()
        
        # Add data to map
        datafiltered <- data[which(data$year == input$yearselected), ]
        ordercounties <- match(map@data$NAME, datafiltered$county)
        map@data <- datafiltered[ordercounties, ]
        
        # Create variableplot
        map$variableplot <- as.numeric(
            map@data[, input$variableselected])
        
        # Create leaflet
        pal <- colorBin("YlOrRd", domain = map$variableplot, bins = 7)
        labels <- sprintf("%s: %g", map$county, map$variableplot) %>%
            lapply(htmltools::HTML)
        
        l <- leaflet(map) %>%
            addTiles() %>%
            addPolygons(
                fillColor = ~ pal(variableplot),
                color = "white",
                dashArray = "3",
                fillOpacity = 0.7,
                label = labels
            ) %>%
            leaflet::addLegend(
                pal = pal, values = ~variableplot,
                opacity = 0.7, title = NULL
            )
    })
}

# shinyApp()
shinyApp(ui = ui, server = server)