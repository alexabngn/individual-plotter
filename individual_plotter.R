# Load necessary libraries
{
library(shiny)
library(plotly)
library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(hms)
library(moveHMM)
library(car)
library(purrr)
library(dunn.test)
library(tidyr)
library(stringr)
library(plyr)
library(REdaS)
library(circular)
library(discr)
library(ggnewscale)
library(ggpubr)
library(gganimate)
library(grid)
library(gifski)
library(gridExtra)
}

# Define UI
{
ui <- fluidPage(
  titlePanel("Individual Trajectory Plotter"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload your data file", accept = c(".RData", ".csv")),
      uiOutput("video_selector"),
      actionButton("plot_button", "Plot Trajectories"),
      textOutput("message")  # Add a text output for messages
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data Review", verbatimTextOutput("data_review")),
        tabPanel("All Trajectories", plotOutput("all_trajectory", height = "800px")),
        tabPanel("Trajectory", plotOutput("individual_trajectory", height = "800px")),
        tabPanel("Heat Map", plotOutput("individual_heat_map", height = "900px")),
        tabPanel("Trajectory (color)", plotOutput("individual_by_color")),
        tabPanel("All Speed", plotOutput("all_speed", height = "800px")),
        tabPanel("Speed", plotOutput("speed_analysis", height = "600px")),
        tabPanel("Speed Distributions (onoff)", plotOutput("speed_distribution_onoff", height = "600px")),
        tabPanel("Speed Distributions (color)", plotOutput("speed_distribution_color", height = "600px")),
        tabPanel("Speed Histogram", plotOutput("speed_histogram", height = "600px")),
        tabPanel("Speed Statistics", verbatimTextOutput("speed_stats")),
        tabPanel("PosRel2Basin Shape", plotOutput("posrel2basin_shape", height = "600px")),
        tabPanel("PosRel2Basin Plot", plotOutput("posrel2basin_plot", height = "600px")),
        tabPanel("PosRel2Basin Stats", verbatimTextOutput("posrel2basin_stats")),
        tabPanel("PosRel2Light Plot", plotOutput("posrel2light_plot", height = "600px")),
        tabPanel("PosRel2Light Stats", verbatimTextOutput("posrel2light_stats")),
        tabPanel("MoveRel2Basin Plot", plotOutput("moverel2basin_plot", height = "600px")),
        tabPanel("MoveRel2Basin Stats", verbatimTextOutput("moverel2basin_stats")),
        tabPanel("MoveRel2Light Plot", plotOutput("moverel2light_plot", height = "600px")),
        tabPanel("MoveRel2Light Stats", verbatimTextOutput("moverel2light_stats")),
        tabPanel("Animation", imageOutput("animation"))

      )
    )
  )
)
}

# Define server logic
{
server <- function(input, output, session) {
  
  data <- reactiveVal(NULL)
  
  ########################################
  # Selection of file
  ########################################
  observeEvent(input$file, {
    file <- input$file
    ext <- tools::file_ext(file$datapath)
    if (ext == "RData") {
      load(file$datapath)
      data(d)
    } else if (ext == "csv") {
      d <- read.csv(file$datapath)
      data(d)
    } else {
      output$message <- renderText("Unsupported file type.")
      return(NULL)
    }
    updateSelectInput(session, "video", choices = unique(d$id))
  })
  
  ########################################
  # Selection of ID
  ########################################
  output$video_selector <- renderUI({
    selectInput("video", "Select ID", choices = unique(data()$id))
  })
  
  observeEvent(input$plot_button, {
    req(data())
    
    ########################################
    # Prepare data
    ########################################
    d <- data()
    indiv_id <- input$video
    indiv_tracks <- dplyr::filter(d, id == indiv_id)
    
    ########################################
    # Prepare color scales
    ########################################
    color_scale_color <- scale_color_manual(values = c("Red" = "red", "Green" = "green", "Blue" = "blue", "White" = "grey"))
    color_scale_fill <- scale_fill_manual(values = c("Red" = "red", "Green" = "green", "Blue" = "blue", "White" = "grey"))
    onoff_scale_fill <- scale_fill_manual(values = c("on" = "#cfc75b", "off" = "#6e6e6e"))
    
    ########################################
    # Organize
    ########################################
    d$onoff <- factor(d$onoff, levels = c('on', 'off'))
    d$color <- factor(d$color, levels = c("Red", "Green", "Blue", "White", "dark"))
    
    ########################################
    # Functions
    ########################################
    calculate_speed <- function(d) {
      d %>%
        dplyr::arrange(time) %>%
        dplyr::mutate(
          dx = c(NA, diff(x)),
          dy = c(NA, diff(y)),
          time_diff = c(NA, diff(as.numeric(time)))
        ) %>%
        dplyr::mutate(
          dist = sqrt(dx^2 + dy^2),
          speed = dist / time_diff
        )
    }
    
    ########################################
    # Calculate speed
    ########################################
    {
    d <- d %>%
      group_by(id) %>%
      group_modify(~ calculate_speed(.x)) %>%
      ungroup() %>%
      dplyr::select(-time_diff, -dist) %>%
      arrange(id, time) %>%
      group_by(id) %>%
      dplyr::mutate(bin = floor((time - min(time)) / 10)) # Define speed bins
    
    d_speed_bin <- d %>%
      group_by(id, bin) %>%
      dplyr::summarise(mean_speed = mean(speed), .groups = 'drop')
    
    d_speed <- d %>%
      dplyr::select(id, time, bin, onoff, color, x, y) %>%
      dplyr::distinct(id, bin, .keep_all = TRUE) %>%
      left_join(d_speed_bin, by = c("id", "bin"))
    
    indiv_tracks_speed_raw <- dplyr::filter(d, id == indiv_id)
    indiv_tracks_speed_binned <- dplyr::filter(d_speed, id == indiv_id)
    }
    
    ########################################
    # Calculate position rel 2 basin
    ########################################
    d$anglerel2basin <- with(d, {
      angle_rad <- atan2(y, x)
      angle_deg <- angle_rad * (180 / pi)
      as.bearing(as.bearing(-as.bearing(angle_deg))+90) 
    })
    
    {
    d_posrel2basin_binned <- dplyr::filter(d, id == indiv_id)
    
    # Binning anglerel2basin
    d_posrel2basin_binned$var <- round(d_posrel2basin_binned$anglerel2basin / 10) * 10
    
    d_posrel2basin_binned <- ddply(d_posrel2basin_binned, ~onoff + color + var, function(x) {
      data.frame(x, count = 1:nrow(x))
    })
    d_posrel2basin_binned <- d_posrel2basin_binned %>% arrange(id, time)
    d_posrel2basin_binned$count <- 50 + d_posrel2basin_binned$count 
    }
    
    ########################################
    # Calculate position rel 2 light
    ########################################
    d$color_angle <- with(d, {
      angle_rad <- atan2(y = y_light, x = x_light)
      angle_deg <- angle_rad * (180 / pi)
      as.bearing(as.bearing(-as.bearing(angle_deg))+90) 
    })
    
    {
    d <- dplyr::mutate(d, color_angle = case_when(color_angle == 90 ~ NA, TRUE ~ color_angle))
    d$anglerel2light <- as.bearing(d$anglerel2basin - d$color_angle)
    
    # Binning anglerel2light
    d_posrel2light_binned <- dplyr::filter(d, id == indiv_id, !is.na(color_angle))
    
    d_posrel2light_binned$var <- round(d_posrel2light_binned$anglerel2light / 10) * 10
    
    d_posrel2light_binned <- ddply(d_posrel2light_binned, ~onoff + color + var, function(x) {
      data.frame(x, count = 1:nrow(x))
    })
    d_posrel2light_binned <- d_posrel2light_binned %>% arrange(id, time)
    d_posrel2light_binned$count <- 50 + d_posrel2light_binned$count 
    }
    
    ########################################
    # Calculate move rel 2 basin
    ########################################
    {
    d <- d %>%
      group_by(id) %>%
      arrange(id, time) %>%
      dplyr::mutate(
        next_x = lead(x), 
        next_y = lead(y), 
        dx_angle = next_x - x, 
        dy_angle = next_y - y
      ) %>%
      ungroup()
    
    d$moverel2basin <- with(d, {
      angle_rad <- atan2(dy_angle, dx_angle)
      angle_deg <- angle_rad * (180 / pi)
      as.bearing(as.bearing(-as.bearing(angle_deg))+90) 
    })
    }
    
    {
      d_moverel2basin_binned <- dplyr::filter(d, id == indiv_id)
      
      # Binning anglerel2basin
      d_moverel2basin_binned$var <- round(d_moverel2basin_binned$moverel2basin / 10) * 10
      
      d_moverel2basin_binned <- ddply(d_moverel2basin_binned, ~onoff + color + var, function(x) {
        data.frame(x, count = 1:nrow(x))
      })
      d_moverel2basin_binned <- d_moverel2basin_binned %>% arrange(id, time)
      d_moverel2basin_binned$count <- 50 + d_moverel2basin_binned$count 
    }
    
    ########################################
    # Calculate move rel 2 light
    ########################################
    {
    d$moverel2light <- as.bearing(d$moverel2basin - d$color_angle)
    d <- dplyr::select(d, -dx_angle, -dy_angle, -type, -next_x, -next_y)
    }
    
    {
      # Binning anglerel2light
      d_moverel2light_binned <- dplyr::filter(d, id == indiv_id, !is.na(color_angle))
      
      d_moverel2light_binned$var <- round(d_moverel2light_binned$moverel2light / 10) * 10
      
      d_moverel2light_binned <- ddply(d_moverel2light_binned, ~onoff + color + var, function(x) {
        data.frame(x, count = 1:nrow(x))
      })
      d_moverel2light_binned <- d_moverel2light_binned %>% arrange(id, time)
      d_moverel2light_binned$count <- 50 + d_moverel2light_binned$count 
    }
    
    ########################################
    # Outputs to show on GUI
    ########################################
    # Data review
    output$data_review <- renderPrint({
      list(
        COLUMN_NAMES = names(indiv_tracks),
        DATA = indiv_tracks,
        VIDEO_COUNT = dplyr::count(indiv_tracks, id),
        ONOFF_COUNT_PER_VIDEO = dplyr::count(indiv_tracks, id, onoff),
        COLOR_COUNT_PER_VIDEO = dplyr::count(indiv_tracks, id, color)
      )
    })
    
    # Plot all trajectories
    output$all_trajectory <- renderPlot({
      ggplot(d, aes(x, y, color = color)) +
        geom_path() +
        facet_wrap(~id+onoff) +
        color_scale_color +
        theme_minimal() 
    })
    
    # Plot individual trajectory
    output$individual_trajectory <- renderPlot({
      ggplot(indiv_tracks, aes(x, y)) +
        geom_path() +
        theme_minimal()
    })
    
    # Plot individual heat map
    output$individual_heat_map <- renderPlot({
      ggplot(indiv_tracks, aes(x = x, y = y)) +
        stat_density_2d(aes(fill = stat(density)), geom = "raster", contour = FALSE, adjust = 1) +
        scale_fill_gradientn(colors = c("White" = "grey", "Blue" = "blue", "Green" = "green", "Red" = "red"), 
                             values = scales::rescale(c(0, 0.25, 0.75, 1)), 
                             name = "Density") +
        geom_point(aes(x_light, y_light, color = color), size = 4) +
        coord_fixed() +
        labs(x = "x coordinate (cm)", y = "y coordinate (cm)") +
        theme_minimal() +
        facet_wrap(~color, ncol = 2) +
        ggtitle("Heat map of fish trajectories per color") +
        theme(strip.text = element_text(size = 12),
              legend.position = "none")
    })
    
    # Plot individual trajectory by color
    output$individual_by_color <- renderPlot({
      ggplot(indiv_tracks, aes(x,y, colour = color)) +
        geom_path() +
        geom_point(aes(x_light, y_light, color = color), size = 8) +
        geom_segment(aes(x = x, y = y, xend = lead(x), yend = lead(y), colour = color),
                     arrow = arrow(length = unit(0.1, "inches"), type = "open")) +
        facet_wrap(~color) +
        color_scale_color 
    })
    
    # Plot all speed bins vs raw
    output$all_speed <- renderPlot({
      ggplot() +
        geom_path(data = d, aes(x = as_hms(time), y = speed), alpha = 0.5, size = 1) +
        geom_path(data = d_speed, aes(x = as_hms(time), y = mean_speed), alpha = 0.5, size = 1, color = "red") +
        xlab("Time (hh:mm:ss)") +
        ylab("Speed (pixel distance/second)") +
        ggtitle("Raw speed vs. Binned speed of all seabass") +
        theme(strip.text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +
        facet_wrap(~id)
    })
    
    # Plot speed bins vs raw
    output$speed_analysis <- renderPlot({
      ggplot() +
        geom_path(data = indiv_tracks_speed_raw, aes(x = as_hms(time), y = speed), alpha = 0.5, size = 1) +
        geom_path(data = indiv_tracks_speed_binned, aes(x = as_hms(time), y = mean_speed), alpha = 0.5, size = 1, color = "red") +
        xlab("Time (hh:mm:ss)") +
        ylab("Speed (pixel distance/second)") +
        ggtitle("Raw speed vs. Binned speed of an individual seabass") +
        theme(strip.text = element_text(size = 24), axis.text = element_text(size = 24), axis.title = element_text(size = 24))
    })
    
    # Plot speed distributions per onoff
    output$speed_distribution_onoff <- renderPlot({
      ggplot(indiv_tracks_speed_binned, aes(x = onoff, y = mean_speed, fill = onoff)) + 
        geom_boxplot() +
        geom_jitter(position = position_jitter(width = 0.1), color = "#5680b0", size = 3, alpha = 0.5) +
        ggtitle("Speed distribution of an individual seabass in absence and presence of light") +
        xlab("Event") +
        ylab("Speed (distance/time)") +
        onoff_scale_fill +
        theme_minimal()
    })
    
    # Plot speed distributions per color
    output$speed_distribution_color <- renderPlot({
      ggplot(indiv_tracks_speed_binned, aes(x = color, y = mean_speed, fill = color)) + 
        geom_boxplot() +
        geom_jitter(position = position_jitter(width = 0.1), color = "#5680b0", size = 3, alpha = 0.5) +
        ggtitle("Speed distribution of an individual seabass in absence and presence of light") +
        xlab("Event") +
        ylab("Speed (distance/time)") +
        color_scale_fill +
        theme_minimal()
    })
    
    # Plot speed histogram
    output$speed_histogram <- renderPlot({
      ggplot(indiv_tracks_speed_binned, aes(x = mean_speed)) +
        geom_histogram(binwidth = 1, aes(fill = color), alpha = 0.7) +
        labs(title = "Histogram of Mean Speed", x = "Mean Speed", y = "Frequency") +
        theme_minimal() +
        facet_wrap(~color, scales = "free_x")
    })
    
    # Speed statistics
    output$speed_stats <- renderPrint({
      kruskal_result <- kruskal.test(mean_speed ~ onoff, data = indiv_tracks_speed_binned)
      if (kruskal_result$p.value < 0.05) {
        dunn_result <- dunn.test(indiv_tracks_speed_binned$mean_speed, indiv_tracks_speed_binned$onoff, method = "bonferroni")
        significant_comparisons <- which(dunn_result$P.adjusted < 0.05)
        for (i in significant_comparisons) {
          comparison <- strsplit(dunn_result$comparisons[i], " - ")[[1]]
          group1 <- comparison[1]
          group2 <- comparison[2]
          mean_speed_group1 <- mean(indiv_tracks_speed_binned$mean_speed[indiv_tracks_speed_binned$onoff == group1], na.rm = TRUE)
          mean_speed_group2 <- mean(indiv_tracks_speed_binned$mean_speed[indiv_tracks_speed_binned$onoff == group2], na.rm = TRUE)
          if (mean_speed_group1 < mean_speed_group2) {
            slower_group <- group1
            faster_group <- group2
          } else {
            slower_group <- group2
            faster_group <- group1
          }
          cat(slower_group, "has a significantly slower mean speed than", faster_group, ".\n")
        }
      } else {
        cat("There is no significant difference between the mean speeds.\n")
      }
      
      fligner_result <- fligner.test(mean_speed ~ onoff, data = indiv_tracks_speed_binned)
      if (fligner_result$p.value < 0.05) {
        cat("There is a significant difference in the variability of speed between 'on' and 'off'.\n")
      } else {
        cat("There is no significant difference in the variability of speed between 'on' and 'off'.\n")
      }
      
      kruskal_result <- kruskal.test(mean_speed ~ color, data = indiv_tracks_speed_binned)
      if (kruskal_result$p.value < 0.05) {
        dunn_result <- dunn.test(indiv_tracks_speed_binned$mean_speed, indiv_tracks_speed_binned$color, method = "bonferroni")
        significant_comparisons <- which(dunn_result$P.adjusted < 0.05)
        for (i in significant_comparisons) {
          comparison <- strsplit(dunn_result$comparisons[i], " - ")[[1]]
          group1 <- comparison[1]
          group2 <- comparison[2]
          mean_speed_group1 <- mean(indiv_tracks_speed_binned$mean_speed[indiv_tracks_speed_binned$color == group1], na.rm = TRUE)
          mean_speed_group2 <- mean(indiv_tracks_speed_binned$mean_speed[indiv_tracks_speed_binned$color == group2], na.rm = TRUE)
          if (mean_speed_group1 < mean_speed_group2) {
            slower_group <- group1
            faster_group <- group2
          } else {
            slower_group <- group2
            faster_group <- group1
          }
          cat(slower_group, "has a significantly slower mean speed than", faster_group, ".\n")
        }
      } else {
        cat("There is no significant difference between the mean speeds among the color groups.\n")
      }
    })
    
    # PosRel2Basin shape
    output$posrel2basin_shape <- renderPlot({
      ggplot(data = dplyr::filter(indiv_tracks, !color == "dark"),
             aes(x = x, y = y, color = color)) +
        geom_path(size = 0.5) +
        geom_point(size = 2.5) +
        geom_point(aes(x = x_light, y = y_light, color = color), size = 8) +
        geom_segment(aes(xend = 0, yend = 0, color = color),
                     linetype = "solid", alpha = 0.5, size = 1) +
        annotate("point", x = 0, y = 0, size = 6, color = "black") +
        labs(x = "x coordinate (cm)", y = "y Coordinate (cm)") +
        facet_wrap(~color) +
        color_scale_color +
        theme_minimal() +
        theme(legend.position = "none",
              panel.background = element_rect(color = "white"),
              panel.grid.major.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_line(color = "grey", size = 0.1),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 24),
              strip.text = element_text(size = 24),
              strip.background = element_rect(fill = "white", colour = "white"),
              plot.background = element_rect(color = "white"))
    })
    
    # PosRel2Basin plot
    output$posrel2basin_plot <- renderPlot({
      stat_rayleigh_posrel2basin <- plyr::ddply(d_posrel2basin_binned, ~onoff + color, function(x) {
        r <- rayleigh.test(as.bearing(x$anglerel2basin))
        if (r$p.value < 0.05) {
          s <- TRUE
        } else {
          s <- FALSE
        }
        mean <- mean.circular(as.bearing(x$anglerel2basin))
        data.frame(id=unique(x$id), mean = mean, R = r$statistic, p.value = r$p.value, signif = s)
      })
      
      ggplot() +
        geom_point(aes(x = var, y = count, color = color),
                   data = dplyr::filter(d_posrel2basin_binned, !color == "dark"), size = 2) +
        coord_polar() +
        geom_segment(aes(x = mean, y = 0, xend = mean, yend = 50, linetype = signif, color = signif),
                     data = dplyr::filter(stat_rayleigh_posrel2basin, !color == "dark"),
                     arrow = arrow(length = unit(0.05, "npc"), type = "closed", ends = "last"), size = 1) +
        color_scale_color +
        scale_x_continuous(
          "Position relative to basin center (°)",
          limits = c(0, 360),
          breaks = c(0, 90, 180, 270)
        ) +
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 10)) +
        scale_linetype_manual(values = c("FALSE" = NA, "TRUE" = "solid")) +
        facet_wrap(~color) +
        theme_void() +
        theme(panel.background = element_rect(color = "white"),
              panel.grid.major.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_line(color = "grey", size = 0.1),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 24),
              strip.text = element_text(size = 24),
              strip.background = element_rect(fill = "white", colour = "white"),
              legend.position = "none",
              plot.background = element_rect(color = "white"))
    })
    
    # PosRel2Basin stats
    output$posrel2basin_stats <- renderPrint({
      
      stat_rayleigh_posrel2basin <- plyr::ddply(d_posrel2basin_binned, ~onoff + color, function(x) {
        r <- rayleigh.test(as.bearing(x$anglerel2basin))
        if (r$p.value < 0.05) {
          s <- TRUE
        } else {
          s <- FALSE
        }
        mean <- mean.circular(as.bearing(x$anglerel2basin))
        data.frame(id=unique(x$id), mean = mean, R = r$statistic, p.value = r$p.value, signif = s)
      })
      
      print(stat_rayleigh_posrel2basin)
    })
    
    # PosRel2Light plot
    output$posrel2light_plot <- renderPlot({
      
      stat_rayleigh_posrel2light <- plyr::ddply(d_posrel2light_binned, ~onoff + color, function(x) {
        r <- rayleigh.test(as.bearing(x$anglerel2light))
        if (r$p.value < 0.05) {
          s <- TRUE
        } else {
          s <- FALSE
        }
        mean <- mean.circular(as.bearing(x$anglerel2light))
        data.frame(id=unique(x$id), mean = mean, R = r$statistic, p.value = r$p.value, signif = s)
      })
      
      ggplot() +
        geom_point(aes(x = var, y = count, color = color),
                   data = dplyr::filter(d_posrel2light_binned, !color == "dark"), size = 2) +
        coord_polar() +
        geom_segment(aes(x = mean, y = 0, xend = mean, yend = 50, linetype = signif, color = signif),
                     data = dplyr::filter(stat_rayleigh_posrel2light, !color == "dark"),
                     arrow = arrow(length = unit(0.05, "npc"), type = "closed", ends = "last"), size = 1) +
        color_scale_color +
        scale_x_continuous(
          "Position relative to light (°)",
          limits = c(0, 360),
          breaks = c(0, 90, 180, 270),
          labels = c("Near Light", "", "Far", "")
        ) +
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 10)) +
        scale_linetype_manual(values = c("FALSE" = NA, "TRUE" = "solid")) +
        facet_wrap(~color) +
        theme_void() +
        theme(panel.background = element_rect(color = "white"),
              panel.grid.major.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_line(color = "grey", size = 0.1),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 24),
              strip.text = element_text(size = 24),
              strip.background = element_rect(fill = "white", colour = "white"),
              legend.position = "none",
              plot.background = element_rect(color = "white"))
    })
    
    # PosRel2Light stats
    output$posrel2light_stats <- renderPrint({
      
      stat_rayleigh_posrel2light <- plyr::ddply(dplyr::filter(d_posrel2light_binned, !is.na(color_angle)), ~onoff + color, function(x) {
        r <- rayleigh.test(as.bearing(x$anglerel2light))
        if (r$p.value < 0.05) {
          s <- TRUE
        } else {
          s <- FALSE
        }
        mean <- mean.circular(as.bearing(x$anglerel2light))
        data.frame(id=unique(x$id), mean = mean, R = r$statistic, p.value = r$p.value, signif = s)
      })
      
      print(stat_rayleigh_posrel2light)
    })
    
    # MoveRel2Basin plot
    output$moverel2basin_plot <- renderPlot({
      stat_rayleigh_moverel2basin<- plyr::ddply(d_moverel2basin_binned, ~onoff + color, function(x) {
        r <- rayleigh.test(as.bearing(x$moverel2basin))
        if (r$p.value < 0.05) {
          s <- TRUE
        } else {
          s <- FALSE
        }
        mean <- mean.circular(as.bearing(x$moverel2basin))
        data.frame(id=unique(x$id), mean = mean, R = r$statistic, p.value = r$p.value, signif = s)
      })
      
      ggplot() +
        geom_point(aes(x = var, y = count, color = color),
                   data = dplyr::filter(d_moverel2basin_binned, !color == "dark"), size = 2) +
        coord_polar() +
        geom_segment(aes(x = mean, y = 0, xend = mean, yend = 50, linetype = signif, color = signif),
                     data = dplyr::filter(stat_rayleigh_moverel2basin, !color == "dark"),
                     arrow = arrow(length = unit(0.05, "npc"), type = "closed", ends = "last"), size = 1) +
        color_scale_color +
        scale_x_continuous(
          "Movement relative to basin center (°)",
          limits = c(0, 360),
          breaks = c(0, 90, 180, 270)
        ) +
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 10)) +
        scale_linetype_manual(values = c("FALSE" = NA, "TRUE" = "solid")) +
        facet_wrap(~color) +
        theme_void() +
        theme(panel.background = element_rect(color = "white"),
              panel.grid.major.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_line(color = "grey", size = 0.1),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 24),
              strip.text = element_text(size = 24),
              strip.background = element_rect(fill = "white", colour = "white"),
              legend.position = "none",
              plot.background = element_rect(color = "white"))
    })
    
    # MoveRel2Basin stats
    output$moverel2basin_stats <- renderPrint({
      
      stat_rayleigh_moverel2basin <- plyr::ddply(d_moverel2basin_binned, ~onoff + color, function(x) {
        r <- rayleigh.test(as.bearing(x$moverel2basin))
        if (r$p.value < 0.05) {
          s <- TRUE
        } else {
          s <- FALSE
        }
        mean <- mean.circular(as.bearing(x$moverel2basin))
        data.frame(id=unique(x$id), mean = mean, R = r$statistic, p.value = r$p.value, signif = s)
      })
      
      print(stat_rayleigh_moverel2basin)
    })
    
    # MoveRel2Light plot
    output$moverel2light_plot <- renderPlot({
      stat_rayleigh_moverel2light<- plyr::ddply(d_moverel2light_binned, ~onoff + color, function(x) {
        r <- rayleigh.test(as.bearing(x$moverel2light))
        if (r$p.value < 0.05) {
          s <- TRUE
        } else {
          s <- FALSE
        }
        mean <- mean.circular(as.bearing(x$moverel2light))
        data.frame(id=unique(x$id), mean = mean, R = r$statistic, p.value = r$p.value, signif = s)
      })
      
      ggplot() +
        geom_point(aes(x = var, y = count, color = color),
                   data = dplyr::filter(d_moverel2light_binned, !color == "dark"), size = 2) +
        coord_polar() +
        geom_segment(aes(x = mean, y = 0, xend = mean, yend = 50, linetype = signif, color = signif),
                     data = dplyr::filter(stat_rayleigh_moverel2light, !color == "dark"),
                     arrow = arrow(length = unit(0.05, "npc"), type = "closed", ends = "last"), size = 1) +
        color_scale_color +
        scale_x_continuous(
          "Movement relative to light (°)",
          limits = c(0, 360),
          breaks = c(0, 90, 180, 270),
          labels = c("Towards Light", "", "Away", "")
        ) +
        scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 10)) +
        scale_linetype_manual(values = c("FALSE" = NA, "TRUE" = "solid")) +
        facet_wrap(~color) +
        theme_void() +
        theme(panel.background = element_rect(color = "white"),
              panel.grid.major.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_line(color = "grey", size = 0.1),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 24),
              strip.text = element_text(size = 24),
              strip.background = element_rect(fill = "white", colour = "white"),
              legend.position = "none",
              plot.background = element_rect(color = "white"))
    })
    
    # MoveRel2Light stats
    output$moverel2light_stats <- renderPrint({
      
      stat_rayleigh_moverel2light <- plyr::ddply(d_moverel2light_binned, ~onoff + color, function(x) {
        r <- rayleigh.test(as.bearing(x$moverel2basin))
        if (r$p.value < 0.05) {
          s <- TRUE
        } else {
          s <- FALSE
        }
        mean <- mean.circular(as.bearing(x$moverel2basin))
        data.frame(id=unique(x$id), mean = mean, R = r$statistic, p.value = r$p.value, signif = s)
      })
      
      print(stat_rayleigh_moverel2light)
    })
    
    # Animate
    output$animation <- renderImage({
      p <- ggplot(indiv_tracks, aes(x, y, group = id, color = color)) +
        geom_point() +
        geom_point(aes(x_light, y_light, color = color), size = 10) +
        theme(legend.position = "none") +
        color_scale_color
      
      p1 <- p + 
        transition_states(states = time, transition_length = 1, state_length = 1, wrap = FALSE) + 
        shadow_wake(wake_length = 1, alpha = 0.7) +
        enter_fade() +
        exit_fade()
      
      anim <- animate(p1, fps = 3, renderer = gifski_renderer())
      
      anim_file <- tempfile(fileext = ".gif")
      anim_save(anim_file, animation = anim)
      
      list(src = anim_file, contentType = 'image/gif')
    }, deleteFile = TRUE)
    
  })
}
}

  # Run the application 
shinyApp(ui = ui, server = server)

