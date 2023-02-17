# to do
# check force curve is data in N?
# pressure curve and area curve


# =============================================================================

# Packages required
library(ggplot2)
library(rgl)
library(sp)
library(spatstat)
library(rgeos)
library(zoo)
library(ggmap)
library(scales)
library(grid)
library(gganimate)


# =============================================================================

#' Load emed data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @param rem_zeros Logical. If TRUE, remove columns and rows from edges which
#'   contain only zeros
#' @param flip Logical. Flips pressure data around its vertical axis.
#' @return Array. A 3D array covering each timepoint of the measurement. z
#'   dimension represents time

load_emed <- function(pressure_filepath, rem_zeros = TRUE, flip = FALSE) {
  # test inputs
  if(is.character(pressure_filepath) == FALSE)
    stop("filepath needs to be a character string")
  
  # Read unformated emed data
  pressure_raw <- readLines(pressure_filepath, warn = FALSE)
  
  # Determine dimensions of active sensor array
  ## determine position breaks, no of frames, y dimension, and frame type
  breaks <- grep("Page", pressure_raw)
  y_dim <- (breaks[2] - breaks[1]) - 13
  frame_type <- pressure_raw[breaks + 8]
  
  # identify and remove any summary frames
  MVP <- which(grepl("MVP", frame_type, fixed = TRUE))
  MPP <- which(grepl("MPP", frame_type, fixed = TRUE))
  if (length(MVP) == 0 && length(MPP) == 0) {
    breaks <- breaks
  } else {
    breaks <- breaks[-c(MVP, MPP)]
  }
  
  # determine x and z dimensions
  z_dim <- length(breaks)
  frame1 <- pressure_raw[(breaks[2] - 3):(breaks[1] + 11)]
  tc_frame1 <- textConnection(frame1)
  frame1 <- read.table(tc_frame1, sep = "\t")
  close(tc_frame1)
  x_dim <- ncol(frame1) - 2
  
  # make empty array to hold data
  pressure_array <- array(NA, dim = c(y_dim, x_dim, z_dim))
  
  # Make a list of individual frames
  for (i in seq_along(breaks)) {
    y <- pressure_raw[(breaks[i] + 11):(breaks[i] + 10 + y_dim)]
    tc_y <- textConnection(y)
    y <- read.table(tc_y, sep = "\t")
    y <- y[2:(x_dim + 1)]
    pressure_array[, , i] <- as.matrix(y)
    close(tc_y)
  }
  
  # if required, remove zero columns and rows
  dims <- dim(pressure_array)
  if (rem_zeros == TRUE) {
    sens_coords <- sensor_coords(pressure_array)
    min_x <- min(sens_coords$x_coord)
    min_x_n <- as.integer((min_x - 0.0025) / 0.005)
    max_x <- max(sens_coords$x_coord)
    max_x_n <- as.integer((max_x + 0.0025) / 0.005)
    min_y <- min(sens_coords$y_coord)
    min_y_n <- as.integer((min_y - 0.0025) / 0.005)
    max_y <- max(sens_coords$y_coord)
    max_y_n <- as.integer((max_y - 0.0025) / 0.005)
    pressure_array <- pressure_array[(dims[1] - max_y_n):(dims[1] - min_y_n), 
                                     min_x_n:max_x_n, ]
  }
  
  # if required, flip array around vertical axis
  dims <- dim(pressure_array)
  if (flip == TRUE) {pressure_array <- pressure_array[, c(dims[2]:1), ]}
  
  # return emed data
  return(pressure_array)
}
  

# =============================================================================

#' Load tekscan i-scan data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @param sensor_type String. Currently only "6900" supported
#' @param sensor_pad Integer. If sensor has multiple sensor pads choose the one
#'   of interest
#' @return Array. A 3D array covering each timepoint of the measurement. z
#'   dimension represents time
load_iscan <- function(pressure_filepath, sensor_type, sensor_pad) {
  # test inputs
  if(is.character(pressure_filepath) == FALSE)
    stop("filepath needs to be a character string")
  
  # Read unformated emed data
  pressure_raw <- readLines(pressure_filepath, warn = FALSE)
  
  # Determine dimensions of active sensor array
  ## determine position breaks, no of frames, y dimension, and frame type
  breaks <- grep("Frame", pressure_raw)
  z_dim <- pressure_raw[grep("END_FRAME", pressure_raw)] %>% 
    str_match_all("[0-9]+") %>% unlist %>% as.numeric
  if (sensor_type == "6900") {
    n_rows <- 11
    n_cols <- 11
    pad_pos <- data.frame(row_pos = c(3, 3, 16, 16), col_pos = c(3, 16, 3, 16))
  }
  
  # make empty array to hold data
  pressure_array <- array(NA, dim = c(n_rows, n_cols, z_dim))
  
  # Make a list of individual frames
  for (i in seq_along(breaks)) {
    y <- pressure_raw[(breaks[i] + pad_pos[sensor_pad, 1]):(breaks[i] + pad_pos[sensor_pad, 1] + n_rows - 1)]
    tc_y <- textConnection(y)
    y <- read.table(tc_y, sep = ",")
    y <- y[, pad_pos[sensor_pad, 2]:(pad_pos[sensor_pad, 2] + n_cols - 1)]
    pressure_array[, , i] <- as.matrix(y)
    close(tc_y)
  }
  
  # return i-Scan data
  return(pressure_array)
}



# =============================================================================

#' Generate force curve with option to plot
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_frames. Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param plot Logical. If TRUE also plots data as line curve
#' @return Numeric vector containing force values

force_curve <- function(pressure_frames) {
  # 
  force <- rep(NA, times = dim(pressure_frames)[3])
  for (i in 1:dim(pressure_frames)[3]) {
    force[i] <- sum(pressure_frames[, , i])
  }
  
  
  
  # return
  return(force)
}


pressure_curve

area_curve

# =============================================================================

#' Find footprint of pressure file
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_frames Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param value String. "max" = footprint of maximum sensors. "mean" = average
#'   value of sensors over time (usually for static analyses). "frame" = an individual frame 
#' @param frame Integer. 
#' @param plot Logical. Display pressure image
#' @return Matrix. Maximum or mean values for all sensors
#' @example

footprint <- function(pressure_frames, value = "max", frame, plot = FALSE) {
  if (value == "max") {
    mat <- apply(simplify2array(pressure_frames), 1:2, max)
  }
  if (value == "mean") {
    mat <- apply(simplify2array(pressure_frames), 1:2, mean)
  }
  if (value == "frame") {
    mat <- pressure_frames[,, frame]
  }
  
  # plot if requested
  if (plot == TRUE) {g <- plot_footprint(mat)}
  
  # return footprint
  return(mat)
}


# =============================================================================

#' Produce plot of pressure data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_frames Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param value String. "max" = footprint of maximum sensors. "mean" = average
#'   value of sensors over time (usually for static analyses). "frame" = an
#'   individual frame
#' @param frame Integer. 
#' @param interp Logical. If data should be interpolated to produce a clearer image
#' @param sens_x Numeric. Dimension of sensor in x direction, equivalent to
#'   column width in pressure matrix
#' @param sens_y Numeric. Dimension of sensor in y direction, equivalent to row
#'   height in pressure matrix
#' @param plot_COP Logical. If TRUE, overlay COP data on plot. Default = FALSE
#' @param plot_outline Logical. If TRUE, overlay convext hull outline on plot 
#' @param plot Logical. If TRUE, plot will be displayed
#' @return ggplot plot object

plot_pressure <- function(pressure_frames, value = "max", frame, interp = FALSE, 
                           sens_x = 0.005, sens_y = 0.005, plot_COP = FALSE, 
                           plot_outline = FALSE, plot = TRUE) {
  # if necessary, generate pressure matrix
  if (length(dim(pressure_frames)) == 3) {
    pressure_matrix <- footprint(pressure_frames, value, frame)
  } else {
    pressure_matrix <- pressure_frames
  }
  
  # interpolate if required
  if (interp == TRUE) {
    f_print2 <- matrix(rep(NA, (nrow(pressure_matrix) * ncol(pressure_matrix) * 5)),
                       nrow = nrow(pressure_matrix), ncol = ncol(pressure_matrix) * 5)
    for (i in 1:nrow(pressure_matrix)) {
      xa <- approx(pressure_matrix[i, ], n = 5 * ncol(pressure_matrix))
      f_print2[i, ] <- xa$y 
    }
    
    #Increase number of rows
    f_print3 <- matrix(rep(NA, (nrow(pressure_matrix) * 5 * ncol(pressure_matrix) * 5)),
                       nrow = nrow(pressure_matrix) * 5, ncol = ncol(pressure_matrix) * 5)
    for (i in 1:ncol(f_print2)) {
      xa <- approx(f_print2[ ,i], n = 5 * nrow(f_print2))
      f_print3[ ,i] <- xa$y 
    }
    pressure_matrix <- f_print3
  }
  
  # generate coordinates for each sensor
  dims <- dim(pressure_matrix)
  x_cor <- seq(from = (sens_x / 2), by = sens_x, length.out = dims[2])
  x_cor <- rep(x_cor, each = dims[1])
  y_cor <- seq(from = (sens_y / 2) + ((dims[1] - 1) * sens_y), 
               by = sens_y * -1, length.out = dims[1])
  y_cor <- rep(y_cor, times = dims[2])
  
  # combine with pressure values
  cor <- cbind(x_cor, y_cor, as.vector(pressure_matrix))
  cor <- as.data.frame(cor)
  colnames(cor) <- c("x", "y", "value")
  
  # add colours
  colour <- c()
  plot_cs <- list(cs_breaks = c(0, 15, 40, 60, 100, 150, 220, 300, 5000),
                  cs_cols = c("white", "grey", "blue", "light blue",
                              "green", "yellow", "red", "deeppink"))
  cols <- unlist(plot_cs[[2]])
  colour <- c()
  for (i in 1:(dims[1] * dims[2])) {
    for (j in seq_along(plot_cs[[1]])) {
      if(cor$value[i] >= plot_cs[[1]][j] & cor$value[i] < plot_cs[[1]][j + 1]) {
        colour = append(colour, j)
      }
    }
  }
  
  # combine with data frame
  cor <- cbind(cor, colour)
  
  ## plot 
  g <- ggplot()
  g <- g + geom_raster(data = cor, aes(x = x, y = y, fill = as.factor(colour)))
  g <- g + scale_fill_manual(values = cols)
  g <- g + scale_x_continuous(expand = c(0, 0))
  g <- g + scale_y_continuous(expand = c(0, 0))
  g <- g + coord_fixed()
  
  # add COP?
  if (plot_COP == TRUE) {
    cop_df <- gen_cop(pressure_frames)
    g <- g + geom_point(data = cop_df, aes(x = x_coord, y = y_coord))
  }
  
  # add outline
  if (plot_outline == TRUE) {
    ch_out <- footprint_outline(pressure_frames)
    g <- g + geom_path(data = ch_out, aes(x = x_coord, y = y_coord), 
                       colour = "black")
    g <- g + geom_point(data = ch_out, aes(x = x_coord, y = y_coord), 
                        colour = "purple")
  }
  
  # formatting
  g <- g + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                 axis.text.y = element_blank(), axis.ticks = element_blank(), 
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = "none")
  
  # display plot immediately if requested
  if (plot == TRUE) {print(g)}
  
  # return ggplot object
  return(g)
}

plot_sensor <- function(pressure_frames, value = "max", frame, interp = FALSE, 
                        sens_x = 0.005, sens_y = 0.005, plot_COP = FALSE, 
                        plot_outline = FALSE, plot = TRUE) {
  # if necessary, generate pressure matrix
  if (length(dim(pressure_frames)) == 3) {
    pressure_matrix <- footprint(pressure_frames, value, frame)
  } else {
    pressure_matrix <- pressure_frames
  }
  
  # interpolate if required
  if (interp == TRUE) {
    f_print2 <- matrix(rep(NA, (nrow(pressure_matrix) * ncol(pressure_matrix) * 5)),
                       nrow = nrow(pressure_matrix), ncol = ncol(pressure_matrix) * 5)
    for (i in 1:nrow(pressure_matrix)) {
      xa <- approx(pressure_matrix[i, ], n = 5 * ncol(pressure_matrix))
      f_print2[i, ] <- xa$y 
    }
    
    #Increase number of rows
    f_print3 <- matrix(rep(NA, (nrow(pressure_matrix) * 5 * ncol(pressure_matrix) * 5)),
                       nrow = nrow(pressure_matrix) * 5, ncol = ncol(pressure_matrix) * 5)
    for (i in 1:ncol(f_print2)) {
      xa <- approx(f_print2[ ,i], n = 5 * nrow(f_print2))
      f_print3[ ,i] <- xa$y 
    }
    pressure_matrix <- f_print3
  }
  
  # generate coordinates for each sensor
  dims <- dim(pressure_matrix)
  x_cor <- seq(from = (sens_x / 2), by = sens_x, length.out = dims[2])
  x_cor <- rep(x_cor, each = dims[1])
  y_cor <- seq(from = (sens_y / 2) + ((dims[1] - 1) * sens_y), 
               by = sens_y * -1, length.out = dims[1])
  y_cor <- rep(y_cor, times = dims[2])
  
  # combine with pressure values
  cor <- cbind(x_cor, y_cor, as.vector(pressure_matrix))
  cor <- as.data.frame(cor)
  colnames(cor) <- c("x", "y", "value")
  
  # add colours
  #colour <- c()
  #plot_cs <- list(cs_breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
  #                cs_cols = c("white", "grey", "blue", "light blue",
  #                            "green", "yellow", "red", "deeppink"))
  #cols <- unlist(plot_cs[[2]])
  #colour <- c()
  #for (i in 1:(dims[1] * dims[2])) {
  #  for (j in seq_along(plot_cs[[1]])) {
  #    if(cor$value[i] >= plot_cs[[1]][j] & cor$value[i] < plot_cs[[1]][j + 1]) {
  #      colour = append(colour, j)
  #    }
  #  }
  #}
  
  # combine with data frame
  #cor <- cbind(cor, colour)
  
  ## plot 
  g <- ggplot()
  g <- g + geom_raster(data = cor, aes(x = x, y = y, fill = value))
  #g <- g + scale_fill_manual(values = cols)
  g <- g + scale_fill_continuous(name = "Pressure (kPa)")
  g <- g + scale_x_continuous(expand = c(0, 0))
  g <- g + scale_y_continuous(expand = c(0, 0))
  g <- g + coord_fixed()
  
  # add COP?
  if (plot_COP == TRUE) {
    cop_df <- gen_cop(pressure_frames)
    g <- g + geom_point(data = cop_df, aes(x = x_coord, y = y_coord))
  }
  
  # formatting
  g <- g + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                 axis.text.y = element_blank(), axis.ticks = element_blank(), 
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())
  
  # display plot immediately if requested
  if (plot == TRUE) {print(g)}
  
  # return ggplot object
  return(g)
}


# =============================================================================

#' Produce animation of pressure footprint
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_frames Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param sens_x Numeric. Dimension of sensor in x direction, equivalent to
#'   column width in pressure matrix
#' @param sens_y Numeric. Dimension of sensor in y direction, equivalent to row
#'   height in pressure matrix
#' @param exp_name Name of export file(s)
#' @param path Location animation files are to be saved
#' @param type The file type for the animation. Currently only "png" supported.
#' @param res Resolution of images in dpi
#' @param width Width of animation
#' @param height Height of animation
#' @return A series of image files that can be combined into an avi file or
#'   similar. For example, ImageJ can be used to create the video

animate_footprint <- function(pressure_frames, sens_x = 0.005,
                              sens_y = 0.005, exp_name, path, type, res, width,
                              height) {
  for (i in 1:(dim(pressure_frames)[3])) {
    filename = paste0(path, "/", exp_name, i, ".png")
    g = plot_footprint(pressure_frames, value = "frame", frame = i, 
                       sens_x = sens_x, sens_y = sens_y)
      ggsave(filename, g, dpi = res, width = width, height = height,
             bg = "transparent")
  }
}

animate_footprintx <- function(pressure_frames, sens_x = 0.005, 
                              sens_y = 0.005) {
  ## helper function
  # generate dataframe with coords for each frame
  press_df <- function(pressure_matrix) {
    # generate coordinates for each sensor
    dims <- dim(pressure_matrix)
    x_cor <- seq(from = (sens_x / 2), by = sens_x, length.out = dims[2])
    x_cor <- rep(x_cor, each = dims[1])
    y_cor <- seq(from = (sens_y / 2) + ((dims[1] - 1) * sens_y), 
                 by = sens_y * -1, length.out = dims[1])
    y_cor <- rep(y_cor, times = dims[2])
    
    # combine with pressure values
    cor <- cbind(x_cor, y_cor, as.vector(pressure_matrix))
    cor <- as.data.frame(cor)
    colnames(cor) <- c("x", "y", "value")
    
    # add colours
    colour <- c()
    plot_cs <- list(cs_breaks = c(0, 15, 40, 60, 100, 150, 220, 300, 5000),
                    cs_cols = c("white", "grey", "blue", "light blue",
                                "green", "yellow", "red", "deeppink"))
    cols <- unlist(plot_cs[[2]])
    colour <- c()
    for (i in 1:(dims[1] * dims[2])) {
      for (j in seq_along(plot_cs[[1]])) {
        if(cor$value[i] >= plot_cs[[1]][j] & cor$value[i] < plot_cs[[1]][j + 1]) {
          colour = append(colour, j)
        }
      }
    }
    
    # combine with data frame
    cor <- cbind(cor, colour)  
  }

  # generate big df with frame number added
  df <- data.frame(x = as.numeric(character()),
                   y = as.numeric(character()),
                   value = as.numeric(character()),
                   colour = as.numeric(character()))
  for (i in 1:(dim(pressure_frames)[3])) {
    x <- press_df(pressure_frames[, , i])
    frame <- rep(i, times = nrow(x))
    x <- cbind(x, frame)
    df <- rbind(df, x)
  }
  
  # plot 
  ## plot 
  g <- ggplot(df, aes(x = x, y = y, fill = as.factor(colour), frame = frame))
  g <- g + geom_raster()
  g <- g + scale_fill_manual(values = cols)
  g <- g + scale_x_continuous(expand = c(0, 0))
  g <- g + scale_y_continuous(expand = c(0, 0))
  g <- g + coord_fixed()
  g <- g + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                 axis.text.y = element_blank(), axis.ticks = element_blank(), 
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = "none")
  gg_animate(g)
}



# =============================================================================

#' Interpolate pressure data. Useful for normalizing to stance phase, for example
#' @author Scott Telfer
#' \email{scott.telfer@gmail.com}
#' @param pressure_frames Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time.
#' @param interp_to Number of frames to interpolate to

pressure_interp <- function(pressure_frames, interp_to) {
  # make new array
  interp_array <- array(NA, dim = c(nrow(pressure_frames), 
                                    ncol(pressure_frames), interp_to))
  
  # interpolation function
  approxP <- function(x, interp_to) {
    y <- approx(x, n = as.numeric(interp_to))
    y$x <- NULL
    y <- unname(unlist(y))
    return(y)
  }
  
  # interpolate array
  array_dim <- dim(pressure_frames)
  for (i in 1:array_dim[1]) {
    for (j in 1:array_dim[2]) {
      interp_array[i, j, ] <- approxP(pressure_frames[i, j, ], interp_to)
    }
  }
  
  # return interpolated array
  return(interp_array)
}


# =============================================================================

#' Manually define orthoganal area
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_frames Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param image_value String. "max" = footprint of maximum sensors. "mean"
#'   average value of sensors over time (usually for static analyses)
#' @return Array. A 3D array covering each timepoint of the measurement for the
#'   selected region. z dimension represents time

select_region <- function(pressure_frames, image_value = "max", 
                          sens_x = 0.005, sens_y = 0.005) {
  # plot footprint
  g <- plot_footprint(pressure_frames)
  print(g)
  
  # interactively select area
  message("Select footprint")
  mask <- gglocator(4)
  mask <- data.frame(x = mask[, 1], y = mask[, 2])
  mask <- mask[c(1:nrow(mask), 1), ]
  
  # which sensors fall in mask
  dims <- dim(pressure_frames)
  x_cor <- seq(from = (sens_x / 2), by = sens_x, length.out = dims[2])
  x_cor <- rep(x_cor, each = dims[1])
  y_cor <- seq(from = (sens_y / 2) + ((dims[1] - 1) * sens_y), 
               by = sens_y * -1, length.out = dims[1])
  y_cor <- rep(y_cor, times = dims[2])
  sensor_def <- point.in.polygon(point.x = x_cor, point.y = y_cor, 
                                 pol.x = mask$x, pol.y = mask$y)
  masked <- which(sensor_def == 1)
  
  # define orthogonal area
  x_max <- ceiling(max(x_cor[masked]) / 0.005)
  x_min <- floor(min(x_cor[masked])  / 0.005)
  y_max <- ceiling(max(y_cor[masked]) / 0.005)
  y_max <- dims[1] - y_max
  y_min <- floor(min(y_cor[masked]) / 0.005)
  y_min <- dims[1] - y_min
  
  # define reduced array
  pressure_array_masked <- pressure_frames[y_max:y_min, x_min:x_max, ]
  
  # return pressure frames for selected area
  return(pressure_array_masked)
}


# =============================================================================

#' Generate center of pressure coordinates
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_frames Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param sens_x_size Numeric. The size of each pressure sensor (x direction)
#' @param sens_y_size Numeric. The size of each pressure sensor (y direction)
#' @param interpol Numeric. Include if the COP is to be interpolated to a set
#'   number of timepoints
#' @return Data frame with x and y coordinates of COP throughout trial

cop <- function(pressure_frames, sens_x = 0.005, sens_y = 0.005,
                    interpol = FALSE) {
  # array dimensions
  x <- pressure_frames
  dims <- dim(x)
  
  # loading totals by column
  col_total <- data.frame(matrix(NA, nrow = dims[2], ncol = dims[3]))
  for (i in 1:dims[3]) {col_total[, i] <- colSums(x[, , i])}
  
  # Loading totals by row
  row_total <- data.frame(matrix(NA, nrow = dims[1], ncol = dims[3]))
  for (i in 1:dims[3]) {row_total[, i] <- rowSums(x[, , i])}
  
  # Sensor spacing in x direction
  sens_spacing_x <- seq(from = sens_x / 2, by = sens_x, length.out = dims[2])
  
  # Sensor spacing in y direction
  sens_spacing_y <- rev(seq(from = sens_y / 2, by = sens_y, 
                            length.out = dims[1]))
  
  # COP coordinates in x direction
  x_coord <- c()
  for (i in 1:dims[3]) {
    p_total <- sum(col_total[, i])
    x_coord[i] <- (sum(sens_spacing_x * col_total[, i])) / p_total 
  }
  
  # COP coordinates in y direction
  y_coord <- c()
  for (i in 1:dims[3]) {
    p_total <- sum(row_total[, i])
    y_coord[i] <- (sum(sens_spacing_y * row_total[, i])) / p_total
  }
  
  # interpolate if required
  if (hasArg(interpol) == TRUE) {
    x_coord <- approx(x_coord, n = interpol)
    x_coord <- x_coord$y
    y_coord <- approx(y_coord, n = interpol)
    y_coord <- y_coord$y
  }
  
  # combine coordinates into dataframe
  COP_df <- data.frame(x_coord, y_coord)
  
  # return COP coordinates
  return(COP_df)
}


# =============================================================================

#' Determine Center of Pressure Excursion Index (CPEI)
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_frames Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param side String. "right" or "left". Required for automatic detection of
#'   points
#' @param plot_result Logical. Plots pressure image with COP and CPEI overlaid
#' @return Numeric. CPEI value

cpei <- function(pressure_frames, side, plot_result = FALSE) {
  ## geometry helper functions 
  # Finds the intersection point of a line from a point perpendicular to 
  # another line
  perp_intersect <- function(line, point) {
    x1 <- line[1, 1]
    x2 <- line[2, 1]
    x3 <- point[1, 1]
    y1 <- line[1, 2]
    y2 <- line[2, 2]
    y3 <- point[1, 2]
    k <- ((y2-y1) * (x3-x1) - (x2-x1) * (y3-y1)) / ((y2-y1)^2 + (x2-x1)^2)
    x_coord <- x3 - k * (y2-y1)
    y_coord <- y3 + k * (x2-x1)
    
    # return
    return(data.frame(x_coord, y_coord))
  }
  
  # Finds intersection point of two lines defined by points
  line_intersect <- function(line1, line2) {
    line2[1, 1] <- line2[1, 1] + 0.000001 # allows vertical lines to be solved 
    m1 <- lm(line1[, 2] ~ line1[, 1])
    m2 <- lm(line2[, 2] ~ line2[, 1])
    cm <- rbind(coef(m1), coef(m2)) # coefficient matrix
    intersect <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
    
    # return
    return(intersect)
  }
  
  # Finds point closest to line
  closest_point <- function(points, line) {
    distances <- c()
    for (i in 1:nrow(points)) {
      perp_point <- perp_intersect(line, points[i, ])
      distance <- as.numeric(dist(rbind(perp_point, points[i, ])))
      distances <- c(distances, distance)
    }
    min_dist <- which.min(distances)
    return(min_dist)
  }
  
  # manual select function
  manually_select <- function(n_points, mess) {
    message(mess)
    x <- gglocator(n_points)
    colnames(x) <- c("x_coord", "y_coord")
    return(x)
  }
  
  ## Automatically identify all required points
  # outline of main part of foot
  toe_cut <- ceiling(nrow(pressure_frames) / 5)
  pressure_frames_trimmed <- pressure_frames[toe_cut:nrow(pressure_frames),
                                             , ]
  fp_out <- footprint_outline(pressure_frames_trimmed)
  
  # find the 2 sets of points with the largest dstance between them
  dis <- as.matrix(dist(fp_out))
  dis <- dis[row(dis) == (col(dis) - 1)]
  max_2 <- order(dis)[(length(dis) - 1):length(dis)]
  
  # identify medial and lateral borders
  bor_1 <- (fp_out[max_2[1], 1] + fp_out[max_2[1] + 1, 1]) / 2
  bor_2 <- (fp_out[max_2[2], 1] + fp_out[max_2[2] + 1, 1]) / 2 
  if (side == "right" & bor_2 > bor_1) {
    m_bor = fp_out[max_2[1]:(max_2[1] + 1), ]
    l_bor = fp_out[max_2[2]:(max_2[2] + 1), ]
  }
  if (side == "right" & bor_1 > bor_2) {
    m_bor = fp_out[max_2[2]:(max_2[2] + 1), ]
    l_bor = fp_out[max_2[1]:(max_2[1] + 1), ]
  }
  if (side == "left" & bor_2 > bor_1) {
    m_bor = fp_out[max_2[2]:(max_2[2] + 1), ]
    l_bor = fp_out[max_2[1]:(max_2[1] + 1), ]
  }
  if (side == "left" & bor_1 > bor_2) {
    m_bor = fp_out[max_2[1]:(max_2[1] + 1), ]
    l_bor = fp_out[max_2[2]:(max_2[2] + 1), ]
  }
  
  # find most proximal heel sensor
  sc <- sensor_coords(pressure_frames)
  sc_min <- which(sc$y_coord == min(sc$y_coord))
  sc_dist <- which.min(dist(rbind(m_bor[1, ], 
                                  sc[sc_min, 2:3]))[1:length(sc_min)])
  heel_coord <- sc[sc_min[sc_dist], 2:3]
  
  # find most distal toe sensor
  sc_min <- which(sc$y_coord == max(sc$y_coord))
  sc_dist <- which.min(dist(rbind(m_bor[1, ],
                                  sc[sc_min, 2:3]))[1:length(sc_min)])
  toe_coord <- sc[sc_min[sc_dist], 2:3]
  
  # Find medial start and end points of COP
  cop <- gen_cop(pressure_frames)
  cop_start <- cop[1:15, ] 
  start_point_n <- closest_point(cop_start, m_bor)
  start_point <- cop_start[start_point_n, ]
  cop_end <- cop[(nrow(cop) - 15):nrow(cop), ]
  end_point_n <- closest_point(cop_end, m_bor)
  end_point <- cop_end[end_point_n, ] 
  
  ## check if automatic identification worked
  # plot max footprint with additional points
  g <- plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)
  g <- g + geom_point(data = heel_coord, aes (x = x_coord, y = y_coord), shape = 6)
  g <- g + geom_point(data = toe_coord, aes (x = x_coord, y = y_coord), shape = 6)
  g <- g + geom_point(data = start_point, aes (x = x_coord, y = y_coord), shape = 2)
  g <- g + geom_point(data = end_point, aes (x = x_coord, y = y_coord), shape = 2)
  g <- g + geom_line(data = m_bor, aes(x = x_coord, y = y_coord), colour = "purple")
  g <- g + geom_line(data = l_bor, aes(x = x_coord, y = y_coord), colour = "orange")
  print(g)
  auto_worked <- readline("Have points been correctly identified? c: manually select cop; a: manually select all")
  
  ## if automatic identification failed, redo manually
  if (auto_worked == "a") {
    # plot footprint
    g <- plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)
    
    # select points
    m_bor <- manually_select(2, "select medial border: proximal point first, then distal")
    l_bor <- manually_select(2, "select lateral border: proximal point first, then distal")
    heel_coord <- manually_select(1, "select most proximal point of heel") 
    toe_coord <- manually_select(1, "select most distal point of toes")
    cop_start <- manually_select(1, "select the most medial point near the start of the COP")
    cop_end <- manually_select(1, "select the most medial point near the end of the COP")
  }
  
  if(auto_worked == "c") {
    # plot footprint
    g <- plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)
    
    # select points
    start_point <- manually_select(1, "select the most medial point near the start of the COP")
    end_point <- manually_select(1, "select the most medial point near the end of the COP")
  }
  
  ## Calculate CPEI
  # determine foot length
  heel_intersect <- perp_intersect(m_bor, heel_coord)
  toe_intersect <- perp_intersect(m_bor, toe_coord)
  foot_length <- dist(rbind(heel_intersect, toe_intersect))
  
  # determine distal tertile
  dist_tert_point <- heel_intersect + ((toe_intersect - heel_intersect) * 0.67)
  dist_tert_point2 <- heel_coord + ((toe_intersect - heel_intersect) * 0.67)
  tert_line <- as.data.frame(rbind(dist_tert_point, dist_tert_point2))

  # intersects with borders
  med_intersect <- line_intersect(tert_line, m_bor)
  lat_intersect <- line_intersect(tert_line, l_bor)
  ff_line <- as.data.frame(rbind(med_intersect, lat_intersect))
  colnames(ff_line) <- c("x_coord", "y_coord")
  
  # forefoot width
  ff_width <- as.numeric(dist(ff_line))
  
  # find intersection between ff width line and COP
  m_ff <- lm(ff_line[, 2] ~ ff_line[, 1])
  m_ff_cop <- (coef(m_ff)[2] * cop$x_coord) + coef(m_ff)[1]
  cross_point <- min(which(cop$y_coord > m_ff_cop))
  ff_cop_int <- line_intersect(ff_line, rbind(cop[cross_point, ], 
                                              cop[cross_point -1, ]))
  
  # find intersection between ff width line and CPEI construction line 
  cpei_con <- as.data.frame(rbind(start_point, end_point))
  colnames(cpei_con) <- c("x_coord", "y_coord")
  ff_con_int <- line_intersect(ff_line, cpei_con)
  
  # calculate CPE
  CPE_df <- as.data.frame(rbind(ff_cop_int, ff_con_int))
  colnames(CPE_df) <- c("x_coord", "y_coord")
  CPE <- as.numeric(dist(CPE_df))
  
  # calculate CPEI
  CPEI <- (CPE / ff_width) * 100
  
  # plot
  if (plot_result == TRUE) {
    g <- plot_footprint(pressure_frames, plot_COP = TRUE)
    g <- g + geom_line(data = ff_line, aes(x_coord, y_coord), 
                       linetype = "dashed")
    g <- g + geom_line(data = CPE_df, aes(x_coord, y_coord),
                       size = 2)
    g <- g + geom_line(data = cpei_con, aes(x_coord, y_coord), colour = "blue", 
                       alpha = 0.8)
    print(g)
  }
  
  # return CPEI
  return(CPEI)
}


# =============================================================================

#' Determine outline (convex hull) of footprint
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_frames Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param sens_x Numeric. Dimension of sensor in x direction, equivalent to
#'   column width in pressure matrix
#' @param sens_y Numeric. Dimension of sensor in y direction, equivalent to row
#'   height in pressure matrix
#' @return Data frame. x and y cordinates of convex hull outline

footprint_outline <- function(pressure_frames, sens_x = 0.005, 
                              sens_y = 0.005) {
  # determine active sensor coordinates
  sens_coords <- sensor_coords(pressure_frames, sens_x, sens_y)
  
  # calculate convex hull
  con_hull <- chull(sens_coords[, 2:3])
  con_hull <- sens_coords[con_hull, 2:3]
  con_hull <- rbind(con_hull, con_hull[1, ])
  rownames(con_hull) <- c()
  
  # return convex hull coordinates
  return(con_hull)
}


# =============================================================================

#' Coordinates of active sensors
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_frames Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param sens_x Numeric. Dimension of sensor in x direction, equivalent to
#'   column width in pressure matrix
#' @param sens_y Numeric. Dimension of sensor in y direction, equivalent to row
#'   height in pressure matrix
#' @param threshold Numeric. Level required to consider sensor active
#' @return Data frame. Includes number, x and y cordinates of sensors

sensor_coords <- function(pressure_frames, sens_x = 0.005, sens_y = 0.005, 
                          threshold = 5) {
  # max pressure footprint
  max_fp <- footprint(pressure_frames)
  
  # data frame with active sensors as coordinates
  P <- c(max_fp)
  x_cor <- seq(from = sens_x / 2, by = sens_x, length.out = ncol(max_fp))
  x_cor <- rep(x_cor, each = nrow(max_fp))
  y_cor <- seq(from = (sens_y / 2) + ((nrow(max_fp) - 1) * sens_y),
               by = (-1 * sens_y), length.out = nrow(max_fp))
  y_cor <- rep(y_cor, times = ncol(max_fp))
  coords <- data.frame(Sens_No = 1:length(x_cor), x_coord = x_cor, 
                       y_coord = y_cor, P = P)
  coords <- coords[which(P >= 5), ]
  coords <- coords[, 1:3]
  
  # return sensor coordinates
  return(coords)
}


# =============================================================================

#' gglocator
#' @param n Integer. Number of points to select
#' @return Data frame. x and y coordinates of selected point

gglocator <- function(n = 1, message = FALSE, xexpand = c(.0, 0), 
                      yexpand = c(.0, 0), mercator = TRUE) {
  
  if(n > 1){
    df <- NULL
    for(k in 1:n){
      df <- rbind(df, gglocator(message = message,
                                xexpand = xexpand, yexpand = yexpand, 
                                mercator = mercator))
    }
    return(df)
  }
  
  object <- last_plot()
  if(is.null(object)){
    stop("no plots available")
  }
  
  # find the correct viewport for the npc coordinates
  x <- unlist(current.vpTree())
  x <- unname(x[grep("\\.name$", names(x))])
  x <- grep("panel", x, fixed = TRUE, value = TRUE)
  n_panels <- length(x)
  if(n_panels == 0){
    stop("ggmap plot not detected in current device")
  }
  if(n_panels > 1){
    x <- x[1]
    warning(gettextf("multiple plots detected, choosing one (\"%s\")",
                     x), domain = NA)
  }
  previous_viewport <- current.vpPath()
  seekViewport(x, recording = FALSE)
  
  # when exiting function, return to previous position in viewport tree
  on.exit(upViewport(0, recording = FALSE))
  if(!is.null(previous_viewport)){
    on.exit(downViewport(previous_viewport, strict = TRUE, recording = FALSE),
            add = TRUE)
  }
  
  # get the position relative to that viewport
  loc <-  as.numeric(grid.locator("npc"))
  
  # scale the position to the plot
  
  # get the x.range and y.range from ggplot
  plot_info <- ggplot_build(object)
  if("layout" %in% names(plot_info)){
    ranges <- plot_info$layout$panel_ranges[[1]]
  } else{
    ranges <- plot_info$panel$ranges[[1]]
  }
  xrng <- ranges$x.range
  yrng <- ranges$y.range
  
  xrng <- expand_range(range = xrng, mul = xexpand[1], add = xexpand[2])
  yrng <- expand_range(range = yrng, mul = yexpand[1], add = yexpand[2])
  
  # format and return
  point <- data.frame(xrng[1] + loc[1]*diff(xrng),
                      yrng[1] + loc[2]*diff(yrng))
  if(isTRUE(mercator)){
    yrng2 <- LonLat2XY(0, yrng, zoom = 0, ypix = 256)$y
    point[[2]] <- XY2LonLat(y = yrng2[1] + loc[2] * diff(yrng2),
                            x = 0, X = 0, Y = 0, zoom = 0, ypix = 256)[[2]]
  }
  names(point) <- with(object, c(deparse(mapping$x), deparse(mapping$y)))
  point
}


# =============================================================================

#file_path <- "C:/Users/telfe/Dropbox/My_Projects/Twinfoot/data/TWINFOOT001A/right/798.lst"
#file_path <- "C:/Users/telfe/Dropbox/My_Projects/Twinfoot/data/TWINFOOT001A/static.lst"
#pressure_frames <- load_emed(file_path, flip = TRUE)
#plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)
#cpei(pressure_frames, side = "right", plot = TRUE)

