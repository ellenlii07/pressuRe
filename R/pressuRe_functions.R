# to do
# add more input tests to throw errors
# pti to pressure curve function
# option to choose different colors in plot_pressure
# fscan import function
# more output variables need to be added to  mask analysis
# edit mask functionality using identity
# automask to work on different sensors (currently just emed) 0.0025s in sensor coords
# allow for multiple masking schemes
# change footprint function name to more generic name?
# create_mask and edit_mask need to work interactively how to check for cran
# select steps > output data frame
# Do we have all pedar insoles? Double check areas
# force curve (and others) to work with insole data
# filepath to fscan example
# can load-iscan be made more generic?
# can cop use sensor_coords to get coordinates?
# add support for pliance
# edit mask needs to be written
# global pressure_import function (leave for V2)
# document pedar insole size data

# data list:
## Array. pressure data
## Character string. data type (usually collection system, e.g. emed)
## Numeric. sens_size. sensor size
## Numeric. Single number time between
## List. Mask list
## events (for example, to define start/end of individual steps for insole data)


# =============================================================================

# Packages required
#library(rgeos)
#library(zoo)


# =============================================================================

#' @title Load emed data
#' @description Imports and formats .lst files collected on emed system and
#'    exported from Novel software
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. Character defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#' @examples
#' pressure_data = load_emed("inst/extdata/emed_test.lst")
#' @importFrom stringr str_extract_all str_detect
#' @export

load_emed <- function(pressure_filepath) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  if (str_split(basename(pressure_filepath), "\\.")[[1]][2] != "lst")
    stop("incorrect file extension, expected .lst")

  # Read unformatted emed data
  pressure_raw <- readLines(pressure_filepath, warn = FALSE)

  # get sensor size
  sens_size_ln <- which(grepl("Sensor size", pressure_raw))[1]
  sens_size <- str_extract_all(pressure_raw[sens_size_ln], "\\d+\\.\\d+")
  sens_size <- as.numeric(unlist(sens_size))
  if (str_detect(pressure_raw[sens_size_ln], "cm") == TRUE) {
    sens_size <- sens_size * 0.01
  }

  # get capture frequency
  time_line <- which(grepl("Time", pressure_raw))[1]
  time_ln <- str_split(pressure_raw[time_line], "picture:")[[1]][2]
  time <- as.numeric(unlist(str_extract_all(time_ln, "\\d+\\.\\d+")))
  if (str_detect(time_ln, "ms") == TRUE) {time <- time / 1000}

  # determine position breaks
  breaks <- grep("Page", pressure_raw)

  # identify and remove any summary frames
  frame_type <- pressure_raw[breaks + 8]
  MVP <- which(grepl("MVP", frame_type, fixed = TRUE))
  MPP <- which(grepl("MPP", frame_type, fixed = TRUE))
  if (length(MVP) > 0 | length(MPP) > 0) {
    breaks <- breaks[1:(min(c(MVP, MPP)) - 1)]
  }

  # get blank lines
  ends <- which(pressure_raw == "\x0C")

  # how many frames in each measurement
  nfs <- sum(str_detect(pressure_raw[breaks + 8], "Pict-No\\.: 1 "))

  # how many measurements
  n_meas <- floor(length(breaks)/nfs)

  # empty array
  pressure_array <- array(0, dim = c(95, 64, n_meas))
  for (i in 1:n_meas) {
    for (j in 1:nfs) {
      # start line
      str <- (nfs * i) - nfs + j

      # load as table
      y <- pressure_raw[(breaks[str] + 10):(ends[which(ends > breaks[str])[1]] - 2)]
      num_col <- unlist(str_split(y[1], "\\s+"))
      wids <- rep(8, times = length(num_col))
      if (nfs > 1) {
        z <- read.fwf(textConnection(y), widths = wids)
      } else {
        z <- read.table(textConnection(y), sep = "\t")
      }

      colnames(z) <- unname(z[1, ])
      z <- z[2:nrow(z), ]

      # remove zeros
      z[is.na(z)] <- 0

      # remove force column
      if (colnames(z)[ncol(z)] == "Force") {z <- z[, 1:(ncol(z) - 1)]}

      # remove force row
      if (z[nrow(z), 1] == "Force") {z <- z[1:(nrow(z) - 1),]}

      # column numbers
      cn <- as.numeric(colnames(z)[2:length(z)])

      # row numbers
      rn <- unname(unlist((z[, 1])))
      z <- as.matrix(z[, 2:ncol(z)])
      z <- matrix(as.numeric(z), ncol = ncol(z))

      # add to array
      pressure_array[rn, cn, i] <- z
    }
  }

  # remove zero columns and rows
  ## make max footprint
  fp <- apply(simplify2array(pressure_array), 1:2, max)

  ## rows
  rsums <- rowSums(fp)
  minr <- min(which(rsums > 0))
  maxr <- max(which(rsums > 0))

  ## columns
  csums <- colSums(fp)
  minc <- min(which(csums > 0))
  maxc <- max(which(csums > 0))

  ## update pressure array
  pressure_array <- pressure_array[minr:maxr, minc:maxc,]

  # return formatted emed data
  return(list(pressure_array = pressure_array, pressure_system = "emed",
              sens_size = sens_size,
              time = time, masks = NULL, events = NULL))
}


# =============================================================================

#' @title Load pedar data
#' @description Imports and formats .asc files collected on pedar system and
#'    exported from Novel software
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. Character defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#' @examples
#' pressure_data = load_pedar("inst/extdata/pedar_example.asc")
#' @importFrom stringr str_split str_trim
#' @export

load_pedar <- function(pressure_filepath) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  if (str_split(basename(pressure_filepath), "\\.")[[1]][2] != "asc")
    stop("incorrect file extension, expected .asc")

  # measurement data
  ## header
  header <- readLines(pressure_filepath, n = 3)

  # Identify insole
  insole_type_ln <- str_split(header[2], "\\:")[[1]][2]
  insole_type <- str_split(insole_type_ln, "\\-")[[1]][1]
  insole_type <- str_trim(insole_type)[[1]][1]

  # Get time between measurements
  time_ln <- str_split(header[3], "\\[Hz\\]\\:")[[1]][2]
  time <- 1 / as.numeric(time_ln)

  # import data
  pressure_data <- as.data.frame(read.table(pressure_filepath,
                                            sep = "", skip = 10,
                                            header = FALSE))

  # make into 2 (left and right insole) x 99 (sensors) x time array
  pressure_array <- array(0, dim = c(2, 99, nrow(pressure_data)))
  for (i in 1:nrow(pressure_data)) {
    pressure_array[1, , i] <- unlist(unname(pressure_data[i, c(101:199)]))
    pressure_array[2, , i] <- unlist(unname(pressure_data[i, c(2:100)]))
  }

  # return
  return(list(pressure_array = pressure_array, pressure_system = "pedar",
              sens_size = insole_type, time = time, masks = NULL,
              events = NULL))
}


# =============================================================================

#' @title Load F-scan data
#' @description Imports and formats .asc files collected on f-scan system and
#'    exported from Tekscan software
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. Character defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#'  @examples
#'  pressure_data <- load_fscan("inst/extdata/)
#'  @importFrom
#'  @export

load_fscan <- function(pressure_filepath) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  if (str_split(basename(pressure_filepath), "\\.")[[1]][2] != "asf")
    stop("incorrect file extension, expected .asf")
}


# =============================================================================

#' @title Load tekscan i-scan data
#' @description Imports and formats .asc files collected on i-scan system and
#'    exported from Tekscan software
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @param sensor_type String. Currently only "6900" supported
#' @param sensor_pad Integer. If sensor has multiple sensor pads choose the one
#'   of interest
#' @return Array. A 3D array covering each timepoint of the measurement. z
#'   dimension represents time
#' @examples
#' pressure_data <- load_iscan()
#' @importFrom stringr str_match_all
#' @export

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

#' @title Interpolate pressure data
#' @description Interpolates data over time. Useful for normalizing to stance
#' phase, for example
#' @author Scott Telfer
#' \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item should be a 3D array covering each
#' timepoint of the measurement. z dimension represents time.
#' @param interp_to Integer. Number of frames to interpolate to
#' @return
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. Character defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#' @examples
#' pressure_interp(pressure_data, interp_to = 101)
#' @export

pressure_interp <- function(pressure_data, interp_to) {
  # check inputs
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")
  if (is.numeric(interp_to) & interp_to < 2)
    stop("interp_to must be a numeric value greater than 2")
  interp_to <- round(interp_to)

  # make new empty array
  dims <- dim(pressure_data[[1]])
  interp_array <- array(NA, dim = c(dims[1], dims[2], interp_to))

  # interpolate array
  pressure_array <- pressure_data[[1]]
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      interp_array[i, j, ] <- approxP(pressure_array[i, j, ], interp_to)
    }
  }

  # timing
  time_seq <- seq(0, by = pressure_data[[3]], length.out = dims[3])
  time_seq_int <- approxP(time_seq, interp_to)
  time_sample_int <- time_seq_int[2] - time_seq_int[1]

  # update pressure data
  pressure_data[[1]] <- interp_array
  pressure_data[[4]] <- time_sample_int

  # return interpolated pressure data
  return(pressure_data)
}


# =============================================================================

#' @title Select steps
#' @description Select steps, usually from insole data, and format for analysis
#' @param pressure_data. List
#' @param threshold_R Numeric. Threshold force to define start and end of step
#' @param threshold_L Numeric. Threshold force to define start and end of step
#' @param min_frames Numeric. Minimum number of frames that need to be in step
#' @param steps_Rn Numeric. Target number of steps for right foot. User will be
#' asked to keep selected steps until this target is reached or they run out of
#' candidate steps
#' @param steps_Ln Numeric. Target number of steps for left foot. User will be
#' asked to keep selected steps until this target is reached or they run out of
#' candidate steps
#' @param skip Numeric. Usually the first few steps of a trial are accelerating
#' and not representative of steady state walking so this removes them
#' @return
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. Character defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   }
#' @examples
#' select_steps(pressure_data)
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab ggtitle
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @export

select_steps <- function (pressure_data, threshold_R = 20,
                          threshold_L = 20, min_frames = 10,
                          steps_Rn = 5, steps_Ln = 5, skip = 2) {
  # check this is pedar (or other suitable) data
  if (!(pressure_data[[2]] == "pedar" || pressure_data[[2]] == "fscan"))
    stop("data should be from pedar or f-scan")

  # if pedar, get insole size
  pedar_insole_type <- pressure_data[[3]]

  # Read in required pedar sensor areas (imported in mm^2 but adjust for kPa)
  load("data/pedar_insole_areas.rda")
  pedarSensorAreas <- as.vector(pedar_insole_areas[[pedar_insole_type]] * 0.001)

  # Convert pressure data to force data (R & L)
  dims <- dim(pressure_data[[1]])
  force_R_mat <- matrix(pressure_data[[1]][1, , ],
                        dims[2], dims[3]) * pedarSensorAreas
  force_L_mat <- matrix(as.vector(pressure_data[[1]][2, , ]),
                        dims[2], dims[3]) * pedarSensorAreas
  force_R <- colSums(force_R_mat)
  force_L <- colSums(force_L_mat)

  # Adjust thresholds to avoid errors
  threshold_R <- threshold_R + 0.01
  threshold_L <- threshold_L + 0.01

  # Get events
  FS_events_R <- which(force_R[-length(force_R)] < threshold_R &
                         force_R[-1] > threshold_R) + 1
  FO_events_R <- which(force_R[-length(force_R)] > threshold_R &
                         force_R[-1] < threshold_R) + 1
  FS_events_L <- which(force_L[-length(force_L)] < threshold_L &
                         force_L[-1] > threshold_L) + 1
  FO_events_L <- which(force_L[-length(force_L)] > threshold_L &
                         force_L[-1] < threshold_L) + 1

  # Create steps
  ## Adjust events to ensure stance phase is used
  if (FO_events_R[1] < FS_events_R[1]) {
    FO_events_R <- FO_events_R[-1]
  }
  if ((length(FS_events_R)) > (length(FO_events_R))) {
    FS_events_R <- FS_events_R[-length(FS_events_R)]
  }
  if (FO_events_L[1] < FS_events_L[1]) {
    FO_events_L <- FO_events_L[-1]
  }
  if (length(FS_events_L) > length(FO_events_L)) {
    FS_events_L <- FS_events_L[-length(FS_events_L)]
  }

  # Make steps into data frame
  ## right
  df_R <- data.frame(step = integer(), frame = integer(), force = double())
  for (i in 1:length(FS_events_R)) {
    if (FO_events_R[i] - FS_events_R[i] > min_frames) {
      force_step <- force_R[FS_events_R[i]:FO_events_R[i]]
      step <- data.frame(step = rep(i, length_out = length(force_step)),
                         frame = c(1:length(force_step)),
                         force = force_step)
      df_R <- rbind(df_R, step)
    }
  }

  ## left
  df_L <- data.frame(step = integer(), frame = integer(), force = double())
  for (i in 1:length(FS_events_L)) {
    if (FO_events_L[i] - FS_events_L[i] > min_frames) {
      force_step <- force_L[FS_events_L[i]:FO_events_L[i]]
      step <- data.frame(step = rep(i, length_out = length(force_step)),
                         frame = c(1:length(force_step)),
                         force = force_step)
      df_L <- rbind(df_L, step)
    }
  }

  # Approve or discard steps
  ## right
  include_stps_R <- c()
  for (stp in 1:length(FS_events_R)) {
    # highlighted step df
    df1 <- df_R %>% filter(step %in% stp)

    # plot
    g <- ggplot()
    g <- g + geom_line(data = df_R, aes(x = frame, y = force,
                                 group = step),
                       linewidth = 1.5, color = "grey")
    g <- g + geom_line(data = df1, aes(x = frame, y = force),
                       linewidth = 1.5, color = "red")
    g <- g + xlab("Frame no") + ylab("Force (N)")
    g <- g + ggtitle(paste0("Right step ", stp))
    print(g)
    Sys.sleep(2.5)

    # get user to approve or reject step
    resp <- menu(c("Y", "N"),
                 title = "Do you want to keep (Y) or discard (N) this step?")
    include_stps_R <- c(include_stps_R, resp)

    # check if number of steps has been reached
    if (sum(include_stps_R == 1) >= steps_Rn) {break}
  }

  # Approve or discard steps
  ## left
  include_stps_L <- c()
  for (stp in 1:length(FS_events_L)) {
    # highlighted step df
    df1 <- df_L %>% filter(step %in% stp)

    # plot
    g <- ggplot()
    g <- g + geom_line(data = df_L, aes(x = frame, y = force,
                                        group = step),
                       linewidth = 1.5, color = "grey")
    g <- g + geom_line(data = df1, aes(x = frame, y = force),
                       linewidth = 1.5, color = "red")
    g <- g + xlab("Frame no") + ylab("Force (N)")
    g <- g + ggtitle(paste0("Left step ", stp))
    print(g)
    Sys.sleep(2.5)

    # get user to approve or reject step
    resp <- menu(c("Y", "N"),
                 title = "Do you want to keep (Y) or discard (N) this step?")
    include_stps_L <- c(include_stps_R, resp)

    # check if number of steps has been reached
    if (sum(include_stps_R == 1) >= steps_Rn) {break}
  }

  # make events df
  event_df <- data.frame(side =  c(rep("RIGHT", length.out = length(include_stps_R)),
                                   rep("LEFT", length.out = length(include_stps_L))),
                         FON = c(FS_events_R[which(include_stps_R == 1)],
                                 FS_events_L[which(include_stps_L == 1)]),
                         FOFF = c(FO_events_R[which(include_stps_R == 1)],
                                  FO_events_L[which(include_stps_L == 1)]))

  # add events to pressure data
  pressure_data[[6]] <- event_df

  # Return
  return(pressure_data)
}


# =============================================================================

#' @title Detect foot side
#' @description Detects which foot plantar pressure data is from (left or
#' right), usually would only be needed for pressure plate data. Generally
#' seems to work but may be thrown off by severe deformities
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data. List. First item should be a 3D array covering each
#' timepoint of the measurement. z dimension represents time
#' @return String. "LEFT" or "RIGHT"
#' @examples
#' auto_detect_side(pressure_data)
#' @importFrom sf st_polygon st_as_sf st_convex_hull st_combine st_intersection st_area
#' @export

auto_detect_side <- function(pressure_data) {
  # max pressure footprint
  fp <- footprint(pressure_data)
  fp_coords <- sensor_coords(pressure_data)
  P <- c(fp)
  sc_df <- data.frame(x = fp_coords$x_coord, y = fp_coords$y_coord, P = P)
  sc_df <- sc_df[which(P >= 5), ]
  sc_df$P <- NULL
  sc_mat <- as.matrix(sc_df)

  # Bounding box
  mbb <- getMinBBox(sc_mat)
  mbb_coords <- mbb$pts

  # Which sides of BBox are longest?
  if (dist(mbb_coords[c(1, 2), ]) > dist(mbb_coords[c(2, 3), ])) {
    side1 <- mbb_coords[c(1, 2), ]
    side2 <- mbb_coords[c(3, 4), ]
  } else {
    side1 <- mbb_coords[c(2, 3), ]
    side2 <- mbb_coords[c(4, 1), ]
  }

  # order so lowest first
  side1 <- side1[order(side1[, 2]), ]
  side2 <- side2[order(side2[, 2]), ]

  # which side is L/R
  if (side1[1, 1] > side2[1, 1]) {
    side1_ <- side1 # side1_ is right
    side2_ <- side2
  } else {
    side1_ <- side2
    side2_ <- side1
  }

  # get midpoints
  midpoint_top <- (side1_[2, ] + side2_[2, ]) / 2
  midpoint_bottom <- (side1_[1, ] + side2_[1, ]) / 2

  # make lat and med box
  side1_pts <- rbind(side1_, midpoint_top, midpoint_bottom, side1_[1, ])
  side2_pts <- rbind(side2_, midpoint_top, midpoint_bottom, side2_[1, ])
  box1 <- st_polygon(list(side1_pts))
  box2 <- st_polygon(list(side2_pts))

  # make chull
  df.sf <- sc_df %>%
    st_as_sf(coords = c( "x", "y" ))
  fp_chull <- st_convex_hull(st_combine(df.sf))

  # area in each half box
  side1_count <- st_area(st_intersection(fp_chull, box1))
  side2_count <- st_area(st_intersection(fp_chull, box2))

  # side
  if (side1_count < side2_count) {
    side <- "RIGHT"
  } else {
    side <- "LEFT"
  }

  # Return side
  return(side)
}


# =============================================================================

#' @title Force curve
#' @description Generates force vector with option to plot
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data. List. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param plot Logical. If TRUE also plots data as line curve
#' @return Numeric vector containing force values
#' @examples
#' force_curve(pressure_data, plot = TRUE)
#' @importFrom ggplot2 aes ggplot geom_line theme_bw xlab ylab
#' @export

force_curve <- function(pressure_data, plot = FALSE) {
  # check input
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")

  # convert to force
  sens_area <- pressure_data[[3]][1] * pressure_data[[3]][2]
  force_array <- pressure_data[[1]] * sens_area * 1000

  # create empty vector
  force <- rep(NA, times = dim(force_array)[3])

  # find total force for each frame and store in vector
  for (i in 1:dim(force_array)[3]) {
    force[i] <- sum(force_array[, , i])
  }

  # plot, if required
  if (plot == TRUE) {
    # make df
    force_df <- data.frame(force = force,
                           time = seq(0, by = pressure_data[[4]],
                                      length.out = length(force)))

    # plot
    g <- ggplot(force_df, aes(x = time, y = force))
    g <- g + geom_line()
    g <- g + theme_bw()
    g <- g + xlab("time (s)") + ylab("Force (N)")
    print(g)
  }

  # return
  return(force)
}


# =============================================================================

#' @title Pressure curve
#' @description Generate pressure curve with option to plot
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data. List. First item should be a 3D array covering each
#' timepoint of the measurement. z dimension represents time
#' @param variable. String. Whether "peak", "mean", or "pti" should be used
#' @param plot Logical. If TRUE also plots data as line curve
#' @return Numeric vector containing force values
#' @examples
#' pressure_curve(pressure_data, variable = "peak", plot = TRUE)
#' @importFrom ggplot2 aes ggplot geom_line theme_bw xlab ylab
#' @export

pressure_curve <- function(pressure_frames, variable = "peak",
                           plot = FALSE) {
  # check input
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")
  if (is.character(variable) == FALSE)
    stop("variable parameter must be a string")

  # pressure array
  pressure_array <- pressure_data[[1]]

  # create empty vector
  pressure <- rep(NA, times = dim(pressure_array)[3])

  # find pressure for each frame and store in vector
  for (i in 1:dim(pressure_array)[3]) {
    if (variable == "peak") {
      pressure[i] <- max(pressure_array[, , i])
    }
    if (variable == "mean") {
      active_sens <- which(pressure_array[, , i] > 0, arr.ind = TRUE)
      pressure[i] <- mean(pressure_array[active_sens[, 1], active_sens[, 2], i])
    }
  }

  # plot, if required
  if (plot == TRUE) {
    # make df
    pressure_df <- data.frame(pressure = pressure,
                              time = seq(0, by = pressure_data[[3]],
                                         length.out = length(pressure)))

    # plot
    g <- ggplot(pressure_df, aes(x = time, y = pressure))
    g <- g + geom_line()
    g <- g + theme_bw()
    g <- g + xlab("time (s)") + ylab("Pressure (kPa)")
    print(g)
  }

  # return
  return(pressure)
}


# =============================================================================

#' @title Contact area curve
#' @description Generate contact area curve with option to plot
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data. List. First item should be a 3D array covering each
#' timepoint of the measurement. z dimension represents time
#' @param threshold. Numeric. The minimum pressure required for a sensor to be
#' considered active
#' @param plot Logical. If TRUE also plots data as line curve
#' @return Numeric vector containing force values
#' @examples
#' area_curve(pressure_data, plot = TRUE)
#' @importFrom ggplot2 aes ggplot geom_line theme_bw xlab ylab
#' @export

area_curve <- function(pressure_data, threshold = 0, plot = FALSE) {
  # check input
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")
  if (is.numeric(threshold) == FALSE)
    stop("threshold must be a numeric value")

  # pressure array
  pressure_array <- pressure_data[[1]]

  # create empty vector
  area <- rep(NA, times = dim(pressure_array)[3])

  # sensor size
  sen_size <- pressure_data[[2]][1] * pressure_data[[2]][2]

  # find active area for each frame and store in vector
  for (i in 1:dim(pressure_array)[3]) {
    area[i] <- (sum(pressure_array[, , i] > threshold)) * sen_size
  }

  # plot if requested
  if (plot == TRUE) {
    # make df
    area_df <- data.frame(area = area,
                          time = seq(0, by = pressure_data[[3]],
                                     length.out = length(area)))

    # plot
    g <- ggplot(area_df, aes(x = time, y = area))
    g <- g + geom_line()
    g <- g + theme_bw()
    g <- g + xlab("time (s)") + ylab("Area (m^2)")
    print(g)
  }

  # return
  return(area)
}


# =============================================================================

#' @title Center of pressure
#' @description Generates center of pressure coordinates
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement. z dimension represents time
#' @return Data frame with x and y coordinates of COP throughout trial
#' @examples
#' cop(pressure_data)
#' @export

cop <- function(pressure_data) {
  # check input
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")

  # array dimensions
  x <- pressure_data[[1]]
  dims <- dim(x)

  # individual sensor dimensions
  sens_dim <- pressure_data[[3]]

  # loading totals by column
  col_total <- data.frame(matrix(NA, nrow = dims[2], ncol = dims[3]))
  for (i in 1:dims[3]) {col_total[, i] <- colSums(x[, , i])}

  # Loading totals by row
  row_total <- data.frame(matrix(NA, nrow = dims[1], ncol = dims[3]))
  for (i in 1:dims[3]) {row_total[, i] <- rowSums(x[, , i])}

  # Sensor spacing in x direction
  sens_spacing_x <- seq(from = sens_dim[1] / 2, by = sens_dim[1],
                        length.out = dims[2])

  # Sensor spacing in y direction
  sens_spacing_y <- rev(seq(from = sens_dim[2] / 2, by = sens_dim[2],
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

  # combine coordinates into dataframe
  COP_df <- data.frame(x_coord, y_coord)

  # return COP coordinates
  return(COP_df)
}


# =============================================================================

#' @title Footprint
#' @description Find footprint of pressure file
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param variable String. "max" = maximum value of each sensor across full
#' dataset. "mean" = average value of sensors over full dataset."frame" = an
#' individual pressure frame
#' @param frame Integer. Only used if variable = "frame".
#' @param plot Logical. Display pressure image
#' @return Matrix. Maximum or mean values for all sensors
#' @examples
#' footprint(pressure_data, plot = TRUE)

footprint <- function(pressure_data, variable = "max", frame,
                      plot = FALSE) {
  # check input
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")
  stopifnot("Unknown Variable" = variable %in% c("max", "mean", "frame"))
  if (variable == "frame" & missing(frame))
    stop("For a single frame, a value is needed for frame variable")

  # calculate footprint for different variables
  if (variable == "max") {
    mat <- apply(simplify2array(pressure_data[[1]]), 1:2, max)
  }
  if (variable == "mean") {
    mat <- apply(simplify2array(pressure_data[[1]]), 1:2, mean)
  }
  if (variable == "frame") {
    mat <- pressure_data[[1]][,, frame]
  }

  pd <- pressure_data
  pd[[1]] <- mat

  # plot if requested
  if (plot == TRUE) {g <- plot_pressure(pd)}

  # return footprint
  return(mat)
}


# =============================================================================

#' @title Plot pressure
#' @description Produce plot of pressure data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param variable String. "max" = footprint of maximum sensors. "mean" = average
#'   value of sensors over time (usually for static analyses). "frame" = an
#'   individual frame
#' @param frame Integer.
#' @param smooth Logical. If TRUE, plot will interpolate between sensors to
#' increase data density
#' @param plot_COP Logical. If TRUE, overlay COP data on plot. Default = FALSE
#' @param plot_outline Logical. If TRUE, overlay convex hull outline on plot
#' @param plot Logical. If TRUE, plot will be displayed
#' @return ggplot plot object
#' @examples
#' plot_pressure(pressure_data, variable = "max", plot_COP = TRUE)
#' @importFrom ggplot2 ggplot aes geom_raster
#' @export

plot_pressure <- function(pressure_data, variable = "max", smooth = FALSE, frame,
                          plot_COP = FALSE, plot_outline = FALSE,
                          plot = TRUE) {
  # check inputs
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure data must contain array")

  # get footprint
  fp <- footprint(pressure_data, variable = variable)
  fp_rows <- nrow(fp)
  fp_cols <- ncol(fp)

  # interpolate if required
  if (smooth == TRUE) {
    f_print2 <- matrix(rep(NA, (fp_rows * fp_cols * 5)),
                       nrow = fp_rows, fp_cols * 5)
    for (i in 1:nrow(fp)) {
      xa <- approx(fp[i, ], n = 5 * fp_cols)
      f_print2[i, ] <- xa$y
    }

    # Increase number of rows
    f_print3 <- matrix(rep(NA, fp_rows * 5 * fp_cols * 5),
                       nrow = fp_rows * 5, fp_cols * 5)
    for (i in 1:ncol(f_print2)) {
      xa <- approx(f_print2[ ,i], n = 5 * nrow(f_print2))
      f_print3[ ,i] <- xa$y
    }
    fp <- f_print3
  }

  # generate coordinates for each sensor
  sens_coords <- sensor_coords(pressure_data)

  # overall dimensions of array
  dims <- dim(pressure_data[[1]])

  # combine with pressure values
  cor <- cbind(sens_coords, as.vector(fp))
  cor <- as.data.frame(cor)
  colnames(cor) <- c("x", "y", "value")

  # add colors
  plot_cs <- list(cs_breaks = c(0, 15, 40, 60, 100, 150, 220, 300, 5000),
                  cs_cols = c("white", "grey", "blue", "light blue",
                              "green", "yellow", "red", "deeppink"))
  cols <- unlist(plot_cs[[2]])
  color <- rep(NA, times = nrow(cor))
  for (i in 1:nrow(cor)) {
    if (cor$value[i] < 15) {color[i] = 1}
    if (cor$value[i] >= 15 && cor$value[i] < 40) {color[i] = 2}
    if (cor$value[i] >= 40 && cor$value[i] < 60) {color[i] = 3}
    if (cor$value[i] >= 60 && cor$value[i] < 100) {color[i] = 4}
    if (cor$value[i] >= 100 && cor$value[i] < 150) {color[i] = 5}
    if (cor$value[i] >= 150 && cor$value[i] < 220) {color[i] = 6}
    if (cor$value[i] >= 220 && cor$value[i] < 300) {color[i] = 7}
    if (cor$value[i] >= 300) {color[i] = 8}
  }

  # combine with data frame
  cor <- cbind(cor, color)

  ## plot
  g <- ggplot()
  g <- g + geom_raster(data = cor, aes(x = x, y = y, fill = as.factor(color)))
  g <- g + scale_fill_manual(values = cols)
  g <- g + scale_x_continuous(expand = c(0, 0))
  g <- g + scale_y_continuous(expand = c(0, 0))
  g <- g + coord_fixed()

  # add COP?
  if (plot_COP == TRUE) {
    cop_df <- cop(pressure_data)
    g <- g + geom_point(data = cop_df, aes(x = x_coord, y = y_coord))
  }

  # add outline
  if (plot_outline == TRUE) {
    ch_out <- pressure_outline(pressure_data)
    g <- g + geom_path(data = ch_out, aes(x = x_coord, y = y_coord),
                       colour = "black")
    g <- g + geom_point(data = ch_out, aes(x = x_coord, y = y_coord),
                        colour = "purple")
  }

  # formatting
  g <- g + theme_void() + theme(legend.position = "none")

  # display plot immediately if requested
  if (plot == TRUE) {print(g)}

  # return ggplot object
  return(g)
}


# =============================================================================

#' @title Get outline of pressure region
#' Determine outline (convex hull) of pressure measurement
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param frame Integer. Frame number to use
#' @return Polygon. sfg object representing convex hull outline
#' @examples
#' pressure_outline(pressure_data)
#' pressure_outline(pressure_data, frame = 50)
#' @importFrom magrittr "%>%"
#' @importFrom sf st_as_sf st_convex_hull st_combine st_buffer
#' @export

pressure_outline <- function(pressure_data, frame) {
  # check inputs
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure data must contain array")

  # check that this is not pedar data
  if (pressure_data[[2]] == "pedar")
    stop("data cannot be from pedar")

  # determine active sensor coordinates
  if (missing(frame) == TRUE) {
    coords <- sensor_coords(pressure_data)
  } else {
    coords <- sensor_coords(pressure_data, "frame", frame)
  }

  # calculate convex hull
  sens_coords_df <- coords %>%
    st_as_sf(coords = c("x_coord", "y_coord"))
  fp_chull <- st_convex_hull(st_combine(sens_coords_df))

  # add buffer here
  buffer_size <- pressure_data[[3]][1] / 2
  fp_chull_buffer <- st_buffer(fp_chull, buffer_size)

  # return convex hull coordinates
  return(fp_chull_buffer)
}


# =============================================================================

#' @title Animate pressure
#' @description Produce animation (gif) of pressure distribution
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data Array. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param fps Numeric. Number of frames per second in animation
#' @param file Name (inlcuding path) of export file
#' @param preview Logical. Whether to play the animation
#' @return Animation in gif format
#' @examples
#' animate_pressure(pressure_data, fps = 10, "testgif.gif")
#' @importFrom stringr str_ends
#' @importFrom magick image_animate image_write
#' @export

animate_pressure <- function(pressure_data, fps, filename, preview = FALSE) {
  # parameter check
  if (str_ends(filename, ".gif") == FALSE)
    stop("filename must end in .gif")

  # colours
  plot_cs <- list(cs_breaks = c(0, 15, 40, 60, 100, 150, 220, 300, 5000),
                  cs_cols = c("white", "grey", "blue", "light blue",
                                     "green", "yellow", "red", "deeppink"))
                                     cols <- unlist(plot_cs[[2]])

  # generate dataframe with coords for each frame (helper function)
  press_df <- function(pressure_data, frame) {
    # generate coordinates for each sensor
    sens_coords <- sensor_coords(pressure_data)

    # combine with pressure values
    dims <- dim(pressure_data[[1]])
    P <- c(footprint(pressure_data, "frame", frame))
    cor <- cbind(sens_coords, P)
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

    # return
    return(cor)
  }

  # plot
  img <- image_graph(600, 340, res = 96)
  for (i in 1:dim(pressure_data[[1]])[3]) {
    df <- press_df(pressure_data, i)
    g <- ggplot(df, aes(x = x, y = y, fill = as.factor(colour)))
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
    print(g)
  }

  # turn off device and clean up data
  dev.off()
  gc()

  # create animation
  animation <- image_animate(img, fps = fps, optimize = TRUE)

  # save animation
  image_write(animation, filename)

  # preview if required
  if (preview == TRUE) {print(animation)}
}


# =============================================================================

#' @title automask footprint
#' @description Automatically create mask of footprint data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement. z dimension represents time
#' @param side Character. "RIGHT" or "LEFT"
#' @param sens Numeric. Number of frames per second in animation
#' @param plot Logical. Whether to play the animation
#' @return List. Contains polygon with each mask
#' @examples
#' automask(pressure_data, sens = 4, plot = TRUE)
#' @export

automask <- function(pressure_data, side,  sens = 4, plot = FALSE) {
  # Helper functions
  vector_to_polygon <- function(x) {
    xa <- c()
    for (i in 1:((length(x)) / 2)) {
      xa <- paste0(xa, x[(2*i) - 1], " ", x[2 * i], ", ")
    }

    xb <- (paste0("POLYGON ((", xa, x[1], " ", x[2], "))"))
    return(xb)
  } #Turns coords into polygon
  extend_line <- function(x, extend = 1) {
    m1 = (x[4] - x[2]) / (x[3] - x[1])
    c1 = x[2] - (m1 * x[1])
    new_x1 = max(x[c(1,3)]) + extend
    new_x2 = min(x[c(1,3)]) - extend
    new_y1 = m1 * new_x1 + c1
    new_y2 = m1 * new_x2 + c1
    ext_line = c(new_x1, new_y1, new_x2, new_y2)
  } #Adds 1 unit to end of line
  line_int <- function(x, y) {
    #correct for vertical vectors
    if(x[3] == x[1]) {x[1] <- x[1] + 0.0001}
    if(y[3] == y[1]) {y[1] <- y[1] + 0.0001}

    #Find equation constants
    m1 <- (x[4] - x[2]) / (x[3] - x[1])
    m2 <- (y[4] - y[2]) / (y[3] - y[1])
    c1 <- x[2] - (m1 * x[1])
    c2 <- y[2] - (m2 * y[1])

    myCoeffMat <- matrix(c(-m1, 1, -m2, 1), nrow = 2, ncol = 2, byrow = TRUE)
    myRhsMat <- matrix(c(c1, c2), nrow = 2, ncol = 1, byrow = TRUE)
    myInverse <- solve(myCoeffMat)
    myResult <- myInverse %*% myRhsMat

    # return
    return(myResult)
  } #Where two lines intercept


  # ===========================================================================

  # Find footprint (max)
  max_df <- footprint(pressure_data)

  # coordinates
  sens_coords <- sensor_coords(pressure_data)

  # make data frame
  P <- c(max_df)
  em_act_df <- data.frame(x = sens_coords$x_coord, y = sens_coords$y_coord,
                          P = P)
  em_act_df <- em_act_df[which(P >= 5), ]
  em_act_df$P <- NULL


  # ===========================================================================

  # Define minimum bounding box
  em_act_m <- as.matrix(em_act_df)
  mbb <- getMinBBox(em_act_m)
  mbb_df <- data.frame(x = mbb$pts[, 1], y = mbb$pts[, 2])

  # Define convex hull, expanding to include all sensors
  chull_elements <- chull(x = em_act_df$x, y = em_act_df$y)
  chull_polygon <- vector_to_polygon(c(t(em_act_df[chull_elements, ])))
  chull_ex <- gBuffer(readWKT(chull_polygon), width = 0.005,
                      joinStyle = "MITRE")
  chull_ex_df <- fortify(chull_ex)
  chull_ex_df <- data.frame(x = chull_ex_df$long, y = chull_ex_df$lat)


  # ===========================================================================

  # Define angles for dividing lines between metatarsals
  ## Find longest vectors (these are the med and lat edges of the footprint)
  z_dist <- Mod(diff(chull_ex_df$x + 1i * chull_ex_df$y))
  vec_1 <- order(z_dist, decreasing = TRUE)[1]
  vec_2 <- order(z_dist, decreasing = TRUE)[2]

  ## Get coords for longest vector, reorder so lowest (y-axis) is first
  vec_1 <- c(chull_ex_df[vec_1, 1], chull_ex_df[vec_1, 2],
             chull_ex_df[vec_1 + 1, 1], chull_ex_df[vec_1 + 1, 2])
  vec_2 <- c(chull_ex_df[vec_2, 1], chull_ex_df[vec_2,2],
             chull_ex_df[vec_2 + 1, 1], chull_ex_df[vec_2 + 1, 2])
  if (vec_1[2] > vec_1[4]) {vec_1 <- c(vec_1[3], vec_1[4], vec_1[1], vec_1[2])}
  if (vec_2[2] > vec_2[4]) {vec_2 <- c(vec_2[3], vec_2[4], vec_2[1], vec_2[2])}

  ## Define as lat or med
  if (side == "RIGHT" & sum(vec_1[c(1,3)]) > sum(vec_2[c(1,3)])) {
    lat_side = vec_1
    med_side = vec_2
  } else if (side == "RIGHT" & sum(vec_1[c(1,3)]) < sum(vec_2[c(1,3)])) {
    lat_side = vec_2
    med_side = vec_1
  } else if (side == "LEFT" & sum(vec_1[c(1,3)]) > sum(vec_2[c(1,3)])) {
    lat_side = vec_2
    med_side = vec_1
  } else if (side == "LEFT" & sum(vec_1[c(1,3)]) < sum(vec_2[c(1,3)])) {
    lat_side = vec_1
    med_side = vec_2
  }

  ## Find angle between vectors
  lat_v_ang <- c(lat_side[3] - lat_side[1], lat_side[4] - lat_side[2])
  med_v_ang <- c(med_side[3] - med_side[1], med_side[4] - med_side[2])
  lat_v_alpha <- atan2(lat_v_ang[2], lat_v_ang[1]) * (180 / pi)
  med_v_alpha <- atan2(med_v_ang[2], med_v_ang[1]) * (180 / pi)
  alpha <- abs(lat_v_alpha - med_v_alpha)

  # get angle vectors
  if (side == "LEFT") {
    MTH_hx_alpha = med_v_alpha + (0.33 * alpha)
    MTH_12_alpha = med_v_alpha + (0.3  * alpha)
    MTH_22_alpha = med_v_alpha + (0.38 * alpha)
    MTH_23_alpha = med_v_alpha + (0.47 * alpha)
    MTH_34_alpha = med_v_alpha + (0.64 * alpha)
    MTH_45_alpha = med_v_alpha + (0.81 * alpha)
  } else if (side == "RIGHT") {
    MTH_hx_alpha = med_v_alpha - (0.33 * alpha)
    MTH_12_alpha = med_v_alpha - (0.3 * alpha)
    MTH_22_alpha = med_v_alpha - (0.38 * alpha)
    MTH_23_alpha = med_v_alpha - (0.47 * alpha)
    MTH_34_alpha = med_v_alpha - (0.64 * alpha)
    MTH_45_alpha = med_v_alpha - (0.81 * alpha)
  }


  # ===========================================================================

  # Define cutting masks for MTH areas
  ## Find coord pairs for MTH and mid 2nd division lines
  crossing_point <- line_int(vec_1, vec_2)
  MT_lines <- c("MTH_hx_alpha", "MTH_12_alpha", "MTH_22_alpha", "MTH_23_alpha",
                "MTH_34_alpha", "MTH_45_alpha")
  MTH_hx_line <- c(NA, NA, NA, NA)
  MTH_12_line <- c(NA, NA, NA, NA)
  MTH_22_line <- c(NA, NA, NA, NA)
  MTH_23_line <- c(NA, NA, NA, NA)
  MTH_34_line <- c(NA, NA, NA, NA)
  MTH_45_line <- c(NA, NA, NA, NA)

  for (i in 1:6) {
    MT_line_name = paste0(substr(MT_lines[i], 1, 7), "line")
    if (get(MT_lines[i]) <= 90) {
      x2 = cos(get(MT_lines[i]) * pi / 180) * 5 + crossing_point[1]
      y2 = sin(get(MT_lines[i]) * pi / 180) * 5 + crossing_point[2]
    } else if (get(MT_lines[i]) > 90) {
      angle = get(MT_lines[i]) - 90
      x2 = crossing_point[1] - sin(angle * pi / 180) * 5
      y2 = cos(angle * pi / 180) * 5 + crossing_point[2]
    }
    assign(MT_line_name, c(crossing_point[1], crossing_point[2], x2, y2))
  }

  if (side == "RIGHT") {
    MT_hx_lat <- readWKT(vector_to_polygon(c(MTH_hx_line, MTH_hx_line[3] + 1,
                                             MTH_hx_line[4], MTH_hx_line[1] +
                                               1, MTH_hx_line[2])))
    MT_hx_med <- readWKT(vector_to_polygon(c(MTH_hx_line, MTH_hx_line[3] - 1,
                                             MTH_hx_line[4], MTH_hx_line[1] -
                                               1, MTH_hx_line[2])))
    MT_12_lat <- readWKT(vector_to_polygon(c(MTH_12_line, MTH_12_line[3] + 1,
                                             MTH_12_line[4], MTH_12_line[1] +
                                               1, MTH_12_line[2])))
    MT_12_med <- readWKT(vector_to_polygon(c(MTH_12_line, MTH_12_line[3] - 1,
                                             MTH_12_line[4], MTH_12_line[1] -
                                               1, MTH_12_line[2])))
    MT_23_lat <- readWKT(vector_to_polygon(c(MTH_23_line, MTH_23_line[3] + 1,
                                             MTH_23_line[4], MTH_23_line[1] +
                                               1, MTH_23_line[2])))
    MT_23_med <- readWKT(vector_to_polygon(c(MTH_23_line, MTH_23_line[3] - 1,
                                             MTH_23_line[4], MTH_23_line[1] -
                                               1, MTH_23_line[2])))
    MT_34_lat <- readWKT(vector_to_polygon(c(MTH_34_line, MTH_34_line[3] + 1,
                                             MTH_34_line[4], MTH_34_line[1] +
                                               1, MTH_34_line[2])))
    MT_34_med <- readWKT(vector_to_polygon(c(MTH_34_line, MTH_34_line[3] - 1,
                                             MTH_34_line[4], MTH_34_line[1] -
                                               1, MTH_34_line[2])))
    MT_45_lat <- readWKT(vector_to_polygon(c(MTH_45_line, MTH_45_line[3] + 1,
                                             MTH_45_line[4], MTH_45_line[1] +
                                               1, MTH_45_line[2])))
    MT_45_med <- readWKT(vector_to_polygon(c(MTH_45_line, MTH_45_line[3] - 1,
                                             MTH_45_line[4], MTH_45_line[1] -
                                               1, MTH_45_line[2])))
  } else if (side == "LEFT") {
    MT_hx_lat <- readWKT(vector_to_polygon(c(MTH_hx_line, MTH_hx_line[3] - 1,
                                             MTH_hx_line[4], MTH_hx_line[1] -
                                               1, MTH_hx_line[2])))
    MT_hx_med <- readWKT(vector_to_polygon(c(MTH_hx_line, MTH_hx_line[3] + 1,
                                             MTH_hx_line[4], MTH_hx_line[1] +
                                               1, MTH_hx_line[2])))
    MT_12_lat <- readWKT(vector_to_polygon(c(MTH_12_line, MTH_12_line[3] - 1,
                                             MTH_12_line[4], MTH_12_line[1] -
                                               1, MTH_12_line[2])))
    MT_12_med <- readWKT(vector_to_polygon(c(MTH_12_line, MTH_12_line[3] + 1,
                                             MTH_12_line[4], MTH_12_line[1] +
                                               1, MTH_12_line[2])))
    MT_23_lat <- readWKT(vector_to_polygon(c(MTH_23_line, MTH_23_line[3] - 1,
                                             MTH_23_line[4], MTH_23_line[1] -
                                               1, MTH_23_line[2])))
    MT_23_med <- readWKT(vector_to_polygon(c(MTH_23_line, MTH_23_line[3] + 1,
                                             MTH_23_line[4], MTH_23_line[1] +
                                               1, MTH_23_line[2])))
    MT_34_lat <- readWKT(vector_to_polygon(c(MTH_34_line, MTH_34_line[3] - 1,
                                             MTH_34_line[4], MTH_34_line[1] -
                                               1, MTH_34_line[2])))
    MT_34_med <- readWKT(vector_to_polygon(c(MTH_34_line, MTH_34_line[3] + 1,
                                             MTH_34_line[4], MTH_34_line[1] +
                                               1, MTH_34_line[2])))
    MT_45_lat <- readWKT(vector_to_polygon(c(MTH_45_line, MTH_45_line[3] - 1,
                                             MTH_45_line[4], MTH_45_line[1] -
                                               1, MTH_45_line[2])))
    MT_45_med <- readWKT(vector_to_polygon(c(MTH_45_line, MTH_45_line[3] + 1,
                                             MTH_45_line[4], MTH_45_line[1] +
                                               1, MTH_45_line[2])))
  }


  #-------------------------------------------------------------------------

  # Define heel and midfoot cutting areas. Length of foot is taken as min
  # bounding box.
  ## Which sides of BBox are longest?
  s1d <- sqrt((mbb$pts[2,1] - mbb$pts[1,1]) ^ 2 +
                (mbb$pts[2,2] - mbb$pts[1,2]) ^ 2)
  s2d <- sqrt((mbb$pts[3,1] - mbb$pts[2,1]) ^ 2 +
                (mbb$pts[3,2] - mbb$pts[2,2]) ^ 2)
  s3d <- sqrt((mbb$pts[4,1] - mbb$pts[3,1]) ^ 2 +
                (mbb$pts[4,2] - mbb$pts[3,2]) ^ 2)
  s4d <- sqrt((mbb$pts[1,1] - mbb$pts[4,1]) ^ 2 +
                (mbb$pts[1,2] - mbb$pts[4,2]) ^ 2)
  longsides <- order(-c(s1d, s2d, s3d, s4d))
  longsides <- longsides[c(1,2)]

  ## Get coordinate pairs for longside lines
  if (is.element(1, longsides) == TRUE & is.element(3, longsides) == TRUE) {
    longside1 <- as.vector(t(mbb$pts[c(1,2), ]))
    longside2 <- as.vector(t(mbb$pts[c(3,4), ]))
  } else if (is.element(2, longsides) == TRUE &
             is.element(4, longsides) == TRUE) {
    longside1 <- as.vector(t(mbb$pts[c(2,3), ]))
    longside2 <- as.vector(t(mbb$pts[c(4,1), ]))
  }

  ## Reorder coordinates to ensure lowest y is first
  if (longside1[2] > longside1[4]) {longside1 = longside1[c(3, 4, 1, 2)]}
  if (longside2[2] > longside2[4]) {longside2 = longside2[c(3, 4, 1, 2)]}

  ## Find coordinate pairs for 27% and 55% BBox lines
  if(longside1[1] >= longside1[3]) {
    p1x_27 = (abs(longside1[1] - longside1[3]) * 0.27) - max(longside1[c(1,3)])
    p1y_27 = ((longside1[4] - longside1[2]) * 0.27) + longside1[2]
    p1x_55 = (abs(longside1[1] - longside1[3]) * 0.55) - max(longside1[c(1,3)])
    p1y_55 = ((longside1[4] - longside1[2]) * 0.55) + longside1[2]
  } else if(longside1[1] < longside1[3]) {
    p1x_27 = (abs(longside1[1] - longside1[3]) * 0.27) + min(longside1[c(1,3)])
    p1y_27 = ((longside1[4] - longside1[2]) * 0.27) + longside1[2]
    p1x_55 = (abs(longside1[1] - longside1[3]) * 0.55) + min(longside1[c(1,3)])
    p1y_55 = ((longside1[4] - longside1[2]) * 0.55) + longside1[2]
  }

  if (longside2[1] >= longside2[3]) {
    p2x_27 = (abs(longside2[1] - longside2[3]) * 0.27) - max(longside2[c(1,3)])
    p2y_27 = ((longside2[4] - longside2[2]) * 0.27) + longside2[2]
    p2x_55 = (abs(longside2[1] - longside2[3]) * 0.55) - max(longside2[c(1,3)])
    p2y_55 = ((longside2[4] - longside2[2]) * 0.55) + longside2[2]
  } else if (longside2[1] < longside2[3]) {
    p2x_27 = (abs(longside2[1] - longside2[3]) * 0.27) + min(longside2[c(1,3)])
    p2y_27 = ((longside2[4] - longside2[2]) * 0.27) + longside2[2]
    p2x_55 = (abs(longside2[1] - longside2[3]) * 0.55) + min(longside2[c(1,3)])
    p2y_55 = ((longside2[4] - longside2[2]) * 0.55) + longside2[2]
  }

  ## Make 27% and 55% line vectors
  l_27 <- c(p1x_27, p1y_27, p2x_27, p2y_27)
  l_55 <- c(p1x_55, p1y_55, p2x_55, p2y_55)

  ## Mid 2nd line angle
  alpha_2 <- (MTH_12_alpha + MTH_23_alpha) / 2

  ## Find where 2nd mid line and 27% and 55% lines cross
  l27_int <- as.vector(line_int(MTH_22_line, l_27))
  l55_int <- as.vector(line_int(MTH_22_line, l_55))

  ## Eqn for perpindicular lines
  m3 <- -1 / ((MTH_22_line[4] - MTH_22_line[2]) /
                (MTH_22_line[3] - MTH_22_line[1]))
  c3 <- l27_int[2] - (m3 * l27_int[1])
  c4 <- l55_int[2] - (m3 * l55_int[1])
  l27 <- c(1, ((m3 * 1) + c3), -1, ((m3 * -1) + c3))
  l55 <- c(1, ((m3 * 1) + c4), -1, ((m3 * -1) + c4))

  ## Form polygons that will make cuts
  heel_cut_dist <- readWKT(vector_to_polygon(c(l27[1], l27[2], l27[3], l27[4],
                                               l27[3], l27[4] + 1, l27[1],
                                               l27[2] + 1)))
  mfoot_cut_prox <- readWKT(vector_to_polygon(c(l27[1], l27[2], l27[3], l27[4],
                                                l27[3], l27[4] - 1, l27[1],
                                                l27[2] - 1)))
  mfoot_cut_dist <- readWKT(vector_to_polygon(c(l55[1], l55[2], l55[3], l55[4],
                                                l55[3], l55[4] + 1, l55[1],
                                                l55[2] + 1)))
  ffoot_cut_prox <- readWKT(vector_to_polygon(c(l55[1], l55[2], l55[3], l55[4],
                                                l55[3], l55[4] - 1,
                                                l55[1], l55[2] - 1)))


  # ===========================================================================

  # Define toe area
  ## Make smaller dataframes
  max_df2 <- max_df[1:round(nrow(max_df) * 0.5), ]
  max_df3 <- max_df2[2:nrow(max_df2), ]

  ## Find front edge of toes
  f_edge <- apply(max_df2, 2, function(x) which(x > 0)[1])

  ## Find side-specific first active column (columns with hallux essentially)
  if(side == "RIGHT") {f_edge_act <- which(!is.na(f_edge))[1]}
  if(side == "LEFT") {
    f_edge_act <- length(f_edge) - which(!is.na(rev(f_edge)))[1]
  }

  ## Prepare empty variables that will be filled using toe line algorithm
  lmin <- c(rep(NA), ncol(max_df2))
  lmin2 <- list()
  lmin3 <- list()
  d_mat <- matrix(NA, nrow(max_df2) - 1, ncol(max_df2))

  ## TRUE/FALSE matrices for local minima and maxima
  TFmat <- rollapply(max_df2, 2, function(x) which.min(x) == 2)
  FTmat <- rollapply(max_df2, 2, function(x) which.max(x) == 2)

  ## Difference matrix
  d_mat <- apply(max_df2, 2, function(x) diff(x))

  ## In TFmat (minima), make row 1 FALSE
  TFmat[1, ] <- FALSE

  ## In TFmat, adjust to ensure only local minima are detected
  TFmat[TFmat == TRUE & d_mat == 0 & !max_df3 == 0] <- FALSE
  TFmat[nrow(TFmat), ] <- FALSE

  ## Produce initial toe line
  for (i in 1:ncol(max_df2)) {
    # Get locations of maxs and mins from TF and FT matrices
    lmin2[[i]] = which(rollapply(TFmat[ ,i], 2, identical, c(TRUE, FALSE)))
    lmin3[[i]] = which(rollapply(FTmat[ ,i], 2, identical, c(TRUE, FALSE)))

    if(length(lmin2[[i]]) >= 2 & lmin2[[i]][2] - lmin2[[i]][1] < sens) {
      lmin2[[i]] = lmin2[[i]][-1]
    }

    # At point where lmin2 is still a list, check that each identified min is
    # at least 10% lower than its two neighbours in the column
    for (j in 1:length(lmin2[[i]])) {
      if(length(lmin2[[i]]) >= 1) {
        if (max_df2[lmin2[[i]][j] + 1,i] * 1.1 >
            max_df2[lmin2[[i]][j] + 2,i] &
            max_df2[lmin2[[i]][j] + 1,i] * 1.1 >
            max_df2[lmin2[[i]][j],i ]) {
          lmin2[[i]][j] = NA
        }
      }
    }

    if(length(lmin3[[i]]) == 1 & lmin3[[i]][1] > 9) {
      lmin2[[i]][1] = NA
    }

    # Remove 1st element if too close to front edge
    if (!is.na(lmin2[[i]][1]) & lmin2[[i]][1] < f_edge[i] + 3) {
      lmin2[[i]][1] = NA
    }


    lmin[i] <- lmin2[[i]][1]
  }

  ## Fix lateral columns (if required)
  h_lmin <- round(length(lmin) / 2)
  if (side == "RIGHT") {
    for (i in (length(lmin) - 4):length(lmin)) {
      lmin2[[i]] = which(rollapply(TFmat[ ,i], 2, identical, c(TRUE, FALSE)))
      lmin3[[i]] = which(rollapply(FTmat[ ,i], 2, identical, c(TRUE, FALSE)))
      if (length(lmin3[[i]]) == 0) {lmin3[[i]] = NA}
      if (length(lmin2[[i]]) == 0) {lmin2[[i]] = NA}
      if (length(lmin3[[i]]) > 1) {lmin3[[i]] = lmin3[[i]][length(lmin3[[i]])]}
      if (lmin2[[i]][1] < lmin3[[i]]) {
        lmin2[[i]] = max(lmin2[[i]][lmin2[[i]] <= lmin3[[i]]])
      }
      lmin[i] = lmin2[[i]]
    }
  }

  if (side == "LEFT") {
    for (i in 1:4) {
      lmin2[[i]] = which(rollapply(TFmat[ ,i], 2, identical, c(TRUE, FALSE)))
      lmin3[[i]] = which(rollapply(FTmat[ ,i], 2, identical, c(TRUE, FALSE)))
      if (length(lmin3[[i]]) == 0) {lmin3[[i]] = NA}
      if (length(lmin2[[i]]) == 0) {lmin2[[i]] = NA}
      if (length(lmin3[[i]]) > 1) {lmin3[[i]] = lmin3[[i]][length(lmin3[[i]])]}
      if (lmin2[[i]][1] < lmin3[[i]]) {
        lmin2[[i]] = max(lmin2[[i]][lmin2[[i]] <= lmin3[[i]]])
      }
      lmin[i] = lmin2[[i]][1]
    }
  }

  ## Remove from lmin points too far from front edge (all)
  lmin[lmin > 18] <- NA
  lmin[lmin > (f_edge + 10)] <- NA

  ## Remove points too close to front edge (medial side) from lmin
  if (side == "RIGHT") {
    for (i in f_edge_act:(f_edge_act + 8)) {
      if (!is.na(lmin[i]) & lmin[i] < f_edge[i] + 3) {lmin[i] = NA}
    }
  }
  if (side == "LEFT") {
    for (i in (f_edge_act - 8):f_edge_act) {
      if (!is.na(lmin[i]) & lmin[i] < f_edge[i] + 3) {lmin[i] = NA}
    }
  }

  ## Remove "kinks" in line
  for (i in 2:(length(lmin) - 1)) {
    if (!any(is.na(c(lmin[i - 1], lmin[i], lmin[i + 1]))) &
        lmin[i] >= lmin[i + 1] + 2 & lmin[i] >= lmin[i - 1] + 1) {
      lmin[i] <- NA
    }
    if (!any(is.na(c(lmin[i - 1], lmin[i], lmin[i + 1]))) &
        lmin[i] >= lmin[i + 1] + 1 & lmin[i] >= lmin[i - 1] + 2) {
      lmin[i] <- NA
    }
    if (!any(is.na(c(lmin[i - 1], lmin[i], lmin[i + 1]))) &
        lmin[i] <= lmin[i + 1] - 2 & lmin[i] <= lmin[i - 1] - 1) {
      lmin[i] <- NA
    }
    if (!any(is.na(c(lmin[i - 1], lmin[i], lmin[i + 1]))) &
        lmin[i] <= lmin[i + 1] - 1 & lmin[i] <= lmin[i - 1] - 2) {
      lmin[i] <- NA
    }
  }

  ## Fix outer edges
  if (side == "LEFT") {
    lmin3a = which(rollapply(FTmat[ ,1], 2, identical, c(TRUE, FALSE)))
    lmin3b = which(rollapply(TFmat[ ,1], 2, identical, c(FALSE, TRUE)))
    if (is.na(lmin[1]) == TRUE & length(lmin3a) >= 3) {
      lmin[1] = lmin3b[2]
    }
  } else if (side == "RIGHT") {
    lmin3a = which(rollapply(FTmat[ ,length(lmin)], 2, identical,
                             c(TRUE, FALSE)))
    lmin3b = which(rollapply(TFmat[ ,length(lmin)], 2, identical,
                             c(FALSE, TRUE)))
    if (is.na(lmin[length(lmin)]) == TRUE & length(lmin3a) >= 3) {
      lmin[length(lmin)] = lmin3b[2]
    }
  }

  # Identify columns to keep i.e. not NA
  cols_keep <- which(!is.na(lmin))

  # Identify columns to skip
  cols_skip <- which(is.na(lmin))

  # Replace NAs in lmin with 1s (just for logic purposes, remove later)
  lmin[is.na(lmin)] <- 1

  # Does lmin have == value in the next rows? Adjust if TRUE
  lmin4 <- lmin
  for (i in cols_keep) {
    if (nrow(max_df2) >= lmin[i] + 2) {
      if (max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 2, i]) {
        lmin4[i] = lmin4[i] + 0.5
      }
    }

    if (nrow(max_df2) >= lmin[i] + 3) {
      if (max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 2, i] &
          max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 3, i]) {
        lmin4[i] = lmin4[i] + 0.5
      }
    }

    if (nrow(max_df2) >= lmin[i] + 4) {
      if (max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 2, i] &
          max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 3, i] &
          max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 4, i]) {
        lmin4[i] = lmin4[i] + 0.5
      }
    }

    if (nrow(max_df2) >= lmin[i] + 5) {
      if (max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 2, i] &
          max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 3, i] &
          max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 4, i] &
          max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 5, i]) {
        lmin4[i] = lmin4[i] + 0.5
      }
    }

    if (nrow(max_df2) >= lmin[i] + 6) {
      if (max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 2, i] &
          max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 3, i] &
          max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 4, i] &
          max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 5, i] &
          max_df2[lmin[i] + 1, i] == max_df2[lmin[i] + 6, i]) {
        lmin4[i] = lmin4[i] + 0.5
      }
    }
  }

  # Turn lmin into coordinates for toe line
  x_cor_lmin <- seq(from = 0.0025, by = 0.005, length.out = ncol(max_df))
  y_cor_lmin <- ((nrow(max_df) * 0.005) - 0.0025) - (lmin4 * 0.005)
  toe_line_df <- data.frame(x = x_cor_lmin, y = y_cor_lmin)
  toe_line_df <- toe_line_df[-cols_skip, ]
  toe_line <- c(t(toe_line_df))

  # Add points to make cut box distal
  toe_cut_dist <- c(toe_line, toe_line[length(toe_line) - 1] + 1,
                    toe_line[length(toe_line)],
                    toe_line[length(toe_line) - 1] + 1,
                    toe_line[length(toe_line)] + 1,
                    toe_line[1] - 1, toe_line[2] + 1,
                    toe_line[1] - 1, toe_line[2])
  toe_cut_dist <- readWKT(vector_to_polygon(toe_cut_dist))

  # Add points to make cut box proximal
  toe_cut_prox <- c(toe_line, toe_line[length(toe_line) - 1] + 1,
                    toe_line[length(toe_line)],
                    toe_line[length(toe_line) - 1] + 1,
                    toe_line[length(toe_line)] - 1, toe_line[1] - 1,
                    toe_line[2]-1, toe_line[1] - 1, toe_line[2])
  toe_cut_prox <- readWKT(vector_to_polygon(toe_cut_prox))


  # ===========================================================================

  # Make all masks
  heel_mask <- gDifference(chull_ex, heel_cut_dist)
  midfoot_mask <- gDifference(chull_ex, mfoot_cut_prox)
  midfoot_mask <- gDifference(midfoot_mask, mfoot_cut_dist)
  forefoot_mask <- gDifference(chull_ex, ffoot_cut_prox)
  forefoot_mask <- gDifference(forefoot_mask, toe_cut_dist)
  hallux_mask <- gDifference(chull_ex, toe_cut_prox)
  hallux_mask <- gDifference(hallux_mask, MT_hx_lat)
  l_toes_mask <- gDifference(chull_ex, toe_cut_prox)
  l_toes_mask <- gDifference(l_toes_mask, MT_hx_med)
  MTH_1_mask <- gDifference(forefoot_mask, MT_12_lat)
  MTH_2_mask <- gDifference(forefoot_mask, MT_12_med)
  MTH_2_mask <- gDifference(MTH_2_mask, MT_23_lat)
  MTH_3_mask <- gDifference(forefoot_mask, MT_23_med)
  MTH_3_mask <- gDifference(MTH_3_mask, MT_34_lat)
  MTH_4_mask <- gDifference(forefoot_mask, MT_34_med)
  MTH_4_mask <- gDifference(MTH_4_mask, MT_45_lat)
  MTH_5_mask <- gDifference(forefoot_mask, MT_45_med)

  heel_mask_ <- st_as_sf(heel_mask)
  midfoot_mask_ <- st_as_sf(midfoot_mask)
  forefoot_mask_ <- st_as_sf(forefoot_mask)
  hallux_mask_ <- st_as_sf(hallux_mask)
  l_toes_mask_ <- st_as_sf(l_toes_mask)
  MTH_1_mask_ <- st_as_sf(MTH_1_mask)
  MTH_2_mask_ <- st_as_sf(MTH_2_mask)
  MTH_3_mask_ <- st_as_sf(MTH_3_mask)
  MTH_4_mask_ <- st_as_sf(MTH_4_mask)
  MTH_5_mask_ <- st_as_sf(MTH_5_mask)

  masks_emed <- list("heel_mask", "midfoot_mask", "forefoot_mask",
                     "hallux_mask", "l_toes_mask", "MTH_1_mask", "MTH_2_mask",
                     "MTH_3_mask", "MTH_4_mask", "MTH_5_mask")
  masks_emed_df <- list("heel_mask_df", "midfoot_mask_df", "forefoot_mask_df",
                        "hallux_mask_df", "l_toes_mask_df", "MTH_1_mask_df",
                        "MTH_2_mask_df", "MTH_3_mask_df", "MTH_4_mask_df",
                        "MTH_5_mask_df")
  mask_list <- list(heel_mask = heel_mask_,
                    midfoot_mask = midfoot_mask_,
                    forefoot_mask = forefoot_mask_,
                    hallux_mask = hallux_mask_,
                    l_toe_mask = l_toes_mask_,
                    MTH_1_mask = MTH_1_mask_,
                    MTH_2_mask = MTH_2_mask_,
                    MTH_3_mask = MTH_3_mask_,
                    MTH_4_mask = MTH_4_mask_,
                    MTH_5_mask = MTH_5_mask_)


  # ===========================================================================

  # Plot Footprint and masks if required
  if (plot == TRUE) {
    # make masks into df format
    mask_df <- masks_2_df(mask_list)

    # Plot footprint and masks
    g <- plot_pressure(pressure_data, "max", plot = FALSE)
    g <- g + geom_path(data = mask_df, aes(x = x, y = y, group = mask),
                       color = "red", size = 1)
    print(g)
  }

  # Return masks for analysis
  return(mask_list)
}



# =============================================================================

#' @title Create mask
#' @description Manually create mask region
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @param image String. "max" = footprint of maximum sensors. "mean"
#'   average value of sensors over time (usually for static analyses)
#' @return List New mask is added to the relevant A 3D array covering each timepoint of the measurement for the
#'   selected region. z dimension represents time
#' @importFrom
#' @export

create_mask <- function(pressure_data, image = "max") {
  # plot footprint
  g <- plot_footprint(pressure_data[[1]])
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

#' @title Edit mask
#' @description Manually edit mask
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @return List.
#' @examples
#' edit_mask(pressure_data)

edit_mask <- function(pressure_data) {

}


# =============================================================================

#' @title CPEI
#' @description Determine Center of Pressure Excursion Index (CPEI) for
#' footprint pressure data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_frames List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @param side String. "right" or "left". Required for automatic detection of
#'   points
#' @param plot_result Logical. Plots pressure image with COP and CPEI overlaid
#' @return Numeric. CPEI value
#' @examples
#' cpei(pressure_data, side = "right", plot = TRUE)

cpei <- function(pressure_data, side, plot_result = FALSE) {
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

  if (auto_worked == "c") {
    # plot footprint
    g <- plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)

    # select points
    start_point <- manually_select(1, "select the most medial point near the start of the COP")
    end_point <- manually_select(1, "select the most medial point near the end of the COP")
  }

  # Calculate CPEI
  ## determine foot length
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

#' Analyze masked regions of pressure data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param masks. List. Masks used to define the regions to be analysed
#' @param partial_sensors Logical Defines how sensors that do not
#'   lie wholly within mask are dealt with. If FALSE, they will be excluded;
#'   if TRUE, for relevant variables their contribution will be weighted by the
#'   proportion of the sensor that falls within the mask border
#' @param variable String. Variable to be determined. "peak_sensor",
#'  "peak_sensor_ts", "peak_mask", "peak_mask_ts", "mean_mask",
#'  "contact_area_peak", "contact_area_ts", "pti_1", "pti_2", force_peak",
#'  "force_ts"
#' @return Data frame.
#' @examples
#'  mask_analysis(pressure_data, masks, TRUE, variable = "force_ts")

mask_analysis <- function(pressure_data, masks, partial_sensors = FALSE,
                          variable = "peak_sensor") {
  # Helper functions
  sensor_to_polygon <- function(act_sens, row_no, max_df_rows, max_df_cols) {
    sensor_x_coord <- (act_sens[row_no, 2] * 0.005) - 0.0025
    sensor_y_coord <- ((max_df_rows * 0.005) - 0.0025) -
      act_sens[row_no, 1] * 0.005
    x1 <- sensor_x_coord - 0.0025
    y1 <- sensor_y_coord + 0.0025
    x2 <- sensor_x_coord + 0.0025
    y2 <- sensor_y_coord + 0.0025
    x3 <- sensor_x_coord + 0.0025
    y3 <- sensor_y_coord - 0.0025
    x4 <- sensor_x_coord - 0.0025
    y4 <- sensor_y_coord - 0.0025
    sens_polygon <- st_polygon(list(matrix(c(x1, x2, x3, x4, x1,
                                             y1, y2, y3, y4, y1), 5, 2)))
    return(sens_polygon)
  }

  # sensor area
  sensor_area <- pressure_data[[3]][1] * pressure_data[[3]][2]

  # Make active sensors into polygons
  ## Find max footprint
  max_df <- apply(simplify2array(pressure_data[[1]]), 1:2, max)

  ## overall dimensions
  dims <- dim(max_df)

  ## Find which sensors are non zero
  act_sens <- which(max_df > 0, arr.ind = TRUE)

  ## make list of polygons
  act_sens_poly <- list()
  for (i in 1:nrow(act_sens)) {
    act_sens_poly[[i]] <- sensor_to_polygon(act_sens, i, dims[1], dims[2])
  }

  # For each region mask, find which polygons intersect
  sens_mask_df <- matrix(rep(0, length.out = (nrow(act_sens) * length(masks))),
                         nrow = nrow(act_sens), ncol = length(masks))
  for (i in 1:length(masks)){
    for (j in 1:length(act_sens_poly)) {
      x <- st_intersects(masks[[i]], act_sens_poly[[j]])
      if (identical(x[[1]], integer(0)) == FALSE) {
        y <- st_intersection(masks[[i]], act_sens_poly[[j]])
        sens_mask_df[j, i] <- st_area(y) / sensor_area
      }
    }
  }

  # combine sensor locations with mask areas
  sens_mask_df <- cbind(act_sens, sens_mask_df)

  # Analyse regions for maximum value of any sensor within region during trial
  if (variable == "press_peak_sensor") {
    peak_sens <- rep(NA, times = length(masks))
    for (i in 1:length(overlap_list)) {
      peak_sens[i] <- max(max_df[act_sens[overlap_list[[i]],]])
    }
  }

  # Analyse regions for maximum value of any sensor for each measurement frame
  if (variable == "press_peak_sensor_ts") {
    mask_pp <- matrix(rep(0, length.out = dim(pressure_data[[1]])[3] * length(masks)),
                      nrow = dim(pressure_data[[1]])[3],
                      ncol = length(masks))
    for (mask in seq_along(masks)) {
      mask_mat <- sens_mask_df[, (mask + 2)]
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- max(pressure_data[[1]][, , i])
        pp <- rep(0, length.out = length(mask_mat))
        for (j in 1:length(mask_mat)) (
          pp[j] <- max(P[sens_mask_df[j, 1], sens_mask_df[j, 2]])
        )
        mask_pp[i, mask] <- pp
      }
    }
  }

  # Analyse regions for peak regional pressure (defined as the maximum mean
  # pressure of active sensors in region during the trial)
  if (variable == "press_peak_region") {

  }

  # Analyse regions for peak regional pressure (defined as the maximum mean
  # pressure of active sensors in region during the trial). Outputs 101 point vector
  if (variable == "press_peak_region_ts") {

  }

  # Analyse regions for overall mean pressure during the trial
  if (variable == "press_mean_region") {

  }

  # Analyse regions for maximum contact area of region during trial
  if (variable == "contact_area_peak") {

  }

  # Analyse regions for contact area throughout the trial (outputs vector)
  if (variable == "contact_area_ts") {

  }

  # Analyse regions for pressure time integral (Novel definition)
  if (variable == "press_ti_1") {

  }

  # Analyse regions for pressure time integral (Melai definition)
  if (variable == "press_ti_2") {

  }

  # Analyse regions for maximum force during the trial
  if (variable == "force_peak") {

  }

  # Analyse regions for force throughout the trial (outputs vector)
  if (variable == "force_ts") {
    mask_force <- matrix(rep(0, length.out = dim(pressure_data[[1]])[3] * length(masks)),
                         nrow = dim(pressure_data[[1]])[3],
                         ncol = length(masks))
    for (mask in seq_along(masks)) {
      mask_mat <- sens_mask_df[, (mask + 2)]
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- pressure_data[[1]][, , i] * sensor_area * 1000
        force <- rep(0, length.out = length(mask_mat))
        for (j in 1:length(mask_mat)) (
          force[j] <- P[sens_mask_df[j, 1], sens_mask_df[j, 2]] * mask_mat[j]
        )
        mask_force[i, mask] <- sum(force)
      }
    }
  }

  # return
  mask_force <- as.data.frame(mask_force)
  colnames(mask_force) <- names(masks)
  return(mask_force)
}


# =============================================================================
# =============================================================================

# helper functions

#' @title Get coordinates of active sensors
#' @description Produces a data frame with coordinates of sensors
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param pressure_image. Character. Which pressure image to use. Options are
#' "all_active", "all", or "frame".
#' @param frame Numeric. If pressure image is frame, the numeric value should be
#' provided here
#' @return Data frame. x and y coordinates of sensors
#' @examples
#' sensor_coords(pressure_data)
#' @noRd

sensor_coords <- function(pressure_data, pressure_image = "all_active", frame) {
  # pressure image
  if (pressure_image == "all_active") {
    sens <- footprint(pressure_data, variable = "max")
  }
  if (pressure_image == "all") {
    dims <- dim(pressure_data[[1]])
    sens <- matrix(rep(10, length.out = dims[1] * dims[2]), nrow = dims[1],
                   ncol = dims[2])
  }
  if (pressure_image == "frame") {
    sens <- footprint(pressure_data, variable = "frame", frame = frame)
  }

  # dimensions
  sens_x <- pressure_data[[3]][1]
  sens_y <- pressure_data[[3]][2]

  # data frame with active sensors as coordinates
  x_cor <- seq(from = sens_x / 2, by = sens_x, length.out = ncol(sens))
  x_cor <- rep(x_cor, each = nrow(sens))
  y_cor <- seq(from = (sens_y / 2) + ((nrow(sens) - 1) * sens_y),
               by = (-1 * sens_y), length.out = nrow(sens))
  y_cor <- rep(y_cor, times = ncol(sens))
  coords <- data.frame(x_coord = x_cor, y_coord = y_cor)

  # remove inactive sensors
  P <- c(sens)
  coords <- coords[which(P > 0), ]

  # return sensor coordinates
  return(coords)
}


#' interpolation function
#' @importFrom stats approx
approxP <- function(x, interp_to) {
  y <- approx(x, n = interp_to)
  y$x <- NULL
  y <- unname(unlist(y))
  return(y)
}

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


#' Get minimum bounding box (angle)
#' @param xy Matrix. n row by 2 column with x-y coords of points to be bounded
getMinBBox <- function(xy) {
  stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) >= 2, ncol(xy) == 2)

  ## rotating calipers algorithm using the convex hull
  H    <- chull(xy)                    # hull indices, vertices ordered clockwise
  n    <- length(H)                    # number of hull vertices
  hull <- xy[H, ]                      # hull vertices

  ## unit basis vectors for all subspaces spanned by the hull edges
  hDir  <- diff(rbind(hull, hull[1,])) # account for circular hull vertices
  hLens <- sqrt(rowSums(hDir^2))       # length of basis vectors
  huDir <- diag(1/hLens) %*% hDir      # scaled to unit length

  ## unit basis vectors for the orthogonal subspaces
  ## rotation by 90 deg -> y' = x, x' = -y
  ouDir <- cbind(-huDir[ , 2], huDir[ , 1])

  ## project hull vertices on the subspaces spanned by the hull edges, and on
  ## the subspaces spanned by their orthogonal complements - in subspace coords
  projMat <- rbind(huDir, ouDir) %*% t(hull)

  ## range of projections and corresponding width/height of bounding rectangle
  rangeH  <- matrix(numeric(n*2), ncol=2)   # hull edge
  rangeO  <- matrix(numeric(n*2), ncol=2)   # orth subspace
  widths  <- numeric(n)
  heights <- numeric(n)
  for(i in seq(along=H)) {
    rangeH[i, ] <- range(projMat[  i, ])
    rangeO[i, ] <- range(projMat[n+i, ])  # orth subspace is in 2nd half
    widths[i]   <- abs(diff(rangeH[i, ]))
    heights[i]  <- abs(diff(rangeO[i, ]))
  }

  ## extreme projections for min-area rect in subspace coordinates
  eMin  <- which.min(widths*heights)   # hull edge leading to minimum-area
  hProj <- rbind(   rangeH[eMin, ], 0)
  oProj <- rbind(0, rangeO[eMin, ])

  ## move projections to rectangle corners
  hPts <- sweep(hProj, 1, oProj[ , 1], "+")
  oPts <- sweep(hProj, 1, oProj[ , 2], "+")

  ## corners in standard coordinates, rows = x,y, columns = corners
  ## in combined (4x2)-matrix: reverse point order to be usable in polygon()
  basis <- cbind(huDir[eMin, ], ouDir[eMin, ])  # basis formed by hull edge and orth
  hCorn <- basis %*% hPts
  oCorn <- basis %*% oPts
  pts   <- t(cbind(hCorn, oCorn[ , c(2, 1)]))

  return(list(pts=pts, width=widths[eMin], height=heights[eMin]))
}


#' sf masks to dataframe
masks_2_df <- function(masks) {
  # no. of masks
  n_masks <- length(masks)

  # empty df
  df <- data.frame(mask = factor(),
                   x = double(),
                   y = double())

  # get coordinates of each mask and store in df
  for (mask in seq_along(masks)) {
    coords <- st_coordinates(masks[[mask]])[, c(1:2)]
    mask_name <- rep(names(masks)[mask], length.out = nrow(coords))
    df_ <- data.frame(mask = mask_name, x = coords[, 1], y = coords[, 2])
    df <- bind_rows(df, df_)
  }

  # return
  return(df)
}
