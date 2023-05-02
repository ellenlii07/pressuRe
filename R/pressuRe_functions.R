# to do (current version)
# CPEI manual edit to be built into function
# Do we have all pedar insoles? Double check areas
# add support for pliance
# edit mask needs to be written
# add more input tests to throw errors
# cop for pedar
# check/make pedar work for a lot of these functions
# automask2: toe line modifications
# fscan processing needs to be checked (work with NA?)
# in edit_mask, make edit_list a vector that works with numbers or names?

# to do (future)
# global pressure_import function (leave for V2)
# create masks for iscan during startup

# data list:
## Array. pressure data
## String. data type (usually collection system, e.g. emed)
## Numeric. sens_size. sensor size
## Numeric. Single number time between samples
## List. Mask list
## events (for example, to define start/end of individual steps for insole data)


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
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' @importFrom stringr str_extract_all str_detect
#' @importFrom utils read.fwf read.table
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
  breaks <- breaks[str_detect(frame_type, "Pict")]

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
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. String with sensor type
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#' @examples
#' pedar_data <- system.file("extdata", "pedar_example.asc", package = "pressuRe")
#' pressure_data <- load_pedar(pedar_data)
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

#' @title Load Tekscan data
#' @description Imports and formats files collected on tekscan systems and
#'    exported from Tekscan software
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#'  @examples
#' tekscan_data <- system.file("extdata", "fscan_testL.asf", package = "pressuRe")
#' tekscan_data <- system.file("extdata", "iscan_test.csv", package = "pressuRe")
#' pressure_data <- load_tekscan(tekscan_data)
#'  @importFrom
#'  @export

load_tekscan <- function(pressure_filepath) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  file_ext <- str_split(basename(pressure_filepath), "\\.")[[1]][2]
  if (file_ext != "asf" & file_ext != "csv")
    stop("incorrect file extension, expected .asf or .csv")

  # Read unformatted data
  pressure_raw <- readLines(pressure_filepath, warn = FALSE)

  # get sensor size
  sens_size <- c(NA, NA)
  sens_width <- which(grepl("ROW_SPACING", pressure_raw))
  sens_width_ <- str_extract_all(pressure_raw[sens_width], "\\d+\\.\\d+")
  sens_size[1] <- as.numeric(unlist(sens_width_))
  sens_height <- which(grepl("COL_SPACING", pressure_raw))
  sens_height_ <- str_extract_all(pressure_raw[sens_height], "\\d+\\.\\d+")
  sens_size[2] <- as.numeric(unlist(sens_height_))
  if (str_detect(pressure_raw[sens_width], "mm") == TRUE) {
    sens_size <- sens_size * 0.001
  }
  if (str_detect(pressure_raw[sens_width], "centimeters") == TRUE) {
    sens_size <- sens_size * 0.01
  }

  # get capture frequency
  time_line <- which(grepl("SECONDS_PER_FRAME", pressure_raw))
  time_ln <- str_split(pressure_raw[time_line], "FRAME ")[[1]][2]
  time <- as.numeric(unlist(str_extract_all(time_ln, "\\d+\\.\\d+")))

  # determine position breaks
  breaks <- grep("Frame", pressure_raw)

  # determine matrix dimensions
  col_n <- length(strsplit(pressure_raw[breaks[1] + 1], ",")[[1]])
  row_n <- breaks[2] - breaks[1] - 2

  # empty array
  pressure_array <- array(0, dim = c(row_n, col_n, length(breaks)))

  # fill array
  for (i in 1:length(breaks)) {
    # frame
    y <- pressure_raw[(breaks[i] + 1):(breaks[i] + row_n)]
    z <- read.table(textConnection(y), sep = ",")
    z[z == "B"] <- -1
    z <- as.data.frame(lapply(z, as.numeric))
    pressure_array[, , i] <- as.matrix(z)
  }

  # return formatted data
  return(list(pressure_array = pressure_array, pressure_system = "tekscan",
              sens_size = sens_size,
              time = time, masks = NULL, events = NULL))
}


# =============================================================================

#' @title Load footscan data
#' @description Imports and formats files collected on footscan systems
#' (formerly RSScan)
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_filepath String. Filepath pointing to emed pressure file
#' @return A list with information about the pressure data.
#' \itemize{
#'   \item pressure_array. 3D array covering each timepoint of the measurement.
#'            z dimension represents time
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#'  @examples
#' footscan_data <- system.file("extdata", "footscan_test.xls", package = "pressuRe")
#' pressure_data <- load_footscan(footscan_data)
#'  @importFrom readxl read_excel
#'  @export

load_footscan <- function(pressure_filepath) {
  # check parameters
  ## file exists
  if (file.exists(pressure_filepath) == FALSE)
    stop("file does not exist")

  ## extension is correct
  if (str_split(basename(pressure_filepath), "\\.")[[1]][2] != "xls")
    stop("incorrect file extension, expected .xls")

  # Read unformatted data
  pressure_raw <- readxl::read_excel(pressure_filepath,
                                     .name_repair = "minimal")

  # sensor size
  sens_size <- c(0.00508, 0.00762)
  sens_area <- sens_size[1] * sens_size[2]

  # determine position breaks
  c1 <- unname(unlist(as.vector(pressure_raw[, 1])))
  breaks <- grep("Frame", c1)

  # get capture frequency
  time <- regmatches(c1[breaks[2]], gregexpr("(?<=\\().*?(?=\\))",
                                             c1[breaks[2]], perl = T))[[1]]
  time <- as.numeric(unlist(str_split(time, " "))[1]) * 0.001

  # determine matrix dimensions
  col_n <- ncol(pressure_raw)
  row_n <- breaks[2] - breaks[1] - 2

  # empty array
  pressure_array <- array(0, dim = c(row_n, col_n, length(breaks)))

  # frame
  for (i in 1:length(breaks)) {
    z <- pressure_raw[(breaks[i] + 1):(breaks[i] + row_n), ]
    z <- as.matrix(as.data.frame(lapply(z, as.numeric)))
    z <- z[nrow(z):1, ]
    pressure_array[, , i] <- z / (sens_area * 1000)
  }

  # return formatted data
  return(list(pressure_array = pressure_array, pressure_system = "footscan",
              sens_size = sens_size,
              time = time, masks = NULL, events = NULL))
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
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'  }
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- pressure_interp(pressure_data, interp_to = 101)
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
  time_seq <- seq(0, by = pressure_data[[4]], length.out = dims[3])
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
#' @param pressure_data List
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
#'   \item pressure_system. String defining pressure system
#'   \item sens_size. Numeric vector with the dimensions of the sensors
#'   \item time. Numeric value for time between measurements
#'   \item masks. List
#'   \item events. List
#'   }
#' @examples
#' \dontrun{
#' pedar_data <- system.file("extdata", "pedar_example.asc", package = "pressuRe")
#' pressure_data <- load_pedar(pedar_data)
#' pressure_data <- select_steps(pressure_data)
#' }
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab ggtitle
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom utils menu
#' @export

select_steps <- function (pressure_data, threshold_R = 30,
                          threshold_L = 30, min_frames = 10,
                          steps_Rn = 5, steps_Ln = 5, skip = 2) {
  # set up global variables
  frame <- NULL

  # check session is interactive
  if (interactive() == FALSE)
    stop("user needs to select suitable steps")

  # check this is pedar (or other suitable) data
  if (!(pressure_data[[2]] == "pedar" || pressure_data[[2]] == "tekscan"))
    stop("data should be from pedar or f-scan")

  # make force vectors
  force_R <- force_pedar(pressure_data)[,1]
  force_L <- force_pedar(pressure_data)[,2]

  # Adjust thresholds to avoid errors
  threshold_R <- threshold_R + 0.01
  threshold_L <- threshold_L + 0.01

  # throw error if threshold is less than min of trial
  min_R <- min(force_R)
  min_L <- min(force_L)
  if (min_R > threshold_R | min_L > threshold_L)
    stop("threshold is less than minimum force recorded in trial")

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
    Sys.sleep(1.5)

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
    Sys.sleep(1.5)

    # get user to approve or reject step
    resp <- menu(c("Y", "N"),
                 title = "Do you want to keep (Y) or discard (N) this step?")
    include_stps_L <- c(include_stps_L, resp)

    # check if number of steps has been reached
    if (sum(include_stps_L == 1) >= steps_Ln) {break}
  }

  # make events df
  len_R <- length(which(include_stps_R == 1))
  len_L <- length(which(include_stps_L == 1))
  event_df <- data.frame(side = c(rep("RIGHT", length.out = len_R),
                                  rep("LEFT", length.out = len_L)),
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
#' @param pressure_data List. First item should be a 3D array covering each
#' timepoint of the measurement. z dimension represents time
#' @return String. "LEFT" or "RIGHT"
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' auto_detect_side(pressure_data)
#' @importFrom sf st_polygon st_as_sf st_convex_hull st_combine
#' st_intersection st_area
#' @importFrom stats dist
#' @export

auto_detect_side <- function(pressure_data) {
  # throw error if pedar data
  if (pressure_data[[2]] == "pedar")
    stop("This function does not work for pedar data")

  # max pressure footprint
  sc_df <- sensor_2_polygon(pressure_data, output = "df")[, c(1, 2)]
  sc_df <- sc_df[!duplicated(sc_df), ]

  # Bounding box
  mbb <- getMinBBox(as.matrix(sc_df))
  side1 <- mbb[c(1, 2), ]
  side2 <- mbb[c(3, 4), ]

  # get midpoints
  midpoint_top <- (side1[2, ] + side2[2, ]) / 2
  midpoint_bottom <- (side1[1, ] + side2[1, ]) / 2

  # make lat and med box
  side1_pts <- rbind(side1, midpoint_top, midpoint_bottom, side1[1, ])
  side2_pts <- rbind(side2, midpoint_top, midpoint_bottom, side2[1, ])
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

#' @title Whole pressure curve
#' @description Generates vectors with option to plot for force, peak/mean
#' pressure and area for complete measurement. Useful for checking data
#' @param pressure_data List. A 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param variable String. "peak_pressure", "force", or "area"
#' @param side For insole data only
#' @param threshold Numeric. Threshold value for sensor to be considered active.
#' Currently only applies to insole data
#' @param plot Logical. If TRUE also plots data as line curve
#' @return Numeric vector containing variable values
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' whole_pressure_curve(pressure_data, variable = "force", plot = TRUE)
#' @importFrom ggplot2 aes ggplot geom_line theme_bw xlab ylab
#' @export

whole_pressure_curve <- function(pressure_data, variable, side, threshold = 10,
                                 plot = FALSE) {
  # set global variables
  time <- value <- NULL

  # check input
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")
  if(pressure_data[[2]] == "pedar" & missing(side) == TRUE)
    stop("pedar data needs to have side defined")

  # create empty vector to store variable
  values <- rep(NA, times = dim(pressure_data[[1]])[3])

  # force
  if (variable == "force") {
    if (pressure_data[[2]] == "pedar") {
        if (side == "RIGHT") {
          values <- force_pedar(pressure_data, variable = "force")[, 1]
        }
        if (side == "LEFT") {
          values <- force_pedar(pressure_data, variable = "force")[, 2]
        }
    } else {
      sens_area <- pressure_data[[3]][1] * pressure_data[[3]][2]
      force_array <- pressure_data[[1]] * sens_area * 1000

      # find total force for each frame and store in vector
      for (i in 1:dim(force_array)[3]) {
        values[i] <- sum(force_array[, , i])
      }
    }
  }

  # peak pressure
  if (variable == "peak_pressure") {
    if (pressure_data[[2]] == "pedar") {
      if (side == "RIGHT") {
        mat_r <- pressure_data[[1]][1, ,]
        values <- apply(mat_r, 2, max, na.rm = TRUE)
      }
      if (side == "LEFT") {
        mat_r <- pressure_data[[1]][2, ,]
        values <- apply(mat_r, 2, max, na.rm = TRUE)
      }
    } else {
      for (i in 1:dim(pressure_data[[1]])[3]) {
        values[i] <- max(pressure_data[[1]][, , i])
      }
    }
  }

  # area
  if (variable == "area") {
    if (pressure_data[[2]] == "pedar") {
      if (side == "RIGHT") {
        values <- force_pedar(pressure_data, variable = "area", threshold)[, 1]
      }
      if (side == "LEFT") {
        values <- force_pedar(pressure_data, variable = "area", threshold)[, 2]
      }
    } else {
      # sensor size
      sen_size <- pressure_data[[3]][1] * pressure_data[[3]][2]

      # find active area for each frame and store in vector
      for (i in 1:dim(pressure_data[[1]])[3]) {
        values[i] <- (sum(pressure_data[[1]][, , i] > 0)) * sen_size
      }
    }
  }

  # plot, if required
  if (plot == TRUE) {
    # make df
    variable_df <- data.frame(value = values,
                              time = seq(0, by = pressure_data[[4]],
                                         length.out = length(values)))

    # plot
    g <- ggplot(variable_df, aes(x = time, y = value))
    g <- g + geom_line()
    g <- g + theme_bw()
    g <- g + xlab("time (s)") + ylab(variable)
    print(g)
  }

  # return
  return(values)
}



# =============================================================================

#' @title Center of pressure
#' @description Generates center of pressure coordinates
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement. z dimension represents time
#' @return Data frame with x and y coordinates of COP throughout trial
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' cop(pressure_data)
#' @export

cop <- function(pressure_data) {
  # check input
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure_frames input must contain an array")

  # check not pedar
  if (pressure_data[[2]] == "pedar")
    stop("pedar data currently not supported for this function")

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
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' footprint(pressure_data, plot = TRUE)
#' @export

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
#' @param step_n Integer. Step number to plot (only for insole data)
#' @param smooth Logical. If TRUE, plot will interpolate between sensors to
#' increase data density
#' @param plot_COP Logical. If TRUE, overlay COP data on plot. Default = FALSE
#' @param plot_outline Logical. If TRUE, overlay convex hull outline on plot
#' @param plot_colors String. "default": novel color scheme; "custom": user
#' supplied
#' @param break_values Vector. If plot_colors is "custom", values to split
#' colors at
#' @param break_colors Vector. If plot_colors is "custom", colors to use.
#' Should be one shorter than break_values
#' @param plot Logical. If TRUE, plot will be displayed
#' @param legend Logical. If TRUE, legend will be added to plot
#' @return ggplot plot object
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' plot_pressure(pressure_data, variable = "max", plot_COP = FALSE)
#' plot_pressure(pressure_data, variable = "frame", frame = 20,
#'               plot_colors = "custom", break_values = c(100, 200, 300, 750),
#'               break_colors = c("blue", "green", "yellow", "red", "pink"))
#' @importFrom ggplot2 ggplot aes geom_raster geom_polygon scale_fill_manual
#' theme geom_point element_rect binned_scale unit
#' @importFrom scales manual_pal
#' @export

plot_pressure <- function(pressure_data, variable = "max", smooth = FALSE, frame,
                          step_n = "max", plot_COP = FALSE, plot_outline = FALSE,
                          plot_colors = "default", break_values, break_colors,
                          plot = TRUE, legend = TRUE) {
  # set global variables
  x <- y <- id <- cols <- x_coord <- y_coord <- value <- NULL

  # check inputs
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure data must contain array")

  # if pedar
  if (pressure_data[[2]] == "pedar") {
    if (step_n == "max") {pedar_var <- "max"; step_no <- NA}
    if (is.numeric(step_n)) {pedar_var <- "step_max"; step_no <- step_n}
    cor <- plot_pedar(pressure_data, pedar_var, step_no, foot_side = "both")
    x_lim <- max(cor$x)
    y_lim <- max(cor$y)
    legend_spacing <- (cor$x[2] - cor$x[1]) * 10
  } else {
    # max size of df
    fp_sens <- sensor_2_polygon(pressure_data, pressure_image = "all_active",
                                output = "df")
    x_lim <- max(fp_sens$x)
    y_lim <- max(fp_sens$y)

    # footprint
    fp <- footprint(pressure_data, variable = variable, frame)

    # get footprint vector
    fp <- as.vector(fp)
    fp <- fp[fp > 0]

    # generate coordinates for each sensor
    sens_poly <- sensor_2_polygon(pressure_data, pressure_image = variable,
                                  frame, output = "df")

    # combine with pressure values
    ids <- c(1:length(as.vector(fp)))
    vals <- data.frame(id = ids, value = as.vector(fp))

    # merge value and coordinate frames
    cor <- merge(sens_poly, vals, by = c("id"))

    # add colors
    cor <- generate_colors(cor, col_type = plot_colors, break_values,
                           break_colors)

    # plot
    legend_spacing <- pressure_data[[3]][1] * 100
  }

  if (plot_colors == "default") {
    break_colors <- c("grey","lightblue", "darkblue","green","yellow",
                      "red", "pink")
    break_values <- c(0, 40, 60, 100, 150, 220, 300)
  }

  # legend range
  range_max <- max(footprint(pressure_data))

  # plot
  g <- ggplot()
  g <- g + geom_polygon(data = cor, aes(x = x, y = y, group = id, fill = value))
  g <- g + binned_scale("fill", "foo",
                        binned_pal(manual_pal(break_colors)),
                        guide = "coloursteps", breaks = break_values,
                        #limits = c(0, max(cor$value)), show.limits = FALSE,
                        limits = c(0, range_max), show.limits = FALSE,
                        name = "Pressure (kPa)")
  g <- g + scale_x_continuous(expand = c(0, 0), limits = c(0, x_lim))
  g <- g + scale_y_continuous(expand = c(0, 0), limits = c(0, y_lim))
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
  g <- g + theme_void()
  g <- g + theme(panel.background = element_rect(fill = "white",
                                                 colour = "white"),
                 legend.box.spacing = unit(legend_spacing, "cm"))
  if (legend == FALSE) {g <- g + theme(legend.position = "none")}

  # display plot if requested
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
#' @param pressure_image String. "max" = footprint of maximum sensors. "frame"
#' = an individual frame
#' @param frame Integer. Frame number to use
#' @return Polygon. sfg object representing convex hull outline
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_outline(pressure_data)
#' pressure_outline(pressure_data, frame = 50)
#' @importFrom magrittr "%>%"
#' @importFrom sf st_as_sf st_convex_hull st_combine st_buffer
#' @export

pressure_outline <- function(pressure_data, pressure_image = "max", frame) {
  # check inputs
  if (is.array(pressure_data[[1]]) == FALSE)
    stop("pressure data must contain array")

  # check that this is not pedar data
  if (pressure_data[[2]] == "pedar")
    stop("data cannot be from pedar")

  # determine active sensor coordinates
  if (pressure_image == "max") {
    coords <- sensor_coords(pressure_data, pressure_image = "max")
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
#' @param plot_colors String
#' @param fps Numeric. Number of frames per second in animation
#' @param dpi Numeric. Resolution of gif
#' @param file_name Name (inlcuding path) of export file
#' @param preview Logical. Whether to play the animation
#' @return Animation in gif format
#' @examples
#' \dontrun{
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' animate_pressure(pressure_data, fps = 10, file_name = "testgif.gif")
#' }
#' @importFrom stringr str_ends str_pad
#' @importFrom magick image_graph image_animate image_write image_info
#' image_read
#' @importFrom ggplot2 ggplot aes geom_polygon scale_x_continuous
#' scale_y_continuous coord_fixed theme_void xlim ylim ggsave
#' @importFrom grDevices dev.off
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export

animate_pressure <- function(pressure_data, plot_colors = "default", fps,
                             dpi = 96, file_name, preview = FALSE) {
  # parameter check
  if (str_ends(file_name, ".gif") == FALSE)
    stop("filename must end in .gif")

  # max size of df
  fp_max <- sensor_2_polygon(pressure_data, pressure_image = "all_active",
                             output = "df")
  x_lim <- max(fp_max$x)
  y_lim <- max(fp_max$y)

  # set up temp directory
  temp_dir <- tempdir()

  # number of frames
  n_frames <- dim(pressure_data[[1]])[3]

  # plot
  img_fns <- rep(NA, length.out = n_frames)
  pb <- txtProgressBar(min = 1, max = n_frames, style = 3)
  print("processing images")
  for (i in 1:n_frames) {
    # make plot
    g <- plot_pressure(pressure_data, variable = "frame", frame = i,
                       plot_colors = plot_colors, plot = FALSE)
    img_fn <- paste0(temp_dir, "/img", str_pad(i, nchar(n_frames),
                                               pad = "0"), ".png")
    ggsave(img_fn, g, width = 6.45, height = 5.77, dpi = dpi)
    img_fns[i] <- img_fn
    setTxtProgressBar(pb, i)
  }

  # update progress
  close(pb)
  print("generating and saving animation")

  # load images back in
  allInfo <- image_info(image_read(img_fns))
  images <- image_read(img_fns)

  # create animation
  animation <- magick::image_animate(images, fps = fps)#, optimize = TRUE)

  # save animation
  magick::image_write(animation, file_name)
  file.remove(img_fns)
  invisible(gc())
}


# =============================================================================

#' @title automask footprint
#' @description Automatically create mask of footprint data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement. z dimension represents time
#' @param foot_side String. "RIGHT", "LEFT", or "auto". Auto uses
#' auto_detect_side function
#' @param mask_scheme List.
#' @param plot Logical. Whether to play the animation
#' @return List. Contains polygon with each mask
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- automask(pressure_data, foot_side = "auto", plot = TRUE)
#' @importFrom zoo rollapply
#' @importFrom sf st_union st_difference
#' @export

automask <- function(pressure_data, foot_side = "auto", mask_scheme,
                     plot = FALSE) {
  # check data isn't from pedar
  if (pressure_data[[2]] == "pedar")
    stop("automask doesn't work with pedar data")

  # global variables
  x <- y <- mask <- x_coord <- y_coord <- me <- heel_cut_dist <-
    mfoot_cut_prox <- ffoot_cut_prox <- NULL

  # Find footprint (max)
  max_df <- footprint(pressure_data)

  # coordinates
  sens_coords <- sensor_coords(pressure_data)

  # side
  if (foot_side == "auto") {
    side <- auto_detect_side(pressure_data)
  } else {
    side <- foot_side
  }

  # Define minimum bounding box
  sens_coords_m <- as.matrix(sens_coords)
  mbb <- getMinBBox(sens_coords_m)

  # Define convex hull, expanding to include all sensors
  df_sf <- sens_coords %>%
    st_as_sf(coords = c( "x_coord", "y_coord" ))
  fp_chull <- st_convex_hull(st_union(df_sf))
  fp_chull <- st_buffer(fp_chull, pressure_data[[3]][1])

  # Simple 3 mask (forefoot, midfoot, and hindfoot)
  ## Bounding box sides
  side1 <- mbb[c(1:2), ]
  side2 <- mbb[c(3:4), ]

  ## Get distal 27% and 55% lines that divide into masks
  side1_27 <- side1[1, ] + ((side1[2, ] - side1[1, ]) * 0.27)
  side2_27 <- side2[1, ] + ((side2[2, ] - side2[1, ]) * 0.27)
  side1_55 <- side1[1, ] + ((side1[2, ] - side1[1, ]) * 0.55)
  side2_55 <- side2[1, ] + ((side2[2, ] - side2[1, ]) * 0.55)
  line_27_df <- rbind(side1_27, side2_27)
  line_55_df <- rbind(side1_55, side2_55)
  line_27 <- st_linestring(as.matrix(line_27_df))
  line_27 <- st_extend_line(line_27, 1)
  line_27_dist_poly <- st_line2polygon(line_27, 1, "+Y")
  line_27_prox_poly <- st_line2polygon(line_27, 1, "-Y")
  line_55 <- st_linestring(as.matrix(line_55_df))
  line_55 <- st_extend_line(line_55, 1)
  line_55_dist_poly <- st_line2polygon(line_55, 1, "+Y")
  line_55_prox_poly <- st_line2polygon(line_55, 1, "-Y")

  ## forefoot, midfoot, and hindfoot masks
  hf_mask <- st_difference(fp_chull, line_27_dist_poly)
  mf_mask <- st_difference(fp_chull, line_55_dist_poly)
  mf_mask <- st_difference(mf_mask, line_27_prox_poly)
  ff_mask <- st_difference(fp_chull, line_55_prox_poly)

  # toe masks
  ## toe line
  toe_line_mat <- toe_line(pressure_data)

  ## toe cut polygons
  toe_poly_prox <- st_line2polygon(as.matrix(toe_line_mat), 1, "-Y")
  toe_poly_dist <- st_line2polygon(as.matrix(toe_line_mat), 1, "+Y")

  ## make masks
  toe_mask <- st_difference(fp_chull, toe_poly_prox)
  ff_mask <- st_difference(ff_mask, toe_poly_dist)

  # Define angles for dividing lines between metatarsals
  ## get edge lines
  edges <- edge_lines(pressure_data, side)

  ## angle between lines
  med_pts <- st_coordinates(edges[[1]])
  lat_pts <- st_coordinates(edges[[2]])
  med_line_angle <- (atan((med_pts[2, 1] - med_pts[1, 1]) /
                                 (med_pts[2, 2] - med_pts[1, 2]))) * 180 / pi
  lat_line_angle <- (atan((lat_pts[2, 1] - lat_pts[1, 1]) /
                                 (lat_pts[2, 2] - lat_pts[1, 2]))) * 180 / pi
  alpha <- lat_line_angle - med_line_angle


  ## intersection point
  med_lat_int <- st_intersection(edges[[1]], edges[[2]])

  ## rotation angles for 2-4 lines
  MT_hx_alpha <- (0.33 * alpha)
  MT_12_alpha <- (0.3  * alpha)
  MT_23_alpha <- (0.47 * alpha)
  MT_34_alpha <- (0.64 * alpha)
  MT_45_alpha <- (0.81 * alpha)

  ## polys for met and hal cuts
  if (side == "RIGHT") {lat_dir <- "+X"; med_dir <- "-X"}
  if (side == "LEFT") {lat_dir <- "-X"; med_dir <- "+X"}
  MT_hx_line <- rot_line(edges[[1]], MT_hx_alpha, med_lat_int)
  MT_hx_poly_lat <- st_line2polygon(st_coordinates(MT_hx_line)[, 1:2], 1, lat_dir)
  MT_hx_poly_med <- st_line2polygon(st_coordinates(MT_hx_line)[, 1:2], 1, med_dir)
  MT_12_line <- rot_line(edges[[1]], MT_12_alpha, med_lat_int)
  MT_12_poly_lat <- st_line2polygon(st_coordinates(MT_12_line)[, 1:2], 1, lat_dir)
  MT_12_poly_med <- st_line2polygon(st_coordinates(MT_12_line)[, 1:2], 1, med_dir)
  MT_23_line <- rot_line(edges[[1]], MT_23_alpha, med_lat_int)
  MT_23_poly_lat <- st_line2polygon(st_coordinates(MT_23_line)[, 1:2], 1, lat_dir)
  MT_23_poly_med <- st_line2polygon(st_coordinates(MT_23_line)[, 1:2], 1, med_dir)
  MT_34_line <- rot_line(edges[[1]], MT_34_alpha, med_lat_int)
  MT_34_poly_lat <- st_line2polygon(st_coordinates(MT_34_line)[, 1:2], 1, lat_dir)
  MT_34_poly_med <- st_line2polygon(st_coordinates(MT_34_line)[, 1:2], 1, med_dir)
  MT_45_line <- rot_line(edges[[1]], MT_45_alpha, med_lat_int)
  MT_45_poly_lat <- st_line2polygon(st_coordinates(MT_45_line)[, 1:2], 1, lat_dir)
  MT_45_poly_med <- st_line2polygon(st_coordinates(MT_45_line)[, 1:2], 1, med_dir)

  ## make met masks
  hal_mask <- st_difference(toe_mask, MT_hx_poly_lat)
  l_toe_mask <- st_difference(toe_mask, MT_hx_poly_med)
  MT1_mask <- st_difference(ff_mask, MT_12_poly_lat)
  MT2_mask <- st_difference(ff_mask, MT_23_poly_lat)
  MT2_mask <- st_difference(MT2_mask, MT_12_poly_med)
  MT3_mask <- st_difference(ff_mask, MT_34_poly_lat)
  MT3_mask <- st_difference(MT3_mask, MT_23_poly_med)
  MT4_mask <- st_difference(ff_mask, MT_45_poly_lat)
  MT4_mask <- st_difference(MT4_mask, MT_34_poly_med)
  MT5_mask <- st_difference(ff_mask, MT_45_poly_med)

  # Make mask list
  mask_list <- list(heel_mask = hf_mask,
                    midfoot_mask = mf_mask,
                    forefoot_mask = ff_mask,
                    hal_mask = hal_mask,
                    l_toe_mask = l_toe_mask,
                    MT1_mask = MT1_mask,
                    MT2_mask = MT2_mask,
                    MT3_mask = MT3_mask,
                    MT4_mask = MT4_mask,
                    MT5_mask = MT5_mask)


  # Plot Footprint and masks if required
  if (plot == TRUE) {
    # make masks into df format
    mask_df <- masks_2_df(mask_list)

    # Plot footprint and masks
    g <- plot_pressure(pressure_data, "max", plot = FALSE)
    g <- g + scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 0.15))
    g <- g + scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 0.30))
    g <- g + geom_path(data = mask_df, aes(x = x, y = y, group = mask),
                       color = "red", linewidth = 1)
    print(g)
  }

  # Return masks for analysis
  pressure_data[[5]] <- mask_list
  return(pressure_data)
}


# =============================================================================

#' @title Create mask
#' @description Manually create mask region
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @param image String. "max" = footprint of maximum sensors. "mean"
#'   average value of sensors over time (usually for static analyses)
#' @param n_verts Numeric. Number of vertices in mask
#' @param n_masks Numeric. Number of masks to add
#' @param threshold Numeric. Distance between adjacent mask vertices before
#' sharing vertex coordinates
#' @param plot_existing_mask Logical. Show exisiting masks
#' @param preview Logical. Show new maks on pressure image
#' @return List New mask is added to the relevant A 3D array covering each
#' timepoint of the measurement for the selected region. z dimension represents
#' time
#' @examples
#' \dontrun{
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- create_mask(pressure_data, 4)
#' }
#' @importFrom grDevices x11
#' @importFrom ggmap gglocator
#' @importFrom ggplot2 aes geom_path
#' @importFrom sf st_polygon
#' @export

create_mask <- function(pressure_data, n_verts = 4, n_masks = 1,
                        threshold = 0.005, plot_existing_mask = TRUE,
                        image = "max", preview = TRUE) {
  # global variables
  x <- y <- NULL

  # check session is interactive
  if (interactive() == FALSE)
    stop("user needs to select mask vertices")

  # plot existing masks or just footprint
  if (plot_existing_mask == TRUE) {
    if (length(pressure_data[[5]]) > 0) {
    g <- plot_masks(pressure_data)
    }
  } else {
    grDevices::x11()
    g <- plot_pressure(pressure_data, image, plot = FALSE)
    g <- g + scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 0.15))
    g <- g + scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 0.30))
    print(g)
  }

  # empty data frame to store vertex coordinates
  mask_vertices <- data.frame(x = double(), y = double())

  for(mask_n in 1:n_masks){
    # interactively select area
    message("Select mask corners")
    mask <- gglocator(n_verts)
    valid_mask <- which(dist(mask) < threshold)

    while(length(valid_mask) > 0){
      message("Some vertices are closer than your designated threshold.
              Reselect mask corners")
      mask <- gglocator(n_verts)
      valid_mask <- which(dist(mask) < threshold)
    }

    # allow closely selected mask vertices to share a vertex
    mask_vertices[(nrow(mask_vertices) + 1):
                    (nrow(mask_vertices) + n_verts), ] <- mask
    if (mask_n > 1) {
      v_distance <- as.matrix(dist(mask_vertices))
      close_v <- which((v_distance < threshold) & (v_distance > 0),
                       arr.ind = TRUE)
      for (change_v in 1:nrow(close_v) / 2) {
        mask[close_v[change_v, 1] - ((mask_n - 1) * n_verts), 1] <-
          mask_vertices[close_v[change_v, 2], 1]
        mask[close_v[change_v, 1] - ((mask_n - 1) * n_verts), 2] <-
          mask_vertices[close_v[change_v, 2], 2]
        mask_vertices[close_v[change_v, 1], ] <-
          mask_vertices[close_v[change_v, 2], ]
      }
    }

    mask <- data.frame(x = mask[, 1], y = mask[, 2])
    mask <- mask[c(1:nrow(mask), 1), ]

    # preview
    if (preview == TRUE) {
      g <- g + geom_path(data = mask, aes(x, y), color = "blue")
      print(g)
    }

    # mask to polygon
    mask_pol <- st_polygon(list(as.matrix(mask)))

    # update mask list
    mask_list <- pressure_data[[5]]
    mask_list[[length(mask_list) + 1]] <- mask_pol
    pressure_data[[5]] <- mask_list
  }

  # return pressure frames for selected area
  return(pressure_data)
}


# =============================================================================

#' @title Edit mask
#' @description Manually edit mask, also simplifies automask vertices to the
#' threshold value
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @param n_edit Numeric. Number of vertices to edit
#' @param threshold Numeric. Distance between point clicked and vertex that is
#' selected (smaller thresholds are recommended when editing automasks becasue
#' they contain more vertices in the polygon)
#' @param edit_list List. Mask numbers that want to be edited. (Default is to
#' load all masks so that adjacent masks with shared coordinates are modified
#' together)
#' @return List. Edited mask is added to the relevant A 3D array covering each
#' timepoint of the measurement for the selected region. z dimension represents
#' time
#' @examples
#' \dontrun{
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- automask(pressure_data, foot_side = "auto", plot = TRUE)
#' edit_mask(pressure_data)
#' }
#' @importFrom grDevices rainbow
#' @export

edit_mask <- function(pressure_data, n_edit, threshold = 0.002,
                      edit_list = seq(1,length(pressure_data[[5]]))) {
  # check session is interactive
  if (interactive() == FALSE) {
    stop("user needs to select mask vertices")
  }

  # global variables
  X <- Y <- NULL

  if(n_edit > 0){
    # plot original mask data
    g <- plot_masks(pressure_data)

    # compile all mask vertices
    mask_org <- pressure_data[[5]]
    mask_vertices_edit <- data.frame(x = double(), y = double())

    for (n_mask in edit_list) {
      mask_vertices <- data.frame(st_coordinates(mask_org[[n_mask]])[,1:2])
      mask_vertices_edit[(nrow(mask_vertices_edit) + 1):
                           (nrow(mask_vertices) + nrow(mask_vertices_edit)), ] <- mask_vertices
    }

    # select point near edit vertex and select new value
    color_v <- rainbow(n_edit + 1)
    for (edit in seq(1, n_edit)){
      message("Select mask vertex to edit, then it's new location")
      mask <- gglocator(2)
      mask_vertices_edit[nrow(mask_vertices_edit) + 1, ] <- mask[1, ]
      edit_distance <- as.matrix(dist(mask_vertices_edit))
      close_v <-
        which((edit_distance[nrow(mask_vertices_edit),] < threshold) &
                (edit_distance[nrow(mask_vertices_edit),] > 0), arr.ind = TRUE)

      while (length(close_v) == 0) {
        message("The point selected isn't close to any existing vertex, try again")
        mask <- gglocator(2)
        mask_vertices_edit[nrow(mask_vertices_edit) + 1, ] <- mask[1,]
        edit_distance <- as.matrix(dist(mask_vertices_edit))
        close_v <-
          which((edit_distance[nrow(mask_vertices_edit),] < threshold) &
                  (edit_distance[nrow(mask_vertices_edit),] > 0), arr.ind = TRUE)
      }

      for (n_mask in edit_list) {
        mask_vertices <- data.frame(st_coordinates(mask_org[[n_mask]])[, 1:2])
        new_mask_v <- mask_vertices
        mask_vertices[nrow(mask_vertices) + 1, ] <- mask[1, ]
        v_distance <- as.matrix(dist(mask_vertices))
        v_distance[upper.tri(v_distance)] <- 0.0
        close_selected <- which((v_distance[nrow(mask_vertices), ] < threshold)
                                 & (v_distance[nrow(mask_vertices), ] > 0),
                                arr.ind = TRUE)

        if (length(close_selected) > 0){
          message(n_mask, ":", close_selected)

          if(close_selected[1] == 1){
            new_mask_v[c(1, nrow(new_mask_v)), ] <- mask[2, ]
            close_selected <- head(close_selected[-1], -1)
          }
          if(length(close_selected) > 0){
            new_mask_v[close_selected[1], ] <- mask[2, ]
            if(length(close_selected) > 1){
              new_mask_v <- new_mask_v[-tail(close_selected, -1), ]
            }
          }

          g <- g + geom_path(data = new_mask_v, aes(X, Y), color = color_v[edit + 1])
          mask_org[[n_mask]] <-
            st_polygon(list(as.matrix(new_mask_v)))
        }
      }

      print(g)
    }

    # return pressure frames for selected area
    pressure_data[[5]] <- mask_org
    return(pressure_data)
  } else {
    message("You must edit at least one vertex, change the value for 'n_edit'")
  }
}


# =============================================================================

#' @title Pedar mask template
#' @description Add a Pedar mask template using 3 in current literature
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @param mask_type String.
#' "mask1" splits the insole into 4 regions using sensel boundaries:
#' rearfoot, midfoot, forefoot, and toes-- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9470545/
#' "mask2" splits the insole into 4 regions using percentages:
#' rearfoot, forefoot, hallux, and lesser toes-- https://jfootankleres.biomedcentral.com/articles/10.1186/1757-1146-7-18
#' "mask3" splits the foot into 9 regions using sensel boundaries:
#'  medial rearfoot, lateral rearfoot, medial midfoot, lateral midfoot, MTPJ1,
#'  MTPJ2-3, MTPJ4-5, hallux, and lesser toes-- https://jfootankleres.biomedcentral.com/articles/10.1186/1757-1146-7-20
#' @param plot Logical. Plot masks if TRUE
#' @return List. Masks are added to the relevant A 3D array covering each
#' timepoint of the measurement for the selected region. z dimension represents
#' time
#' @examples
#' pedar_data <- system.file("extdata", "pedar_example.asc", package = "pressuRe")
#' pressure_data <- load_pedar(pedar_data)
#' pressure_data <- pedar_mask(pressure_data, "mask1")
#' @importFrom sf st_union st_difference st_bbox
#' @export

pedar_mask <- function(pressure_data, mask_type, plot = TRUE) {

  # set global variables
  pedar_insole_grid <- x <- y <- id <- NULL

  # check this is pedar (or other suitable) data
  if (pressure_data[[2]] != "pedar")
    stop("data should be from pedar")

  # load pedar coords
  load("data/pedar_insole_grid.rda")

  # add coordinates to df
  xs <- c()
  for (i in c(101:199, 1:99)) {
    for (j in c(1, 3, 5, 7)) {xs = append(xs, pedar_insole_grid[i, j])}
  }
  ys <- c()
  for (i in c(101:199, 1:99)) {
    for (j in c(2, 4, 6, 8)) {ys = append(ys, pedar_insole_grid[i, j])}
  }

  # make ids
  ids <- c()
  for (i in 1:99) {ids = append(ids, i)} #left
  for (i in 1:99) {ids = append(ids, i)} #right

  position <- data.frame(id = rep(ids, each = 4), x = xs, y = ys)

  # defines the sensels or percentages for masking
  if(mask_type == "mask1"){
    rearfoot_L <- pedar_polygon(position, 1:26, "L")
    midfoot_L <-pedar_polygon(position, 27:54, "L")
    forefoot_L <- pedar_polygon(position, 55:82, "L")
    toes_L <- pedar_polygon(position, 83:99, "L")

    rearfoot_R <- pedar_polygon(position, 1:26, "R")
    midfoot_R <-pedar_polygon(position, 27:54, "R")
    forefoot_R <- pedar_polygon(position, 55:82, "R")
    toes_R <- pedar_polygon(position, 83:99, "R")

    mask_list <- list(L_rearfoot_mask = rearfoot_L,
                        L_midfoot_mask = midfoot_L,
                        L_forefoot_mask = forefoot_L,
                        L_toes_mask = toes_L,
                        R_rearfoot_mask = rearfoot_R,
                        R_midfoot_mask = midfoot_R,
                        R_forefoot_mask = forefoot_R,
                        R_toes_mask = toes_R)

      pressure_data[[5]] <- mask_list

  } else if (mask_type == "mask2"){
    bbox_L <- sf::st_bbox(pedar_polygon(position, 1:99, "L"))
    outline_mask <- pedar_polygon(position, 1:99, "L")

    rearfoot_line_Ly <- (bbox_L$ymax-bbox_L$ymin)/2 + bbox_L$ymin
    forefoot_line_Ly <- (bbox_L$ymax-bbox_L$ymin)*.85 + bbox_L$ymin
    hallux_line_Lx <- bbox_L$xmax - (bbox_L$xmax-bbox_L$xmin)*.35

    rearfoot_line_L <- data.frame(x_coord = c(bbox_L$xmin,bbox_L$xmax),
                                y_coord =  c(rearfoot_line_Ly,rearfoot_line_Ly))
    rearfoot_line_L <- st_extend_line(st_linestring(as.matrix(rearfoot_line_L)),1)
    rearfoot_line_L_dist <- st_line2polygon(rearfoot_line_L, 1, "+Y")
    rearfoot_L <- st_difference(outline_mask, rearfoot_line_L_dist)

    forefoot_line_L <- data.frame(x_coord = c(bbox_L$xmin,bbox_L$xmax),
                                  y_coord =  c(forefoot_line_Ly,forefoot_line_Ly))
    forefoot_line_L <- st_extend_line(st_linestring(as.matrix(forefoot_line_L)),1)
    forefoot_line_L_dist <- st_line2polygon(forefoot_line_L, 1, "+Y")
    forefoot_line_L_prox <- st_line2polygon(forefoot_line_L, 1, "-Y")
    forefoot_L <- st_difference(outline_mask, forefoot_line_L_dist)

    toes_L <- st_difference(outline_mask, forefoot_line_L_prox)
    hallux_line_L <- data.frame(x_coord = c(hallux_line_Lx,hallux_line_Lx),
                                  y_coord =  c(bbox_L$ymin,bbox_L$ymax))
    hallux_line_L <- st_extend_line(st_linestring(as.matrix(hallux_line_L)),1)
    hallux_line_L_lat <- st_line2polygon(hallux_line_L, 1, "-X")
    hallux_line_L_med <- st_line2polygon(hallux_line_L, 1, "+X")
    hallux_L <- st_difference(toes_L, hallux_line_L_lat)
    lesser_toes_L <- st_difference(toes_L, hallux_line_L_med)

    bbox_R <- sf::st_bbox(pedar_polygon(position, 1:99, "R"))
    outline_mask <- pedar_polygon(position, 1:99, "R")

    rearfoot_line_Ry <- (bbox_R$ymax-bbox_R$ymin)/2 + bbox_R$ymin
    forefoot_line_Ry <- (bbox_R$ymax-bbox_R$ymin)*.85 + bbox_R$ymin
    hallux_line_Rx <- bbox_R$xmin + (bbox_R$xmax-bbox_R$xmin)*.35

    rearfoot_line_R <- data.frame(x_coord = c(bbox_R$xmin,bbox_R$xmax),
                                  y_coord =  c(rearfoot_line_Ry,rearfoot_line_Ry))
    rearfoot_line_R <- st_extend_line(st_linestring(as.matrix(rearfoot_line_R)),1)
    rearfoot_line_R_dist <- st_line2polygon(rearfoot_line_R, 1, "+Y")
    rearfoot_R <- st_difference(outline_mask, rearfoot_line_R_dist)

    forefoot_line_R <- data.frame(x_coord = c(bbox_R$xmin,bbox_R$xmax),
                                  y_coord =  c(forefoot_line_Ry,forefoot_line_Ry))
    forefoot_line_R <- st_extend_line(st_linestring(as.matrix(forefoot_line_R)),1)
    forefoot_line_R_dist <- st_line2polygon(forefoot_line_R, 1, "+Y")
    forefoot_line_R_prox <- st_line2polygon(forefoot_line_R, 1, "-Y")
    forefoot_R <- st_difference(outline_mask, forefoot_line_R_dist)

    toes_R <- st_difference(outline_mask, forefoot_line_R_prox)
    hallux_line_R <- data.frame(x_coord = c(hallux_line_Rx,hallux_line_Rx),
                                y_coord =  c(bbox_R$ymin,bbox_R$ymax))
    hallux_line_R <- st_extend_line(st_linestring(as.matrix(hallux_line_R)),1)
    hallux_line_R_lat <- st_line2polygon(hallux_line_R, 1, "+X")
    hallux_line_R_med <- st_line2polygon(hallux_line_R, 1, "-X")
    hallux_R <- st_difference(toes_R, hallux_line_R_lat)
    lesser_toes_R <- st_difference(toes_R, hallux_line_R_med)

    mask_list <- list(L_rearfoot_mask = rearfoot_L,
                      L_forefoot_mask = forefoot_L,
                      L_hallux_mask = hallux_L,
                      L_lesser_toes_mask = lesser_toes_L,
                      R_rearfoot_mask = rearfoot_R,
                      R_forefoot_mask = forefoot_R,
                      R_hallux_mask = hallux_R,
                      R_lesser_toes_mask = lesser_toes_R)

    pressure_data[[5]] <- mask_list
  } else if (mask_type == "mask3"){
    med_rf_L <- st_union(pedar_polygon(position, c(1:2), "L"),
                         st_union(pedar_polygon(position, c(6:8), "L"),
                         st_union(pedar_polygon(position, c(13:15), "L"),
                                  pedar_polygon(position, c(20:22), "L"))))
    lat_rf_L <- st_union(pedar_polygon(position, c(3:5), "L"),
                         st_union(pedar_polygon(position, c(9:12), "L"),
                                  st_union(pedar_polygon(position, c(16:19), "L"),
                                  pedar_polygon(position, c(23:26), "L"))))
    med_mf_L <- st_union(pedar_polygon(position, c(27:29), "L"),
                         st_union(pedar_polygon(position, c(34:36), "L"),
                                  st_union(pedar_polygon(position, c(41:43), "L"),
                                           st_union(pedar_polygon(position, c(48:50), "L"),
                                 pedar_polygon(position, c(55:57), "L")))))
    lat_mf_L <- st_union(pedar_polygon(position, c(30:33), "L"),
                         st_union(pedar_polygon(position, c(37:40), "L"),
                                  st_union(pedar_polygon(position, c(44:47), "L"),
                                           st_union(pedar_polygon(position, c(51:54), "L"),
                                 pedar_polygon(position, c(58:59), "L")))))
    MTPJ1_L<- st_union(pedar_polygon(position, c(62:63), "L"),
                       st_union(pedar_polygon(position, c(69:70), "L"),
                     pedar_polygon(position, c(76:77), "L")))
    MTPJ23_L<- st_union(pedar_polygon(position, c(64:66), "L"),
                        st_union(pedar_polygon(position, c(71:73), "L"),
                     pedar_polygon(position, c(78:80), "L")))
    MTPJ45_L<- st_union(pedar_polygon(position, c(60:61), "L"),
                        st_union(pedar_polygon(position, c(67:68), "L"),
                        pedar_polygon(position, c(74:75), "L")))
    hallux_L<- st_union(pedar_polygon(position, c(83:84), "L"),
                        st_union(pedar_polygon(position, c(90:91), "L"),
                        pedar_polygon(position, c(96), "L")))
    lesser_toes_L <- st_union(pedar_polygon(position, c(81:82), "L"),
                              st_union(pedar_polygon(position, c(85:89), "L"),
                                       st_union(pedar_polygon(position, c(92:95), "L"),
                              pedar_polygon(position, c(97:99), "L"))))


    med_rf_R <- st_union(pedar_polygon(position, c(1:2), "R"),
                         st_union(pedar_polygon(position, c(6:8), "R"),
                                  st_union(pedar_polygon(position, c(13:15), "R"),
                                           pedar_polygon(position, c(20:22), "R"))))
    lat_rf_R <- st_union(pedar_polygon(position, c(3:5), "R"),
                         st_union(pedar_polygon(position, c(9:12), "R"),
                                  st_union(pedar_polygon(position, c(16:19), "R"),
                                           pedar_polygon(position, c(23:26), "R"))))
    med_mf_R <- st_union(pedar_polygon(position, c(27:29), "R"),
                         st_union(pedar_polygon(position, c(34:36), "R"),
                                  st_union(pedar_polygon(position, c(41:43), "R"),
                                           st_union(pedar_polygon(position, c(48:50), "R"),
                                                    pedar_polygon(position, c(55:57), "R")))))
    lat_mf_R <- st_union(pedar_polygon(position, c(30:33), "R"),
                         st_union(pedar_polygon(position, c(37:40), "R"),
                                  st_union(pedar_polygon(position, c(44:47), "R"),
                                           st_union(pedar_polygon(position, c(51:54), "R"),
                                                    pedar_polygon(position, c(58:59), "R")))))
    MTPJ1_R<- st_union(pedar_polygon(position, c(62:63), "R"),
                       st_union(pedar_polygon(position, c(69:70), "R"),
                                pedar_polygon(position, c(76:77), "R")))
    MTPJ23_R<- st_union(pedar_polygon(position, c(64:66), "R"),
                        st_union(pedar_polygon(position, c(71:73), "R"),
                                 pedar_polygon(position, c(78:80), "R")))
    MTPJ45_R<- st_union(pedar_polygon(position, c(60:61), "R"),
                        st_union(pedar_polygon(position, c(67:68), "R"),
                                 pedar_polygon(position, c(74:75), "R")))
    hallux_R<- st_union(pedar_polygon(position, c(83:84), "R"),
                        st_union(pedar_polygon(position, c(90:91), "R"),
                                 pedar_polygon(position, c(96), "R")))
    lesser_toes_R <- st_union(pedar_polygon(position, c(81:82), "R"),
                              st_union(pedar_polygon(position, c(85:89), "R"),
                                       st_union(pedar_polygon(position, c(92:95), "R"),
                                                pedar_polygon(position, c(97:99), "R"))))
    mask_list <- list(L_medial_rearfoot_mask = med_rf_L,
                      L_lateral_rearfoot_mask = lat_rf_L,
                      L_medial_midfoot_mask = med_mf_L,
                      L_lateral_midfoot_mask = lat_mf_L,
                      L_MTPJ1_mask = MTPJ1_L,
                      L_MTPJ23_mask = MTPJ23_L,
                      L_MTPJ45_mask = MTPJ45_L,
                      L_hallux_mask = hallux_L,
                      L_lesser_toes_mask = lesser_toes_L,
                      R_medial_rearfoot_mask = med_rf_R,
                      R_lateral_rearfoot_mask = lat_rf_R,
                      R_medial_midfoot_mask = med_mf_R,
                      R_lateral_midfoot_mask = lat_mf_R,
                      R_MTPJ1_mask = MTPJ1_R,
                      R_MTPJ23_mask = MTPJ23_R,
                      R_MTPJ45_mask = MTPJ45_R,
                      R_hallux_mask = hallux_R,
                      R_lesser_toes_mask = lesser_toes_R)

    pressure_data[[5]] <- mask_list
  } else {
    stop("Please select an existing mask.")
  }

  # plot
  if (plot == TRUE) {
    plot_masks(pressure_data)
  }

  # return
  return(pressure_data)
}


# =============================================================================

#' @title CPEI
#' @description Determine Center of Pressure Excursion Index (CPEI) for
#' footprint pressure data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @param foot_side String. "right" or "left". Required for automatic detection of
#'   points
#' @param plot_result Logical. Plots pressure image with COP and CPEI overlaid
#' @return Numeric. CPEI value
#' @examples
#' \dontrun{
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' cpei(pressure_data, foot_side = "auto", plot = TRUE)
#' }
#' @importFrom sf st_convex_hull st_linestring st_distance st_coordinates
#' @importFrom dplyr pull summarise
#' @export

cpei <- function(pressure_data, foot_side, plot_result = TRUE) {
  # check set up
  if (interactive() == FALSE)
    stop("we recommend user reviews each measurement")

  # global variables
  x <- y <- me <- X <- Y <- NULL

  # side
  if (foot_side == "auto") {
    side <- auto_detect_side(pressure_data)
  } else {
    side <- foot_side
  }

  # footprint coordinates
  sc_df <- sensor_2_polygon(pressure_data, output = "df")[, c(1, 2)]
  sc_df <- sc_df[!duplicated(sc_df), ]

  # determine medial line
  ## medial edge coords
  unq_y <- unique(sc_df$y)
  med_edge <- data.frame(x = rep(NA, length.out = length(unq_y)),
                         y = unq_y)
  for (i in 1:length(unq_y)) {
    if (side == "LEFT") {
      med_edge[i, 1] <- sc_df %>% filter(y == unq_y[i]) %>%
        summarise(me = max(x)) %>% pull(me)
    } else {
      med_edge[i, 1] <- sc_df %>% filter(y == unq_y[i]) %>%
        summarise(me = min(x)) %>% pull(me)
    }
  }

  ## convex hull of medial edge
  med_edge_sf <- med_edge %>% st_as_sf(coords = c("x", "y"))
  med_edge_chull <- st_convex_hull(st_combine(med_edge_sf))

  ## longest med edge line
  me_dis <- as.matrix(dist(st_coordinates(med_edge_chull)))
  me_dis <- me_dis[row(me_dis) == (col(me_dis) - 1)]
  me_max <- order(me_dis)[length(me_dis)]
  med_side <- st_coordinates(med_edge_chull)[c(me_max, me_max + 1), c(1, 2)]
  med_side_line <- st_linestring(med_side)
  med_side_line <- st_extend_line(med_side_line, 0.1)

  ## confirm no points to medial side of line


  # bounding box for whole foot
  ## Bounding box
  mbb <- getMinBBox(as.matrix(sc_df))
  side1 <- mbb[c(1:2), ]
  side2 <- mbb[c(3:4), ]

  ## Get distal trisection points and make line
  trisection_left <- side2[1, ] + ((side2[2, ] - side2[1, ]) * 0.667)
  trisection_right <- side1[1, ] + ((side1[2, ] - side1[1, ]) * 0.667)
  trisection_df <- data.frame(x = c(trisection_left[1], trisection_right[1]),
                              y = c(trisection_left[2], trisection_right[2]))
  tri_line <- st_linestring(as.matrix(trisection_df))
  tri_line <- st_extend_line(tri_line, 1)
  #tri_int <- st_intersection(st_cast(fp_chull, "MULTILINESTRING"), tri_line)

  # trisection of medial line
  med_edge_tri_pt <- st_intersection(med_side_line, tri_line)

  # perp med line
  med_slope <- (med_side_line[2, 2] - med_side_line[1, 2]) /
                 (med_side_line[2, 1] - med_side_line[1, 1])
  perp_med_line_slope <- -1 / med_slope
  tri_med_int <- med_edge_tri_pt[2] - (perp_med_line_slope * med_edge_tri_pt[1])
  perp_med_line <- rbind(med_edge_tri_pt, c(0, tri_med_int))
  perp_med_line <- st_extend_line(st_linestring(perp_med_line), 0.2)

  # ff width
  ## convex hull
  df.sf <- sc_df %>%
    st_as_sf(coords = c( "x", "y" ))
  fp_chull <- st_convex_hull(st_combine(df.sf))
  ff_lat_int <- st_intersection(fp_chull, perp_med_line)
  ff_lat_int <- st_coordinates(ff_lat_int)[, c(1, 2)]
  if (side == "RIGHT") {
    ff_lat_int <- ff_lat_int[which.max(ff_lat_int[, 1]), ]
  } else {
    ff_lat_int <- ff_lat_int[which.min(ff_lat_int[, 1]), ]
  }
  ff_width_pts <- rbind(ff_lat_int, st_coordinates(med_edge_tri_pt))
  ff_width <- as.numeric(unname(dist(ff_width_pts)))

  # cop_line
  cop_df <- cop(pressure_data)
  cop_sf <- cop_df %>%
    st_as_sf(coords = c("x_coord", "y_coord"))
  cop_chull <- st_convex_hull(st_combine(cop_sf))

  # cop_straight
  cop_dis <- as.matrix(dist(st_coordinates(cop_chull)))
  cop_dis <- cop_dis[row(cop_dis) == (col(cop_dis) - 1)]
  cop_max <- order(cop_dis)[length(cop_dis)]
  cop_side <- st_coordinates(cop_chull)[c(cop_max, cop_max + 1), c(1, 2)]
  cop_side_line <- st_linestring(cop_side)
  cop_side_line <- st_extend_line(cop_side_line, 0.1)

  # cross point cop and per_med_line
  cop1_per_med_pt <- st_intersection(st_linestring(as.matrix(cop_df)),
                                     perp_med_line)

  # cross point cop straight and per_med_line
  cop2_per_med_pt <- st_intersection(cop_side_line, perp_med_line)

  # CPE distance
  CPE <- as.numeric(unname(st_distance(cop1_per_med_pt, cop2_per_med_pt)))

  # calculate CPEI
  CPEI <- (CPE / ff_width) * 100

  # plot
  if (plot_result == TRUE) {
    # make CPEI dfs for plotting
    med_side_df <- as.data.frame(med_side)
    cop_side_df <- as.data.frame(cop_side)
    cpe_df <- as.data.frame(rbind(cop1_per_med_pt, cop2_per_med_pt))
    colnames(cpe_df) <- c("X", "Y")

    # plot
    g <- plot_pressure(pressure_data, plot_colors = "custom",
                       break_values = c(100, 750),
                       break_colors = c("lightgrey", "grey", "darkgrey"),
                       plot_COP = TRUE)
    g <- g + geom_line(data = med_side_df, aes(X, Y), linetype = "dashed",
                       color = "black", size = 1.5)
    g <- g + geom_line(data = cpe_df, aes(X, Y), color = "black",
                       size = 2)
    g <- g + geom_line(data = cop_side_df, aes(X, Y), colour = "blue",
                       alpha = 0.8)
    print(g)
  }

  # manual select function
  #manually_select <- function(n_points, mess) {
  #  message(mess)
  #  x <- gglocator(n_points)
  #  colnames(x) <- c("x_coord", "y_coord")
  #  return(x)
  #}

  ## check if automatic identification worked
  # plot max footprint with additional points
  #g <- plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)
  #g <- g + geom_point(data = heel_coord, aes (x = x_coord, y = y_coord), shape = 6)
  #g <- g + geom_point(data = toe_coord, aes (x = x_coord, y = y_coord), shape = 6)
  #g <- g + geom_point(data = start_point, aes (x = x_coord, y = y_coord), shape = 2)
  #g <- g + geom_point(data = end_point, aes (x = x_coord, y = y_coord), shape = 2)
  #g <- g + geom_line(data = m_bor, aes(x = x_coord, y = y_coord), colour = "purple")
  #g <- g + geom_line(data = l_bor, aes(x = x_coord, y = y_coord), colour = "orange")
  print(g)
  #auto_worked <- readline("Have points been correctly identified? c: manually select cop; a: manually select all")

  ## if automatic identification failed, redo manually
  #if (auto_worked == "a") {
  #  # plot footprint
  #  g <- plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)

    # select points
  #  m_bor <- manually_select(2, "select medial border: proximal point first, then distal")
  #  l_bor <- manually_select(2, "select lateral border: proximal point first, then distal")
  #  heel_coord <- manually_select(1, "select most proximal point of heel")
  #  toe_coord <- manually_select(1, "select most distal point of toes")
  #  cop_start <- manually_select(1, "select the most medial point near the start of the COP")
  #  cop_end <- manually_select(1, "select the most medial point near the end of the COP")
  #}

  #if (auto_worked == "c") {
  #  # plot footprint
  #  g <- plot_footprint(pressure_frames, plot_COP = TRUE, plot_outline = TRUE)

    # select points
  #  start_point <- manually_select(1, "select the most medial point near the start of the COP")
  #  end_point <- manually_select(1, "select the most medial point near the end of the COP")
  #}

  # return CPEI
  return(CPEI)
}


# =============================================================================

#' Analyze masked regions of pressure data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param partial_sensors Logical Defines how sensors that do not
#'   lie wholly within mask are dealt with. If FALSE, they will be excluded;
#'   if TRUE, for relevant variables their contribution will be weighted by the
#'   proportion of the sensor that falls within the mask border
#' @param variable String. Variable to be determined. "press_peak_sensor",
#'  "press_peak_sensor_ts", "press_peak_mask", "press_peak_mask_ts",
#'  "contact_area_peak", "contact_area_ts", "pti_1", "pti_2",
#'  "force_time_integral", "force_peak", "force_ts"
#' @param pressure_units String. Default "kPa". Other options: "MPa", "Ncm2"
#' (Newtons per square centimeter)
#' @param area_units String. Default is "cm2" (square centimeters). Other
#' options "m2" (square meters); "mm2" (square millimeters)
#' @return Data frame.
#' @examples
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- automask(pressure_data)
#' mask_analysis(pressure_data, TRUE, variable = "force_ts")
#' @importFrom sf st_intersects st_geometry st_area
#' @importFrom pracma trapz
#' @export

mask_analysis <- function(pressure_data, partial_sensors = FALSE,
                          variable = "press_peak_sensor",
                          pressure_units = "kPa", area_units = "cm2") {
  # set global variables
  sens_poly <- act_sens <- area <- max_df <- overlap_list <- NULL

  # masks
  masks <- pressure_data[[5]]

  # set up mask/sensor areas
  ## sensor area
  sensor_area <- pressure_data[[3]][1] * pressure_data[[3]][2]

  # active sensors
  act_sens <- which(footprint(pressure_data, "max") > 0)

  ## Make active sensors into polygons
  sens_coords <- sensor_coords(pressure_data, pressure_image = "all_active")
  sens_poly <- sensor_2_polygon(pressure_data)

  ## For each region mask, find which polygons intersect
  sens_mask_df <- matrix(rep(0, length.out = (length(sens_poly) * length(masks))),
                         nrow = length(sens_poly), ncol = length(masks))
  sens_mask_df <- as.data.frame(sens_mask_df)
  colnames(sens_mask_df) <- names(masks)
  for (i in 1:length(masks)){
    for (j in 1:length(sens_poly)) {
      x <- st_intersects(masks[[i]], sens_poly[[j]])
      if (identical(x[[1]], integer(0)) == FALSE) {
        y <- st_intersection(masks[[i]], sens_poly[[j]])
        sens_mask_df[j, i] <- st_area(y) / sensor_area
      }
    }
  }

  # calculate output variables. If variable name ends "_ts" output is vector
  # (per mask) equal to length of measurement, otherwise output is a single
  # value per mask
  ## storage matrices
  if (str_ends(variable, "_ts")) {
    output_mat <- matrix(rep(0, length.out = dim(pressure_data[[1]])[3] *
                               length(masks)),
                         nrow = dim(pressure_data[[1]])[3],
                         ncol = length(masks))
  } else {
    output_mat <- matrix(rep(NA, times = length(masks)), nrow = 1)
  }

  ## Analyse regions for maximum value of any sensor within region during trial
  if (variable == "press_peak_sensor") {
    P <- c(footprint(pressure_data))
    P <- P[act_sens]
    for (mask in seq_along(masks)) {
      output_mat[, mask] <- max(P[which(sens_mask_df[, mask] > 0)])
    }
    if (pressure_units == "MPa") {output_mat <- output_mat * 0.001}
    if (pressure_units == "Ncm2") {output_mat <- output_mat * 0.1}
  }

  # Analyse regions for maximum value of any sensor for each measurement frame
  if (variable == "press_peak_sensor_ts") {
    for (mask in seq_along(masks)) {
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- c(pressure_data[[1]][, , i])
        P <- P[act_sens]
        output_mat[i, mask] <- max(P[which(sens_mask_df[, mask] > 0)])
      }
    }
    if (pressure_units == "MPa") {output_mat <- output_mat * 0.001}
    if (pressure_units == "Ncm2") {output_mat <- output_mat * 0.1}
  }

  # Analyse regions for peak regional pressure (defined as the pressure
  # averaged over the mask)
  if (variable == "press_peak_mask") {
    P <- c(footprint(pressure_data))
    P <- P[act_sens] * sensor_area * 1000
    CA <- rep(sensor_area, length.out = length(P))
    for (mask in seq_along(masks)) {
      force <- sum(P * sens_mask_df[, mask])
      contact_area <- sum(CA * sens_mask_df[, mask])
      output_mat[, mask] <- force / contact_area / 1000
    }
    if (pressure_units == "MPa") {output_mat <- output_mat * 0.001}
    if (pressure_units == "Ncm2") {output_mat <- output_mat * 0.1}
  }

  # Analyse regions for peak mask pressure (defined as the maximum mean pressure
  # of active sensors in region during the trial). Outputs 101 point vector (kPa)
  if (variable == "press_peak_mask_ts") {
    for (mask in seq_along(masks)) {
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- c(pressure_data[[1]][, , i])
        P <- P[act_sens] * sensor_area * 1000
        CA <- rep(sensor_area, length.out = length(P))
        force <- sum(P * sens_mask_df[, mask])
        contact_area <- sum(CA * sens_mask_df[, mask])
        output_mat[i, mask] <- force / contact_area / 1000
      }
    }
    if (pressure_units == "MPa") {output_mat <- output_mat * 0.001}
    if (pressure_units == "Ncm2") {output_mat <- output_mat * 0.1}
  }

  # Analyse regions for maximum contact area of region during trial
  if (variable == "contact_area_peak") {
    CA <- rep(sensor_area, length.out = length(act_sens))
    for (mask in seq_along(masks)) {
      cArea <- sum(CA * sens_mask_df[, mask])
      if (area_units == "cm2") {cArea <- cArea * 10000}
      if (area_units == "mm2") {cArea <- cArea * 1000000}
      output_mat[, mask] <- cArea
    }
  }

  # Analyse regions for contact area throughout the trial (outputs vector)
  if (variable == "contact_area_ts") {
    for (mask in seq_along(masks)) {
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- c(pressure_data[[1]][, , i])
        P <- P[act_sens]
        P[P > 0] <- sensor_area
        cArea <- sum(P * sens_mask_df[, mask])
        if (area_units == "cm2") {cArea <- cArea * 10000}
        if (area_units == "mm2") {cArea <- cArea * 1000000}
        output_mat[i, mask] <- cArea
      }
    }
  }

  # Analyse regions for pressure time integral (Novel definition)
  if (variable == "pti_1") {
    for (mask in seq_along(masks)) {
      mask_pp <- rep(NA, length.out = dim(pressure_data[[1]])[3])
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- c(pressure_data[[1]][, , i])
        P <- P[act_sens]
        mask_pp[i] <- max(P[which(sens_mask_df[, mask] > 0)]) *
          pressure_data[[4]] * 1000
      }
      output_mat[, mask] <- sum(mask_pp)
    }
    if (area_units == "cm2") {output_mat <- output_mat * 1e-04}
    if (area_units == "mm2") {output_mat <- output_mat * 1e-06}
  }

  # Analyse regions for pressure time integral (Melai definition)
  if (variable == "pti_2") {
    time_seq <- seq(0, by = pressure_data[[4]],
                    length.out = dim(pressure_data[[1]])[3])
    for (mask in seq_along(masks)) {
      force <- rep(NA, length.out = dim(pressure_data[[1]])[3])
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- c(pressure_data[[1]][, , i])
        P <- P[act_sens] * sensor_area * 1000
        force[i] <- sum(P * sens_mask_df[, mask])
      }
      CA <- rep(sensor_area, length.out = length(act_sens))
      CA <- sum(CA * sens_mask_df[, mask])
      output_mat[, mask] <- pracma::trapz(time_seq, force) / CA / 10000
    }
    if (area_units == "cm2") {output_mat <- output_mat * 1e-04}
    if (area_units == "mm2") {output_mat <- output_mat * 1e-06}
  }

  # Analyse regions for force time integral
  if (variable == "force_time_integral") {
    time_seq <- seq(0, by = pressure_data[[4]],
                    length.out = dim(pressure_data[[1]])[3])
    for (mask in seq_along(masks)) {
      force <- rep(NA, length.out = dim(pressure_data[[1]])[3])
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- c(pressure_data[[1]][, , i])
        P <- P[act_sens] * sensor_area * 1000
        force[i] <- sum(P * sens_mask_df[, mask])
      }
      output_mat[, mask] <- pracma::trapz(time_seq, force)
    }
    if (area_units == "cm2") {output_mat <- output_mat * 1e-04}
    if (area_units == "mm2") {output_mat <- output_mat * 1e-06}
  }

  # Analyse regions for maximum force during the trial
  if (variable == "force_peak") {
    P <- c(footprint(pressure_data))
    P <- P[act_sens] * sensor_area * 1000
    for (mask in seq_along(masks)) {
      output_mat[, mask] <- sum(P * sens_mask_df[, mask])
    }
  }

  # Analyse regions for force throughout the trial (outputs vector)
  if (variable == "force_ts") {
    for (mask in seq_along(masks)) {
      for (i in 1:(dim(pressure_data[[1]])[3])) {
        P <- c(pressure_data[[1]][, , i])
        P <- P[act_sens] * sensor_area * 1000
        output_mat[i, mask] <- sum(P * sens_mask_df[, mask])
      }
    }
  }

  # return
  output_df <- as.data.frame(output_mat)
  colnames(output_df) <- names(masks)
  return(output_df)
}


# =============================================================================
# =============================================================================

# helper functions

#' interpolation function
#' @importFrom stats approx
#' @noRd
approxP <- function(x, interp_to) {
  y <- approx(x, n = interp_to)
  y$x <- NULL
  y <- unname(unlist(y))
  return(y)
}


#' @title Get minimum bounding box (angle)
#' @description Create a bounding box around sensor coordinate array. Adapted
#' from function in shotGroups package
#' @param xy Matrix. n row by 2 column with x-y coords of points to be
#' bounded
#' @return List
#' @importFrom grDevices chull
#' @noRd
getMinBBox <- function(xy) {
  stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) >= 2, ncol(xy) == 2)

  ## rotating calipers algorithm using the convex hull
  H    <- chull(xy)           # hull indices, vertices ordered clockwise
  n    <- length(H)           # number of hull vertices
  hull <- xy[H, ]             # hull vertices

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
  rangeH  <- matrix(numeric(n * 2), ncol = 2)   # hull edge
  rangeO  <- matrix(numeric(n * 2), ncol = 2)   # orth subspace
  widths  <- numeric(n)
  heights <- numeric(n)
  for(i in seq(along = H)) {
    rangeH[i, ] <- range(projMat[  i, ])
    rangeO[i, ] <- range(projMat[n + i, ])  # orth subspace is in 2nd half
    widths[i]   <- abs(diff(rangeH[i, ]))
    heights[i]  <- abs(diff(rangeO[i, ]))
  }

  ## extreme projections for min-area rect in subspace coordinates
  eMin  <- which.min(widths * heights)   # hull edge leading to minimum-area
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

  # Which sides of BBox are longest?
  if (dist(pts[c(1, 2), ]) > dist(pts[c(2, 3), ])) {
    side1 <- pts[c(1, 2), ]
    side2 <- pts[c(3, 4), ]
  } else {
    side1 <- pts[c(2, 3), ]
    side2 <- pts[c(4, 1), ]
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

  mat <- rbind(side1_, side2_)

  #return(list(pts=pts, width=widths[eMin], height=heights[eMin]))
  return(mat)
}


#' @title masks (sf format) to dataframe
#' @description converts mask coordinates into a data frame
#' @param masks List
#' @returns Data frame
#' @importFrom dplyr bind_rows
#' @noRd
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

#' @title Get coordinates of active sensors
#' @description Produces a data frame with coordinates of sensors
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. Includes a 3D array covering each timepoint of the
#'   measurement. z dimension represents time
#' @param pressure_image. String. Which pressure image to use. Options are
#' "all_active", "all", or "frame".
#' @param frame Numeric. If pressure image is frame, the numeric value should be
#' provided here
#' @return Data frame. x and y coordinates of sensors
#' @examples
#' coords <- sensor_coords(pressure_data)
#' @noRd
sensor_coords <- function(pressure_data, pressure_image = "all_active", frame) {
  # pressure image
  if (pressure_image == "all_active" | pressure_image == "max") {
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
  P <- as.vector(sens)
  if (pressure_image != "all") {
    coords <- coords[which(P > 0), ]
  }

  # return sensor coordinates
  return(coords)
}



#' @title sensor array to polygon
#' @description converts sensor positions (coordinates) to polygons
#' @param pressure_data List
#' @param pressure_image String
#' @param frame Numeric
#' @param output String. "df" dataframe for plotting or "sf" shape poly for analysis
#' @importFrom sf st_polygon st_coordinates
#' @noRd
sensor_2_polygon <- function(pressure_data, pressure_image = "all_active",
                             frame = NA, output = "sf") {
  # set global variables
  pedar_insole_grid <- NULL

  # sensor coordinates
  sens_polygons <- list()
  if (pressure_data[[2]] != "pedar") {
    sens_coords <- sensor_coords(pressure_data, pressure_image, frame)

    # sensor dimensions
    width <- pressure_data[[3]][1] / 2
    height <- pressure_data[[3]][2] / 2

    # get corner points and make into polygon
    for (sens in 1:nrow(sens_coords)) {
      x1 <- sens_coords[sens, 1] - width
      y1 <- sens_coords[sens, 2] + height
      x2 <- sens_coords[sens, 1] + width
      y2 <- sens_coords[sens, 2] + height
      x3 <- sens_coords[sens, 1] + width
      y3 <- sens_coords[sens, 2] - height
      x4 <- sens_coords[sens, 1] - width
      y4 <- sens_coords[sens, 2] - height
      if (y1 < 0) {y1 <- 0}
      if (y2 < 0) {y2 <- 0}
      if (y3 < 0) {y3 <- 0}
      if (y4 < 0) {y4 <- 0}
      sens_polygons[[sens]] <- st_polygon(list(matrix(c(x1, x2, x3, x4, x1,
                                                        y1, y2, y3, y4, y1),
                                                      5, 2)))
    }
  } else {
    load("data/pedar_insole_grid.rda")
    rs <- c(1:99, 101:199)
    for (sens in 1:length(rs)) {
      x1 <- pedar_insole_grid[rs[sens], 1]
      y1 <- pedar_insole_grid[rs[sens], 2]
      x2 <- pedar_insole_grid[rs[sens], 3]
      y2 <- pedar_insole_grid[rs[sens], 4]
      x3 <- pedar_insole_grid[rs[sens], 5]
      y3 <- pedar_insole_grid[rs[sens], 6]
      x4 <- pedar_insole_grid[rs[sens], 7]
      y4 <- pedar_insole_grid[rs[sens], 8]
      sens_polygons[[sens]] <- st_polygon(list(matrix(c(x1, x2, x3, x4, x1,
                                                        y1, y2, y3, y4, y1),
                                                      5, 2)))
    }
  }

  # make df if required
  if (output == "df") {
    id_df <- data.frame(x = double(),
                        y = double(),
                        z = integer())
    # make into identity df
    for (i in 1:length(sens_polygons)) {
      mat <- st_coordinates(sens_polygons[[i]])[, c(1, 2)]
      mat_df <- data.frame(mat)
      id <- rep(i, length.out = nrow(mat_df))
      mat_df <- cbind(mat_df, id)
      colnames(mat_df) <- c("x", "y", "id")
      id_df <- rbind(id_df, mat_df)
    }
  }

  # check for minus values

  # return
  if (output == "sf") {return(sens_polygons)}
  if (output == "df") {return(id_df)}
}

#' @title pedar force
#' @description generates force curve from pedar data
#' @param pressure_data List.
#' @param variable "force", "area"
#' @param threshold Numeric. Threshold value for sensor to be considered active.
#' Currently only applies to insole data
#' @param threshold Numeric. For area calculations, minimum threshold for sensor to
#' be active
#' @return Vector
#' @noRd
force_pedar <- function(pressure_data, variable = "force", threshold = 10) {
  # set global variables
  pedar_insole_areas <- NULL

  # check this is pedar data
  if (pressure_data[[2]] != "pedar")
    stop("must be pedar data")

  # get insole sizing
  pedar_insole_type <- pressure_data[[3]]

  # get pedar sensor areas
  load("data/pedar_insole_areas.rda")
  pedarSensorAreas <- as.vector(pedar_insole_areas[[pedar_insole_type]] *
                                  0.001)

  # change array direction
  force_array <- aperm(pressure_data[[1]], c(3, 2, 1))

  # output
  output_df = data.frame()

  # calculate force
  if (variable == "force") {
    force_right <- rowSums(force_array[, , 1] * pedarSensorAreas)
    force_left <- rowSums(force_array[, , 2] * pedarSensorAreas)
    output_df <- cbind(force_right, force_left)
  }

  # calculate area
  if (variable == "area") {
    area_right <- force_array[, , 1] > threshold
    area_right <- rowSums(area_right * pedarSensorAreas)
    area_left <- force_array[, , 2] > threshold
    area_left <- rowSums(area_left * pedarSensorAreas)
    output_df <- cbind(area_right, area_left)
  }

  # return
  return(output_df)
}

#' @title plot pedar
#' @description produces dataframe that plot_pressure uses to plot pedar data
#' @param pressure_data List.
#' @param pressure_image String. "max" = maximum values for full trial;
#' "mean" = mean value across trial; "step_max"
#' @param foot_side String. "both", "left", or "right"
#' @param step_n Numeric. If pressure_image is "step_max", the step number to
#' analyze
#' @return df
#' @noRd
plot_pedar <- function(pressure_data, pressure_image = "max",
                       foot_side = "left", step_n) {
  # set global variables
  pedar_insole_grid <- x <- y <- id <- NULL

  # check this is pedar (or other suitable) data
  if (!(pressure_data[[2]] == "pedar"))
    stop("data should be from pedar")

  # load pedar coords
  load("data/pedar_insole_grid.rda")

  # get L+ R data frames
  dims <- dim(pressure_data[[1]])
  pressure_R_mat <- matrix(pressure_data[[1]][1, , ],
                           dims[2], dims[3])
  pressure_L_mat <- matrix(as.vector(pressure_data[[1]][2, , ]),
                           dims[2], dims[3])

  # pressure image
  if (pressure_image == "max") {
    R_data <- apply(pressure_R_mat, 1, max)
    L_data <- apply(pressure_L_mat, 1, max)
  }

  # separate into steps
  if (pressure_image == "step_max") {
    start_end_R <- pressure_data[[6]]
    start_end_L <- pressure_data[[6]]
    #pressure_R_mat <- pressure_R_mat[start_end_R[],]
  }

  # combine data
  comb_data <- c(L_data, R_data)

  # make ids
  ids <- c()
  for (i in 1:99) {ids = append(ids, paste0("L", i))}
  for (i in 1:99) {ids = append(ids, paste0("R", i))}

  # make df
  df <- data.frame(id = ids, value = comb_data)

  # add coordinates to df
  xs <- c()
  for (i in c(101:199, 1:99)) {
    for (j in c(1, 3, 5, 7)) {xs = append(xs, pedar_insole_grid[i, j])}
  }
  ys <- c()
  for (i in c(101:199, 1:99)) {
    for (j in c(2, 4, 6, 8)) {ys = append(ys, pedar_insole_grid[i, j])}
  }
  position <- data.frame(id = rep(ids, each = 4), x = xs, y = ys)
  df <- merge(df, position, by = c("id"))

  # add color
  df <- generate_colors(df)

  # side
  if (foot_side == "left") {
    df <- df[397:792,]
  }
  if (foot_side == "right") {
    df <- df[1:396,]
  }

  # return
  return(df)
}


#' @title Generate colors
#' @description Let's the user prescribe the color scale for plots
#' @param df Dataframe.
#' @param col_type String. "default": novel color scheme; "custom": user
#' supplied
#' @param break_values Vector. Vector with break values to be used. Should be one
#' shorter than break_values
#' @param break_colors Vector. Vector with colors to be used. Should be one
#' longer than break_values
#' @return Data frame.
#' @noRd
generate_colors <- function(df, col_type = "default", break_values,
                            break_colors) {
  if (col_type == "default") {
    break_values <- c(0, 40, 60, 100, 150, 220, 300, 1000000)
    break_colors <- c("grey","lightblue", "darkblue","green","yellow",
                      "red","pink")
  } else {
    # check break_values and break_colors
    if (break_values[1] <= 0)
      stop("first break value should be >0")
    if(length(break_values) != (length(break_colors) - 1))
      stop("break_values should be one shorter than break_colors")
    if ("white" %in% break_colors)
      stop("break_colors cannot contain white")

    # add high final value to break_values
    break_values <- c(-1, 0, break_values, 1000000)
    break_colors <- c("white", break_colors)
  }

  # add col column
  df$cols <- cut(df$value, breaks = break_values,
                 labels = break_colors)

  #return
  return (df)
}

#' @title extend st line
#' @param line sf linestring
#' @param distance Numeric.
#' @param end String
#' @importFrom sf st_sfc st_crs st_coordinates
#' @noRd
st_extend_line <- function(line, distance, end = "BOTH") {
  if (!(end %in% c("BOTH", "HEAD", "TAIL")) | length(end) != 1)
    stop("'end' must be 'BOTH', 'HEAD' or 'TAIL'")

  M <- st_coordinates(line)[,1:2]
  keep <- !(end == c("TAIL", "HEAD"))

  ends <- c(1, nrow(M))[keep]
  i <- c(2, nrow(M) - 1)
  j <- c(1, -1)
  headings <- mapply(i, j, FUN = function(i, j) {
    Ax <- M[i-j,1]
    Ay <- M[i-j,2]
    Bx <- M[i,1]
    By <- M[i,2]
    unname(atan2(Ay-By, Ax-Bx))
  })
  #headings <- st_ends_heading(line)[keep]
  distances <- if (length(distance) == 1) rep(distance, 2) else rev(distance[1:2])

  M[ends,] <- M[ends,] + distances[keep] * c(cos(headings), sin(headings))
  newline <- st_linestring(M)

  # If input is sfc_LINESTRING and not sfg_LINESTRING
  if (is.list(line)) newline <- st_sfc(newline, crs = st_crs(line))

  return(newline)
}

#' @title line to polygon
#' @description Extrude line to form polygon
#' @param mat Matrix. xy points describing line or polyline
#' @param distance Numeric. Distance to extrude line
#' @param direction String. "+X" "-X", "+y", or "-Y"
#' @return Polygon. sf polygon object
#' @importFrom sf st_polygon
#' @noRd
st_line2polygon <- function(mat, distance, direction) {
  # get ends
  if (mat[1, 1] >= mat[nrow(mat), 1]) {dir1 <- "RL"} else {dir1 <- "LR"}
  if (mat[1, 2] >= mat[nrow(mat), 2]) {dir2 <- "TB"} else {dir2 <- "BT"}

  # +Y
  if (direction == "+Y" ) {
    mat_pts <- rbind(mat, c(mat[nrow(mat), 1], mat[nrow(mat), 2] + distance),
                     c(mat[1, 1], mat[1, 2] + distance), mat[1, ])
    if (dir1 == "RL") {mat <- mat[nrow(mat):1, ]}
  }
  # -Y
  if (direction == "-Y") {
    mat_pts <- rbind(mat, c(mat[nrow(mat), 1], mat[nrow(mat), 2] - distance),
                     c(mat[1, 1], mat[1, 2] - distance), mat[1, ])
    if (dir1 == "LR") {mat <- mat[nrow(mat):1, ]}
  }
  # +X
  if (direction == "+X") {
    mat_pts <- rbind(mat, c(mat[nrow(mat), 1] + distance, mat[nrow(mat), 2]),
                     c(mat[1, 1] + distance, mat[1, 2]), mat[1, ])
    if (dir2 == "TB") {mat <- mat[nrow(mat):1, ]}
  }

  # -X
  if (direction == "-X") {
    mat_pts <- rbind(mat, c(mat[nrow(mat), 1] - distance, mat[nrow(mat), 2]),
                     c(mat[1, 1] - distance, mat[1, 2]), mat[1, ])
    if (dir2 == "BT") {mat <- mat[nrow(mat):1, ]}
  }

  # make polygon
  mat_poly <- st_polygon(list(mat_pts))
}

#' @title Get toe line
#' @description Calculates line proximal to toes
#' @param pressure_data List. First item should be a 3D array covering each
#' timepoint of the measurement. z dimension represents time
#' @importFrom zoo as.zoo rollapply
#' @importFrom stats na.omit
#' @importFrom sf st_cast st_join st_sf
#' @noRd
toe_line <- function(pressure_data) {
  # global variables
  clusterID <- NULL

  # get max footprint and take top half
  pf_max <- footprint(pressure_data)
  pf_max_top <- pf_max[1:(round(nrow(pf_max)) / 2), ]
  pressure_data2 <- pressure_data
  pressure_data2[[1]] <- pf_max_top
  act_sens_vec <- which(pf_max_top > 0)

  # remove small islands (toes)
  ## make polygon df
  polygons <- sensor_2_polygon(pressure_data2, output = "sf")
  pg_df <- data.frame(sens_id = 1:length(polygons))
  pg_df$geometry <- st_sfc(polygons)
  pg_df <- st_as_sf(pg_df)

  ## find islands
  geometries <- st_cast(st_union(st_buffer(pg_df, 0.0001)), "POLYGON")
  dissolved <- st_sf(geometries)
  dissolved$clusterID = 1:length(geometries)
  pg_df <- st_join(pg_df, dissolved)

  ## remove small islands
  pg_df <- pg_df %>% filter(clusterID == which.max(table(pg_df$clusterID)))
  sens_keep <- pg_df[["sens_id"]]
  sens_keep <- act_sens_vec[sens_keep]
  pf_max_top_vec <- c(pf_max_top)
  zeros <- c(1:length(pf_max_top_vec))
  zeros <- zeros[!zeros %in% sens_keep]
  pf_max_top_vec[zeros] <- 0
  pf_max_top2 <- matrix(pf_max_top_vec, nrow = nrow(pf_max_top), ncol = ncol(pf_max_top))

  # which cols had toes removed
  good_cols <- rep(NA, length.out = ncol(pf_max))
  for (i in 1:ncol(pf_max)) {
    og <- which.max(pf_max_top[, i] > 0)
    og_rm <- which.max(pf_max_top2[, i] > 0)
    if (og != og_rm) {good_cols[i] <- og_rm - 2}
  }

  # remove negative values
  good_cols[good_cols < 1] <- NA

  # which of the remaining cols have minima between two maxima
  remaining_cols <- which(is.na(good_cols))
  for (i in seq_along(remaining_cols)) {
    mins <- which(rollapply(as.zoo(pf_max_top[, remaining_cols[i]]), 3,
                            function(x) which.min(x) == 2) == TRUE)
    mins_rev <- which(rollapply(as.zoo(rev(pf_max_top[, remaining_cols[i]])), 3,
                                function(x) which.min(x) == 2) == TRUE)
    maxs <- which(rollapply(as.zoo(pf_max_top[, remaining_cols[i]]), 3,
                            function(x) which.max(x) == 2) == TRUE)
    if (length(maxs) >= 2) {
      if (mins[1] > maxs[1] & mins[1] < maxs[2]) {
        reps <- rle(pf_max_top[(mins[1] + 1):ncol(pf_max_top), remaining_cols[i]])$lengths[1]
        if (reps > 1) {
          good_cols[i] <- mins[1] + (reps / 2)
        } else {good_cols[remaining_cols[i]] <- mins[1]}
      }
    }
  }

  # check points are less than 1/3 from front of foot but > 15%
  ff_dist <- round(nrow(pf_max) / 6)
  ff_prox <- round(nrow(pf_max) / 12)
  remaining_cols2 <- which(!is.na(good_cols))
  for (fp_col in seq_along(remaining_cols2)) {
    if ((good_cols[remaining_cols2[fp_col]] - (which(pf_max_top[, remaining_cols2[fp_col]] > 0)[1])) >= ff_dist) {
      good_cols[remaining_cols2[fp_col]] <- NA
    }
  }
  remaining_cols3 <- which(!is.na(good_cols))
  for (fp_col in seq_along(remaining_cols3)) {
    if ((good_cols[remaining_cols3[fp_col]] - (which(pf_max_top[, remaining_cols3[fp_col]] > 0)[1])) < ff_prox) {
      good_cols[remaining_cols3[fp_col]] <- NA
    }
  }

  # coords
  toe_line_mat <- matrix(NA, length(good_cols), 2)
  dims <- pressure_data[[3]]
  for (i in 1:length(good_cols)) {
    if (is.na(good_cols[i]) != TRUE) {
      xy <- c(-(dims[2] / 2) + (dims[2] * i),
              (-(dims[1] / 2) + (nrow(pf_max) * dims[1]) - dims[1] * good_cols[i]))
      toe_line_mat[i, ] <- xy
    }
  }
  toe_line_mat <- na.omit(toe_line_mat)

  # extend ends
  toe_line_mat <- rbind(c(toe_line_mat[1, 1] - 1, toe_line_mat[1, 2]), toe_line_mat,
                        c(toe_line_mat[nrow(toe_line_mat), 1] + 1,
                          toe_line_mat[nrow(toe_line_mat), 2]))

  #toe_path_df <- as.data.frame(cbind(toe_line_mat[, 1], toe_line_mat[, 2]))
  #colnames(toe_path_df) <- c("x", "y")
  #g <- plot_pressure(pressure_data)
  #g <- g + geom_path(data = toe_path_df, aes(x = x, y = y), color = "red")
  #g

  # return
  return(toe_line_mat)
}

#' @title edge lines
#' @description create edge liens for footprint data
#' @param mat Matrix
#' @param side String
#' @importFrom sf st_as_sfc st_combine st_coordinates st_convex_hull
#' @importFrom utils head tail
#' @noRd
edge_lines <- function(pressure_data, side) {
  # global variables
  x <- y <- x_coord <- y_coord <- me <- sc_df <- NULL

  # max pressure image
  max_fp <- footprint(pressure_data, "max")

  # coordinates
  sens_coords <- sensor_coords(pressure_data)

  # Find longest vectors (these are the med and lat edges of the footprint)
  ## unique y coordinates of sensors
  unq_y <- unique(sens_coords$y)

  ## how much to shorten by
  shorten_n <- ceiling(ncol(max_fp) / 10)
  unq_y_short <- unq_y[-match(tail(sort(unq_y), shorten_n), unq_y)]
  unq_y_short <- unq_y_short[-match(head(sort(unq_y_short), shorten_n),
                                    unq_y_short)]

  ## make edges
  med_edge <- data.frame(x = rep(NA, length.out = length(unq_y_short)),
                         y = unq_y_short)
  lat_edge <- med_edge
  for (i in 1:length(unq_y_short)) {
    med_edge[i, 1] <- sens_coords %>% filter(y_coord == unq_y_short[i]) %>%
      summarise(me = max(x_coord)) %>% pull(me)
    lat_edge[i, 1] <- sens_coords %>% filter(y_coord == unq_y_short[i]) %>%
      summarise(me = min(x_coord)) %>% pull(me)
  }
  if (side == "RIGHT") {
    x <- med_edge
    med_edge <- lat_edge
    lat_edge <- x
  }

  ### convex hulls of edges
  med_edge_sf <- med_edge %>% st_as_sf(coords = c("x", "y"))
  med_edge_chull <- st_convex_hull(st_combine(med_edge_sf))
  lat_edge_sf <- lat_edge %>% st_as_sf(coords = c("x", "y"))
  lat_edge_chull <- st_convex_hull(st_combine(lat_edge_sf))

  ## longest med and lateral edge line
  me_dis <- as.matrix(dist(st_coordinates(med_edge_chull)))
  me_dis <- me_dis[row(me_dis) == (col(me_dis) - 1)]
  me_max <- order(me_dis)[length(me_dis)]
  med_side <- st_coordinates(med_edge_chull)[c(me_max, me_max + 1), c(1, 2)]
  med_side_line <- st_linestring(med_side)
  med_side_line <- st_extend_line(med_side_line, 1)
  le_dis <- as.matrix(dist(st_coordinates(lat_edge_chull)))
  le_dis <- le_dis[row(le_dis) == (col(le_dis) - 1)]
  le_max <- order(le_dis)[length(le_dis) - 1]
  lat_side <- st_coordinates(lat_edge_chull)[c(le_max, le_max + 1), c(1, 2)]
  lat_side_line <- st_linestring(lat_side)
  lat_side_line <- st_extend_line(lat_side_line, 1)

  # check that no points fall outside of line
  if (side == "RIGHT") {
    med_poly <- st_line2polygon(med_side_line, 1, "-X")  + c(-0.001, 0)
    lat_poly <- st_line2polygon(lat_side_line, 1, "+X")  + c(0.001, 0)
  }
  #med_int <- st_intersects(st_as_sfc(med_edge_sf), med_poly)
  #lat_int <- st_intersects(st_as_sfc(lat_edge_sf), lat_poly)

  # if points do fall in poly, use fitted line approach

  # return lines
  return(list(med_side_line, lat_side_line))
}

#' Color binning function
#' taken from ggplot2 codebase
#' @noRd
binned_pal <- function(palette) {
  function(x) {
    palette(length(x))
  }
}


rot_line <- function(line, ang, cnt) {
  ang_ <- ang * pi /180
  new_line <- (line - cnt) * matrix(c(cos(ang_), sin(ang_),
                                      -sin(ang_), cos(ang_)), 2, 2) + cnt
  return(new_line)
}

#' @title Visual masks
#' @description Visualize the existing masks
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_data List. First item is a 3D array covering each timepoint
#' of the measurement.
#' @param visual_list List. Mask numbers that want to be viewed. (Default is
#' all exisiting masks)
#' @return ggplot plot object
#' @examples
#' \dontrun{
#' emed_data <- system.file("extdata", "emed_test.lst", package = "pressuRe")
#' pressure_data <- load_emed(emed_data)
#' pressure_data <- automask(pressure_data, foot_side = "auto", plot = TRUE)
#' plot_masks(pressure_data)
#' }
#' @noRd

plot_masks <- function(pressure_data, visual_list = seq(1, length(pressure_data[[5]]))){
  # global variables
  X <- Y <- NULL

  # plot footprint
  grDevices::x11()
  g <- plot_pressure(pressure_data, step_n = 2, plot = FALSE)
  if (pressure_data[[2]] != "pedar"){
    g <- g + scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 0.15))
    g <- g + scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 0.30))
  }
  print(g)

  # plot original mask data
  for (n_mask in visual_list) {
    g <-
      g + geom_path(data = data.frame(st_coordinates(pressure_data[[5]][[n_mask]])[,1:2]), aes(X, Y), color = "red", linewidth = 1)
  }
  print(g)

  return(g)
}

#' @title Pedar mask regions
#' @description Creates a mask region from pedar sensel coordinates
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param position Dataframe. A n x 3 dataframe of sensel coordinates
#' [sensel id, x,y]
#' @param sensel_list List. List of sensels to include
#' @param foot_side String. "LEFT" for left foot, "RIGHT" for right foot
#' @return polygon. A mask polygon
#' @noRd

pedar_polygon <- function(position, sensel_list, foot_side){
  # global variables
  y <- id <- NULL

  # mask value based on insole side
  if(foot_side == "LEFT"){
    coord <- subset(position[1:(99 * 4), ], id %in% sensel_list)
  } else if (foot_side == "RIGHT"){
    coord <- subset(position[(99 * 4 + 1):nrow(position), ], id %in% sensel_list)
  }

  # mask from convex hull of sensel vertices
  sens_coords<-data.frame(x_coord = double(), y_coord = double())
  sens_coords[1:(length(sensel_list)*4),] <- coord[,2:3]
  df_sf <- sens_coords %>%
    st_as_sf(coords = c( "x_coord", "y_coord" ))
  fp_chull <- st_convex_hull(st_union(df_sf))

  # sensel 48 is slanted, so this will account for it in the mask
  if(length(which(sensel_list %in% 48)) > 0 ){
    if(length(which(sensel_list %in% c(27,61))) == 2){

     tri_coords <- data.frame(x_coord= double(), y_coord = double())
     if(foot_side == "L"){
       edge_sensels <- c(6,13,20,27,34,41,48,55,62,69)
       tri_coord <- subset(position[1:(99*4),], id %in% edge_sensels)
       tri_coords <- tri_coord[(seq(1,nrow(tri_coord/4), by = 4)+1),2:3]
       tri_coords[nrow(tri_coords)+1,] <- tri_coord[2,2:3]
     }else if(foot_side == "R"){
       edge_sensels <- c(13,20,27,34,41,48,55,62,69,76)
       tri_coord <- subset(position[(99*4 + 1):nrow(position),], id %in%  edge_sensels)
       tri_coords <- tri_coord[(seq(1,nrow(tri_coord/4), by = 4)+2),2:3]
       tri_coords[nrow(tri_coords)+1,] <- tri_coord[3,2:3]
      }
    } else {
      mask_y_low <- min(sens_coords$y_coord)
      mask_y_high <- max(sens_coords$y_coord)
      if (foot_side == "LEFT"){
        tri_coord <- subset(position[1:(99*4),], id %in% 48)
        mask_x_low <- max(subset(coord, y %in% mask_y_low)[,2])
        mask_x_high <- max(subset(coord, y %in% mask_y_high)[,2])
      } else if (foot_side == "RIGHT"){
        coord <- subset(position[(99*4 + 1):nrow(position),], id %in% sensel_list)
        tri_coord <- subset(position[(99*4 + 1):nrow(position),], id %in%  48)
        mask_x_low <- min(subset(coord, y %in% mask_y_low)[,2])
        mask_x_high <- min(subset(coord, y %in% mask_y_high)[,2])
      }
      tri_coords <- data.frame(x_coord = double(), y_coord = double())
      tri_coords[1,] <- c(mask_x_high, mask_y_high)
      tri_coords[2,] <- tri_coord[3, 2:3]
      tri_coords[3,] <- c(mask_x_low, mask_y_low)
      tri_coords[4,] <- c(mask_x_high, mask_y_high)
    }

      tri_sf <-  st_polygon(list(as.matrix(tri_coords)))
      fp_chull <- st_difference(fp_chull, tri_sf)
  }

  # return
  return(fp_chull)
}
