#' Perform a fit based analysis of MRS data.
#' 
#' Note that TARQUIN and LCModel require these packages to be installed, and 
#' the functions set_tqn_cmd and set_lcm_cmd (respectively) need to be used to 
#' specify the location of these software packages.
#'
#' Fitting approaches described in the following references:
#' VARPRO
#' van der Veen JW, de Beer R, Luyten PR, van Ormondt D. Accurate quantification
#' of in vivo 31P NMR signals using the variable projection method and prior 
#' knowledge. Magn Reson Med 1988;6:92-98.
#' 
#' TARQUIN
#' Wilson, M., Reynolds, G., Kauppinen, R. A., Arvanitis, T. N. & Peet, A. C. 
#' A constrained least-squares approach to the automated quantitation of in vivo 
#' 1H magnetic resonance spectroscopy data. Magn Reson Med 2011;65:1-12.
#' 
#' LCModel
#' Provencher SW. Estimation of metabolite concentrations from localized in vivo
#' proton NMR spectra. Magn Reson Med 1993;30:672-679.
#' 
#' @param metab metabolite data.
#' @param basis basis class object or character vector to basis file in 
#' LCModel .basis format.
#' @param method 'ABFIT' (default), 'VARPRO', 'VARPRO_3P', 'TARQUIN' or 
#' 'LCMODEL'.
#' @param w_ref water reference data for concentration scaling (optional).
#' @param opts options to pass to the analysis method.
#' @param parallel perform analyses in parallel (TRUE or FALSE).
#' @param time measure the time taken for the analysis to complete
#' (TRUE or FALSE).
#' @param progress option is passed to plyr::alply function to display a
#' progress bar during fitting. Default value is "text", set to "none" to
#' disable.
#' @return MRS analysis object.
#' @examples
#' fname <- system.file("extdata","philips_spar_sdat_WS.SDAT",package="spant")
#' svs <- read_mrs(fname, format="spar_sdat")
#' \dontrun{
#' basis <- sim_basis_1h_brain_press(svs)
#' fit_result <- fit_mrs(svs, basis)
#' }
#' @export
fit_mrs <- function(metab, basis = NULL, method = 'ABFIT', w_ref = NULL,
                    opts = NULL, parallel = FALSE, time = TRUE,
                    progress = "text") {
  
  # start the clock
  if (time) ptm <- proc.time()
  
  # force uppercase for matching
  METHOD <- toupper(method)
  
  # let's assume a PRESS basis and simulate if one isn't specified
  if (is.null(basis)) {
    if (METHOD == "LCMODEL") {
      lcm_compat = TRUE
    } else {
      lcm_compat = FALSE
    }
    TE1 = 0.01
    TE2 = metab$meta$EchoTime - TE1
    warning("Basis set not specified, so simulating default PRESS brain basis.")
    basis <- sim_basis_1h_brain_press(metab, lcm_compat = lcm_compat, TE1 = TE1, 
                                      TE2 = TE2) 
  }
  
  if (class(basis) == "basis_set") {
    # TODO check basis matches mrs_data 
  } else if (is.character(basis)) {
    if (!file.exists(basis)) {
      stop("Error, basis file not found.")
    }
  } else {
    stop("Error, specified basis is not the correct data type.")
  }
  
  # put data into TD
  if (is_fd(metab)) {
    metab <- fd2td(metab)
  }
  if (!is.null(w_ref)) {
    if (is_fd(w_ref)) {
      w_ref <- fd2td(w_ref)
    }
  }
  
  if (METHOD == "ABFIT") {
    # read basis into memory if a file
    if (is.character(basis)) {
      basis <- read_basis(basis)
    }
    
    # use default fitting opts if not specified 
    if (is.null(opts)) opts <- abfit_opts()
    
    acq_paras <- get_acq_paras(metab)
    
    result_list <- plyr::alply(metab$data, c(2, 3, 4, 5, 6), abfit, 
                               acq_paras, basis, opts, 
                               .parallel = parallel, 
                               .paropts = list(.inorder = TRUE,
                                               .packages = "spant"),
                               .progress = progress, .inform = FALSE)
    
  } else if (METHOD == "VARPRO") {
    # read basis into memory if a file
    if (is.character(basis)) {
      basis <- read_basis(basis)
    }
    
    acq_paras <- get_acq_paras(metab)
    
    result_list <- plyr::alply(metab$data, c(2, 3, 4, 5, 6), varpro,
                               acq_paras, basis, opts, 
                               .parallel = parallel, 
                               .paropts = list(.inorder = TRUE,
                                               .packages = "spant"),
                               .progress = progress, .inform = FALSE)
    
  } else if (METHOD == "VARPRO_3P") {
    # read basis into memory if a file
    if (is.character(basis)) {
      basis <- read_basis(basis)
    }
    
    acq_paras <- get_acq_paras(metab)
    
    result_list <- plyr::alply(metab$data, c(2, 3, 4, 5, 6),
                               varpro_3_para, acq_paras, basis, opts, 
                               .parallel = parallel, 
                               .paropts = list(.inorder = TRUE,
                                               .packages = "spant"),
                               .progress = progress, .inform = FALSE)
    
  } else if (METHOD == "TARQUIN") {
    if (!is.null(w_ref)) { 
      if (Ndyns(w_ref) > 1) {
        w_ref <- mean_dyns(w_ref)
        warning("Using the mean reference signal for water scaling.")
      }
      # repeat the refernce signal to match the number of dynamics
      if (Ndyns(metab) > 1) {
        w_ref <- rep_dyn(w_ref, Ndyns(metab))
      }
    }
    
    # combine metab and w_ref data together
    if (!is.null(w_ref)) {
      metab <- comb_metab_ref(metab, w_ref)
    }
    
    # write basis object (if specified) to file
    basis_file <- tempfile(fileext = ".basis")
    write_basis(basis, basis_file)
    
    temp_mrs <- metab
    temp_mrs$data = temp_mrs$data[1, 1, 1, 1, 1, 1,]
    dim(temp_mrs$data) <- c(1, 1, 1, 1, 1, 1, length(temp_mrs$data))
  
    #result_list <- apply(metab$data, c(2,3,4,5,6), tarquin_fit, temp_mrs, 
    #                     basis_file, opts)
    
    result_list <- plyr::alply(metab$data, c(2, 3, 4, 5, 6), tarquin_fit, 
                               temp_mrs, basis_file, opts, 
                               .parallel = parallel, 
                               .paropts = list(.inorder = TRUE,
                                               .packages = "spant"),
                               .progress = progress, .inform = FALSE)
    
                         #.paropts = list(.options.snow=snowopts),
                         #.paropts = list(.export="N",.packages="spant"),
  } else if (METHOD == "LCMODEL") {
    if (!is.null(w_ref)) { 
      if (Ndyns(w_ref) > 1) {
        w_ref <- mean_dyns(w_ref)
        warning("Using the mean reference signal for water scaling.")
      }
      
      # repeat the reference signal to match the number of dynamics
      if (Ndyns(metab) > 1) {
        w_ref <- rep_dyn(w_ref, Ndyns(metab))
      }
    }
    
    # combine metab and w_ref data together
    if (!is.null(w_ref)) {
      metab <- comb_metab_ref(metab, w_ref)
    }
  
    if (is.character(basis)) {
      basis_file <- basis  
    } else {
      # write basis object (if specified) to file
      basis_file <- tempfile(fileext = ".basis")
      write_basis(basis, basis_file)
    }
    
    temp_mrs <- metab
    temp_mrs$data = temp_mrs$data[1, 1, 1, 1, 1, 1,]
    dim(temp_mrs$data) <- c(1, 1, 1, 1, 1, 1, length(temp_mrs$data))
    
    result_list <- plyr::alply(metab$data, c(2,3,4,5,6), lcmodel_fit, 
                               temp_mrs, basis_file, opts,
                               .parallel = parallel, 
                               .paropts = list(.inorder = TRUE,
                                               .packages = "spant"),
                               .progress = progress, .inform = FALSE)
  } else if (exists(method)) {
    message(paste("Using external fit method :", method))
    
    # read basis into memory if a file
    if (is.character(basis)) {
      basis <- read_basis(basis)
    }
    
    acq_paras <- get_acq_paras(metab)
    
    result_list <- plyr::alply(metab$data, c(2, 3, 4, 5, 6),
                               get(method), acq_paras, basis, opts, 
                               .parallel = parallel, 
                               .paropts = list(.inorder = TRUE,
                                               .packages = "spant"),
                               .progress = progress, .inform = FALSE)
    
  } else {
    stop(paste('Fit method not found : ', method))
  }
  
  # combine fit results into dataframe and fits into a list
  fit_num <- length(result_list)
  labs <- attr(result_list, "split_labels")
  colnames(labs) <- c("X", "Y", "Z", "Dynamic", "Coil")
  
  # first 4 elements are amps, crlbs, diags, fit_table,
  # extra elements are temp files
  res_n <- length(result_list[[1]])
  
  result_list <- unlist(result_list, recursive = FALSE)
  
  if (res_n > 4) { # looks like temp files were used, so check for any colisions
    file_list_vec = vector()
    for (n in (5:res_n)) {
      file_list_vec <- c(file_list_vec, 
                         result_list[seq(from = n, by = res_n, 
                                         length.out = fit_num)])
    }
    # check that all temp i/o files are unique
    file_list_vec <- stats::na.omit(file_list_vec)
    file_list_vec <- as.character(file_list_vec)
    file_list_vec <- grep("NA", file_list_vec, value = TRUE, invert = TRUE)
    if (sum(duplicated(file_list_vec)) > 0 ) {
      stop("Error, duplicate temp file detected.")
    }
  }
  
  fits <- result_list[seq(from = 4, by = res_n, length.out = fit_num)]
 
  # check for any masked voxels 
  if (anyNA(fits)) {
    na_res <- TRUE
    na_res_vec <- is.na(fits)
    na_res_vec_ind <- which(is.na(fits))
    good_row <- which(!na_res_vec)[[1]]
    rowN <- length(fits)
  } else {
    na_res <- FALSE
  }
  
  df_list_amps  <- result_list[seq(from = 1, by = res_n, length.out = fit_num)]
  df_list_crlbs <- result_list[seq(from = 2, by = res_n, length.out = fit_num)]
  df_list_diags <- result_list[seq(from = 3, by = res_n, length.out = fit_num)]
 
  if (na_res) {
    amp_colnames <- colnames(df_list_amps[[good_row]])
    amps <- as.data.frame(matrix(NA, rowN, length(amp_colnames)))
    colnames(amps) <- amp_colnames
    amps_not_na <- plyr::ldply(df_list_amps[!na_res_vec], data.frame)[-1]
    amps[(!na_res_vec),] <- amps_not_na
    
    crlbs <- as.data.frame(matrix(NA, rowN, length(amp_colnames)))
    crlbs_not_na <- plyr::ldply(df_list_crlbs[!na_res_vec], data.frame)[-1]
    crlbs[(!na_res_vec),] <- crlbs_not_na
    colnames(crlbs) <- paste(colnames(amps), ".sd", sep = "")
    
    diags_colnames <- colnames(df_list_diags[[good_row]])
    diags <- as.data.frame(matrix(NA, rowN, length(diags_colnames)))
    colnames(diags) <- diags_colnames
    diags_not_na <- plyr::ldply(df_list_diags[!na_res_vec], data.frame)[-1]
    diags[(!na_res_vec),] <- diags_not_na
  } else {
    amps <- plyr::ldply(df_list_amps, data.frame)[-1]
  
    crlbs <- plyr::ldply(df_list_crlbs, data.frame)[-1]
    colnames(crlbs) <- paste(colnames(amps), ".sd", sep = "")
  
    diags <- plyr::ldply(df_list_diags, data.frame)[-1]
  }
  
  res_tab <- cbind(labs, amps, crlbs, diags)
  res_tab[, 1:5] <- sapply(res_tab[, 1:5], as.numeric)
  
  # stop the clock
  if (time) {
    proc_time <- proc.time() - ptm
  } else {
    proc_time <- NULL
  }
  
  out <- list(res_tab = res_tab, fits = fits, 
              data = metab, basis = basis, amp_cols = ncol(amps), 
              proc_time = proc_time, method = method, opts = opts)
  
  class(out) <- "fit_result"
  return(out)
}

tarquin_fit <- function(element, temp_mrs, basis_file, opts) {
  metab <- temp_mrs
  # write the metabolite data to file
  if (length(dim(element)) == 2) {
    metab$data[1, 1, 1, 1, 1, 1,] <- element[1,]
  } else if (length(dim(element)) == 0) {
    metab$data[1, 1, 1, 1, 1, 1,] <- element
    ref_file <- NA
  } else {
    stop('Unexpected data dimensions.')
  }
  
  metab_file <- tempfile()
  write_mrs_dpt_v2(metab_file, metab)
  
  if (length(dim(element)) == 2) {
    # looks like there is reference data so let's write it to file
    ref <- temp_mrs
    ref$data[1, 1, 1, 1, 1, 1,] <- element[2,]
    ref_file <- tempfile()
    write_mrs_dpt_v2(ref_file, ref)
    opts = c(opts, "--input_w", ref_file)
  }
  
  # output files
  result_csv_f <- tempfile()
  result_fit_f <- tempfile()
  
  # run TARQUIN
  cmd = paste(getOption("spant.tqn_cmd"), "--input", metab_file, "--format dpt",
              "--basis_lcm", basis_file, "--output_csv", result_csv_f,
              "--output_fit", result_fit_f, paste(opts, collapse = " "))
  
  
  # This one doesn't work on windows for some reason.
  #res = system(cmd, intern = TRUE, ignore.stderr = TRUE, ignore.stdout = TRUE)
  
  res = system(cmd, intern = TRUE, ignore.stderr = TRUE)
  
  if (!file.exists(result_csv_f)) {
    print(res)
    print(cmd)
    stop("Error loading data with above TARQUIN command.")
  }
  
  res_tab <- read_tqn_result(result_csv_f)
  fit <- read_tqn_fit(result_fit_f)
  
  return(append(res_tab,list(fit = fit, metab_file = metab_file, 
                             ref_file = ref_file, result_csv_f = result_csv_f, 
                             result_fit_f = result_fit_f)))  
}

lcmodel_fit <- function(element, temp_mrs, basis_file, opts) {
  metab <- temp_mrs
  # write the metabolite data to file
  if (length(dim(element)) == 2) {
    metab$data[1, 1, 1, 1, 1, 1,] <- element[1,]
  } else if (length(dim(element)) == 0) {
    metab$data[1, 1, 1, 1, 1, 1,] <- element
    ref_file <- NA
  } else {
    stop('Unexpected data dimensions.')
  }
  
  metab_file <- tempfile()
  write_mrs_lcm_raw(metab_file, metab)
  
  if (length(dim(element)) == 2) {
    # looks like there is reference data so let's write it to file
    ref <- temp_mrs
    ref$data[1, 1, 1, 1, 1, 1,] <- element[2,]
    ref_file <- tempfile()
    write_mrs_lcm_raw(ref_file, ref)
    opts = c(opts, "dows=T", paste0("filh2o='", ref_file, "'"))
  }
  #print(opts)
  
  # output files
  coord_f <- tempfile()
  #coord_f <- "test.coord"
  # control file
  control_f <- tempfile()
  #print(coord_f)
  
  sink(control_f)
  cat("$LCMODL\n")
  cat(paste0("nunfil=", dim(metab$data)[7], "\n"))
  cat(paste0("deltat=", metab$resolution[7],"\n"))
  cat(paste0("hzpppm=", metab$ft / 1e6, "\n"))
  cat(paste0("filbas='", basis_file, "'\n"))
  cat(paste0("filraw='", metab_file, "'\n"))
  cat("lcoord=9\n")                        # output a coord file
  cat(paste0("filcoo='", coord_f, "'\n"))  # coord filename
  cat("lps=0\n")                           # don't output a ps plot
  #lcm_control.write("lps=7\n")            # output a ps plot
  #lcm_control.write("filps='"+ps_f+"'\n") # ps plot filename
  cat("neach=999\n")                       # plot all metabs
  #cat("nsimul=11\n")                       # default is 13 which includes Gua 
                                           # and the -CrCH2 signal
  for (opt in opts) { # append any extra options
    cat(paste0(opt, "\n"))
  }
  
  cat("$END")
  sink()
  
  # used for debugging
  #file.copy(control_f, "~/control.file")
  
  # run LCModel
  cmd <- paste(getOption("spant.lcm_cmd"), "<", control_f)
  res = system(cmd, intern = TRUE, ignore.stderr = TRUE, ignore.stdout = TRUE)
  
  #print(cmd)
  if (!file.exists(coord_f)) {
    print(res)
    print(cmd)
    stop("Error with above LCMODEL command.")
  }
  
  # used for debugging
  #file.copy(coord_f, "~/coord.file")
  
  coord_res <- read_lcm_coord(coord_f)
  res_tab <- coord_res$res_tab
  fit <- coord_res$fit
  
  return(append(res_tab, list(fit = fit, metab_file = metab_file, 
                             ref_file = ref_file, coord_f = coord_f, 
                             control_f = control_f)))  
}

#' Reader for csv results generated by TARQUIN.
#' @param result_f TARQUIN result file.
#' @param remove_rcs omit row, column and slice ids from output.
#' @return list of amplitudes, crlbs and diagnostics.
#' @examples
#' \dontrun{
#' result <- read_tqn_result(system.file("extdata","result.csv",package="spant"))
#' }
#' @export
read_tqn_result <- function(result_f, remove_rcs=TRUE) {
  amps <- utils::read.csv(result_f, skip = 1, nrow = 1)
  crlbs <- utils::read.csv(result_f, skip = 4, nrow = 1)
  diags <- utils::read.csv(result_f, skip = 7, nrow = 1)
  if (remove_rcs) {
    amps <- amps[,c(-1, -2, -3)]
    crlbs <- crlbs[,c(-1, -2, -3)]
    diags <- diags[,c(-1, -2, -3)]
  }
  return(list(amps = amps, crlbs = crlbs, diags = diags))
}

#' Reader for csv fit results generated by TARQUIN.
#' @param fit_f TARQUIN fit file.
#' @return A data frame of the fit data points.
#' @examples
#' \dontrun{
#' fit <- read_tqn_fit(system.file("extdata","fit.csv",package="spant"))
#' }
#' @export
read_tqn_fit <- function(fit_f) {
  fit <- utils::read.csv(fit_f, skip = 1)
  class(fit) <- c("fit_table", "data.frame")
  return(fit)
}

#' Read an LCModel formatted coord file containing fit information.
#' @param coord_f path to the coord file.
#' @return list containing a table of fit point and results structure containing
#' signal amplitudes, errors and fitting diagnostics.
#' @export
read_lcm_coord <- function(coord_f) {
  line_reader <- readLines(coord_f)
  n <- 0
  for (line in line_reader) {
    n <- n + 1
    if (endsWith(line, "lines in following concentration table = NCONC+1")) {
      signals <- as.integer(strsplit(trimws(line)," ")[[1]][1]) - 1
      #print(signals)
    } else if (endsWith(line, "lines in following misc. output table")) {
      misc <- as.integer(strsplit(trimws(line)," ")[[1]][1])
      #print(misc)
    } else if (endsWith(line,"points on ppm-axis = NY")) {
      data_start = n
      points <- as.integer(strsplit(trimws(line)," ")[[1]][1])
      #print(data_start)
    }
  }
  
  FWHM <- as.double(strsplit(trimws(line_reader[signals + 8]),"  *")[[1]][3])
  SNR <- as.double(strsplit(trimws(line_reader[signals + 8]),"  *")[[1]][7])
  diags <- data.frame(FWHM = FWHM, SNR = SNR)
  
  #print(coord_f)  
  # -1 width needed to avoid issues when the metab name is
  # prefixed with a + or -
  metab_table <- utils::read.fwf(coord_f, widths = c(9, 5, 8, -1, 40), skip = 6, 
                                 n = signals, header = FALSE, 
                                 col.names = c("amp", "SD", "TCr_ratio", "Metab")
                                 , row.names = "Metab")
  
  row.names(metab_table) <- trimws(row.names(metab_table))
  
  metab_table$SD <- as.double(sub("%", "", metab_table$SD)) # convert to doubles
  metab_table <- t(metab_table)
  amps <- as.data.frame(t(as.data.frame(metab_table[1,])))
  row.names(amps) <- "1"
  crlbs <- as.data.frame(t(as.data.frame(metab_table[2,])))
  max_crlbs <- crlbs == 999
  crlbs <- amps*crlbs/100
  crlbs[max_crlbs] = Inf
  row.names(crlbs) <- "1"
  res_tab <- list(amps = amps, crlbs = crlbs, diags = diags)
  
  data_lines <- ceiling(points/10)
  #print(data_lines)
  n <- 0
  colnames <- vector()
  fit_tab_list <- list()
  repeat {
    header_line <- line_reader[data_start + n * (data_lines + 1)]
    if ( endsWith(header_line,"lines in following diagnostic table:") ) {break}
    name <- strsplit(trimws(header_line),"  *")[[1]][1]
    skip_pt <- data_start + (n * (data_lines + 1))
    data <- as.vector(t(as.matrix(utils::read.table(coord_f, skip = skip_pt, 
                                                    nrows = data_lines,
                                                    fill = T))))
    
    fit_tab_list <- c(fit_tab_list, list(data))
    colnames <- c(colnames, name)
    n = n + 1
  }
  colnames[1] = "PPMScale"
  colnames[2] = "Data"
  colnames[3] = "Fit"
  colnames[4] = "Baseline"
  names(fit_tab_list) <- colnames
  fit_tab <- stats::na.omit(as.data.frame(fit_tab_list))
  fit_tab$Fit <- fit_tab$Fit - fit_tab$Baseline
  fit_tab[5:ncol(fit_tab)] <- fit_tab[5:ncol(fit_tab)] - fit_tab$Baseline 
  class(fit_tab) <- c("fit_table", "data.frame")
  
  return(list(fit = fit_tab, res_tab = res_tab))
}

test_fit <- function(method = "VARPRO_3P", n = 10, preproc = TRUE) {
  # real data 
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
  mrs_data <- read_mrs(fname, format = "spar_sdat")
  if (preproc) {
    mrs_data <- hsvd_filt(mrs_data)
    mrs_data <- align(mrs_data, 2.01)
  }
  basis <- sim_basis_1h_brain_press(mrs_data)
  fit <- fit_mrs(mrs_data, basis, method = method)
  graphics::plot(fit, xlim = c(4,0.5))
  system.time(replicate(n, fit_mrs(mrs_data, basis, method = method)))
}

dim.fit_result <- function(x) {
  dim(x$data$data)
}

#' Calculate diagnostic information for object of class \code{fit_result}.
#' @param x \code{fit_result} object.
#' @param amps known metabolite amplitudes.
#' @return a dataframe of diagnostic information.
#' @export
fit_diags <- function(x, amps = NULL) {
  time <- as.numeric(x$proc_time[3])
  mean_iters <- mean(x$res_tab$res.niter)
  mean_res <- mean(x$res_tab$res.deviance)
  if (!is.null(amps)) {
    basis_n <- length(x$basis$names)
    amps <- matrix(amps, nrow = nrow(x$res_tab), ncol = basis_n, byrow = TRUE)
    mean_error <- mean(rowSums((amps - x$res_tab[6:(5 + basis_n)])^2))
    data.frame(duration = time, mean_iters = mean_iters, mean_res = mean_res,
               mean_error = mean_error)
  } else {
    data.frame(duration = time, mean_iters = mean_iters, mean_res = mean_res)
  }
}

#' Extract the fit amplitudes from an object of class \code{fit_result}.
#' @param x \code{fit_result} object.
#' @param inc_index include columns for the voxel index.
#' @param sort_names sort the basis set names alphabetically.
#' @param append_common_1h_comb append commonly used 1H metabolite combinations
#' eg tNAA = NAA + NAAG.
#' @return a dataframe of amplitudes.
#' @export
fit_amps <- function(x, inc_index = FALSE, sort_names = FALSE,
                     append_common_1h_comb = TRUE) {
  
  basis_n <- length(x$basis$names)
  out <- x$res_tab[, 6:(5 + basis_n)]
  if (sort_names) out <- out[, order(colnames(out))]
  if (inc_index)  out <- cbind(x$res_tab[, 1:5], out)
  
  if (append_common_1h_comb) {
    # create some common metabolite combinations
    if (("NAA" %in% colnames(out)) & ("NAAG" %in% colnames(out))) {
      out['tNAA'] <- out['NAA'] + out['NAAG']
    }
    
    if (("Cr" %in% colnames(out)) & ("PCr" %in% colnames(out))) {
      out['tCr'] <- out['Cr'] + out['PCr']
    }
    
    if (("PCh" %in% colnames(out)) & ("GPC" %in% colnames(out))) {
      out['tCho'] <- out['PCh'] + out['GPC']
    }
    
    if (("Glu" %in% colnames(out)) & ("Gln" %in% colnames(out))) {
      out['Glx'] <- out['Glu'] + out['Gln']
    }
    
    if (("Lip09" %in% colnames(out)) & ("MM09" %in% colnames(out))) {
      out['tLM09'] <- out['Lip09'] + out['MM09']
    }
    
    if (("Lip13a" %in% colnames(out)) & ("Lip13b" %in% colnames(out)) & 
        ("MM12" %in% colnames(out)) & ("MM14" %in% colnames(out))) {
      out["tLM13"] <- out["Lip13a"] + out["Lip13b"] + out["MM12"] + out["MM14"]
    }
    
    if (("Lip20" %in% colnames(out)) & ("MM20" %in% colnames(out))) {
      out['tLM20'] <- out['Lip20'] + out['MM20']
    }
  }
  out
}

#' Combine a list of fit_result objects into a single fit object containing a 
#' dynamic series.
#' @param fit_list list of fit_result objects.
#' @return fit_result object.
#' @export
comb_fits <- function(fit_list) {
  res_tab <- do.call(rbind, lapply(fit_list, `[[`, "res_tab"))
  res_tab$Dynamic <- 1:length(fit_list)
  fits <- unlist(lapply(fit_list, `[[`, "fits"), FALSE)
  data <- append_dyns(lapply(fit_list, `[[`, "data"))
  proc_time <- fit_list[[1]]$proc_time
  for (n in 2:length(fit_list)) {
    proc_time <- proc_time + fit_list[[n]]$proc_time
  }
  
  fit_out <- list(res_tab = res_tab, fits = fits, data = data,
                  basis = fit_list[[1]]$basis,
                  amp_cols = fit_list[[1]]$amp_cols, proc_time = proc_time)
  
  class(fit_out) <- "fit_result"
  fit_out
}

#' Combine all fitting data points into a single dataframe.
#' @param fit_res a single fit_result object.
#' @return a dataframe containing the fit data points.
#' @export
get_fit_table <- function(fit_res) {
  
  # search for masked spectra
  na_fits <- is.na(fit_res$fits)
  
  if (sum(na_fits) > 0) {
    fit_n <- which(!na_fits)[1]
    # create a df with consistent dimensions containing only NA's
    na_frame <- fit_res$fits[[fit_n]]
    na_frame[] <- NA
    fits_w_na <- fit_res$fits
    fits_w_na[na_fits] <- list(na_frame)
    frows <- nrow(fits_w_na[[1]])
    full_fit_df <- do.call("rbind", fits_w_na)
  } else {
    frows <- nrow(fit_res$fits[[1]])
    full_fit_df <- do.call("rbind", fit_res$fits)
  }
  
  # add some ID columns
  full_fit_df$id <- rep(1:length(na_fits), each = frows)
  full_fit_df$X  <- rep(fit_res$res_tab$X, each = frows)
  full_fit_df$Y  <- rep(fit_res$res_tab$Y, each = frows)
  full_fit_df$Z  <- rep(fit_res$res_tab$Z, each = frows)
  full_fit_df$Dynamic <- rep(fit_res$res_tab$Dynamic, each = frows)
  full_fit_df$Coil <- rep(fit_res$res_tab$Coil, each = frows)
  full_fit_df
}