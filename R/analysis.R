#' Perform a fit based analysis of MRS data.
#' 
#' Note that TARQUIN and LCModel require these packages to be installed, and 
#' the functions set_tqn_path and set_lcm_path (respectivly) need to be used 
#' to specify the location of these software packages.
#'
#' Fitting approaches decribed in the following references:
#' VARPRO
#' van der Veen JW, de Beer R, Luyten PR, van Ormondt D. Accurate quantification
#' of in vivo 31P NMR signals using the variable projection method and prior 
#' knowledge. Magn Reson Med 1988;6:92-98
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
#' @param mrs_data An MRS data object.
#' @param basis A basis class object or character vector to basis file in 
#' LCModel .basis format.
#' @param method 'VARPRO', 'TARQUIN' or 'LCMODEL'
#' @param opts Options to pass to the analysis method.
#' @param parallel Perform analysis in parallel (TRUE or FALSE)
#' @param cores Number of cores to use for parallel analysis.
#' @return MRS analysis object.
#' @examples
#' fname <- system.file("extdata","philips_spar_sdat_WS.SDAT",package="spant")
#' svs <- read_mrs(fname, format="spar_sdat")
#' \dontrun{
#' fit_result <- fit_mrs(svs, basis)
#' }
#' @export
fit_mrs <- function(mrs_data, basis, method = 'VARPRO', opts = NULL, 
                    parallel = FALSE, cores = 4) {
  
  if (class(basis) == "basis_set") {
    # TODO check basis matches mrs_data 
  } else if (is.character(basis)) {
    if (!file.exists(basis)) {
      stop("Error, basis file not found.")
    }
  } else {
    stop("Error, specified basis is not the correct data type.")
  }
  
  #if (parallel) { registerDoSNOW(makeCluster(cores, type="SOCK")) }
  #if (parallel) { doMC::registerDoMC(cores=cores) }
  
  # force uppercase for matching
  method <- toupper(method)
  
  if (method == "VARPRO") {
    # put data into TD
    if (is_fd(mrs_data)) {
      mrs_data <- fd2td(mrs_data)
    }
    
    # read basis into memory if a file
    if (is.character(basis)) {
      basis <- read_basis(basis)
    }
    
    if (is.null(opts)) {
      opts <- varpro_opts()
    }
    
    temp_mrs <- mrs_data
    temp_mrs$data = temp_mrs$data[1, 1, 1, 1, 1, 1,]
    dim(temp_mrs$data) <- c(1, 1, 1, 1, 1, 1, length(temp_mrs$data))
  
    #result_list <- apply(mrs_data$data, c(2,3,4,5,6), varpro_fit, temp_mrs, 
    #                     basis, opts)
    
    result_list <- plyr::alply(mrs_data$data, c(2, 3, 4, 5, 6), varpro_fit, 
                               temp_mrs, basis, opts, 
                               .parallel = parallel, 
                               .paropts = list(.inorder = TRUE),
                               .progress = "text", .inform = FALSE)
    
  } else if (method == "VARPRO_3P") {
    # put data into TD
    if (is_fd(mrs_data)) {
      mrs_data <- fd2td(mrs_data)
    }
    
    # read basis into memory if a file
    if (is.character(basis)) {
      basis <- read_basis(basis)
    }
    
    if (is.null(opts)) {
      opts <- varpro_3_para_opts()
    }
    
    temp_mrs <- mrs_data
    temp_mrs$data = temp_mrs$data[1, 1, 1, 1, 1, 1,]
    dim(temp_mrs$data) <- c(1, 1, 1, 1, 1, 1, length(temp_mrs$data))
  
    #result_list <- apply(mrs_data$data, c(2,3,4,5,6), varpro_fit, temp_mrs, 
    #                     basis, opts)
    
    result_list <- plyr::alply(mrs_data$data, c(2, 3, 4, 5, 6),
                               varpro_3_para_fit, temp_mrs, basis, opts, 
                               .parallel = parallel, 
                               .paropts = list(.inorder = TRUE),
                               .progress = "text", .inform = FALSE)
    
    
  } else if (method == "TARQUIN") {
    # write basis object (if specified) to file
    basis_file <- tempfile(fileext = ".basis")
    write_basis(basis, basis_file)
    
    temp_mrs <- mrs_data
    temp_mrs$data = temp_mrs$data[1, 1, 1, 1, 1, 1,]
    dim(temp_mrs$data) <- c(1, 1, 1, 1, 1, 1, length(temp_mrs$data))
  
    #result_list <- apply(mrs_data$data, c(2,3,4,5,6), tarquin_fit, temp_mrs, 
    #                     basis_file, opts)
    
    result_list <- plyr::alply(mrs_data$data, c(2, 3, 4, 5, 6), tarquin_fit, 
                               temp_mrs, basis_file, opts, 
                               .parallel = parallel, 
                               .paropts = list(.inorder = TRUE),
                               .progress = "text", .inform = FALSE)
    
                         #.paropts = list(.options.snow=snowopts),
                         #.paropts = list(.export="N",.packages="spant"),
  } else if (method == "LCMODEL") {
    # write basis object (if specified) to file
    basis_file <- tempfile(fileext = ".basis")
    write_basis(basis, basis_file)
    
    temp_mrs <- mrs_data
    temp_mrs$data = temp_mrs$data[1, 1, 1, 1, 1, 1,]
    dim(temp_mrs$data) <- c(1, 1, 1, 1, 1, 1, length(temp_mrs$data))
    
    result_list <- plyr::alply(mrs_data$data, c(2,3,4,5,6), lcmodel_fit, 
                               temp_mrs, basis_file, opts,
                               .parallel = parallel, 
                               .paropts = list(.inorder = TRUE),
                               .progress = "text", .inform = FALSE)
  } else {
    stop(paste('Error, unsupported fit method : ', method))
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
    if (sum(duplicated(file_list_vec)) > 0 ) {
      stop("Error, duplicate temp file detected.")
    }
    #print(file_list_vec)
  }
  
  df_list_amps <- result_list[seq(from = 1, by = res_n,length.out = fit_num)]
  amps = plyr::ldply(df_list_amps, data.frame)[-1]
  df_list_crlbs <- result_list[seq(from = 2, by = res_n, length.out = fit_num)]
  crlbs = plyr::ldply(df_list_crlbs, data.frame)[-1]
  colnames(crlbs) <- paste(colnames(crlbs), ".sd", sep = "")
  df_list_diags <- result_list[seq(from = 3, by = res_n, length.out = fit_num)]
  diags = plyr::ldply(df_list_diags, data.frame)[-1]
  fits <- result_list[seq(from = 4, by = res_n, length.out = fit_num)]
  
  out <- list(results = cbind(labs, amps, crlbs, diags), fits = fits, 
              data = mrs_data, amp_cols = ncol(amps))
  
  class(out) <- "analysis_results"
  return(out)
}

varpro_fit <- function(element, temp_mrs, basis, opts) {
  metab <- temp_mrs
  metab$data[1, 1, 1, 1, 1, 1,] <- element
  varpro(metab, basis, opts)
}

varpro_3_para_fit <- function(element, temp_mrs, basis, opts) {
  metab <- temp_mrs
  metab$data[1, 1, 1, 1, 1, 1,] <- element
  varpro_3_para(metab, basis, opts)
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
  
  res = system(cmd, intern = TRUE, ignore.stderr = TRUE, ignore.stdout = TRUE)
  
  if (!file.exists(result_csv_f)) {
    print(res)
    print(cmd)
    stop("Error loading data with above TARQUIN command.")
  }
  
  results <- read_tqn_result(result_csv_f)
  fit <- read_tqn_fit(result_fit_f)
  
  return(append(results,list(fit = fit, metab_file = metab_file, 
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
  cat("nsimul=11\n")                       # default is 13 which includes Gua 
                                           # and the -CrCH2 signal
  for (opt in opts) { # append any extra options
    cat(paste0(opt, "\n"))
  }
  
  cat("$END")
  sink()
  
  # run LCModel
  cmd <- paste(getOption("spant.lcm_cmd"), "<", control_f)
  res = system(cmd, intern = TRUE, ignore.stderr = TRUE, ignore.stdout = TRUE)
  
  #print(cmd)
  if (!file.exists(coord_f)) {
    print(res)
    print(cmd)
    stop("Error with above LCMODEL command.")
  }
  
  coord_res <- read_lcm_coord(coord_f)
  results <- coord_res$results
  fit <- coord_res$fit
  
  return(append(results,list(fit = fit, metab_file = metab_file, 
                             ref_file = ref_file, coord_f = coord_f, 
                             control_f = control_f)))  
}

#' Reader for csv results generated by TARQUIN.
#' @param result_f TARQIUN result file.
#' @param remove_rcs Omit row, column and slice ids from output.
#' @return A list of amplitudies, crlbs and diagnostics.
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
#' @param fit_f TARQUIN fit file
#' @return A data frame of the fit data points.
#' @examples
#' \dontrun{
#' fit <- read_tqn_fit(system.file("extdata","fit.csv",package="spant"))
#' }
#' @export
read_tqn_fit <- function(fit_f) {
  fit <- utils::read.csv(fit_f, skip = 1)
  class(fit) <- "fit_table"
  return(fit)
}

#' Read an LCModel formatted coord file contaning fit information.
#' @param coord_f Path to the coord file.
#' @return A list containg a table of fit point and results structure contaning
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
  results <- list(amps = amps, crlbs = crlbs, diags = diags)
  
  data_lines <- ceiling(points/10)
  #print(data_lines)
  n <- 0
  colnames <- vector()
  fit_table_list <- list()
  repeat {
    header_line <- line_reader[data_start + n * (data_lines + 1)]
    if ( endsWith(header_line,"lines in following diagnostic table:") ) {break}
    name <- strsplit(trimws(header_line),"  *")[[1]][1]
    skip_pt <- data_start + (n * (data_lines + 1))
    data <- as.vector(t(as.matrix(utils::read.table(coord_f, skip = skip_pt, 
                                                    nrows = data_lines,
                                                    fill = T))))
    
    fit_table_list <- c(fit_table_list, list(data))
    colnames <- c(colnames, name)
    n = n + 1
  }
  colnames[1] = "PPMScale"
  colnames[2] = "Data"
  colnames[3] = "Fit"
  colnames[4] = "Baseline"
  names(fit_table_list) <- colnames
  fit_table <- stats::na.omit(as.data.frame(fit_table_list))
  fit_table$Fit <- fit_table$Fit - fit_table$Baseline
  fit_table[5:ncol(fit_table)] <- fit_table[5:ncol(fit_table)] - fit_table$Baseline 
  class(fit_table) <- "fit_table"
  
  return(list(fit = fit_table, results = results))
}

get_corr_factor <- function(te, tr, B0, gm_vol, wm_vol, csf_vol) {
  # Correction factor calcualted according to the method of Gasparovic et al (MRM 55:1219-1226 2006)
  # NOTE - gives concs as Mol/kg of water NOT Mol/liter of tissue like default LCM/TQN analysis.
  if (B0 == 3.0) {
    # Wanasapura values given in Harris paper
    t1_gm    = 1.331
    t2_gm    = 0.110
    t1_wm    = 0.832
    t2_wm    = 0.0792
    t1_csf   = 3.817
    t2_csf   = 0.503
    t1_metab = 1.15
    t2_metab = 0.3
  } else if (B0 == 1.5) {
    # values from Gasparovic 2006 MRM paper
    t1_gm    = 1.304
    t2_gm    = 0.093
    t1_wm    = 0.660
    t2_wm    = 0.073
    t1_csf   = 2.93
    t2_csf   = 0.23
    t1_metab = 1.15
    t2_metab = 0.3
  } else {
    stop("Error. Relaxation values not available for this field strength.")
  }
  
  # MR-visable water densities
  gm_vis  = 0.78
  wm_vis  = 0.65
  csf_vis = 0.97
  
  # molal concentration (moles/gram) of MR-visible water
  water_conc = 55510.0
  
  # fractions of water attributable to GM, WM and CSF
  f_gm  = gm_vol * gm_vis / 
          (gm_vol * gm_vis + wm_vol * wm_vis + csf_vol * csf_vis)
  
  f_wm  = wm_vol * wm_vis / 
          (gm_vol * gm_vis + wm_vol * wm_vis + csf_vol * csf_vis)
  
  f_csf = csf_vol * csf_vis / 
          (gm_vol * gm_vis + wm_vol * wm_vis + csf_vol * csf_vis)
  
  #This might give the result in Mol/kg?
  #f_gm  = gm_vol  * gm_vis   / ( gm_vol + wm_vol + csf_vol )
  #f_wm  = wm_vol  * wm_vis   / ( gm_vol + wm_vol + csf_vol )
  #f_csf = csf_vol * csf_vis  / ( gm_vol + wm_vol + csf_vol )
  
  # Relaxtion attenuation factors
  R_h2o_gm  = exp(-te / t2_gm) * (1.0 - exp(-tr / t1_gm))
  R_h2o_wm  = exp(-te / t2_wm) * (1.0 - exp(-tr / t1_wm))
  R_h2o_csf = exp(-te / t2_csf) * (1.0 - exp(-tr / t1_csf))
  R_metab   = exp(-te / t2_metab) * (1.0 - exp(-tr / t1_metab))
  
  corr_factor = ((f_gm * R_h2o_gm + f_wm * R_h2o_wm + f_csf * R_h2o_csf) / 
                ((1 - f_csf) * R_metab)) * water_conc
  
  return(corr_factor)
}

#' @export
int_spec <- function(mrs_data, xlim = NULL, scale = "ppm", mode = "real") {
  
  if (!is_fd(mrs_data)) {
    mrs_data <- td2fd(mrs_data)
  }
    
  if ( scale == "ppm" ) {
    x_scale <- ppm(mrs_data)
  } else if (scale == "hz") {
    x_scale <- hz(mrs_data)
  } else if (scale == "points") {
    x_scale <- pts(mrs_data)
  }
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[N(mrs_data)])
  }
  
  x_inds <- get_seg_ind(x_scale, xlim[1], xlim[2])
  subset <- x_inds[1]:x_inds[2]
  
  data_arr <- mrs_data$data[,,,,,, subset, drop = F]
  
  if (mode == "real") {
    data_arr <- Re(data_arr)
  } else if (mode == "imag") {
    data_arr <- Im(data_arr)
  } else if (mode == "abs") {
    data_arr <- Mod(data_arr)
  }
  
  apply(data_arr, c(1, 2, 3, 4, 5, 6), sum)
}

crop_range <- function(map, lower, upper) {
  # TODO check upper > lower and both are >0<1
  map_range <- range(map, na.rm = TRUE)
  #upper_lim <- map_range[1]+upper/100*(map_range[2] - map_range[1])
  #lower_lim <- map_range[1]+lower/100*(map_range[2] - map_range[1])
  upper_lim <- upper
  lower_lim <- lower
  map <- ifelse(map > upper_lim, upper_lim,map)  
  ifelse(map < lower_lim, lower_lim, map)  
}

test_varpro <- function() {
  # real data 
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
  mrs_data <- read_mrs(fname,format = "spar_sdat")
  mrs_data <- hsvd_filt(mrs_data)
  mrs_data <- align(mrs_data, 2.01)
  acq_paras <- get_acq_paras(mrs_data)
  basis <- sim_basis_1h_brain_press(acq_paras, xlim = c(4,0))
  fit <- fit_mrs(mrs_data, basis)
  plot(fit, xlim = c(6,0))
  system.time(replicate(10, fit_mrs(mrs_data, basis)))
}

varpro_3_para <- function(mrs_data, basis, opts = NULL) {
  # use default fitting opts if not specified 
  if (is.null(opts)) {
      opts <- varpro_opts()
  }
  
  # convert basis from FD to TD
  basis_td <- apply(basis$data, 2, ift_shift)
  
  y <- drop(mrs_data$data)
  Npts <- length(y)
  Nbasis <- dim(basis$data)[2]
  
  # phase, global damping, global shift
  par <- c(0, opts$init_damping, 0)
           
  t <- seconds(mrs_data)  
  # lm control options
  ctrl <- minpack.lm::nls.lm.control()
  ctrl$maxiter = opts$maxiters
  # do the fit
  lower <- c(-pi, 0, -opts$max_shift)
  upper <- c(pi, opts$max_damping, opts$max_shift)
  
  if (opts$anal_jac) {
    res <- minpack.lm::nls.lm(par, lower, upper, varpro_3_para_fn,
                              varpro_3_para_anal_jac, ctrl, y, basis_td, t,
                              opts$nstart)
  } else {
    res <- minpack.lm::nls.lm(par, lower, upper, varpro_3_para_fn, NULL, ctrl, y,
                              basis_td, t, opts$nstart)
  }
  
  # apply phase to y
  y <- y * exp(1i * (res$par[1]))
  
  # apply global broadening term to basis
  basis_mod <- basis_td * matrix(exp(-t * t * lw2beta(res$par[2])),
                                 ncol = ncol(basis_td), nrow = nrow(basis_td),
                                 byrow = FALSE)
  
  # apply shift terms to basis
  t_mat <- matrix(t, nrow = Npts, ncol = Nbasis)
  freq_vec <- 2i * pi * rep(res$par[3], Nbasis)
  freq_mat <- matrix(freq_vec, nrow = Npts, ncol = Nbasis, byrow = TRUE)
 
  basis_mod <- basis_mod * exp(t_mat * freq_mat)
  
  # get ahat
  y_real <- c(Re(y[opts$nstart:Npts]), Im(y[opts$nstart:Npts]))
  basis_real <- rbind(Re(basis_mod[opts$nstart:Npts,]),
                      Im(basis_mod[opts$nstart:Npts,]))
  
  ahat <- nnls::nnls(basis_real, y_real)$x
  
  yhat <- basis_mod %*% ahat
  
  # zero pad
  yhat <- c(yhat, rep(0, Npts))
  YHAT <- ft_shift(as.vector(yhat))
  
  # zero pad
  y <- c(y, rep(0, Npts))
  Y <- ft_shift(y)
  resid <- Y - YHAT
  
  BL <- smoother::smth.gaussian(Re(resid), opts$bl_smth_pts, tails = TRUE) + 
        1i * smoother::smth.gaussian(Im(resid), opts$bl_smth_pts, tails = TRUE)
  
  RESID <- Y - YHAT
  
  offset <- max(Re(Y)) - min(Re(RESID))
  resid <- vec2mrs_data(RESID + offset, fd = TRUE)
  
  amps <- data.frame(t(ahat))
  colnames(amps) <- basis$names
  
  fit <- data.frame(PPMScale = ppm(mrs_data, N = Npts * 2), Data = Re(Y),
                    Fit = Re(YHAT), Baseline = Re(BL))
  
  class(fit) <- "fit_table"
  
  diags <- data.frame(phase = res$par[1] * 180 / pi, damping = res$par[2],
                      shift = res$par[3], res$deviance, res$niter, res$info, 
                      res$deviance, res$message)
  
  list(amps = amps, crlbs = t(rep(NA, length(amps))), diags = diags, fit = fit)
}

varpro <- function(mrs_data, basis, opts = NULL) {
  # use default fitting opts if not specified 
  if (is.null(opts)) {
      opts <- varpro_opts()
  }
  
  # convert basis from FD to TD
  basis_td <- apply(basis$data, 2, ift_shift)
  
  y <- drop(mrs_data$data)
  Npts <- length(y)
  Nbasis <- dim(basis$data)[2]
  
  # phase, global damping, basis shifts, basis dampings
  par <- c(0, opts$init_g_damping, rep(0, Nbasis), rep(0, Nbasis))
           
  t <- seconds(mrs_data)  
  # lm control options
  ctrl <- minpack.lm::nls.lm.control()
  ctrl$maxiter = opts$maxiters
  # do the fit
  lower <- c(-pi, 0, rep(-opts$max_shift, Nbasis), rep(0, Nbasis))
  upper <- c(pi, opts$max_g_damping, rep(opts$max_shift, Nbasis),
             rep(opts$max_ind_damping, Nbasis))
  
  if (opts$anal_jac) {
    res <- minpack.lm::nls.lm(par, lower, upper, varpro_fn, varpro_anal_jac,
                              ctrl, y, basis_td, t, opts$nstart)
  } else {
    res <- minpack.lm::nls.lm(par, lower, upper, varpro_fn, NULL, ctrl, y,
                              basis_td, t, opts$nstart)
  }
  
  # apply phase to y
  y <- y * exp(1i * (res$par[1]))
  
  # apply global broadening term to basis
  basis_mod <- basis_td * matrix(exp(-t * t * lw2beta(res$par[2])),
                                 ncol = ncol(basis_td), nrow = nrow(basis_td),
                                 byrow = FALSE)
  
  # apply shift and lb terms to basis
  t_mat <- matrix(t, nrow = Npts, ncol = Nbasis)
  freq_vec <- 2i * pi * res$par[3:(2 + Nbasis)]
  lb_vec <- lw2alpha(res$par[(3 + Nbasis):(2 + 2 * Nbasis)])
  freq_lb_mat <- matrix(freq_vec - lb_vec, nrow = Npts, ncol = Nbasis,
                        byrow = TRUE) 
  
  basis_mod <- basis_mod * exp(t_mat * freq_lb_mat)
  
  # get ahat
  y_real <- c(Re(y[opts$nstart:Npts]), Im(y[opts$nstart:Npts]))
  basis_real <- rbind(Re(basis_mod[opts$nstart:Npts,]),
                      Im(basis_mod[opts$nstart:Npts,]))
  
  ahat <- nnls::nnls(basis_real, y_real)$x
  
  yhat <- basis_mod %*% ahat
  
  # zero pad
  yhat <- c(yhat, rep(0, Npts))
  YHAT <- ft_shift(as.vector(yhat))
  
  # zero pad
  y <- c(y, rep(0, Npts))
  Y <- ft_shift(y)
  resid <- Y - YHAT
  
  BL <- smoother::smth.gaussian(Re(resid), opts$bl_smth_pts, tails = TRUE) + 
        1i * smoother::smth.gaussian(Im(resid), opts$bl_smth_pts, tails = TRUE)
  
  RESID <- Y - YHAT
  
  offset <- max(Re(Y)) - min(Re(RESID))
  resid <- vec2mrs_data(RESID + offset, fd = TRUE)
  
  amps <- data.frame(t(ahat))
  colnames(amps) <- basis$names
  
  fit <- data.frame(PPMScale = ppm(mrs_data, N = Npts * 2), Data = Re(Y),
                    Fit = Re(YHAT), Baseline = Re(BL))
  
  class(fit) <- "fit_table"
  
  diags <- data.frame(res$deviance, res$niter, res$info, res$deviance,
                      res$message)
  
  list(amps = amps, crlbs = t(rep(NA, length(amps))), diags = diags, fit = fit)
}

varpro_anal_jac <- function(par, y, basis, t, nstart) {
  Npts <- length(y)
  Nbasis <- dim(basis)[2]
  
  # apply phase to y
  y <- y * exp(1i * (par[1]))
  
  # apply global broadening term to basis
  basis_mod <- basis * matrix(exp(-t * t * lw2beta(par[2])), ncol = ncol(basis),
                              nrow = nrow(basis), byrow = F)
  
  # apply shift and lb terms to basis
  t_mat <- matrix(t, nrow = Npts, ncol = Nbasis)
  freq_vec <- 2i * pi * par[3:(2 + Nbasis)]
  lb_vec <- lw2alpha(par[(3 + Nbasis):(2 + 2 * Nbasis)])
  freq_lb_mat <- matrix(freq_vec - lb_vec, nrow = Npts, ncol = Nbasis,
                        byrow = TRUE) 
  
  basis_mod <- basis_mod * exp(t_mat * freq_lb_mat)
  
  t_cut <- c(t[nstart:Npts], t[nstart:Npts])
  y_real <- c(Re(y[nstart:Npts]), Im(y[nstart:Npts]))
  basis_real <- rbind(Re(basis_mod[nstart:Npts,]), Im(basis_mod[nstart:Npts,]))
  
  #ahat <- nnls::nnls(basis_real, y_real)$x
  ahat <- MASS::ginv(basis_real) %*% y_real
  
  unmod_basis_real <- basis_real %*% ahat
  unmod_basis_cplx <- basis_mod %*% ahat
  
  # find phase and global lw grads 
  phase_jac <- unmod_basis_cplx * 1i
  phase_jac_real <- c(Re(phase_jac[nstart:Npts]), Im(phase_jac[nstart:Npts]))
  g_lw_jac_real <- (unmod_basis_real * t_cut * t_cut * 2 * par[2] * pi ^ 2 / 4
                   / (-log(0.5)))
  
  res_size <- 2 * (Npts - nstart + 1)
  
  shift_jac <- vector(length = res_size * Nbasis)
  lw_jac <- vector(length = res_size * Nbasis)
  
  # individual shifts and dampings
  for (n in 1:Nbasis) {
    #adj_basis <- basis_mod
    #adj_basis[,n] <- -adj_basis[,n]*2i*pi*t
    #adj_basis_comb <- adj_basis %*% ahat 
    #adj_basis_comb <- unmod_basis_cplx - ahat[n]*basis_mod[,n] - ahat[n]*basis_mod[,n]*2i*pi*t
    
    adj_basis_comb <- unmod_basis_cplx - ahat[n] * basis_mod[,n] *
                      (1 + 2i * pi * t)
    
    shift_jac_real_temp <- c(Re(adj_basis_comb[nstart:Npts]),
                             Im(adj_basis_comb[nstart:Npts]))
    
    shift_jac[((n - 1) * res_size + 1):(n * res_size)] <- shift_jac_real_temp
    
    adj_basis <- basis_mod
    adj_basis[,n] <- -adj_basis[,n]*pi*t
    adj_basis_comb <- adj_basis %*% ahat 
    
    #adj_basis_comb <- unmod_basis_cplx - ahat[n]*basis_mod[,n] - ahat[n]*basis_mod[,n]*pi*t
    #adj_basis_comb <- unmod_basis_cplx - ahat[n] * basis_mod[,n] * (1 + pi * t)
    
    lw_jac_real_temp <- c(Re(adj_basis_comb[nstart:Npts]),
                          Im(adj_basis_comb[nstart:Npts]))
    
    lw_jac[((n - 1) * res_size + 1):(n * res_size)] <- lw_jac_real_temp
  }
  
  c(phase_jac_real, g_lw_jac_real, shift_jac, lw_jac)
}
    
# varpro_num_jac <- function(par, y, basis, t, nstart) {
#   numDeriv::jacobian(varpro_num_jac_eval, par, method="simple", side=NULL, method.args=list(), y, basis, t, nstart)
# }

# varpro_num_jac_eval <- function(par, y, basis, t, nstart) {
#   Npts <- length(y)
#   Nbasis <- dim(basis)[2]
#   
#   # apply phase to y
#   y <- y*exp(1i*(par[1]))
#   
#   # apply global broadening term to basis
#   basis_mod <- basis * matrix(exp(-t*t*lw2beta(par[2])),ncol=ncol(basis),nrow=nrow(basis),byrow=F)
#   
#   # apply shift and lb terms to basis
#   t_mat <- matrix(t, nrow=Npts, ncol=Nbasis)
#   freq_vec <- 2i*pi*par[3:(2+Nbasis)]
#   lb_vec <- lw2alpha(par[(3+Nbasis):(2+2*Nbasis)])
#   freq_lb_mat <- matrix(freq_vec-lb_vec, nrow=Npts, ncol=Nbasis, byrow=TRUE) 
#   basis_mod <- basis_mod * exp(t_mat * freq_lb_mat)
#   
#   y_real <- c(Re(y[nstart:Npts]),Im(y[nstart:Npts]))
#   basis_real <- rbind(Re(basis_mod[nstart:Npts,]),Im(basis_mod[nstart:Npts,]))
#   ahat <- nnls::nnls(basis_real, y_real)$x
#   
#   res <- y_real - basis_real %*% ahat
#   
#   res
# }
  
varpro_fn <- function(par, y, basis, t, nstart, sc_res = FALSE) {
  Npts <- length(y)
  Nbasis <- dim(basis)[2]
  
  # apply phase to y
  y <- y * exp(1i * (par[1]))
  
  # apply global broadening term to basis
  basis_mod <- basis * matrix(exp(-t * t * lw2beta(par[2])), ncol = ncol(basis),
                              nrow = nrow(basis), byrow = F)
  
  # apply shift and lb terms to basis
  t_mat <- matrix(t, nrow = Npts, ncol = Nbasis)
  freq_vec <- 2i * pi * par[3:(2 + Nbasis)]
  lb_vec <- lw2alpha(par[(3 + Nbasis):(2 + 2 * Nbasis)])
  freq_lb_mat <- matrix(freq_vec - lb_vec, nrow = Npts, ncol = Nbasis,
                        byrow = TRUE) 
  
  basis_mod <- basis_mod * exp(t_mat * freq_lb_mat)
  
  y_real <- c(Re(y[nstart:Npts]), Im(y[nstart:Npts]))
  basis_real <- rbind(Re(basis_mod[nstart:Npts,]), Im(basis_mod[nstart:Npts,]))
  
  ahat <- nnls::nnls(basis_real, y_real)$x
  res <- y_real - basis_real %*% ahat
  
  if ( sc_res ) {
    res <- sum(res ^ 2)
  }
  res
}

varpro_3_para_fn <- function(par, y, basis, t, nstart, sc_res = FALSE) {
  Npts <- length(y)
  Nbasis <- dim(basis)[2]
  
  # apply phase to y
  y <- y * exp(1i * (par[1]))
  
  # apply global broadening term to basis
  basis_mod <- basis * matrix(exp(-t * t * lw2beta(par[2])), ncol = ncol(basis),
                              nrow = nrow(basis), byrow = F)
  
  # apply global shift to basis
  t_mat <- matrix(t, nrow = Npts, ncol = Nbasis)
  freq_vec <- 2i * pi * rep(par[3], Nbasis)
  freq_mat <- matrix(freq_vec, nrow = Npts, ncol = Nbasis, byrow = TRUE)
 
  basis_mod <- basis_mod * exp(t_mat * freq_mat)
  
  y_real <- c(Re(y[nstart:Npts]), Im(y[nstart:Npts]))
  basis_real <- rbind(Re(basis_mod[nstart:Npts,]), Im(basis_mod[nstart:Npts,]))
  
  ahat <- nnls::nnls(basis_real, y_real)$x
  res <- y_real - basis_real %*% ahat
  
  if ( sc_res ) {
    res <- sum(res ^ 2)
  }
  res
}

#' Return a list of options for VARPRO based fitting.
#' @param nstart Position in the time-domain to start fitting, units of data
#' points.
#' @param init_g_damping Starting value for the global Gaussian line-broadening
#' term - measured in Hz.
#' @param maxiters Maximum number of lemar iterations to perform.
#' @param max_shift Maximum shift allowed to each element in the basis set, 
#' measured in Hz.
#' @param max_g_damping Maximum permitted global Gaussian line-broadening.
#' @param max_ind_damping Maximum permitted Lorentzian line-broadening for each
#' element in the basis set, measured in Hz.
#' @param anal_jac Option to use the analytic or numerical Jacobian (logical).
#' @param bl_smth_pts Number of data points to use in the baseline smoothing
#' calculation.
#' @return List of options.
#' @examples
#' varpro_opts(nstart = 20)
#' @export
varpro_opts <- function(nstart = 10, init_g_damping = 2, maxiters = 200,
                        max_shift = 5, max_g_damping = 5, max_ind_damping = 5,
                        anal_jac = TRUE, bl_smth_pts = 80) {
  
  list(nstart = nstart, init_g_damping = init_g_damping, maxiters = maxiters,
       max_shift = max_shift, max_g_damping = max_g_damping,
       max_ind_damping = max_ind_damping, anal_jac = anal_jac,
       bl_smth_pts = bl_smth_pts)
}

#' Return a list of options for VARPRO based fitting with 3 free parameters:
#' * zero'th order phase correction
#' * global damping
#' * global frequency shift.
#' @param nstart Position in the time-domain to start fitting, units of data
#' points.
#' @param init_damping Starting value for the global Gaussian line-broadening
#' term - measured in Hz.
#' @param maxiters Maximum number of lemar iterations to perform.
#' @param max_shift Maximum global shift allowed, measured in Hz.
#' @param max_damping Maximum damping allowed, FWHM measured in Hz.
#' @param anal_jac Option to use the analytic or numerical Jacobian (logical).
#' @param bl_smth_pts Number of data points to use in the baseline smoothing
#' calculation.
#' @return List of options.
#' @examples
#' varpro_opts(nstart = 20)
#' @export
varpro_3_para_opts <- function(nstart = 10, init_damping = 2, maxiters = 200,
                        max_shift = 5, max_damping = 5, anal_jac = FALSE,
                        bl_smth_pts = 80) {
  
  list(nstart = nstart, init_damping = init_damping, maxiters = maxiters,
       max_shift = max_shift, max_damping = max_damping, anal_jac = anal_jac,
       bl_smth_pts = bl_smth_pts)
}

varpro_3_para_anal_jac <- function(par, y, basis, t, nstart) {
  Npts <- length(y)
  Nbasis <- dim(basis)[2]

  # apply phase to y
  y <- y*exp(1i * (par[1]))

  # apply global broadening term to basis
  basis_mod <- basis * matrix(exp(-t * t * lw2beta(par[2])), ncol = ncol(basis),
                              nrow = nrow(basis), byrow = F)

  # apply global shift to basis
  t_mat <- matrix(t, nrow = Npts, ncol = Nbasis)
  freq_vec <- 2i * pi * rep(par[3], Nbasis)
  freq_mat <- matrix(freq_vec, nrow = Npts, ncol = Nbasis, byrow = TRUE)
  basis_mod <- basis_mod * exp(t_mat * freq_mat)

  t_cut <- c(t[nstart:Npts], t[nstart:Npts])
  y_real <- c(Re(y[nstart:Npts]), Im(y[nstart:Npts]))
  basis_real <- rbind(Re(basis_mod[nstart:Npts,]), Im(basis_mod[nstart:Npts,]))
  ahat <- nnls::nnls(basis_real, y_real)$x
  #ahat <- MASS::ginv(basis_real) %*% y_real
  
  

  unmod_basis_real <- basis_real %*% ahat
  unmod_basis_cplx <- basis_mod %*% ahat

  phase_jac <- unmod_basis_cplx * 1i
  phase_jac_real <- c(Re(phase_jac[nstart:Npts]), Im(phase_jac[nstart:Npts]))
  g_lw_jac_real <- (unmod_basis_real * t_cut * t_cut * 2 * par[2] * (pi ^ 2) / 4
                    / (-log(0.5)))
  
  shift_jac <- -unmod_basis_cplx * 2i * pi * t
  shift_jac_real <- c(Re(shift_jac[nstart:Npts]), Im(shift_jac[nstart:Npts]))

  c(phase_jac_real, g_lw_jac_real, shift_jac_real)
}