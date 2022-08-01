calc_siemens_paras <- function(vars, is_ima) {
  
  # correct fs if the FID has been decimated for ima data
  if (is_ima) {
    if (vars$rm_oversampling) vars$fs <- vars$fs / 2
  }
  
  res <- c(NA, vars$x_dim / vars$x_pts, vars$y_dim / vars$y_pts,
           vars$z_dim / vars$z_pts, 1, NA, 1 / vars$fs)
  
  ima_norm <- c(vars$norm_sag, vars$norm_cor, vars$norm_tra)
  ima_norm <- l2_norm_vec(ima_norm)
  ima_pos  <- c(vars$pos_sag,  vars$pos_cor,  vars$pos_tra)
  rotation <- vars$ip_rot
  
  x <- fGSLCalcPRS(ima_norm, rotation)
  col_vec <- x$dGp
  
  # to get agreement with spec2nii we need to make SVS and MRSI special cases
  if (vars$x_pts * vars$y_pts * vars$z_pts == 1) {
    # SVS
    row_vec <- x$dGr
  } else {
    # MRSI
    row_vec <- -x$dGr
  }
  
  # sli_vec <- ima_norm
  sli_vec <- cross(row_vec, col_vec)
  ima_pos <- ima_pos - row_vec * (vars$x_pts / 2 - 0.5) * vars$x_dim /
                       vars$x_pts - col_vec * (vars$y_pts / 2 - 0.5) *
                       vars$y_dim / vars$y_pts
  
  affine <- cbind(c(row_vec * res[2], 0),
                  c(col_vec * res[3], 0),
                  c(sli_vec * res[4], 0),
                  c(ima_pos, 1))
  
  # NIfTI sign swap
  affine[1:2,] <- -affine[1:2,]
  
  # TODO parse from the data file and use sensible ref based on nuc
  nuc <- def_nuc()
  ref <- def_ref()
  
  return(list(res = res, nuc = nuc, ref = ref,
              affine = affine))
}

read_twix <- function(fname, verbose, full_fid = FALSE,
                      omit_svs_ref_scans = TRUE, extra) {
  
  # check the file size
  fbytes <- file.size(fname)
  
  # read the binary header
  con <- file(fname, "rb")
  first_int <- read_uint32(con)
  second_int <- read_uint32(con)
  #print(first_int)
  #print(second_int)
  if ((first_int < 10000) && (second_int <= 64)) {
    if (verbose) cat("TWIX file is VD format.\n")
    version <- "vd"
    Nscans <- second_int
    measID <- read_uint32(con)
    #print(measID)
    fileID <- read_uint32(con)
    #print(fileID)
    measOffset <- read_uint64(con) # usually 10240 bytes
    #print(measOffset)
    measLength <- read_uint64(con)
    seek(con, measOffset)
    hdrLength <- read_uint32(con)
    dataStart <- measOffset + hdrLength
    #print(dataStart)
  } else {
    if (verbose) cat("TWIX file is VB format.\n")
    version <- "vb"
    dataStart <- first_int
    Nscans <- 1
  }
  close(con)
  vars <- read_siemens_txt_hdr(fname, version, verbose)
  
  # read data points 
  con <- file(fname, "rb")
  seek(con, dataStart)
 
  # scan ids
  inds <- NULL
  
  # ima_echoes <- 0
  cPos <- dataStart
  
  # overestimate the number of data points and pre-allocate
  expected_pts <- as.integer((fbytes - dataStart) / 4)  
  
  if (verbose) cat(paste("Scans           :", Nscans, "\n"))
  
  for (scans in 0:(Nscans - 1)) {
    # the final scan is the one we are interested in, so clear the last one 
    raw_pts <- c(NA)
    inds <- NULL
    length(raw_pts) <- expected_pts
    raw_pt_start <- 1
    
    n <- 0
    while (TRUE) {
      if (verbose) cat(".")
      seek(con, cPos, "start")
      n <- n + 1
      ulDMALength_bin <- intToBits(read_int32(con))
      ulDMALength <- sum(2^.subset(0:24, as.logical(ulDMALength_bin[1:25])))
      #print("ulDMALenth")
      #print(ulDMALength)
      seek(con, cPos, "start")
      if (version == "vb") {
        seek(con, 20, "current") # move ahead 20 bytes
      } else {
        seek(con, 40, "current") # move ahead 40 bytes
      }
      
      # read next 64 bits
      info_bits <- c(as.logical(intToBits(read_int32(con))),
                     as.logical(intToBits(read_int32(con))))
      samples_in_scan       <- read_uint16(con) # NCol
      used_channels         <- read_uint16(con) # NCha
      Lin                   <- read_uint16(con)
      Ave                   <- read_uint16(con)
      Sli                   <- read_uint16(con)
      Par                   <- read_uint16(con)
      Eco                   <- read_uint16(con)
      Phs                   <- read_uint16(con)
      Rep                   <- read_uint16(con)
      Set                   <- read_uint16(con)
      Seg                   <- read_uint16(con)
      Ida                   <- read_uint16(con)
      Idb                   <- read_uint16(con)
      Idc                   <- read_uint16(con)
      Idd                   <- read_uint16(con)
      Ide                   <- read_uint16(con)
      
      # not sure if this next part is right as I don't have any test data
      # that uses it
      seek(con, 4, "current")
      kspace_center_column  <- read_uint16(con) # centerCol
      centerLin             <- read_uint16(con)
      centerPar             <- read_uint16(con)
      cutOff                <- read_uint16(con)
      coilSelect            <- read_uint16(con)
      
      MDH_ACQEND            <- info_bits[1]
      MDH_RTFEEDBACK        <- info_bits[2]
      MDH_HPFEEDBACK        <- info_bits[3]
      MDH_SYNCDATA          <- info_bits[6]
      MDH_RAWDATACORRECTION <- info_bits[11]
      MDH_REFPHASESTABSCAN  <- info_bits[15]
      MDH_PHASESTABSCAN     <- info_bits[16]
      MDH_SIGNREV           <- info_bits[18]
      MDH_PHASCOR           <- info_bits[22]
      MDH_PATREFSCAN        <- info_bits[23]
      MDH_PATREFANDIMASCAN  <- info_bits[24]
      MDH_REFLECT           <- info_bits[25]
      MDH_NOISEADJSCAN      <- info_bits[26]
      MDH_IMASCAN           <- TRUE
      
      #print(info_bits)
      #print(samples_in_scan)
      #print(used_channels)
      
      if (MDH_ACQEND || MDH_RTFEEDBACK || MDH_HPFEEDBACK || MDH_PHASCOR || MDH_NOISEADJSCAN || MDH_SYNCDATA) {
        MDH_IMASCAN <- FALSE
      }
      
      if (MDH_PHASESTABSCAN || MDH_REFPHASESTABSCAN) {
        MDH_PATREFSCAN <- FALSE
        MDH_PATREFANDIMASCAN <- FALSE
        MDH_IMASCAN <- FALSE
      }
      
      if (MDH_PATREFSCAN && !(MDH_PATREFANDIMASCAN)) {
        MDH_IMASCAN <- FALSE
      }
      
      if (version == "vb") {
        if (!(MDH_SYNCDATA) && !(MDH_ACQEND)) {
          ulDMALength <- (2 * 4 * samples_in_scan + 128) * used_channels;
        }
      } else if (version == "vd") {
        if (!(MDH_SYNCDATA) && !(MDH_ACQEND) && (ulDMALength != 0)) {
          ulDMALength <- 192 + (2 * 4 * samples_in_scan + 32) * used_channels;
        }
      }
      
      if (MDH_IMASCAN) {
        # this chunk of data is from all coils
        inds <- c(inds, c(Lin, Ave, Sli, Par, Eco, Phs, Rep, Set, Seg, Ida, Idb,
                          Idc, Idd, Ide))
        
        ima_coils <- used_channels
        ima_samples <- samples_in_scan
        # if (Necho > ima_echoes) ima_echoes <- Necho
        ima_kspace_center_column <- kspace_center_column 
        # read in the data points
        for (x in 0:(used_channels - 1)) {
          if (version == "vb") {
            seek(con, 128 + cPos + x * (2 * 4 * samples_in_scan + 128), "start")
          } else if (version == "vd") {
            seek(con, 192 + 32 + cPos + x * (2 * 4 * samples_in_scan + 32),
                 "start")
          }
          # read (samples_in_scan * 2) floats
          fid <- readBin(con, "numeric", size = 4L, n = (samples_in_scan * 2),
                         endian = "little")
          
          raw_pt_end <- raw_pt_start + (samples_in_scan * 2) - 1
          raw_pts[raw_pt_start:raw_pt_end] <- fid
          expt <- length(raw_pts[raw_pt_start:raw_pt_end])
          actu <- length(fid)
          if (expt != actu) {
            print(expt)
            print(actu)
            stop("Unexpected number of data points")
          }
          raw_pt_start <- raw_pt_end + 1
        }
      }
      
      if (MDH_ACQEND || (ulDMALength == 0)) { # break out to the next scan
        if (scans < (Nscans - 1)) {
          cPos <- cPos + ulDMALength
          cPos <- cPos + 512 - (cPos %% 512)
          seek(con, cPos, "start")
          hdrLength <- read_uint32(con)
          cPos <- cPos + hdrLength
          seek(con, cPos, "start")
        }
        seek(con, cPos, "start")
        break
      }
      
      if (MDH_SYNCDATA) {
        cPos <- cPos + ulDMALength
        seek(con, cPos, "start")
        next
      }
      
      cPos <- cPos + ulDMALength
      seek(con, cPos, "start")
      
      if (seek(con) > fbytes) {
        stop("Read past the end of file.")
        break
      }
    }
  }
  close(con)
  
  if (verbose) cat("\n")
  raw_pts <- raw_pts[1:raw_pt_end]
  fid_offset <- floor(ima_kspace_center_column / 2) + 1
  dynamics <- length(raw_pts) / ima_coils / (ima_samples * 2)
  if (verbose) cat(paste("Raw data points :", length(raw_pts), "\n"))
  if (verbose) cat(paste("Coils           :", ima_coils, "\n"))
  if (verbose) cat(paste("Complex pts     :", ima_samples, "\n"))
  if (verbose) cat(paste("Dynamics        :", dynamics, "\n"))
  if (verbose) cat(paste("kspace center   :", ima_kspace_center_column, "\n"))
  #if (verbose) cat(paste("FID offset pts  :", fid_offset, "\n"))
  
  # make complex
  data <- raw_pts[c(TRUE, FALSE)] - 1i * raw_pts[c(FALSE, TRUE)]
  
  data <- array(data, dim = c(ima_samples, ima_coils, dynamics, 1, 1, 1, 1))
  data <- aperm(data, c(7,6,5,4,3,2,1))
   
  if (!full_fid & (floor(ima_kspace_center_column / 2) > 0)) {
    #data <- data[,,,,,,(fid_offset + 1):ima_samples, drop = FALSE]
    data <- data[,,,,,,(ima_kspace_center_column + 1):ima_samples, drop = FALSE]
  }
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)
  
  # get the resolution and geom info
  paras <- calc_siemens_paras(vars, FALSE)
  
  meta <- list(EchoTime = vars$te)

  mrs_data <- mrs_data(data = data, ft = vars$ft, resolution = paras$res,
                       ref = paras$ref, nuc = paras$nuc,
                       freq_domain = freq_domain, affine = paras$affine,
                       meta = meta, extra = extra)
  
  # some extra info specific to twix data
  mrs_data$twix_inds <- as.data.frame(matrix(inds, ncol = 14, byrow = TRUE))
  
  ind_names <- c("Lin", "Ave", "Sli", "Par", "Eco", "Phs", "Rep", "Set", "Seg",
                 "Ida", "Idb", "Idc", "Idd", "Ide")
  
  names(mrs_data$twix_inds) <- ind_names
  
  # remove SVS ref scans if required
  if (omit_svs_ref_scans & (max(mrs_data$twix_inds$Phs) == 1) &
      max(mrs_data$twix_inds$Ave) > 0) {
    
    mrs_data <- get_dyns(mrs_data, mrs_data$twix_inds$Phs == 1)
  }
  
  # correct the voxel dimensions for MRSI
  # this is needed because the Siemens software automatically zero-pads k-space
  # to a power of two, but we don't want this when dealing with raw twix data
  x_pts <- max(mrs_data$twix_inds$Seg) + 1
  y_pts <- max(mrs_data$twix_inds$Lin) + 1
  
  if (x_pts > 1 || y_pts > 1) {
    mrs_data$resolution[2] <- vars$x_dim / x_pts
    mrs_data$resolution[3] <- vars$y_dim / y_pts
    
    # fix the affine
    mrs_data$affine[,1] <- mrs_data$affine[,1] * vars$x_dim / x_pts
    mrs_data$affine[,2] <- mrs_data$affine[,2] * vars$y_dim / y_pts
  }
  
  return(mrs_data)
}

#' Read the text format header found in Siemens IMA and TWIX data files.
#' @param input file name to read or raw data.
#' @param version software version, can be "vb" or "vd".
#' @param verbose print information to the console.
#' @return a list of parameter values
#' @export
read_siemens_txt_hdr <- function(input, version = "vd", verbose) {
  if (is.character(input)) {
    con <- file(input, 'rb', encoding = "UTF-8")
  } else {
    # assume binary
    con <- rawConnection(input, "rb")
  }
  while (TRUE) {
    line <- readLines(con, n = 1, skipNul = TRUE, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "ulVersion") && version == 'vb') break

    if (startsWith(line, "ulVersion") && version == 'vd') {
      line_no_sp <- gsub(" ", "", line)
      line_no_sp_no_tab <- gsub("\t", "", line_no_sp)
      # skip if there isn't an equal sign once spaces have been removed
      if (!startsWith(line_no_sp_no_tab, "ulVersion=")) next
      tSequenceFilename <- readLines(con, n = 1)
      tProtocolName <- readLines(con, n = 1)
      last_ulVersion_pos <- seek(con)
      #if ((tProtocolName != "tProtocolName\t = \t\"AdjCoilSens\"") && (tProtocolName != "tProtocolName\t = \t\"CBU_MPRAGE_32chn\"")) {
        #print(tProtocolName)
      #  break
      #}
    }
  }
  
  # TODO reads the whole file looking for "ulVersion" and then
  # tracks back to the last one. Must be a better way of finding the
  # last one without reading all the spectral data. This is particularly
  # bad for ima data when there will only ever be one set.
  if (version == 'vd') seek(con, where = last_ulVersion_pos)
  
  vars <- list(averages = NA,
               fs = NA,
               ft = NA,
               te = NA,
               N = NA,
               x_pts = 1,
               y_pts = 1,
               z_pts = 1,
               z_dim = NA,
               x_dim = NA,
               y_dim = NA,
               ip_rot = 0,
               pos_sag = 0,
               pos_cor = 0,
               pos_tra = 0,
               norm_sag = 0,
               norm_cor = 0,
               norm_tra = 0)
  
  # when a parameter is missing from an ima file it means it's zero (I think)
  slice_dPhaseFOV    <- 0
  slice_dReadoutFOV  <- 0
  slice_dThickness   <- 0
  slice_dSag         <- 0
  slice_dCor         <- 0
  slice_dTra         <- 0
  slice_norm_sag     <- 0
  slice_norm_cor     <- 0
  slice_norm_tra     <- 0
  slice_ip_rot       <- 0
  voi_dPhaseFOV      <- 0
  voi_dReadoutFOV    <- 0
  voi_dThickness     <- 0
  voi_dSag           <- 0
  voi_dCor           <- 0
  voi_dTra           <- 0
  voi_norm_sag       <- 0
  voi_norm_cor       <- 0
  voi_norm_tra       <- 0
  voi_ip_rot         <- 0
  scan_reg_pos_tra   <- 0
  
  while (TRUE) {
    line <- readLines(con, n = 1, skipNul = TRUE)
    if (grepl("### ASCCONV END ###", line, fixed = TRUE, useBytes = TRUE)) {
      break
    } else if (startsWith(line, "lAverages")) {
      vars$averages <- as.integer(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sRXSPEC.alDwellTime[0]")) {
      vars$fs <- 1e9 / (as.numeric(strsplit(line, "=")[[1]][2]))
    } else if (startsWith(line, "sTXSPEC.asNucleusInfo[0].lFrequency")) {
      vars$ft <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "alTE[0]")) {
      vars$te <- (as.numeric(strsplit(line, "=")[[1]][2])) / 1e6
    } else if (startsWith(line, "sSpecPara.lVectorSize")) {
      vars$N <- as.integer(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.lFinalMatrixSizePhase")) {
      vars$x_pts <- as.integer(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.lFinalMatrixSizeRead")) {
      vars$y_pts <- as.integer(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.lFinalMatrixSizeSlice")) {
      vars$z_pts <- as.integer(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSliceArray.asSlice[0].dThickness")) {
      slice_dThickness <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSliceArray.asSlice[0].dPhaseFOV")) {
      slice_dPhaseFOV <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSliceArray.asSlice[0].dReadoutFOV")) {
      slice_dReadoutFOV <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.sVoI.dThickness")) {
      voi_dThickness <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.sVoI.dPhaseFOV")) {
      voi_dPhaseFOV <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.sVoI.dReadoutFOV")) {
      voi_dReadoutFOV <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.sVoI.dInPlaneRot")) {
      voi_ip_rot <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.sVoI.sPosition.dSag")) {
      voi_dSag <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.sVoI.sPosition.dCor")) {
      voi_dCor <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.sVoI.sPosition.dTra")) {
      voi_dTra <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "lScanRegionPosTra")) {
      scan_reg_pos_tra <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.sVoI.sNormal.dSag")) {
      voi_norm_sag <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.sVoI.sNormal.dCor")) {
      voi_norm_cor <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.sVoI.sNormal.dTra")) {
      voi_norm_tra <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSliceArray.asSlice[0].dInPlaneRot")) {
      slice_ip_rot <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSliceArray.asSlice[0].sPosition.dSag")) {
      slice_dSag <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSliceArray.asSlice[0].sPosition.dCor")) {
      slice_dCor <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSliceArray.asSlice[0].sPosition.dTra")) {
      slice_dTra <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSliceArray.asSlice[0].sNormal.dSag")) {
      slice_norm_sag <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSliceArray.asSlice[0].sNormal.dCor")) {
      slice_norm_cor <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSliceArray.asSlice[0].sNormal.dTra")) {
      slice_norm_tra <- as.numeric(strsplit(line, "=")[[1]][2])
    } else if (startsWith(line, "sSpecPara.ucRemoveOversampling")) {
      vars$rm_oversampling <- as.numeric(strsplit(line, "=")[[1]][2])
    }
  }
  
  if (verbose) cat(paste("Table position  :", scan_reg_pos_tra, "mm\n"))
 
  # how many voxels do we expect?
  Nvoxels <- vars$x_pts * vars$y_pts * vars$z_pts
  
  if (Nvoxels > 1) {
    # looks like MRSI
    vars$x_dim    <- slice_dReadoutFOV 
    vars$y_dim    <- slice_dPhaseFOV 
    vars$z_dim    <- slice_dThickness
    vars$pos_sag  <- slice_dSag 
    vars$pos_cor  <- slice_dCor
    vars$pos_tra  <- slice_dTra + scan_reg_pos_tra # don't ask
    vars$norm_sag <- slice_norm_sag
    vars$norm_cor <- slice_norm_cor
    vars$norm_tra <- slice_norm_tra
    vars$ip_rot   <- slice_ip_rot
  } else if (Nvoxels == 1) {
    # looks like SVS
    vars$x_dim    <- voi_dReadoutFOV 
    vars$y_dim    <- voi_dPhaseFOV 
    vars$z_dim    <- voi_dThickness
    vars$pos_sag  <- voi_dSag 
    vars$pos_cor  <- voi_dCor
    vars$pos_tra  <- voi_dTra + scan_reg_pos_tra # don't ask
    vars$norm_sag <- voi_norm_sag
    vars$norm_cor <- voi_norm_cor
    vars$norm_tra <- voi_norm_tra
    vars$ip_rot   <- voi_ip_rot
  } else {
    stop("Unexpected number of voxels found.")
  }
  
  close(con)
  vars
}

read_char <- function(con) {
  readBin(con, "character", size = 1L, n = 1, endian = "little")
}

read_float32 <- function(con) {
  readBin(con, "numeric", size = 4L, n = 1, endian = "little")
}

read_int16 <- function(con) {
  readBin(con, "integer", size = 2L, n = 1, endian = "little")
}

read_int32 <- function(con) {
  readBin(con, "integer", size = 4L, n = 1, endian = "little")
}

read_uint16 <- function(con) {
  readBin(con, "integer", size = 2L, n = 1, endian = "little", signed = FALSE)
}

read_uint32 <- function(con) {
  int <- readBin(con, "integer", size = 4L, n = 1, endian = "little")
  raw_bits <- intToBits(int)
  # warning - R doesn't nativly support unsigned 32 bit integers - so this will
  # be an approximation by converting to double
  sum(2^.subset(0:31, as.logical(raw_bits)))
}

read_uint64 <- function(con) {
  intvec <- readBin(con, "integer", size = 4L, n = 2, endian = "little")
  raw_bits <- c(intToBits(intvec[1]), intToBits(intvec[2]))
  # warning - R doesn't nativly support 64 bit integers - so this will be 
  # an approximation by converting to double
  sum(2^.subset(0:63, as.logical(raw_bits)))
}