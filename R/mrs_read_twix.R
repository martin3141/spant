read_twix <- function(fname, verbose) {
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
  vars <- read_siemens_txt_hdr(fname, version)
  
  # read data points 
  con <- file(fname, "rb")
  seek(con, dataStart)
  
  ima_echoes <- 0
  cPos <- dataStart
  
  # overestimate the number of data points and pre-allocate
  expected_pts <- as.integer((fbytes - dataStart) / 4)  
  raw_pts <- c(NA)
  length(raw_pts) <- expected_pts
  raw_pt_start <- 1
  
  for (scans in 0:(Nscans - 1)) {
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
      samples_in_scan       <- read_uint16(con)
      used_channels         <- read_uint16(con)
      seek(con, 8, "current")
      Necho                 <- read_uint16(con)
      seek(con, 22, "current")
      kspace_center_column  <- read_uint16(con)
      
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
        ima_coils <- used_channels
        ima_samples <- samples_in_scan
        if (Necho > ima_echoes) ima_echoes <- Necho
        ima_kspace_center_column <- kspace_center_column 
        # read in the data points
        for (x in 0:(used_channels - 1)) {
          if (version == "vb") {
            seek(con, 128 + cPos + x * (2 * 4 * samples_in_scan + 128), "start")
          } else if (version == "vd") {
            seek(con, 192 + 32 + cPos + x * (2 * 4 * samples_in_scan + 32), "start")
          }
          # read (samples_in_scan * 2) floats
          fid <- readBin(con, "numeric", size = 4L, n = (samples_in_scan * 2), endian = "little")
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
  if (verbose) cat(paste("Coils          :", ima_coils, "\n"))
  if (verbose) cat(paste("Complex pts    :", ima_samples, "\n"))
  if (verbose) cat(paste("Dynamics       :", dynamics, "\n"))
  if (verbose) cat(paste("FID offset pts :", fid_offset, "\n"))
  
  # make complex
  data <- raw_pts[c(TRUE, FALSE)] - 1i * raw_pts[c(FALSE, TRUE)]
  
  data <- array(data, dim = c(ima_samples, ima_coils, dynamics, 1, 1, 1, 1))
  data <- aperm(data, c(7,6,5,4,3,2,1))
   
  if (floor(ima_kspace_center_column / 2) > 0) {
    data <- data[,,,,,,(fid_offset + 1):ima_samples, drop = FALSE]
  }
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / vars$fs)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  ref <- def_acq_paras()$ref
  
  mrs_data <- list(ft = vars$ft, data = data, resolution = res,
                   te = vars$te, ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}

read_siemens_txt_hdr <- function(fname, version) {
  con <- file(fname, 'rb', encoding = "UTF-8")
  while (TRUE) {
    line <- readLines(con, n = 1 ,skipNul = TRUE)
    if (startsWith(line, "ulVersion") && version == 'vb') {
      break
    }
    
    if (startsWith(line, "ulVersion") && version == 'vd') {
      tSequenceFilename <- readLines(con, n = 1)
      tProtocolName <- readLines(con, n = 1)
      if ((tProtocolName != "tProtocolName\t = \t\"AdjCoilSens\"") && (tProtocolName != "tProtocolName\t = \t\"CBU_MPRAGE_32chn\"")) {
        #print(tProtocolName)
        break
      }
    }
  }
  
  vars <- vector(mode = "list", length = 5)
  names(vars) <- c("averages", "fs", "ft", "te", "N")
  
  while (TRUE) {
    line <- readLines(con, n = 1, skipNul = TRUE)
    if (grepl("### ASCCONV END ###", line, fixed = TRUE)) {
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
    }
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