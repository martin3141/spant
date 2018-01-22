#' @export
read_twix <- function(fname) {
  # read the binary header
  con <- file(fname, "rb")
  first_int <- read_uint32(con)
  second_int <- read_uint32(con)
  #print(first_int)
  #print(second_int)
  if ((first_int < 10000) && (second_int <= 64)) {
    cat("Twix file is VD format.\n")
    version <- "vd"
    Nscans <- second_int
    measID <- read_uint32(con)
    #print(measID)
    fileID <- read_uint32(con)
    #print(fileID)
    measOffset <- read_uint64(con) # usually 10240 bytes
    #print(measOffset)
    measLength <- read_uint64(con)
    hdrLength <- read_uint32(con)
    dataStart <- measOffset + hdrLength
    #print(dataStart)
  } else {
    cat("Twix file is VB format.\n")
    version <- "vb"
    dataStart <- first_int
    Nscans <- 1
  }
  close(con)
  vars <- read_twix_txt_hdr(fname, version)
}

read_twix_txt_hdr <- function(fname, version) {
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

read_float32 <- function(con) {
  readBin(con, "integer", size = 4L, n = 1, endian = "little")
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