# stolen from oro.dicom to convert unsigned int tags to character arrays
rawToHex <- function(bytes) toupper(paste(rev(bytes), collapse = ""))

read_raw <- function(fname, endian = "little") {
  # read the full file as raw
  fsize <- file.info(fname)$size
  return(readBin(fname, "raw", n = as.integer(fsize), endian = endian))
}

#' A very simple DICOM reader.
#' 
#' Note this reader is very basic and does not use a DICOM dictionary or try to
#' convert the data to the correct datatype. For a more robust and sophisticated
#' reader use the oro.dicom package.
#' @param input either a file path or raw binary object.
#' @param tags a named list of tags to be extracted from the file.
#' eg tags <- list(spec_data = "7FE1,1010", pat_name = "0010,0010")
#' @param endian can be "little" or "big".
#' @param debug print out some debugging information, can be "little" or "big".
#' @return a list with the same structure as the input, but with tag codes
#' replaced with the corresponding data in a raw format.
#' @export
dicom_reader <- function(input, tags = list(sop_class_uid = "0008,0016"),
                         endian = "little", debug = FALSE) {
  
  # assume as file path if input is character, otherwise assume raw
  if (is.character(input)) {
    fsize <- file.info(input)$size
    fraw  <- readBin(input, "raw", n = as.integer(fsize), endian = endian)
  } else {
    fraw <- input 
  }
  
  # skip the 128 byte preamble and read the next 4 bytes (should be "DICM")
  dicom_prefix <- rawToChar(fraw[129:132])
  
  if (dicom_prefix != "DICM") {
    stop("Could not find DICM prefix - input is probally not DICOM. Try manually setting the correct 'format' argument of read_mrs.")
  }
  
  # strip the preamble
  fraw <- fraw[133:length(fraw)]
  
  # prepare the output structure
  out <- tags
  out[] <- list(NULL)
  
  short_vrs <- c("OB","OW","SQ","UN","UT")
  long_vrs  <- c("AE","AS","AT","CS","DA","DS","DT","FL","FD","IS","LO",
                 "LT","OF","PN","SH","SL","SS","ST","TM","UI","UL","US")
  
  # for implicit dicom we don't know when we have hit a sequence without
  # prior knowledge - this is a list that should work for Philips DICOM MRS
  seq_tags <- c("5200,9229", "0020,9116", "5200,9230", "0018,9114",
                "0020,9113", "2005,140F", "0028,9110", "0018,9103")
  
  # current position in the data
  pos <- 0
  
  if (debug) debug_table <- NULL
  
  while (pos < length(fraw)) {
    
    # read the tag ID
    group   <- rawToHex(fraw[pos + 1:2])
    element <- rawToHex(fraw[pos + 3:4])
    tag_str <- paste(group, element, sep = ",") 
    
    # tentatively read in the VR - only used for explicit VRs
    vr <- try(rawToChar(fraw[pos + 5:6]), TRUE)
    
    # MRS is a special case where we need to check for zero padding
    # 5600,0020 - standard dicom MRS data tag
    # 2005,1270 - Philips private dicom MRS data tag
    if (tag_str %in% c("5600,0020", "2005,1270", "0065,FF07", "0065,FF08")) {
      # check for zero padding
      zp <- is.na(pos[7]) & is.na(pos[8])
      if ((vr %in% c(long_vrs, short_vrs)) & zp) {
        # explicit VR with two bytes of zero padding
        length_raw <- fraw[pos + 9:12]
        length <- read_uint32(length_raw)
        # move pos to start of the value
        pos <- pos + 12
      } else if (vr %in% c(long_vrs, short_vrs)) {
        # explicit VR without zero padding
        length_raw <- fraw[pos + 7:8]
        length <- readBin(length_raw, "integer", size = 2, signed = FALSE)
        pos <- pos + 8
      } else {
        # looks like implicit VR
        length_raw <- fraw[pos + 5:8]
        length <- read_uint32(length_raw)
        pos <- pos + 8
      }
    } else { 
      if (vr %in% c(short_vrs)) {
        # explicit VR with two bytes of zero padding
        length_raw <- fraw[pos + 9:12]
        length <- read_uint32(length_raw)
        # move pos to start of the value
        pos <- pos + 12
      } else if (vr %in% long_vrs) {
        # explicit VR without zero padding
        length_raw <- fraw[pos + 7:8]
        length <- readBin(length_raw, "integer", size = 2, signed = FALSE)
        pos <- pos + 8
      } else {
        # looks like implicit VR
        length_raw <- fraw[pos + 5:8]
        length <- read_uint32(length_raw)
        pos <- pos + 8
      }
    }
    
    if (rawToHex(length_raw) %in% c("FFFFFFFF", "00000A00")) length <- 0 
    
    # don't skip over sequences
    if (vr %in% c("SQ")) length <- 0
    
    # don't skip over items
    if (tag_str %in% c("FFFE,E000")) length <- 0
    
    # only needed for implicit DICOM
    if (tag_str %in% seq_tags) length <- 0
      
    if (length == 0) {
      start_byte <- pos + 1
      end_byte   <- NA
    } else {
      start_byte <- pos + 1
      end_byte   <- pos + length
    }
  
    if (debug) {
      tag <- paste("(", tag_str, ")", sep = "")
      x <- data.frame(tag, start_byte, end_byte, length, vr)
      debug_table <- rbind(debug_table, x)
    }
    
    # check if current tag is the in list of those to be returned
    search_res <- grep(tag_str, tags) 
    
    if (length(search_res) == 1) out[[search_res]] <- fraw[start_byte:end_byte]
    
    pos <- pos + length
  }
  
  if (debug) print(debug_table)
  
  return(out)
}